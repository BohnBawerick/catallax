#!/usr/bin/env python3
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time


toolchain = os.environ.get("RUST_TOOLCHAIN", "stable")
step_name = os.environ.get("RUST_STEP_NAME", "Set up Rust toolchain")
idle_seconds = int(os.environ.get("RUSTUP_INSTALL_TIMEOUT_SECONDS", "120"))
max_install_seconds = int(os.environ.get("RUSTUP_INSTALL_MAX_SECONDS", "900"))
probe_seconds = int(os.environ.get("RUSTUP_PROBE_TIMEOUT_SECONDS", "30"))
attempts = int(os.environ.get("RUSTUP_INSTALL_ATTEMPTS", "3"))


def requested_items(variable):
    return [name for item in os.environ.get(variable, "").split(",") if (name := item.strip())]


components = requested_items("RUST_COMPONENTS")
targets = requested_items("RUST_TARGETS")


def output(args):
    # Keep downloads in the explicitly monitored install command.
    env = {**os.environ, "RUSTUP_AUTO_INSTALL": "0"}
    try:
        return subprocess.run(args, check=True, text=True, stdout=subprocess.PIPE,
                              env=env, timeout=probe_seconds).stdout.strip()
    except subprocess.TimeoutExpired:
        sys.exit(f"::error::{step_name}: {' '.join(args)} timed out after {probe_seconds}s")


def installed():
    try:
        listed = output(["rustup", "toolchain", "list"]).splitlines()
        names = [line.split()[0] for line in listed if line.strip()]
        if not any(name == toolchain or name.startswith(f"{toolchain}-") for name in names):
            return False
        # Rustup resolves the requested name to its configured host for this query.
        installed_components = output(["rustup", "component", "list", "--toolchain", toolchain,
                                       "--installed"]).splitlines()
        installed_components = [line.split()[0] for line in installed_components if line.strip()]
        if any(not any(name == component or name.startswith(f"{component}-")
                       for name in installed_components) for component in components):
            return False
        installed_targets = output(["rustup", "target", "list", "--toolchain", toolchain,
                                    "--installed"]).splitlines()
    except subprocess.CalledProcessError as error:
        print(f"{step_name}: {' '.join(error.cmd)} exited with status {error.returncode}", file=sys.stderr)
        return False
    return all(target in installed_targets for target in targets)


def partial_download_sizes(downloads):
    sizes = {}
    for archive in downloads.glob("*.partial"):
        try:
            sizes[archive] = archive.stat().st_size
        except FileNotFoundError:
            # Rustup renames a complete archive while the watchdog is polling.
            pass
    return sizes


def install_attempt(command, deadline):
    home = os.environ.get("RUSTUP_HOME")
    downloads = (Path(home) if home else Path.home() / ".rustup") / "downloads"
    sizes = partial_download_sizes(downloads)
    # A spinner can redraw while Rustup is stalled. Watch the archive bytes instead.
    env = {**os.environ, "RUSTUP_TERM_PROGRESS_WHEN": "never"}
    with tempfile.TemporaryDirectory(prefix="rustup-install-") as directory:
        log = Path(directory) / "rustup.log"
        # Separate file handles let us stream output on Windows as well as Unix.
        with log.open("wb") as writer, log.open("rb") as reader, \
                subprocess.Popen(command, stdout=writer, stderr=subprocess.STDOUT, env=env) as process:
            last_progress = time.monotonic()
            try:
                while process.poll() is None:
                    now = time.monotonic()
                    if now >= deadline:
                        return f"reached the overall install limit of {max_install_seconds}s"
                    chunk = reader.read(65536)
                    if chunk:
                        sys.stdout.buffer.write(chunk)
                        sys.stdout.buffer.flush()
                    current_sizes = partial_download_sizes(downloads)
                    if chunk or any(size > sizes.get(path, 0) for path, size in current_sizes.items()):
                        last_progress = now
                    sizes = current_sizes
                    remaining = min(deadline, last_progress + idle_seconds) - now
                    if remaining <= 0:
                        return f"made no install progress for {idle_seconds}s"
                    time.sleep(min(0.1, remaining))
            finally:
                if process.poll() is None:
                    process.kill()
                    process.wait()
                while chunk := reader.read(65536):
                    sys.stdout.buffer.write(chunk)
                sys.stdout.buffer.flush()
            if process.returncode:
                return f"rustup exited with status {process.returncode}"


def install():
    command = ["rustup", "toolchain", "install", toolchain, "--profile", "minimal",
               "--no-self-update", "--force-non-host"]
    for component in components:
        command.extend(["--component", component])
    for target in targets:
        command.extend(["--target", target])
    deadline = time.monotonic() + max_install_seconds
    attempt = 0
    while attempt < attempts and time.monotonic() < deadline:
        attempt += 1
        print(f"{step_name}: installing Rust toolchain {toolchain} "
              f"(attempt {attempt}/{attempts}, idle limit {idle_seconds}s)", flush=True)
        error = install_attempt(command, deadline)
        if error is None:
            return
        print(f"{step_name}: {error}", file=sys.stderr)
    sys.exit(f"::error::{step_name}: failed to install Rust toolchain {toolchain} after {attempt} "
             f"attempt{'s' if attempt != 1 else ''}")


if installed():
    print(f"{step_name}: reusing installed Rust toolchain {toolchain}")
else:
    install()

github_env = os.environ.get("GITHUB_ENV")
github_path = os.environ.get("GITHUB_PATH")
if github_path:
    try:
        cargo = output(["rustup", "which", "--toolchain", toolchain, "cargo"])
    except subprocess.CalledProcessError as error:
        sys.exit(f"::error::{step_name}: {' '.join(error.cmd)} exited with status {error.returncode}")
    with open(github_path, "a", encoding="utf-8") as stream:
        stream.write(f"{os.path.dirname(cargo)}\n")
if github_env:
    with open(github_env, "a", encoding="utf-8") as stream:
        stream.write(f"RUSTUP_TOOLCHAIN={toolchain}\n")
