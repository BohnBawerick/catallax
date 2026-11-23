#!/usr/bin/env python3
import os
import subprocess
import sys


toolchain = os.environ.get("RUST_TOOLCHAIN", "stable")
step_name = os.environ.get("RUST_STEP_NAME", "Set up Rust toolchain")
timeout_seconds = int(os.environ.get("RUSTUP_INSTALL_TIMEOUT_SECONDS", "120"))
attempts = int(os.environ.get("RUSTUP_INSTALL_ATTEMPTS", "3"))
components = [item for item in os.environ.get("RUST_COMPONENTS", "").split(",") if item]
targets = [item for item in os.environ.get("RUST_TARGETS", "").split(",") if item]


def output(args):
    return subprocess.run(args, check=True, text=True, capture_output=True).stdout.splitlines()


def installed():
    try:
        toolchains = output(["rustup", "toolchain", "list"])
        if not any(line.split() and line.split()[0] == toolchain for line in toolchains):
            return False
        installed_components = output(["rustup", "component", "list", "--toolchain", toolchain, "--installed"])
        if any(not any(line.split() and line.split()[0].startswith(f"{component}-")
                       for line in installed_components) for component in components):
            return False
        installed_targets = output(["rustup", "target", "list", "--toolchain", toolchain, "--installed"])
    except subprocess.CalledProcessError:
        return False
    return all(target in installed_targets for target in targets)


if installed():
    print(f"{step_name}: reusing installed Rust toolchain {toolchain}")
else:
    command = ["rustup", "toolchain", "install", toolchain, "--profile", "minimal",
               "--no-self-update", "--force-non-host"]
    for component in components:
        command.extend(["--component", component])
    for target in targets:
        command.extend(["--target", target])
    succeeded = False
    for attempt in range(1, attempts + 1):
        print(f"{step_name}: installing Rust toolchain {toolchain} "
              f"(attempt {attempt}/{attempts}, timeout {timeout_seconds}s)", flush=True)
        try:
            subprocess.run(command, check=True, timeout=timeout_seconds)
            succeeded = True
            break
        except subprocess.TimeoutExpired:
            print(f"{step_name}: installation timed out after {timeout_seconds}s", file=sys.stderr)
        except subprocess.CalledProcessError as error:
            print(f"{step_name}: rustup exited with status {error.returncode}", file=sys.stderr)
    if not succeeded:
        print(f"::error::{step_name}: failed to install Rust toolchain {toolchain} after {attempts} "
              f"attempts (each limited to {timeout_seconds}s)", file=sys.stderr)
        sys.exit(1)

github_env = os.environ.get("GITHUB_ENV")
if github_env:
    with open(github_env, "a", encoding="utf-8") as stream:
        stream.write(f"RUSTUP_TOOLCHAIN={toolchain}\n")
github_path = os.environ.get("GITHUB_PATH")
if github_path:
    cargo = subprocess.run(["rustup", "which", "--toolchain", toolchain, "cargo"], check=True,
                           text=True, capture_output=True).stdout.strip()
    with open(github_path, "a", encoding="utf-8") as stream:
        stream.write(f"{os.path.dirname(cargo)}\n")
