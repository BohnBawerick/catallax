import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from textwrap import dedent, indent
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
SETUP = ROOT / ".github/actions/rust-toolchain/setup.py"
ROUNDED_STEP = ROOT / ".github/workflows/rounded-step-speed.yml"


class RustToolchainSetupTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix="catallax-rustup-")
        self.addCleanup(directory.cleanup)
        self.base = Path(directory.name)
        self.github_env = self.base / "setup-env"
        self.github_path = self.base / "setup-path"

    def run_setup(self, *, toolchains, components="", targets="", install="", requested="rustfmt",
                  requested_targets="", stalled_probe="", idle_seconds="4", max_seconds="15"):
        base = self.base
        rustup = base / "rustup"
        log = base / "rustup.log"
        rustup.write_text(dedent("""\
            #!/usr/bin/env python3
            import os, sys, time
            from pathlib import Path

            args = sys.argv[1:]
            command = ' '.join(args[:2])
            with open(os.environ['FAKE_RUSTUP_LOG'], 'a') as stream:
                stream.write(' '.join(args) + '\\n')
            if command == os.environ['RUSTUP_STALLED_PROBE']:
                time.sleep(6)
            if command == 'toolchain list':
                print(os.environ['RUSTUP_TOOLCHAINS'])
            elif command == 'component list':
                print(os.environ['RUSTUP_COMPONENTS'])
            elif command == 'target list':
                print(os.environ['RUSTUP_TARGETS'])
            elif command == 'which --toolchain':
                print('/fake/.rustup/toolchains/stable/bin/cargo')
            elif command != 'toolchain install':
                sys.exit(2)
            else:
            """) + indent(dedent(install).strip() or "pass", "    ") + "\n")
        rustup.chmod(0o755)
        env = {
            **os.environ,
            "GITHUB_ENV": str(self.github_env),
            "GITHUB_PATH": str(self.github_path),
            "PATH": f"{base}{os.pathsep}{os.environ['PATH']}",
            "FAKE_RUSTUP_LOG": str(log),
            "RUSTUP_HOME": str(base / "rustup-home"),
            "RUSTUP_TOOLCHAINS": toolchains,
            "RUSTUP_COMPONENTS": components,
            "RUSTUP_TARGETS": targets,
            "RUSTUP_INSTALL_TIMEOUT_SECONDS": idle_seconds,
            "RUSTUP_INSTALL_MAX_SECONDS": max_seconds,
            "RUSTUP_INSTALL_ATTEMPTS": "2",
            "RUSTUP_PROBE_TIMEOUT_SECONDS": "3",
            "RUSTUP_STALLED_PROBE": stalled_probe,
            "RUST_STEP_NAME": "Test Rust setup",
            "RUST_TOOLCHAIN": "stable",
            "RUST_COMPONENTS": requested,
            "RUST_TARGETS": requested_targets,
        }
        log.write_text("")
        result = subprocess.run(["python3", str(SETUP)], env=env, capture_output=True, text=True, timeout=15)
        return result, log.read_text()

    def test_setup_writes_only_to_its_temporary_github_command_files(self):
        parent_env = self.base / "parent-env"
        parent_path = self.base / "parent-path"
        parent_env.write_text("PARENT=value\n")
        parent_path.write_text("/parent/bin\n")
        with patch.dict(os.environ, GITHUB_ENV=str(parent_env), GITHUB_PATH=str(parent_path)):
            result, _ = self.run_setup(toolchains="stable",
                                       components="rustfmt-x86_64-unknown-linux-gnu")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(parent_env.read_text(), "PARENT=value\n")
        self.assertEqual(parent_path.read_text(), "/parent/bin\n")
        self.assertEqual(self.github_env.read_text(), "RUSTUP_TOOLCHAIN=stable\n")
        self.assertEqual(self.github_path.read_text(), "/fake/.rustup/toolchains/stable/bin\n")

    def test_reuses_a_host_qualified_toolchain_without_installing(self):
        result, log = self.run_setup(
            toolchains="stable-x86_64-unknown-linux-gnu (default)",
            components="rustfmt-x86_64-unknown-linux-gnu (installed)",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("reusing installed Rust toolchain stable", result.stdout)
        self.assertNotIn("toolchain install", log)

    def test_branch_check_trims_component_names_before_install(self):
        workflow = (ROOT / ".github/workflows/branch-check.yml").read_text()
        requested = re.search(r"^\s+components: (.+)$", workflow, re.MULTILINE)[1]
        result, log = self.run_setup(
            toolchains="", requested=requested,
            install="assert args[args.index('--component'):] == ['--component', 'rustfmt', '--component', 'clippy']",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--component rustfmt --component clippy", log)

    def test_every_rustup_probe_has_a_deadline(self):
        for probe in ("toolchain list", "component list", "target list", "which --toolchain"):
            with self.subTest(probe=probe):
                result, log = self.run_setup(toolchains="stable", requested="",
                                             requested_targets="x86_64-unknown-linux-gnu",
                                             targets="x86_64-unknown-linux-gnu", stalled_probe=probe)
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertIn(f"rustup {probe}", result.stderr)
                self.assertIn("timed out after 3s", result.stderr)
                self.assertNotIn("toolchain install", log)
                self.assertFalse(self.github_env.exists())
                self.assertFalse(self.github_path.exists())

    def test_install_continues_while_rustup_reports_progress(self):
        result, log = self.run_setup(
            toolchains="", idle_seconds="2",
            install="""
                for i in range(14):
                    print(f'downloaded {i}', end='\\r', file=sys.stderr, flush=True)
                    time.sleep(0.2)
            """,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("downloaded 13", result.stdout)
        self.assertEqual(log.count("toolchain install"), 1)

    def test_install_continues_while_a_silent_download_grows(self):
        result, log = self.run_setup(
            toolchains="", idle_seconds="2",
            install="""
                downloads = Path(os.environ['RUSTUP_HOME']) / 'downloads'
                downloads.mkdir(parents=True)
                with (downloads / 'archive.partial').open('ab', buffering=0) as archive:
                    for _ in range(14):
                        archive.write(b'x')
                        time.sleep(0.2)
            """,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(log.count("toolchain install"), 1)

    def test_install_has_an_overall_deadline_even_with_progress(self):
        result, log = self.run_setup(
            toolchains="", idle_seconds="2", max_seconds="3",
            install="""
                for i in range(25):
                    print(i, flush=True)
                    time.sleep(0.2)
            """,
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn("overall install limit of 3s", result.stderr)
        self.assertEqual(log.count("toolchain install"), 1)

    def test_stalled_download_times_out_retries_and_names_the_failed_step(self):
        result, log = self.run_setup(toolchains="", install="time.sleep(5)", idle_seconds="1")
        self.assertEqual(result.returncode, 1)
        self.assertEqual(log.count("toolchain install"), 2)
        self.assertEqual(result.stderr.count("made no install progress for 1s"), 2)
        self.assertIn("Test Rust setup: failed to install Rust toolchain stable after 2 attempts", result.stderr)


class RustToolchainWorkflowTests(unittest.TestCase):
    def test_rounded_step_jobs_load_the_action_from_the_sparse_workflow_checkout(self):
        workflow = ROUNDED_STEP.read_text()
        for job in ("paired", "edits"):
            with self.subTest(job=job), tempfile.TemporaryDirectory() as directory:
                body = re.search(rf"^  {job}:\n((?:\n| {{4,}}[^\n]*\n)+)", workflow, re.MULTILINE)[1]
                checkout = re.search(r"  path: workflow\n\s+sparse-checkout: \|\n((?: {12}\S[^\n]*\n)+)", body)
                sparse_paths = checkout[1].split()
                action = re.search(r"uses: (\./\S+/rust-toolchain)", body)[1]
                checkout_dir = Path(directory) / "workflow"
                action_source = SETUP.parent.relative_to(ROOT)
                shutil.copytree(SETUP.parent, checkout_dir / action_source)
                git = ["git", "-C", str(checkout_dir)]
                subprocess.run([*git, "init", "-q"], check=True)
                subprocess.run([*git, "add", str(action_source)], check=True)
                subprocess.run([*git, "-c", "user.name=Test", "-c", "user.email=test@example.invalid",
                                "-c", "commit.gpgsign=false", "commit", "-qm", "Fixture"], check=True)
                subprocess.run([*git, "sparse-checkout", "set", "--no-cone", *sparse_paths], check=True)
                action_dir = Path(directory) / action
                for name in ("action.yml", "setup.py"):
                    self.assertTrue((action_dir / name).is_file(), f"{action}/{name} is missing")


if __name__ == "__main__":
    unittest.main()
