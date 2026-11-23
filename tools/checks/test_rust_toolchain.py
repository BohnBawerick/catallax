import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SETUP = ROOT / ".github/actions/rust-toolchain/setup.py"


class RustToolchainSetupTests(unittest.TestCase):
    def run_setup(self, *, toolchains, components="", targets="", install="", timeout="1"):
        directory = tempfile.TemporaryDirectory(prefix="catallax-rustup-")
        self.addCleanup(directory.cleanup)
        base = Path(directory.name)
        rustup = base / "rustup"
        log = base / "rustup.log"
        rustup.write_text(
            "#!/usr/bin/env python3\n"
            "import os, sys, time\n"
            "with open(os.environ['RUSTUP_LOG'], 'a') as stream: stream.write(' '.join(sys.argv[1:]) + '\\n')\n"
            "args = sys.argv[1:]\n"
            "if args[:2] == ['toolchain', 'list']: print(os.environ['RUSTUP_TOOLCHAINS'])\n"
            "elif args[:2] == ['component', 'list']: print(os.environ['RUSTUP_COMPONENTS'])\n"
            "elif args[:2] == ['target', 'list']: print(os.environ['RUSTUP_TARGETS'])\n"
            "elif args[:2] == ['toolchain', 'install']: " + (install or "pass") + "\n"
            "elif args[:1] == ['which']: print('/fake/.rustup/toolchains/stable/bin/cargo')\n"
            "else: sys.exit(2)\n"
        )
        rustup.chmod(0o755)
        env = {
            **os.environ,
            "PATH": f"{base}{os.pathsep}{os.environ['PATH']}",
            "RUSTUP_LOG": str(log),
            "RUSTUP_TOOLCHAINS": toolchains,
            "RUSTUP_COMPONENTS": components,
            "RUSTUP_TARGETS": targets,
            "RUSTUP_INSTALL_TIMEOUT_SECONDS": timeout,
            "RUSTUP_INSTALL_ATTEMPTS": "2",
            "RUST_STEP_NAME": "Test Rust setup",
            "RUST_TOOLCHAIN": "stable",
            "RUST_COMPONENTS": "rustfmt",
        }
        result = subprocess.run(["python3", str(SETUP)], env=env, capture_output=True, text=True)
        return result, log.read_text() if log.exists() else ""

    def test_reuses_a_complete_toolchain_without_installing(self):
        result, log = self.run_setup(
            toolchains="stable-x86_64-unknown-linux-gnu\nstable",
            components="rustfmt-x86_64-unknown-linux-gnu (installed)",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("reusing installed Rust toolchain stable", result.stdout)
        self.assertNotIn("toolchain install", log)

    def test_stalled_download_times_out_retries_and_names_the_failed_step(self):
        result, log = self.run_setup(toolchains="", install="time.sleep(5)")
        self.assertEqual(result.returncode, 1)
        self.assertEqual(log.count("toolchain install"), 2)
        self.assertIn("Test Rust setup: failed to install Rust toolchain stable after 2 attempts", result.stderr)


if __name__ == "__main__":
    unittest.main()
