#!/usr/bin/env python3
"""Join cargo-dist's two macOS archives with Apple's lipo."""
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

from smoke import exercise


def main(directory):
    with tempfile.TemporaryDirectory(prefix="catallax-universal-") as temporary:
        root = Path(temporary)
        for target in ("aarch64", "x86_64"):
            with tarfile.open(directory / f"catallax-cli-{target}-apple-darwin.tar.xz") as archive:
                archive.extractall(root / target, filter="data")
        sources = {target: next((root / target).rglob("catallax-verify")).parent for target in ("aarch64", "x86_64")}
        output = root / "catallax-cli-universal2-apple-darwin"
        output.mkdir()
        for binary in ("catallax", "catallax-verify", "exact-brep"):
            subprocess.run(["lipo", "-create", str(sources["aarch64"] / binary), str(sources["x86_64"] / binary), "-output", str(output / binary)], check=True)
            subprocess.run(["lipo", "-verify_arch", "arm64", "x86_64", str(output / binary)], check=True)
            # Ad-hoc signing is required for arm64 execution, not an Apple Developer signature.
            subprocess.run(["codesign", "--force", "--sign", "-", str(output / binary)], check=True)
        for name in ("LICENSE", "README.md"):
            shutil.copy2(sources["aarch64"] / name, output / name)
        exercise(output, ("arch", "-arm64"))
        exercise(output, ("arch", "-x86_64"))
        destination = directory / "catallax-cli-universal2-apple-darwin.tar.xz"
        with tarfile.open(destination, "w:xz") as archive:
            archive.add(output, arcname=output.name)
        print(f"Universal archive includes arm64 and x86_64 executables: {destination.name}")


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
