#!/usr/bin/env python3
"""Join cargo-dist's two macOS archives with Apple's lipo."""
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

from smoke import exercise


def _copy_static_assets(source, destination):
    for name in ("LICENSE", "NOTICE", "README.md"):
        shutil.copy2(source / name, destination / name)
    shutil.copytree(source / "licenses", destination / "licenses")


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
            subprocess.run(["lipo", str(output / binary), "-verify_arch", "arm64", "x86_64"], check=True)
            # Ad-hoc signing is required for arm64 execution, not an Apple Developer signature.
            subprocess.run(["codesign", "--force", "--sign", "-", str(output / binary)], check=True)
        _copy_static_assets(sources["aarch64"], output)
        exercise(output, ("arch", "-arm64"))
        exercise(output, ("arch", "-x86_64"))
        destination = directory / "catallax-cli-universal2-apple-darwin.tar.xz"
        with tarfile.open(destination, "w:xz") as archive:
            archive.add(output, arcname=output.name)
        print(f"Universal archive includes arm64 and x86_64 executables: {destination.name}")


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
