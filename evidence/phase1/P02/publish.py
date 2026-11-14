"""Publish portable P02 copies from a private capture, without running geometry."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("capture", type=Path)
parser.add_argument("--checkout-root", required=True)
parser.add_argument("--binary-root", required=True)
args = parser.parse_args()
destination = Path(__file__).resolve().parent
capture = args.capture.resolve()
assert capture != destination and capture not in destination.parents
replacements = {
    args.checkout_root: "${P02_CHECKOUT}",
    args.binary_root: "${P02_BINARY_ROOT}",
}
sanitized = []
for source in sorted(capture.rglob("*")):
    if not source.is_file():
        continue
    relative = source.relative_to(capture)
    if relative.as_posix() in {"SHA256SUMS", "comparison.json", "publication.json", "publish.py"}:
        continue
    original = source.read_bytes()
    portable = original
    for location, token in replacements.items():
        portable = portable.replace(location.encode(), token.encode())
    target = destination / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(portable)
    if portable != original:
        sanitized.append(relative.as_posix())
(destination / "publication.json").write_text(json.dumps({
    "format": "p02-portable-publication/1",
    "label": "Sanitized portable copies, not original raw logs",
    "tokens": {
        "${P02_CHECKOUT}": "Checkout containing this evidence directory",
        "${P02_BINARY_ROOT}": "Checkout containing the observed pre-existing release binary",
    },
    "sanitized_files": sanitized,
}, indent=2) + "\n")
subprocess.run([sys.executable, "-B", str(destination / "analyze.py")], check=True)
(destination / "SHA256SUMS").write_text("".join(
    f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(destination).as_posix()}\n"
    for path in sorted(destination.rglob("*"))
    if path.is_file() and path.name != "SHA256SUMS"
))
