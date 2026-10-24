#!/usr/bin/env python3
"""Exercise the actual programs extracted from a cargo-dist archive."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import zipfile


def exercise(directory, prefix=()):
    suffix = ".exe" if os.name == "nt" else ""
    cli = directory / ("catallax" + suffix)
    verify = directory / ("catallax-verify" + suffix)
    for binary in (cli, verify):
        subprocess.run([*prefix, str(binary), "about"], check=True)
    question = {"format": "exact-brep-question/2", "id": "release-smoke", "operands": {
        "a": {"shape": "block", "size": ["1", "2", "3"]}}, "steps": [], "result": "a", "limits": {"deadline_ms": 0}}
    with tempfile.TemporaryDirectory(prefix="catallax-release-smoke-") as temporary:
        out = Path(temporary)
        request = {"format": "catallax-request/1", "source": "question", "problem": question,
                   "out": str(out / "original"), "max_work": 1_000_000}
        result = subprocess.run([*prefix, str(verify), "stream"], input=json.dumps(request) + "\n", text=True, capture_output=True, check=True)
        response = json.loads(result.stdout)
        assert response["manifest"]["status"] == "answered", response
        if os.environ.get("GITHUB_SHA"):
            assert response["manifest"]["kernel"]["commit"] == os.environ["GITHUB_SHA"], "packaged source commit differs"
        step = (out / "original/answer.step").read_bytes()
        assert step.startswith(b"ISO-10303-21;"), "missing STEP exchange"
        subprocess.run([*prefix, str(verify), "replay", str(out / "original/question.json"), "--out", str(out / "replay")], check=True)
        assert (out / "replay/answer.step").read_bytes() == step, "portable replay changed STEP bytes"
    print(f"Packaged CLI and verifier built a block and replayed identical STEP bytes{': ' + ' '.join(prefix) if prefix else ''}")


def archive(path):
    with tempfile.TemporaryDirectory(prefix="catallax-release-unpack-") as temporary:
        root = Path(temporary)
        if path.suffix == ".zip":
            with zipfile.ZipFile(path) as bundle:
                bundle.extractall(root)
        else:
            with tarfile.open(path) as bundle:
                bundle.extractall(root, filter="data")
        suffix = ".exe" if os.name == "nt" else ""
        binaries = list(root.rglob("catallax-verify" + suffix))
        assert len(binaries) == 1, "archive must contain exactly one verifier"
        exercise(binaries[0].parent)


if __name__ == "__main__":
    archives = [Path(p) for p in sys.argv[1:] if p.endswith((".tar.xz", ".zip"))]
    assert len(archives) == 1, "expected one platform archive"
    archive(archives[0])
