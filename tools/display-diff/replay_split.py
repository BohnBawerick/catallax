#!/usr/bin/env python3
"""Replay a green pre-split display-captures artifact through both layouts.

For example, download https://github.com/BohnBawerick/catallax/actions/runs/45898498941
and pass exceptions.jsonl from the recorded head and base checkouts.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from zipfile import ZipFile

from check import load


def verdict(root, flags):
    check = Path(__file__).with_name("check.py")
    result = subprocess.run(
        [sys.executable, "-B", check, "compare", root / "display-main",
         root / "display-head", "--gate", *flags],
        capture_output=True, text=True, check=False,
    )
    return result.returncode, result.stdout + result.stderr


def replay(recorded, split, flags, expected, diagnostic):
    # Each upload contains one file at the archive root. The two head
    # downloads go into the same directory without overwriting each other.
    for directory, filename in (
        ("display-main", "1.jsonl"),
        ("display-head", "1.jsonl"),
        ("display-head", "4.jsonl"),
    ):
        archive = split / f"{directory}-{filename}.zip"
        with ZipFile(archive, "w") as upload:
            upload.write(recorded / directory / filename, filename)
        with ZipFile(archive) as download:
            download.extractall(split / directory)

    original = verdict(recorded, flags)
    assembled = verdict(split, flags)
    if original != assembled or assembled[0] != expected or diagnostic not in assembled[1]:
        raise SystemExit(f"Expected exit {expected} with {diagnostic!r}\n"
                         f"Unsplit: {original}\nSplit: {assembled}")
    print(f"Both layouts exit {expected} ({diagnostic}): {assembled[1].splitlines()[-1]}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path)
    parser.add_argument("--exceptions", type=Path)
    parser.add_argument("--baseline-exceptions", type=Path)
    args = parser.parse_args()
    flags = []
    if args.exceptions:
        flags.extend(("--exceptions", args.exceptions))
    if args.baseline_exceptions:
        flags.extend(("--baseline-exceptions", args.baseline_exceptions))

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        recorded, split = root / "unsplit", root / "split"
        shutil.copytree(args.artifact, recorded)
        split.mkdir()
        replay(recorded, split, flags, 0, "0 unexpected changes")

        main = load(recorded / "display-main" / "1.jsonl")
        head = load(recorded / "display-head" / "1.jsonl")
        key = next(key for key, row in head.items()
                   if key in main and row["sha256"] == main[key]["sha256"] and not row["refused"])
        changed = head[key]["sha256"]
        head[key]["sha256"] = ("1" if changed[0] == "0" else "0") + changed[1:]
        rows = "".join(json.dumps(row) + "\n" for row in head.values())
        (recorded / "display-head" / "4.jsonl").write_text(rows)
        replay(recorded, split, flags, 1, "1 worker differences")
        (recorded / "display-head" / "1.jsonl").write_text(rows)
        replay(recorded, split, flags, 1, "1 unexpected changes, 0 new calls, 0 worker differences")


if __name__ == "__main__":
    main()
