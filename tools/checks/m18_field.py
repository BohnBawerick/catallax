#!/usr/bin/env python3
"""Capture native-field display and independent P1-P5 construction evidence."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys


def command(args, output):
    result = subprocess.run([str(a) for a in args], capture_output=True, check=False)
    output.write_bytes(result.stdout)
    output.with_suffix(output.suffix + ".stderr").write_bytes(result.stderr)
    if result.returncode:
        raise RuntimeError(f"{args[0]} exited {result.returncode}; see {output.name}")
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--display-only", action="store_true")
    parser.add_argument("--replay-inputs", type=Path)
    parser.add_argument("--baseline-binary", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    binary = args.binary.resolve()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    fixtures = out / "fixtures"
    fixtures.mkdir(exist_ok=True)
    source = root / "crates/catallax-features/tests/fixtures"
    for name in ["straight-sweep", "round-turn-sweep", "loft-natural-circles", "loft-spline-prism", "spline-sweep-planar", "spline-sweep-spatial"]:
        (fixtures / f"{name}.json").write_bytes((source / f"{name}.json").read_bytes())
    straight = json.loads((fixtures / "straight-sweep.json").read_text())
    profile = straight["steps"][0]
    translated = json.loads(json.dumps(profile))
    translated["id"] = "end"
    translated["plane"]["translation"] = ["9", "-2", "14"]
    loft = {"format": "catallax-recipe/4", "name": "translated loft", "steps": [profile, translated, {"id": "loft", "op": "loft", "sections": [{"profile": "profile", "start": ["a"]}, {"profile": "end", "start": ["a"]}], "shape": "ruled"}], "result": "loft"}
    (fixtures / "translated-loft.json").write_text(json.dumps(loft))
    circle = json.loads((fixtures / "loft-natural-circles.json").read_text())
    circle["steps"].pop(2)
    circle["steps"][-1]["sections"].pop()
    (fixtures / "two-circle-loft.json").write_text(json.dumps(circle))
    pose = {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]}
    smooth = {"format": "catallax-recipe/4", "name": "circular smooth blend", "steps": [{"id": "a", "op": "sphere", "radius": "2", "pose": pose}, {"id": "b", "op": "sphere", "radius": "1", "pose": pose}, {"id": "smooth", "op": "smooth_blend", "target": "a", "tool": "b", "kind": "union", "radius": "1"}], "result": "smooth"}
    (fixtures / "smooth.json").write_text(json.dumps(smooth))
    arc = {"format": "catallax-recipe/4", "name": "cardinal arc", "steps": [
        {"id": "profile", "op": "sketch", "plane": pose, "loops": [[{"kind": "circle", "id": "tube", "centre": ["4", "0"], "radius": "1"}]]},
        {"id": "path", "op": "path", "start": ["4", "0", "0"], "segments": [{"kind": "arc", "id": "turn", "centre": ["0", "0", "0"], "to": ["0", "0", "-4"], "axis": ["0", "1", "0"]}]},
        {"id": "sweep", "op": "sweep", "profile": "profile", "path": "path", "orientation": "rotation_minimising"}], "result": "sweep"}
    (fixtures / "cardinal-arc.json").write_text(json.dumps(arc))
    driver = root / ".agents/skills/verify-catallax/scripts/control-catallax.py"
    command([sys.executable, driver, "doctor", "--binary", binary, "--json"], out / "doctor.json")
    for workers in [1, 2, 4]:
        command([sys.executable, driver, "run", fixtures, "--source", "field", "--binary", binary, "--max-work", "20000000", "--workers", str(workers), "--fresh", "--out", out / f"field-{workers}", "--json"], out / f"field-{workers}.json")
    for workers in [2, 4]:
        result = json.loads(command([sys.executable, driver, "compare", out / "field-1/run.json", out / f"field-{workers}/run.json", "--json"], out / f"field-1-vs-{workers}.json"))
        if not result["equal"]:
            raise RuntimeError(f"field outcomes differ at {workers} workers")
    native = out / "native-steps"
    native.mkdir(exist_ok=True)
    if not args.display_only:
        for fixture in sorted(fixtures.glob("*.json")):
            command([binary, "steps", fixture, "--max-work", "20000000"], native / fixture.name)
    if args.display_only:
        return
    parts = ["P1-bracket", "P2-collar", "P3-filleted-boss", "P4-shelled-tray", "P5-flanged-housing"]
    previous = None
    summary = {}
    for workers in [1, 2, 4]:
        folder = out / f"steps-{workers}"
        folder.mkdir(exist_ok=True)

        def build(name):
            data = command([binary, "--threads", str(workers), "steps", root / f"parts/{name}.json", "--max-work", "20000000"], folder / f"{name}.json")
            return name, json.loads(data)

        with ThreadPoolExecutor(max_workers=workers) as pool:
            rows = dict(pool.map(build, parts))
        if previous is not None and rows != previous:
            raise RuntimeError(f"independent prefix outcomes/work differ at {workers} workers/threads")
        previous = rows
        summary[str(workers)] = {name: report["steps"][-1] for name, report in rows.items()}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if args.replay_inputs:
        if not args.baseline_binary:
            parser.error("--replay-inputs requires --baseline-binary")
        inputs = json.loads(args.replay_inputs.read_text())
        replay = out / "replay.jsonl"
        replay.write_text("".join(json.dumps(row["problem"]) + "\n" for row in inputs))
        for label, executable in [("baseline", args.baseline_binary.resolve()), ("candidate", binary)]:
            for workers in [1, 2, 4]:
                command([sys.executable, driver, "run", replay, "--binary", executable, "--workers", str(workers), "--fresh", "--out", out / f"{label}-{workers}", "--json"], out / f"{label}-{workers}.json")
        for label, workers in [("baseline", 2), ("baseline", 4), ("candidate", 1), ("candidate", 2), ("candidate", 4)]:
            result = json.loads(command([sys.executable, driver, "compare", out / "baseline-1/run.json", out / f"{label}-{workers}/run.json", "--json"], out / f"replay-vs-{label}-{workers}.json"))
            if not result["equal"]:
                raise RuntimeError(f"saved replay changed for {label} at {workers} workers")


if __name__ == "__main__":
    main()
