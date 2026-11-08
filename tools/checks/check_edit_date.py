#!/usr/bin/env python3
"""Exercise the frozen title-date producer across midnight, keeping comparison strict."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess

import edit_speed


CLOCK = """
mod web_time {
    pub const UNIX_EPOCH: std::time::SystemTime = std::time::UNIX_EPOCH;
    pub struct SystemTime;
    impl SystemTime {
        pub fn now() -> std::time::SystemTime {
            let secs = std::env::var("EDIT_CHECK_WALL_EPOCH").unwrap().parse::<u64>().unwrap();
            UNIX_EPOCH + std::time::Duration::from_secs(secs)
        }
    }
}
"""
MAIN = """
fn main() {
    let wall = web_time::SystemTime::now().duration_since(web_time::UNIX_EPOCH).unwrap().as_secs();
    println!("{}|{}", wall, today());
}
"""


def bundle(path, date):
    command = path / "commands/0001"
    command.mkdir(parents=True)
    (path / "scenario.json").write_text(json.dumps({"steps": [{"tool": "make_drawing"}]}))
    (path / "result.json").write_text(json.dumps({"steps": [{"pass": True}]}))
    (path / "protocol.json").write_text(json.dumps({"drawing": {"title": {"text": "Date control", "date": date}}}))
    (command / "command.json").write_text(json.dumps({
        "summary": {"drawing": {"date": date}},
        "rebuild": {"meshes_rebuilt": 0, "background_meshes_rebuilt": 0},
    }))


def check(app, out):
    out.mkdir(parents=True)
    epoch = int(os.environ["SOURCE_DATE_EPOCH"])
    expected = datetime.fromtimestamp(epoch, timezone.utc).date().isoformat()
    source = (app / "crates/soliddraft-app/src/drawing/mod.rs").read_text()
    functions = re.findall(r"^pub fn today\(\) -> String \{\n.*?^\}", source, re.M | re.S)
    assert len(functions) == 1, "unexpected released title-date function"
    # Compile the actual released producer after setup, not a copy of its date
    # algorithm. The shim controls only its wall-clock dependency.
    probe = out / "date-probe.rs"
    probe.write_text(CLOCK + functions[0] + MAIN)
    binary = out / "date-probe"
    compiled = subprocess.run(["rustc", "--edition", "2024", "--crate-name", "released_edit_date_probe",
                               str(probe), "-o", str(binary)], text=True, capture_output=True)
    (out / "rustc.log").write_text(compiled.stdout + compiled.stderr)
    compiled.check_returncode()

    midnight = (epoch // 86_400 + 1) * 86_400
    walls = [midnight - 5, midnight + 40]
    dates = []
    for wall in walls:
        observed = subprocess.run([str(binary)], env={**os.environ, "EDIT_CHECK_WALL_EPOCH": str(wall)},
                                  check=True, text=True, capture_output=True)
        seen_wall, date = observed.stdout.strip().split("|", 1)
        assert int(seen_wall) == wall, "probe did not cross the requested wall-clock boundary"
        assert date == expected, f"title date changed with the wall clock: {date} != {expected}"
        dates.append(date)

    left, right = out / "before", out / "after"
    bundle(left, dates[0])
    bundle(right, dates[1])
    comparator = edit_speed.load(app / "ci/compare-scenario-bundles.py", "released_date_comparator")
    pins = ("a" * 40, "b" * 40)
    edit_speed.compare_application(comparator, left, right, pins)

    # A genuine date difference must still fail the existing strict comparison.
    protocol = json.loads((right / "protocol.json").read_text())
    protocol["drawing"]["title"]["date"] = datetime.fromtimestamp(midnight, timezone.utc).date().isoformat()
    (right / "protocol.json").write_text(json.dumps(protocol))
    try:
        edit_speed.compare_application(comparator, left, right, pins)
    except AssertionError:
        pass
    else:
        raise AssertionError("strict comparison accepted a changed title date")
    report = {"epoch": epoch, "wall_epochs": walls, "dates": dates,
              "changed_date_rejected": True}
    (out / "date-check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--app", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    check(args.app.resolve(), args.out.resolve())
