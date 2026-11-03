#!/usr/bin/env python3
"""Pair the five released SolidDraft edits and compare their exact native answers."""
import argparse
from collections import defaultdict
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import sys

TARGETS = {"bracket": {27: "edit_dimension", 43: "set_value"},
           "drawing": {24: "edit_extrude"},
           "references": {17: "edit_extrude", 18: "edit_dimension"}}
FILES = ("answer.step", "answer.txt", "display.obj", "mesh.txt",
         "names.txt", "operations.txt")
MEASURED_PAIRS = 9


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def pin(app, sha):
    source = app / "crates/soliddraft-kernel-catallax/src/lib.rs"
    text, count = re.subn(r'(pub const CATALLAX_REVISION: &str = ")[0-9a-f]{40}(";)',
                          rf'\g<1>{sha}\2', source.read_text())
    assert count == 1
    source.write_text(text)


def build(args):
    manifest = args.app / "crates/soliddraft-kernel-catallax/Cargo.toml"
    original = manifest.read_text()
    native = args.app / "crates/soliddraft-kernel-catallax/src/native.rs"
    source = native.read_text()
    old = "Ok(FaceRef { face: face_pattern, touching: vec![], near })"
    assert source.count(old) == 1, "unexpected released face-selector constructor"
    native.write_text(source.replace(
        old, "Ok(FaceRef { face: Some(face_pattern), touching: vec![], near })"))
    args.binaries.mkdir(parents=True)
    env = {**os.environ, "CARGO_TARGET_DIR": str(args.binaries.parent / "edit-build")}
    for revision, root, sha in [("baseline", args.baseline, os.environ["BASE_SHA"]),
                                ("candidate", args.candidate, os.environ["HEAD_SHA"])]:
        text, count = re.subn(
            r'^(catallax-[a-z-]+) = \{ git = "https://github.com/BohnBawerick/catallax", rev = "[0-9a-f]{40}" \}',
            lambda m: f'{m[1]} = {{ path = "{root / "crates" / m[1]}" }}', original, flags=re.M)
        assert count == 13, f"unexpected adapter dependency inventory: {count}"
        manifest.write_text(text)
        pin(args.app, sha)
        subprocess.run(["cargo", "build", "--release", "-p", "soliddraft-cli", "-j", "2"],
                       cwd=args.app, env=env, check=True)
        binary = Path(env["CARGO_TARGET_DIR"]) / "release/soliddraft"
        (args.binaries / revision).write_bytes(binary.read_bytes())
        (args.binaries / revision).chmod(0o755)


def without_clock(value):
    if isinstance(value, dict):
        return {k: without_clock(v) for k, v in value.items() if k != "wall_ms"}
    if isinstance(value, list):
        return [without_clock(v) for v in value]
    return value


def call_answer(call):
    value = without_clock(call)
    if value["work"] is not None:
        value["work"].pop("used")
    if value.get("operation_log") is not None:
        # Gate request work separately; dispatch counts and cumulative per-step
        # prefixes are diagnostics, not answers or certificates.
        for operation in value["operation_log"]["operations"]:
            operation.pop("work")
            operation.pop("decisions")
    return value


def accepted_work(value, where):
    if type(value) is not int or value < 0:
        raise ValueError(f"accepted work must be a nonnegative integer: {value!r}: {where}")
    return value


def read_native_work(path):
    return accepted_work(int(path.read_text()), path)


def compare_calls(left, right, where):
    def groups(calls):
        answers = defaultdict(list)
        for call in calls:
            cost = None if call["work"] is None else accepted_work(call["work"]["used"], where)
            answers[json.dumps(call_answer(call), sort_keys=True)].append(cost)
        return answers
    base, head = groups(left), groups(right)
    assert base.keys() == head.keys(), f"operation answers, marks or refusals changed: {where}"
    for answer, before in base.items():
        after = head[answer]
        assert len(before) == len(after), f"operation multiplicity changed: {where}"
        if before[0] is None:
            continue
        before, after = sorted(before), sorted(after)
        assert all(b <= a for a, b in zip(before, after)), \
            f"accepted call work rose from {before} to {after}: {where}"
    print(json.dumps(dict(comparison=str(where),
                          base_work=sum(c["work"]["used"] for c in left
                                        if c["work"] is not None and not c.get("nested", False)),
                          head_work=sum(c["work"]["used"] for c in right
                                        if c["work"] is not None and not c.get("nested", False)),
                          unmeasured_calls=sum(c["work"] is None for c in right),
                          matched_calls=len(right))), flush=True)


def compare_work(left, right):
    for name in ("kernel.json", "display.kernel.json", "view.kernel.json", "measurement.kernel.json"):
        a, b = left / name, right / name
        assert a.exists() == b.exists(), f"missing operation evidence: {a}, {b}"
        if a.exists():
            observed = [json.loads(path.read_text()) for path in (a, b)]
            foreground = [[call for call in value["calls"] if call["on_capture_thread"]]
                          for value in observed]
            for value, calls in zip(observed, foreground):
                value["calls"] = [call_answer(call) for call in calls]
            assert without_clock(observed[0]) == without_clock(observed[1]), \
                f"foreground questions, answers, marks or certificates changed: {a}, {b}"
            compare_calls(*foreground, b)


def compare_application(comparator, left, right, pins):
    for difference in comparator.problems(left, right, pins):
        if re.fullmatch(r"\d{4}/command.json", difference):
            values = [json.loads(comparator.text(root / "commands" / difference, root, pins))
                      for root in (left, right)]
            # Completion delivery may cross a capture boundary. Keep the raw
            # counter, and compare actual background work over the whole run.
            for value in values:
                value["rebuild"].pop("background_meshes_rebuilt")
            if values[0] == values[1]:
                continue
        raise AssertionError(f"application answer differences in {right}: {difference}")
    def calls(root):
        return [call for path in root.glob("commands/*/*kernel.json")
                for call in json.loads(path.read_text())["calls"]]
    compare_calls(calls(left), calls(right), right)


def run(args):
    args.out.mkdir(parents=True)
    core = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {core})
    comparator = load(args.app / "ci/compare-scenario-bundles.py", "compare_bundles")
    pins = os.environ["BASE_SHA"], os.environ["HEAD_SHA"]
    samples, references, native_references = [], {}, {}
    scenarios = args.out / "inputs"
    scenarios.mkdir()
    for name, targets in TARGETS.items():
        source = json.loads((args.app / f"ci/scenarios/{name}.json").read_text())
        for number, tool in targets.items():
            assert source["steps"][number - 1]["tool"] == tool, f"wrong edit target: {name} {number}"
        source["steps"] = source["steps"][:max(targets)]
        # File-export acceptance is outside the timed clicks and belongs to the
        # existing application suite. Keep every command, view and measurement.
        for step in source["steps"]:
            step.pop("checked_step_export", None)
            step.pop("drawing", None)
        (scenarios / f"{name}.json").write_text(json.dumps(source, indent=2) + "\n")
    for pair in range(MEASURED_PAIRS + 1):
        order = ["baseline", "candidate"] if pair % 2 == 0 else ["candidate", "baseline"]
        for revision in order:
            sha = pins[revision == "candidate"]
            root = args.baseline if revision == "baseline" else args.candidate
            pin(args.app, sha)
            for name, targets in TARGETS.items():
                dest = args.out / f"{pair}-{revision}" / name
                load_before = os.getloadavg()
                command = [sys.executable, str(args.app / "ci/run-scenarios.py"),
                           str(scenarios / f"{name}.json"), "--out", str(dest),
                           "--soliddraft", str(args.binaries / revision),
                           "--catallax-verify", str(root / "target/release/catallax-verify")]
                dest.parent.mkdir(exist_ok=True)
                with (dest.parent / f"{name}.log").open("w") as log:
                    result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
                assert result.returncode == 0, f"scenario failed: {dest.parent / (name + '.log')}"
                report = json.loads((dest / "result.json").read_text())
                assert report["pass"], f"incomplete or failed scenario: {dest}"
                if name in references:
                    compare_application(comparator, references[name], dest, pins)
                else:
                    references[name] = dest
                for number in targets:
                    folder = dest / "commands" / f"{number:04}"
                    compare_work(references[name] / "commands" / f"{number:04}", folder)
                    step = report["steps"][number - 1]
                    timing = step["timing"]
                    calls = []
                    for log_name in ("kernel.json", "display.kernel.json"):
                        calls.extend(json.loads((folder / log_name).read_text())["calls"])
                    foreground = [c for c in calls if c["on_capture_thread"] and not c.get("nested", False)]
                    native = json.loads((folder / "kernel.json").read_text())["questions"][-1]["problem"]
                    assert native["format"].startswith("catallax-recipe/")
                    proof = dest.parent / "native" / name / f"{number:04}"
                    proof.mkdir(parents=True)
                    native_input = proof / "native.json"
                    native_input.write_text(json.dumps(native, indent=2) + "\n")
                    captured = subprocess.run([str(root / "target/release/examples/edit_answer"),
                                               str(native_input), str(proof)],
                                              text=True, capture_output=True, check=True)
                    (proof / "timing.json").write_text(captured.stdout)
                    key = name, number
                    answer = {f: hashlib.sha256((proof / f).read_bytes()).hexdigest() for f in FILES}
                    work = read_native_work(proof / "work.txt")
                    if key not in native_references:
                        native_references[key] = dict(answer=answer, work=work)
                    reference = native_references[key]
                    assert answer == reference["answer"], f"exact answer/certificate changed: {proof}"
                    assert work <= reference["work"], \
                        f"accepted native work rose from {reference['work']} to {work}: {proof}"
                    row = dict(pair=pair, warmup=pair == 0, revision=revision, scenario=name, step=number,
                               tool=step["tool"], click_ms=timing["elapsed_ms"], split_ms=timing["split_ms"],
                               kernel_ms=sum(c["wall_ms"] for c in foreground),
                               work=timing["work"], certified=json.loads(captured.stdout),
                               base_work=reference["work"],
                               head_work=work if revision == "candidate" else None,
                               accepted_work=work,
                               answer_sha256=answer, load_before=load_before, load_after=os.getloadavg())
                    samples.append(row)
                    print(json.dumps(row), flush=True)
                (args.out / "samples.json").write_text(json.dumps(samples, indent=2) + "\n")
    medians = []
    for name, targets in TARGETS.items():
        for number in targets:
            for revision in ("baseline", "candidate"):
                rows = [s for s in samples if s["scenario"] == name and s["step"] == number
                        and s["revision"] == revision and not s["warmup"]]
                medians.append(dict(scenario=name, step=number, revision=revision,
                                    click_ms=statistics.median(s["click_ms"] for s in rows),
                                    kernel_ms=statistics.median(s["kernel_ms"] for s in rows),
                                    application_ms=statistics.median(s["click_ms"] - s["kernel_ms"] for s in rows),
                                    split_ms={phase: statistics.median(s["split_ms"][phase] for s in rows)
                                              for phase in ("edit", "display", "measurement")},
                                    certified_build_ms=statistics.median(s["certified"]["timing"]["build_ms"] for s in rows),
                                    certified_display_ms=statistics.median(s["certified"]["timing"]["display_ms"] for s in rows)))
    report = dict(base=pins[0], head=pins[1], soliddraft="04f45bf40c119d5387836173727a3905db094e4a",
                  core=core, measured_pairs=MEASURED_PAIRS,
                  binary_sha256={r: hashlib.sha256((args.binaries / r).read_bytes()).hexdigest()
                                            for r in ("baseline", "candidate")},
                  comparison="identical application answers and native STEP, mesh, certificates, names and operation outcomes; accepted work never increases",
                  medians=medians, samples=samples)
    (args.out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


def profile(args):
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    binary = args.baseline / "target/release/examples/edit_answer"
    metadata = dict(binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    valgrind=subprocess.check_output(["valgrind", "--version"], text=True).strip(),
                    metric="Callgrind instruction costs, not wall-time measurements")
    for name, targets in TARGETS.items():
        for number in targets:
            source = args.out / "0-baseline/native" / name / f"{number:04}/native.json"
            dest = args.out / "profiles" / name / f"{number:04}"
            dest.mkdir(parents=True)
            with (dest / "instrumented.stdout").open("w") as stdout, \
                    (dest / "callgrind.log").open("w") as stderr:
                subprocess.run(["valgrind", "--tool=callgrind", "--error-exitcode=1",
                                f"--callgrind-out-file={dest / 'callgrind.out'}",
                                str(binary), str(source), str(dest / "answer")],
                               stdout=stdout, stderr=stderr, check=True)
            with (dest / "inclusive.txt").open("w") as stdout:
                subprocess.run(["callgrind_annotate", "--inclusive=yes", "--tree=both",
                                "--auto=no", "--threshold=99", str(dest / "callgrind.out")],
                               stdout=stdout, check=True)
            metadata.setdefault("inputs", []).append(dict(
                scenario=name, step=number, sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    (args.out / "profiles/metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("build", "run", "profile"))
    for name in ("app", "baseline", "candidate", "binaries", "out"):
        parser.add_argument("--" + name, type=lambda value: Path(value).resolve(), required=True)
    arguments = parser.parse_args()
    {"build": build, "run": run, "profile": profile}[arguments.mode](arguments)
