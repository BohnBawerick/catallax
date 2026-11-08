#!/usr/bin/env python3
"""Pair the five released SolidDraft edits and compare their exact native answers."""
import argparse
from collections import Counter, defaultdict
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
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


def pin_date(app, epoch):
    source = app / "crates/soliddraft-app/src/drawing/mod.rs"
    text = source.read_text()
    old = "let secs = web_time::SystemTime::now().duration_since(web_time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0);"
    assert text.count(old) == 1, "unexpected released title-date clock"
    source.write_text(text.replace(old, f"let secs = {epoch}u64;"))


def build(args):
    manifest = args.app / "crates/soliddraft-kernel-catallax/Cargo.toml"
    original = manifest.read_text()
    pin_date(args.app, int(os.environ["SOURCE_DATE_EPOCH"]))
    native = args.app / "crates/soliddraft-kernel-catallax/src/native.rs"
    source = native.read_text()
    old = "Ok(FaceRef { face: face_pattern, touching: vec![], near })"
    assert source.count(old) == 1, "unexpected released face-selector constructor"
    native.write_text(source.replace(
        old, "Ok(FaceRef { face: Some(face_pattern), touching: vec![], near })"))
    adapter = args.app / "crates/soliddraft-kernel-catallax/src/lib.rs"
    text = adapter.read_text()
    old = '                radius: native::number(radius, "fillet radius")?,'
    assert text.count(old) == 1, "unexpected released fillet constructor"
    adapter.write_text(text.replace(
        old, old + "\n                overflow: None,\n                corner: None,"))
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


def read_native_work(path, reference_work=None):
    work = accepted_work(int(path.read_text()), path)
    if reference_work is not None:
        assert work <= reference_work, \
            f"accepted work rose from {reference_work} to {work}: {path}"
    return work


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


def phases(args):
    root = args.out / "boolean-phases"
    root.mkdir()
    compiled = subprocess.run(
        ["cargo", "build", "--locked", "--release", "-p", "catallax-features",
         "--bench", "rebuild", "--message-format=json", "-j", "2"],
        cwd=args.baseline, stdout=subprocess.PIPE, text=True, check=True)
    (root / "build.jsonl").write_text(compiled.stdout)
    artifacts = [json.loads(line) for line in compiled.stdout.splitlines()]
    executables = [item["executable"] for item in artifacts
                   if item.get("reason") == "compiler-artifact"
                   and item["target"]["name"] == "rebuild" and item.get("executable")]
    assert len(executables) == 1, "expected one rebuild benchmark executable"
    binary = Path(executables[0])
    core = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {core})
    stages = ("intersect", "imprint", "classify", "keep", "sew")
    samples, medians = [], []
    for name, targets in TARGETS.items():
        for number in targets:
            source = args.out / "0-baseline/native" / name / f"{number:04}/native.json"
            expected = {step["id"] for step in json.loads(source.read_text())["steps"]
                        if step["op"] == "boolean"}
            reference = None
            case = []
            for repetition in range(4):
                dest = root / name / f"{number:04}" / str(repetition)
                dest.mkdir(parents=True)
                load_before = os.getloadavg()
                process = subprocess.run(
                    [str(binary), "--input", str(source), "--repeat", "1", "--stages",
                     "--step-output", str(dest / "answer.step")],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                load_after = os.getloadavg()
                (dest / "stdout.jsonl").write_text(process.stdout)
                (dest / "stderr.txt").write_text(process.stderr)
                process.check_returncode()
                records = [json.loads(line) for line in process.stdout.splitlines()]
                rows = [row for row in records if "stage" in row]
                assert Counter((row["step"], row["stage"]) for row in rows) == Counter(
                    (step, stage) for step in expected for stage in stages), \
                    f"missing Boolean phase evidence: {dest}"
                normalized = without_clock(rows)
                if reference is None:
                    reference = normalized
                assert normalized == reference, f"phase decisions or work varied: {dest}"
                assert (dest / "answer.step").read_bytes() == \
                    (source.parent / "answer.step").read_bytes(), \
                    f"phase benchmark's complete build changed STEP bytes: {dest}"
                row = dict(scenario=name, step=number, repetition=repetition,
                           warmup=repetition == 0, load_before=load_before, load_after=load_after,
                           stages_ms={stage: sum(value["wall_ms"] for value in rows
                                                 if value["stage"] == stage) for stage in stages},
                           operations=rows)
                samples.append(row)
                if repetition:
                    case.append(row)
            medians.append(dict(
                scenario=name, step=number,
                stages_ms={stage: statistics.median(row["stages_ms"][stage] for row in case)
                           for stage in stages},
                operations=[dict(
                    step=step,
                    stages_ms={stage: statistics.median(
                        value["wall_ms"] for row in case for value in row["operations"]
                        if value["step"] == step and value["stage"] == stage) for stage in stages})
                    for step in sorted(expected)]))
    report = dict(
        base=os.environ["BASE_SHA"], core=core,
        binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        measured_repetitions=3,
        scope="Independent native prefix phase replay with the existing benchmark's raw scope; "
              "not an application click or its aggregate work scope",
        medians=medians, samples=samples)
    (root / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


def check_profile_answer(source, dest):
    for name in (*FILES, "stats.txt", "work.txt"):
        assert (dest / "answer" / name).read_bytes() == (source.parent / name).read_bytes(), \
            f"profiling changed an exact answer, certificate, statistic or work: {dest} {name}"


def capture_cpu_profile(binary, source, dest):
    # Use the installed tool directly: Ubuntu's /usr/bin/perf wrapper requires
    # a package matching the hosted VM's exact kernel, unlike the perf ABI.
    perf = next(path for path in sorted(Path("/usr/lib/linux-tools").glob("*/perf"))
                if path.is_file())
    dest.mkdir(parents=True)
    shutil.copyfile(binary, dest / "edit_answer")
    metadata = dict(tool=str(perf), version=subprocess.check_output(
        [str(perf), "--version"], text=True).strip(), frequency=997, stack_bytes=16384,
        metric="Sampled user CPU-clock periods; not instrumented elapsed time",
        core=min(os.sched_getaffinity(0)), load_average=os.getloadavg())
    metadata["inline_expansion"] = False
    # Raw stacks contain process memory: never inherit runner credentials.
    profile_env = dict(PATH=os.defpath, LANG="C.UTF-8", HOME=str(dest),
                       DEBUGINFOD_URLS="")
    for repetition in range(5):
        capture = dest / str(repetition)
        capture.mkdir()
        print(f"CPU capture {repetition}: recording native probe", flush=True)
        with (capture / "probe.stdout").open("w") as stdout, \
                (capture / "record.log").open("w") as stderr:
            subprocess.run([str(perf), "record", "--no-buildid-cache", "-e", "cpu-clock:u",
                            "-F", "997", "--call-graph", "dwarf,16384",
                            "-o", str(capture / "perf.data"), "--",
                            str(binary), str(source), str(capture / "answer")],
                           stdout=stdout, stderr=stderr, check=True, env=profile_env)
        print(f"CPU capture {repetition}: resolving DWARF call chains", flush=True)
        with (capture / "stacks.txt").open("w") as stdout, \
                (capture / "decode.log").open("w") as stderr:
            subprocess.run([str(perf), "script", "--demangle", "--no-inline",
                            "--ns", "-F", "+period", "-i", str(capture / "perf.data")],
                           stdout=stdout, stderr=stderr, check=True, env=profile_env)
        assert "catallax_boolean" in (capture / "stacks.txt").read_text(), \
            "CPU samples did not resolve Boolean call stacks"
        check_profile_answer(source, capture)
    (dest / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def capture_profile(binary, source, dest, options=()):
    dest.mkdir(parents=True)
    with (dest / "instrumented.stdout").open("w") as stdout, \
            (dest / "callgrind.log").open("w") as stderr:
        subprocess.run(["valgrind", "--tool=callgrind", "--error-exitcode=1", *options,
                        f"--callgrind-out-file={dest / 'callgrind.out'}",
                        str(binary), str(source), str(dest / "answer")],
                       stdout=stdout, stderr=stderr, check=True)
    with (dest / "inclusive.txt").open("w") as stdout:
        subprocess.run(["callgrind_annotate", "--inclusive=yes", "--tree=both",
                        "--auto=no", "--threshold=99", str(dest / "callgrind.out")],
                       stdout=stdout, check=True)
    check_profile_answer(source, dest)


def profile(args):
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    binary = args.baseline / "target/release/examples/edit_answer"
    metadata = dict(binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    valgrind=subprocess.check_output(["valgrind", "--version"], text=True).strip(),
                    metric="Callgrind instruction costs, not wall-time measurements")
    for name, targets in TARGETS.items():
        for number in targets:
            source = args.out / "0-baseline/native" / name / f"{number:04}/native.json"
            capture_profile(binary, source, args.out / "profiles" / name / f"{number:04}")
            if name == "bracket" and number == 27:
                functions = ("catallax_boolean::graph::section",
                             "catallax_boolean::graph::carrier_hits",
                             "catallax_boolean::imprinting::imprint")
                capture_profile(
                    binary, source, args.out / "profiles/non-carrier/bracket/0027",
                    ["--collect-atstart=no",
                     *(f"--toggle-collect={function}" for function in functions)])
                metadata["non_carrier"] = dict(
                    included=[functions[0], functions[2]], excluded=[functions[1]],
                    scope="Rest of intersect plus imprint; carrier collection toggles off/on")
                capture_cpu_profile(binary, source, args.out / "profiles/cpu/bracket/0027")
            metadata.setdefault("inputs", []).append(dict(
                scenario=name, step=number, sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    (args.out / "profiles/metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("build", "run", "phases", "profile"))
    for name in ("app", "baseline", "candidate", "binaries", "out"):
        parser.add_argument("--" + name, type=lambda value: Path(value).resolve(), required=True)
    arguments = parser.parse_args()
    {"build": build, "run": run, "phases": phases, "profile": profile}[arguments.mode](arguments)
