#!/usr/bin/env python3
"""Compare public display results on two revisions.

Run capture inside each disposable checkout with the same copy of this tool and
its fixtures. Capture writes one JSONL file per worker count.

Two collectors produce the rows. The interface collector, the default, builds
the tests with the `display-capture` feature: mesh and field display record
an explicit, versioned encoding of each result and its work. The probe
collector rewrites four function bodies and hashes `Debug` text; keep it for
revisions that predate the interface. Compare captures from one collector
only, and pin exceptions against that collector's rows. `migrate` checks that
both collectors on one revision record the same calls, results and work.

The CI gate is `compare --gate`: it compares interface digests with main
and requires the candidate's 1- and 4-worker captures to match exactly.
"""
import argparse
import collections
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

PROBES = {
    "catallax-mesh/src/lib.rs": ("mesh_body_impl",),
    "catallax-mesh/src/display.rs": ("display_mesh_body",),
    "catallax-mesh/src/preview.rs": ("preview_mesh_body",),
    "catallax-field/src/display.rs": ("mesh_plan_inner",),
}
THREAD_COUNTS = (1, 2, 4, 24)
GATE_THREADS = (1, 4)
# Cargo features that turn on the interface collector in each tested package.
CAPTURE_FEATURES = {
    "catallax-mesh": "display-capture",
    "catallax-field": "display-capture",
    "catallax-app": "catallax-field/display-capture",
    "catallax-features": "catallax-field/display-capture",
}
# Formats 1 and 2 share a digest encoding; only 2 exposes the outcome.
INTERFACE_FORMAT = 2
SUITES = {
    "catallax-mesh": "native_boundaries preview_boundary",
    "catallax-app": (
        "coaxial_fillet_display coil_display ellipse_profiles fillet_display fillet_over_fillet_display "
        "application fine_mesh generator_fillet_display handover open_spherical_cap "
        "p2_display p5_display preview_mesh rebuild_work "
        "rational_lofts rational_patches rolling_fillets spherical_corner_mesh "
        "spline_display spline_sweeps straight_sweep_display"
    ),
    "catallax-features": (
        "algebraic_extrusions algebraic_revolutions bezier_profiles certified_prism "
        "coaxial_fillets display_diff_ellipse draft_annular draft_annular_clearance "
        "draft_oblique draft_preserved_bore "
        "draft_preserved_boss equal_crossing_fillets fillet_over_fillet "
        "fillet_sphere_ring g2_fillets m11_spline_fillet oblique_cone_fillets "
        "oblique_fillets ruled_edge_rows saddle_fillets slot_fillets strip_ends "
        "variable_fillets bored_shell boolean_face_thicken circular_sphere_end_conditions "
        "concave_draft cone_distance disjoint_bezier drawing_basic_shapes nonconvex_shell "
        "nonconvex_thicken polygon_end_conditions torus_sections trihedral_corner"
    ),
    "catallax-field": (
        "algebraic_extrusions chamfer_display curved_end_conditions display display_preview "
        "lattice lattice_stage m18_native nearest_face parallel patch_shell "
        "rolling_conics rolling_ends rolling_over rolling_quadrics "
        "rolling_tori spline_spine unproved_preview"
    ),
}
SHARED_CASES = {
    "catallax-app/tests/spline_sweeps.rs": "spline_sweeps.rs",
    "catallax-features/tests/certified_prism.rs": "certified_prism.rs",
    "catallax-features/tests/equal_crossing_fillets.rs": "equal_crossing_fillets.rs",
    "catallax-features/tests/variable_fillets.rs": "variable_fillets.rs",
    "catallax-field/tests/patch_shell.rs": "patch_shell.rs",
}

DISPLAY_PATHS = (
    ".github/actions/rust-toolchain/",
    ".github/workflows/display-diff.yml",
    "tools/display-diff/",
    "crates/catallax-mesh/",
    "crates/catallax-app/Cargo.toml",
    "crates/catallax-app/src/certified.rs",
    "crates/catallax-app/src/lib.rs",
    "crates/catallax-app/src/lifecycle.rs",
    "crates/catallax-app/src/numerical_sweeps.rs",
    "crates/catallax-app/src/preview.rs",
    "crates/catallax-app/src/views.rs",
    "crates/catallax-field/Cargo.toml",
    "crates/catallax-field/src/display.rs",
    "crates/catallax-field/src/display/",
    "crates/catallax-field/src/mesh.rs",
    "crates/catallax-field/src/patch_shell.rs",
    "crates/catallax-field/src/patch_shell/",
    "crates/catallax-field/src/unproved_batch.rs",
    "crates/catallax-field/src/unproved_wgsl.rs",
    "crates/catallax-field/src/workers.rs",
)
DISPLAY_TESTS = frozenset(
    f"crates/{package}/tests/{name}.rs"
    for package, suites in SUITES.items() for name in suites.split()
)


def affects_display(paths):
    return any(
        path in DISPLAY_TESTS or any(
            path.startswith(rule) if rule.endswith("/") else path == rule
            for rule in DISPLAY_PATHS
        )
        for path in paths
    )

RUST_PROBE = r'''
fn differential_record(
    function: &str,
    before: Option<catallax_numbers::OperationWork>,
    tags: &str,
    options: &impl std::fmt::Debug,
    result: &impl std::fmt::Debug,
) {
    use std::io::Write;
    use std::process::{Command, Stdio};
    let work = before.map_or_else(|| "null".to_owned(), |before| {
        let after = catallax_numbers::current_operation_work();
        let counters = [
            after.used.wrapping_sub(before.used),
            after.refunded.wrapping_sub(before.refunded),
            after.uncharged.wrapping_sub(before.uncharged),
            u64::from(after.exhausted),
        ];
        format!("{counters:?}")
    });
    let mut child = Command::new("sha256sum")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("start sha256sum");
    {
        let mut input = std::io::BufWriter::new(child.stdin.take().unwrap());
        write!(input, "{result:?}").unwrap();
        input.flush().unwrap();
    }
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    let digest = std::str::from_utf8(&output.stdout)
        .unwrap()
        .split_whitespace()
        .next()
        .unwrap();
    let thread = std::thread::current();
    eprintln!(
        "\nmeshdiff\t{}\t{function}\t{tags}\t{options:?}\t{work}\t{digest}",
        thread.name().unwrap_or("unnamed"),
    );
}
'''

CARRIERS = r'''
    let tags = model.faces_of_body(body).into_iter().filter_map(|face| {
        match model.geom().surface(model.face(face).surface) {
            catallax_shapes::Surface::BSpline(support) => Some(match support {
                catallax_shapes::SplineSupport::BezierExtrusion(_) => "bezier",
                catallax_shapes::SplineSupport::EllipticOffset(_) => "ellipse-offset",
                catallax_shapes::SplineSupport::LinearFillet(_) => "linear-fillet",
                catallax_shapes::SplineSupport::ObliqueDraft(_) => "oblique-draft",
                catallax_shapes::SplineSupport::RationalAtlas(_) => "atlas",
                catallax_shapes::SplineSupport::RationalPatch(_) => "rational-patch",
                catallax_shapes::SplineSupport::RollingEnvelope(_) => "envelope",
                catallax_shapes::SplineSupport::RollingMitre(_) => "mitre",
                catallax_shapes::SplineSupport::RollingPipe(_) => "pipe",
                catallax_shapes::SplineSupport::RollingSaddle(_) => "saddle",
                catallax_shapes::SplineSupport::SplineFillet(_) => "spline-fillet",
                _ => "other-spline",
            }),
            _ => None,
        }
    }).collect::<std::collections::BTreeSet<_>>().into_iter().collect::<Vec<_>>().join(",");
'''


def instrument(source, names):
    for name in names:
        # Revisions with the capture interface run display_mesh_body's work
        # inside display_within_allowance. Record that function whole.
        code_name = name
        if name == "display_mesh_body" and "fn display_within_allowance(" in source:
            code_name = "display_within_allowance"
        # Anchor at the function's own declaration. Its closing brace is unindented.
        pattern = rf"(^[a-z() ]*fn {code_name}\([^{{]*\{{\n)(.*?)(^\}})"
        match = re.search(pattern, source, re.MULTILINE | re.DOTALL)
        if match is None:
            raise ValueError(f"no declaration for {name}")
        tags = '    let tags = String::new();\n' if name == "mesh_plan_inner" else CARRIERS
        first, last = match.span(2)
        old = match.group(2)
        # display_mesh_body establishes its own allowance on standalone calls.
        # Record the inner invocation within that allowance exactly once.
        prefix = ""
        if code_name == "display_mesh_body":
            boundary = "    let faces = model.faces_of_body(body);"
            start = old.index(boundary)
            prefix, old = old[:start], old[start:]
        wrapper = (
            prefix
            + tags
            + "    let before = catallax_numbers::operation_budget_active()\n"
            + "        .then(catallax_numbers::current_operation_work);\n"
            + "    let result = (|| {\n"
            + old
            + "    })();\n"
            + f'    differential_record("{name}", before, &tags, options, &result);\n'
            + "    result\n"
        )
        source = source[:first] + wrapper + source[last:]
    return source + RUST_PROBE


ELLIPSE_TEST = r"""
use catallax_features::evaluate;
use catallax_mesh::{MeshOptions, PreviewOptions, display_mesh_body, preview_mesh_body};
use catallax_numbers::{Decider, OperationResult, Rational, with_operation_budget};
use catallax_recipe::Recipe;

#[test]
fn ellipse_rim_at_both_chords() {
    let recipe = Recipe::parse(include_str!(".display_diff_ellipse.json")).unwrap();
    let mut dec = Decider::default();
    let built = evaluate(&recipe, &mut dec).unwrap();
    for chord in ["1/100", "1/200"] {
        let chord = Rational::parse(chord).unwrap();
        let options = MeshOptions { chord_tolerance: chord.clone(), ..MeshOptions::default() };
        let (result, work) = with_operation_budget(100_000_000, || {
            display_mesh_body(&built.model, built.result, &options, "ellipse-rim", &mut dec)
        });
        assert!(!work.exhausted);
        assert!(
            matches!(result, OperationResult::Proven(_) | OperationResult::Bounded(_)),
            "{result:?}"
        );
        let options = PreviewOptions { chord_tolerance: chord, ..PreviewOptions::default() };
        let (result, work) = with_operation_budget(100_000_000, || {
            preview_mesh_body(&built.model, built.result, &options, &mut dec)
        });
        assert!(!work.exhausted);
        result.unwrap();
    }
}
"""


@contextmanager
def probes(root, collector):
    originals = {}
    try:
        for filename, names in PROBES.items() if collector == "probe" else ():
            path = root / "crates" / filename
            originals[path] = path.read_text()
            path.write_text(instrument(originals[path], names))
        for filename, fixture in SHARED_CASES.items():
            path = root / "crates" / filename
            originals[path] = path.read_text()
            cases = Path(__file__).with_name("cases") / fixture
            path.write_text(originals[path] + "\n" + cases.read_text())
        yield
    finally:
        for path, source in originals.items():
            path.write_text(source)


def capture(root, args, output, threads):
    with output.open("w") as out:
        for package, names in SUITES.items():
            selected = [
                name for name in names.split()
                if not args.suite or name in args.suite
            ]
            missing = [name for name in selected if not (root / "crates" / package / "tests" / f"{name}.rs").is_file()]
            if missing:
                raise ValueError(f"missing {package} display test targets: {missing}")
            if not selected:
                continue
            cmd = ["cargo", "test", "--no-fail-fast", "--locked", "--release", "-p", package]
            if args.collector == "interface":
                cmd += ["--features", CAPTURE_FEATURES[package]]
            for name in selected:
                cmd += ["--test", name]
            cmd += ["--", "--test-threads=1", "--nocapture"]
            counts = collections.Counter()
            suite = None
            print(f"{package}: {' '.join(selected)}", file=sys.stderr, flush=True)
            with tempfile.TemporaryFile(mode="w+t") as log:
                with subprocess.Popen(
                    cmd, cwd=root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, env={**os.environ, "CATALLAX_THREADS": str(threads)},
                ) as run:
                    for line in run.stdout:
                        log.write(line)
                        target = re.match(r"\s*Running tests/(\w+)\.rs ", line)
                        if target:
                            suite = target[1]
                        row = parse(line, args.collector)
                        if row is None:
                            continue
                        assert suite in selected, (suite, line)
                        key = (package, suite, *row.pop("call"))
                        counts[key] += 1
                        row["key"] = [*key, counts[key]]
                        print(json.dumps(row, sort_keys=True), file=out, flush=True)
                    code = run.wait()
                if code:
                    log.seek(0)
                    for line in log:
                        if parse(line, args.collector) is None:
                            sys.stderr.write(line)
                print(json.dumps({"package": package, "exit": code, "calls": sum(counts.values())}),
                      file=sys.stderr, flush=True)
                if code:
                    raise SystemExit(code)


def parse(line, collector):
    """One record line as a row with its call key, or None for other output."""
    fields = line.rstrip("\n").split("\t")
    if collector == "probe" and fields[0] == "meshdiff":
        case, function, tags, options, work, digest = fields[1:]
        return {"call": [case, function, options], "carriers": tags,
                "work": json.loads(work), "sha256": digest}
    if collector == "interface" and fields[0] == "displaycapture":
        version, case, function, tags, options, work, digest, *status = fields[1:]
        version = int(version)
        if version not in (1, INTERFACE_FORMAT):
            raise ValueError(f"unsupported display record format {version}")
        row = {"call": [case, function, options], "carriers": tags, "format": version,
               "work": json.loads(work), "sha256": digest}
        if version == INTERFACE_FORMAT and status in (["ok"], ["refused"]):
            row["refused"] = status == ["refused"]
        elif version != 1 or status:
            raise ValueError(f"invalid outcome for display record format {version}: {status}")
        return row
    return None


def load(path):
    result = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        row = json.loads(line)
        if "work" not in row:
            raise ValueError(f"{path}:{number}: missing work field, use null for unmeasured work")
        key = tuple(row["key"])
        if key in result:
            raise ValueError(f"{path}:{number}: duplicate case: {key}")
        result[key] = row
    if not result:
        raise ValueError(f"empty capture: {path}")
    return result


def result_hash(row):
    # Pin both the result digest and its accounting and carrier metadata.
    return hashlib.sha256(json.dumps(row, sort_keys=True).encode()).hexdigest()


def load_exceptions(path):
    exceptions = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        row = json.loads(line)
        fields = {"key", "reason", "main_sha256", "candidate_sha256"}
        if (
            not isinstance(row, dict) or row.keys() != fields
            or not isinstance(row["key"], list) or len(row["key"]) != 6
            or any(not isinstance(value, str) for value in row["key"][:-1])
            or type(row["key"][-1]) is not int or row["key"][-1] < 1
            or not isinstance(row["reason"], str) or not row["reason"].strip()
            or any(not isinstance(row[key], str) or not re.fullmatch(r"[a-f0-9]{64}", row[key])
                   for key in ("main_sha256", "candidate_sha256"))
        ):
            raise ValueError(f"{path}:{number}: invalid exception")
        key = tuple(row["key"])
        if key in exceptions:
            raise ValueError(f"{path}:{number}: duplicate exception: {key}")
        exceptions[key] = row
    return exceptions


def collector(rows, path):
    formats = {row.get("format") for row in rows.values()}
    if len(formats) != 1:
        raise ValueError(f"{path}: rows from more than one collector")
    return formats.pop()


def compare(baseline, candidate, exception_path=None):
    paths = (baseline, candidate)
    baseline, candidate = map(load, paths)
    if collector(baseline, paths[0]) != collector(candidate, paths[1]):
        raise ValueError("captures come from different collectors; use migrate to match them")
    exceptions = load_exceptions(exception_path) if exception_path else {}
    exact = recorded = unexpected = 0
    for key in sorted(baseline.keys() | candidate.keys()):
        old, new = baseline.get(key), candidate.get(key)
        if old == new:
            exact += 1
            continue
        exception = exceptions.get(key)
        if (exception is not None and exception["main_sha256"] == result_hash(old)
                and exception["candidate_sha256"] == result_hash(new)):
            del exceptions[key]
            recorded += 1
        else:
            unexpected += 1
            print(json.dumps({"key": key, "main": old, "candidate": new}, sort_keys=True))
    for key in sorted(exceptions):
        print(f"Stale or unused exception: {key}")
    print(f"{exact} exact matches, {recorded} recorded changes, {unexpected} unexpected changes")
    return bool(unexpected or exceptions)


def gate_rows(path, versions):
    rows = load(path)
    version = collector(rows, path)
    if version not in versions:
        raise ValueError(f"{path}: expected interface record format {versions}, got {version}")
    for key, row in rows.items():
        if not re.fullmatch(r"[a-f0-9]{64}", row["sha256"]):
            raise ValueError(f"{path}: invalid digest: {key}")
        if version == INTERFACE_FORMAT and type(row.get("refused")) is not bool:
            raise ValueError(f"{path}: missing or invalid refusal status: {key}")
    return rows


def compare_gate(baseline, candidate, exception_path, baseline_exception_path):
    one, four = (candidate / f"{count}.jsonl" for count in GATE_THREADS)
    main = gate_rows(baseline / one.name, (1, INTERFACE_FORMAT))
    head = gate_rows(one, (INTERFACE_FORMAT,))
    parallel = gate_rows(four, (INTERFACE_FORMAT,))
    worker_differences = 0
    for key in sorted(head.keys() | parallel.keys()):
        if head.get(key) != parallel.get(key):
            worker_differences += 1
            print(json.dumps({"problem": "worker difference", "key": key,
                              "1 worker": head.get(key), "4 workers": parallel.get(key)}))

    exceptions = load_exceptions(exception_path) if exception_path else {}
    old_exceptions = load_exceptions(baseline_exception_path) if baseline_exception_path else {}
    exceptions = {key: pin for key, pin in exceptions.items() if pin != old_exceptions.get(key)}
    exact = recorded = unexpected = 0
    for key, before in sorted(main.items()):
        after = head.get(key)
        pin = exceptions.get(key)
        if after is None:
            problem = "missing display call"
        elif before["sha256"] == after["sha256"]:
            if before.get("refused", after["refused"]) == after["refused"]:
                exact += 1
                continue
            problem = "inconsistent refusal status"
        elif after["refused"] and not before.get("refused", False):
            # Format 1 cannot prove that a different refusal already existed.
            problem = "new refusal"
        elif pin and (pin["main_sha256"], pin["candidate_sha256"]) == (
            before["sha256"], after["sha256"]
        ):
            del exceptions[key]
            recorded += 1
            continue
        else:
            problem = "unlisted digest change"
        unexpected += 1
        print(json.dumps({"problem": problem, "key": key, "main": before, "candidate": after}))

    for key in sorted(exceptions):
        print(f"Stale or unused exception: {key}")
    print(f"{exact} unchanged digests, {recorded} listed changes, {unexpected} unexpected changes, "
          f"{len(head.keys() - main.keys())} new calls, {worker_differences} worker differences")
    return bool(unexpected or exceptions or worker_differences)


def migrate(probe_path, interface_path):
    """Match a probe capture with an interface capture of the same revision.

    Calls pair in record order within each test and function. The pairs must
    agree on carriers and work, and options and digests must correspond one
    to one, so every probe result maps to exactly one interface result.
    """
    probe, interface = load(probe_path), load(interface_path)
    if collector(probe, probe_path) is not None or collector(interface, interface_path) is None:
        raise ValueError("migrate takes a probe capture, then an interface capture")

    def calls(rows):
        groups = collections.defaultdict(list)
        for key, row in rows.items():
            groups[key[:4]].append(row)
        return groups

    old, new = calls(probe), calls(interface)
    problems = []
    maps = {"options": ({}, {}), "sha256": ({}, {})}
    for call in sorted(old.keys() | new.keys()):
        before, after = old.get(call, []), new.get(call, [])
        if len(before) != len(after):
            problems.append(f"{call}: {len(before)} probe and {len(after)} interface records")
            continue
        for first, second in zip(before, after):
            pair = {"options": (first["key"][4], second["key"][4]),
                    "sha256": (first["sha256"], second["sha256"])}
            for field in ("carriers", "work"):
                if first[field] != second[field]:
                    problems.append(f"{first['key']}: {field} {first[field]} != {second[field]}")
            for field, (a, b) in pair.items():
                forward, backward = maps[field]
                if forward.setdefault(a, b) != b or backward.setdefault(b, a) != a:
                    problems.append(f"{first['key']}: {field} does not correspond one to one")
    for problem in problems:
        print(problem)
    print(f"{len(probe)} probe and {len(interface)} interface records, "
          f"{len(maps['sha256'][0])} distinct results, {len(problems)} problems")
    return bool(problems)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(required=True, dest="command")
    run = commands.add_parser("capture")
    run.add_argument("output", type=Path)
    run.add_argument("--collector", choices=("interface", "probe"), default="interface")
    run.add_argument("--threads", type=int, action="append", choices=THREAD_COUNTS,
                     help="capture only these worker counts; compare --gate uses 1 and 4")
    run.add_argument("--suite", action="append", choices=sorted(" ".join(SUITES.values()).split()))
    diff = commands.add_parser("compare")
    diff.add_argument("baseline", type=Path)
    diff.add_argument("candidate", type=Path)
    diff.add_argument("--exceptions", type=Path)
    diff.add_argument("--gate", action="store_true", help="require CI interface and worker checks")
    diff.add_argument("--baseline-exceptions", type=Path,
                      help="with --gate, ignore exception records already on main")
    commands.add_parser("affected", help="read NUL-delimited changed paths on stdin")
    match = commands.add_parser("migrate")
    match.add_argument("probe", type=Path)
    match.add_argument("interface", type=Path)
    args = parser.parse_args()

    if args.command == "affected":
        print(str(affects_display(os.fsdecode(sys.stdin.buffer.read()).split("\0"))).lower())
        return 0

    if args.command == "capture":
        args.output.mkdir(parents=True, exist_ok=True)
        root = Path.cwd()
        tests = root / "crates/catallax-features/tests"
        fixture = Path(__file__).with_name("ellipse-rim.json")
        with (tests / ".display_diff_ellipse.json").open("xb") as copy:
            copy.write(fixture.read_bytes())
        try:
            test = tests / "display_diff_ellipse.rs"
            with test.open("x") as source:
                source.write(ELLIPSE_TEST)
            try:
                with probes(root, args.collector):
                    for threads in args.threads or THREAD_COUNTS:
                        capture(root, args, args.output / f"{threads}.jsonl", threads)
            finally:
                test.unlink()
        finally:
            (tests / ".display_diff_ellipse.json").unlink()
        return 0

    if args.command == "migrate":
        return int(any([
            migrate(args.probe / f"{count}.jsonl", args.interface / f"{count}.jsonl")
            for count in THREAD_COUNTS
        ]))

    if args.gate:
        return int(compare_gate(args.baseline, args.candidate, args.exceptions, args.baseline_exceptions))
    if args.baseline_exceptions:
        parser.error("--baseline-exceptions requires --gate")

    first, *rest = (f"{count}.jsonl" for count in THREAD_COUNTS)
    failed = compare(args.baseline / first, args.candidate / first, args.exceptions)
    for side in (args.baseline, args.candidate):
        for name in rest:
            print(f"{side} {first} vs {name}")
            failed |= compare(side / first, side / name)
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
