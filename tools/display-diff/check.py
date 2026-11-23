#!/usr/bin/env python3
"""Compare public display results on two revisions.

Run capture inside each disposable checkout with the same copy of this tool and
its fixtures. Capture writes one JSONL file per worker count.
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
        # Anchor at the function's own declaration. Its closing brace is unindented.
        pattern = rf"(^[a-z() ]*fn {name}\([^{{]*\{{\n)(.*?)(^\}})"
        match = re.search(pattern, source, re.MULTILINE | re.DOTALL)
        if match is None:
            raise ValueError(f"no declaration for {name}")
        tags = '    let tags = String::new();\n' if name == "mesh_plan_inner" else CARRIERS
        first, last = match.span(2)
        old = match.group(2)
        # display_mesh_body establishes its own allowance on standalone calls.
        # Record the inner invocation within that allowance exactly once.
        prefix = ""
        if name == "display_mesh_body":
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
def probes(root):
    originals = {}
    try:
        for filename, names in PROBES.items():
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
                        if not line.startswith("meshdiff\t"):
                            continue
                        case, function, tags, options, work, digest = line.rstrip().split("\t")[1:]
                        assert suite in selected, (suite, line)
                        key = (package, suite, case, function, options)
                        counts[key] += 1
                        row = {"key": [*key, counts[key]], "carriers": tags,
                               "work": json.loads(work), "sha256": digest}
                        print(json.dumps(row, sort_keys=True), file=out, flush=True)
                    code = run.wait()
                if code:
                    log.seek(0)
                    for line in log:
                        if not line.startswith("meshdiff\t"):
                            sys.stderr.write(line)
                print(json.dumps({"package": package, "exit": code, "calls": sum(counts.values())}),
                      file=sys.stderr, flush=True)
                if code:
                    raise SystemExit(code)


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


def compare(baseline, candidate, exception_path=None):
    baseline, candidate = map(load, (baseline, candidate))
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(required=True, dest="command")
    run = commands.add_parser("capture")
    run.add_argument("output", type=Path)
    run.add_argument("--suite", action="append", choices=sorted(" ".join(SUITES.values()).split()))
    diff = commands.add_parser("compare")
    diff.add_argument("baseline", type=Path)
    diff.add_argument("candidate", type=Path)
    diff.add_argument("--exceptions", type=Path)
    args = parser.parse_args()

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
                with probes(root):
                    for threads in THREAD_COUNTS:
                        capture(root, args, args.output / f"{threads}.jsonl", threads)
            finally:
                test.unlink()
        finally:
            (tests / ".display_diff_ellipse.json").unlink()
        return 0

    first, *rest = (f"{count}.jsonl" for count in THREAD_COUNTS)
    failed = compare(args.baseline / first, args.candidate / first, args.exceptions)
    for side in (args.baseline, args.candidate):
        for name in rest:
            print(f"{side} {first} vs {name}")
            failed |= compare(side / first, side / name)
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
