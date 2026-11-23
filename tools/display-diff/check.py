#!/usr/bin/env python3
"""Compare public display results on two revisions.

Run capture inside each disposable checkout with the tool and its fixture from
the same revision. Capture writes one JSONL file per worker count.
"""
import argparse
import collections
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
        "fine_mesh generator_fillet_display p5_display preview_mesh "
        "rational_lofts rational_patches rolling_fillets spherical_corner_mesh "
        "spline_display spline_sweeps straight_sweep_display"
    ),
    "catallax-features": (
        "algebraic_extrusions algebraic_revolutions bezier_profiles certified_prism "
        "coaxial_fillets display_diff_ellipse draft_annular draft_oblique draft_preserved_bore "
        "draft_preserved_boss equal_crossing_fillets fillet_over_fillet "
        "fillet_sphere_ring g2_fillets m11_spline_fillet oblique_cone_fillets "
        "oblique_fillets ruled_edge_rows saddle_fillets slot_fillets strip_ends "
        "variable_fillets"
    ),
    "catallax-field": (
        "patch_shell rolling_conics rolling_ends rolling_over rolling_quadrics "
        "rolling_tori spline_spine unproved_preview"
    ),
}

RUST_PROBE = r'''
fn differential_record(
    function: &str,
    before: catallax_numbers::OperationWork,
    tags: &str,
    options: &impl std::fmt::Debug,
    result: &impl std::fmt::Debug,
) {
    use std::io::Write;
    use std::process::{Command, Stdio};
    let after = catallax_numbers::current_operation_work();
    let work = [
        after.used.wrapping_sub(before.used),
        after.refunded.wrapping_sub(before.refunded),
        after.uncharged.wrapping_sub(before.uncharged),
        u64::from(after.exhausted),
    ];
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
        "\nmeshdiff\t{}\t{function}\t{tags}\t{options:?}\t{work:?}\t{digest}",
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
            + "    let before = catallax_numbers::current_operation_work();\n"
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


def capture(root, args, output, threads):
    originals = {}
    try:
        for filename, names in PROBES.items():
            path = root / "crates" / filename
            originals[path] = path.read_text()
            path.write_text(instrument(originals[path], names))
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
    finally:
        for path, source in originals.items():
            path.write_text(source)


def load(path):
    result = {}
    for line in path.read_text().splitlines():
        row = json.loads(line)
        key = tuple(row["key"])
        assert key not in result, key
        result[key] = row
    if not result:
        raise ValueError(f"empty capture: {path}")
    return result


def compare(baseline, candidate, exception_path=None):
    baseline, candidate = map(load, (baseline, candidate))
    exceptions = {}
    if exception_path:
        for line in exception_path.read_text().splitlines():
            row = json.loads(line)
            key = tuple(row["key"])
            assert key not in exceptions and row["reason"], row
            exceptions[key] = row
    differences = set()
    unexpected = 0
    for key in sorted(baseline.keys() | candidate.keys()):
        if baseline.get(key) == candidate.get(key):
            continue
        differences.add(key)
        row = {"key": key, "main": baseline.get(key), "candidate": candidate.get(key)}
        exception = exceptions.get(key)
        if exception is None or any(exception[field] != row[field] for field in ("main", "candidate")):
            unexpected += 1
            print(json.dumps(row, sort_keys=True))
    stale = exceptions.keys() - differences
    for key in sorted(stale):
        print(f"Unused exception: {key}")
    print(f"{len(baseline.keys() | candidate.keys()) - len(differences)} exact matches, "
          f"{len(differences) - unexpected} recorded changes, {unexpected} unexpected changes")
    return bool(unexpected or stale)


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
