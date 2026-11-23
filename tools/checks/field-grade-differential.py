#!/usr/bin/env python3
"""Capture and compare existing lattice and gyroid output tests on two revisions.

Run capture with this same script on each revision, using the same Rust toolchain,
target architecture and thread count. Capture temporarily adds probes to the source
in the supplied checkout and restores it on exit. Use a disposable build checkout.
On macOS, export SDKROOT with xcrun --sdk macosx --show-sdk-path before capture
so the Python launcher and Cargo use the selected Xcode SDK.

    field-grade-differential.py capture CHECKOUT BASELINE
    field-grade-differential.py capture CHECKOUT CANDIDATE
    field-grade-differential.py compare BASELINE CANDIDATE
"""

import argparse
from collections import defaultdict
import filecmp
import hashlib
import os
from pathlib import Path
import subprocess
import tempfile


FIELD_TESTS = (
    "graded_fields", "lattice_stage", "lattice", "printing",
    "adaptive_print", "mass", "cubic_moments",
)
CLI_TESTS = ("graded_fields", "m19_print_verdict", "verification")
FIELD = "crates/catallax-field/src/"
CLI = "crates/catallax-cli/src/"

# The CLI suite also checks unknown grade syntax. Its parser diagnostics list
# supported variants, so record only the existing constant and linear inputs.
LEGACY_REQUEST = """request.problem.get("recipe").unwrap_or(&request.problem)
        .get("steps").and_then(serde_json::Value::as_array).is_none_or(|steps| {
            steps.iter().all(|step| {
                let key = match step["op"].as_str() {
                    Some("lattice") => "radius",
                    Some("gyroid") => "level",
                    _ => return true,
                };
                step[key].get("constant").is_some() || step[key].get("linear").is_some()
            })
        })"""

PROBES = {
    FIELD + "lattice_stage.rs": {
        "    pub fn export": ("stage-export", "write!(out, \"{result:?}\")?;"),
        "    pub fn labelled_stl": (
            "stage-stl",
            """match result {
                Ok((bytes, output)) => {
                    writeln!(out, "{output:?}")?;
                    out.write_all(bytes)?;
                }
                Err(error) => write!(out, "{error:?}")?,
            }""",
        ),
        "pub fn preview_plan": ("preview", "write!(out, \"{result:?}\")?;"),
    },
    FIELD + "printing.rs": {
        "pub fn mesh_plan(": ("print-mesh", "write!(out, \"{result:?}\")?;"),
        "pub fn mesh_plan_connected": ("connected-mesh", "write!(out, \"{result:?}\")?;"),
    },
    FIELD + "printing/brep.rs": {
        "pub fn field(": ("qualified-field-print", "write!(out, \"{result:?}\")?;"),
    },
    FIELD + "printing/adaptive.rs": {
        "pub fn mesh_with_workers": ("adaptive-print", "write!(out, \"{result:?}\")?;"),
    },
    FIELD + "mass.rs": {
        "pub fn properties(": ("mass", "write!(out, \"{result:?}\")?;"),
    },
    CLI + "verification.rs": {
        "pub fn answer(": (
            "request",
            """match result {
                Ok(response) => {
                    fn stable(value: &mut serde_json::Value) {
                        if let Some(object) = value.as_object_mut() {
                            object.remove("wall_ms");
                            object.values_mut().for_each(stable);
                        } else if let Some(array) = value.as_array_mut() {
                            array.iter_mut().for_each(stable);
                        }
                    }
                    let mut response = response.clone();
                    stable(&mut response);
                    if let Some(kernel) = response["manifest"]["kernel"].as_object_mut() {
                        kernel.remove("commit");
                    }
                    if let Some(display) = response.pointer_mut("/manifest/face_measurements/display")
                        .and_then(serde_json::Value::as_object_mut) {
                        display.remove("first_mesh_ms");
                    }
                    writeln!(out, "{response}")?;
                }
                Err(error) => writeln!(out, "{error:?}")?,
            }
            for name in ["print.stl", "display.obj", "answer.step"] {
                match std::fs::read(request.out.join(name)) {
                    Ok(bytes) => {
                        writeln!(out, "{name} {}", bytes.len())?;
                        out.write_all(&bytes)?;
                    }
                    Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                        writeln!(out, "{name} absent")?;
                    }
                    Err(error) => return Err(error),
                }
            }""",
            LEGACY_REQUEST,
        ),
    },
}


def instrument(source, signature, name, write_result, include="true"):
    start = source.index("{", source.index("\n" + signature))
    indent = signature[:len(signature) - len(signature.lstrip())]
    end = source.index("\n" + indent + "}", start)
    body = source[start + 1:end]
    return source[:start + 1] + f"""
        let result = (|| {{ {body}
        }})();
        if let Some(dir) = std::env::var_os("CATALLAX_GRADE_CAPTURE") && {include} {{
            use std::io::Write as _;
            static NEXT: std::sync::atomic::AtomicUsize =
                std::sync::atomic::AtomicUsize::new(0);
            let thread = std::thread::current();
            let test = thread.name().unwrap_or("unnamed");
            let index = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let process = std::path::Path::new(&dir).join(std::process::id().to_string());
            std::fs::create_dir_all(&process).expect("create field grade capture directory");
            let file = process.join(format!("{{test}}-{name}-{{index:05}}"));
            let save = || -> std::io::Result<()> {{
                let mut out = std::io::BufWriter::new(std::fs::File::create(file)?);
                let result = &result;
                {write_result}
                out.flush()
            }};
            save().expect("write field grade baseline capture");
        }}
        result""" + source[end:]


def capture(checkout, output):
    checkout = checkout.resolve()
    output = output.resolve()
    output.mkdir()
    originals = {}
    try:
        for name, probes in PROBES.items():
            path = checkout / name
            source = path.read_text()
            originals[path] = source
            for signature, probe in probes.items():
                source = instrument(source, signature, *probe)
            path.write_text(source)

        with tempfile.TemporaryDirectory(prefix="grade-differential-") as scratch:
            for crate, tests in (("field", FIELD_TESTS), ("cli", CLI_TESTS)):
                for test in tests:
                    case = output / f"{crate}-{test}"
                    case.mkdir()
                    filters = (
                        ["the_field_adapter_shows_a_lattice_as_an_unproved_preview"]
                        if crate == "cli" and test == "verification" else []
                    )
                    subprocess.run(
                        ["cargo", "test", "--locked", "--release", "-p", f"catallax-{crate}",
                         "--test", test, "--", "--test-threads=1", "--include-ignored", *filters],
                        cwd=checkout,
                        env={**os.environ, "CATALLAX_GRADE_CAPTURE": str(case),
                             "TMPDIR": scratch},
                        check=True,
                    )
                    if not any(case.iterdir()):
                        raise ValueError(f"No results in {case}")
    finally:
        for path, source in originals.items():
            path.write_text(source)


def captures(root):
    cases = defaultdict(list)
    for process in root.glob("*/*"):
        files = sorted(process.iterdir())
        signature = []
        for path in files:
            digest = hashlib.sha256()
            with path.open("rb") as file:
                for block in iter(lambda: file.read(1 << 20), b""):
                    digest.update(block)
            signature.append((path.name, digest.digest()))
        cases[process.parent.name].append((signature, files))
    return cases


def compare(baseline, candidate):
    expected = captures(baseline)
    actual = captures(candidate)
    if not expected or expected.keys() != actual.keys():
        raise SystemExit(f"Capture suites differ: {sorted(expected)} vs {sorted(actual)}")
    total = 0
    for case in sorted(expected):
        left = sorted(expected[case], key=lambda process: process[0])
        right = sorted(actual[case], key=lambda process: process[0])
        if len(left) != len(right):
            raise SystemExit(f"{case}: capture process counts differ")
        count = 0
        for (_, baseline_files), (_, candidate_files) in zip(left, right):
            if len(baseline_files) != len(candidate_files):
                raise SystemExit(f"{case}: capture counts differ")
            for a, b in zip(baseline_files, candidate_files):
                if a.name != b.name or not filecmp.cmp(a, b, shallow=False):
                    raise SystemExit(f"{case}: changed result {a.name} vs {b.name}")
                count += 1
        total += count
        print(f"{case}: {count} exact matches")
    print(f"{total} exact matches")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("capture", "compare"))
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    args = parser.parse_args()
    if args.mode == "capture":
        capture(args.left, args.right)
    else:
        compare(args.left, args.right)


if __name__ == "__main__":
    main()
