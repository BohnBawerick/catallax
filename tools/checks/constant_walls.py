#!/usr/bin/env python3
"""Replay the scalar Shell and Thicken corpus against two catallax-verify builds.

On an offload checkout at the baseline, use capture with the existing Rust test
command. Then use corpus on the baseline checkout (or before grading changes)
and run the verifier driver on that same corpus with each binary:

    constant_walls.py capture CAPTURE -- cargo test [test target arguments]
    constant_walls.py corpus CAPTURE CORPUS
    control-catallax.py run CORPUS --binary BASELINE --out BASE_RUN --fresh
    control-catallax.py run CORPUS --binary CANDIDATE --out NEW_RUN --fresh
    constant_walls.py compare BASE_RUN/run.json NEW_RUN/run.json

The corpus includes result recipes and each scalar wall prefix. It requests
both B-rep and field outcomes, at the existing 100M request allowance.
Capture keeps import files through the last scalar wall. Corpus exits with
status 1 if any variant lacks a replayable import source.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
from textwrap import indent

ROOT = Path(__file__).resolve().parents[2]
DRIVER = ROOT / ".claude/skills/verify-catallax/scripts/control-catallax.py"
OPS = {"shell", "thicken"}
IMPORT_OPS = {"import_step", "import_iges", "import_mesh"}
FORMATS = {f"catallax-recipe/{version}" for version in range(3, 6)}

# These probes are only compiled in the disposable offload test checkout.
RECIPE_PROBE = """
        if let Some(out) = std::env::var_os("CATALLAX_CONSTANT_WALL_CAPTURE")
            && self.steps.iter().any(|s| matches!(
                s.op, Operation::Shell { .. } | Operation::Thicken { .. }
            ))
        {
            use std::hash::{Hash, Hasher};
            let mut hash = std::collections::hash_map::DefaultHasher::new();
            self.hash(&mut hash);
            let name = format!("recipe-{:016x}.json", hash.finish());
            std::fs::write(std::path::Path::new(&out).join(name), self.to_text())
                .expect("write captured wall recipe");
        }
"""
FILES_PROBE = """\
if let Some(out) = std::env::var_os("CATALLAX_CONSTANT_WALL_CAPTURE")
    && let Some(last) = recipe.steps.iter().rposition(|s| matches!(
        s.op, Operation::Shell { .. } | Operation::Thicken { .. }
    ))
{
    let mut seen = std::collections::BTreeSet::new();
    for step in &recipe.steps[..=last] {
        let file = match &step.op {
            Operation::ImportStep { file, .. }
            | Operation::ImportIges { file, .. }
            | Operation::ImportMesh { file, .. } => file,
            _ => continue,
        };
        if seen.insert(&file.sha256)
            && let Some(bytes) = files.get(&file.sha256)
            && catallax_import::sha256_hex(bytes) == file.sha256
        {
            let name = format!("source-{}", file.sha256);
            std::fs::write(std::path::Path::new(&out).join(name), bytes)
                .expect("write captured wall source");
        }
    }
}
"""


def capture(out, command):
    if not command:
        raise ValueError("capture needs a test command")
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    probes = (
        ("crates/catallax-recipe/src/recipe.rs", "        crate::check::check(self)", RECIPE_PROBE),
        ("crates/catallax-features/src/evaluate.rs", "    recipe.check().map_err(BuildError::Recipe)?;\n    catallax_numbers::reuse_operation_enclosures();", indent(FILES_PROBE, "    ")),
        ("crates/catallax-field/src/lib.rs", "        recipe.check().map_err(|e| {", indent(FILES_PROBE, "        ")),
    )
    saved = {}
    try:
        for filename, anchor, probe in probes:
            path = ROOT / filename
            subprocess.run(["git", "-C", ROOT, "diff", "--quiet", "HEAD", "--", filename], check=True)
            before = path.read_text()
            if before.count(anchor) != 1:
                raise ValueError(f"{filename}: expected one instrumentation anchor")
            saved[path] = before
            path.write_text(before.replace(anchor, probe + anchor))
        subprocess.run(command, cwd=ROOT, check=True,
                       env=os.environ | {"CATALLAX_CONSTANT_WALL_CAPTURE": str(out)})
    finally:
        for path, before in saved.items():
            path.write_text(before)


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def scalar_recipes(value, files=None):
    if isinstance(value, dict):
        if value.get("format") == "catallax-import-recipe/1":
            yield from scalar_recipes(value["recipe"], value["files"])
        elif value.get("format") in FORMATS:
            for i, step in enumerate(value.get("steps", []), 1):
                if step.get("op") in OPS:
                    yield value, files, i
        else:
            for child in value.values():
                yield from scalar_recipes(child, files)
    elif isinstance(value, list):
        for child in value:
            yield from scalar_recipes(child, files)


def corpus(captured, out):
    cases = {}
    skipped = set()
    sources = {p.name.removeprefix("source-"): p
               for p in captured.glob("source-*")}
    tracked = subprocess.check_output(["git", "-C", ROOT, "ls-files", "-z", "*.json", "*.jsonl"])
    captured_recipes = sorted(captured.glob("recipe-*.json"))
    paths = [ROOT / os.fsdecode(name) for name in tracked.split(b"\0") if name]
    for path in paths + captured_recipes:
        raw = path.read_text()
        if "catallax-recipe/" not in raw or not any(f'"{op}"' in raw for op in OPS):
            continue
        values = [json.loads(line) for line in raw.splitlines() if line.strip()] if path.suffix == ".jsonl" else [json.loads(raw)]
        for value in values:
            for recipe, files, prefix in scalar_recipes(value):
                steps = recipe["steps"]
                ids = {step["id"] for step in steps[:prefix]}
                if recipe["result"] == steps[prefix - 1]["id"] and prefix == len(steps):
                    variants = [recipe]
                else:
                    wall = recipe | {"steps": steps[:prefix], "result": steps[prefix - 1]["id"]}
                    if "materials" in recipe:
                        wall["materials"] = [m for m in recipe["materials"] if m["region"] in ids]
                    variants = [recipe, wall]
                for variant in variants:
                    required = {step["file"]["sha256"] for step in variant["steps"]
                                if step.get("op") in IMPORT_OPS}
                    source_files = {}
                    missing = []
                    for digest in sorted(required):
                        text = (files or {}).get(digest)
                        if text is None and digest in sources:
                            try:
                                text = sources[digest].read_bytes().decode()
                            except UnicodeDecodeError:
                                pass  # The verifier's import wrapper accepts text sources.
                        if isinstance(text, str) and hashlib.sha256(text.encode()).hexdigest() == digest:
                            source_files[digest] = text
                        else:
                            missing.append(digest)
                    if missing:
                        key = (path, variant["name"], variant["result"], tuple(missing))
                        if key not in skipped:
                            skipped.add(key)
                            print(f"{path}: {variant['name']} result {variant['result']} missing "
                                  f"replayable imports {', '.join(missing)}", file=sys.stderr)
                        continue
                    case = ({"format": "catallax-import-recipe/1", "recipe": variant,
                             "files": source_files} if source_files else variant)
                    cases[encoded(case)] = case
    if not cases and not skipped:
        raise ValueError("no scalar wall recipes found")
    out.mkdir(parents=True, exist_ok=False)
    for key, case in sorted(cases.items()):
        digest = hashlib.sha256(key.encode()).hexdigest()
        for route in ("recipe", "field"):
            question = {"format": "catallax-question/1", "source": route,
                        "problem": case, "max_work": 100_000_000}
            (out / f"{digest}-{route}.json").write_text(encoded(question) + "\n")
    print(f"{len(cases)} scalar recipe variants, {len(captured_recipes)} captured inputs, {len(sources)} captured import files")
    return int(bool(skipped))


def compare(baseline, candidate):
    spec = importlib.util.spec_from_file_location("driver", DRIVER)
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    result = driver.compare(baseline, candidate)
    a, b = [json.loads(path.read_text()) for path in (baseline, candidate)]
    for i, (x, y) in enumerate(zip(a["cases"], b["cases"])):
        paths = [p.parent / c["path"] / "operations.json"
                 for p, c in ((baseline, x), (candidate, y))]
        logs = [driver.without_clocks(json.loads(path.read_text())) for path in paths]
        if x["response"]["work"] != y["response"]["work"] or logs[0] != logs[1]:
            result["failures"].append({"index": i, "id": x["id"], "cause": "work-or-operation-log"})
    result["equal"] = not result["failures"]
    print(json.dumps(result, indent=2))
    return 0 if result["equal"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="mode", required=True)
    cmd = subs.add_parser("capture")
    cmd.add_argument("out", type=Path)
    cmd.add_argument("command", nargs=argparse.REMAINDER)
    cmd = subs.add_parser("corpus")
    cmd.add_argument("captured", type=Path)
    cmd.add_argument("out", type=Path)
    cmd = subs.add_parser("compare")
    cmd.add_argument("baseline", type=Path)
    cmd.add_argument("candidate", type=Path)
    args = parser.parse_args()
    if args.mode == "capture":
        command = args.command[1:] if args.command[:1] == ["--"] else args.command
        return capture(args.out, command)
    if args.mode == "corpus":
        return corpus(args.captured, args.out)
    return compare(args.baseline, args.candidate)


if __name__ == "__main__":
    sys.exit(main())
