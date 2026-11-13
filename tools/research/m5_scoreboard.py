#!/usr/bin/env python3
"""Prepare private M5 inputs from the pinned R8 corpus and unchanged P0 recipes."""

import argparse
import copy
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from r8_step_census import ARCHIVE_BYTES, ARCHIVE_MD5, SEED, choose, digest, restrict


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE_SHA256 = "52e6dd1b6fa38e3cd99af59b662370829129540030975919c3e256dce6ad1dbe"
# First eleven complete admissions below 1M import work in the published R8 census.
MODEL_IDS = ("00000131", "00000173", "00000190", "00000405", "00000511",
             "00000549", "00000605", "00000632", "00000633", "00000695", "00000821")
LENGTHS = {"block": ("size",), "cylinder": ("radius", "height"),
           "sphere": ("radius",), "cone": ("base_radius", "top_radius", "height"),
           "torus": ("major", "minor")}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def quaternion_product(a, b):
    w, x, y, z = map(Fraction, a)
    v, i, j, k = map(Fraction, b)
    return [str(n) for n in (w*v-x*i-y*j-z*k, w*i+x*v+y*k-z*j,
                             w*j-x*k+y*v+z*i, w*k+x*j-y*i+z*v)]


def p0_copy(recipe, variant):
    """Change geometry exactly, or move only the periodic charts of round stock."""
    result = copy.deepcopy(recipe)
    result["name"] += f" / {variant}"
    result.pop("description", None)
    for step in result["steps"]:
        if step["op"] == "boolean":
            continue
        dimensions = LENGTHS[step["op"]]
        pose = step["pose"]
        if variant == "translated":
            pose["translation"] = [str(Fraction(x) + delta) for x, delta in
                                   zip(pose["translation"], (101, -37, 19))]
        elif variant == "scaled":
            for key in dimensions:
                value = step[key]
                step[key] = ([str(2 * Fraction(x)) for x in value] if isinstance(value, list)
                             else str(2 * Fraction(value)))
            pose["translation"] = [str(2 * Fraction(x)) for x in pose["translation"]]
        elif variant == "seam-shifted":
            if step["op"] != "block":
                pose["rotation"] = quaternion_product(pose["rotation"], (3, 0, 0, 4))
        else:
            raise ValueError(f"unknown P0 copy: {variant}")
    return result


def prepare_p0(scratch):
    inputs = scratch / "p0-inputs"
    inputs.mkdir()
    records = []
    for path in sorted((ROOT / "parts/P0").glob("*.json")):
        source = json.loads(path.read_text())
        for variant in ("base", "translated", "scaled", "seam-shifted"):
            name = f"{path.stem}-{variant}.json"
            if variant == "base":
                shutil.copyfile(path, inputs / name)
            else:
                write_json(inputs / name, p0_copy(source, variant))
            records.append({"input": name, "source": path.name,
                            "variant": variant, "sha256": digest(inputs / name)})
    write_json(scratch / "p0-selection.json", records)


def import_wrapper(data, model_id):
    sha = hashlib.sha256(data).hexdigest()
    return {"format": "catallax-import-recipe/1", "files": {sha: data.decode("utf-8")},
            "recipe": {"format": "catallax-recipe/3", "name": f"ABC {model_id}",
                       "steps": [{"id": "imported", "op": "import_step",
                                  "file": {"name": "part.step", "sha256": sha},
                                  "pose": {"rotation": ["1", "0", "0", "0"],
                                           "translation": ["0", "0", "0"]}}],
                       "result": "imported"}}


def prepare_real(archive, scratch):
    if (archive.stat().st_size != ARCHIVE_BYTES or digest(archive, "md5") != ARCHIVE_MD5
            or digest(archive) != ARCHIVE_SHA256):
        raise ValueError("pinned R8 archive integrity mismatch")
    entries = subprocess.check_output(["bsdtar", "-tf", str(archive)], text=True).splitlines()
    selected, total = choose(entries, 3000, SEED)
    selected = restrict(selected, MODEL_IDS)
    extracted = scratch / "extracted"
    inputs = scratch / "real-inputs"
    extracted.mkdir()
    inputs.mkdir()
    members = scratch / "extract-list-private.txt"
    members.write_text("".join(entry + "\n" for _, entry in selected))
    records = []
    try:
        subprocess.run(["bsdtar", "-xf", str(archive), "-C", str(extracted),
                        "-T", str(members)], check=True)
        for model_id, member in selected:
            path = extracted / member
            if not path.is_file() or path.is_symlink():
                raise ValueError("missing or non-regular corpus model")
            data = path.read_bytes()
            wrapper = import_wrapper(data, model_id)
            write_json(inputs / f"{model_id}.json", wrapper)
            records.append({"model_id": model_id, "source_sha256": hashlib.sha256(data).hexdigest()})
    finally:
        shutil.rmtree(extracted)
        members.unlink()
    write_json(scratch / "real-selection.json", {"archive_sha256": ARCHIVE_SHA256,
               "sample_seed": SEED, "chunk_models": total, "models": records})


def publish(evidence, output):
    driver = ROOT / ".claude/skills/verify-catallax/scripts/control-catallax.py"
    sys.path.insert(0, str(driver.parent))
    spec = importlib.util.spec_from_file_location("m5_driver", driver)
    control = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(control)

    def compare_runs(left, right):
        # Collected artifacts can move hosts without changing their saved requests.
        with tempfile.TemporaryDirectory(prefix="publish-", dir=evidence) as folder:
            paths = []
            for index, path in enumerate((left, right)):
                value = json.loads(path.read_text())
                for case in value["cases"]:
                    case["path"] = str(path.parent / "cases" / case["key"])
                target = Path(folder) / f"{index}.json"
                write_json(target, value)
                paths.append(str(target))
            return control.compare(*paths)

    report = {"format": "catallax-m5-scoreboard/1", "max_work": 1_000_000,
              "comparison": "not-requested", "timing": "shared machine, not an M26 timing pass",
              "source_archive_sha256": ARCHIVE_SHA256, "real_part_ids": list(MODEL_IDS),
              "unmeasured": ["selectable-preview", "display", "parameter-edit", "background-proof",
                             "requested-measurement", "independent-step-check", "field-comparison"],
              "groups": {}}
    for group, count in (("p0", 32), ("real", len(MODEL_IDS))):
        baseline = evidence / f"{group}-baseline/run.json"
        first = evidence / f"{group}-1/run.json"
        run = json.loads(first.read_text())
        assert len(run["cases"]) == count
        assert compare_runs(first, baseline)["equal"], "baseline answers changed"
        for workers in (1, 2, 4):
            path = evidence / f"{group}-{workers}/run.json"
            config = json.loads(path.with_name("config.json").read_text())
            current = json.loads(path.read_text())
            assert current["execution"]["fresh"] is True
            assert current["execution"]["workers"] == workers
            assert str(config["kernel_threads"]) == str(workers)
            assert config["max_work"] == 1_000_000
            assert current["config"]["binary_sha256"] == run["config"]["binary_sha256"]
            assert all(code == 0 for code in current["execution"]["worker_exits"].values())
            assert compare_runs(first, path)["equal"], "worker answers or artifacts differ"
            assert all(control.without_clocks(a["response"]["operation_log"]) ==
                       control.without_clocks(b["response"]["operation_log"])
                       for a, b in zip(run["cases"], current["cases"])), "operation logs differ"
        rows = []
        for case in run["cases"]:
            directory = first.parent / "cases" / case["key"]
            request = json.loads((directory / "request.json").read_text())
            problem = request["problem"]
            response = json.loads((directory / "response.json").read_text())
            manifest = json.loads((directory / "manifest.json").read_text())
            assert response["manifest"] == manifest
            assert manifest["status"] in ("answered", "refused")
            assert request["max_work"] == response["work"]["limit"] == 1_000_000
            row = {"name": problem.get("name") or problem["recipe"]["name"],
                   "status": manifest["status"], "exact": manifest["exact"],
                   "result": manifest["result"], "work": response["work"],
                   "request_wall_ms": response["wall_ms"],
                   "answer_sha256": hashlib.sha256(control.canonical(control.comparable(manifest)).encode()).hexdigest(),
                   "refusal": ({key: value for key, value in manifest["refusal"].items()
                                if key in ("code", "step", "stage")} if manifest["refusal"] else None),
                   "operations": [{key: op[key] for key in ("id", "operation", "status", "work", "wall_ms")}
                                  for op in response["operation_log"]["operations"]]}
            row["blocked_operation"] = next((op for op in reversed(row["operations"])
                                             if op["status"] != "answered"), None)
            step = directory / "answer.step"
            if manifest["step"]:
                assert manifest["status"] == "answered" and step.is_file() and step.stat().st_size
                row["step_sha256"] = digest(step)
                row["step_bytes"] = step.stat().st_size
            else:
                assert not step.exists(), "refused request left STEP output"
            rows.append(row)
        report["production_revision"] = run["cases"][0]["response"]["manifest"]["kernel"]["commit"]
        report["baseline_revision"] = json.loads(baseline.read_text())["cases"][0]["response"]["manifest"]["kernel"]["commit"]
        report["binary_sha256"] = run["config"]["binary_sha256"]
        report["groups"][group] = {"cases": rows, "baseline_equal": True,
                                  "workers_and_threads": [1, 2, 4], "answers_work_and_artifacts_identical": True}
    assert sum(row["status"] == "answered" and row["exact"]
               for row in report["groups"]["real"]["cases"]) == 10
    write_json(output.with_suffix(".json"), report)
    lines = ["# M5 P0 and real-part scoreboard", "",
             "Generated by `tools/research/m5_scoreboard.py publish` from retained CLI artifacts.", "",
             f"Kernel source: `{report['production_revision']}`. Baseline main: `{report['baseline_revision']}`. Each complete request has the unchanged 1,000,000-work scoreboard allowance.", "",
             "All answers and artifacts match the pre-change main snapshot. Fresh runs with 1, 2 and 4 workers, each using the same number of kernel threads, have identical manifests, deterministic work and STEP bytes.", "",
             "Clocks below are end-to-end proof-request observations on a shared machine. They do not qualify M26. Selectable preview, display, parameter edits, separate background proof, requested measurements and independent STEP validation were not requested. Field comparison remains on demand and was not requested.", "",
             "## P0", "",
             "The eight authored recipes are unchanged. Each has three generated copies: translation by (101, -37, 19) mm, uniform scale 2, and a local z-axis chart rotation (3, 0, 0, 4) on round stock. The last preserves geometry and moves periodic seams. Operation ids and Boolean steps are unchanged.", "",
             "| Case | Result | Work | Request ms |", "| --- | --- | ---: | ---: |"]
    for row in report["groups"]["p0"]["cases"]:
        outcome = "answered, exact" if row["status"] == "answered" and row["exact"] else row["status"]
        if row["refusal"]:
            block = row["blocked_operation"]
            location = f"{block['id']}/{block['operation']}" if block else "unknown operation"
            outcome = f"`{row['refusal']['code']}` at `{location}`"
        lines.append(f"| {row['name']} | {outcome} | {row['work']['used']:,} | {row['request_wall_ms']:.3f} |")
    lines += ["", "## Ten real parts", "",
              "The first ten complete admissions below 1M import work in the published R8 census were selected by public model id. Nine complete requests answered; 00000549 exhausted the same allowance in the complete recipe import step, before export. The next prior admission, 00000821, was added to reach ten answering parts. All eleven candidates remain in the table. These are capability controls from ABC chunk 0000, not a representative admission-rate sample or evidence that a particular user makes them. The original census remains unchanged.", "",
              "The source and licence scope are the [R8 census](R8-public-step-census.md#source-and-licence-scope). Per-model source exceptions remain unresolved. Only ids, counts, outcomes, work and artifact hashes are published. Raw evidence retained privately.", "",
              "| ABC model | Result | Faces | Work | Request ms |", "| --- | --- | ---: | ---: | ---: |"]
    for row in report["groups"]["real"]["cases"]:
        outcome = "answered, exact" if row["status"] == "answered" and row["exact"] else row["status"]
        if row["refusal"]:
            block = row["blocked_operation"]
            location = f"{block['id']}/{block['operation']}" if block else "unknown operation"
            outcome = f"`{row['refusal']['code']}` at `{location}`"
        faces = row["result"]["faces"] if row["result"] else "unmeasured"
        lines.append(f"| {row['name']} | {outcome} | {faces} | {row['work']['used']:,} | {row['request_wall_ms']:.3f} |")
    lines += ["", "[Machine-readable results](M5-scoreboard.json) retain result counts, per-operation work and clocks, exact-answer hashes and hashes of the actual STEP files. An answered recipe includes checked kernel construction and STEP export. It does not establish independent import validity, mass properties or SolidDraft acceptance.", "",
              "## Reproduce", "",
              "Use the checksum-pinned archive from the R8 reproduction instructions. Set `M5_PROOF` to fresh private storage outside the checkout, `M5_ARCHIVE` to that archive, and `M5_BINARY` and `M5_BASELINE` to preserved candidate and baseline verifiers. Build heavy binaries on the work PC. The preparation tool verifies archive size, MD5 and SHA-256 and deletes extracted model copies. It keeps private, content-addressed request wrappers for the CLI. Delete the archive and disposable input wrappers after collecting evidence.", "",
              "```sh",
              'python3 -B tools/research/m5_scoreboard.py p0 --scratch "$M5_PROOF/inputs"',
              'python3 -B tools/research/m5_scoreboard.py real --archive "$M5_ARCHIVE" --scratch "$M5_PROOF/inputs"',
              'for group in p0 real; do',
              '  CATALLAX_THREADS=1 python3 -B .claude/skills/verify-catallax/scripts/control-catallax.py run "$M5_PROOF/inputs/$group-inputs" --binary "$M5_BASELINE" --out "$M5_PROOF/$group-baseline" --source recipe --max-work 1000000 --workers 1 --fresh --json',
              '  for workers in 1 2 4; do',
              '    CATALLAX_THREADS="$workers" python3 -B .claude/skills/verify-catallax/scripts/control-catallax.py run "$M5_PROOF/inputs/$group-inputs" --binary "$M5_BINARY" --out "$M5_PROOF/$group-$workers" --source recipe --max-work 1000000 --workers "$workers" --fresh --json',
              '  done', 'done',
              'python3 -B tools/research/m5_scoreboard.py publish --scratch "$M5_PROOF" --output docs/milestones/M5-scoreboard',
              '```', ""]
    output.with_suffix(".md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("p0", "real", "publish"))
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--output", type=Path, help="publication path stem for publish")
    args = parser.parse_args()
    scratch = args.scratch.resolve()
    if scratch.is_relative_to(ROOT):
        parser.error("private inputs must stay outside the checkout")
    scratch.mkdir(parents=True, exist_ok=True)
    if args.kind == "p0":
        prepare_p0(scratch)
    elif args.kind == "publish":
        if args.output is None:
            parser.error("publish requires --output")
        publish(scratch, args.output)
    elif args.archive is None:
        parser.error("real parts require --archive")
    else:
        prepare_real(args.archive.resolve(), scratch)


if __name__ == "__main__":
    main()
