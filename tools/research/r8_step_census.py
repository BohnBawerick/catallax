#!/usr/bin/env python3
"""Sample one pinned ABC STEP chunk, observe import, and remove model copies.

Only counts and ABC model ids enter the public CSV. Full import diagnostics
and archive member names stay in the supplied private scratch directory.
"""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import random
import re
import shutil
import subprocess
import time


ARCHIVE = "abc_0000_step_v00.7z"
ARCHIVE_BYTES = 1_594_129_754
ARCHIVE_MD5 = "695388be7a278c7798c8c8ae239772ac"
SEED = 20261111
KINDS = ("plane", "cylinder", "cone", "sphere", "torus", "bspline",
         "rational_bspline", "other", "unresolved")


def digest(path, algorithm="sha256"):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, algorithm).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def choose(entries, count, seed):
    models = {}
    for entry in entries:
        path = PurePosixPath(entry)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("unsafe archive member")
        if not entry.lower().endswith((".step", ".stp")):
            continue
        match = re.fullmatch(r"(\d{8})_[a-zA-Z0-9_]+\.(?:step|stp)", path.name)
        if not match:
            raise ValueError("unexpected STEP member naming")
        model_id = match[1]
        if model_id in models:
            raise ValueError("more than one STEP file for a model id")
        models[model_id] = entry
    ids = sorted(models)
    if count > len(ids):
        raise ValueError("sample exceeds available model ids")
    selected = sorted(random.Random(seed).sample(ids, count))
    return [(model_id, models[model_id]) for model_id in selected], len(ids)


def observe(binary, path, model_id, raw, timeout):
    stdout = raw / f"{model_id}.jsonl"
    stderr = raw / f"{model_id}.stderr"
    started = time.monotonic()
    with stdout.open("wb") as out, stderr.open("wb") as err:
        try:
            completed = subprocess.run([str(binary), str(path)], stdout=out,
                                       stderr=err, timeout=timeout, check=False)
            outcome = "process_error" if completed.returncode else "missing_result"
            returncode = completed.returncode
        except subprocess.TimeoutExpired:
            outcome, returncode = "observation_timeout", None
    record = {"model_id": model_id, "outcome": outcome,
              "coverage_complete": False, "face_kinds": {},
              "all_face_kinds_buildable": None, "needs_spline_faces": None,
              "fillet_count": None, "chamfer_count": None,
              "feature_count_status": "unknown_no_feature_history_in_step"}
    for line in stdout.read_text().splitlines():
        try:
            stage = json.loads(line)
        except json.JSONDecodeError:
            continue
        if stage.pop("stage") == "faces":
            record.update(stage)
        elif returncode == 0:
            record.update(stage)
    write_json(raw / f"{model_id}.process.json", {
        "returncode": returncode, "wall_seconds": time.monotonic() - started,
        "observation_timeout_seconds": timeout,
    })
    return record


def summarize(records):
    complete = [row for row in records if row["coverage_complete"]]
    kinds = Counter()
    for row in records:
        kinds.update(row["face_kinds"])
    def label(value):
        return "unknown" if value is None else str(value).lower()
    all_kinds = Counter(label(row["all_face_kinds_buildable"]) for row in records)
    splines = Counter(label(row["needs_spline_faces"]) for row in records)
    return {
        "evaluated_model_ids": len(records), "complete_source_traversals": len(complete),
        "face_kind_counts": dict(sorted(kinds.items())),
        "face_kind_part_counts": {kind: sum(row["face_kinds"].get(kind, 0) > 0
                                            for row in records) for kind in KINDS},
        "all_face_kinds_buildable": dict(sorted(all_kinds.items())),
        "needs_spline_faces": dict(sorted(splines.items())),
        "all_face_kinds_buildable_share_of_complete": (
            sum(row["all_face_kinds_buildable"] is True for row in complete) / len(complete)
            if complete else None),
        "outcome_counts": dict(sorted(Counter(row["outcome"] for row in records).items())),
        "refusal_counts": dict(sorted(Counter(row["refusal_kind"] for row in records
                                               if row["outcome"] == "refused").items())),
        "fillet_count": None, "chamfer_count": None,
        "feature_count_unknown_model_ids": len(records),
    }


def table(path, records):
    columns = ["model_id", "coverage_complete", "source_parse", "referenced_faces",
               *KINDS, "solid_entities", "all_face_kinds_buildable", "needs_spline_faces",
               "unreferenced_faces", "missing_references", "unresolved_faces",
               "assembly_or_mapping_entities", "unknown_representation_items",
               "fillet_count", "chamfer_count",
               "feature_count_status", "outcome", "refusal_code", "refusal_kind", "bodies",
               "work_used", "work_limit"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, columns, lineterminator="\n")
        writer.writeheader()
        for record in records:
            row = {key: record.get(key) for key in columns}
            row.update({kind: record["face_kinds"].get(kind, 0) for kind in KINDS})
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", required=True, type=Path)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--count", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--workers", type=int, default=len(os.sched_getaffinity(0)))
    parser.add_argument("--timeout", type=int, default=60,
                        help="external observation seconds, not a proof allowance")
    parser.add_argument("--validate-workers", type=int, default=24,
                        help="repeat this many selected ids with 1, 2 and 4 workers")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    scratch, binary = args.scratch.resolve(), args.binary.resolve()
    if scratch.is_relative_to(repo):
        parser.error("scratch must be outside the source checkout")
    if min(args.count, args.workers, args.timeout) <= 0 or args.validate_workers < 0:
        parser.error("count, workers and timeout must be positive")
    archive = scratch / ARCHIVE
    if shutil.disk_usage(scratch).free < 30_000_000_000:
        raise SystemExit("less than 30 GB free; no extraction started")
    if archive.stat().st_size != ARCHIVE_BYTES or digest(archive, "md5") != ARCHIVE_MD5:
        raise SystemExit("archive size or published MD5 mismatch")
    provenance = {
        "schema": "catallax-r8-step-census/1",
        "production_revision": subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip(),
        "binary_sha256": digest(binary), "collector_sha256": digest(Path(__file__)),
        "observer_sha256": digest(repo / "crates/catallax-import/examples/r8_census.rs"),
        "syntax_parser_sha256": digest(repo / "crates/catallax-import/src/part21.rs"),
        "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
        "archive_sha256": digest(archive), "archive_md5": ARCHIVE_MD5,
        "archive_bytes": ARCHIVE_BYTES, "seed": args.seed, "workers": args.workers,
        "observation_timeout_seconds": args.timeout,
        "import_policy": "import_step_bodies with unchanged decider defaults and existing 20000000-unit application build allowance; no repair or tolerance changes",
        "licence_scope": "aggregate analysis under dataset MIT grant; per-model source exceptions not individually verified",
        "sample_unit": "distinct public ABC model ids; exact byte duplicates retained and counted",
    }
    with (scratch / "step-entries.txt").open("w") as stream:
        subprocess.run(["bsdtar", "-tf", str(archive)], stdout=stream, check=True)
    selected, available = choose((scratch / "step-entries.txt").read_text().splitlines(),
                                 args.count, args.seed)
    provenance["chunk_model_ids"] = available
    write_json(scratch / "selection-private.json", selected)
    (scratch / "extract-list.txt").write_text("".join(entry + "\n" for _, entry in selected))
    extracted = scratch / "selected-step"
    extracted.mkdir(exist_ok=False)
    raw = scratch / "raw"
    raw.mkdir(exist_ok=False)
    records = []
    try:
        subprocess.run(["bsdtar", "-xf", str(archive), "-C", str(extracted),
                        "-T", str(scratch / "extract-list.txt")], check=True)
        paths = [(model_id, extracted / entry) for model_id, entry in selected]
        if any(not path.is_file() or path.is_symlink() for _, path in paths):
            raise RuntimeError("missing or non-regular selected input")
        hashes = Counter(digest(path) for _, path in paths)
        provenance.update({"distinct_content_hashes": len(hashes),
                           "duplicate_content_groups": sum(n > 1 for n in hashes.values()),
                           "duplicate_model_ids_beyond_first": sum(n - 1 for n in hashes.values())})
        def run(item, folder=raw):
            model_id, path = item
            return observe(binary, path, model_id, folder, args.timeout)
        with ThreadPoolExecutor(max_workers=args.workers) as pool, \
                (scratch / "records-private.jsonl").open("w") as checkpoint:
            for record in pool.map(run, paths):
                records.append(record)
                checkpoint.write(json.dumps(record, sort_keys=True) + "\n")
                checkpoint.flush()
                if len(records) % 100 == 0:
                    print(f"observed {len(records)}/{len(paths)}", flush=True)
        baseline = {record["model_id"]: record for record in records}
        checks = []
        for workers in (1, 2, 4):
            check_paths = paths[:args.validate_workers]
            folder = scratch / f"workers-{workers}"
            folder.mkdir(exist_ok=False)
            with ThreadPoolExecutor(max_workers=workers) as pool:
                repeated = list(pool.map(lambda item: run(item, folder), check_paths))
            identical = all(row == baseline[row["model_id"]] for row in repeated)
            raw_identical = all((folder / f"{model_id}{suffix}").read_bytes()
                                == (raw / f"{model_id}{suffix}").read_bytes()
                                for model_id, _ in check_paths
                                for suffix in (".jsonl", ".stderr"))
            checks.append({"workers": workers, "models": len(repeated),
                           "identical": identical and raw_identical})
        provenance["worker_checks"] = checks
        provenance.update(summarize(records))
        assert len(records) == len(selected)
        assert sum(provenance["outcome_counts"].values()) == len(records)
        assert digest(binary) == provenance["binary_sha256"]
        assert all(check["identical"] for check in checks), "worker repeat mismatch"
        table(scratch / "R8-public-step-census.csv", records)
        write_json(scratch / "R8-public-step-census.json", provenance)
    finally:
        shutil.rmtree(extracted)
        archive.unlink()
        write_json(scratch / "cleanup.json", {
            "archive_removed": not archive.exists(),
            "extracted_models_removed": not extracted.exists(),
        })


if __name__ == "__main__":
    main()
