#!/usr/bin/env python3
"""Rerun the pinned 745 planar refusals through recorded import repair."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
from decimal import Decimal
import json
import os
from pathlib import Path
import random
import subprocess

from planar_self_intersection import cohort
from r8_step_census import digest, write_json


def distribution(values):
    values = sorted(Decimal(value) for value in values)
    if not values:
        return {"count": 0}
    # Nearest-rank percentiles retain the observed certified decimal bounds.
    return {"count": len(values), "min_mm": str(values[0]),
            **{f"p{p}_mm": str(values[max(0, (len(values) * p + 99) // 100 - 1)])
               for p in (50, 95, 99)}, "max_mm": str(values[-1])}


def public_row(row):
    faces = row.get("source_incidence", [])
    error = row.get("refusal", {})
    result = {key: row.get(key) for key in
              ("model_id", "source_sha256", "sampled", "outcome", "complete_coverage_admission",
               "bodies", "work_used", "work_limit")}
    result.update(refusal_code=error.get("code"), refusal_stage=error.get("stage"),
                  refusal_measurement=error.get("measurement"),
                  refusal_lower_mm=error.get("deviation_lower_mm"),
                  refusal_upper_mm=error.get("deviation_upper_mm"), repaired_faces=len(faces),
                  refusal_feature_size_squared_mm2=error.get("feature_size_squared_mm2"),
                  repair_methods=";".join(sorted({face["method"] for face in faces})))
    for field in ("source_discrepancy_upper_mm", "geometric_deviation_upper_mm",
                  "source_to_result_upper_mm", "vertex_movement_upper_mm"):
        result["max_" + field] = str(max(Decimal(face[field]) for face in faces)) if faces else None
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("models", "selection", "pair-details", "baseline", "binary", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--workers", type=int, default=len(os.sched_getaffinity(0)))
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    out = args.out.resolve()
    if out.is_relative_to(repo) or args.workers < 1:
        parser.error("evidence must stay outside the checkout; workers must be positive")
    paths = cohort(json.loads(args.pair_details.read_text()),
                   json.loads(args.selection.read_text()), args.models.resolve())
    previous = {row["model_id"]: row for row in
                map(json.loads, args.baseline.read_text().splitlines())}
    if len(paths) != 745 or set(paths) != set(previous):
        parser.error("expected the same pinned 745-model cohort")
    sampled = sorted(random.Random(20261111).sample(sorted(paths), 30))
    out.mkdir(exist_ok=False)
    raw = out / "raw"
    raw.mkdir()
    binary = args.binary.resolve()

    def run(model_id, folder):
        source_sha = digest(paths[model_id])
        if source_sha != previous[model_id]["source_sha256"]:
            raise ValueError("source hash differs from the original census")
        result = subprocess.run([str(binary), str(paths[model_id])], capture_output=True, check=False)
        stdout, stderr = result.stdout, result.stderr
        row = json.loads(stdout) if not result.returncode else {"outcome": "process_error"}
        (folder / (model_id + ".json")).write_bytes(stdout)
        (folder / (model_id + ".stderr")).write_bytes(stderr)
        if row.get("work_limit") not in (None, 20_000_000):
            raise ValueError("the import work allowance changed")
        if row.get("source_sha256") not in (None, source_sha):
            raise ValueError("observer source hash differs")
        admitted = row["outcome"] == "admitted"
        source = previous[model_id]
        complete = admitted and source["coverage_complete"] and source["all_face_kinds_buildable"]
        complete &= (row.get("bodies") == source["solid_entities"] and row.get("left_open") == 0
                     and row.get("topology_valid") is True)
        if admitted:
            for face in row["source_incidence"]:
                if Decimal(face["source_to_result_upper_mm"]) > Decimal("0.01"):
                    raise ValueError("admitted source bound exceeds the unchanged import cap")
        return {**row, "model_id": model_id, "source_sha256": source_sha,
                "sampled": model_id in sampled, "complete_coverage_admission": bool(complete)}

    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool, (out / "records-private.jsonl").open("w") as log:
        for row in pool.map(lambda i: run(i, raw), sorted(paths)):
            rows.append(row)
            log.write(json.dumps(row, sort_keys=True) + "\n")
            log.flush()
    original = {row["model_id"]: row for row in rows}
    # Retain the seeded sample, both constructions and a charged-work refusal.
    representatives = {}
    for row in rows:
        for face in row.get("source_incidence", []):
            representatives.setdefault(face["method"], row["model_id"])
        if row.get("refusal", {}).get("code") == "resource_limit":
            representatives.setdefault("resource_limit", row["model_id"])
    checked_ids = sorted(set(sampled) | set(representatives.values()))
    checks = []
    for workers in (1, 2, 4):
        folder = out / f"workers-{workers}"
        folder.mkdir()
        old_threads = os.environ.get("CATALLAX_THREADS")
        os.environ["CATALLAX_THREADS"] = str(workers)
        try:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                repeated = list(pool.map(lambda i: run(i, folder), checked_ids))
        finally:
            if old_threads is None:
                os.environ.pop("CATALLAX_THREADS", None)
            else:
                os.environ["CATALLAX_THREADS"] = old_threads
        identical = all(row == original[row["model_id"]] for row in repeated)
        identical &= all((folder / (i + suffix)).read_bytes() == (raw / (i + suffix)).read_bytes()
                         for i in checked_ids for suffix in (".json", ".stderr"))
        checks.append({"workers": workers, "kernel_threads": workers,
                       "models": len(checked_ids), "identical": identical})
    faces = [face for row in rows for face in row.get("source_incidence", [])]
    manifest = {
        "schema": "catallax-planar-incidence-repair/1", "models": len(rows),
        "seed": 20261111, "sampled_ids": sampled, "checked_ids": checked_ids, "worker_checks": checks,
        "outcomes": dict(Counter(row["outcome"] for row in rows)),
        "complete_coverage_admissions": sum(row["complete_coverage_admission"] for row in rows),
        "first_refusal_codes": dict(Counter(row["refusal"]["code"] for row in rows if "refusal" in row)),
        "planar_refusal_stages": dict(Counter(row["refusal"]["stage"] for row in rows
                                             if row.get("refusal", {}).get("stage"))),
        "nearest_feature_refusals": dict(Counter(
            "bound_exceeds_feature" if row["refusal"].get("feature_size_squared_mm2") is not None
            else "missing_complete_cover" for row in rows
            if row.get("refusal", {}).get("stage") == "nearest_feature")),
        "repair_methods_by_face": dict(Counter(face["method"] for face in faces)),
        "admitted_face_bounds": {field: distribution(face[field] for face in faces) for field in
                                 ("source_discrepancy_upper_mm", "geometric_deviation_upper_mm",
                                  "source_to_result_upper_mm", "vertex_movement_upper_mm")},
        "refusal_bound_scopes": dict(Counter(row["refusal"].get("measurement", "unavailable")
                                             for row in rows if "refusal" in row)),
        "import_work_limit": 20_000_000, "observation_timeout_seconds": None,
        "workers": args.workers, "kernel_threads": os.environ.get("CATALLAX_THREADS", "auto"),
        "base_revision": subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip(),
        "binary_sha256": digest(binary), "collector_sha256": digest(Path(__file__)),
        "observer_sha256": digest(repo / "crates/catallax-import/examples/r8_planar_repair.rs"),
        "repair_sources_sha256": {name: digest(repo / "crates/catallax-import/src" / name)
                                  for name in ("planar_incidence.rs", "planar_arithmetic.rs", "planar_deviation.rs",
                                               "incidence.rs", "incidence_arithmetic.rs", "incidence_repair.rs",
                                               "incidence_features.rs", "incidence_straight.rs", "deviation.rs",
                                               "repair.rs", "gap.rs", "intersect.rs", "wire.rs", "build.rs")},
        "cargo_lock_sha256": digest(repo / "Cargo.lock"),
        "baseline_sha256": digest(args.baseline), "selection_sha256": digest(args.selection),
        "source_archive_sha256": "52e6dd1b6fa38e3cd99af59b662370829129540030975919c3e256dce6ad1dbe",
        "evidence": "raw evidence retained privately",
    }
    table = [public_row(row) for row in rows]
    with (out / "R8-planar-incidence-repair.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)
    write_json(out / "R8-planar-incidence-repair.json", manifest)
    assert all(check["identical"] for check in checks), "worker comparison failed"
    print(json.dumps({key: manifest[key] for key in
                      ("models", "outcomes", "complete_coverage_admissions", "first_refusal_codes",
                       "planar_refusal_stages", "admitted_face_bounds", "worker_checks")}))


if __name__ == "__main__":
    main()
