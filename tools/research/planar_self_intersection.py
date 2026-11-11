#!/usr/bin/env python3
"""Audit the saved planar cohort and rerun its public import entry."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import json
import os
from pathlib import Path, PurePosixPath
import random
import subprocess

from r8_step_census import digest, observe, write_json


def cohort(details, selection, models):
    ids = sorted(i for i, row in details.items()
                 if row["intersection_path"] == "planar_polygon_pass")
    members = dict(selection)
    paths = {}
    for i in ids:
        member = PurePosixPath(members[i])
        if member.is_absolute() or ".." in member.parts:
            raise ValueError("unsafe archive member")
        path = models / member
        if not path.is_file() or path.is_symlink() or not path.resolve().is_relative_to(models.resolve()):
            raise ValueError("missing or unsafe model input")
        paths[i] = path
    return paths


def public_table(path, rows):
    fields = ["model_id", "source_sha256", "sampled", "classification", "shared_edges",
              "shared_vertices", "audit_status", "off_plane_vertices_a", "off_plane_vertices_b",
              "outcome", "refusal_code", "refusal_kind", "bodies", "work_used", "work_limit"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fields})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("models", "selection", "saved-refusals", "pair-details", "audit-binary", "import-binary", "out"):
        p.add_argument("--" + name, required=True, type=Path)
    p.add_argument("--workers", type=int, default=len(os.sched_getaffinity(0)))
    p.add_argument("--seed", type=int, default=20261111)
    p.add_argument("--sample", type=int, default=30)
    args = p.parse_args()
    repo = Path(__file__).resolve().parents[2]
    out = args.out.resolve()
    if out.is_relative_to(repo):
        p.error("evidence must remain outside the checkout")
    paths = cohort(json.loads(args.pair_details.read_text()),
                   json.loads(args.selection.read_text()), args.models.resolve())
    if len(paths) != 745 or not 30 <= args.sample <= len(paths) or args.workers < 1:
        p.error("expected the pinned 745-model cohort, at least 30 samples and positive workers")
    sampled = sorted(random.Random(args.seed).sample(sorted(paths), args.sample))
    out.mkdir(exist_ok=False)
    raw = out / "raw"
    raw.mkdir()
    audit_binary, import_binary = args.audit_binary.resolve(), args.import_binary.resolve()

    def run(i, folder):
        try:
            completed = subprocess.run(
                [str(audit_binary), str(paths[i]), str(args.saved_refusals / (i + ".stderr"))],
                capture_output=True, timeout=60, check=False)
            stdout, stderr = completed.stdout, completed.stderr
            status = "process_error" if completed.returncode else "observed"
        except subprocess.TimeoutExpired as error:
            stdout, stderr = error.stdout or b"", error.stderr or b""
            status = "observation_timeout"
        (folder / (i + ".audit.json")).write_bytes(stdout)
        (folder / (i + ".audit.stderr")).write_bytes(stderr)
        audit = {"classification": "not_verified", "off_plane_vertices": [None, None],
                 "shared_edges": None, "shared_vertices": None, "geometry_relation_verified": False}
        if status == "observed":
            try:
                audit = json.loads(stdout)
            except ValueError:
                status = "invalid_output"
        audit["audit_status"] = status
        record = observe(import_binary, paths[i], i, folder, 60)
        if record.get("work_limit") is not None and record["work_limit"] != 20_000_000:
            raise RuntimeError("import allowance changed")
        return {**record, **audit, "source_sha256": digest(paths[i]), "sampled": i in sampled,
                "off_plane_vertices_a": audit["off_plane_vertices"][0],
                "off_plane_vertices_b": audit["off_plane_vertices"][1]}

    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool, (out / "records-private.jsonl").open("w") as log:
        for row in pool.map(lambda i: run(i, raw), paths):
            rows.append(row)
            log.write(json.dumps(row, sort_keys=True) + "\n")
            log.flush()
    baseline = {row["model_id"]: row for row in rows}
    checks = []
    original_threads = os.environ.get("CATALLAX_THREADS")
    try:
        for workers in (1, 2, 4):
            os.environ["CATALLAX_THREADS"] = str(workers)
            folder = out / f"workers-{workers}"
            folder.mkdir()
            with ThreadPoolExecutor(max_workers=workers) as pool:
                repeated = list(pool.map(lambda i: run(i, folder), sampled))
            equal = all(row == baseline[row["model_id"]] for row in repeated)
            equal &= all((folder / (i + suffix)).read_bytes() == (raw / (i + suffix)).read_bytes()
                         for i in sampled for suffix in (".audit.json", ".audit.stderr", ".jsonl", ".stderr"))
            checks.append({"workers": workers, "kernel_threads": workers,
                           "models": len(repeated), "identical": equal})
    finally:
        if original_threads is None:
            os.environ.pop("CATALLAX_THREADS", None)
        else:
            os.environ["CATALLAX_THREADS"] = original_threads
    manifest = {
        "schema": "catallax-planar-self-intersection/1", "models": len(rows),
        "seed": args.seed, "sampled_ids": sampled, "workers": args.workers,
        "kernel_threads": original_threads or "auto",
        "base_revision": subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip(),
        "audit_binary_sha256": digest(audit_binary), "import_binary_sha256": digest(import_binary),
        "collector_sha256": digest(Path(__file__)),
        "pair_details_sha256": digest(args.pair_details),
        "selection_sha256": digest(args.selection),
        "source_archive_sha256": "52e6dd1b6fa38e3cd99af59b662370829129540030975919c3e256dce6ad1dbe",
        "observer_sha256": digest(repo / "crates/catallax-import/examples/r8_planar_contact.rs"),
        "intersection_source_sha256": digest(repo / "crates/catallax-import/src/intersect.rs"),
        "classifications": dict(Counter(row["classification"] for row in rows)),
        "audit_statuses": dict(Counter(row["audit_status"] for row in rows)),
        "outcomes": dict(Counter(row["outcome"] for row in rows)),
        "first_refusals": dict(Counter(row.get("refusal_kind") for row in rows if row["outcome"] == "refused")),
        "worker_checks": checks, "import_work_limit": 20_000_000,
        "observation_timeout_seconds": 60, "evidence": "raw evidence retained privately",
    }
    public_table(out / "R8-planar-self-intersection.csv", rows)
    write_json(out / "R8-planar-self-intersection.json", manifest)
    assert all(row["identical"] for row in checks), "worker comparison failed"
    print(json.dumps({key: manifest[key] for key in ("models", "classifications", "outcomes", "first_refusals", "worker_checks")}))


if __name__ == "__main__":
    main()
