#!/usr/bin/env python3
"""Tabulate saved first refusals and audit their face references without importing."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess

from r8_step_census import ARCHIVE, ARCHIVE_BYTES, ARCHIVE_MD5, digest, write_json


def counterfactual(count):
    return {"status": "unmeasured_first_refusal_only", "no_other_refusal_models": None,
            "confirmed_models": 0, "unknown_models": count, "possible_range": [0, count]}


def summarize(records, errors, details):
    unsupported, pairs = Counter(), Counter()
    for record in records:
        if record["outcome"] != "refused":
            continue
        model_id = record["model_id"]
        error = errors[model_id]
        assert error["code"] == record["refusal_code"]
        if error["code"] == "unsupported":
            unsupported[error["type_name"]] += 1
        elif error["code"] == "self_intersecting":
            detail = details.get(model_id)
            if detail and detail.get("source_reference_audit") is True:
                assert detail["import_stage"] == "build_solid.check_solid_before_topology"
                assert detail["intersection_path"] in (
                    "planar_polygon_pass", "curved_or_circular_edge_pass")
                assert len(detail["face_kinds"]) == 2
                assert set(detail["face_kinds"]) <= {"plane", "cylinder", "cone", "sphere", "torus"}
                assert set(detail["face_kinds"]) <= {kind for kind, count in record["face_kinds"].items() if count}
                key = (detail["intersection_path"], tuple(sorted(detail["face_kinds"])))
            else:
                key = ("unresolved_metadata", ("unknown", "unknown"))
            pairs[key] += 1
    return {
        "unsupported_models": sum(unsupported.values()),
        "unsupported": [{"reported_type": name, "entity_type": name.split(" (", 1)[0],
                         "models": count, "counterfactual": counterfactual(count)}
                        for name, count in sorted(unsupported.items(), key=lambda item: (-item[1], item[0]))],
        "self_intersection_models": sum(pairs.values()),
        "self_intersections": [{"import_stage": "build_solid.check_solid_before_topology",
                                "intersection_path": path, "face_kinds": list(kinds),
                                "models": count, "counterfactual": counterfactual(count)}
                               for (path, kinds), count in sorted(pairs.items(), key=lambda item: (-item[1], item[0]))],
        "counterfactual_status": "unmeasured_first_refusal_only",
        "counterfactual_unknown_models": sum(unsupported.values()) + sum(pairs.values()),
    }


def render(summary):
    lines = ["## Refusal triage from retained diagnostics", "",
             "Unsupported rows use the exact STEP type and diagnostic variant saved by the importer. "
             "Self-intersection rows use the actual failed face references, not the model's aggregate face kinds.", "",
             "No candidate fix was evaluated. The records stop at the first refusal, so later refusals are unobserved. "
             "`unknown [0, N]` is a bound on models with no other refusal if that cause were fixed, not a predicted gain. "
             "Zero models are confirmed to have no other refusal in every row. Even a fully identified pair does not establish isolation.", "",
             "| Exact unsupported STEP type / diagnostic variant | Models | No other refusal after fixing this cause |",
             "|---|---:|---|"]
    for row in summary["unsupported"]:
        n = row["models"]
        lines.append(f"| `{row['reported_type']}` | {n:,} | unknown [0, {n:,}] |")
    lines += ["", "Every self-intersection below was returned by `build_solid`'s pre-topology "
              "`check_solid` stage. The pass is recovered from the recorded pair's surface and boundary-curve references. "
              "Two planes can still use the curved pass when a boundary is circular or elliptical.", "",
              "| Import substage | Failed face kinds | Models | No other refusal after fixing this cause |",
              "|---|---|---:|---|"]
    for row in summary["self_intersections"]:
        n = row["models"]
        lines.append(f"| `{row['intersection_path']}` | {' / '.join(row['face_kinds'])} | {n:,} | unknown [0, {n:,}] |")
    lines += ["", "The read-only audit verifies that both saved face ids belong to the saved solid, follows "
              "`SURFACE_CURVE`/`SEAM_CURVE` aliases, and reconstructs only the planar/curved dispatch predicate. "
              "It does not run, suppress or relax an intersection check, repeat a proof allowance, or construct an accepted body.", "",
              "The raw records cannot establish counterfactual admissions: an unsupported carrier must be implemented "
              "before its later checks can run, and a reported intersection must be resolved before the rest of that "
              "solid and later solids can be certified. A diagnostic bypass would not demonstrate that a model passes. "
              "The existing 982-placement cohort illustrates this limit: removing its first refusal produced no admissions.", "",
              "The [triage manifest](R8-public-step-refusal-triage.json) contains the counts, explicit unknowns and audit identities. "
              "No additional import answers were produced or changed.", ""]
    return "\n".join(lines)


def audit(binary, model_id, path, saved_error, raw, timeout):
    with (raw / f"{model_id}.json").open("wb") as out, (raw / f"{model_id}.stderr").open("wb") as err:
        try:
            completed = subprocess.run([str(binary), "--refusal-detail", str(path), str(saved_error)],
                                       stdout=out, stderr=err, timeout=timeout, check=False)
            reason = "process_error" if completed.returncode else None
        except subprocess.TimeoutExpired:
            reason = "observation_timeout"
    if reason is None:
        try:
            return json.loads((raw / f"{model_id}.json").read_text())
        except ValueError:
            reason = "invalid_audit_output"
    return {"source_reference_audit": False, "reason": reason}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--scratch", required=True, type=Path)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--render-only", action="store_true")
    parser.add_argument("--report", type=Path)
    parser.add_argument("--workers", type=int, default=len(os.sched_getaffinity(0)))
    parser.add_argument("--timeout", type=int, default=60)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    evidence, scratch = args.evidence.resolve(), args.scratch.resolve()
    if scratch.is_relative_to(repo):
        parser.error("scratch must be outside the source checkout")
    records = [json.loads(line) for line in (evidence / "records-private.jsonl").read_text().splitlines()]
    errors = {row["model_id"]: json.loads((evidence / "raw" / (row["model_id"] + ".stderr")).read_text())
              for row in records if row["outcome"] == "refused"}
    original = json.loads((evidence / "R8-public-step-census.json").read_text())
    if args.render_only:
        config = json.loads((scratch / "audit-config.json").read_text())
        details = json.loads((scratch / "pair-details.json").read_text())
    else:
        if args.binary is None:
            parser.error("collection requires --binary")
        binary = args.binary.resolve()
        archive = scratch / ARCHIVE
        if shutil.disk_usage(scratch).free < 30_000_000_000:
            raise SystemExit("less than 30 GB free")
        if (archive.stat().st_size != ARCHIVE_BYTES or digest(archive, "md5") != ARCHIVE_MD5
                or digest(archive) != original["archive_sha256"]):
            raise SystemExit("archive identity mismatch")
        wanted = {row["model_id"] for row in records if row.get("refusal_code") == "self_intersecting"}
        selected = [(model_id, entry) for model_id, entry in
                    json.loads((evidence / "selection-private.json").read_text()) if model_id in wanted]
        assert len(selected) == len(wanted)
        for _, entry in selected:
            path = PurePosixPath(entry)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError("unsafe source member")
        config = {"schema": "catallax-r8-refusal-triage/1",
                  "original_production_revision": original["production_revision"],
                  "audit_source_revision": subprocess.check_output(
                      ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip(),
                  "binary_sha256": digest(binary), "archive_sha256": original["archive_sha256"],
                  "original_records_sha256": digest(evidence / "records-private.jsonl"),
                  "observer_sha256": digest(repo / "crates/catallax-import/examples/r8_census.rs"),
                  "pair_auditor_sha256": digest(repo / "crates/catallax-import/examples/r8_census/refusal_detail.rs"),
                  "workers": args.workers, "observation_timeout_seconds": args.timeout,
                  "import_attempts_repeated": 0, "geometry_checks_bypassed": 0}
        write_json(scratch / "audit-config.json", config)
        (scratch / "extract-list.txt").write_text("".join(entry + "\n" for _, entry in selected))
        extracted, raw = scratch / "selected-step", scratch / "pair-audits"
        extracted.mkdir(exist_ok=False)
        raw.mkdir(exist_ok=False)
        details = {}
        try:
            subprocess.run(["bsdtar", "-xf", str(archive), "-C", str(extracted),
                            "-T", str(scratch / "extract-list.txt")], check=True)
            def run(item):
                model_id, entry = item
                path = extracted / entry
                if not path.is_file() or path.is_symlink():
                    raise ValueError("non-regular source member")
                return model_id, audit(binary, model_id, path, evidence / "raw" / (model_id + ".stderr"), raw, args.timeout)
            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                for model_id, detail in pool.map(run, selected):
                    details[model_id] = detail
                    write_json(scratch / "pair-details.json", details)
        finally:
            shutil.rmtree(extracted)
            archive.unlink()
            write_json(scratch / "cleanup.json", {"archive_removed": not archive.exists(),
                       "extracted_models_removed": not extracted.exists()})
    summary = {**config, **summarize(records, errors, details), "generator_sha256": digest(Path(__file__))}
    assert summary["unsupported_models"] == 1314 and summary["self_intersection_models"] == 822
    write_json(scratch / "R8-public-step-refusal-triage.json", summary)
    section = render(summary)
    (scratch / "report-section.md").write_text(section)
    if args.report:
        start, end = "<!-- r8-refusal-triage:start -->", "<!-- r8-refusal-triage:end -->"
        report = args.report.read_text()
        block = start + "\n" + section + end
        if start in report:
            before, remainder = report.split(start, 1)
            _, after = remainder.split(end, 1)
            report = before + block + after
        else:
            report = report.replace("## Source coverage and face kinds", block + "\n\n## Source coverage and face kinds", 1)
        args.report.write_text(report)


if __name__ == "__main__":
    main()
