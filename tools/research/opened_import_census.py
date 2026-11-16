"""Rerun all 3,000 frozen sources through the tier consumer and actual STEP export.

Retained source surfaces, admitted solids and independently valid exports are
separate measures. Captures stay outside the checkout. A partial run has no
complete percentage or acceptance credit. OCCT may process/heal inside its
builder; these observations establish neither unhealed validity nor fidelity.
"""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import time

from export_validity_census import (
    ADMISSIONS_SHA,
    ARCHIVE_SHA,
    HASHES_SHA,
    RAW_HELPER_SHA,
    SELECTION_SHA,
    WORK,
    digest,
    git,
    public_build,
    save,
)

ROOT = Path(__file__).resolve().parents[2]
OBSERVER = ROOT / "crates/catallax-import/examples/opened_import.rs"
BASELINE = ROOT / "docs/research/step-census-periodic/results.csv"
VALIDITY = ROOT / "docs/research/step-census-circular-pcurves/results.csv"
RETAINED = ROOT / "docs/research/step-census/retained-controls.json"
CONE_REFERENCE = ROOT / "docs/research/step-census-opened/cone-reference.json"


def build(args):
    assert not args.record.resolve().is_relative_to(ROOT)
    assert not git("status", "--porcelain")
    args.record.parent.mkdir(parents=True, exist_ok=True)
    cargo_json = args.record.with_suffix(".cargo.jsonl")
    command = [
        "cargo",
        "build",
        "--release",
        "-p",
        "catallax-import",
        "--example",
        "opened_import",
        "--message-format=json-render-diagnostics",
    ]
    with cargo_json.open("w") as stream:
        subprocess.run(command, cwd=ROOT, stdout=stream, check=True, timeout=3600)
    artifacts = [
        a
        for line in cargo_json.read_text().splitlines()
        if (a := json.loads(line))["reason"] == "compiler-artifact"
        and a["target"]["name"] == "opened_import"
        and a.get("executable")
    ]
    assert len(artifacts) == 1
    artifact = artifacts[0]
    assert artifact["target"]["kind"] == ["example"] and not artifact["profile"]["test"]
    revision = git("rev-parse", "HEAD")
    record = dict(
        format="catallax-opened-observer-build/1",
        source_revision=revision,
        kernel_revision=revision,
        git_tree=git("rev-parse", "HEAD^{tree}"),
        artifact_path=artifact["executable"],
        cargo_artifact=artifact,
        binary_sha256=digest(Path(artifact["executable"])),
        cargo_json_sha256=digest(cargo_json),
        cargo_lock_sha256=digest(ROOT / "Cargo.lock"),
        observer_sha256=digest(OBSERVER),
        collector_sha256=digest(Path(__file__)),
        raw_helper_sha256=RAW_HELPER_SHA,
        build_command=command,
        rustc=subprocess.check_output(["rustc", "--version"], text=True).strip(),
    )
    assert not git("status", "--porcelain")
    save(args.record, record)
    print(
        json.dumps(dict(source_revision=revision, binary_sha256=record["binary_sha256"])),
        flush=True,
    )


def raw_surface(path):
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.Interface import Interface_Static
    from OCP.STEPControl import STEPControl_Reader
    from OCP.StepData import StepData_Factors
    from OCP.StepShape import StepShape_ShellBasedSurfaceModel
    from OCP.StepToTopoDS import StepToTopoDS_Builder, StepToTopoDS_NMTool
    from OCP.Transfer import Transfer_TransientProcess

    reader = STEPControl_Reader()
    settings = {"read.stdsameparameter.mode": 0, "read.maxprecision.mode": 0}
    previous = {k: Interface_Static.IVal_s(k) for k in settings}
    try:
        for key, value in settings.items():
            assert Interface_Static.SetIVal_s(key, value)
        assert reader.ReadFile(str(path)) == IFSelect_RetDone
        model = reader.StepModel()
        roots = [
            model.Value(i)
            for i in range(1, model.NbEntities() + 1)
            if isinstance(model.Value(i), StepShape_ShellBasedSurfaceModel)
        ]
        assert len(roots) == 1
        process = Transfer_TransientProcess()
        process.SetModel(model)
        factors = StepData_Factors()
        factors.InitializeFactors(1.0, 1.0, 1.0)
        builder = StepToTopoDS_Builder()
        builder.SetPrecision(1e-7)
        builder.SetMaxTol(1e-7)
        builder.Init(roots[0], process, StepToTopoDS_NMTool(), factors)
        assert builder.IsDone()
        return builder.Value()
    finally:
        for key, value in previous.items():
            Interface_Static.SetIVal_s(key, value)


def raw(args):
    helper = ROOT / "tools/step-regressions/check_numerical_step.py"
    assert digest(helper) == RAW_HELPER_SHA
    sys.path.insert(0, str(helper.parent))
    import OCP
    from OCP.BRepCheck import BRepCheck_Analyzer
    from check_numerical_step import raw_solid, topology

    assert OCP.__version__ == "8.0.1.1"
    response = json.loads((args.capture / "response.json").read_text())
    assert response["source_sha256"] == args.source_sha256
    results = []
    for entry in response["files"]:
        assert entry["file"] == f"body-{entry['body']}.step"
        path = args.capture / entry["file"]
        assert digest(path) == entry["sha256"] and path.stat().st_size == entry["bytes"]
        before = path.read_bytes()
        try:
            assert b"SI_UNIT(.MILLI.,.METRE.)" in b"".join(before.split())
            shape = raw_solid(path) if entry["kind"] == "solid" else raw_surface(path)
            counts = topology(shape)
            assert counts["solids"] == (1 if entry["kind"] == "solid" else 0)
            outcome = (
                "valid_after_builder"
                if BRepCheck_Analyzer(shape).IsValid()
                else "invalid_after_builder"
            )
        except Exception as error:
            outcome, counts = "reader_refused", None
            print(type(error).__name__, str(error), file=sys.stderr)
        assert path.read_bytes() == before
        results.append(
            dict(
                body=entry["body"],
                kind=entry["kind"],
                sha256=entry["sha256"],
                outcome=outcome,
                topology=counts,
            )
        )
    save(
        args.out,
        dict(
            source_sha256=args.source_sha256,
            ocp_version=OCP.__version__,
            files=results,
            unhealed_validity_established=False,
            complete_source_body_bound_established=False,
        ),
    )


def cause(refusal):
    if not refusal:
        return ""
    code = refusal.get("code", "unknown")
    return code + (":" + refusal.get("type_name", "unknown") if code == "unsupported" else "")


def opening(row, checkpoint):
    bodies = checkpoint["bodies"]
    tiers = Counter(b["tier"] for b in bodies)
    solids = tiers["valid_solid"] + tiers["repaired_solid"]
    row.update(
        opened=bool(bodies),
        outcome=(
            "refused"
            if not bodies
            else (
                "full"
                if solids == len(bodies) and row["source_coverage_complete"]
                else "partial" if solids else "surface"
            )
        ),
        observed_body_count=len(bodies),
        valid_solid_body_count=tiers["valid_solid"],
        repaired_solid_body_count=tiers["repaired_solid"],
        surface_body_count=tiers["surface_body"],
        open_work_used=checkpoint["work_used"],
        open_work_exhausted=checkpoint["work_exhausted"],
        read_work_used=checkpoint["read_work_used"],
        file_admission_work_used=checkpoint["file_admission_work_used"],
        first_refusal=cause(checkpoint["first_refusal"]),
    )
    row["proofs"] = dict(Counter(b["proof"] for b in bodies))
    if not bodies:
        row["first_refusal"] = "no_supported_bodies"
    after = []
    if not row["source_coverage_complete"]:
        after.append("source_coverage_incomplete")
    for body in bodies:
        if body["tier"] != "surface_body":
            continue
        if body["healing_refusal"]:
            after.append(cause(body["healing_refusal"]["refusal"]))
        elif body["refusal"]:
            after.append(cause(body["refusal"]))
        elif body["proof"] == "retained_source_surfaces" and row["first_refusal"]:
            after.append(row["first_refusal"])
        after.extend(i["code"] for i in body["healing_issues"] or [])
    row["after_fix_causes"] = sorted(set(after))
    row["unresolved_checks"] = sorted(
        {i["code"] for b in bodies for i in b["healing_issues"] or []}
    )


def observe(binary, source, baseline, folder, retain):
    key, expected = baseline["model_id"], baseline["source_sha256"]
    row = dict(
        model_id=key,
        source_sha256=expected,
        outcome="process_error",
        opened=False,
        baseline_outcome=baseline["baseline_outcome"],
        baseline_first_refusal=baseline["baseline_first_refusal"],
        previous_outcome=baseline["outcome"],
        previous_after_fix_cause=baseline["after_fix_cause"],
        source_coverage_complete=baseline["coverage_complete"] == "True",
    )
    assert digest(source) == expected
    folder.mkdir(parents=True, exist_ok=True)
    folder = folder / f"attempt-{len(list(folder.glob('attempt-*')))}"
    folder.mkdir()
    capture = folder / "capture"
    start = time.monotonic()
    command = [str(binary), str(source), expected, str(capture)] + (
        ["--retain-controls"] if retain else []
    )
    with (folder / "producer.stdout").open("wb") as out, (folder / "producer.stderr").open(
        "wb"
    ) as err:
        try:
            process = subprocess.run(
                command,
                stdout=out,
                stderr=err,
                timeout=60,
                check=False,
                env={**os.environ, "CATALLAX_THREADS": "1"},
            )
        except subprocess.TimeoutExpired:
            checkpoints = []
            for line in (folder / "producer.stdout").read_text().splitlines():
                try:
                    value = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if value.get("stage") == "open" and value.get("outcome") == "opened":
                    checkpoints.append(value)
            if checkpoints:
                assert len(checkpoints) == 1 and checkpoints[0]["source_sha256"] == expected
                opening(row, checkpoints[0])
            row.update(timeout_stage="export" if checkpoints else "open", export_validity="timeout")
            if not checkpoints:
                row["outcome"] = "timeout"
            return row
        finally:
            assert digest(source) == expected
    if process.returncode or not (capture / "response.json").is_file():
        return row | dict(cause="producer_process_or_missing_response")
    response = json.loads((capture / "response.json").read_text())
    assert response["source_sha256"] == expected and response["work_limit"] == WORK
    if response["stage"] == "open":
        assert response["outcome"] == "refused"
        return row | dict(
            outcome="refused",
            first_refusal=cause(response["refusal"]),
            after_fix_causes=[cause(response["refusal"])],
            open_work_used=response["work_used"],
            open_work_exhausted=response["work_exhausted"],
            export_validity="not_opened",
        )
    opening(row, response["open"])
    row.update(
        export_work_used=response["work_used"],
        exported_body_count=len(response["files"]),
        exported_step_hashes=[f["sha256"] for f in response["files"]],
    )
    if response["outcome"] == "refused":
        row["export_refusal"] = cause(response["refusal"]) or "resource_limit"
    if not response["files"]:
        row["export_validity"] = (
            "source_surfaces_without_checked_brep"
            if not row.get("export_refusal")
            else "export_refused"
        )
        return row
    remaining = 60 - (time.monotonic() - start)
    if remaining <= 0:
        return row | dict(export_validity="timeout", timeout_stage="raw_reader")
    with (folder / "raw.stdout").open("wb") as out, (folder / "raw.stderr").open("wb") as err:
        try:
            process = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(Path(__file__).resolve()),
                    "raw",
                    "--capture",
                    str(capture),
                    "--source-sha256",
                    expected,
                    "--out",
                    str(folder / "raw.json"),
                ],
                stdout=out,
                stderr=err,
                timeout=remaining,
                check=False,
            )
        except subprocess.TimeoutExpired:
            return row | dict(export_validity="timeout", timeout_stage="raw_reader")
    if process.returncode or not (folder / "raw.json").is_file():
        return row | dict(export_validity="process_error", cause="raw_process_or_missing_response")
    checked = json.loads((folder / "raw.json").read_text())
    assert checked["source_sha256"] == expected and checked["ocp_version"] == "8.0.1.1"
    assert [f["sha256"] for f in checked["files"]] == row["exported_step_hashes"]
    counts = Counter((f["kind"], f["outcome"]) for f in checked["files"])
    row.update(
        valid_exported_solid_count=counts["solid", "valid_after_builder"],
        valid_exported_surface_count=counts["surface", "valid_after_builder"],
        invalid_exported_body_count=sum(
            v for (k, o), v in counts.items() if o == "invalid_after_builder"
        ),
        reader_refused_body_count=sum(v for (k, o), v in counts.items() if o == "reader_refused"),
    )
    row["all_solids_export_valid"] = (
        row["outcome"] == "full" and row["valid_exported_solid_count"] == row["observed_body_count"]
    )
    row["export_validity"] = (
        "export_refused"
        if row.get("export_refusal")
        else (
            "reader_refused"
            if row["reader_refused_body_count"]
            else (
                "invalid_after_builder"
                if row["invalid_exported_body_count"]
                else (
                    "all_solids_valid_after_builder"
                    if row["all_solids_export_valid"]
                    else "valid_exported_subset"
                )
            )
        )
    )
    return row


def collect(args):
    assert args.workers == 2 and not args.out.resolve().is_relative_to(ROOT)
    assert digest(BASELINE) == ADMISSIONS_SHA and digest(args.source_hashes) == HASHES_SHA
    assert digest(args.selection) == SELECTION_SHA
    frozen_manifest = json.loads(BASELINE.with_suffix(".json").read_text())
    assert (
        frozen_manifest["complete"]
        and frozen_manifest["evaluated"] == frozen_manifest["denominator"] == 3000
    )
    assert frozen_manifest["seed"] == 20261111 and frozen_manifest["archive_sha256"] == ARCHIVE_SHA
    assert (
        frozen_manifest["selection_sha256"] == SELECTION_SHA
        and frozen_manifest["source_hash_manifest_sha256"] == HASHES_SHA
    )
    record = json.loads(args.build_record.read_text())
    assert record["format"] == "catallax-opened-observer-build/1"
    assert record["source_revision"] == git("rev-parse", "HEAD") and not git(
        "status", "--porcelain"
    )
    assert record["observer_sha256"] == digest(OBSERVER) and record["collector_sha256"] == digest(
        Path(__file__)
    )
    binary = Path(record["artifact_path"])
    assert digest(binary) == record["binary_sha256"]
    csv.field_size_limit(64 * 1024 * 1024)
    with BASELINE.open(newline="") as stream:
        baseline = {r["model_id"]: r for r in csv.DictReader(stream)}
    with VALIDITY.open(newline="") as stream:
        validity = {r["model_id"]: r for r in csv.DictReader(stream)}
    hashes = json.loads(args.source_hashes.read_text())
    selection = dict(json.loads(args.selection.read_text()))
    assert len(hashes) == len(selection) == len(baseline) == 3000 and set(hashes) == set(
        selection
    ) == set(baseline)
    assert all(
        baseline[k]["source_sha256"] == hashes[k] == validity[k]["source_sha256"] for k in hashes
    )
    retained = {r["model_id"] for r in json.loads(RETAINED.read_text())["rows"]}
    cone_reference = json.loads(CONE_REFERENCE.read_text())
    assert cone_reference["source_revision"] == "88bb7851a00f0a8531fd85618a2cf083fcda510c"
    assert (
        len(cone_reference["results"]) == 5
        and sum(r["expected_body_count"] for r in cone_reference["results"]) == 7
    )
    assert all(
        r["outcome"] == "raw_valid" and r["source_sha256"] == hashes[r["model_id"]]
        for r in cone_reference["results"]
    )
    files = {}
    for key, member in selection.items():
        path = PurePosixPath(member)
        assert not path.is_absolute() and ".." not in path.parts
        source = (args.sources / str(path)).resolve()
        assert source.is_relative_to(args.sources.resolve()) and digest(source) == hashes[key]
        files[key] = source
    identity = dict(
        format="catallax-opened-import-census/1",
        seed=20261111,
        denominator=3000,
        source_revision=record["source_revision"],
        binary_sha256=record["binary_sha256"],
        build_record=public_build(record),
        archive_sha256=ARCHIVE_SHA,
        source_hash_manifest_sha256=HASHES_SHA,
        selection_sha256=SELECTION_SHA,
        admission_csv_sha256=ADMISSIONS_SHA,
        raw_helper_sha256=RAW_HELPER_SHA,
        open_work_limit=WORK,
        cone_reference_sha256=digest(CONE_REFERENCE),
        aggregate_export_work_limit=WORK,
        observation_timeout_seconds=60,
        workers=2,
        per_import_workers=1,
        ocp_version="8.0.1.1",
        unhealed_validity_established=False,
        complete_source_body_bound_established=False,
    )
    args.out.mkdir(parents=True, exist_ok=args.resume)
    if args.resume:
        assert (
            json.loads((args.out / "identity.json").read_text()) == identity
        ), "resume producer changed"
    save(args.out / "identity.json", identity)
    rows_dir = args.out / "rows"
    rows_dir.mkdir(exist_ok=True)
    rows = (
        {p.stem: json.loads(p.read_text()) for p in rows_dir.glob("*.json")} if args.resume else {}
    )
    assert set(rows).issubset(hashes)
    observations = list(rows.values())
    save(args.out / "results.json", identity | dict(complete=False, evaluated=len(observations)))
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending = {
            pool.submit(
                observe, binary, files[k], baseline[k], args.out / "parts" / k, k in retained
            ): k
            for k in sorted(hashes)
            if k not in rows
        }
        for future in as_completed(pending):
            key = pending[future]
            row = future.result()
            save(rows_dir / f"{key}.json", row)
            observations.append(row)
            if len(observations) % 50 == 0 or len(observations) == 3000:
                save(
                    args.out / "results.json",
                    identity
                    | dict(
                        complete=False,
                        evaluated=len(observations),
                        outcomes=dict(Counter(r["outcome"] for r in observations)),
                    ),
                )
                print(
                    json.dumps(
                        dict(
                            evaluated=len(observations),
                            outcomes=dict(Counter(r["outcome"] for r in observations)),
                        )
                    ),
                    flush=True,
                )
    assert len(observations) == 3000 and {r["model_id"] for r in observations} == set(hashes)
    assert (
        all(digest(files[k]) == hashes[k] for k in hashes)
        and digest(binary) == record["binary_sha256"]
    )
    observations.sort(key=lambda r: r["model_id"])
    table = args.out / "results.csv"
    fields = sorted(set().union(*(r.keys() for r in observations)))
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            {
                k: json.dumps(v, separators=(",", ":")) if isinstance(v, (list, dict)) else v
                for k, v in r.items()
            }
            for r in observations
        )
    by_id = {r["model_id"]: r for r in observations}
    expected_full = {k: int(r["bodies"]) for k, r in baseline.items() if r["outcome"] == "full"}
    expected_full.update(
        {r["model_id"]: r["expected_body_count"] for r in cone_reference["results"]}
    )
    expected_valid = {k for k, r in validity.items() if r["validity_outcome"] == "raw_valid"}
    expected_valid.update(r["model_id"] for r in cone_reference["results"])
    assert len(expected_full) == 405 and len(expected_valid) == 374
    tier_losses = [
        k
        for k, count in expected_full.items()
        if by_id[k]["outcome"] != "full" or by_id[k].get("observed_body_count") != count
    ]
    export_losses = [k for k in expected_valid if not by_id[k].get("all_solids_export_valid")]
    outcomes = Counter(r["outcome"] for r in observations)
    exports = Counter(r.get("export_validity", "not_reached") for r in observations)
    opened = sum(r["opened"] for r in observations)
    qualified = sum(bool(r.get("all_solids_export_valid")) for r in observations)
    preserved_9338 = bool(
        by_id["00009338"].get("all_solids_export_valid")
        and by_id["00009338"].get("proofs", {}).get("numerical")
        and not by_id["00009338"].get("first_refusal")
    )
    for key in retained:
        captures = list((args.out / "parts" / key).glob("attempt-*/capture/source.step"))
        assert captures and digest(captures[-1]) == hashes[key]
        assert (captures[-1].parent / "controls.json").is_file()
    summary = identity | dict(
        complete=True,
        evaluated=3000,
        outcomes=dict(outcomes),
        export_outcomes=dict(exports),
        opened=opened,
        opened_share_percent=f"{opened/30:.2f}",
        full_solid_files=outcomes["full"],
        full_solid_share_percent=f"{outcomes['full']/30:.2f}",
        all_solids_export_valid_files=qualified,
        all_solids_export_valid_share_percent=f"{qualified/30:.2f}",
        body_tiers={
            name: sum(r.get(name + "_body_count", 0) for r in observations)
            for name in ["valid_solid", "repaired_solid", "surface"]
        },
        tier_losses=tier_losses,
        frozen_valid_export_losses=export_losses,
        baseline_full_solid_files=405,
        baseline_valid_export_files=374,
        integrated_00009338_preserved=preserved_9338,
        acceptance_passed=not tier_losses
        and not export_losses
        and not outcomes["process_error"]
        and not exports["process_error"]
        and preserved_9338,
        retained_controls=[by_id[k] for k in sorted(retained)],
        integrated_00009338=by_id["00009338"],
        refusal_causes=dict(
            Counter(r.get("first_refusal", "") for r in observations if r["outcome"] == "refused")
        ),
        per_model_csv_sha256=digest(table),
        scope="readable source surfaces versus admitted solids versus exports valid after builder processing; no common source-fidelity proof",
    )
    save(args.out / "results.json", summary)
    print(
        json.dumps(
            {
                k: summary[k]
                for k in [
                    "complete",
                    "evaluated",
                    "outcomes",
                    "opened",
                    "full_solid_files",
                    "all_solids_export_valid_files",
                    "tier_losses",
                    "frozen_valid_export_losses",
                    "acceptance_passed",
                ]
            }
        ),
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("build")
    command.add_argument("--record", type=Path, required=True)
    command = commands.add_parser("raw")
    command.add_argument("--capture", type=Path, required=True)
    command.add_argument("--source-sha256", required=True)
    command.add_argument("--out", type=Path, required=True)
    command = commands.add_parser("collect")
    for field in ["build-record", "sources", "selection", "source-hashes", "out"]:
        command.add_argument("--" + field, type=Path, required=True)
    command.add_argument("--workers", type=int, default=2)
    command.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    {"build": build, "raw": raw, "collect": collect}[args.command](args)


if __name__ == "__main__":
    main()
