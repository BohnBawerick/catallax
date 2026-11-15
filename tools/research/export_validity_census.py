"""Measure raw STEP validity for the frozen 400 admissions, without changing admission.

Build records bind the actual Cargo artifact to a clean observer revision whose
production kernel matches the final reviewed exporter. Captures and paths stay
private. Publication retains all 3,000 original IDs and requires every frozen
admission to have a terminal observation. No partial export qualifies a file.
"""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
ADMISSIONS_SHA = "b08224a6f919405de46d35896d8d682febbd10c36495d34759f4958c80be3bd7"
HASHES_SHA = "3323e104c377023846073a05d7e164ca1082dc3dcf32cc046eee897401454cc7"
SELECTION_SHA = "1b11410b220d5e8a3156794f6c8f9cc09091a43dcac9af801215fb3a38d571a5"
RAW_HELPER_SHA = "5633cdcf16efc13f328ed85b33a738f6a9296e2afb925f4c64bbed7c91025376"
ARCHIVE_SHA = "52e6dd1b6fa38e3cd99af59b662370829129540030975919c3e256dce6ad1dbe"
TERMINAL = {
    "raw_valid",
    "raw_invalid",
    "raw_reader_refused",
    "export_refused",
    "reimport_refused",
    "admission_changed",
    "timeout",
    "process_error",
}
WORK = 20_000_000


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()


def population(admissions, hashes):
    rows = {r["model_id"]: r for r in admissions}
    assert len(rows) == len(admissions) == len(hashes) == 3000, "missing/duplicate original ID"
    assert set(rows) == set(hashes), "different original IDs"
    assert all(r["source_sha256"] == hashes[k] for k, r in rows.items()), "changed source hash"
    assert Counter(r["outcome"] for r in rows.values()) == Counter(
        full=400, partial=3, refused=2595, timeout=2
    )
    frozen = {k: r for k, r in rows.items() if r["outcome"] == "full"}
    assert all(int(r["bodies"]) > 0 and r["coverage_complete"] == "True" for r in frozen.values())
    return rows, frozen


def public_build(record):
    fields = "source_revision git_tree kernel_revision binary_sha256 cargo_json_sha256 cargo_lock_sha256 observer_sha256 collector_sha256 raw_helper_sha256 build_command rustc".split()
    result = {k: record[k] for k in fields}
    artifact = record["cargo_artifact"]
    result["cargo_artifact"] = {k: artifact[k] for k in ["profile", "features", "fresh"]}
    result["cargo_artifact"]["target"] = {
        k: artifact["target"][k] for k in ["name", "kind", "crate_types", "edition"]
    }
    return result


def build(args):
    assert not args.record.resolve().is_relative_to(ROOT)
    assert not git("status", "--porcelain"), "build requires a clean checkout"
    revision = git("rev-parse", "HEAD")
    assert len(args.kernel_revision) == 40 and all(
        c in "0123456789abcdef" for c in args.kernel_revision
    )
    # Observer/example changes are allowed; every production implementation and
    # dependency must remain the final reviewed kernel for this comparison.
    subprocess.run(
        [
            "git",
            "-C",
            str(ROOT),
            "diff",
            "--exit-code",
            args.kernel_revision,
            "HEAD",
            "--",
            ":(glob)crates/*/src/**",
            ":(glob)crates/*/Cargo.toml",
            "Cargo.toml",
            "Cargo.lock",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    command = [
        "cargo",
        "build",
        "--release",
        "-p",
        "catallax-import",
        "--example",
        "export_import",
        "--message-format=json",
    ]
    log = args.record.with_suffix(".cargo.jsonl")
    with log.open("xb") as out:
        subprocess.run(command, cwd=ROOT, stdout=out, check=True, timeout=3600)
    artifacts = [
        r
        for line in log.read_text().splitlines()
        if (r := json.loads(line))["reason"] == "compiler-artifact"
        and r["target"]["name"] == "export_import"
        and r["target"]["kind"] == ["example"]
        and r.get("executable")
    ]
    assert len(artifacts) == 1 and artifacts[0]["profile"]["test"] is False
    assert revision == git("rev-parse", "HEAD") and not git("status", "--porcelain")
    artifact = artifacts[0]
    binary = Path(artifact["executable"]).resolve()
    record = dict(
        format="catallax-export-observer-build/1",
        source_revision=revision,
        git_tree=git("rev-parse", "HEAD^{tree}"),
        kernel_revision=args.kernel_revision,
        artifact_path=str(binary),
        binary_sha256=digest(binary),
        cargo_artifact=artifact,
        cargo_json_sha256=digest(log),
        cargo_lock_sha256=digest(ROOT / "Cargo.lock"),
        observer_sha256=digest(ROOT / "crates/catallax-import/examples/export_import.rs"),
        collector_sha256=digest(Path(__file__)),
        raw_helper_sha256=RAW_HELPER_SHA,
        build_command=command,
        rustc=subprocess.check_output(["rustc", "--version"], text=True).strip(),
    )
    assert not args.record.exists()
    save(args.record, record)
    print(
        json.dumps(dict(source_revision=revision, binary_sha256=record["binary_sha256"])),
        flush=True,
    )


def raw(args):
    helper = ROOT / "tools/step-regressions/check_numerical_step.py"
    assert digest(helper) == RAW_HELPER_SHA
    sys.path.insert(0, str(helper.parent))
    import OCP
    from OCP.BRepCheck import BRepCheck_Analyzer
    from check_numerical_step import raw_solid, topology

    assert OCP.__version__ == "8.0.1.1"
    response = json.loads((args.capture / "response.json").read_text())
    assert response["stage"] == "export" and response["outcome"] == "complete"
    assert response["source_sha256"] == args.source_sha256
    files = response["files"]
    assert len(files) == response["bodies"] > 0
    assert [entry["body"] for entry in files] == list(range(len(files)))
    results = []
    for entry in files:
        assert entry["file"] == f"body-{entry['body']}.step"
        path = args.capture / entry["file"]
        before = path.read_bytes()
        assert digest(path) == entry["sha256"] and len(before) == entry["bytes"]
        try:
            compact = b"".join(before.split())
            assert b"SI_UNIT(.MILLI.,.METRE.)" in compact, "output unit not millimetres"
            shape = raw_solid(path)
            counts = topology(shape)
            assert counts["solids"] == 1, "one explicit solid per body required"
            valid = BRepCheck_Analyzer(shape).IsValid()
            status = "raw_valid" if valid else "raw_invalid"
        except Exception as error:
            status = "raw_reader_refused"
            counts = None
            print(type(error).__name__, str(error), file=sys.stderr)
        assert path.read_bytes() == before
        results.append(
            dict(body=entry["body"], sha256=entry["sha256"], outcome=status, topology=counts)
        )
    save(
        args.out,
        dict(
            format="catallax-raw-export-files/1",
            source_sha256=args.source_sha256,
            ocp_version=OCP.__version__,
            files=results,
        ),
    )


def observe(binary, source, baseline, folder, timeout):
    key, expected = baseline["model_id"], baseline["source_sha256"]
    row = dict(
        model_id=key,
        source_sha256=expected,
        outcome="process_error",
        expected_body_count=int(baseline["bodies"]),
    )
    assert digest(source) == expected, (key, "changed source bytes")
    folder.mkdir(exist_ok=True)
    attempt = folder / f"attempt-{len(list(folder.glob('attempt-*')))}"
    attempt.mkdir()
    capture = attempt / "capture"
    start = time.monotonic()
    with (attempt / "producer.stdout").open("wb") as out, (attempt / "producer.stderr").open(
        "wb"
    ) as err:
        try:
            process = subprocess.run(
                [str(binary), str(source), expected, str(capture)],
                stdout=out,
                stderr=err,
                timeout=timeout,
                check=False,
                env={**os.environ, "CATALLAX_THREADS": "1"},
            )
        except subprocess.TimeoutExpired:
            checkpoints = []
            for line in (attempt / "producer.stdout").read_text().splitlines():
                try:
                    value = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if value.get("stage") == "import" and value.get("outcome") == "admitted":
                    assert value["source_sha256"] == expected
                    checkpoints.append(value)
            assert len(checkpoints) <= 1
            if checkpoints:
                return row | dict(
                    outcome="timeout",
                    timeout_stage="export",
                    import_work_used=checkpoints[0]["work_used"],
                    import_category=checkpoints[0]["import_category"],
                    observed_body_count=checkpoints[0]["bodies"],
                )
            return row | dict(outcome="timeout", timeout_stage="import")
        finally:
            assert digest(source) == expected, (key, "source changed during observation")
    if process.returncode or not (capture / "response.json").is_file():
        return row | dict(cause="producer_process_or_missing_response")
    response = json.loads((capture / "response.json").read_text())
    assert response["source_sha256"] == expected
    if response["stage"] == "import":
        assert response["outcome"] == "refused"
        return row | dict(
            outcome="reimport_refused",
            cause=response["refusal"].get("code", "unknown"),
            import_work_used=response["work_used"],
        )
    admission = response["import"]
    assert admission["source_sha256"] == expected and admission["outcome"] == "admitted"
    assert response["work_limit"] == admission["work_limit"] == WORK
    assert 0 <= admission["work_used"] <= WORK and 0 <= response["work_used"] <= WORK
    row |= dict(
        import_work_used=admission["work_used"],
        export_work_used=response["work_used"],
        import_category=admission["import_category"],
        exported_body_count=len(response["files"]),
        expected_body_count=int(baseline["bodies"]),
        observed_body_count=admission["bodies"],
    )
    if admission["bodies"] != int(baseline["bodies"]):
        return row | dict(
            outcome="admission_changed", cause="body_count_differs_from_frozen_admission"
        )
    if response["outcome"] == "refused":
        cause = response["refusal"]["code"]
        return row | dict(
            outcome="process_error" if cause.startswith("observer.") else "export_refused",
            cause=cause,
        )
    assert response["outcome"] == "complete" and len(response["files"]) == admission["bodies"]
    remaining = timeout - (time.monotonic() - start)
    if remaining <= 0:
        return row | dict(outcome="timeout", timeout_stage="raw_reader")
    with (attempt / "raw.stdout").open("wb") as out, (attempt / "raw.stderr").open("wb") as err:
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
                    str(attempt / "raw.json"),
                ],
                stdout=out,
                stderr=err,
                timeout=remaining,
                check=False,
            )
        except subprocess.TimeoutExpired:
            return row | dict(outcome="timeout", timeout_stage="raw_reader")
    if process.returncode or not (attempt / "raw.json").is_file():
        return row | dict(cause="raw_process_or_missing_response")
    checked = json.loads((attempt / "raw.json").read_text())
    assert checked["source_sha256"] == expected and len(checked["files"]) == admission["bodies"]
    assert checked["ocp_version"] == "8.0.1.1"
    assert [f["sha256"] for f in checked["files"]] == [f["sha256"] for f in response["files"]]
    assert all(
        f["outcome"] in {"raw_valid", "raw_invalid", "raw_reader_refused"} for f in checked["files"]
    )
    outcomes = Counter(r["outcome"] for r in checked["files"])
    outcome = (
        "raw_reader_refused"
        if outcomes["raw_reader_refused"]
        else "raw_invalid" if outcomes["raw_invalid"] else "raw_valid"
    )
    return row | dict(
        outcome=outcome,
        raw_valid_solid_count=outcomes["raw_valid"],
        raw_invalid_solid_count=outcomes["raw_invalid"],
        raw_reader_refused_body_count=outcomes["raw_reader_refused"],
        exported_step_hashes=[r["sha256"] for r in checked["files"]],
    )


def summary(admissions, observations):
    expected = {k: r for k, r in admissions.items() if r["outcome"] == "full"}
    actual = {r["model_id"]: r for r in observations}
    assert len(actual) == len(observations), "duplicate validity observation"
    assert set(actual) <= set(expected), "non-frozen input measured"
    for key, row in actual.items():
        assert (
            row["source_sha256"] == expected[key]["source_sha256"]
        ), "changed validity source hash"
        assert row["outcome"] in TERMINAL, "unknown validity outcome"
        if row["outcome"] == "raw_valid":
            count = int(expected[key]["bodies"])
            assert (
                row["exported_body_count"] == row["raw_valid_solid_count"] == count > 0
            ), "partial export cannot qualify a file"
            assert len(row["exported_step_hashes"]) == count
    counts = dict(Counter(r["outcome"] for r in observations))
    complete = set(actual) == set(expected) and not counts.get("process_error", 0)
    return dict(
        complete=complete,
        evaluated=len(actual),
        frozen_admitted_inputs=len(expected),
        denominator=3000,
        seed=20261111,
        export_outcomes=counts,
        raw_valid_files=counts.get("raw_valid", 0),
        raw_valid_share_percent=f"{counts.get('raw_valid', 0) / 30:.2f}" if complete else None,
    )


def collect(args):
    assert args.timeout == 60 and args.workers > 0
    assert not args.out.resolve().is_relative_to(ROOT)
    assert digest(args.admissions) == ADMISSIONS_SHA and digest(args.source_hashes) == HASHES_SHA
    assert digest(args.selection) == SELECTION_SHA
    manifest = json.loads(args.admission_manifest.read_text())
    assert manifest["complete"] and manifest["denominator"] == manifest["evaluated"] == 3000
    assert manifest["seed"] == 20261111 and manifest["per_model_csv_sha256"] == ADMISSIONS_SHA
    assert (
        manifest["source_hash_manifest_sha256"] == HASHES_SHA
        and manifest["selection_sha256"] == SELECTION_SHA
    )
    assert manifest["archive_sha256"] == ARCHIVE_SHA
    assert manifest["source_revision"] == "51e9c0ffa13efcbcfb5c004b152e4c6014c6b5c1"
    csv.field_size_limit(64 * 1024 * 1024)
    with args.admissions.open(newline="") as stream:
        admissions, frozen = population(
            list(csv.DictReader(stream)), json.loads(args.source_hashes.read_text())
        )
    selection = dict(json.loads(args.selection.read_text()))
    assert set(selection) == set(admissions)
    files = {}
    for key in frozen:
        member = PurePosixPath(selection[key])
        assert not member.is_absolute() and ".." not in member.parts
        path = (args.sources / str(member)).resolve()
        assert (
            path.is_relative_to(args.sources.resolve())
            and digest(path) == frozen[key]["source_sha256"]
        )
        files[key] = path
    record = json.loads(args.build_record.read_text())
    assert record["format"] == "catallax-export-observer-build/1"
    assert record["source_revision"] == git("rev-parse", "HEAD") and not git(
        "status", "--porcelain"
    )
    assert record["git_tree"] == git("rev-parse", "HEAD^{tree}")
    binary = Path(record["artifact_path"])
    assert (
        digest(binary) == record["binary_sha256"]
        and binary.resolve() == Path(record["cargo_artifact"]["executable"]).resolve()
    )
    assert record["observer_sha256"] == digest(
        ROOT / "crates/catallax-import/examples/export_import.rs"
    )
    assert record["collector_sha256"] == digest(Path(__file__))
    assert record["cargo_lock_sha256"] == digest(ROOT / "Cargo.lock")
    assert (
        record["raw_helper_sha256"]
        == digest(ROOT / "tools/step-regressions/check_numerical_step.py")
        == RAW_HELPER_SHA
    )
    assert record["cargo_artifact"]["target"]["name"] == "export_import" and record[
        "cargo_artifact"
    ]["target"]["kind"] == ["example"]
    assert record["cargo_artifact"]["profile"]["test"] is False
    identity = dict(
        format="catallax-frozen-admission-export-validity/1",
        source_revision=record["source_revision"],
        kernel_revision=record["kernel_revision"],
        binary_sha256=record["binary_sha256"],
        build_record=public_build(record),
        build_record_sha256=digest(args.build_record),
        admission_csv_sha256=ADMISSIONS_SHA,
        source_hash_manifest_sha256=HASHES_SHA,
        admission_source_revision=manifest["source_revision"],
        admission_manifest_sha256=digest(args.admission_manifest),
        archive_sha256=ARCHIVE_SHA,
        selection_sha256=SELECTION_SHA,
        raw_helper_sha256=RAW_HELPER_SHA,
        per_import_work_limit=WORK,
        aggregate_export_work_limit=WORK,
        observation_timeout_seconds=60,
        workers=args.workers,
        per_import_workers=1,
        ocp_version="8.0.1.1",
        raw_precision_mm="1e-7",
        raw_max_tolerance_mm="1e-7",
        output_unit_factors=[1, 1, 1],
        exported_step_healing=False,
        exported_step_same_parameter=False,
        source_coverage="frozen_full_admissions_only",
        scope="raw exported-solid validity; not source fidelity, arbitrary edit, material or full-body deviation qualification",
    )
    args.out.mkdir(exist_ok=args.resume)
    identity_path = args.out / "identity.json"
    if args.resume:
        assert json.loads(identity_path.read_text()) == identity
    else:
        save(identity_path, identity)
    rows_dir = args.out / "rows"
    rows_dir.mkdir(exist_ok=True)
    parts = args.out / "parts"
    parts.mkdir(exist_ok=True)
    observations = {}
    for path in rows_dir.glob("*.json"):
        row = json.loads(path.read_text())
        assert path.stem == row["model_id"] and row["model_id"] not in observations
        if row["outcome"] != "process_error":
            observations[row["model_id"]] = row
    summary(admissions, list(observations.values()))
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs = {
            pool.submit(observe, binary, files[key], frozen[key], parts / key, args.timeout): key
            for key in sorted(frozen)
            if key not in observations
        }
        for job in as_completed(jobs):
            key = jobs[job]
            row = job.result()
            save(rows_dir / (key + ".json"), row)
            observations[key] = row
            if len(observations) % 25 == 0:
                print(
                    json.dumps(
                        dict(
                            evaluated=len(observations),
                            outcomes=dict(Counter(r["outcome"] for r in observations.values())),
                        )
                    ),
                    flush=True,
                )
    assert digest(binary) == record["binary_sha256"] and not git("status", "--porcelain")
    assert all(
        digest(files[key]) == frozen[key]["source_sha256"] for key in frozen
    ), "source set changed during capture"
    report = identity | summary(admissions, list(observations.values()))
    save(args.out / "public-summary.json", report)
    fields = "model_id source_sha256 admission_outcome baseline_first_refusal admission_after_fix_cause validity_outcome import_category expected_body_count observed_body_count exported_body_count raw_valid_solid_count raw_invalid_solid_count raw_reader_refused_body_count import_work_used export_work_used cause timeout_stage exported_step_hashes".split()
    with (args.out / "results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for key, admitted in sorted(admissions.items()):
            observed = observations.get(key, {})
            row = {k: observed.get(k, "") for k in fields}
            row.update(
                model_id=key,
                source_sha256=admitted["source_sha256"],
                admission_outcome=admitted["outcome"],
                baseline_first_refusal=admitted["baseline_first_refusal"],
                admission_after_fix_cause=admitted["after_fix_cause"],
                validity_outcome=observed.get("outcome", "outside_frozen_full_population"),
            )
            row["exported_step_hashes"] = json.dumps(
                observed.get("exported_step_hashes", []), separators=(",", ":")
            )
            writer.writerow(row)
    report["per_model_csv_sha256"] = digest(args.out / "results.csv")
    save(args.out / "public-summary.json", report)
    print(json.dumps(report, sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("build")
    command.add_argument("--record", type=Path, required=True)
    command.add_argument("--kernel-revision", required=True)
    command = commands.add_parser("raw")
    command.add_argument("--capture", type=Path, required=True)
    command.add_argument("--source-sha256", required=True)
    command.add_argument("--out", type=Path, required=True)
    command = commands.add_parser("collect")
    for field in ["build-record", "sources", "selection", "source-hashes", "out"]:
        command.add_argument("--" + field, type=Path, required=True)
    command.add_argument(
        "--admissions", type=Path, default=ROOT / "docs/research/step-census-periodic/results.csv"
    )
    command.add_argument(
        "--admission-manifest",
        type=Path,
        default=ROOT / "docs/research/step-census-periodic/results.json",
    )
    command.add_argument("--workers", type=int, default=2)
    command.add_argument("--timeout", type=int, default=60)
    command.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    {"build": build, "raw": raw, "collect": collect}[args.command](args)


if __name__ == "__main__":
    main()
