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
import rss_guard

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


def step_reference_graph(data):
    """Read serialized references without interpreting strings or report comments."""
    import re
    graph, current = {}, None
    tokens = re.finditer(rb"'(?:''|[^'])*'|/\*.*?\*/|#([0-9]+)(\s*=)?|;", data, re.S)
    for token in tokens:
        if token[1] is not None:
            target = int(token[1])
            if token[2] is not None:
                current = target
                graph[current] = set()
            elif current is not None:
                graph[current].add(target)
        elif token[0] == b";":
            current = None
    return graph


def product_reachable_entities(model, data, identifiers):
    """Bind OCCT product roots to references in the actual emitted STEP bytes."""
    from OCP.StepShape import StepShape_ShapeDefinitionRepresentation

    # Keep wrappers alive while binding typed references. OCP 8.0.1.1 returns
    # zero for both model.Number and IdentLabel, including model.Value entries.
    entities = [model.Value(index) for index in range(1, model.NbEntities() + 1)]
    labels = dict(zip(entities, identifiers))
    references = step_reference_graph(data)
    assert set(references) == set(identifiers)

    def product_bound(property_definition):
        pending, seen = [property_definition], set()
        while pending:
            value = pending.pop()
            if value is None or labels[value] in seen:
                continue
            seen.add(labels[value])
            definition = value.Definition()
            if (definition.ProductDefinition() is not None
                    or definition.ProductDefinitionRelationship() is not None):
                return True
            nested = definition.ProductDefinitionShape()
            if nested is not None:
                pending.append(nested)
        return False

    pending = []
    for entity in entities:
        if isinstance(entity, StepShape_ShapeDefinitionRepresentation):
            definition = entity.Definition().PropertyDefinition()
            if definition is not None and product_bound(definition):
                pending.append(labels[entity.UsedRepresentation()])
    reached = set()
    while pending:
        identifier = pending.pop()
        if identifier in reached:
            continue
        reached.add(identifier)
        pending.extend(references[identifier])
    return reached


def raw_native_document(path, entry):
    """Check every named native root before default occurrence transfer.

    Prototype checks use raw StepToTopoDS, explicit source dimensional factors
    and the stated maximum local tolerance of that source body. Default
    TransferRoots is a separate occurrence-graph check, not a raw-body rescue.
    """
    import math
    from OCP.BRep import BRep_Builder
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.Interface import Interface_Static
    from OCP.STEPControl import STEPControl_Reader
    from OCP.StepData import StepData_Factors
    from OCP.StepShape import StepShape_ManifoldSolidBrep
    from OCP.StepToTopoDS import StepToTopoDS_Builder
    from OCP.TopoDS import TopoDS_Compound
    from OCP.Transfer import Transfer_TransientProcess
    from check_numerical_step import topology

    reader = STEPControl_Reader()
    settings = {"read.stdsameparameter.mode": 0, "read.maxprecision.mode": 0}
    previous = {key: Interface_Static.IVal_s(key) for key in settings}
    try:
        for key, value in settings.items():
            assert Interface_Static.SetIVal_s(key, value), key
        assert reader.ReadFile(str(path)) == IFSelect_RetDone
        model = reader.StepModel()
        # This writer emits one DATA instance per line, in entity order.
        # OCP 8.0.1.1 leaves IdentLabel unset (zero), so bind the actual P21
        # instance numbers to reader record ordinals, not to optional labels.
        import re
        source_bytes = path.read_bytes()
        identifiers = [int(value) for value in re.findall(rb"^#([0-9]+)=", source_bytes, re.M)]
        assert len(identifiers) == model.NbEntities() and len(set(identifiers)) == len(identifiers)
        reached_identifiers = product_reachable_entities(model, source_bytes, identifiers)
        labelled = {}
        for index, identifier in enumerate(identifiers, 1):
            entity = model.Value(index)
            if isinstance(entity, StepShape_ManifoldSolidBrep):
                labelled[identifier] = entity
        roots = []
        for root in entry["source_roots"]:
            factors_values = root.get("unit_factors", entry.get("unit_factors"))
            assert len(factors_values) == 3 and all(math.isfinite(v) and v > 0 for v in factors_values)
            counts, valid_solids, outcome, detail = None, 0, "reader_refused", ""
            inspected_parts = 0
            try:
                entities = [labelled[entity] for entity in root["entities"]]
                assert len(entities) == root["parts"], (root["name"], len(entities), root["parts"])
                tolerance = root["tolerance_mm"]
                assert math.isfinite(tolerance) and tolerance > 0
                factors = StepData_Factors()
                factors.InitializeFactors(*factors_values)
                process = Transfer_TransientProcess()
                process.SetModel(model)
                compound = TopoDS_Compound()
                aggregate = BRep_Builder()
                aggregate.MakeCompound(compound)
                valid = True
                for entity in entities:
                    builder = StepToTopoDS_Builder()
                    builder.SetPrecision(min(tolerance, 1e-7))
                    builder.SetMaxTol(tolerance)
                    builder.Init(entity, process, factors)
                    assert builder.IsDone(), root["name"]
                    shape = builder.Value()
                    assert topology(shape)["solids"] == 1
                    part_valid = BRepCheck_Analyzer(shape).IsValid()
                    valid_solids += int(part_valid)
                    valid &= part_valid
                    inspected_parts += 1
                    aggregate.Add(compound, shape)
                counts = topology(compound)
                coverage = all(counts[key] == expected for key, expected in root["counts"].items())
                outcome = "valid_after_builder" if valid and coverage else "invalid_after_builder"
                if not coverage:
                    detail = "native topology counts differ from the emitted body"
            except Exception as error:
                detail = f"{type(error).__name__}: {error}"
            roots.append(dict(body=root["body"], source_root=root["name"], topology=counts,
                              outcome=outcome, valid_solid_count=valid_solids, detail=detail,
                              occurrence_reachable=all(entity in reached_identifiers for entity in root["entities"]),
                              parts_expected=root["parts"], parts_inspected=inspected_parts,
                              observation_complete=inspected_parts == root["parts"]))
    finally:
        for key, value in previous.items():
            Interface_Static.SetIVal_s(key, value)
    occurrence_outcome, counts = "reader_refused", None
    try:
        expected_roots = reader.NbRootsForTransfer()
        assert expected_roots > 0 and reader.TransferRoots() == expected_roots
        shape = reader.OneShape()
        assert not shape.IsNull() and topology(shape)["solids"] > 0
        counts = topology(shape)
        occurrence_outcome = "valid_after_builder" if BRepCheck_Analyzer(shape).IsValid() else "invalid_after_builder"
    except Exception as error:
        print(type(error).__name__, str(error), file=sys.stderr)
    valid = occurrence_outcome == "valid_after_builder" and all(
        root["outcome"] == "valid_after_builder" for root in roots
    )
    return counts, roots, occurrence_outcome, valid


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
        bodies = entry["bodies"]
        assert bodies and len(set(bodies)) == len(bodies)
        native = entry["kind"] in ("native_document", "native_body", "selected_document")
        assert entry["file"] == ("document.step" if entry["kind"] in ("native_document", "selected_document") else f"body-{bodies[0]}.step")
        path = args.capture / entry["file"]
        assert digest(path) == entry["sha256"] and path.stat().st_size == entry["bytes"]
        before = path.read_bytes()
        valid_solids, native_roots, occurrence_outcome = 0, [], None
        try:
            if native:
                counts, native_roots, occurrence_outcome, valid = raw_native_document(path, entry)
                assert sorted(root["body"] for root in native_roots) == sorted(bodies)
                valid_solids = sum(root["valid_solid_count"] for root in native_roots)
                outcome = "valid_after_builder" if valid else "invalid_after_builder"
            else:
                assert len(bodies) == 1
                assert b"SI_UNIT(.MILLI.,.METRE.)" in b"".join(before.split())
                shape = raw_solid(path) if entry["kind"] == "solid" else raw_surface(path)
                counts = topology(shape)
                assert counts["solids"] == (1 if entry["kind"] == "solid" else 0)
                outcome = (
                    "valid_after_builder"
                    if BRepCheck_Analyzer(shape).IsValid()
                    else "invalid_after_builder"
                )
                valid_solids = int(entry["kind"] == "solid" and outcome == "valid_after_builder")
        except Exception as error:
            outcome, counts = "reader_refused", None
            print(type(error).__name__, str(error), file=sys.stderr)
        assert path.read_bytes() == before
        results.append(
            dict(
                bodies=bodies,
                kind=entry["kind"],
                sha256=entry["sha256"],
                outcome=outcome,
                topology=counts,
                valid_solid_count=valid_solids,
                raw_source_roots=native_roots,
                occurrence_outcome=occurrence_outcome,
                occurrence_graph_complete=bool(native_roots) and all(root["occurrence_reachable"] for root in native_roots),
                observation_complete=(
                    outcome != "reader_refused"
                    and (not native or bool(native_roots)
                         and all(root["observation_complete"] for root in native_roots))),
                validation_method=(
                    "raw_source_roots_at_stated_tolerance_then_default_occurrence_transfer"
                    if native else "raw_identity_unit_builder_at_fixed_precision"
                ),
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
    solids = tiers["valid_solid"] + tiers["repaired_solid"] + tiers["tolerant_solid"]
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
        tolerant_solid_body_count=tiers["tolerant_solid"],
        surface_body_count=tiers["surface_body"],
        open_work_used=checkpoint["work_used"],
        open_work_exhausted=checkpoint["work_exhausted"],
        read_work_used=checkpoint["read_work_used"],
        file_admission_work_used=checkpoint["file_admission_work_used"],
        native_admission_work_used=checkpoint.get("native_admission_work_used"),
        checked_upgrade_work_used=checkpoint.get("checked_upgrade_work_used"),
        first_refusal=cause(checkpoint["first_refusal"]),
    )
    row["proofs"] = dict(Counter(b["proof"] for b in bodies))
    row["body_first_blockers"] = [
        dict(body=body["body"], source_body=body["source_body"],
             stage=blocker["stage"], refusal=cause(blocker["refusal"]),
             resource_limited=blocker["resource_limited"], work_used=blocker["work_used"])
        for body in bodies
        if (blocker := body.get("first_blocker"))
    ]
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


def retain_opening_checkpoint(row, stdout, expected):
    checkpoints = []
    for line in stdout.read_text().splitlines():
        try:
            value = json.loads(line)
        except json.JSONDecodeError:
            continue
        if value.get("stage") == "open" and value.get("outcome") == "opened":
            checkpoints.append(value)
    if not checkpoints:
        return False
    assert len(checkpoints) == 1 and checkpoints[0]["source_sha256"] == expected
    opening(row, checkpoints[0])
    return True


def observe(binary, source, baseline, folder, retain, max_rss_mb=6144, timeout=60, work=WORK,
            raw_timeout=None):
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
        has_valid_exported_solid=None,
        any_valid_exported_solid_count=0,
        any_valid_exported_solid_count_is_lower_bound=True,
        export_observation_complete=False,
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
            process = rss_guard.run(
                command,
                stdout=out,
                stderr=err,
                timeout=timeout,
                max_rss_mb=max_rss_mb,
                env={**os.environ, "CATALLAX_THREADS": "1", "CATALLAX_OBSERVER_WORK": str(work)},
            )
        except (subprocess.TimeoutExpired, rss_guard.MemoryExceeded) as error:
            termination = "out_of_memory" if isinstance(error, rss_guard.MemoryExceeded) else "timeout"
            opened = retain_opening_checkpoint(row, folder / "producer.stdout", expected)
            row.update(
                timeout_stage="export" if opened else "open",
                export_validity=termination,
            )
            if isinstance(error, rss_guard.MemoryExceeded):
                row["observer_peak_rss_mb"] = error.peak_mb
            if not opened:
                row["outcome"] = termination
            return row
        finally:
            assert digest(source) == expected
    if process.returncode or not (capture / "response.json").is_file():
        if retain_opening_checkpoint(row, folder / "producer.stdout", expected):
            row["export_validity"] = "process_error"
        return row | dict(cause="producer_process_or_missing_response")
    response = json.loads((capture / "response.json").read_text())
    assert response["source_sha256"] == expected and response["work_limit"] == work
    if response["stage"] == "open":
        assert response["outcome"] == "refused"
        return row | dict(
            outcome="refused",
            first_refusal=cause(response["refusal"]),
            after_fix_causes=[cause(response["refusal"])],
            open_work_used=response["work_used"],
            open_work_exhausted=response["work_exhausted"],
            export_validity="not_opened",
            has_valid_exported_solid=False,
            any_valid_exported_solid_count=0,
            any_valid_exported_solid_count_is_lower_bound=False,
            export_observation_complete=True,
            source_document_exported=False,
        )
    opening(row, response["open"])
    row.update(
        export_work_used=response["work_used"],
        exported_body_count=len({body for entry in response["files"] for body in entry["bodies"]}),
        exported_step_hashes=[f["sha256"] for f in response["files"]],
    )
    if response["outcome"] == "refused":
        row["export_refusal"] = cause(response["refusal"]) or "resource_limit"
        row["non_budget_whole_document_export_refusal"] = document_export_refusal(response)
    if not response["files"]:
        row.update(export_inventory(
            [], row["observed_body_count"], row["outcome"] == "full",
            response["outcome"] == "complete"))
        row["export_validity"] = (
            "source_surfaces_without_checked_brep"
            if not row.get("export_refusal")
            else "export_refused"
        )
        return row
    remaining = raw_timeout if raw_timeout is not None else timeout - (time.monotonic() - start)
    if remaining <= 0:
        return row | dict(export_validity="timeout", timeout_stage="raw_reader")
    with (folder / "raw.stdout").open("wb") as out, (folder / "raw.stderr").open("wb") as err:
        try:
            process = rss_guard.run(
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
                max_rss_mb=max_rss_mb,
            )
        except (subprocess.TimeoutExpired, rss_guard.MemoryExceeded) as error:
            termination = "out_of_memory" if isinstance(error, rss_guard.MemoryExceeded) else "timeout"
            return row | dict(export_validity=termination, timeout_stage="raw_reader")
    if process.returncode or not (folder / "raw.json").is_file():
        return row | dict(export_validity="process_error", cause="raw_process_or_missing_response")
    checked = json.loads((folder / "raw.json").read_text())
    assert checked["source_sha256"] == expected and checked["ocp_version"] == "8.0.1.1"
    assert [f["sha256"] for f in checked["files"]] == row["exported_step_hashes"]
    row.update(export_inventory(
        checked["files"], row["observed_body_count"], row["outcome"] == "full",
        response["outcome"] == "complete"))
    row["export_validity"] = (
        "export_refused" if row.get("export_refusal")
        else "raw_observation_incomplete" if not row["export_observation_complete"]
        else "reader_refused" if row["reader_refused_body_count"]
        else "invalid_after_builder" if row["invalid_exported_body_count"]
        else "all_solids_valid_after_builder" if row["all_solids_export_valid"]
        else "valid_exported_subset"
    )
    return row


def document_export_refusal(response):
    refusal = response.get("refusal")
    if (response.get("stage") == "export" and refusal
            and response["failed_body"] is None and not response["files"]
            and response["open"]["bodies"] and not response["work_exhausted"]
            and refusal["code"] != "numbers.operation-budget"):
        return cause(refusal)
    return None


def export_inventory(files, body_count, all_bodies_solid, observation_complete):
    observation_complete = observation_complete and all(
        entry.get("observation_complete", entry["outcome"] != "reader_refused")
        and all(root.get("observation_complete", root["outcome"] != "reader_refused")
                for root in entry.get("raw_source_roots", []))
        for entry in files
    )
    counts = Counter((entry["kind"], entry["outcome"])
                     for entry in files for _ in entry["bodies"])
    solid = sum(counts[kind, "valid_after_builder"]
                for kind in ("solid", "native_document", "native_body", "selected_document"))
    any_solid = sum(entry["valid_solid_count"] for entry in files)
    document = any(entry["kind"] in ("native_document", "selected_document")
                   and entry.get("occurrence_graph_complete") is True for entry in files)
    prototypes_valid = all_bodies_solid and solid == body_count
    return dict(
        valid_exported_solid_count=solid,
        valid_exported_surface_count=counts["surface", "valid_after_builder"],
        invalid_exported_body_count=sum(v for (_, outcome), v in counts.items()
                                        if outcome == "invalid_after_builder"),
        reader_refused_body_count=sum(v for (_, outcome), v in counts.items()
                                      if outcome == "reader_refused"),
        any_valid_exported_solid_count=any_solid,
        has_valid_exported_solid=(
            any_solid > 0 if observation_complete or any_solid > 0 else None),
        any_valid_exported_solid_count_is_lower_bound=not observation_complete,
        export_observation_complete=observation_complete,
        source_document_exported=document,
        prototype_solids_export_valid=prototypes_valid,
        all_solids_export_valid=(
            prototypes_valid and observation_complete and document),
    )


def recorded_tiers(row):
    """Read an actual per-file tier inventory, never infer it from 'full'."""
    try:
        total = int(row["observed_body_count"])
        tiers = {
            "valid_solid": int(row["valid_solid_body_count"]),
            "repaired_solid": int(row["repaired_solid_body_count"]),
            "tolerant_solid": int(row.get("tolerant_solid_body_count") or 0),
            "surface_body": int(row["surface_body_count"]),
        }
    except (KeyError, TypeError, ValueError):
        return None
    return tiers if min(tiers.values()) >= 0 and sum(tiers.values()) == total else None


def comparison_outcomes(by_id, expected_full, expected_valid, prior_rows):
    unknown = {"timeout", "out_of_memory", "process_error", "raw_observation_incomplete"}
    opening_unknown = {key for key, row in by_id.items() if row["outcome"] in unknown}
    export_unknown = opening_unknown | {
        key for key, row in by_id.items() if row.get("export_validity") in unknown
    }
    transitions = {}
    downgrades, cooccurrences = [], {}
    for key, row in by_id.items():
        before = recorded_tiers(prior_rows.get(key, {}))
        after = recorded_tiers(row)
        if before is None:
            transition = "unknown_prior_tiers"
        elif after is None or key in opening_unknown:
            transition = "unknown_current_tiers"
        elif sum(before.values()) != sum(after.values()):
            transition = "unknown_body_inventory"
        else:
            checked_before = before["valid_solid"] + before["repaired_solid"]
            checked_after = after["valid_solid"] + after["repaired_solid"]
            checked_drop = checked_before - checked_after
            native_gain = after["tolerant_solid"] - before["tolerant_solid"]
            lowered = (
                before["valid_solid"] > after["valid_solid"]
                or checked_drop > 0
                or checked_before + before["tolerant_solid"]
                   > checked_after + after["tolerant_solid"]
            )
            if lowered:
                downgrades.append(key)
            if checked_drop > 0 and native_gain > 0:
                transition = "checked_loss_with_tolerant_gain"
                cooccurrences[key] = dict(checked_loss=checked_drop, tolerant_gain=native_gain)
            else:
                transition = (
                    "downgraded" if lowered else
                    "same_recorded_tier_totals" if before == after else "changed"
                )
        transitions[key] = dict(
            prior_body_tiers=before, current_body_tiers=after, tier_transition=transition
        )
    coverage_losses = {
        key for key, count in expected_full.items()
        if key not in opening_unknown
        and (by_id[key]["outcome"] != "full"
             or by_id[key].get("observed_body_count") != count)
    }
    unknown_tiers = sorted(
        key for key, value in transitions.items()
        if value["tier_transition"].startswith("unknown_")
    )
    return dict(
        solid_coverage_losses=sorted(coverage_losses),
        tier_losses=sorted(coverage_losses | set(downgrades)),
        recorded_tier_downgrades=sorted(downgrades),
        checked_loss_tolerant_gain=dict(sorted(cooccurrences.items())),
        checked_loss_tolerant_gain_files=len(cooccurrences),
        tier_comparison_unknown=unknown_tiers,
        tier_comparison_unknown_count=len(unknown_tiers),
        tier_transition_counts=dict(Counter(v["tier_transition"] for v in transitions.values())),
        tier_comparison_basis="aggregate per-file tier counts; co-occurrence does not establish same-body transitions",
        per_file_tier_transitions=transitions,
        frozen_valid_export_losses=sorted(
            key for key in expected_valid
            if key not in export_unknown and not by_id[key].get("prototype_solids_export_valid")
        ),
        not_measured=sorted(export_unknown),
        baseline_not_measured=sorted(
            (set(expected_full) & opening_unknown) | (expected_valid & export_unknown)
        ),
    )


def collect(args):
    assert args.workers == 2 and not args.out.resolve().is_relative_to(ROOT)
    assert args.max_rss_mb > 0
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
    tier_baseline = baseline
    if args.tier_baseline is not None:
        with args.tier_baseline.open(newline="") as stream:
            tier_baseline = {row["model_id"]: row for row in csv.DictReader(stream)}
        assert tier_baseline.keys() == hashes.keys()
        assert all(tier_baseline[key]["source_sha256"] == hashes[key] for key in hashes)
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
        max_rss_mb_per_observer=args.max_rss_mb,
        ocp_version="8.0.1.1",
        tier_baseline_csv_sha256=digest(args.tier_baseline or BASELINE),
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
                observe, binary, files[k], baseline[k], args.out / "parts" / k, k in retained,
                args.max_rss_mb,
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
    by_id = {r["model_id"]: r for r in observations}
    expected_full = {k: int(r["bodies"]) for k, r in baseline.items() if r["outcome"] == "full"}
    expected_full.update(
        {r["model_id"]: r["expected_body_count"] for r in cone_reference["results"]}
    )
    expected_valid = {k for k, r in validity.items() if r["validity_outcome"] == "raw_valid"}
    expected_valid.update(r["model_id"] for r in cone_reference["results"])
    assert len(expected_full) == 405 and len(expected_valid) == 374
    expected_full.update({
        key: int(row["observed_body_count"])
        for key, row in tier_baseline.items()
        if row.get("outcome") == "full" and row.get("observed_body_count")
    })
    expected_valid.update(
        key for key, row in tier_baseline.items()
        if row.get("all_solids_export_valid") in (True, "True")
    )
    comparison = comparison_outcomes(by_id, expected_full, expected_valid, tier_baseline)
    for key, transition in comparison.pop("per_file_tier_transitions").items():
        by_id[key].update(transition)
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
    blockers = {}
    for row in observations:
        for blocker in row.get("body_first_blockers", []):
            key = blocker["stage"], blocker["refusal"], blocker["resource_limited"]
            count, files = blockers.setdefault(key, [0, set()])
            blockers[key][0] = count + 1
            files.add(row["model_id"])
    summary = identity | dict(
        complete=True,
        evaluated=3000,
        outcomes=dict(outcomes),
        export_outcomes=dict(exports),
        opened=opened,
        opened_share_percent=f"{opened/30:.2f}",
        full_solid_files=outcomes["full"],
        full_solid_share_percent=f"{outcomes['full']/30:.2f}",
        files_with_any_valid_exported_solid=sum(r.get("has_valid_exported_solid") is True for r in observations),
        any_valid_export_unmeasured_files=[r["model_id"] for r in observations if r.get("has_valid_exported_solid") is None],
        any_valid_export_lower_bound_files=[r["model_id"] for r in observations if r.get("any_valid_exported_solid_count_is_lower_bound", r.get("has_valid_exported_solid") is None)],
        valid_solid_count_basis="individual observed raw-builder-valid exported solids; incomplete export observation is a lower bound, not a negative verdict on skipped siblings",
        prototype_solids_export_valid_files=sum(bool(r.get("prototype_solids_export_valid")) for r in observations),
        source_document_exported_files=sum(bool(r.get("source_document_exported")) for r in observations),
        strict_source_export_basis="complete source-body solid coverage, occurrence-preserving document export and all exported solids valid; isolated prototype exports are excluded",
        non_budget_whole_document_export_failures=dict(Counter(
            r["non_budget_whole_document_export_refusal"] for r in observations
            if r.get("non_budget_whole_document_export_refusal"))),
        non_budget_whole_document_export_failure_files=[
            r["model_id"] for r in observations if r.get("non_budget_whole_document_export_refusal")],
        body_first_blockers=[
            dict(stage=stage, refusal=refusal, resource_limited=resource, bodies=count, files=len(files))
            for (stage, refusal, resource), (count, files) in sorted(blockers.items())
        ],
        body_first_blockers_unmeasured_files=[r["model_id"] for r in observations if "body_first_blockers" not in r],
        all_solids_export_valid_files=qualified,
        all_solids_export_valid_share_percent=f"{qualified/30:.2f}",
        body_tiers={
            name: sum(r.get(name + "_body_count", 0) for r in observations)
            for name in ["valid_solid", "repaired_solid", "tolerant_solid", "surface"]
        },
        **comparison,
        baseline_full_solid_files=len(expected_full),
        baseline_valid_export_files=len(expected_valid),
        integrated_00009338_preserved=preserved_9338,
        acceptance_passed=not comparison["tier_losses"]
        and not comparison["frozen_valid_export_losses"]
        and not comparison["baseline_not_measured"]
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
                    "checked_loss_tolerant_gain_files",
                    "tier_comparison_unknown_count",
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
    command.add_argument("--tier-baseline", type=Path)
    command.add_argument("--workers", type=int, default=2)
    command.add_argument("--max-rss-mb", type=int, default=6144)
    command.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    {"build": build, "raw": raw, "collect": collect}[args.command](args)


if __name__ == "__main__":
    main()
