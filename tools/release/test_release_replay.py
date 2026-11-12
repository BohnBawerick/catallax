"""Replay policy applies to saved manifests, inputs and exact output bytes."""
import copy
import importlib.util
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import release


class ReplayTransitionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location(
            "release_test_driver", release.ROOT / ".claude/skills/verify-catallax/scripts/control-catallax.py")
        cls.driver = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.driver)

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def case(self, name, status="answered", step=b"exact STEP bytes", value=1, step_fields=None):
        directory = self.root / name
        directory.mkdir()
        manifest = {"status": status, "result": value if status == "answered" else None,
                    "refusal": {"code": "numbers.operation-budget"} if status == "refused" else None,
                    "kernel": {"commit": name}, "stats": {"wall_ms": len(name)}}
        if step_fields is not None:
            manifest["step"] = {"file": "answer.step", "added_edges": 0, "added_vertices": 0,
                                "max_deviation": "0", **step_fields}
        release.write_json(directory / "manifest.json", manifest)
        release.write_json(directory / "request.json", {"problem": {"value": 1}, "max_work": 100,
                                                       "out": str(directory)})
        if status == "answered":
            (directory / "answer.step").write_bytes(step)
        return {"path": str(directory), "id": "synthetic-case", "response": {"work": {"used": 1}}}

    def report(self, name, cases):
        path = self.root / (name + ".json")
        release.write_json(path, {"cases": cases, "config": {"binary_sha256": name}})
        return path

    def approval(self, before, after, index=0):
        request = release.read_json(Path(before["path"]) / "request.json")
        request.pop("out")

        def answer(case):
            directory = Path(case["path"])
            return self.driver.digest({
                "manifest": self.driver.comparable(release.read_json(directory / "manifest.json")),
                "files": {name: release.sha256(directory / name) if (directory / name).exists() else None
                          for name in release.OUTPUTS}})

        return {"index": index, "request_sha256": self.driver.digest(request),
                "refusal_sha256": answer(before), "answer_sha256": answer(after)}

    def compare(self, before, after, upgrades=(), corrections=()):
        return release.compare_baseline(self.driver, self.report("base", before),
                                        self.report("candidate", after), list(upgrades), list(corrections))

    def test_preserves_earlier_answers_and_unchanged_refusals(self):
        result = self.compare([self.case("a"), self.case("b", "refused")],
                              [self.case("c"), self.case("d", "refused")])
        self.assertEqual(result, {"cases": 2, "equal": True, "failures": [], "upgrades": 0, "corrections": 0})

    def test_listed_refusal_upgrade_pins_the_new_answer(self):
        before, after = self.case("before", "refused"), self.case("after")
        result = self.compare([before], [after], [self.approval(before, after)])
        self.assertEqual(result, {"cases": 1, "equal": False, "failures": [], "upgrades": 1, "corrections": 0})

    def test_changed_answer_is_rejected_even_if_listed(self):
        before, after = self.case("before"), self.case("after", value=2)
        with self.assertRaisesRegex(ValueError, "replay differs"):
            self.compare([before], [after], [self.approval(before, after)])

    def test_answer_now_refuses_is_rejected_even_if_listed(self):
        before, after = self.case("before"), self.case("after", "refused")
        with self.assertRaisesRegex(ValueError, "replay differs"):
            self.compare([before], [after], [self.approval(before, after)])

    def test_unlisted_upgrade_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "unlisted upgrade"):
            self.compare([self.case("before", "refused")], [self.case("after")])

    def test_upgrade_cannot_hide_input_or_allowance_changes(self):
        for field in ("problem", "max_work"):
            with self.subTest(field=field):
                before = self.case(field + "before", "refused")
                after = self.case(field + "after")
                approval = self.approval(before, after)
                path = Path(after["path"]) / "request.json"
                request = release.read_json(path)
                request[field] = 200
                release.write_json(path, request)
                with self.assertRaisesRegex(ValueError, "replay differs"):
                    self.compare([before], [after], [approval])

    def test_approved_answer_must_match_manifest_and_all_output_bytes(self):
        for field in ("manifest.json", *release.OUTPUTS):
            with self.subTest(field=field):
                before = self.case(field + "before", "refused")
                after = self.case(field + "after")
                approval = self.approval(before, after)
                path = Path(after["path"]) / field
                if field == "manifest.json":
                    manifest = release.read_json(path)
                    manifest["result"] = 2
                    release.write_json(path, manifest)
                else:
                    path.write_bytes(b"changed bytes")
                with self.assertRaisesRegex(ValueError, "approved upgrade differs"):
                    self.compare([before], [after], [approval])

    def test_earlier_answer_bytes_cannot_change(self):
        with self.assertRaisesRegex(ValueError, "replay differs"):
            self.compare([self.case("before")], [self.case("after", step=b"different")])

    def test_unused_approval_is_rejected(self):
        before, after = self.case("before"), self.case("after")
        with self.assertRaisesRegex(ValueError, "unused upgrade"):
            self.compare([before], [after], [self.approval(before, after)])

    def test_pins_request_and_original_refusal(self):
        before, after = self.case("before", "refused"), self.case("after")
        approval = self.approval(before, after)
        for field in ("request_sha256", "refusal_sha256"):
            with self.subTest(field=field):
                wrong = copy.deepcopy(approval)
                wrong[field] = "0" * 64
                with self.assertRaisesRegex(ValueError, "approved upgrade differs"):
                    self.compare([before], [after], [wrong])


    def correction(self, before, after, fields=()):
        old, new = [Path(case["path"]) for case in (before, after)]
        approval = self.approval(before, after)
        return {"index": 0, "class": "invalid-to-valid-step",
                "request_sha256": approval["request_sha256"],
                "before_sha256": approval["refusal_sha256"], "after_sha256": approval["answer_sha256"],
                "manifest_fields": list(fields),
                "sources": [{"pr": "https://github.com/BohnBawerick/catallax/pull/123", "commit": "e" * 40}],
                "evidence": {"sha256": "4" * 64, "before_step_sha256": release.sha256(old / "answer.step"),
                             "after_step_sha256": release.sha256(new / "answer.step")}}

    def corrected_pair(self):
        # These bytes exercise approval binding, not a real geometry-validity claim.
        return (self.case("before", step=b"invalid synthetic STEP", step_fields={}),
                self.case("after", step=b"corrected synthetic STEP", step_fields={"added_edges": 1}))

    def test_reviewed_step_correction_does_not_claim_strict_equality(self):
        before, after = self.corrected_pair()
        approval = self.correction(before, after, ("step.added_edges",))
        result = self.compare([before], [after], corrections=[approval])
        self.assertEqual(result, {"cases": 1, "equal": False, "failures": [], "upgrades": 0, "corrections": 1})
        self.assertFalse(self.driver.compare(self.report("left", [before]),
                                             self.report("right", [after]))["equal"])

    def test_byte_only_correction_can_be_enumerated(self):
        before = self.case("before", step=b"invalid synthetic STEP", step_fields={})
        after = self.case("after", step=b"corrected synthetic STEP", step_fields={})
        approval = self.correction(before, after)
        self.assertEqual(self.compare([before], [after], corrections=[approval])["corrections"], 1)

    def test_reviewed_fingerprints_cannot_hide_unrelated_manifest_or_output_changes(self):
        for field in ("result", "properties", "marks", "display.obj", "print.stl", "answer.igs"):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as directory:
                previous_root, self.root = self.root, Path(directory)
                try:
                    before, after = self.corrected_pair()
                    if field.endswith((".obj", ".stl", ".igs")):
                        (Path(after["path"]) / field).write_bytes(b"unrelated output")
                    else:
                        path = Path(after["path"]) / "manifest.json"
                        manifest = release.read_json(path)
                        manifest[field] = 2
                        release.write_json(path, manifest)
                    approval = self.correction(before, after, ("step.added_edges",))
                    with self.assertRaisesRegex(ValueError, "unapproved manifest scope|another output"):
                        self.compare([before], [after], corrections=[approval])
                finally:
                    self.root = previous_root

    def test_correction_pins_request_answers_scope_and_evidence_subjects(self):
        before, after = self.corrected_pair()
        approval = self.correction(before, after, ("step.added_edges",))
        for field in ("request_sha256", "before_sha256", "after_sha256", "manifest_fields",
                      "before_step_sha256", "after_step_sha256"):
            with self.subTest(field=field):
                changed = copy.deepcopy(approval)
                if field in changed["evidence"]:
                    changed["evidence"][field] = "0" * 64
                else:
                    changed[field] = [] if field == "manifest_fields" else "0" * 64
                with self.assertRaisesRegex(ValueError, "correction differs|manifest scope|STEP evidence"):
                    self.compare([before], [after], corrections=[changed])

    def test_correction_cannot_authorize_changed_requests_or_statuses(self):
        before, after = self.corrected_pair()
        approval = self.correction(before, after, ("step.added_edges",))
        request_path = Path(after["path"]) / "request.json"
        request = release.read_json(request_path)
        request["max_work"] = 2_000_000
        release.write_json(request_path, request)
        with self.assertRaisesRegex(ValueError, "replay differs"):
            self.compare([before], [after], corrections=[approval])
        old_request = release.read_json(Path(before["path"]) / "request.json")
        release.write_json(request_path, old_request)
        manifest_path = Path(after["path"]) / "manifest.json"
        manifest = release.read_json(manifest_path)
        manifest["status"] = "refused"
        release.write_json(manifest_path, manifest)
        with self.assertRaisesRegex(ValueError, "unsupported status transition"):
            self.compare([before], [after], corrections=[approval])

    def test_correction_requires_both_step_files_and_all_listed_rows_used(self):
        before, after = self.corrected_pair()
        approval = self.correction(before, after, ("step.added_edges",))
        with self.assertRaisesRegex(ValueError, "unused correction"):
            self.compare([before], [before], corrections=[approval])
        (Path(after["path"]) / "answer.step").unlink()
        approval["after_sha256"] = release.answer_digest(self.driver, after["path"])
        with self.assertRaisesRegex(ValueError, "both STEP files"):
            self.compare([before], [after], corrections=[approval])

    def test_reviewed_correction_may_change_approximated_curve_count(self):
        before = self.case("before", step=b"invalid synthetic STEP", step_fields={"approximated_curves": 0})
        after = self.case("after", step=b"corrected synthetic STEP", step_fields={"approximated_curves": 2})
        approval = self.correction(before, after, ("step.approximated_curves",))
        self.assertEqual(self.compare([before], [after], corrections=[approval])["corrections"], 1)

    def refusal(self, name, code="numbers.operation-budget", stage="carrier", detail="ran out"):
        case = self.case(name, status="refused")
        path = Path(case["path"]) / "manifest.json"
        manifest = release.read_json(path)
        manifest["refusal"] = {"code": code, "step": "f2", "stage": stage, "detail": detail}
        release.write_json(path, manifest)
        return case

    def refusal_move(self, before, after, fields):
        approval = self.approval(before, after)
        return {"index": 0, "class": "same-code-refusal-stage",
                "request_sha256": approval["request_sha256"],
                "before_sha256": approval["refusal_sha256"], "after_sha256": approval["answer_sha256"],
                "manifest_fields": list(fields),
                "sources": [{"pr": "https://github.com/BohnBawerick/catallax/pull/162", "commit": "e" * 40}]}

    def test_same_code_refusal_may_move_to_a_later_stage(self):
        before = self.refusal("before")
        after = self.refusal("after", stage="membership", detail="ran out later")
        approval = self.refusal_move(before, after, ("refusal.stage", "refusal.detail"))
        result = self.compare([before], [after], corrections=[approval])
        self.assertEqual(result, {"cases": 1, "equal": False, "failures": [], "upgrades": 0, "corrections": 1})
        with self.assertRaisesRegex(ValueError, "unlisted refusal change"):
            self.compare([before], [after])
        for field in ("request_sha256", "before_sha256", "after_sha256", "manifest_fields"):
            with self.subTest(field=field):
                changed = copy.deepcopy(approval)
                changed[field] = ["refusal.stage"] if field == "manifest_fields" else "0" * 64
                with self.assertRaisesRegex(ValueError, "correction differs|manifest scope"):
                    self.compare([before], [after], corrections=[changed])

    def test_refusal_stage_class_cannot_hide_a_new_code_or_other_changes(self):
        before = self.refusal("before")
        after = self.refusal("after", code="boolean.self-touching-shell", stage="membership")
        approval = self.refusal_move(before, after, ("refusal.stage",))
        with self.assertRaisesRegex(ValueError, "refusal code changed"):
            self.compare([before], [after], corrections=[approval])
        for change in ("marks", "answer.igs"):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as directory:
                previous_root, self.root = self.root, Path(directory)
                try:
                    before = self.refusal("before")
                    after = self.refusal("after", stage="membership")
                    if change == "answer.igs":
                        (Path(after["path"]) / change).write_bytes(b"unrelated output")
                    else:
                        path = Path(after["path"]) / "manifest.json"
                        manifest = release.read_json(path)
                        manifest[change] = ["new"]
                        release.write_json(path, manifest)
                    approval = self.refusal_move(before, after, ("refusal.stage",))
                    with self.assertRaisesRegex(ValueError, "unapproved manifest scope|another output"):
                        self.compare([before], [after], corrections=[approval])
                finally:
                    self.root = previous_root

    def test_correction_classes_follow_the_status_transition(self):
        before, after = self.corrected_pair()
        step_row = self.correction(before, after, ("step.added_edges",))
        wrong = {key: value for key, value in step_row.items() if key != "evidence"}
        wrong.update(**{"class": "same-code-refusal-stage", "manifest_fields": ["refusal.stage"]})
        with self.assertRaisesRegex(ValueError, "unsupported correction class"):
            self.compare([before], [after], corrections=[wrong])
        old, new = self.refusal("old"), self.refusal("new", stage="membership")
        approval = self.refusal_move(old, new, ("refusal.stage",))
        as_step = approval | {"class": "invalid-to-valid-step", "manifest_fields": [],
                              "evidence": step_row["evidence"]}
        with self.assertRaisesRegex(ValueError, "unsupported correction class"):
            self.compare([old], [new], corrections=[as_step])


    def test_reviewed_refusal_change_binds_a_new_code(self):
        before = self.refusal("before")
        after = self.refusal("after", code="boolean.self-touching-shell", stage="membership")
        path = Path(after["path"]) / "manifest.json"
        manifest = release.read_json(path)
        manifest["question_id"] = "read later"
        release.write_json(path, manifest)
        fields = ("question_id", "refusal.code", "refusal.stage")
        approval = self.refusal_move(before, after, fields) | {"class": "reviewed-refusal-change"}
        self.assertEqual(self.compare([before], [after], corrections=[approval])["corrections"], 1)
        with self.assertRaisesRegex(ValueError, "manifest scope"):
            self.compare([before], [after], corrections=[approval | {"manifest_fields": fields[1:]}])
        with self.assertRaisesRegex(ValueError, "refusal code changed"):
            self.compare([before], [after], corrections=[approval | {"class": "same-code-refusal-stage"}])
        same = self.refusal("same", stage="membership")
        reviewed = self.refusal_move(before, same, ("refusal.stage",)) | {"class": "reviewed-refusal-change"}
        with self.assertRaisesRegex(ValueError, "keeps its code"):
            self.compare([before], [same], corrections=[reviewed])

    def test_reviewed_refusal_change_cannot_hide_other_changes(self):
        before = self.refusal("before")
        after = self.refusal("after", code="boolean.self-touching-shell")
        path = Path(after["path"]) / "manifest.json"
        manifest = release.read_json(path)
        manifest["marks"] = ["new"]
        release.write_json(path, manifest)
        approval = self.refusal_move(before, after, ("refusal.code",)) | {"class": "reviewed-refusal-change"}
        with self.assertRaisesRegex(ValueError, "unapproved manifest scope"):
            self.compare([before], [after], corrections=[approval])


class RebuildRecordTests(unittest.TestCase):
    """The wider replay record: corpus rows, recorded answers and requalification."""

    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location(
            "release_test_driver", release.ROOT / ".claude/skills/verify-catallax/scripts/control-catallax.py")
        cls.driver = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.driver)

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.corpus = self.root / "rebuild-corpus.jsonl"
        self.answers = self.root / "rebuild-answers.json"
        self.recipe = {"format": "catallax-recipe/3", "name": "block",
                       "steps": [{"id": "base", "op": "block", "size": ["4", "4", "4"]}],
                       "result": "base"}
        self.rows = [
            {"format": "catallax-question/1", "source": "recipe",
             "problem": self.recipe, "max_work": 100000000},
            {"format": "catallax-question/1", "source": "recipe",
             "problem": dict(self.recipe, name="second"), "max_work": 50},
        ]
        self.corpus.write_text("".join(self.driver.canonical(row) + "\n" for row in self.rows))
        release.REBUILD_CORPUS = self.corpus
        release.REBUILD_ANSWERS = self.answers
        self.addCleanup(self.restore)
        self.original_committed = release.committed_bytes
        release.committed_bytes = lambda path: path.read_bytes()
        self.addCleanup(setattr, release, "committed_bytes", self.original_committed)
        self.runs = self.root / "runs"
        self.serial = itertools.count()

    def restore(self):
        release.REBUILD_CORPUS = release.ROOT / "tools/release/rebuild-corpus.jsonl"
        release.REBUILD_ANSWERS = release.ROOT / "tools/release/rebuild-answers.json"

    def case(self, name, status="answered", value=2, max_work=100000000, step=b"exact STEP bytes"):
        directory = self.runs / name
        directory.mkdir(parents=True)
        manifest = {"status": status, "result": ({"faces": value} if status == "answered" else None),
                    "refusal": {"code": "intersect.unsupported-pair"} if status == "refused" else None,
                    "kernel": {"commit": "a" * 40}, "stats": {"wall_ms": len(name)}}
        release.write_json(directory / "manifest.json", manifest)
        request = {"format": "catallax-request/1", "max_work": max_work,
                   "source": "recipe", "problem": self.recipe, "out": str(directory)}
        release.write_json(directory / "request.json", request)
        if status == "answered":
            (directory / "answer.step").write_bytes(step)
        files = {p.name: self.driver.file_hash(p) for p in directory.iterdir() if p.is_file()}
        release.write_json(directory / "artifact.json", {"key": "key-" + name, "files": files})
        return {"key": "key-" + name, "id": name, "path": str(directory), "cache_hit": False,
                "response": {"work": {"used": 1, "limit": max_work, "exhausted": False},
                             "manifest": manifest}}

    def report(self, name, cases, workers=4, binary=None):
        # Real runs keep run.json and config.json in one owned directory.
        directory = self.runs / name
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / "run.json"
        release.write_json(path, {
            "format": "catallax-run/1",
            "config": {"binary_sha256": binary or hashlib.sha256(name.encode()).hexdigest(), "workers": workers},
            "execution": {"fresh": True, "workers": workers,
                          "worker_exits": {str(i): 0 for i in range(workers)}},
            "elapsed_ms": 1, "cases": cases})
        release.write_json(path.parent / "config.json", {"workers": workers})
        return path

    def rehash(self, case):
        """Recompute the artifact hash after editing files in a synthetic case."""
        directory = Path(case["path"])
        # The driver hashes the artifact before writing artifact.json itself.
        (directory / "artifact.json").unlink()
        files = {p.name: self.driver.file_hash(p) for p in directory.iterdir() if p.is_file()}
        release.write_json(directory / "artifact.json", {"key": case["key"], "files": files})

    def run_case(self, index, tag="run", status="answered", value=2, step=b"exact STEP bytes", cache=False):
        row = self.rows[index]
        case = self.case(f"case{index}-{tag}-{next(self.serial)}", status, value,
                         max_work=row["max_work"], step=step)
        directory = self.runs / case["id"]
        request = release.read_json(directory / "request.json")
        request["problem"] = row["problem"]
        release.write_json(directory / "request.json", request)
        self.rehash(case)
        case["cache_hit"] = cache
        return case

    def record(self):
        return release.read_json(self.answers)

    def test_record_then_check_reproduces_the_run(self):
        path = self.report("w4", [self.run_case(0), self.run_case(1, value=1)])
        args = type("Args", (), {"run": str(path), "workers": 4, "out": str(self.answers)})()
        release.rebuild_record(args)
        document = self.record()
        self.assertEqual(document["format"], "catallax-rebuild-answers/1")
        self.assertEqual(document["corpus_sha256"], release.sha256(self.corpus))
        self.assertEqual([row["index"] for row in document["cases"]], [0, 1])
        self.assertEqual(document["cases"][0]["faces"], 2)
        self.assertEqual(document["cases"][1]["faces"], 1)
        for row in document["cases"]:
            self.assertIsNone(row["step_faces"])
            self.assertIsNone(row["names"])
            self.assertNotEqual(row["answer_sha256"], row["manifest_sha256"])
        check = type("Args", (), {"run": str(path), "workers": 4, "allow_uncommitted": True})()
        release.check_rebuild_record(check)

    def test_check_rejects_a_changed_answer(self):
        path = self.report("w4", [self.run_case(0), self.run_case(1, value=1)])
        args = type("Args", (), {"run": str(path), "workers": 4, "out": str(self.answers)})()
        release.rebuild_record(args)
        changed = self.report("w4b", [self.run_case(0), self.run_case(1, step=b"changed bytes")])
        with self.assertRaisesRegex(ValueError, "answer_sha256 differs from the record"):
            release.check_rebuild_record(
                type("Args", (), {"run": str(changed), "workers": 4, "allow_uncommitted": True})())

    def test_check_rejects_a_changed_request_or_allowance(self):
        for field, value in (("max_work", 40), ("problem", {"other": 1})):
            with self.subTest(field=field):
                path = self.report("w4", [self.run_case(0), self.run_case(1, value=1)])
                args = type("Args", (), {"run": str(path), "workers": 4, "out": str(self.answers)})()
                release.rebuild_record(args)
                case = self.run_case(1, value=1)
                request = release.read_json(Path(case["path"]) / "request.json")
                request[field] = value
                release.write_json(Path(case["path"]) / "request.json", request)
                self.rehash(case)
                changed = self.report("w4b", [self.run_case(0), case])
                with self.assertRaisesRegex(ValueError, "request differs from the corpus row"):
                    release.check_rebuild_record(
                        type("Args", (), {"run": str(changed), "workers": 4, "allow_uncommitted": True})())

    def test_check_rejects_cached_internal_error_and_extra_work(self):
        path = self.report("w4", [self.run_case(0), self.run_case(1, value=1)])
        args = type("Args", (), {"run": str(path), "workers": 4, "out": str(self.answers)})()
        release.rebuild_record(args)
        cached = self.report("w4c", [self.run_case(0), self.run_case(1, value=1, cache=True)])
        with self.assertRaisesRegex(ValueError, "cached rebuild"):
            release.check_rebuild_record(
                type("Args", (), {"run": str(cached), "workers": 4, "allow_uncommitted": True})())
        case = self.run_case(1, value=1)
        case["response"]["work"] = {"used": 60, "limit": 50, "exhausted": True}
        exhausted = self.report("w4d", [self.run_case(0), case])
        with self.assertRaisesRegex(ValueError, "work outside the request allowance"):
            release.check_rebuild_record(
                type("Args", (), {"run": str(exhausted), "workers": 4, "allow_uncommitted": True})())

    def test_record_rejects_a_run_that_does_not_match_the_corpus(self):
        path = self.report("w4", [self.run_case(0)])
        with self.assertRaisesRegex(ValueError, "does not match the committed corpus"):
            release.rebuild_record(type("Args", (), {"run": str(path), "workers": 4,
                                                     "out": str(self.answers)})())

    def candidate_runs(self, baseline_refusal=False, candidate_change=False):
        """One baseline and three candidate runs over the same two cases."""
        paths = []
        candidate = hashlib.sha256(b"candidate").hexdigest()
        for name, workers in zip(("base", "one", "two", "four"), (1, 1, 2, 4)):
            if name == "base" and baseline_refusal:
                cases = [self.run_case(0, tag=name), self.run_case(1, tag=name, status="refused")]
            elif name == "two" and candidate_change:
                cases = [self.run_case(0, tag=name), self.run_case(1, tag=name, value=3)]
            else:
                cases = [self.run_case(0, tag=name), self.run_case(1, tag=name, value=1)]
            binary = candidate if name != "base" else None
            paths.append(self.report(name, cases, workers=workers, binary=binary))
        return paths

    def test_summary_validates_candidate_equality_and_binds_the_files(self):
        paths = self.candidate_runs(baseline_refusal=True)
        record = type("Args", (), {"run": str(paths[1]), "workers": 1, "out": str(self.answers)})()
        release.rebuild_record(record)
        summary_path = self.root / "summary.json"
        args = type("Args", (), {"runs": [str(path) for path in paths], "out": str(summary_path)})()
        release.rebuild_summary(args)
        summary = release.read_json(summary_path)
        self.assertEqual(summary["format"], release.REBUILD_FORMAT)
        self.assertEqual(summary["corpus_sha256"], release.sha256(self.corpus))
        self.assertEqual(summary["answers_sha256"], release.sha256(self.answers))
        self.assertEqual(summary["changes"], [{"index": 1, "cause": "manifest"}])
        self.assertFalse(summary["comparisons"]["baseline-w1 vs candidate-w1"]["equal"])
        for name in release.COMPARISONS[1:]:
            self.assertTrue(summary["comparisons"][name]["equal"])
        check = type("Args", (), {"summary": str(summary_path), "commit": "a" * 40,
                                  "allow_historical": False})()
        release.check_rebuild(check)

    def test_summary_rejects_a_worker_difference(self):
        paths = self.candidate_runs(candidate_change=True)
        record = type("Args", (), {"run": str(paths[1]), "workers": 1, "out": str(self.answers)})()
        release.rebuild_record(record)
        with self.assertRaisesRegex(ValueError, "differs from the record"):
            release.rebuild_summary(type("Args", (), {"runs": [str(path) for path in paths],
                                                      "out": str(self.root / "summary.json")})())
        self.assertFalse((self.root / "summary.json").exists())

    def test_check_rejects_a_wrong_commit_or_tampered_summary(self):
        paths = self.candidate_runs()
        record = type("Args", (), {"run": str(paths[1]), "workers": 1, "out": str(self.answers)})()
        release.rebuild_record(record)
        summary_path = self.root / "summary.json"
        release.rebuild_summary(type("Args", (), {"runs": [str(path) for path in paths],
                                                  "out": str(summary_path)})())
        check = type("Args", (), {"summary": str(summary_path), "commit": "other",
                                  "allow_historical": False})()
        with self.assertRaisesRegex(ValueError, "not from the selected commit"):
            release.check_rebuild(check)
        tampered = release.read_json(summary_path)
        tampered["runs"]["candidate-w1"]["binary_sha256"] = "0" * 64
        path = self.root / "tampered.json"
        release.write_json(path, tampered)
        check = type("Args", (), {"summary": str(path), "commit": release.read_json(summary_path)["candidate_sha"],
                                  "allow_historical": False})()
        with self.assertRaisesRegex(ValueError, "different candidate binaries"):
            release.check_rebuild(check)


if __name__ == "__main__":
    unittest.main()
