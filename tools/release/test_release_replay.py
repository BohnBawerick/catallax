"""Replay policy applies to saved manifests, inputs and exact output bytes."""
import copy
import importlib.util
from pathlib import Path
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

if __name__ == "__main__":
    unittest.main()
