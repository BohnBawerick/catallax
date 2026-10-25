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

    def case(self, name, status="answered", step=b"exact STEP bytes", value=1):
        directory = self.root / name
        directory.mkdir()
        manifest = {"status": status, "result": value if status == "answered" else None,
                    "refusal": {"code": "numbers.operation-budget"} if status == "refused" else None,
                    "kernel": {"commit": name}, "stats": {"wall_ms": len(name)}}
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
                          for name in ("answer.step", "display.obj", "print.stl")}})

        return {"index": index, "request_sha256": self.driver.digest(request),
                "refusal_sha256": answer(before), "answer_sha256": answer(after)}

    def compare(self, before, after, upgrades=()):
        return release.compare_baseline(self.driver, self.report("base", before),
                                        self.report("candidate", after), list(upgrades))

    def test_preserves_earlier_answers_and_unchanged_refusals(self):
        result = self.compare([self.case("a"), self.case("b", "refused")],
                              [self.case("c"), self.case("d", "refused")])
        self.assertEqual(result, {"cases": 2, "equal": True, "failures": [], "upgrades": 0})

    def test_listed_refusal_upgrade_pins_the_new_answer(self):
        before, after = self.case("before", "refused"), self.case("after")
        result = self.compare([before], [after], [self.approval(before, after)])
        self.assertEqual(result, {"cases": 1, "equal": False, "failures": [], "upgrades": 1})

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
        for field in ("manifest.json", "answer.step", "display.obj", "print.stl"):
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


if __name__ == "__main__":
    unittest.main()
