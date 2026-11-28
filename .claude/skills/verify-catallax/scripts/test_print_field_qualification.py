"""Public print and field qualification through the streamed verifier."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest

from test_control import BINARY, load_driver


FIXTURES = Path(__file__).with_name("fixtures")
CORPUS = FIXTURES / "print-field-qualification.jsonl"
EXPECTED = FIXTURES / "print-field-expected.json"


def recipe_name(problem):
    return problem.get("recipe", problem)["name"]


@unittest.skipUnless(BINARY.exists(), "build catallax-verify first")
class PrintFieldQualificationTests(unittest.TestCase):
    def test_print_and_field_corpus(self):
        driver = load_driver()
        inputs = driver.load([CORPUS], "auto")
        expected = json.loads(EXPECTED.read_text())
        names = [recipe_name(item["problem"]) for item in inputs]
        self.assertCountEqual(names, expected.keys())

        with tempfile.TemporaryDirectory(prefix="catallax-print-field-") as temporary:
            root = Path(temporary)
            for workers in (1, 2, 4):
                out = root / str(workers)
                driver.run(argparse.Namespace(
                    inputs=[CORPUS], source="auto", adapter=None, binary=BINARY,
                    out=out, fresh=True, workers=workers, cores=None,
                    cache=None, max_work=100_000_000, dry_run=False))
                report = json.loads((out / "run.json").read_text())
                self.assertEqual(report["execution"]["worker_exits"],
                                 {str(i): 0 for i in range(workers)})
                self.assertEqual(len(report["cases"]), len(inputs))
                for item, name, case in zip(inputs, names, report["cases"]):
                    with self.subTest(workers=workers, case=name):
                        self.check_case(driver, item, case, expected[name])
                if workers > 1:
                    comparison = driver.compare(root / "1/run.json", out / "run.json")
                    self.assertTrue(comparison["equal"], comparison)

    def check_case(self, driver, item, case, expected):
        directory = Path(case["path"])
        request = json.loads((directory / "request.json").read_text())
        request.pop("out")
        self.assertEqual(request, dict(item, format=driver.CONTRACT))
        self.assertEqual(request["source"], expected["source"])
        self.assertFalse(case["cache_hit"])
        self.assertTrue(driver.valid_artifact(directory, case["key"]))
        response = case["response"]
        manifest = json.loads((directory / "manifest.json").read_text())
        self.assertEqual(manifest, response["manifest"])
        self.assertEqual(manifest["status"], expected["status"])
        self.assertEqual((manifest["refusal"] or {}).get("code"), expected["refusal"])
        self.assertEqual(manifest.get("field"), expected.get("field"))

        problem = item["problem"]
        printing_request = problem["format"] == "catallax-print/1"
        operation = "print" if printing_request else "display"
        self.assertEqual(response["work"]["limit"], item["max_work"])
        self.assertFalse(response["work"]["exhausted"])
        operations = json.loads((directory / "operations.json").read_text())
        self.assertEqual(operations, response["operation_log"])
        self.assertEqual(
            [(op["kernel"], op["status"]) for op in operations["operations"]
             if op["operation"] == operation],
            [(expected["route"], expected["status"])],
        )

        files = ("answer.step", "answer.igs", "display.obj", "print.stl")
        self.assertCountEqual([name for name in files if (directory / name).exists()],
                              expected["files"])
        if not printing_request:
            self.assertNotIn("printing", manifest)
            self.assertEqual(manifest["face_measurements"]["display"]["file"], "display.obj")
            return

        printing = manifest["printing"]
        self.assertIs(manifest["exact"], False)
        self.assertEqual(printing["route"], expected["route"])
        self.assertIs(printing["cell_cap_applicable"], expected["route"] == "field")
        tolerance = Fraction(problem["tolerance"])
        self.assertEqual(Fraction(printing["tolerance_mm"]), tolerance)
        if "max_cells" in problem:
            self.assertEqual(printing["max_cells"], problem["max_cells"])
        self.assertCountEqual(manifest.get("files", {}), expected["files"])
        if expected["status"] == "refused":
            self.assertIsNone(printing["file"])
            self.assertIsNone(printing["deviation_bound_mm"])
            return

        self.assertEqual(printing["file"], "print.stl")
        self.assertEqual(printing["completeness"], "proven-two-sided")
        self.assertIs(printing["watertight"], True)
        self.assertGreater(printing["triangles"], 0)
        self.assertGreater(printing["vertices"], 0)
        bound = Fraction(printing["deviation_bound_mm"])
        self.assertGreaterEqual(bound, 0)
        self.assertLessEqual(bound, tolerance)
        self.assertEqual(printing["sha256"], manifest["files"]["print.stl"]["sha256"])
        for name in expected["files"]:
            path = directory / name
            self.assertEqual(manifest["files"][name]["sha256"], driver.file_hash(path))
            self.assertEqual(manifest["files"][name]["bytes"], path.stat().st_size)
