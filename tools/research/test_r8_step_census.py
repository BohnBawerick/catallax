"""Selection and denominator checks for the R8 collector."""

import unittest
from pathlib import Path
import sys
import tempfile

from r8_step_census import choose, observe, summarize


class CensusTests(unittest.TestCase):
    def test_selection_does_not_depend_on_archive_order(self):
        entries = [f"{i:08d}/{i:08d}_abc_step_000.step" for i in range(100)]
        a, total = choose(entries, 30, 20261111)
        b, _ = choose(list(reversed(entries)), 30, 20261111)
        self.assertEqual(a, b)
        self.assertEqual(total, 100)
        self.assertEqual(len(set(model_id for model_id, _ in a)), 30)
        self.assertNotEqual(a, choose(entries, 30, 1)[0])

    def test_selection_rejects_unsafe_and_ambiguous_inputs(self):
        for entries in [["../00000000_a.step"], ["/00000000_a.step"],
                        ["00000000_a.step", "00000000_b.step"]]:
            with self.assertRaises(ValueError):
                choose(entries, 1, 0)

    def test_unknown_coverage_and_timeouts_stay_in_denominators(self):
        rows = [
            {"coverage_complete": True, "face_kinds": {"plane": 6},
             "all_face_kinds_buildable": True, "needs_spline_faces": False,
             "outcome": "refused", "refusal_kind": "not_solid"},
            {"coverage_complete": False, "face_kinds": {},
             "all_face_kinds_buildable": None, "needs_spline_faces": None,
             "outcome": "observation_timeout"},
        ]
        result = summarize(rows)
        self.assertEqual(result["evaluated_model_ids"], 2)
        self.assertEqual(result["complete_source_traversals"], 1)
        self.assertEqual(result["all_face_kinds_buildable"], {"true": 1, "unknown": 1})
        self.assertEqual(result["all_face_kinds_buildable_share_of_complete"], 1)
        self.assertIsNone(result["fillet_count"])
        self.assertIsNone(result["chamfer_count"])

    def test_timeout_preserves_the_face_checkpoint(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            script = root / "observer.py"
            script.write_text('import time\nprint(\'{"stage":"faces",'
                              '"face_kinds":{"plane":6},"coverage_complete":true}\','
                              'flush=True)\ntime.sleep(10)\n')
            row = observe(Path(sys.executable), script, "00000000", root, 0.2)
            self.assertEqual(row["outcome"], "observation_timeout")
            self.assertEqual(row["face_kinds"], {"plane": 6})


if __name__ == "__main__":
    unittest.main()
