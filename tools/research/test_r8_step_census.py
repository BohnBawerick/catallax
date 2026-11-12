"""Selection and denominator checks for the R8 collector."""

import unittest
from pathlib import Path
import sys
import tempfile
import subprocess
from unittest.mock import patch

from r8_step_census import choose, observe, qualify_coverage, restrict, summarize


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

    def test_filter_preserves_seeded_order_and_rejects_unknown_or_duplicate_ids(self):
        entries = [f"{i:08d}/{i:08d}_abc_step_000.step" for i in range(100)]
        selected, _ = choose(entries, 30, 20261111)
        subset = selected[::3]
        self.assertEqual(restrict(selected, reversed([i for i, _ in subset])), subset)
        for ids in [[], ["99999999"], [subset[0][0]] * 2, ["../00000000"]]:
            with self.assertRaises(ValueError):
                restrict(selected, ids)

    def test_unknown_coverage_and_timeouts_stay_in_denominators(self):
        rows = [
            {"coverage_complete": True, "face_kinds": {"plane": 6},
             "data_section_scan_complete": True, "data_sections": 1,
             "unsupported_data_sections": 0,
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
            def timeout(command, *, stdout, **kwargs):
                stdout.write(b'{"stage":"faces","face_kinds":{"plane":6},'
                             b'"coverage_complete":true}\n')
                stdout.flush()
                raise subprocess.TimeoutExpired(command, kwargs["timeout"])

            with patch("r8_step_census.subprocess.run", side_effect=timeout):
                row = observe(Path(sys.executable), root / "model.step",
                              "00000000", root, 0.2)
            self.assertEqual(row["outcome"], "observation_timeout")
            self.assertEqual(row["face_kinds"], {"plane": 6})

    def test_legacy_or_multiple_sections_cannot_establish_absence(self):
        original = {"coverage_complete": True, "face_kinds": {"plane": 6},
                    "all_face_kinds_buildable": True, "needs_spline_faces": False,
                    "outcome": "refused", "refusal_kind": "unknown_unit"}
        for extra in [{}, {"data_section_scan_complete": True, "data_sections": 2,
                           "unsupported_data_sections": 0},
                      {"data_section_scan_complete": True, "data_sections": 1,
                       "unsupported_data_sections": 1}]:
            row = qualify_coverage({**original, **extra})
            self.assertFalse(row["coverage_complete"])
            self.assertIsNone(row["all_face_kinds_buildable"])
            self.assertIsNone(row["needs_spline_faces"])
            self.assertEqual(row["face_kinds"], {"plane": 6})
        self.assertTrue(original["coverage_complete"])
        summary = summarize([original])
        self.assertEqual(summary["complete_source_traversals"], 0)
        self.assertIsNone(summary["all_face_kinds_buildable_share_of_complete"])
        self.assertEqual(summary["all_face_kinds_buildable"], {"unknown": 1})
        self.assertEqual(summary["needs_spline_faces"], {"unknown": 1})

    def test_known_spline_presence_survives_incomplete_section_coverage(self):
        row = qualify_coverage({"coverage_complete": True,
                                "all_face_kinds_buildable": False,
                                "needs_spline_faces": True})
        self.assertFalse(row["coverage_complete"])
        self.assertIsNone(row["all_face_kinds_buildable"])
        self.assertIs(row["needs_spline_faces"], True)


if __name__ == "__main__":
    unittest.main()
