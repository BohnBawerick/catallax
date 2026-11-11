"""Controls for the planar cohort's selection and public evidence."""

import csv
from pathlib import Path
import tempfile
import unittest

from planar_self_intersection import cohort, public_table


class PlanarCohortTests(unittest.TestCase):
    def test_only_the_saved_planar_dispatch_enters_the_cohort(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "one.step").write_text("synthetic control")
            details = {"one": {"intersection_path": "planar_polygon_pass"},
                       "two": {"intersection_path": "curved_or_circular_edge_pass"}}
            self.assertEqual(cohort(details, [["one", "one.step"]], root),
                             {"one": root / "one.step"})

    def test_member_cannot_escape_the_model_directory(self):
        details = {"one": {"intersection_path": "planar_polygon_pass"}}
        for member in ("../one.step", "/one.step"):
            with self.assertRaisesRegex(ValueError, "unsafe archive member"):
                cohort(details, [["one", member]], Path("/unused"))

    def test_public_rows_exclude_source_paths_and_diagnostics(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "public.csv"
            public_table(path, [{"model_id": "one", "classification": "nonincident_boundary",
                                 "private_path": "private source", "diagnostic": "private detail"}])
            with path.open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(rows[0]["classification"], "nonincident_boundary")
            self.assertNotIn("private", path.read_text())
            self.assertNotIn(b"\r", path.read_bytes())


if __name__ == "__main__":
    unittest.main()
