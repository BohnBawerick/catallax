"""A missing export verdict must not erase a completed opening verdict."""
import unittest

from opened_import_census import comparison_outcomes


class StageComparisonTests(unittest.TestCase):
    def test_unmeasured_stages_and_measured_inventory_losses_stay_separate(self):
        rows = {
            "export_timeout": dict(outcome="full", observed_body_count=1, export_validity="timeout"),
            "export_oom": dict(outcome="full", observed_body_count=1, export_validity="out_of_memory"),
            "open_oom": dict(outcome="out_of_memory"),
            "missing_body": dict(outcome="full", observed_body_count=1, export_validity="timeout"),
            "surface": dict(outcome="surface", export_validity="source_surfaces_without_checked_brep"),
        }
        expected_full = {key: 1 for key in rows}
        expected_full["missing_body"] = 2
        result = comparison_outcomes(rows, expected_full, {"export_oom", "surface"})
        self.assertEqual(result["tier_losses"], ["missing_body", "surface"])
        self.assertEqual(result["frozen_valid_export_losses"], ["surface"])
        self.assertEqual(result["baseline_not_measured"], ["export_oom", "open_oom"])
        self.assertEqual(result["not_measured"], ["export_oom", "export_timeout", "missing_body", "open_oom"])


if __name__ == "__main__":
    unittest.main()
