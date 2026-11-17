"""A missing export verdict must not erase a completed opening verdict."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

from opened_import_census import comparison_outcomes, observe, export_inventory, document_export_refusal


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
        result = comparison_outcomes(rows, expected_full, {"export_oom", "surface"}, {})
        self.assertEqual(result["tier_losses"], ["missing_body", "surface"])
        self.assertEqual(result["frozen_valid_export_losses"], ["surface"])
        self.assertEqual(result["baseline_not_measured"], ["export_oom", "open_oom"])
        self.assertEqual(result["not_measured"], ["export_oom", "export_timeout", "missing_body", "open_oom"])

    def test_a_producer_crash_after_opening_preserves_its_completed_solid_verdict(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "source.step"
            source.write_bytes(b"immutable source")
            expected = hashlib.sha256(source.read_bytes()).hexdigest()
            checkpoint = dict(
                stage="open", outcome="opened", source_sha256=expected,
                work_used=64, work_exhausted=False, read_work_used=16,
                file_admission_work_used=32, first_refusal=None,
                bodies=[dict(tier="valid_solid", proof="certified", healing_issues=[])],
            )
            producer = directory / "crashing-producer"
            producer.write_text(
                f"#!{sys.executable}\nimport os\n"
                f"print({json.dumps(checkpoint)!r}, flush=True)\nos._exit(23)\n"
            )
            producer.chmod(0o700)
            baseline = dict(
                model_id="source", source_sha256=expected,
                baseline_outcome="full", baseline_first_refusal="",
                outcome="full", after_fix_cause="", coverage_complete="True",
            )
            row = observe(producer, source, baseline, directory / "capture", False)
            self.assertEqual(row["outcome"], "full")
            self.assertEqual(row["observed_body_count"], 1)
            self.assertEqual(row["export_validity"], "process_error")
            self.assertIsNone(row["has_valid_exported_solid"])
            self.assertTrue(row["any_valid_exported_solid_count_is_lower_bound"])
            comparison = comparison_outcomes({"source": row}, {"source": 1}, {"source"}, {})
            self.assertEqual(comparison["tier_losses"], [])
            self.assertEqual(comparison["frozen_valid_export_losses"], [])
            self.assertEqual(comparison["baseline_not_measured"], ["source"])

    def test_checked_loss_and_tolerant_gain_are_aggregate_observations(self):
        current = {
            key: dict(
                outcome="full", observed_body_count=1, all_solids_export_valid=True,
                valid_solid_body_count=0, repaired_solid_body_count=0,
                tolerant_solid_body_count=1, surface_body_count=0,
            ) for key in ["demoted", "improved", "unknown"]
        }
        prior = {
            "demoted": dict(
                observed_body_count="1", valid_solid_body_count="0",
                repaired_solid_body_count="1", surface_body_count="0",
            ),
            "improved": dict(
                observed_body_count="1", valid_solid_body_count="0",
                repaired_solid_body_count="0", surface_body_count="1",
            ),
            "unknown": dict(outcome="full", bodies="1"),
        }
        current["swapped"] = dict(
            outcome="full", observed_body_count=2, valid_solid_body_count=0,
            repaired_solid_body_count=0, tolerant_solid_body_count=1, surface_body_count=1,
        )
        prior["swapped"] = dict(
            observed_body_count=2, valid_solid_body_count=0,
            repaired_solid_body_count=1, surface_body_count=1,
        )
        result = comparison_outcomes(current, {key: row["observed_body_count"] for key, row in current.items()}, set(), prior)
        self.assertEqual(result["solid_coverage_losses"], [])
        self.assertEqual(result["tier_losses"], ["demoted", "swapped"])
        self.assertEqual(result["checked_loss_tolerant_gain_files"], 2)
        self.assertEqual(result["checked_loss_tolerant_gain"]["swapped"],
                         dict(checked_loss=1, tolerant_gain=1))
        self.assertEqual(
            result["per_file_tier_transitions"]["swapped"]["tier_transition"],
            "checked_loss_with_tolerant_gain",
        )
        self.assertEqual(result["tier_comparison_unknown"], ["unknown"])
        self.assertEqual(
            result["per_file_tier_transitions"]["unknown"]["tier_transition"],
            "unknown_prior_tiers",
        )

    def test_a_valid_native_sibling_counts_without_claiming_complete_coverage(self):
        files = [
            dict(kind="native_document", bodies=[0, 1], outcome="invalid_after_builder", valid_solid_count=1),
            dict(kind="surface", bodies=[2], outcome="valid_after_builder", valid_solid_count=0),
        ]
        inventory = export_inventory(files, 3, False, True)
        self.assertTrue(inventory["has_valid_exported_solid"])
        self.assertEqual(inventory["any_valid_exported_solid_count"], 1)
        self.assertEqual(inventory["valid_exported_solid_count"], 0)
        self.assertEqual(inventory["invalid_exported_body_count"], 2)

        files[0].update(outcome="valid_after_builder", valid_solid_count=3)
        inventory = export_inventory(files, 3, False, True)
        self.assertEqual(inventory["any_valid_exported_solid_count"], 3)
        self.assertEqual(inventory["valid_exported_solid_count"], 2)


    def test_valid_prototypes_do_not_establish_source_occurrence_coverage(self):
        files = [dict(kind="solid", bodies=[0], outcome="valid_after_builder", valid_solid_count=1)]
        inventory = export_inventory(files, 1, True, True)
        self.assertTrue(inventory["prototype_solids_export_valid"])
        self.assertFalse(inventory["all_solids_export_valid"])
        self.assertFalse(inventory["source_document_exported"])
        row = dict(inventory, outcome="full", observed_body_count=1)
        comparison = comparison_outcomes({"source": row}, {"source": 1}, {"source"}, {})
        self.assertEqual(comparison["frozen_valid_export_losses"], [])
        files[0]["kind"] = "native_document"
        self.assertFalse(export_inventory(files, 1, True, True)["all_solids_export_valid"])
        files[0]["occurrence_graph_complete"] = True
        self.assertTrue(export_inventory(files, 1, True, True)["all_solids_export_valid"])

    def test_a_detached_valid_body_is_not_complete_source_coverage(self):
        files = [dict(kind="selected_document", bodies=[0, 1], outcome="valid_after_builder",
                      valid_solid_count=2, occurrence_graph_complete=False)]
        inventory = export_inventory(files, 2, True, True)
        self.assertTrue(inventory["prototype_solids_export_valid"])
        self.assertTrue(inventory["has_valid_exported_solid"])
        self.assertFalse(inventory["all_solids_export_valid"])
        self.assertFalse(inventory["source_document_exported"])

    def test_stopped_sibling_export_is_an_observed_lower_bound(self):
        files = [dict(kind="native_body", bodies=[0], outcome="invalid_after_builder", valid_solid_count=0)]
        inventory = export_inventory(files, 2, True, False)
        self.assertIsNone(inventory["has_valid_exported_solid"])
        self.assertEqual(inventory["any_valid_exported_solid_count"], 0)
        self.assertTrue(inventory["any_valid_exported_solid_count_is_lower_bound"])
        self.assertFalse(inventory["export_observation_complete"])
        files[0].update(outcome="valid_after_builder", valid_solid_count=1)
        self.assertTrue(export_inventory(files, 2, True, False)["has_valid_exported_solid"])

    def test_a_refused_part_leaves_other_native_parts_unmeasured(self):
        for root in [
            dict(outcome="reader_refused", parts_expected=2, parts_inspected=0, observation_complete=False),
            dict(outcome="reader_refused"),
        ]:
            files = [dict(kind="native_document", bodies=[0], outcome="invalid_after_builder",
                          valid_solid_count=0, raw_source_roots=[root])]
            inventory = export_inventory(files, 1, True, True)
            self.assertIsNone(inventory["has_valid_exported_solid"])
            self.assertFalse(inventory["export_observation_complete"])
            self.assertTrue(inventory["any_valid_exported_solid_count_is_lower_bound"])
            row = dict(inventory, outcome="full", observed_body_count=1,
                       export_validity="raw_observation_incomplete")
            comparison = comparison_outcomes({"source": row}, {"source": 1}, {"source"}, {})
            self.assertEqual(comparison["not_measured"], ["source"])
            self.assertEqual(comparison["baseline_not_measured"], ["source"])
            self.assertEqual(comparison["frozen_valid_export_losses"], [])
            files[0]["valid_solid_count"] = 1
            inventory = export_inventory(files, 1, True, True)
            self.assertTrue(inventory["has_valid_exported_solid"])
            self.assertEqual(inventory["any_valid_exported_solid_count"], 1)
            self.assertTrue(inventory["any_valid_exported_solid_count_is_lower_bound"])

    def test_document_refusals_are_counted_before_any_sibling_fallback(self):
        response = dict(stage="export", failed_body=None, files=[], open=dict(bodies=[{}]),
                        work_exhausted=False, refusal=dict(code="import.unsupported", detail="refused"))
        self.assertEqual(document_export_refusal(response), "import.unsupported")
        response["failed_body"] = 0
        self.assertIsNone(document_export_refusal(response))
        response["failed_body"] = None
        response["work_exhausted"] = True
        self.assertIsNone(document_export_refusal(response))
        response["work_exhausted"] = False
        response["refusal"]["code"] = "numbers.operation-budget"
        self.assertIsNone(document_export_refusal(response))

if __name__ == "__main__":
    unittest.main()
