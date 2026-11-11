"""Checks that triage preserves exact causes and never invents isolated gains."""

import unittest

from r8_refusal_triage import render, summarize


class TriageTests(unittest.TestCase):
    def test_exact_unsupported_variants_remain_separate(self):
        rows = [{"model_id": str(i), "outcome": "refused", "refusal_code": "unsupported"} for i in range(3)]
        errors = {"0": {"code": "unsupported", "type_name": "CONICAL_SURFACE (axis length is not rational)"},
                  "1": {"code": "unsupported", "type_name": "CONICAL_SURFACE (bound is not one complete rim circle)"},
                  "2": {"code": "unsupported", "type_name": "CONICAL_SURFACE (axis length is not rational)"}}
        result = summarize(rows, errors, {})
        self.assertEqual([x["models"] for x in result["unsupported"]], [2, 1])
        self.assertTrue(all(x["entity_type"] == "CONICAL_SURFACE" for x in result["unsupported"]))

    def test_mixed_model_counts_do_not_identify_the_failed_pair(self):
        rows = [{"model_id": "0", "outcome": "refused", "refusal_code": "self_intersecting",
                 "face_kinds": {"plane": 10, "cylinder": 2}}]
        errors = {"0": {"code": "self_intersecting"}}
        unknown = summarize(rows, errors, {})["self_intersections"][0]
        self.assertEqual(unknown["face_kinds"], ["unknown", "unknown"])
        detail = {"0": {"source_reference_audit": True,
                        "import_stage": "build_solid.check_solid_before_topology",
                        "intersection_path": "planar_polygon_pass", "face_kinds": ["plane", "plane"]}}
        known = summarize(rows, errors, detail)["self_intersections"][0]
        self.assertEqual(known["face_kinds"], ["plane", "plane"])

    def test_one_recorded_error_does_not_prove_no_other_refusal(self):
        rows = [{"model_id": "0", "outcome": "refused", "refusal_code": "unsupported"}]
        result = summarize(rows, {"0": {"code": "unsupported", "type_name": "OFFSET_SURFACE"}}, {})
        value = result["unsupported"][0]["counterfactual"]
        self.assertIsNone(value["no_other_refusal_models"])
        self.assertEqual(value["confirmed_models"], 0)
        self.assertEqual(value["unknown_models"], 1)
        self.assertEqual(value["possible_range"], [0, 1])
        self.assertIn("unknown [0, 1]", render(result))


if __name__ == "__main__":
    unittest.main()
