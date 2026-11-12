import unittest

from planar_incidence_repair import distribution, public_row


class RepairEvidenceTests(unittest.TestCase):
    def test_distribution_keeps_tiny_decimal_bounds_and_nearest_ranks(self):
        result = distribution(["0.01", "0.00000000000000000001", "0", "0.002"])
        self.assertEqual(result["count"], 4)
        self.assertEqual(result["min_mm"], "0")
        self.assertEqual(result["p50_mm"], "1E-20")
        self.assertEqual(result["p95_mm"], "0.01")
        self.assertEqual(distribution([]), {"count": 0})

    def test_refusal_measurements_keep_their_scope(self):
        row = public_row({"outcome": "refused", "refusal": {
            "code": "planar_repair_not_verified", "stage": "curved_dependency",
            "measurement": "vertex_plane_residual", "deviation_upper_mm": "0.000001"}})
        self.assertEqual(row["refusal_measurement"], "vertex_plane_residual")
        self.assertEqual(row["refusal_upper_mm"], "0.000001")
        self.assertIsNone(row["max_source_to_result_upper_mm"])


if __name__ == "__main__":
    unittest.main()
