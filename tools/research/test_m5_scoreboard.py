import copy
from fractions import Fraction
import hashlib
import json
import unittest

from m5_scoreboard import ROOT, import_wrapper, p0_copy, quaternion_product


class ScoreboardTests(unittest.TestCase):
    def test_copies_preserve_operations_and_scale_the_gap_exactly(self):
        source = json.loads((ROOT / "parts/P0/near-miss.json").read_text())
        saved = copy.deepcopy(source)
        scaled = p0_copy(source, "scaled")
        gap = Fraction(scaled["steps"][1]["pose"]["translation"][0]) - 2 * Fraction(scaled["steps"][0]["radius"])
        self.assertEqual(gap, Fraction("0.0000002"))
        translated = p0_copy(source, "translated")
        for a, b in zip(source["steps"], translated["steps"]):
            if "pose" in a:
                self.assertEqual([Fraction(y)-Fraction(x) for x, y in
                                  zip(a["pose"]["translation"], b["pose"]["translation"])],
                                 [101, -37, 19])
        self.assertEqual(source, saved)
        for variant in ("translated", "scaled", "seam-shifted"):
            result = p0_copy(source, variant)
            self.assertEqual(result["steps"][-1], source["steps"][-1])
            self.assertEqual(result["result"], source["result"])

    def test_seam_rotation_preserves_the_round_stock_axis(self):
        for rotation in (["1", "0", "0", "0"], ["2", "1", "1", "0"], ["3", "0", "0", "4"]):
            def axis(q):
                w, x, y, z = map(Fraction, q)
                n = w*w+x*x+y*y+z*z
                return [(2*x*z+2*w*y)/n, (2*y*z-2*w*x)/n, (w*w-x*x-y*y+z*z)/n]
            self.assertEqual(axis(rotation), axis(quaternion_product(rotation, [3, 0, 0, 4])))

    def test_wrapper_keeps_original_bytes_and_content_address(self):
        data = b"ISO-10303-21;\nEND-ISO-10303-21;\n"
        result = import_wrapper(data, "00000190")
        sha = hashlib.sha256(data).hexdigest()
        self.assertEqual(result["files"], {sha: data.decode()})
        self.assertEqual(result["recipe"]["steps"][0]["file"]["sha256"], sha)


if __name__ == "__main__":
    unittest.main()
