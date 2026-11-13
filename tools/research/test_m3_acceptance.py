"""The M3 check keeps the authored recipes, edits and request limits fixed."""

import json
from pathlib import Path
import struct
import tempfile
import unittest

from m3_acceptance import (APP_WORK, PARTS, ROOT, SCOREBOARD_WORK, audit_near_gap,
                           audit_stl, digest, prepare)


class M3Inputs(unittest.TestCase):
    def test_stl_audit_rejects_empty_file_even_with_matching_hash(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "print.stl"
            path.write_bytes(bytes(80) + struct.pack("<I", 0))
            printing = {"triangles": 0, "vertices": 0, "sha256": digest(path),
                        "watertight": True, "route": "brep", "deviation_bound_mm": "0",
                        "completeness": "proven-two-sided"}
            with self.assertRaises(AssertionError):
                audit_stl(path, printing)

    def test_canonical_states_and_limits(self):
        with tempfile.TemporaryDirectory() as temporary:
            inputs = Path(temporary) / "inputs"
            prepare(inputs)
            requests = {p.stem: json.loads(p.read_text()) for p in inputs.glob("*.json")}
        self.assertEqual(len(requests), 14)
        for name, request in requests.items():
            self.assertEqual(request["max_work"], SCOREBOARD_WORK if name.startswith("P0") else APP_WORK)
        for part, (file, step, parameter, value, _) in PARTS.items():
            original = json.loads((ROOT / "parts" / file).read_text())
            self.assertEqual(requests[f"{part}-original-step"]["problem"], original)
            edited = requests[f"{part}-edited-step"]["problem"]
            self.assertNotEqual(original, edited)
            self.assertEqual(next(s for s in edited["steps"] if s["id"] == step)[parameter], value)
            next(s for s in edited["steps"] if s["id"] == step)[parameter] = next(
                s for s in original["steps"] if s["id"] == step)[parameter]
            self.assertEqual(edited, original)

    def test_stl_audit_rejects_a_reversed_triangle_even_with_matching_hash(self):
        points = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
        triangles = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "print.stl"
            for reversed_face in (False, True):
                rows = triangles.copy()
                if reversed_face:
                    rows[0] = tuple(reversed(rows[0]))
                data = bytes(80) + struct.pack("<I", len(rows))
                for row in rows:
                    data += struct.pack("<12fH", 0, 0, 0,
                                        *(x for index in row for x in points[index]), 0)
                path.write_bytes(data)
                printing = {"triangles": 4, "vertices": 4, "sha256": digest(path),
                            "watertight": True, "route": "brep", "deviation_bound_mm": "0.001",
                            "completeness": "proven-two-sided"}
                if reversed_face:
                    with self.assertRaises(AssertionError):
                        audit_stl(path, printing)
                else:
                    self.assertEqual(audit_stl(path, printing)["status"], "passed")

    def test_gap_audit_uses_file_decimals_and_rejects_contact_or_overlap(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "supports.step"
            for center in ("10.0000001", "10.", "9.9999999"):
                path.write_text(f"""#1=CYLINDRICAL_SURFACE('',#2,5.);
#2=AXIS2_PLACEMENT_3D('',#3,#4,#5);
#3=CARTESIAN_POINT('',(0.,0.,0.));
#4=DIRECTION('',(0.,0.,1.));
#5=DIRECTION('',(1.,0.,0.));
#6=CYLINDRICAL_SURFACE('',#7,5.);
#7=AXIS2_PLACEMENT_3D('',#8,#4,#5);
#8=CARTESIAN_POINT('',({center},0.,0.));
#9=PLANE('',#2);
#10=PLANE('',#2);
#11=PLANE('',#7);
#12=PLANE('',#7);
""")
                if center == "10.0000001":
                    self.assertEqual(audit_near_gap(path)["gap_mm"], "1/10000000")
                else:
                    with self.assertRaises(AssertionError):
                        audit_near_gap(path)


if __name__ == "__main__":
    unittest.main()
