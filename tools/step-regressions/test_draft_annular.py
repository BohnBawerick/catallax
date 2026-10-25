"""Independent STEP imports, analytic volume and hole membership; no repair."""
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import runpy
import subprocess
import tempfile
import unittest

from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.GeomAbs import GeomAbs_Cone, GeomAbs_Cylinder, GeomAbs_Plane
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_FACE, TopAbs_IN, TopAbs_OUT
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
CHECK = runpy.run_path(str(DRIVER.with_name('check-step.py')))['check']
FIXTURES = ROOT / 'crates/catallax-features/tests/fixtures/draft-annular'
FACE = getattr(TopoDS, 'Face_s', None) or TopoDS.Face


class AnnularDraftImports(unittest.TestCase):
    def test_wall_choices_signed_angle_and_identity_keep_the_hole(self):
        with tempfile.TemporaryDirectory(prefix='catallax-draft-annular-') as temporary:
            for name, volume, cones in [('both', Fraction(80), 2), ('outer', Fraction(272, 3), 1),
                                       ('inner', Fraction(328, 3), 1), ('negative', Fraction(120), 2),
                                       ('zero', Fraction(120), 0)]:
                with self.subTest(name=name):
                    out = Path(temporary) / name
                    completed = subprocess.run([str(DRIVER), 'run', str(FIXTURES / (name + '.json')),
                                                '--source', 'recipe', '--binary', str(BINARY), '--out', str(out),
                                                '--fresh', '--workers', '1', '--json'], check=True,
                                               capture_output=True, text=True)
                    self.assertEqual(json.loads(completed.stdout)['outcomes'], {'answered': 1})
                    case = json.loads((out / 'run.json').read_text())['cases'][0]
                    manifest = case['response']['manifest']
                    self.assertTrue(manifest['exact'])
                    self.assertEqual(manifest['marks'], [])
                    step = Path(case['path']) / 'answer.step'
                    report = CHECK(str(step), 4)
                    self.assertEqual(report['status'], 'passed', report)
                    self.assertAlmostEqual(report['volume_mm3'], float(volume) * math.pi, delta=1e-7)
                    reader = STEPControl_Reader()
                    reader.ReadFile(str(step))
                    reader.TransferRoots()
                    shape = reader.OneShape()
                    for xyz, expected in [((0, 0, 1), TopAbs_OUT), ((5, 0, 1), TopAbs_IN),
                                          ((5, 0, -1), TopAbs_OUT), ((5, 0, 3), TopAbs_OUT)]:
                        state = BRepClass3d_SolidClassifier(shape, gp_Pnt(*xyz), 1e-7).State()
                        self.assertEqual(state, expected)
                    explorer = TopExp_Explorer(shape, TopAbs_FACE)
                    counts = {GeomAbs_Plane: 0, GeomAbs_Cylinder: 0, GeomAbs_Cone: 0}
                    while explorer.More():
                        kind = BRepAdaptor_Surface(FACE(explorer.Current())).GetType()
                        self.assertIn(kind, counts)
                        counts[kind] += 1
                        explorer.Next()
                    self.assertEqual(counts[GeomAbs_Plane], 2)
                    self.assertEqual(counts[GeomAbs_Cone], cones)
                    self.assertEqual(counts[GeomAbs_Cylinder], 2 - cones)


if __name__ == '__main__':
    unittest.main()
