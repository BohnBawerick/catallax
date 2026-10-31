"""Independent raw and ordinary imports of a planar draft, at its 3- and 5-degree angles."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.GeomAbs import GeomAbs_Plane
from OCP.TopAbs import TopAbs_FACE, TopAbs_IN, TopAbs_OUT
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

from test_p5_flanged_housing import (APPLICATION, ROOT, integrated, ordinary_shape, solids,
                                     surface_area)
from test_partial_periodic_trims import BINARY, DRIVER
from raw_step import raw_shapes

FIXTURES = ROOT / 'crates/catallax-cli/tests/fixtures/planar-draft'
# Integrated by height from 0 <= z <= 20, |y| <= 15, |x| <= 20 - z tan(a).
CASES = {
    3: ('fd3a4dbafc1cb6d24e5df5e68895ef5d17f34c00cadd71812badd82d8a5fa50d',
        Fraction('23371.106648603505551534330106')),
    5: ('fd4acf64213aded64a199e3122ef6293fa6db8af847ae776e5287d0524203586',
        Fraction('22950.136037688911937335775967')),
}
# (19, 0, 19) separates the angles: 19 tan(3 degrees) < 1 < 19 tan(5 degrees).
MATERIAL = {(0, 0, 10): {3: TopAbs_IN, 5: TopAbs_IN},
            (19, 0, 19): {3: TopAbs_IN, 5: TopAbs_OUT},
            (21, 0, 10): {3: TopAbs_OUT, 5: TopAbs_OUT},
            (0, 0, 21): {3: TopAbs_OUT, 5: TopAbs_OUT}}


class PlanarDraftImports(unittest.TestCase):
    def assert_solid(self, shape, volume, angle):
        found = solids(shape)
        self.assertEqual(len(found), 1)
        self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
        planes = 0
        explorer = TopExp_Explorer(shape, TopAbs_FACE)
        while explorer.More():
            self.assertTrue(BRepCheck_Analyzer(explorer.Current()).IsValid())
            self.assertEqual(BRepAdaptor_Surface(TopoDS.Face(explorer.Current())).GetType(),
                             GeomAbs_Plane)
            planes += 1
            explorer.Next()
        self.assertEqual(planes, 6)
        values = integrated(shape)
        spread = abs(values[2] - values[1])
        self.assertLessEqual(spread, 1e-10 * max(1, abs(values[2])))
        self.assertLessEqual(abs(values[2] - float(volume)),
                             max(1e-6, surface_area(shape) * 1e-7 + spread))
        for point, states in MATERIAL.items():
            state = BRepClass3d_SolidClassifier(found[0], gp_Pnt(*point), 1e-7).State()
            self.assertEqual(state, states[angle], point)

    def test_raw_and_ordinary_imports_match_the_independent_solid(self):
        with tempfile.TemporaryDirectory(prefix='catallax-planar-draft-') as temporary:
            for angle, (digest, volume) in CASES.items():
                with self.subTest(angle=angle):
                    source = FIXTURES / f'draft-{angle}.json'
                    self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), digest)
                    out = Path(temporary) / str(angle)
                    completed = subprocess.run(
                        [str(DRIVER), 'run', str(source), '--source', 'recipe', '--binary',
                         str(BINARY), '--max-work', str(APPLICATION), '--out', str(out),
                         '--fresh', '--workers', '1', '--json'],
                        check=True, capture_output=True, text=True)
                    self.assertEqual(json.loads(completed.stdout)['outcomes'], {'answered': 1})
                    case = json.loads((out / 'run.json').read_text())['cases'][0]
                    step = Path(case['path']) / 'answer.step'
                    self.assertGreater(step.stat().st_size, 0)
                    shapes = raw_shapes(step)
                    self.assertEqual(len(shapes), 1)
                    self.assert_solid(shapes[0], volume, angle)
                    self.assert_solid(ordinary_shape(step), volume, angle)


if __name__ == '__main__':
    unittest.main()
