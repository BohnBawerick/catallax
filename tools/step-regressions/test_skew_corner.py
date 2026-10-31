"""A placed all-edges triangular prism: exact skew sphere sectors through the public driver."""
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane, GeomAbs_Sphere
from OCP.TopAbs import TopAbs_FACE, TopAbs_IN, TopAbs_OUT
from OCP.TopExp import TopExp_Explorer
from OCP.gp import gp_Pnt

from test_p5_flanged_housing import (APPLICATION, FILE_TO_MODEL, REPLAY, ROOT, RQ, T, comparable,
                                     faces, integrated, moved, ordinary_shape, solids,
                                     surface_area)
from test_partial_periodic_trims import BINARY, DRIVER
from raw_step import raw_shapes

SOURCE = ROOT / 'crates/catallax-cli/tests/fixtures/skew-corner/H-R.json'
DIGEST = 'a1b3cb43adf7b14f7f7e4a98036fad157636359162335193d7f27e8452636853'
# Parallel body of the inward-offset prism at r1: 156 + 46 pi / 3.
VOLUME = 156 + 46 * math.pi / 3
FACES = {GeomAbs_Plane: 5, GeomAbs_Cylinder: 9, GeomAbs_Sphere: 6}
# Unplaced points. The fourth is outside only by the sphere about (5,1,1); the
# fifth is on that sphere's support, outside its sector, and stays material.
MATERIAL = (((2, 2, 5), TopAbs_IN), ((8, 0, 0), TopAbs_OUT),
            ((Fraction(11, 2), Fraction(1, 2), Fraction(1, 2)), TopAbs_IN),
            ((Fraction(29, 5), Fraction(2, 5), Fraction(1, 5)), TopAbs_OUT),
            ((Fraction(13, 3), Fraction(1, 3), Fraction(2, 3)), TopAbs_IN))


class SkewCorners(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='catallax-skew-corner-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def run_case(self, problem, max_work):
        root = Path(tempfile.mkdtemp(dir=self.root))
        source = root / 'problem.json'
        source.write_text(json.dumps(problem))
        out = root / 'run'
        subprocess.run(
            [str(DRIVER), 'run', str(source), '--source', 'recipe', '--binary', str(BINARY),
             '--max-work', str(max_work), '--out', str(out), '--fresh', '--workers', '1',
             '--json'],
            capture_output=True, text=True, check=True)
        case = json.loads((out / 'run.json').read_text())['cases'][0]
        return Path(case['path']), case['response']['manifest']

    def assert_solid(self, shape):
        found = solids(shape)
        self.assertEqual(len(found), 1)
        self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
        explorer = TopExp_Explorer(shape, TopAbs_FACE)
        while explorer.More():
            self.assertTrue(BRepCheck_Analyzer(explorer.Current()).IsValid())
            explorer.Next()
        self.assertEqual(faces(shape), FACES)
        values = integrated(shape)
        self.assertTrue(all(math.isfinite(value) for value in values), values)
        spread = abs(values[2] - values[1])
        self.assertLessEqual(spread, 1e-10 * max(1, abs(values[2])))
        self.assertLessEqual(abs(values[2] - VOLUME),
                             max(1e-6, surface_area(shape) * 1e-7 + spread))
        for point, state in MATERIAL:
            probe = gp_Pnt(*(float(value) for value in moved(RQ, point, T)))
            self.assertEqual(BRepClass3d_SolidClassifier(found[0], probe, 1e-7).State(), state,
                             point)

    def test_placed_prism_exports_valid_sectors_and_replays(self):
        self.assertEqual(hashlib.sha256(SOURCE.read_bytes()).hexdigest(), DIGEST)
        recipe = json.loads(SOURCE.read_text())
        directory, manifest = self.run_case(recipe, APPLICATION)
        self.assertEqual(manifest['format'], 'catallax-manifest/3')
        self.assertEqual(manifest['status'], 'answered')
        self.assertTrue(manifest['exact'])
        self.assertEqual(manifest['marks'], [])
        summary = manifest['result']
        self.assertEqual(summary['faces_by_kind'], {'plane': 5, 'cylinder': 9, 'sphere': 6})
        self.assertEqual(summary['edges_by_kind'], {'line': 18, 'circle': 18})
        self.assertEqual((summary['lumps'], summary['shells'], summary['vertices']), (1, 1, 18))
        written = manifest['step']
        self.assertLessEqual(Fraction(written['max_deviation']), FILE_TO_MODEL)
        self.assertEqual(written['approximated_curves'], 0)
        step = directory / 'answer.step'
        shapes = raw_shapes(step)
        self.assertEqual(len(shapes), 1)
        self.assert_solid(shapes[0])
        self.assert_solid(ordinary_shape(step))
        replay, again = self.run_case(json.loads((directory / 'question.json').read_text()),
                                      APPLICATION)
        self.assertEqual(comparable(again), comparable(manifest))
        self.assertEqual((replay / 'answer.step').read_bytes(), step.read_bytes())
        wider, other = self.run_case(recipe, REPLAY)
        self.assertEqual(comparable(other), comparable(manifest))
        self.assertEqual((wider / 'answer.step').read_bytes(), step.read_bytes())


if __name__ == '__main__':
    unittest.main()
