"""A placed drafted part keeps its unselected through bore on the public route."""
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepGProp import BRepGProp
from OCP.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_FACE, TopAbs_IN, TopAbs_OUT
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

from test_p5_flanged_housing import (APPLICATION, FILE_TO_MODEL, REPLAY, ROOT, RQ, T, comparable,
                                     faces, integrated, moved, ordinary_shape, solids,
                                     surface_area)
from test_partial_periodic_trims import BINARY, DRIVER
from raw_step import raw_shapes

FIXTURES = ROOT / 'crates/catallax-features/tests/fixtures/draft-preserved-bore'
SOURCE = FIXTURES / 'K-R.json'
DIGEST = 'bb1581e3fb960bed5896d5ac2fa93f702223ea332d07eb651c0a53a5d1838df6'
# Removing only the draft, or only the bore, leaves these files untouched.
CONTROLS = {'K-U': 'ab250fae028aa7cb67098276638eb502fec9239c52a0145508e57b0adbe25664',
            'K-P': 'c132c3ef6cc74bf72bdb80f9e29bc6e2c182783010eee6c7996c4741fc27a755'}
# The section at height z is 30(40-2zt)-25pi over 0 <= z <= 20.
TANGENT = math.tan(math.radians(3))
VOLUME = 24000 - 12000 * TANGENT - 500 * math.pi
BORE_AREA = 200 * math.pi
# Unplaced strict points; the last lies on the cap's plane, in its hole.
MATERIAL = (((10, 0, 10), TopAbs_IN), ((0, 0, 10), TopAbs_OUT), ((19, 0, 19), TopAbs_IN),
            ((21, 0, 10), TopAbs_OUT), ((0, 0, 20), TopAbs_OUT))


def bore_area(shape):
    """Bore coverage in the file. A periodic seam may split the support into patches."""
    total = 0
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = TopoDS.Face(explorer.Current())
        surface = BRepAdaptor_Surface(face)
        if surface.GetType() == GeomAbs_Cylinder:
            if abs(surface.Cylinder().Radius() - 5) > 1e-9:
                raise AssertionError('the bore radius changed')
            properties = GProp_GProps()
            BRepGProp.SurfaceProperties_s(face, properties)
            total += properties.Mass()
        explorer.Next()
    return total


class PreservedBore(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='catallax-preserved-bore-')
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
        return Path(case['path']), case['response']

    def assert_solid(self, shape):
        found = solids(shape)
        self.assertEqual(len(found), 1)
        self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
        explorer = TopExp_Explorer(shape, TopAbs_FACE)
        while explorer.More():
            self.assertTrue(BRepCheck_Analyzer(explorer.Current()).IsValid())
            explorer.Next()
        kinds = faces(shape)
        self.assertEqual(kinds.pop(GeomAbs_Plane), 6)
        self.assertGreaterEqual(kinds.pop(GeomAbs_Cylinder), 1)
        self.assertEqual(kinds, {})
        area = surface_area(shape)
        self.assertAlmostEqual(bore_area(shape), BORE_AREA, delta=area * 1e-7)
        values = integrated(shape)
        self.assertTrue(all(math.isfinite(value) for value in values), values)
        spread = abs(values[2] - values[1])
        self.assertLessEqual(spread, 1e-10 * max(1, abs(values[2])))
        self.assertLessEqual(abs(values[2] - VOLUME), max(1e-6, area * 1e-7 + spread))
        for point, state in MATERIAL:
            probe = gp_Pnt(*(float(value) for value in moved(RQ, point, T)))
            self.assertEqual(BRepClass3d_SolidClassifier(found[0], probe, 1e-7).State(), state,
                             point)

    def test_placed_part_answers_raw_valid_and_replays_at_both_allowances(self):
        self.assertEqual(hashlib.sha256(SOURCE.read_bytes()).hexdigest(), DIGEST)
        recipe = json.loads(SOURCE.read_text())
        directory, response = self.run_case(recipe, APPLICATION)
        self.assertEqual(response['work']['limit'], APPLICATION)
        self.assertFalse(response['work']['exhausted'])
        manifest = response['manifest']
        self.assertEqual(manifest['format'], 'catallax-manifest/3')
        self.assertEqual(manifest['status'], 'answered')
        self.assertTrue(manifest['exact'])
        self.assertEqual(manifest['marks'], [])
        summary = manifest['result']
        self.assertEqual(summary['faces_by_kind'], {'plane': 6, 'cylinder': 1})
        self.assertEqual(summary['edges_by_kind'], {'line': 12, 'circle': 2})
        self.assertEqual((summary['lumps'], summary['shells'], summary['vertices']), (1, 1, 8))
        written = manifest['step']
        self.assertLessEqual(Fraction(written['max_deviation']), FILE_TO_MODEL)
        self.assertEqual(written['approximated_curves'], 0)
        step = directory / 'answer.step'
        self.assertGreater(step.stat().st_size, 0)
        shapes = raw_shapes(step)
        self.assertEqual(len(shapes), 1)
        self.assert_solid(shapes[0])
        self.assert_solid(ordinary_shape(step))
        question = json.loads((directory / 'question.json').read_text())
        self.assertEqual(question['max_work'], APPLICATION)
        replay = subprocess.run(
            [str(BINARY), 'replay', str(directory / 'question.json'), '--out',
             str(self.root / 'fresh')],
            capture_output=True, text=True, check=True)
        again = json.loads(replay.stdout)
        self.assertEqual(again['work']['limit'], APPLICATION)
        self.assertEqual(comparable(again['manifest']), comparable(manifest))
        self.assertEqual((self.root / 'fresh/answer.step').read_bytes(), step.read_bytes())
        wider, other = self.run_case(recipe, REPLAY)
        self.assertEqual(comparable(other['manifest']), comparable(manifest))
        self.assertEqual((wider / 'answer.step').read_bytes(), step.read_bytes())

    def test_undrafted_and_unbored_controls_keep_their_bytes(self):
        for name, digest in CONTROLS.items():
            with self.subTest(name=name):
                recipe = json.loads((FIXTURES / (name + '.json')).read_text())
                directory, _ = self.run_case(recipe, APPLICATION)
                written = (directory / 'answer.step').read_bytes()
                self.assertEqual(hashlib.sha256(written).hexdigest(), digest)


if __name__ == '__main__':
    unittest.main()
