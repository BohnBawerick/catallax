"""Independent raw and ordinary imports of a planar draft, at its 3- and 5-degree angles."""
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
from OCP.GeomAbs import GeomAbs_Plane
from OCP.TopAbs import (TopAbs_EDGE, TopAbs_FACE, TopAbs_FORWARD, TopAbs_IN, TopAbs_OUT,
                        TopAbs_VERTEX)
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

from test_draft_preserved_boss import cap_areas
from test_loopless_faces import unique
from test_p5_flanged_housing import (APPLICATION, FILE_TO_MODEL, REPLAY, ROOT, RQ, T, comparable,
                                     faces, integrated, moved, ordinary_shape, solids,
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

BRACKET = ROOT / 'crates/catallax-features/tests/fixtures/concave-draft'
BRACKET_INPUT = '06b794a651ad6bcf026219588baed13a8bef228f9ac34ddaf295ff06474ccec5'
BRACKET_CONTROLS = {'R-PU': 'f8b8b5eb6ee9301e16d92c1b9f19889b2f69cfb627d33c0afef87b060e02346e',
                    'R-P45': '95731a0c75674a8b5f49f106da30841162e054bb99a0d5a0e348ab49b7e20c90'}
# The L at 5 degrees has section 800 - 60 z t over 0 <= z <= 20. Half scale divides its
# volume by 8 and its area by 4, so the unmoved cap of 800 and the drafted end have these
# areas at local heights 0 and 10.
TANGENT = math.tan(math.radians(5))
BRACKET_VOLUME = (16000 - 12000 * TANGENT) / 8
BRACKET_AREA = (3200 - 2000 * TANGENT + 1200 / math.cos(math.radians(5))) / 4
BRACKET_CAPS = {0: 800 / 4, 10: (800 - 1200 * TANGENT) / 4}
# Definition points. The third is on the open seam between the two arms. The last two
# are inside at 3 degrees only: 19 tan(3 degrees) < 1 < 19 tan(5 degrees).
BRACKET_MATERIAL = (((-10, -10, 10), TopAbs_IN), ((-10, 5, 10), TopAbs_IN),
                    ((-10, -5, 10), TopAbs_IN), ((5, 5, 10), TopAbs_OUT),
                    ((21, -10, 10), TopAbs_OUT), ((19, -10, 19), TopAbs_OUT),
                    ((-1, 0, 19), TopAbs_OUT))


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


class ConcaveBracket(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='catallax-concave-draft-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def run_case(self, source, max_work):
        out = Path(tempfile.mkdtemp(dir=self.root)) / 'run'
        subprocess.run(
            [str(DRIVER), 'run', str(source), '--source', 'recipe', '--binary', str(BINARY),
             '--max-work', str(max_work), '--out', str(out), '--fresh', '--workers', '1',
             '--json'],
            capture_output=True, text=True, check=True)
        case = json.loads((out / 'run.json').read_text())['cases'][0]
        return Path(case['path']), case['response']

    def assert_bracket(self, shape):
        found = solids(shape)
        self.assertEqual(len(found), 1)
        self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
        planes = unique(shape, TopAbs_FACE)
        edges = unique(shape, TopAbs_EDGE)
        for face in planes:
            self.assertTrue(BRepCheck_Analyzer(face).IsValid())
        self.assertEqual(faces(shape), {GeomAbs_Plane: 8})
        self.assertEqual((len(edges), len(unique(shape, TopAbs_VERTEX))), (18, 12))
        for edge in edges:
            senses = []
            for face in planes:
                explorer = TopExp_Explorer(face, TopAbs_EDGE)
                while explorer.More():
                    if explorer.Current().IsSame(edge):
                        senses.append(explorer.Current().Orientation() == TopAbs_FORWARD)
                    explorer.Next()
            self.assertEqual(sorted(senses), [False, True], 'each edge joins two opposed uses')
        area = surface_area(shape)
        tolerance = area * 1e-7
        self.assertAlmostEqual(area, BRACKET_AREA, delta=tolerance)
        caps = cap_areas(shape)
        self.assertEqual(sorted(caps), sorted(BRACKET_CAPS))
        for level, expected in BRACKET_CAPS.items():
            self.assertAlmostEqual(caps[level], expected, delta=tolerance)
        values = integrated(shape)
        self.assertTrue(all(math.isfinite(value) for value in values), values)
        spread = abs(values[2] - values[1])
        self.assertLessEqual(spread, 1e-10 * max(1, abs(values[2])))
        self.assertLessEqual(abs(values[2] - BRACKET_VOLUME), max(1e-6, area * 1e-7 + spread))
        for point, state in BRACKET_MATERIAL:
            half = [Fraction(value, 2) for value in point]
            probe = gp_Pnt(*(float(value) for value in moved(RQ, half, T)))
            self.assertEqual(BRepClass3d_SolidClassifier(found[0], probe, 1e-7).State(), state,
                             point)

    def test_placed_bracket_answers_raw_valid_and_replays_at_both_allowances(self):
        source = BRACKET / 'R-PRS.json'
        self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), BRACKET_INPUT)
        directory, response = self.run_case(source, APPLICATION)
        self.assertEqual(response['work']['limit'], APPLICATION)
        self.assertFalse(response['work']['exhausted'])
        manifest = response['manifest']
        self.assertEqual(manifest['format'], 'catallax-manifest/3')
        self.assertEqual(manifest['status'], 'answered')
        self.assertTrue(manifest['exact'])
        self.assertEqual(manifest['marks'], [])
        summary = manifest['result']
        self.assertEqual(summary['faces_by_kind'], {'plane': 8})
        self.assertEqual(summary['edges_by_kind'], {'line': 18})
        self.assertEqual((summary['lumps'], summary['shells'], summary['vertices']), (1, 1, 12))
        written = manifest['step']
        self.assertLessEqual(Fraction(written['max_deviation']), FILE_TO_MODEL)
        self.assertEqual(written['approximated_curves'], 0)
        step = directory / 'answer.step'
        self.assertGreater(step.stat().st_size, 0)
        shapes = raw_shapes(step)
        self.assertEqual(len(shapes), 1)
        self.assert_bracket(shapes[0])
        self.assert_bracket(ordinary_shape(step))
        question = json.loads((directory / 'question.json').read_text())
        self.assertEqual(question['max_work'], APPLICATION)
        fresh = self.root / 'fresh'
        replay = subprocess.run(
            [str(BINARY), 'replay', str(directory / 'question.json'), '--out', str(fresh)],
            capture_output=True, text=True, check=True)
        again = json.loads(replay.stdout)
        self.assertEqual(again['work']['limit'], APPLICATION)
        self.assertEqual(comparable(again['manifest']), comparable(manifest))
        self.assertEqual((fresh / 'answer.step').read_bytes(), step.read_bytes())
        wider, other = self.run_case(source, REPLAY)
        self.assertEqual(comparable(other['manifest']), comparable(manifest))
        self.assertEqual((wider / 'answer.step').read_bytes(), step.read_bytes())

    def test_undrafted_and_exact_angle_controls_keep_their_bytes(self):
        for name, digest in BRACKET_CONTROLS.items():
            with self.subTest(name=name):
                directory, _ = self.run_case(BRACKET / (name + '.json'), APPLICATION)
                written = (directory / 'answer.step').read_bytes()
                self.assertEqual(hashlib.sha256(written).hexdigest(), digest)


if __name__ == '__main__':
    unittest.main()
