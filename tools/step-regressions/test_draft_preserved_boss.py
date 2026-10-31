"""A placed drafted part keeps its unselected raised boss on the public route."""
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepGProp import BRepGProp
from OCP.GeomAbs import GeomAbs_Circle, GeomAbs_Cylinder, GeomAbs_Plane
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_IN, TopAbs_OUT
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

from test_p5_flanged_housing import (APPLICATION, FILE_TO_MODEL, REPLAY, ROOT, RQ, T, comparable,
                                     faces, integrated, moved, ordinary_shape, solids,
                                     surface_area)
from test_partial_periodic_trims import BINARY, DRIVER
from raw_step import raw_shapes

FIXTURES = ROOT / 'crates/catallax-features/tests/fixtures/draft-preserved-boss'
DIGEST = '91a5192fe48d89b87e7f4782b126d41aa2396bfeedfc0fc9a44608d7e8277615'
# Removing only the draft, or drafting before the boss, leaves these files untouched.
CONTROLS = {'M-U': '49c5c3133675ea6f00793768e5383fbf0393dcc95a65cfd49150dee3bbf3909e',
            'M-A': '5f254fadca235d061b4c3fe3ebc44ad8c6e82ed048f18b399a928f1294f038df'}
TANGENT = math.tan(math.radians(3))
SECANT = 1 / math.cos(math.radians(3))
BOSS = 80 * math.pi
# The base's section is 30(40-2zt) over 0 <= z <= 20. The boss's end replaces the
# shoulder disk it covers, so only its side adds area.
VOLUME = 24000 - 12000 * TANGENT + BOSS
AREA = 4000 - 2000 * TANGENT + 1200 * SECANT + 40 * math.pi
# Unplaced strict points: base, boss, past the boss, the joint, the slope, past the wall.
MATERIAL = (((0, 0, 10), TopAbs_IN), ((6, 2, 22), TopAbs_IN), ((11, 2, 22), TopAbs_OUT),
            ((6, 2, 20), TopAbs_IN), ((19, 0, 19), TopAbs_IN), ((21, 0, 10), TopAbs_OUT))
# Plane area across the axis by local height: the bottom, the shoulder less its joint
# and the boss's end.
CAPS = {0: 1200, 20: 1200 - 1200 * TANGENT - 16 * math.pi, 25: 16 * math.pi}
HEIGHT_AXIS = tuple(float(row[2]) for row in RQ)
TRANSLATION = tuple(float(value) for value in T)


def height(point):
    return sum(axis * (value - offset) for axis, value, offset
               in zip(HEIGHT_AXIS, (point.X(), point.Y(), point.Z()), TRANSLATION))


def each(shape, kind):
    explorer = TopExp_Explorer(shape, kind)
    while explorer.More():
        yield explorer.Current()
        explorer.Next()


def area_of(face):
    properties = GProp_GProps()
    BRepGProp.SurfaceProperties_s(face, properties)
    return properties.Mass()


def wall_area(shape):
    """Wall coverage. A periodic seam may split the one support into patches."""
    return sum(area_of(face) for face in each(shape, TopAbs_FACE)
               if BRepAdaptor_Surface(TopoDS.Face(face)).GetType() == GeomAbs_Cylinder)


def cap_areas(shape):
    totals = {}
    for face in each(shape, TopAbs_FACE):
        surface = BRepAdaptor_Surface(TopoDS.Face(face))
        if surface.GetType() != GeomAbs_Plane:
            continue
        plane = surface.Plane()
        normal = plane.Axis().Direction()
        along = normal.X() * HEIGHT_AXIS[0] + normal.Y() * HEIGHT_AXIS[1] \
            + normal.Z() * HEIGHT_AXIS[2]
        if abs(abs(along) - 1) < 1e-9:
            level = round(height(plane.Location()))
            totals[level] = totals.get(level, 0) + area_of(face)
    return totals


def ring_lengths(shape):
    """Circular edge length by local height. Each edge joins exactly two faces."""
    rings = []
    for face in each(shape, TopAbs_FACE):
        for found in each(face, TopAbs_EDGE):
            if BRepAdaptor_Curve(TopoDS.Edge(found)).GetType() != GeomAbs_Circle:
                continue
            for ring in rings:
                if ring[0].IsSame(found):
                    ring[1] += 1
                    break
            else:
                rings.append([found, 1])
    totals = {}
    for edge, users in rings:
        if users != 2:
            raise AssertionError(f'ring edge has {users} faces')
        properties = GProp_GProps()
        BRepGProp.LinearProperties_s(edge, properties)
        level = round(height(BRepAdaptor_Curve(TopoDS.Edge(edge)).Circle().Location()))
        totals[level] = totals.get(level, 0) + properties.Mass()
    return totals


class PreservedBoss(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='catallax-preserved-boss-')
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
        for face in each(shape, TopAbs_FACE):
            self.assertTrue(BRepCheck_Analyzer(face).IsValid())
        kinds = faces(shape)
        self.assertEqual(kinds.pop(GeomAbs_Plane), 7)
        self.assertGreaterEqual(kinds.pop(GeomAbs_Cylinder), 1)
        self.assertEqual(kinds, {})
        area = surface_area(shape)
        tolerance = area * 1e-7
        self.assertAlmostEqual(area, AREA, delta=tolerance)
        self.assertAlmostEqual(wall_area(shape), 40 * math.pi, delta=tolerance)
        covered = cap_areas(shape)
        self.assertEqual(sorted(covered), sorted(CAPS))
        for level, expected in CAPS.items():
            self.assertAlmostEqual(covered[level], expected, delta=tolerance)
        rings = ring_lengths(shape)
        self.assertEqual(sorted(rings), [20, 25], 'both rings stay shared')
        for length in rings.values():
            self.assertAlmostEqual(length, 8 * math.pi, delta=1e-6)
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
        source = FIXTURES / 'M-R.json'
        self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), DIGEST)
        recipe = json.loads(source.read_text())
        directory, response = self.run_case(recipe, APPLICATION)
        self.assertEqual(response['work']['limit'], APPLICATION)
        self.assertFalse(response['work']['exhausted'])
        manifest = response['manifest']
        self.assertEqual(manifest['format'], 'catallax-manifest/3')
        self.assertEqual(manifest['status'], 'answered')
        self.assertTrue(manifest['exact'])
        self.assertEqual(manifest['marks'], [])
        summary = manifest['result']
        self.assertEqual(summary['faces_by_kind'], {'plane': 7, 'cylinder': 1})
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
        fresh = self.root / 'fresh'
        replay = subprocess.run(
            [str(BINARY), 'replay', str(directory / 'question.json'), '--out', str(fresh)],
            capture_output=True, text=True, check=True)
        again = json.loads(replay.stdout)
        self.assertEqual(again['work']['limit'], APPLICATION)
        self.assertEqual(comparable(again['manifest']), comparable(manifest))
        self.assertEqual((fresh / 'answer.step').read_bytes(), step.read_bytes())
        wider, other = self.run_case(recipe, REPLAY)
        self.assertEqual(comparable(other['manifest']), comparable(manifest))
        self.assertEqual((wider / 'answer.step').read_bytes(), step.read_bytes())

    def test_undrafted_and_boss_after_draft_controls_keep_their_bytes(self):
        for name, digest in CONTROLS.items():
            with self.subTest(name=name):
                recipe = json.loads((FIXTURES / (name + '.json')).read_text())
                directory, _ = self.run_case(recipe, APPLICATION)
                written = (directory / 'answer.step').read_bytes()
                self.assertEqual(hashlib.sha256(written).hexdigest(), digest)


if __name__ == '__main__':
    unittest.main()
