"""P5 through the ordinary driver: a certified mixed body and a raw-valid STEP file."""
from fractions import Fraction
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepGProp import BRepGProp
from OCP.GeomAbs import (GeomAbs_BSplineSurface, GeomAbs_Cone, GeomAbs_Cylinder, GeomAbs_Plane,
                         GeomAbs_Torus)
from OCP.GProp import GProp_GProps
from OCP.gp import gp_Pnt
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_FACE, TopAbs_IN, TopAbs_OUT, TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS

from test_loopless_faces import raw_solid
from test_partial_periodic_trims import BINARY, DRIVER

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools/checks'))

from raw_step import raw_shapes  # noqa: E402
APPLICATION = 20_000_000
REPLAY = 100_000_000
FACES = {GeomAbs_Plane: 3, GeomAbs_Cylinder: 7, GeomAbs_Cone: 4, GeomAbs_Torus: 4,
         GeomAbs_BSplineSurface: 4}
SUMMARY = {'plane': 3, 'cylinder': 7, 'cone': 4, 'torus': 4, 'bspline': 4}
FILE_TO_MODEL = Fraction(1, 10_000_000)
RZ = ((Fraction(3, 5), Fraction(-4, 5), 0), (Fraction(4, 5), Fraction(3, 5), 0), (0, 0, 1))
RQ = ((Fraction(-2, 3), Fraction(2, 15), Fraction(11, 15)),
      (Fraction(2, 3), Fraction(-1, 3), Fraction(2, 3)),
      (Fraction(1, 3), Fraction(14, 15), Fraction(2, 15)))
T = (7, 11, -13)
# Definition points strictly away from the boundary: inside at r2, inside at r5/2.
MATERIAL = (((10, 15, 0), True, True), ((0, 15, 0), False, False), ((17, 2, 0), False, False),
            ((0, 2, 20), True, True), ((Fraction(15, 2), Fraction(11, 2), 10), True, True),
            ((Fraction(15, 2), Fraction(29, 5), 10), False, True))


def canonical():
    return json.loads((ROOT / 'parts/P5-flanged-housing.json').read_text())


def step(recipe, name):
    return next(item for item in recipe['steps'] if item['id'] == name)


def radius(value):
    recipe = canonical()
    step(recipe, 'filleted')['radius'] = value
    return recipe


def placed(recipe, edit):
    """Every definition point of the sketch plane, both tool poses and both axes."""
    for point in (step(recipe, 'section')['plane']['translation'],
                  step(recipe, 'housing')['axis']['point'],
                  step(recipe, 'drill')['pose']['translation'],
                  step(recipe, 'sink')['pose']['translation'],
                  step(recipe, 'hole-pattern')['layout']['circular']['axis']['point']):
        point[:] = edit(point)
    return recipe


def translated(by):
    return placed(canonical(), lambda point: [
        str(Fraction(value) + offset) for value, offset in zip(point, by)])


def scaled(factor):
    k = Fraction(factor)
    times = lambda value: str(Fraction(value) * k)
    recipe = placed(canonical(), lambda point: [times(value) for value in point])
    for segment in step(recipe, 'section')['loops'][0]:
        segment['to'] = [times(value) for value in segment['to']]
    for name, lengths in (('drill', ('radius', 'height')),
                          ('sink', ('base_radius', 'top_radius', 'height')),
                          ('filleted', ('radius',))):
        for length in lengths:
            step(recipe, name)[length] = times(step(recipe, name)[length])
    return recipe


def moved(matrix, point, by=(0, 0, 0)):
    return [sum(Fraction(value) * entry for value, entry in zip(point, row)) + offset
            for row, offset in zip(matrix, by)]


def rigid(matrix, quaternion, by=(0, 0, 0), value='2'):
    """The rotation on every frame, point and axis direction, before evaluation."""
    turn = lambda point, by=(0, 0, 0): [str(part) for part in moved(matrix, point, by)]
    def compose(rotation):
        a, b, c, d = quaternion
        w, x, y, z = map(int, rotation)
        return [str(a * w - b * x - c * y - d * z), str(a * x + b * w + c * z - d * y),
                str(a * y - b * z + c * w + d * x), str(a * z + b * y - c * x + d * w)]
    recipe = placed(radius(value), lambda point: turn(point, by))
    for name, key in (('section', 'plane'), ('drill', 'pose'), ('sink', 'pose')):
        step(recipe, name)[key]['rotation'] = compose(step(recipe, name)[key]['rotation'])
    axis = step(recipe, 'housing')['axis']
    axis['direction'] = turn(axis['direction'])
    pattern = step(recipe, 'hole-pattern')['layout']['circular']['axis']
    pattern['direction'] = turn(pattern['direction'])
    return recipe


def faces(shape):
    found = {}
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = TopoDS.Face(explorer.Current())
        kind = BRepAdaptor_Surface(face).GetType()
        found[kind] = found.get(kind, 0) + 1
        explorer.Next()
    return found


def integrated(shape):
    values = []
    for tolerance in (1e-6, 1e-7, 1e-8):
        properties = GProp_GProps()
        error = BRepGProp.VolumePropertiesGK_s(shape, properties, tolerance, True, False, False,
                                               False, False)
        if error < 0:
            raise AssertionError(f'volume integration failed: {error}')
        values.append(properties.Mass())
    return values


def surface_area(shape):
    properties = GProp_GProps()
    BRepGProp.SurfaceProperties_s(shape, properties)
    return properties.Mass()


def solids(shape):
    found = []
    explorer = TopExp_Explorer(shape, TopAbs_SOLID)
    while explorer.More():
        found.append(TopoDS.Solid(explorer.Current()))
        explorer.Next()
    return found


def ordinary_shape(path):
    reader = STEPControl_Reader()
    if reader.ReadFile(str(path)) != IFSelect_RetDone or reader.TransferRoots() < 1:
        raise AssertionError(f'cannot read {path}')
    return reader.OneShape()


def comparable(manifest):
    """A manifest without the build identity, the one part a rerun may change."""
    return {key: value for key, value in manifest.items() if key != 'kernel'}


class FlangedHousing(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='catallax-p5-')
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

    def assert_body(self, directory, manifest, reference, margin):
        self.assertEqual(manifest['format'], 'catallax-manifest/4')
        self.assertEqual(manifest['status'], 'answered')
        self.assertFalse(manifest['exact'])
        self.assertEqual(manifest['result']['faces_by_kind'], SUMMARY)
        regions = manifest['certificate']['regions']
        self.assertEqual(len(regions), 22)
        fitted = [region for region in regions if region['exception'] is not None]
        self.assertEqual({region['role'] for region in fitted}, {'rolling-envelope'})
        self.assertEqual(len(fitted), 4)
        for region in fitted:
            self.assertLessEqual(Fraction(region['exception']['enclosure']['hi']),
                                 Fraction(1, 10_000))
        written = manifest['step']
        self.assertEqual(written['status'], 'answered')
        self.assertLessEqual(Fraction(written['file_to_model_bound_mm']), FILE_TO_MODEL)
        self.assertEqual(len(written['regions']), 22)
        volume = manifest['properties']['values']['volume']
        self.assertEqual(manifest['properties']['scope'], 'original-definition')
        self.assertLess(abs((Fraction(volume['lo']) + Fraction(volume['hi'])) / 2
                            - Fraction(reference)), Fraction(1, 100_000))
        shape = raw_solid(directory / 'answer.step')
        self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
        self.assertEqual(faces(shape), FACES)
        values = integrated(shape)
        self.assertLessEqual(abs(values[2] - values[1]), 1e-10 * max(1, abs(values[2])))
        spread = abs(values[2] - values[1])
        self.assertLessEqual(abs(values[2] - reference),
                             margin + surface_area(shape) * 1e-7 + spread)

    def assert_placed_file(self, path, manifest, reference, margin, matrix, by, edited):
        """The file itself: raw and ordinary solids, volume and material at the placed pose."""
        written = manifest['step']['regions']
        self.assertEqual({region['model_region'] for region in written},
                         {region['id'] for region in manifest['certificate']['regions']})
        self.assertEqual(len({region['file_region'] for region in written}), 22)
        carried = {}
        for region in manifest['certificate']['regions']:
            bound = (region['exception'] or {}).get('enclosure', {}).get('hi', 0)
            carried[region['id']] = max(
                Fraction(bound), Fraction(region['domain'].get('fitted_trim_bound_mm', 0)))
        for region in written:
            model = carried[region['model_region']]
            self.assertLessEqual(Fraction(region['file_to_model_bound_mm']), FILE_TO_MODEL)
            self.assertLessEqual(Fraction(region['file_to_definition_bound_mm']),
                                 FILE_TO_MODEL + model)
        shapes = raw_shapes(path)
        self.assertEqual(len(shapes), 1)
        for shape in (shapes[0], ordinary_shape(path)):
            found = solids(shape)
            self.assertEqual(len(found), 1)
            self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
            explorer = TopExp_Explorer(shape, TopAbs_FACE)
            while explorer.More():
                self.assertTrue(BRepCheck_Analyzer(explorer.Current()).IsValid())
                explorer.Next()
            values = integrated(shape)
            self.assertTrue(all(math.isfinite(value) for value in values), values)
            spread = abs(values[2] - values[1])
            self.assertLessEqual(spread, 1e-10 * max(1, abs(values[2])))
            self.assertLessEqual(abs(values[2] - reference),
                                 margin + surface_area(shape) * 1e-7 + spread)
            for point, *inside in MATERIAL:
                probe = gp_Pnt(*(float(value) for value in moved(matrix, point, by)))
                state = BRepClass3d_SolidClassifier(found[0], probe, 1e-7).State()
                self.assertEqual(state, TopAbs_IN if inside[edited] else TopAbs_OUT, point)

    def test_rows_build_export_and_replay(self):
        rows = (
            ('canonical', canonical(), 12558.681699123754, 0.02, None),
            ('edited', radius('5/2'), 12591.905604508610, 0.03, None),
            ('translated', translated(T), 12558.681699123754, 0.02, None),
            ('half', scaled('1/2'), 1569.8352123904692, 0.005, None),
            ('double', scaled('2'), 100469.45359299003, 0.08, None),
            ('rotated', rigid(RQ, (1, 2, 3, 4)), 12558.681699123754, 0.02,
             (RQ, (0, 0, 0), False)),
            ('moved edit', rigid(RZ, (2, 0, 0, 1), T, '5/2'), 12591.905604508610, 0.03,
             (RZ, T, True)),
        )
        for name, recipe, reference, margin, pose in rows:
            with self.subTest(name):
                directory, manifest = self.run_case(recipe, APPLICATION)
                self.assert_body(directory, manifest, reference, margin)
                if pose:
                    self.assert_placed_file(directory / 'answer.step', manifest, reference,
                                            margin, *pose)
                wider, other = self.run_case(recipe, REPLAY)
                self.assertEqual(other['status'], 'answered')
                self.assertEqual((wider / 'answer.step').read_bytes(),
                                 (directory / 'answer.step').read_bytes())
                self.assertEqual(comparable(other), comparable(manifest))
                replay, again = self.run_case(
                    json.loads((directory / 'question.json').read_text()), APPLICATION)
                self.assertEqual(comparable(again), comparable(manifest))
                self.assertEqual((replay / 'answer.step').read_bytes(),
                                 (directory / 'answer.step').read_bytes())

    def test_uncovered_strata_keep_their_refusals(self):
        rows = (
            ('near-miss', radius('999/1000'), 'fillet.unsupported-overflow'),
            ('tangent', radius('1'), 'fillet.unsupported-overflow'),
            ('too-large', radius('7/2'), 'shapes.rolling-envelope'),
        )
        for name, recipe, code in rows:
            with self.subTest(name):
                directory, manifest = self.run_case(recipe, APPLICATION)
                self.assertEqual(manifest['status'], 'refused')
                self.assertEqual(manifest['refusal']['code'], code)
                self.assertFalse((directory / 'answer.step').exists())

    def test_a_small_allowance_claims_nothing(self):
        directory, manifest = self.run_case(canonical(), 1_000_000)
        self.assertEqual(manifest['status'], 'refused')
        self.assertEqual(manifest['refusal']['code'], 'numbers.operation-budget')
        self.assertIsNone(manifest.get('result'))
        self.assertFalse((directory / 'answer.step').exists())

    def test_an_impossible_export_writes_no_file(self):
        problem = {'format': 'catallax-certified-recipe/1', 'recipe': canonical(),
                   'outputs': {'step': {'file_tolerance': '0', 'decimal_places': 15}}}
        directory, manifest = self.run_case(problem, APPLICATION)
        self.assertEqual(manifest['status'], 'answered')
        self.assertEqual(manifest['step']['status'], 'refused')
        self.assertEqual(manifest['step']['refusal']['code'], 'step.tolerance-not-reached')
        self.assertFalse((directory / 'answer.step').exists())


if __name__ == '__main__':
    unittest.main()
