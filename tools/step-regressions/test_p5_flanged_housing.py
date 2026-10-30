"""P5 through the ordinary driver: a certified mixed body and a raw-valid STEP file."""
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GeomAbs import (GeomAbs_BSplineSurface, GeomAbs_Cone, GeomAbs_Cylinder, GeomAbs_Plane,
                         GeomAbs_Torus)
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_FACE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS

from test_loopless_faces import raw_solid
from test_partial_periodic_trims import BINARY, DRIVER

ROOT = Path(__file__).resolve().parents[2]
APPLICATION = 20_000_000
REPLAY = 100_000_000
FACES = {GeomAbs_Plane: 3, GeomAbs_Cylinder: 7, GeomAbs_Cone: 4, GeomAbs_Torus: 4,
         GeomAbs_BSplineSurface: 4}
SUMMARY = {'plane': 3, 'cylinder': 7, 'cone': 4, 'torus': 4, 'bspline': 4}
FILE_TO_MODEL = Fraction(1, 10_000_000)


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


def rotated():
    """Quaternion (1,2,3,4) on every frame, point and axis direction."""
    matrix = ((-20, 4, 22), (20, -10, 20), (10, 28, 4))
    turn = lambda vector: [
        str(sum(Fraction(value) * Fraction(entry, 30) for value, entry in zip(vector, row)))
        for row in matrix]
    def compose(rotation):
        w, x, y, z = map(int, rotation)
        return [str(w - 2 * x - 3 * y - 4 * z), str(x + 2 * w + 3 * z - 4 * y),
                str(y + 3 * w + 4 * x - 2 * z), str(z + 4 * w + 2 * y - 3 * x)]
    recipe = placed(canonical(), turn)
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
        BRepGProp.VolumePropertiesGK_s(shape, properties, tolerance, True, False, False, False,
                                       False)
        values.append(properties.Mass())
    return values


def surface_area(shape):
    properties = GProp_GProps()
    BRepGProp.SurfaceProperties_s(shape, properties)
    return properties.Mass()


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

    def test_rows_build_export_and_replay(self):
        rows = (
            ('canonical', canonical(), 12558.681699123754, 0.02),
            ('edited', radius('5/2'), 12591.905604508610, 0.03),
            ('translated', translated((7, 11, -13)), 12558.681699123754, 0.02),
            ('half', scaled('1/2'), 1569.8352123904692, 0.005),
            ('double', scaled('2'), 100469.45359299003, 0.08),
        )
        for name, recipe, reference, margin in rows:
            with self.subTest(name):
                directory, manifest = self.run_case(recipe, APPLICATION)
                self.assert_body(directory, manifest, reference, margin)
                wider, other = self.run_case(recipe, REPLAY)
                self.assertEqual(other['status'], 'answered')
                self.assertEqual((wider / 'answer.step').read_bytes(),
                                 (directory / 'answer.step').read_bytes())
                self.assertEqual(other['certificate'], manifest['certificate'])
                replay, again = self.run_case(
                    json.loads((directory / 'question.json').read_text()), APPLICATION)
                self.assertEqual(again['certificate'], manifest['certificate'])
                self.assertEqual((replay / 'answer.step').read_bytes(),
                                 (directory / 'answer.step').read_bytes())

    def test_uncovered_strata_keep_their_refusals(self):
        rows = (
            ('rotated', rotated(), 'fillet.unsupported-overflow'),
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
