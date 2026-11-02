"""A rounded ring between nonparallel cylinders, read without ShapeFix or sewing."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_VERTEX

from test_loopless_faces import BINARY, DRIVER, raw_solid, unique

LEVEL = ['1', '0', '0', '0']
# A quarter turn about x, which lays a cylinder along y.
FLAT = ['1', '1', '0', '0']
# Leans a cylinder's axis to (0, 3/5, 4/5).
LEAN = ['3', '-1', '0', '0']
FIT_MARK = 'features.rolling-envelope-fit'


def cylinder(name, radius, height, rotation, at):
    return {'id': name, 'op': 'cylinder', 'radius': radius, 'height': height,
            'pose': {'rotation': rotation, 'translation': at}}


def ring(kind, narrow, rounded):
    steps = [cylinder('wide', '10', '60', FLAT, ['0', '30', '0']), narrow,
             {'id': 'joined', 'op': 'boolean', 'kind': kind, 'target': 'wide', 'tool': 'wall'}]
    if rounded:
        side = {'side': {'segment': 'wall'}}
        steps.append({'id': 'round', 'op': 'fillet', 'body': 'joined', 'radius': '1',
                      'edges': [{'between': [{'step': 'wide', 'role': side},
                                             {'step': 'wall', 'role': side}]}]})
    return {'format': 'catallax-recipe/3', 'name': 'ring', 'steps': steps,
            'result': 'round' if rounded else 'joined'}


WALL = cylinder('wall', '4', '30', LEVEL, ['1', '0', '0'])
LEANING = cylinder('wall', '4', '30', LEAN, ['1', '0', '0'])
THROUGH = cylinder('wall', '4', '60', LEVEL, ['1', '0', '-30'])
# The change in volume is the kernel's, which matches an independent section
# integral. The counts are faces, edges and vertices. The band and each wall add
# one seam, and each closed edge gains one vertex.
FIXTURES = [
    ('boss', 'union', WALL, 3.53786, (6, 8, 5)),
    ('bore', 'subtract', WALL, -11.64886, (6, 8, 5)),
    ('leaning-boss', 'union', LEANING, 6.46506, (6, 8, 5)),
]


class RingFilletImports(unittest.TestCase):
    def read(self, root, name, recipe):
        source = root / f'{name}.json'
        source.write_text(json.dumps(recipe))
        out = root / name
        completed = subprocess.run(
            [str(DRIVER), 'run', str(source), '--source', 'recipe', '--binary',
             str(BINARY), '--out', str(out), '--fresh', '--workers', '1', '--json'],
            capture_output=True, text=True, check=True)
        self.assertEqual(json.loads(completed.stdout)['outcomes'], {'answered': 1})
        case = json.loads((out / 'run.json').read_text())['cases'][0]
        manifest = case['response']['manifest']
        self.assertEqual(manifest['status'], 'answered')
        return manifest, raw_solid(Path(case['path']) / 'answer.step')

    def rounded(self, root, name, recipe):
        manifest, shape = self.read(root, name, recipe)
        self.assertEqual(manifest['format'], 'catallax-manifest/4')
        self.assertFalse(manifest['exact'])
        self.assertEqual([mark['code'] for mark in manifest['marks']], [FIT_MARK])
        self.assertEqual(manifest['step']['status'], 'answered')
        self.assertTrue(BRepCheck_Analyzer(shape).IsValid(), name)
        return manifest, shape

    def volume(self, shape):
        # The adaptive form: the fixed rule misses a closed spline band.
        properties = GProp_GProps()
        BRepGProp.VolumePropertiesGK_s(shape, properties, 1e-6, False, False)
        return properties.Mass()

    def counts(self, shape):
        return tuple(len(unique(shape, kind))
                     for kind in (TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX))

    def test_rings(self):
        with tempfile.TemporaryDirectory(prefix='catallax-ring-fillet-') as temporary:
            root = Path(temporary)
            for name, kind, narrow, change, counts in FIXTURES:
                with self.subTest(name):
                    manifest, shape = self.rounded(root, name, ring(kind, narrow, True))
                    self.assertEqual(self.counts(shape), counts)
                    _, plain = self.read(root, f'plain-{name}', ring(kind, narrow, False))
                    self.assertAlmostEqual(self.volume(shape) - self.volume(plain), change,
                                           delta=2e-3)
                    roles = [region['role'] for region in manifest['certificate']['regions']]
                    self.assertEqual(roles.count('rolling-saddle'), 1)

    def test_a_through_bore(self):
        with tempfile.TemporaryDirectory(prefix='catallax-ring-fillet-') as temporary:
            root = Path(temporary)
            manifest, shape = self.rounded(root, 'through', ring('subtract', THROUGH, True))
            self.assertEqual(self.counts(shape)[0], 6)
            _, plain = self.read(root, 'plain-through', ring('subtract', THROUGH, False))
            self.assertAlmostEqual(self.volume(shape) - self.volume(plain), -2 * 11.64886,
                                   delta=4e-3)
            kinds = manifest['certificate']['retained_target']
            self.assertEqual(kinds['kind'], 'rolling-saddle-junction')
            self.assertEqual(len(kinds['saddle_faces']), 2)

    def test_a_mirrored_and_a_turned_ring(self):
        recipe = ring('union', WALL, True)
        first = {'id': 'mirrored', 'op': 'mirror', 'body': 'round', 'keep_original': False,
                 'plane': {'point': ['20', '0', '0'], 'normal': ['3', '4', '0']}}
        # Two mirrors make a turn, which keeps the chart's handedness in a general pose.
        second = {'id': 'turned', 'op': 'mirror', 'body': 'mirrored', 'keep_original': False,
                  'plane': {'point': ['1', '2', '3'], 'normal': ['1', '2', '2']}}
        with tempfile.TemporaryDirectory(prefix='catallax-ring-fillet-') as temporary:
            root = Path(temporary)
            _, original = self.rounded(root, 'original', recipe)
            for result, steps in (('mirrored', [first]), ('turned', [first, second])):
                with self.subTest(result):
                    _, shape = self.rounded(root, result, {
                        **recipe, 'name': result, 'result': result,
                        'steps': recipe['steps'] + steps})
                    self.assertEqual(self.counts(shape), FIXTURES[0][-1])
                    self.assertAlmostEqual(self.volume(shape), self.volume(original), delta=1e-3)


if __name__ == '__main__':
    unittest.main()
