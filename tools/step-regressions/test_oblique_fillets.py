"""A rolling-pipe fillet on an oblique ellipse, read without ShapeFix or sewing."""
import json
import math
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
# A half turn about y of 2 atan(7/24): the plane normal is (336, 0, 527) / 625.
TILT = ['24/25', '0', '7/25', '0']


def cylinder(name, height, z):
    return {'id': name, 'op': 'cylinder', 'radius': '5', 'height': height,
            'pose': {'rotation': LEVEL, 'translation': ['0', '0', z]}}


def slab(at, distance):
    corners = [('s', '30', '-30'), ('e', '30', '30'), ('n', '-30', '30'), ('w', '-30', '-30')]
    return [
        {'id': 'cut-section', 'op': 'sketch',
         'plane': {'rotation': TILT, 'translation': ['0', '0', at]},
         'loops': [[{'kind': 'line', 'id': name, 'to': [x, y]} for name, x, y in corners]]},
        {'id': 'cut', 'op': 'extrude', 'profile': 'cut-section', 'distance': distance},
    ]


def joined(kind, target, tool):
    return {'id': 'joined', 'op': 'boolean', 'kind': kind, 'target': target, 'tool': tool}


def rounded(steps, radius, role):
    edge = {'between': [{'step': 'wall', 'role': {'side': {'segment': 'wall'}}},
                        {'step': 'cut', 'role': role}]}
    return {'format': 'catallax-recipe/3', 'name': 'oblique round', 'result': 'round',
            'steps': steps + [{'id': 'round', 'op': 'fillet', 'body': 'joined',
                               'radius': radius, 'edges': [edge]}]}


RIM = [cylinder('wall', '16', '0')] + slab('8', '30') + [joined('subtract', 'wall', 'cut')]
BOSS = [cylinder('wall', '22', '-10')] + slab('-25', '527/25') + [joined('union', 'wall', 'cut')]
BORE = [cylinder('wall', '100', '-60')] + slab('-25', '527/25') + [joined('subtract', 'cut', 'wall')]
# The volume the blend moves, from the section integral in the kernel's tests.
FIXTURES = [
    ('rim-1', RIM, '1', 'start', (4, 5, 3), 8 * 25 * math.pi - 9.8242),
    ('rim-2', RIM, '2', 'start', (4, 5, 3), 8 * 25 * math.pi - 37.6296),
]


class ObliqueFilletImports(unittest.TestCase):
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
        self.assertEqual(manifest['marks'], [])
        self.assertEqual(manifest['result']['faces_by_kind'].get('rolling-pipe', 0),
                         int(recipe['result'] != 'joined'))
        shape = raw_solid(Path(case['path']) / 'answer.step')
        self.assertTrue(BRepCheck_Analyzer(shape).IsValid(), case['path'])
        return shape

    def volume(self, shape):
        # The adaptive form: the fixed rule misses a closed spline band by 2 here.
        properties = GProp_GProps()
        BRepGProp.VolumePropertiesGK_s(shape, properties, 1e-6, False, False)
        return properties.Mass()

    def test_rims(self):
        with tempfile.TemporaryDirectory(prefix='catallax-oblique-fillet-') as temporary:
            root = Path(temporary)
            for name, steps, radius, role, counts, volume in FIXTURES:
                with self.subTest(name):
                    shape = self.read(root, name, rounded(steps, radius, role))
                    self.assertEqual(tuple(len(unique(shape, kind)) for kind in
                                           (TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX)), counts)
                    self.assertAlmostEqual(self.volume(shape), volume, delta=2e-4)

    def test_a_mirrored_and_a_turned_rim(self):
        rounded_rim = rounded(RIM, '1', 'start')
        first = {'id': 'mirrored', 'op': 'mirror', 'body': 'round', 'keep_original': False,
                 'plane': {'point': ['20', '0', '0'], 'normal': ['3', '4', '0']}}
        # Two mirrors make a turn, which keeps the chart's handedness in a general pose.
        second = {'id': 'turned', 'op': 'mirror', 'body': 'mirrored', 'keep_original': False,
                  'plane': {'point': ['1', '2', '3'], 'normal': ['1', '2', '2']}}
        with tempfile.TemporaryDirectory(prefix='catallax-oblique-fillet-') as temporary:
            root = Path(temporary)
            for result, steps in (('mirrored', [first]), ('turned', [first, second])):
                with self.subTest(result):
                    recipe = {**rounded_rim, 'name': result, 'result': result,
                              'steps': rounded_rim['steps'] + steps}
                    shape = self.read(root, result, recipe)
                    self.assertEqual(tuple(len(unique(shape, kind)) for kind in
                                           (TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX)), (4, 5, 3))
                    self.assertAlmostEqual(self.volume(shape), FIXTURES[0][-1], delta=2e-4)

    def test_a_boss_and_its_bore_change_the_volume_oppositely(self):
        with tempfile.TemporaryDirectory(prefix='catallax-oblique-fillet-') as temporary:
            root = Path(temporary)
            boss = self.read(root, 'boss', rounded(BOSS, '1', 'end'))
            bore = self.read(root, 'bore', rounded(BORE, '1', 'end'))
            self.assertEqual(tuple(len(unique(boss, kind)) for kind in
                                   (TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX)), (9, 17, 11))
            self.assertEqual(tuple(len(unique(bore, kind)) for kind in
                                   (TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX)), (8, 17, 11))
            plain_boss = self.read(root, 'plain-boss', {
                'format': 'catallax-recipe/3', 'name': 'boss', 'result': 'joined', 'steps': BOSS})
            plain_bore = self.read(root, 'plain-bore', {
                'format': 'catallax-recipe/3', 'name': 'bore', 'result': 'joined', 'steps': BORE})
            self.assertAlmostEqual(self.volume(boss) - self.volume(plain_boss), 10.6578, delta=2e-4)
            self.assertAlmostEqual(self.volume(bore) - self.volume(plain_bore), -10.6578, delta=2e-4)


if __name__ == '__main__':
    unittest.main()
