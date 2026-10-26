"""Complete periodic bands, including a reversed cylindrical hole wall."""
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


class PeriodicBandImports(unittest.TestCase):
    def check_recipe(self, recipe, expected, volume, root, name):
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
        self.assertTrue(manifest['exact'], manifest)
        self.assertEqual(manifest['marks'], [])
        shape = raw_solid(Path(case['path']) / 'answer.step')
        self.assertTrue(BRepCheck_Analyzer(shape).IsValid(), case['path'])
        self.assertEqual(tuple(len(unique(shape, kind)) for kind in
                               (TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX)), expected)
        properties = GProp_GProps()
        BRepGProp.VolumeProperties_s(shape, properties)
        self.assertAlmostEqual(properties.Mass(), volume, delta=volume * 1e-10)

    def test_cylinders_and_both_frustum_directions(self):
        poses = [
            {'rotation': ['1', '0', '0', '0'], 'translation': ['0', '0', '0']},
            {'rotation': ['1', '2', '3', '4'], 'translation': ['7/3', '-11/7', '13/5']},
        ]
        fixtures = [
            ('cylinder', {'radius': '1', 'height': '2'}, 2 * math.pi),
            ('cone', {'base_radius': '3', 'top_radius': '1', 'height': '4'}, 52 * math.pi / 3),
            ('cone', {'base_radius': '1', 'top_radius': '3', 'height': '4'}, 52 * math.pi / 3),
        ]
        with tempfile.TemporaryDirectory(prefix='catallax-periodic-band-') as temporary:
            root = Path(temporary)
            for index, (operation, dimensions, volume) in enumerate(fixtures):
                for posed, pose in enumerate(poses):
                    with self.subTest(dimensions=dimensions, pose=posed):
                        recipe = {'format': 'catallax-recipe/3', 'name': 'full band',
                                  'steps': [{'id': 'body', 'op': operation, **dimensions, 'pose': pose}],
                                  'result': 'body'}
                        self.check_recipe(recipe, (3, 3, 2), volume, root, f'band-{index}-{posed}')

    def test_complete_cylindrical_hole_wall(self):
        pose = lambda x, y, z: {'rotation': ['1', '0', '0', '0'], 'translation': [x, y, z]}
        recipe = {'format': 'catallax-recipe/3', 'name': 'through hole', 'steps': [
            {'id': 'base', 'op': 'block', 'size': ['4', '4', '2'], 'pose': pose('0', '0', '0')},
            {'id': 'tool', 'op': 'cylinder', 'radius': '1', 'height': '4', 'pose': pose('2', '2', '-1')},
            {'id': 'hole', 'op': 'boolean', 'kind': 'subtract', 'target': 'base', 'tool': 'tool'}],
            'result': 'hole'}
        with tempfile.TemporaryDirectory(prefix='catallax-periodic-hole-') as temporary:
            self.check_recipe(recipe, (7, 15, 10), 32 - 2 * math.pi, Path(temporary), 'hole')

    def test_complete_band_on_an_enclosed_void_shell(self):
        pose = lambda z: {'rotation': ['1', '0', '0', '0'], 'translation': ['0', '0', z]}
        recipe = {'format': 'catallax-recipe/3', 'name': 'closed cylindrical void', 'steps': [
            {'id': 'base', 'op': 'cylinder', 'radius': '3', 'height': '4', 'pose': pose('0')},
            {'id': 'tool', 'op': 'cylinder', 'radius': '1', 'height': '2', 'pose': pose('1')},
            {'id': 'void', 'op': 'boolean', 'kind': 'subtract', 'target': 'base', 'tool': 'tool'}],
            'result': 'void'}
        with tempfile.TemporaryDirectory(prefix='catallax-periodic-void-') as temporary:
            self.check_recipe(recipe, (6, 6, 4), 34 * math.pi, Path(temporary), 'void')

    def test_algebraic_countersink_bands_at_existing_allowance(self):
        with tempfile.TemporaryDirectory(prefix='catallax-periodic-countersink-') as temporary:
            for angle in (82, 100):
                with self.subTest(angle=angle):
                    recipe = {'format': 'catallax-recipe/3', 'name': 'algebraic countersink',
                              'steps': [
                        {'id': 'stock', 'op': 'block', 'size': ['10', '10', '10'],
                         'pose': {'rotation': ['1', '0', '0', '0'],
                                  'translation': ['-5', '-5', '0']}},
                        {'id': 'drill', 'op': 'hole', 'target': 'stock',
                         'at': [{'id': 'centre', 'point': ['0', '0', '0'],
                                 'direction': ['0', '0', '1']}],
                         'diameter': '2', 'depth': {'blind': '4'}, 'tip': 'flat',
                         'head': {'countersink': {'diameter': '8', 'angle': str(angle)}}}],
                              'result': 'drill'}
                    height = 3 / math.tan(math.radians(angle / 2))
                    self.check_recipe(recipe, (9, 17, 11),
                                      1000 - math.pi * (4 + 6 * height),
                                      Path(temporary), f'countersink-{angle}')


if __name__ == '__main__':
    unittest.main()
