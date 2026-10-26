"""Public regular periodic trims; raw translation remains non-repairing."""
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX, TopAbs_SOLID

from test_loopless_faces import BINARY, DRIVER, raw_solid, unique


def pose(x='0', y='0', z='0', rotation=('1', '0', '0', '0')):
    return {'rotation': list(rotation), 'translation': [x, y, z]}


def converged_volume(shape):
    values = []
    for epsilon in (1e-10, 1e-12, 1e-14):
        properties = GProp_GProps()
        error = BRepGProp.VolumePropertiesGK_s(
            shape, properties, epsilon, True, False, False, False, False)
        if error < 0:
            raise AssertionError(f'volume integration failed: {error}')
        values.append(properties.Mass())
    if abs(values[-1] - values[-2]) > max(1, abs(values[-1])) * 1e-10:
        raise AssertionError(f'volume integration did not converge: {values}')
    return values[-1], max(values) - min(values)


class PartialPeriodicTrimImports(unittest.TestCase):
    def check_intersection(self, primitives, topology, expected_volume):
        recipe = {'format': 'catallax-recipe/3', 'name': 'public periodic trim',
                  'steps': [dict(primitives[0], id='a'), dict(primitives[1], id='b'),
                            {'id': 'body', 'op': 'boolean', 'kind': 'intersect',
                             'target': 'a', 'tool': 'b'}], 'result': 'body'}
        with tempfile.TemporaryDirectory(prefix='catallax-partial-periodic-') as temporary:
            root = Path(temporary)
            source = root / 'recipe.json'
            source.write_text(json.dumps(recipe))
            out = root / 'run'
            completed = subprocess.run(
                [str(DRIVER), 'run', str(source), '--source', 'recipe', '--binary',
                 str(BINARY), '--out', str(out), '--fresh', '--workers', '1', '--json'],
                capture_output=True, text=True, check=True)
            self.assertEqual(json.loads(completed.stdout)['outcomes'], {'answered': 1})
            case = json.loads((out / 'run.json').read_text())['cases'][0]
            manifest = case['response']['manifest']
            self.assertTrue(manifest['exact'], manifest)
            self.assertEqual(manifest['marks'], [])
            path = Path(case['path']) / 'answer.step'
            raw = raw_solid(path)
            reader = STEPControl_Reader()
            self.assertEqual(reader.ReadFile(str(path)), IFSelect_RetDone)
            self.assertGreater(reader.TransferRoots(), 0)
            # This second procedure intentionally includes standard reader
            # processing. Only raw_solid is the no-healing assertion.
            for mode, shape in [('raw', raw), ('standard', reader.OneShape())]:
                with self.subTest(reader=mode):
                    self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
                    self.assertEqual(len(unique(shape, TopAbs_SOLID)), 1)
                    self.assertEqual(tuple(len(unique(shape, kind)) for kind in
                                           (TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX)), topology)
                    volume, spread = converged_volume(shape)
                    area = GProp_GProps()
                    BRepGProp.SurfaceProperties_s(shape, area, 1e-12, False)
                    # Diagnostic tube comparison, not a certification or a
                    # change to any kernel/importer tolerance.
                    self.assertLessEqual(abs(volume - expected_volume),
                                         area.Mass() * 1e-7 + spread)

    def test_partial_cylinder_and_two_spherical_patches(self):
        # Independent adaptive OCCT construction of the same public request.
        self.check_intersection([
            {'op': 'cylinder', 'radius': '1', 'height': '11/4', 'pose': pose(z='-11/8')},
            {'op': 'sphere', 'radius': '7/4', 'pose': pose(x='1/8')}],
            (5, 9, 6), 8.638033145352587)

    def test_touching_branch_vertices_preserve_source_loop_incidence(self):
        # Certified public cylinder/cylinder volume, rounded to binary64.
        self.check_intersection([
            {'op': 'cylinder', 'radius': '2', 'height': '10', 'pose': pose(z='-5')},
            {'op': 'cylinder', 'radius': '1', 'height': '10',
             'pose': pose('-5', '1', '0', ('1', '0', '1', '0'))}],
            (4, 5, 3), 10.165315870570877)

    def test_closed_branch_has_bounded_parameter_coupled_pcurves(self):
        self.check_intersection([
            {'op': 'cylinder', 'radius': '1', 'height': '6', 'pose': pose(z='-3')},
            {'op': 'sphere', 'radius': '2', 'pose': pose(x='2')}],
            (2, 1, 1), 3.100408279750527)

    def test_spherical_lens_at_identity_and_rational_pose(self):
        for first, second in [(pose(), pose(x='6')),
                              (pose('7/3', '-11/7', '13/5', ('1', '2', '3', '4')),
                               pose('-5/3', '17/7', '23/5', ('1', '2', '3', '4')))]:
            with self.subTest(pose=first):
                self.check_intersection([
                    {'op': 'sphere', 'radius': '5', 'pose': first},
                    {'op': 'sphere', 'radius': '5', 'pose': second}],
                    (2, 1, 1), 104 * math.pi / 3)


if __name__ == '__main__':
    unittest.main()
