"""Public regular periodic trims; raw translation bypasses caller-side processing."""
from fractions import Fraction
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


def intersection(primitives, kind='intersect'):
    return {'format': 'catallax-recipe/3', 'name': 'public periodic trim',
            'steps': [dict(primitives[0], id='a'), dict(primitives[1], id='b'),
                      {'id': 'body', 'op': 'boolean', 'kind': kind,
                       'target': 'a', 'tool': 'b'}], 'result': 'body'}


class PartialPeriodicTrimImports(unittest.TestCase):
    def answer(self, recipe, root):
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
        return Path(case['path']), manifest

    def check_intersection(self, primitives, topology, expected_volume, kind='intersect'):
        self.check_recipe(intersection(primitives, kind), topology, expected_volume)

    def check_recipe(self, recipe, topology, expected_volume):
        with tempfile.TemporaryDirectory(prefix='catallax-partial-periodic-') as temporary:
            directory, _ = self.answer(recipe, Path(temporary))
            path = directory / 'answer.step'
            raw = raw_solid(path)
            reader = STEPControl_Reader()
            self.assertEqual(reader.ReadFile(str(path)), IFSelect_RetDone)
            self.assertGreater(reader.TransferRoots(), 0)
            # This second procedure intentionally includes standard reader
            # processing. Only raw_solid bypasses caller-side processing;
            # internal repair inside the reader remains unproven.
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

    def test_cut_vertices_share_the_seam_parameter(self):
        # Public block/cylinder union with a cylinder band on each side. Each
        # band's seam meets two cut vertices. Inverting their rounded points
        # gave u values one binary64 step apart, and raw translation paired
        # that seam's pcurves the wrong way round. Independent OCCT Boolean volume.
        self.check_intersection([
            {'op': 'block', 'size': ['64.834', '68.884', '35.02'],
             'pose': pose('60829/7000', '75771/7000', '462887/7000', ('-5', '-1', '-1', '-8'))},
            {'op': 'cylinder', 'radius': '3.934', 'height': '68.225',
             'pose': pose('-244605379/7000000', '907370213/45500000', '3159971679/45500000',
                          ('-14', '5', '-8', '-25'))}],
            (10, 18, 12), 157588.97339949687, kind='union')

    def test_uncertified_pole_keeps_the_answer(self):
        # The ball's frame puts its poles on the y axis. Its membership search
        # finds no certificate, so the scene keeps its legacy encoding inside
        # the unchanged request allowance instead of refusing.
        with tempfile.TemporaryDirectory(prefix='catallax-partial-periodic-') as temporary:
            _, manifest = self.answer(intersection([
                {'op': 'cylinder', 'radius': '1', 'height': '11/4', 'pose': pose(z='-11/8')},
                {'op': 'sphere', 'radius': '7/4',
                 'pose': pose(x='1/8', rotation=('1', '-1', '0', '0'))}]), Path(temporary))
            self.assertEqual(manifest['step']['added_edges'], 0)

    def test_touching_spherical_trims_keep_legacy_encoding(self):
        # The first intersection leaves spherical trims touching at one
        # vertex. Its partial charts crossed themselves in the parameter
        # plane and failed both readers; the legacy file reads under the
        # standard reader, so this scene keeps its legacy encoding.
        tilt = ('1', '1', '1', '1')
        recipe = {'format': 'catallax-recipe/3', 'name': 'touching spherical trims', 'steps': [
            {'id': 'a', 'op': 'cylinder', 'radius': '3', 'height': '6',
             'pose': pose('-3', rotation=tilt)},
            {'id': 'b', 'op': 'sphere', 'radius': '1', 'pose': pose('0', '3', '1', tilt)},
            {'id': 's1', 'op': 'boolean', 'kind': 'intersect', 'target': 'a', 'tool': 'b'},
            {'id': 'c2', 'op': 'block', 'size': ['2', '2', '2'],
             'pose': pose('-4/5', '-1/2', '-2/3', tilt)},
            {'id': 's2', 'op': 'boolean', 'kind': 'union', 'target': 's1', 'tool': 'c2'},
            {'id': 'c3', 'op': 'cylinder', 'radius': '2', 'height': '4',
             'pose': pose('-9/5', '1', '0', tilt)},
            {'id': 's3', 'op': 'boolean', 'kind': 'subtract', 'target': 's2', 'tool': 'c3'}],
            'result': 's3'}
        with tempfile.TemporaryDirectory(prefix='catallax-partial-periodic-') as temporary:
            directory, _ = self.answer(recipe, Path(temporary))
            self.assertNotIn('PCURVE', (directory / 'answer.step').read_text())

    def test_interior_sphere_pole_moves_onto_an_export_cap(self):
        # An L-shaped prism minus a unit ball centred on its top face. The
        # dimple contains the ball's lower pole, so no chart of its own frame
        # is pole-free. An export-only latitude ring moves that pole onto a cap
        # written with a perpendicular frame. Volume: 22 - 2*pi/3.
        for rotation in [('1', '0', '0', '0'), ('1', '2', '3', '4')]:
            with self.subTest(rotation=rotation):
                self.check_recipe({
                    'format': 'catallax-recipe/3', 'name': 'spherical dimple', 'steps': [
                        {'id': 'a', 'op': 'block', 'size': ['4', '2', '2'], 'pose': pose()},
                        {'id': 'b', 'op': 'block', 'size': ['1', '3', '2'], 'pose': pose(y='2')},
                        {'id': 'l', 'op': 'boolean', 'kind': 'union', 'target': 'a', 'tool': 'b'},
                        {'id': 's', 'op': 'sphere', 'radius': '1',
                         'pose': pose('2', '1', '2', rotation)},
                        {'id': 'd', 'op': 'boolean', 'kind': 'subtract', 'target': 'l',
                         'tool': 's'}],
                    'result': 'd'}, (11, 25, 16), 22 - 2 * math.pi / 3)

    def test_trim_close_to_a_pole_uses_a_perpendicular_frame(self):
        # A radius-two ball cut by the plane 3x + 4y = 5e-15. The trim misses
        # both poles of the ball's frame by 1e-15, so charts on that frame
        # need more pcurve spans than reach the tolerance, and the answer must
        # not lose its exactness. The retry on the perpendicular frame moves the
        # old poles away from the trim. Volume: half the ball, 16*pi/3.
        offset = Fraction(1, 10**15)
        self.check_recipe({
            'format': 'catallax-recipe/3', 'name': 'trim close to a pole', 'steps': [
                {'id': 'a', 'op': 'sphere', 'radius': '2', 'pose': pose()},
                {'id': 'b', 'op': 'block', 'size': ['10', '10', '10'],
                 'pose': pose(str(Fraction(3, 5) * offset + 4), str(Fraction(4, 5) * offset - 3),
                              '-5', ('2', '0', '0', '1'))},
                {'id': 'body', 'op': 'boolean', 'kind': 'intersect', 'target': 'a', 'tool': 'b'}],
            'result': 'body'}, (3, 3, 2), 16 * math.pi / 3)

    def test_pole_on_a_trim_writes_another_sphere_frame(self):
        # One octant of a ball: its spherical face has a corner at the ball's
        # own pole, which no chart or cap of that frame admits. The face is
        # written with another exact frame whose poles lie outside it.
        # Volume: 4*pi/3.
        for rotation in [('1', '0', '0', '0'), ('1', '2', '3', '4')]:
            with self.subTest(rotation=rotation):
                corner = pose('7/3', '-11/7', '13/5', rotation)
                self.check_intersection([
                    {'op': 'block', 'size': ['3', '3', '3'], 'pose': corner},
                    {'op': 'sphere', 'radius': '2', 'pose': corner}],
                    (4, 6, 4), 4 * math.pi / 3)

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
