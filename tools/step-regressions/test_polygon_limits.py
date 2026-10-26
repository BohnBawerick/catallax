"""Independent STEP imports of finite convex-polygon first-hit caps."""
import copy
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import re
import runpy
import subprocess
import tempfile
import unittest

from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Curve2d, BRepAdaptor_Surface
from OCP.BRepGProp import BRepGProp
from OCP.GeomAbs import GeomAbs_Circle, GeomAbs_Sphere
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
CHECK = runpy.run_path(str(DRIVER.with_name('check-step.py')))['check']
FIXTURES = ROOT / 'crates/catallax-features/tests/fixtures/polygon-limits'


def spherical_circle_bounds(step_text):
    """Exact radial carrier bounds; neither trim nor import mass bounds."""
    entities = {int(i): value for i, value in re.findall(r'^#(\d+)=(.*);$', step_text, re.M)}

    def refs(value):
        return [int(i) for i in re.findall(r'#(\d+)', value)]

    def xyz(value):
        return [Fraction(x) for x in re.search(r',\(([^()]*)\)\)', value)[1].split(',')]

    def dot(a, b):
        return sum(x * y for x, y in zip(a, b))

    def sqrt_upper(value):
        scale = 1 << 160
        return Fraction(math.isqrt(value.numerator * scale * scale // value.denominator) + 1, scale)

    sphere = next(value for value in entities.values() if value.startswith('SPHERICAL_SURFACE('))
    radius = Fraction(re.search(r',([^,]+)\)$', sphere)[1])
    center = xyz(entities[refs(entities[refs(sphere)[0]])[0]])
    bounds = []
    for circle in entities.values():
        if not circle.startswith('CIRCLE('):
            continue
        r = Fraction(re.search(r',([^,]+)\)$', circle)[1])
        placement = entities[refs(circle)[0]]
        c, a = [xyz(entities[i]) for i in refs(placement)[:2]]
        d = [x - y for x, y in zip(c, center)]
        transverse_sq = dot(d, d) - dot(d, a) ** 2 / dot(a, a)
        # For every point on this exact circle, |F| is bounded by the
        # constant residual plus 2*r*|d_perpendicular|. Divide by R because
        # |distance-R| = |F|/(distance+R). All operations are rational.
        bounds.append((abs(dot(d, d) + r * r - radius * radius)
                       + 2 * r * sqrt_upper(transverse_sq)) / radius)
    return bounds


class PolygonLimitImports(unittest.TestCase):
    def test_both_modes_and_branches_have_valid_closed_solids_and_exact_volume_identity(self):
        with tempfile.TemporaryDirectory(prefix='catallax-polygon-limits-') as temporary:
            directory = Path(temporary)
            for fixture in sorted(FIXTURES.glob('*.json')):
                source = json.loads(fixture.read_text())
                points = [[Fraction(x) for x in s['to']] for s in source['steps'][0]['loops'][0]]
                area = abs(sum(a[0] * b[1] - a[1] * b[0]
                               for a, b in zip(points, points[1:] + points[:1]))) / 2
                axis = 2 if source['steps'][1]['op'] == 'sphere' else 0
                for mode in ['extrude_next', 'extrude_to']:
                    volumes = []
                    for branch in ['entry', 'exit']:
                        with self.subTest(fixture=fixture.stem, mode=mode, branch=branch):
                            recipe = copy.deepcopy(source)
                            if branch == 'exit':
                                recipe['steps'][0]['plane']['translation'][axis] = '5'
                            recipe['steps'][2]['op'] = mode
                            direction = ['0', '0', '0']
                            direction[axis] = '2'
                            recipe['steps'][2]['direction'] = direction
                            if mode == 'extrude_to':
                                recipe['steps'][2]['limit'] = {
                                    'face': {'step': 'target', 'role': {'side': {'segment': 'wall'}}}}
                            name = f'{fixture.stem}-{mode}-{branch}'
                            path = directory / (name + '.json')
                            path.write_text(json.dumps(recipe))
                            out = directory / name
                            completed = subprocess.run(
                                [str(DRIVER), 'run', str(path), '--source', 'recipe',
                                 '--binary', str(BINARY), '--out', str(out), '--fresh', '--json'],
                                check=True, capture_output=True, text=True)
                            self.assertEqual(json.loads(completed.stdout)['outcomes'], {'answered': 1})
                            case = json.loads((out / 'run.json').read_text())['cases'][0]
                            checked = CHECK(str(Path(case['path']) / 'answer.step'), len(points) + 2)
                            self.assertEqual(checked['status'], 'passed', checked)
                            # The shared import checker reports nonadaptive
                            # diagnostic properties. Request adaptive native
                            # integration here for the independent volume oracle.
                            reader = STEPControl_Reader()
                            self.assertEqual(reader.ReadFile(str(Path(case['path']) / 'answer.step')),
                                             IFSelect_RetDone)
                            self.assertGreater(reader.TransferRoots(), 0)
                            properties = GProp_GProps()
                            error = BRepGProp.VolumePropertiesGK_s(reader.OneShape(), properties,
                                                                 1e-12, True, False, False, False, False)
                            self.assertGreaterEqual(error, 0)
                            self.assertLessEqual(error, 1e-10)
                            volume = properties.Mass()
                            volumes.append(volume)
                            if fixture.stem == 'cylinder-triangle' and branch == 'entry':
                                exact = 2.5 - (11 * math.sqrt(15) / 4 + 4 * math.asin(0.25) - 32 / 3)
                                self.assertAlmostEqual(volume, exact, delta=1e-8)
                    # Independent import properties are a numerical sanity
                    # check; certified exact enclosures are tested in Rust.
                    with self.subTest(fixture=fixture.stem, mode=mode, check='volume_identity'):
                        self.assertAlmostEqual(sum(volumes), float(5 * area), delta=1e-7)

    def test_spherical_triangle_and_rotated_square_keep_exact_carriers_and_imported_trims(self):
        with tempfile.TemporaryDirectory(prefix='catallax-spherical-polygon-trims-') as temporary:
            directory = Path(temporary)
            for profile in ['triangle', 'rotated']:
                source = json.loads((FIXTURES / f'sphere-{profile}.json').read_text())
                sides = len(source['steps'][0]['loops'][0])
                for mode in ['extrude_next', 'extrude_to']:
                    for branch in ['entry', 'exit']:
                        with self.subTest(profile=profile, mode=mode, branch=branch):
                            recipe = copy.deepcopy(source)
                            recipe['steps'][2]['op'] = mode
                            recipe['steps'][2]['direction'] = ['0', '0', '2']
                            if branch == 'exit':
                                recipe['steps'][0]['plane']['translation'][2] = '5'
                            if mode == 'extrude_to':
                                recipe['steps'][2]['limit'] = {
                                    'face': {'step': 'target', 'role': {'side': {'segment': 'wall'}}}}
                            name = f'{profile}-{mode}-{branch}'
                            path = directory / (name + '.json')
                            path.write_text(json.dumps(recipe))
                            out = directory / name
                            subprocess.run([str(DRIVER), 'run', str(path), '--source', 'recipe',
                                            '--binary', str(BINARY), '--out', str(out), '--fresh', '--json'],
                                           check=True, capture_output=True, text=True)
                            case = json.loads((out / 'run.json').read_text())['cases'][0]
                            self.assertEqual(case['response']['manifest']['status'], 'answered')
                            step = Path(case['path']) / 'answer.step'
                            bounds = spherical_circle_bounds(step.read_text())
                            self.assertEqual(len(bounds), sides)
                            self.assertTrue(all(b < Fraction(1, 10**15) for b in bounds))
                            self.assertEqual(CHECK(str(step), sides + 2)['status'], 'passed')
                            reader = STEPControl_Reader()
                            self.assertEqual(reader.ReadFile(str(step)), IFSelect_RetDone)
                            self.assertGreater(reader.TransferRoots(), 0)
                            faces = TopExp_Explorer(reader.OneShape(), TopAbs_FACE)
                            spherical_faces = 0
                            while faces.More():
                                face = TopoDS.Face(faces.Current())
                                surface = BRepAdaptor_Surface(face)
                                if surface.GetType() == GeomAbs_Sphere:
                                    spherical_faces += 1
                                    edges = TopExp_Explorer(face, TopAbs_EDGE)
                                    edge_count = 0
                                    while edges.More():
                                        edge = TopoDS.Edge(edges.Current())
                                        self.assertFalse(BRep_Tool.Degenerated_s(edge))
                                        curve = BRepAdaptor_Curve(edge)
                                        self.assertEqual(curve.GetType(), GeomAbs_Circle)
                                        pcurve = BRepAdaptor_Curve2d(edge, face)
                                        # Sampled importer sanity at the existing file
                                        # uncertainty, not a certified trim-distance bound.
                                        for i in range(129):
                                            t = curve.FirstParameter() + (curve.LastParameter()
                                                - curve.FirstParameter()) * i / 128
                                            uv = pcurve.Value(t)
                                            lifted = surface.Value(uv.X(), uv.Y())
                                            self.assertLessEqual(curve.Value(t).Distance(lifted), 1e-7)
                                        edge_count += 1
                                        edges.Next()
                                    self.assertEqual(edge_count, sides)
                                faces.Next()
                            self.assertEqual(spherical_faces, 1)


if __name__ == '__main__':
    unittest.main()
