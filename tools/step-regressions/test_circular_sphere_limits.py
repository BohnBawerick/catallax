"""Unrepaired circular sphere caps: finite material and independent integrals."""
import copy
import json
import math
import os
from pathlib import Path
import runpy
import subprocess
import tempfile
import unittest

from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepGProp import BRepGProp
from OCP.GeomAbs import GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Sphere
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_FACE, TopAbs_IN, TopAbs_OUT, TopAbs_ON
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
CHECK = runpy.run_path(str(DRIVER.with_name('check-step.py')))['check']
SOURCE = ROOT / 'crates/catallax-features/tests/fixtures/circular-sphere-limits/sphere-entry-next.json'
FACE = getattr(TopoDS, 'Face_s', None) or TopoDS.Face


class CircularSphereImports(unittest.TestCase):
    def test_both_modes_branches_edits_reversal_and_placement_without_repair(self):
        with tempfile.TemporaryDirectory(prefix='catallax-circular-sphere-') as temporary:
            for mode in ['extrude_next', 'extrude_to']:
                for exit_cap in [False, True]:
                    for configuration in ['base', 'radius', 'sphere', 'height', 'reverse', 'placed']:
                        with self.subTest(mode=mode, exit=exit_cap, configuration=configuration):
                            recipe = copy.deepcopy(json.loads(SOURCE.read_text()))
                            radius = 4 if configuration == 'radius' else 3
                            sphere = 6 if configuration == 'sphere' else 5
                            start = (11 if configuration == 'height' else 10) if exit_cap else (
                                -1 if configuration == 'height' else 0)
                            forward = -1 if configuration == 'reverse' else 1
                            if forward == -1:
                                start = 20 - start
                            recipe['steps'][0]['loops'][0][0]['radius'] = str(radius)
                            recipe['steps'][0]['plane']['translation'][2] = str(start)
                            recipe['steps'][1]['radius'] = str(sphere)
                            operation = recipe['steps'][2]
                            operation['op'] = mode
                            operation['direction'] = ['0', '0', str(2 * forward)]
                            if mode == 'extrude_to':
                                operation['limit'] = {'face': {'step': 'target', 'role': {
                                    'side': {'segment': 'wall'}}}}

                            def world(x, y, z):
                                if configuration == 'placed':
                                    return (7+x, 11+3*y/5-4*z/5, -3+4*y/5+3*z/5)
                                return (x, y, z)

                            if configuration == 'placed':
                                for index, key, height in [(0, 'plane', start), (1, 'pose', 10)]:
                                    recipe['steps'][index][key] = {
                                        'rotation': ['2', '1', '0', '0'],
                                        'translation': [str(x) for x in world(0, 0, height)]}
                                operation['direction'] = ['0', '-8/5', '6/5']
                            name = f'{mode}-{exit_cap}-{configuration}'
                            directory = Path(temporary) / name
                            directory.mkdir()
                            source = directory / 'recipe.json'
                            source.write_text(json.dumps(recipe))
                            out = directory / 'run'
                            completed = subprocess.run(
                                [str(DRIVER), 'run', str(source), '--source', 'recipe',
                                 '--binary', str(BINARY), '--out', str(out), '--fresh', '--json'],
                                check=True, capture_output=True, text=True)
                            self.assertEqual(json.loads(completed.stdout)['outcomes'], {'answered': 1})
                            case = json.loads((out / 'run.json').read_text())['cases'][0]
                            manifest = case['response']['manifest']
                            self.assertTrue(manifest['exact'])
                            self.assertEqual(manifest['marks'], [])
                            path = Path(case['path']) / 'answer.step'
                            # The spherical face contains its pole, which moves
                            # onto an export-only cap: three model faces, four
                            # exported faces.
                            self.assertEqual(CHECK(str(path), 4)['status'], 'passed')
                            reader = STEPControl_Reader()
                            self.assertEqual(reader.ReadFile(str(path)), IFSelect_RetDone)
                            self.assertGreater(reader.TransferRoots(), 0)
                            shape = reader.OneShape()
                            properties = GProp_GProps()
                            error = BRepGProp.VolumePropertiesGK_s(shape, properties, 1e-12,
                                                                 True, False, False, False, False)
                            self.assertGreaterEqual(error, 0)
                            self.assertLessEqual(error, 1e-10)
                            h = math.sqrt(sphere*sphere-radius*radius)
                            integral = 2*math.pi*(sphere**3-h**3)/3
                            start_axial = forward*(start-10)
                            expected = -start_axial*math.pi*radius*radius + (
                                integral if exit_cap else -integral)
                            self.assertAlmostEqual(properties.Mass(), expected, delta=1e-8)
                            rim_height = -start_axial + (h if exit_cap else -h)
                            for (u, v, t), state in [((0, 0, 1), TopAbs_IN),
                                                     ((0, 0, -0.1), TopAbs_OUT),
                                                     ((radius, 0, 1), TopAbs_ON),
                                                     ((radius+0.1, 0, 1), TopAbs_OUT),
                                                     ((radius, radius, 1), TopAbs_OUT),
                                                     ((0, 0, -start_axial+sphere+1), TopAbs_OUT),
                                                     ((radius, 0, rim_height), TopAbs_ON),
                                                     ((0, 0, -start_axial+(
                                                         sphere if exit_cap else -sphere)), TopAbs_ON)]:
                                p = world(u, v, start + forward*t)
                                self.assertEqual(BRepClass3d_SolidClassifier(shape, gp_Pnt(*p), 1e-7).State(), state, p)
                            kinds = []
                            explorer = TopExp_Explorer(shape, TopAbs_FACE)
                            while explorer.More():
                                kinds.append(BRepAdaptor_Surface(FACE(explorer.Current())).GetType())
                                explorer.Next()
                            # The pole cap is a second spherical face.
                            self.assertCountEqual(kinds, [GeomAbs_Plane, GeomAbs_Cylinder,
                                                          GeomAbs_Sphere, GeomAbs_Sphere])


if __name__ == '__main__':
    unittest.main()
