"""Unrepaired STEP imports of the finite open-base cylinder/sphere shell."""
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
from OCP.GProp import GProp_GProps
from OCP.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane, GeomAbs_Sphere
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_FACE, TopAbs_IN, TopAbs_OUT
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get('CATALLAX_VERIFY_BINARY', ROOT / 'target/release/catallax-verify')).resolve()
DRIVER = ROOT / '.claude/skills/verify-catallax/scripts/control-catallax.py'
CHECK = runpy.run_path(str(DRIVER.with_name('check-step.py')))['check']
FIXTURE = ROOT / 'crates/catallax-features/tests/fixtures/open-spherical-cap/base-open.json'
FACE = getattr(TopoDS, 'Face_s', None) or TopoDS.Face


class OpenSphericalCapImports(unittest.TestCase):
    def test_opening_edit_finite_faces_volume_area_and_material_without_repair(self):
        with tempfile.TemporaryDirectory(prefix='catallax-open-spherical-cap-') as temporary:
            for thickness in [0.5, 0.75]:
                with self.subTest(thickness=thickness):
                    recipe = copy.deepcopy(json.loads(FIXTURE.read_text()))
                    recipe['steps'][3]['thickness'] = str(thickness)
                    directory = Path(temporary) / str(thickness)
                    directory.mkdir()
                    source = directory / 'recipe.json'
                    source.write_text(json.dumps(recipe))
                    out = directory / 'run'
                    completed = subprocess.run(
                        [str(DRIVER), 'run', str(source), '--source', 'recipe',
                         '--binary', str(BINARY), '--max-work', '20000000',
                         '--out', str(out), '--fresh', '--json'],
                        check=True, capture_output=True, text=True)
                    self.assertEqual(json.loads(completed.stdout)['outcomes'], {'answered': 1})
                    case = json.loads((out / 'run.json').read_text())['cases'][0]
                    manifest = case['response']['manifest']
                    self.assertTrue(manifest['exact'])
                    self.assertEqual(manifest['marks'], [])
                    step = Path(case['path']) / 'answer.step'
                    self.assertTrue(step.is_file())
                    checked = CHECK(str(step), 5)
                    self.assertEqual(checked['status'], 'passed', checked)
                    reader = STEPControl_Reader()
                    self.assertEqual(reader.ReadFile(str(step)), IFSelect_RetDone)
                    self.assertGreater(reader.TransferRoots(), 0)
                    shape = reader.OneShape()
                    properties = GProp_GProps()
                    error = BRepGProp.VolumePropertiesGK_s(shape, properties, 1e-12,
                                                         True, False, False, False, False)
                    self.assertGreaterEqual(error, 0)
                    self.assertLessEqual(error, 1e-10)
                    a, radius = 3 - thickness, 5 - thickness
                    join = math.sqrt(radius * radius - a * a)
                    # Independent disk integrals of the source and cavity.
                    expected_volume = math.pi * (122 / 3 - 2 * (radius**3 - join**3) / 3)
                    self.assertAlmostEqual(properties.Mass(), expected_volume, delta=1e-8)
                    area = GProp_GProps()
                    BRepGProp.SurfaceProperties_s(shape, area)
                    expected_area = math.pi * (34 + 9 - a*a + 2*a*join
                                               + 2*radius*(radius-join))
                    self.assertAlmostEqual(area.Mass(), expected_area, delta=1e-7)
                    for xyz, expected in [((0, 0, 0.1), TopAbs_OUT),
                                          ((2, 0, 0.1), TopAbs_OUT),
                                          ((2.75, 0, 0.1), TopAbs_IN),
                                          ((2.75, 0, -0.1), TopAbs_OUT),
                                          ((2.75, 0, 4), TopAbs_IN),
                                          ((0, 0, 4.75), TopAbs_IN),
                                          ((0, 0, 5.1), TopAbs_OUT)]:
                        self.assertEqual(BRepClass3d_SolidClassifier(
                            shape, gp_Pnt(*xyz), 1e-7).State(), expected, xyz)
                    counts = {GeomAbs_Plane: 0, GeomAbs_Cylinder: 0, GeomAbs_Sphere: 0}
                    explorer = TopExp_Explorer(shape, TopAbs_FACE)
                    while explorer.More():
                        kind = BRepAdaptor_Surface(FACE(explorer.Current())).GetType()
                        self.assertIn(kind, counts)
                        counts[kind] += 1
                        explorer.Next()
                    self.assertEqual(counts, {GeomAbs_Plane: 1,
                                              GeomAbs_Cylinder: 2,
                                              GeomAbs_Sphere: 2})


if __name__ == '__main__':
    unittest.main()
