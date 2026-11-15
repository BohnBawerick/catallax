"""Native numerical cap joins against direct-builder STEP and cylinder material."""
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_FACE, TopAbs_IN, TopAbs_OUT, TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.gp import gp_Pnt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools/checks'))
from raw_step import raw_shapes


def count(shape, kind):
    walk = TopExp_Explorer(shape, kind)
    total = 0
    while walk.More():
        total += 1
        walk.Next()
    return total


class HoleCapOracles(unittest.TestCase):
    def test_skew_cap_shell_and_diameter_edit_keep_material(self):
        with tempfile.TemporaryDirectory(prefix='catallax-hole-caps-') as directory:
            proof = Path(directory)
            environment = dict(os.environ)
            environment['CATALLAX_CAP_PROOF'] = str(proof)
            subprocess.run([
                'cargo', 'test', '--release', '-p', 'catallax-features', '--test', 'hole_caps',
                'shell::skew_caps_sew_edit_transfer_and_export_as_a_numerical_shell',
                '--', '--exact', '--nocapture',
            ], cwd=ROOT, env=environment, check=True, timeout=600)
            center = [23.5, 5, -7]

            def radial(distance):
                return [center[0] + distance / math.sqrt(2),
                        center[1] + distance / math.sqrt(2), center[2]]

            for label, radius in [('initial', 1), ('edited', .75)]:
                with self.subTest(label=label):
                    shapes = raw_shapes(proof / (label + '.step'))
                    self.assertEqual(len(shapes), 1)
                    shape = shapes[0]
                    self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
                    self.assertEqual(count(shape, TopAbs_SOLID), 1)
                    self.assertEqual(count(shape, TopAbs_FACE), 6)
                    properties = GProp_GProps()
                    integration_error = BRepGProp.VolumePropertiesGK_s(
                        shape, properties, Eps=1e-12, OnlyClosed=True, CGFlag=True,
                    )
                    self.assertTrue(math.isfinite(integration_error))
                    self.assertGreaterEqual(integration_error, 0)
                    self.assertLessEqual(integration_error, 1e-10)
                    expected_volume = 12 * math.sqrt(3) * math.pi * radius ** 2
                    self.assertAlmostEqual(properties.Mass(), expected_volume, delta=1e-5)
                    for actual, expected in zip(properties.CentreOfMass().Coord(), center):
                        self.assertAlmostEqual(actual, expected, delta=1e-6)
                    for point, expected in [
                        (center, TopAbs_IN), (radial(.6), TopAbs_IN),
                        (radial(.85), TopAbs_IN if label == 'initial' else TopAbs_OUT),
                        (radial(1.1), TopAbs_OUT), ([16.5, 12, 0], TopAbs_OUT),
                        ([30.5, -2, -14], TopAbs_OUT),
                    ]:
                        actual = BRepClass3d_SolidClassifier(shape, gp_Pnt(*point), 1e-7).State()
                        self.assertEqual(actual, expected, (label, point))
