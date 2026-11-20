"""Bend and twist STEP read raw by OpenCASCADE, healing off.

The Rust test `bend_twist` writes the files. A twist keeps volume exactly
(its Jacobian determinant is one), and so does this bend, because the plate's
mid-plane is the neutral plane: the determinant u / radius averages to one
across a thickness centred on it. The probes sit where only the warped shape
has material, or where only the unwarped one has it.
"""
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


def twisted(x, y, z, degrees):
    turn = math.radians(degrees) * z / 100
    return [x * math.cos(turn) - y * math.sin(turn), x * math.sin(turn) + y * math.cos(turn), z]


def bent(x, y, z):
    """The test's bend: axis z through the origin, toward +y, radius 20."""
    turn = -x / 20
    return [-y * math.sin(turn), y * math.cos(turn), z]


CASES = {
    'twisted.step': (10 * 10 * 100, [
        ([0, 0, 50], TopAbs_IN),
        # Half way along, the 90 degree twist has turned the square by 45
        # degrees: its corner diagonal lies on the x axis.
        ([6.5, 0, 50], TopAbs_IN),
        (twisted(4.5, 4.5, 50, 90), TopAbs_IN),
        ([4.5, 4.5, 50], TopAbs_OUT),
        ([8, 0, 50], TopAbs_OUT),
    ]),
    'twisted-edited.step': (10 * 10 * 100, [
        ([0, 0, 50], TopAbs_IN),
        (twisted(4.5, 4.5, 50, 120), TopAbs_IN),
        (twisted(4.5, 4.5, 50, 90), TopAbs_OUT),
    ]),
    'bent.step': (30 * 2 * 20, [
        ([0, 20, 10], TopAbs_IN),
        (bent(-10, 20, 10), TopAbs_IN),
        (bent(14, 20.9, 1), TopAbs_IN),
        # The unbent plate's end lies outside the bend's 21 mm outer radius.
        ([-10, 20, 10], TopAbs_OUT),
        ([0, 22, 10], TopAbs_OUT),
    ]),
}


class BendTwistOracles(unittest.TestCase):
    def test_bent_and_twisted_solids_read_raw_as_valid_with_their_volume(self):
        with tempfile.TemporaryDirectory(prefix='catallax-bend-twist-') as directory:
            proof = Path(directory)
            environment = dict(os.environ)
            environment['CATALLAX_DEFORM_PROOF'] = str(proof)
            subprocess.run([
                'cargo', 'test', '--release', '-p', 'catallax-features', '--test', 'bend_twist',
                '--', '--nocapture',
            ], cwd=ROOT, env=environment, check=True, timeout=1200)
            for name, (volume, probes) in CASES.items():
                with self.subTest(file=name):
                    shapes = raw_shapes(proof / name)
                    self.assertEqual(len(shapes), 1)
                    shape = shapes[0]
                    self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
                    self.assertEqual(count(shape, TopAbs_SOLID), 1)
                    self.assertEqual(count(shape, TopAbs_FACE), 6)
                    properties = GProp_GProps()
                    error = BRepGProp.VolumePropertiesGK_s(
                        shape, properties, Eps=1e-10, OnlyClosed=True,
                    )
                    self.assertTrue(math.isfinite(error))
                    # Fitted faces lie within 1e-6 mm of the warp: under 1e-3
                    # mm^3 over these areas.
                    print(f'{name}: volume {properties.Mass()!r} mm^3, expected {volume}')
                    self.assertAlmostEqual(properties.Mass(), volume, delta=1e-3)
                    for point, expected in probes:
                        actual = BRepClass3d_SolidClassifier(shape, gp_Pnt(*point), 1e-7).State()
                        self.assertEqual(actual, expected, (name, point))


if __name__ == '__main__':
    unittest.main()
