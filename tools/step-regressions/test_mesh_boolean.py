"""Mesh Booleans and decimation in STEP, read raw by OpenCASCADE, healing off.

The Rust test `mesh_boolean_decimate` writes each file beside its exact
volume. The meshes are a radius 10 sphere at the origin and a radius 9 sphere
at (7.5, 1, 0.5), each with 8 rings and 16 sectors for the Booleans and 16 and
32 for decimation. The probes sit where only one sphere, both, or neither has
material, at least 2 mm inside or outside every facet.
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
from OCP.TopAbs import TopAbs_IN, TopAbs_OUT, TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.gp import gp_Pnt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools/checks'))
from raw_step import raw_shapes

ONLY_A = [-6, 0, 0]
BOTH = [3.5, 0.5, 0.25]
ONLY_B = [13, 1.5, 0.5]
NEITHER = [0, 0, 15]

CASES = {
    'union': [(ONLY_A, TopAbs_IN), (BOTH, TopAbs_IN), (ONLY_B, TopAbs_IN), (NEITHER, TopAbs_OUT)],
    'subtract': [(ONLY_A, TopAbs_IN), (BOTH, TopAbs_OUT), (ONLY_B, TopAbs_OUT), (NEITHER, TopAbs_OUT)],
    'intersect': [(ONLY_A, TopAbs_OUT), (BOTH, TopAbs_IN), (ONLY_B, TopAbs_OUT), (NEITHER, TopAbs_OUT)],
    'decimated': [([0, 0, 0], TopAbs_IN), ([4, 4, 0], TopAbs_IN), (NEITHER, TopAbs_OUT)],
    'union-decimated': [(ONLY_A, TopAbs_IN), (BOTH, TopAbs_IN), (ONLY_B, TopAbs_IN), (NEITHER, TopAbs_OUT)],
}


def count(shape, kind):
    walk = TopExp_Explorer(shape, kind)
    total = 0
    while walk.More():
        total += 1
        walk.Next()
    return total


class MeshBooleanOracles(unittest.TestCase):
    def test_mesh_results_read_raw_as_valid_solids_with_their_exact_volume(self):
        with tempfile.TemporaryDirectory(prefix='catallax-mesh-boolean-') as directory:
            proof = Path(directory)
            environment = dict(os.environ)
            environment['CATALLAX_MESH_PROOF'] = str(proof)
            subprocess.run([
                'cargo', 'test', '--release', '-p', 'catallax-features', '--test',
                'mesh_boolean_decimate', '--', '--nocapture',
            ], cwd=ROOT, env=environment, check=True, timeout=1800)
            for name, probes in CASES.items():
                with self.subTest(file=name):
                    shapes = raw_shapes(proof / f'{name}.step')
                    self.assertEqual(len(shapes), 1)
                    shape = shapes[0]
                    self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
                    self.assertEqual(count(shape, TopAbs_SOLID), 1)
                    properties = GProp_GProps()
                    error = BRepGProp.VolumePropertiesGK_s(
                        shape, properties, Eps=1e-10, OnlyClosed=True,
                    )
                    self.assertTrue(math.isfinite(error))
                    volume = float((proof / f'{name}.volume').read_text())
                    print(f'{name}: volume {properties.Mass()!r} mm^3, exact {volume!r}')
                    # Planar faces: the file rounds each vertex to binary64,
                    # so the volumes agree to rounding.
                    self.assertAlmostEqual(properties.Mass(), volume, delta=1e-9 * volume)
                    for point, expected in probes:
                        actual = BRepClass3d_SolidClassifier(shape, gp_Pnt(*point), 1e-7).State()
                        self.assertEqual(actual, expected, (name, point))


if __name__ == '__main__':
    unittest.main()
