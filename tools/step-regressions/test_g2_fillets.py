"""Curvature-continuous fillets: STEP read raw by OpenCASCADE, healing off.

The Rust test `g2_fillets` writes the files. Each is a 20 x 20 x 10 block
with a G2 fillet on top edges along y. The removed volume is 20 mm times the
area between the quintic section and the edge, integrated here in the power
basis, a different method from the Rust test's quadrature.

The G2 check reads OpenCASCADE's own curvature of the blend surface at its
rails. At a rail on a plane both principal curvatures are zero, as the
plane's are, and the surface normal is the plane's. OpenCASCADE's circular
fillet of the same radius, built here for contrast, has curvature 1/r there:
it is tangent to the plane (G1) but its curvature jumps.
"""
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from OCP.BRep import BRep_Tool
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Cylinder
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepLProp import BRepLProp_SLProps
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_IN, TopAbs_OUT, TopAbs_SOLID
from OCP.TopExp import TopExp, TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools/checks'))
from raw_step import raw_shapes


def items(shape, kind):
    walk = TopExp_Explorer(shape, kind)
    while walk.More():
        yield walk.Current()
        walk.Next()


def step_fraction(fullness):
    """The quintic's control step over the contact distance for a right
    angle: 6 w sin(45) / (5 (1 + 2 sin(45))), SMLib's fairing rule."""
    sine = math.sqrt(0.5)
    return 6 * fullness * sine / (5 * (1 + 2 * sine))


def poly_mul(a, b):
    out = [0.0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def power_basis(controls):
    """Power-basis coefficients of a Bezier curve coordinate."""
    n = len(controls) - 1
    coefficients = [0.0] * (n + 1)
    for j, c in enumerate(controls):
        # B_j^n(t) = C(n, j) t^j (1 - t)^(n - j)
        for k in range(n - j + 1):
            coefficients[j + k] += math.comb(n, j) * math.comb(n - j, k) * (-1) ** k * c
    return coefficients


def removed_area(s, fullness):
    f = step_fraction(fullness)
    xs = [-s, -s * (1 - f), -s * (1 - 2 * f), 0.0, 0.0, 0.0]
    ys = [0.0, 0.0, 0.0, -s * (1 - 2 * f), -s * (1 - f), -s]
    x, y = power_basis(xs), power_basis(ys)
    dx = [k * c for k, c in enumerate(x)][1:]
    dy = [k * c for k, c in enumerate(y)][1:]
    integrand = [a - b for a, b in zip(poly_mul(x, dy), poly_mul(y, dx))]
    twice = sum(c / (k + 1) for k, c in enumerate(integrand))
    return abs(twice) / 2


CASES = {
    'g2-box.step': (7, 4000 - 20 * removed_area(2, 1)),
    'g2-box-radius-3.step': (7, 4000 - 20 * removed_area(3, 1)),
    'g2-box-thin.step': (7, 4000 - 20 * removed_area(3, 0.5)),
    'g2-two-edges.step': (8, 4000 - 40 * removed_area(2, 1.25)),
}


def rail_curvatures(face):
    """Largest principal curvature magnitude and the normal at the middle of
    each rail, and the curvature at the middle of the section. The section
    runs along the degree-5 parameter, the rails along the degree-1 one."""
    adaptor = BRepAdaptor_Surface(face)
    u0, u1 = adaptor.FirstUParameter(), adaptor.LastUParameter()
    v0, v1 = adaptor.FirstVParameter(), adaptor.LastVParameter()
    um, vm = (u0 + u1) / 2, (v0 + v1) / 2
    if adaptor.UDegree() == 5:
        points = {'rail0': (u0, vm), 'rail1': (u1, vm), 'middle': (um, vm)}
    else:
        assert adaptor.VDegree() == 5
        points = {'rail0': (um, v0), 'rail1': (um, v1), 'middle': (um, vm)}
    found = {}
    for name, (u, v) in points.items():
        props = BRepLProp_SLProps(adaptor, u, v, 2, 1e-9)
        assert props.IsCurvatureDefined()
        normal = props.Normal()
        found[name] = (max(abs(props.MaxCurvature()), abs(props.MinCurvature())),
                       (normal.X(), normal.Y(), normal.Z()))
    return found


class G2FilletOracles(unittest.TestCase):
    def test_g2_fillets_read_raw_as_valid_with_volume_and_zero_rail_curvature(self):
        with tempfile.TemporaryDirectory(prefix='catallax-g2-') as directory:
            proof = Path(directory)
            environment = dict(os.environ)
            environment['CATALLAX_G2_PROOF'] = str(proof)
            subprocess.run([
                'cargo', 'test', '--release', '-p', 'catallax-features', '--test', 'g2_fillets',
                '--', '--nocapture',
            ], cwd=ROOT, env=environment, check=True, timeout=1200)
            for name, (faces, volume) in CASES.items():
                with self.subTest(file=name):
                    shapes = raw_shapes(proof / name)
                    self.assertEqual(len(shapes), 1)
                    shape = shapes[0]
                    self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
                    self.assertEqual(len(list(items(shape, TopAbs_SOLID))), 1)
                    self.assertEqual(len(list(items(shape, TopAbs_FACE))), faces)
                    properties = GProp_GProps()
                    error = BRepGProp.VolumePropertiesGK_s(
                        shape, properties, Eps=1e-10, OnlyClosed=True,
                    )
                    self.assertTrue(math.isfinite(error))
                    print(f'{name}: volume {properties.Mass()!r} mm^3, expected {volume!r}')
                    self.assertAlmostEqual(properties.Mass(), volume, delta=1e-6)
                    # A corner of the block's x_max top edge, nearer the edge
                    # than any section passes, and the middle of the block.
                    corner = BRepClass3d_SolidClassifier(shape, gp_Pnt(19.95, 10, 9.95), 1e-7)
                    self.assertEqual(corner.State(), TopAbs_OUT)
                    middle = BRepClass3d_SolidClassifier(shape, gp_Pnt(10, 10, 5), 1e-7)
                    self.assertEqual(middle.State(), TopAbs_IN)
                    blends = [TopoDS.Face(f) for f in items(shape, TopAbs_FACE)
                              if BRepAdaptor_Surface(TopoDS.Face(f)).GetType()
                              == GeomAbs_BSplineSurface]
                    self.assertEqual(len(blends), faces - 6)
                    for blend in blends:
                        found = rail_curvatures(blend)
                        rails = [found[k][0] for k in ('rail0', 'rail1')]
                        print(f'{name}: rail curvature {rails}, middle {found["middle"][0]}')
                        for k in rails:
                            self.assertLess(k, 1e-6)
                        # The section turns between its flat ends. A thin
                        # section turns hardest near its ends, not mid-way.
                        self.assertGreater(found['middle'][0], 0.1)
                        # G1: at each rail the normal is the adjoining plane's,
                        # +z on the top face or +-x on a side face.
                        for key in ('rail0', 'rail1'):
                            n = found[key][1]
                            axis = max(range(3), key=lambda i: abs(n[i]))
                            self.assertIn(axis, (0, 2))
                            self.assertAlmostEqual(abs(n[axis]), 1.0, delta=1e-9)

    def test_opencascade_circular_fillet_jumps_in_curvature_at_its_rail(self):
        box = BRepPrimAPI_MakeBox(20, 20, 10).Shape()
        fillet = BRepFilletAPI_MakeFillet(box)
        for edge in items(box, TopAbs_EDGE):
            edge = TopoDS.Edge(edge)
            ends = [BRep_Tool.Pnt_s(TopExp.FirstVertex_s(edge)),
                    BRep_Tool.Pnt_s(TopExp.LastVertex_s(edge))]
            if all(abs(p.X() - 20) < 1e-9 and abs(p.Z() - 10) < 1e-9 for p in ends):
                fillet.Add(2.0, edge)
        fillet.Build()
        self.assertTrue(fillet.IsDone())
        rounds = [TopoDS.Face(f) for f in items(fillet.Shape(), TopAbs_FACE)
                  if BRepAdaptor_Surface(TopoDS.Face(f)).GetType() == GeomAbs_Cylinder]
        self.assertEqual(len(rounds), 1)
        adaptor = BRepAdaptor_Surface(rounds[0])
        props = BRepLProp_SLProps(adaptor, adaptor.FirstUParameter(),
                                  (adaptor.FirstVParameter() + adaptor.LastVParameter()) / 2,
                                  2, 1e-9)
        self.assertTrue(props.IsCurvatureDefined())
        curvature = max(abs(props.MaxCurvature()), abs(props.MinCurvature()))
        print(f'OpenCASCADE circular fillet: rail curvature {curvature}')
        self.assertAlmostEqual(curvature, 0.5, delta=1e-9)


if __name__ == '__main__':
    unittest.main()
