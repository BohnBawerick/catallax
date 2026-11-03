"""Fillets rolling over earlier fillets, read without ShapeFix or sewing.

A later fillet's tangent chain continues over earlier cylinder rounds as
torus or sphere strips, every remaining edge rounds after an earlier round,
and uprights rounded after a top loop swallow its mitres. Each answer is one
raw-valid solid. Its volume matches OpenCASCADE's own fillets of the same
part where those build, and the closed form otherwise: OpenCASCADE does not
round the uprights under a mitred top loop.
"""
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.GProp import GProp_GProps
from OCP.GeomAbs import GeomAbs_Line
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

from test_loopless_faces import BINARY, DRIVER, raw_solid, unique

LEVEL = ['1', '0', '0', '0']
PI = math.pi
SIDES = ['x_min', 'y_min', 'x_max', 'y_max']


def pose(at):
    return {'rotation': LEVEL, 'translation': at}


def side(step, segment):
    return {'step': step, 'role': {'side': {'segment': segment}}}


def cap(step, role):
    return {'step': step, 'role': role}


def block(name, size, at):
    return {'id': name, 'op': 'block', 'size': size, 'pose': pose(at)}


def fillet(name, body, radius, edges):
    if edges != 'all':
        edges = [{'between': pair} for pair in edges]
    return {'id': name, 'op': 'fillet', 'body': body, 'radius': radius, 'edges': edges}


def uprights(step):
    return [[side(step, SIDES[i]), side(step, SIDES[(i + 1) % 4])] for i in range(4)]


def recipe(steps):
    return {'format': 'catallax-recipe/4', 'name': 'fillet over fillet',
            'result': steps[-1]['id'], 'steps': steps}


BLOCK = block('block', ['20', '20', '10'], ['0', '0', '0'])


def boss(upright, top):
    return recipe([BLOCK, fillet('upright', 'block', upright, uprights('block')),
                   fillet('top', 'upright', top, [[cap('block', 'end'), side('block', 'x_max')]])])


def pocket(upright, floor):
    return recipe([BLOCK, block('tool', ['10', '10', '10'], ['5', '5', '4']),
                   {'id': 'pocket', 'op': 'boolean', 'kind': 'subtract', 'target': 'block',
                    'tool': 'tool'},
                   fillet('upright', 'pocket', upright, uprights('tool')),
                   fillet('floor', 'upright', floor,
                          [[cap('tool', 'start'), side('tool', 'x_max')]])])


def corner(front, right):
    return recipe([BLOCK,
                   fillet('front', 'block', front, [[cap('block', 'end'), side('block', 'y_min')]]),
                   fillet('right', 'front', right, [[cap('block', 'end'), side('block', 'x_max')]])])


def bracket(inner, right):
    return recipe([block('base', ['20', '20', '5'], ['0', '0', '0']),
                   block('wall', ['20', '5', '10'], ['0', '0', '0']),
                   {'id': 'l', 'op': 'boolean', 'kind': 'union', 'target': 'base', 'tool': 'wall'},
                   fillet('inner', 'l', inner, [[cap('base', 'end'), side('wall', 'y_max')]]),
                   fillet('right', 'inner', right, [[cap('base', 'end'), side('base', 'x_max')]])])


def based(upright, base):
    return recipe([block('plate', ['30', '30', '5'], ['0', '0', '0']),
                   block('boss', ['10', '10', '10'], ['10', '10', '5']),
                   {'id': 'joined', 'op': 'boolean', 'kind': 'union', 'target': 'plate',
                    'tool': 'boss'},
                   fillet('upright', 'joined', upright, uprights('boss')),
                   fillet('base', 'upright', base, [[cap('plate', 'end'), side('boss', 'x_max')]])])


def every(front, rest):
    return recipe([BLOCK,
                   fillet('front', 'block', front, [[cap('block', 'end'), side('block', 'y_min')]]),
                   fillet('rest', 'front', rest, 'all')])


def mitred(radius):
    return recipe([BLOCK,
                   fillet('top', 'block', radius,
                          [[cap('block', 'end'), side('block', s)] for s in SIDES]),
                   fillet('upright', 'top', radius, uprights('block'))])


def volume(shape):
    properties = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(shape, properties, 1e-8, False, False)
    return properties.Mass()


def box(x, y, z, dx, dy, dz):
    return BRepPrimAPI_MakeBox(gp_Pnt(x, y, z), dx, dy, dz).Shape()


def straight(shape, **fixed):
    """The straight edges whose ends both have these coordinates."""
    axes = {'x': 0, 'y': 1, 'z': 2}
    found = []
    for edge in unique(shape, TopAbs_EDGE):
        curve = BRepAdaptor_Curve(TopoDS.Edge(edge))
        if curve.GetType() != GeomAbs_Line:
            continue
        ends = [curve.Value(curve.FirstParameter()), curve.Value(curve.LastParameter())]
        ends = [(p.X(), p.Y(), p.Z()) for p in ends]
        if all(abs(p[axes[k]] - v) < 1e-9 for k, v in fixed.items() for p in ends):
            found.append(TopoDS.Edge(edge))
    return found


def rounded(shape, radius, edges):
    """OpenCASCADE's own fillet, which continues each edge over its tangent chain."""
    operation = BRepFilletAPI_MakeFillet(shape)
    for edge in edges:
        operation.Add(radius, edge)
    operation.Build()
    if not operation.IsDone():
        raise AssertionError('OpenCASCADE fillet failed')
    return operation.Shape()


def occ_boss(upright, top):
    shape = box(0, 0, 0, 20, 20, 10)
    columns = [e for x in (0, 20) for y in (0, 20) for e in straight(shape, x=x, y=y)]
    shape = rounded(shape, upright, columns)
    return volume(rounded(shape, top, straight(shape, x=20, z=10)))


def occ_pocket(upright, floor):
    cut = BRepAlgoAPI_Cut(box(0, 0, 0, 20, 20, 10), box(5, 5, 4, 10, 10, 10))
    cut.Build()
    shape = cut.Shape()
    columns = [e for x in (5, 15) for y in (5, 15) for e in straight(shape, x=x, y=y)]
    shape = rounded(shape, upright, columns)
    return volume(rounded(shape, floor, straight(shape, x=15, z=4)))


def occ_corner(front, right):
    shape = box(0, 0, 0, 20, 20, 10)
    shape = rounded(shape, front, straight(shape, y=0, z=10))
    return volume(rounded(shape, right, straight(shape, x=20, z=10)))


def occ_bracket(inner, right):
    fuse = BRepAlgoAPI_Fuse(box(0, 0, 0, 20, 20, 5), box(0, 0, 0, 20, 5, 10))
    fuse.Build()
    unify = ShapeUpgrade_UnifySameDomain(fuse.Shape(), True, True, False)
    unify.Build()
    shape = rounded(unify.Shape(), inner, straight(unify.Shape(), y=5, z=5))
    return volume(rounded(shape, right, straight(shape, x=20, z=5)))


def occ_based(upright, base):
    fuse = BRepAlgoAPI_Fuse(box(0, 0, 0, 30, 30, 5), box(10, 10, 5, 10, 10, 10))
    fuse.Build()
    unify = ShapeUpgrade_UnifySameDomain(fuse.Shape(), True, True, False)
    unify.Build()
    shape = unify.Shape()
    columns = [e for x in (10, 20) for y in (10, 20) for e in straight(shape, x=x, y=y)]
    shape = rounded(shape, upright, columns)
    return volume(rounded(shape, base, straight(shape, x=20, z=5)))


def rounded_box(h, x, y, radii, r):
    """A box `h` long with an `x` by `y` section whose four long edges are
    rounded by `radii` and every edge of both ends by `r`: the middle section
    over `h - 2r`, and at each end a section offset inward by
    `r - sqrt(r^2 - s^2)` at depth `s`."""
    area = x * y - (1 - PI / 4) * sum(radius * radius for radius in radii)
    x1, y1 = x - 2 * r, y - 2 * r
    end = (x1 * y1 * r + (x1 + y1) * PI * r * r / 2 + 8 * r ** 3 / 3
           - (1 - PI / 4) * sum((radius - r) ** 2 * r + (radius - r) * PI * r * r / 2
                                + 2 * r ** 3 / 3 for radius in radii))
    return (h - 2 * r) * area + 2 * end


class FilletOverFilletImports(unittest.TestCase):
    def read(self, root, name, recipe):
        source = root / f'{name}.json'
        source.write_text(json.dumps(recipe))
        out = root / name
        completed = subprocess.run(
            [str(DRIVER), 'run', str(source), '--source', 'recipe', '--binary',
             str(BINARY), '--out', str(out), '--fresh', '--workers', '1', '--json'],
            capture_output=True, text=True, check=True)
        self.assertEqual(json.loads(completed.stdout)['outcomes'], {'answered': 1})
        case = json.loads((out / 'run.json').read_text())['cases'][0]
        self.assertEqual(case['response']['manifest']['marks'], [])
        shape = raw_solid(Path(case['path']) / 'answer.step')
        self.assertTrue(BRepCheck_Analyzer(shape).IsValid(), case['path'])
        for face in unique(shape, TopAbs_FACE):
            self.assertTrue(BRepCheck_Analyzer(face).IsValid(), case['path'])
        return shape

    def test_rolled_chains_match_opencascade_fillets(self):
        fixtures = [
            ('boss-torus', boss('3', '1'), occ_boss(3, 1)),
            ('boss-torus-sheet', boss('3', '2'), occ_boss(3, 2)),
            ('boss-sphere', boss('3', '3'), occ_boss(3, 3)),
            ('pocket-torus', pocket('3', '1'), occ_pocket(3, 1)),
            ('pocket-sphere', pocket('3', '3'), occ_pocket(3, 3)),
            ('corner-torus', corner('3', '1'), occ_corner(3, 1)),
            ('corner-sphere', corner('2', '2'), occ_corner(2, 2)),
            ('bracket', bracket('3', '1'), occ_bracket(3, 1)),
            ('boss-base', based('2', '1'), occ_based(2, 1)),
        ]
        with tempfile.TemporaryDirectory(prefix='catallax-fillet-over-fillet-') as temporary:
            root = Path(temporary)
            for name, source, expected in fixtures:
                with self.subTest(name):
                    shape = self.read(root, name, source)
                    self.assertAlmostEqual(volume(shape), expected, delta=1e-6)

    def test_every_edge_and_mitred_corners_match_their_closed_forms(self):
        fixtures = [
            ('every-torus', every('3', '1'), rounded_box(20, 20, 10, [3, 1, 1, 1], 1)),
            ('every-sphere', every('2', '2'), rounded_box(20, 20, 10, [2, 2, 2, 2], 2)),
            # The top loop and four uprights, each run to its corner, and four
            # corner cubes less their ball octants.
            ('mitred', mitred('3'),
             4000 - 9 * (1 - PI / 4) * (4 * 14 + 4 * 7) - 4 * (27 - 27 * PI / 6)),
        ]
        with tempfile.TemporaryDirectory(prefix='catallax-fillet-over-fillet-') as temporary:
            root = Path(temporary)
            for name, source, expected in fixtures:
                with self.subTest(name):
                    shape = self.read(root, name, source)
                    self.assertAlmostEqual(volume(shape), expected, delta=1e-6)


if __name__ == '__main__':
    unittest.main()
