"""Fillet strip ends and corners, read without ShapeFix or sewing.

Pocket corners close with a concave sphere sector, two strips of unequal
reach meet on a step, an arc strip ends at a meridian, a straight strip ends
on a sphere, and ring rounds on a cylinder's bottom edge and a blind pocket's
floor write readable tori. Each answer is one raw-valid solid of the expected
volume.
"""
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace,
                                BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeWire)
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCP.GC import GC_MakeArcOfCircle
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_FACE
from OCP.gp import gp_Pnt, gp_Vec

from test_loopless_faces import BINARY, DRIVER, raw_solid, unique

LEVEL = ['1', '0', '0', '0']
PI = math.pi


def pose(at):
    return {'rotation': LEVEL, 'translation': at}


def side(step, segment):
    return {'step': step, 'role': {'side': {'segment': segment}}}


def cap(step, role):
    return {'step': step, 'role': role}


def rounded(steps, body, radius, edges):
    return {'format': 'catallax-recipe/3', 'name': 'strip ends', 'result': 'round',
            'steps': steps + [{'id': 'round', 'op': 'fillet', 'body': body, 'radius': radius,
                               'edges': [{'between': pair} for pair in edges]}]}


POCKET = [
    {'id': 'block', 'op': 'block', 'size': ['20', '20', '10'], 'pose': pose(['0', '0', '0'])},
    {'id': 'tool', 'op': 'block', 'size': ['10', '10', '10'], 'pose': pose(['5', '5', '4'])},
    {'id': 'pocket', 'op': 'boolean', 'kind': 'subtract', 'target': 'block', 'tool': 'tool'},
]
WALLS = ['x_min', 'y_min', 'x_max', 'y_max']
FLOOR = [[cap('tool', 'start'), side('tool', wall)] for wall in WALLS]
UPRIGHTS = [[side('tool', WALLS[i]), side('tool', WALLS[(i + 1) % 4])] for i in range(4)]
SLOPE = [
    {'id': 'profile', 'op': 'sketch', 'plane': pose(['0', '0', '0']),
     'loops': [[{'kind': 'line', 'id': 'bottom', 'to': ['12', '0']},
                {'kind': 'line', 'id': 'right', 'to': ['12', '5']},
                {'kind': 'line', 'id': 'slope', 'to': ['0', '10']},
                {'kind': 'line', 'id': 'left', 'to': ['0', '0']}]]},
    {'id': 'body', 'op': 'extrude', 'profile': 'profile', 'distance': '10'},
]


def half_disc(at):
    return {'id': 'profile', 'op': 'sketch', 'plane': pose(['0', '0', at]),
            'loops': [[{'kind': 'arc', 'id': 'round', 'centre': ['0', '0'], 'to': ['0', '5'],
                        'turn': 'ccw'},
                       {'kind': 'line', 'id': 'flat', 'to': ['0', '-5']}]]}


D_BOSS = [half_disc('0'), {'id': 'body', 'op': 'extrude', 'profile': 'profile', 'distance': '10'}]
D_POCKET = [
    {'id': 'block', 'op': 'block', 'size': ['20', '20', '10'], 'pose': pose(['-10', '-10', '0'])},
    half_disc('4'),
    {'id': 'tool', 'op': 'extrude', 'profile': 'profile', 'distance': '10'},
    {'id': 'pocket', 'op': 'boolean', 'kind': 'subtract', 'target': 'block', 'tool': 'tool'},
]
BOSS = [{'id': 'boss', 'op': 'cylinder', 'radius': '5', 'height': '10',
         'pose': pose(['0', '0', '0'])}]
BLIND = [
    {'id': 'block', 'op': 'block', 'size': ['20', '20', '10'], 'pose': pose(['-10', '-10', '0'])},
    {'id': 'tool', 'op': 'cylinder', 'radius': '5', 'height': '10', 'pose': pose(['0', '0', '4'])},
    {'id': 'pocket', 'op': 'boolean', 'kind': 'subtract', 'target': 'block', 'tool': 'tool'},
]
PILLOW = [
    {'id': 'block', 'op': 'block', 'size': ['10', '10', '10'], 'pose': pose(['-5', '-5', '0'])},
    {'id': 'ball', 'op': 'sphere', 'radius': '12', 'pose': pose(['0', '0', '-2'])},
    {'id': 'cap', 'op': 'boolean', 'kind': 'intersect', 'target': 'block', 'tool': 'ball'},
]
# Volumes in closed form from the kernel's tests, except the sphere wall's,
# which is a quadrature of the ball's height over the square less the kite.
FIXTURES = [
    ('pocket-corner', rounded(POCKET, 'pocket', '1', FLOOR[:2] + UPRIGHTS[:1]),
     3424 - 71 * PI / 12),
    ('whole-pocket', rounded(POCKET, 'pocket', '1', FLOOR + UPRIGHTS), 3456 - 41 * PI / 3),
    ('half-disc-top', rounded(D_BOSS, 'body', '1', [[cap('body', 'end'), side('body', 'round')]]),
     725 * PI / 6 + PI ** 2),
    ('half-disc-pocket-floor',
     rounded(D_POCKET, 'pocket', '1', [[cap('tool', 'start'), side('tool', 'round')]]),
     4000 - 425 * PI / 6 - PI ** 2),
    ('cylinder-bottom', rounded(BOSS, 'boss', '1', [[cap('boss', 'start'), side('boss', 'wall')]]),
     725 * PI / 3 + 2 * PI ** 2),
    ('blind-pocket-floor',
     rounded(BLIND, 'pocket', '1', [[cap('tool', 'start'), side('tool', 'wall')]]),
     4000 - 425 * PI / 3 - 2 * PI ** 2),
    ('sphere-wall',
     rounded(PILLOW, 'cap', '1', [[side('block', 'x_max'), side('block', 'y_max')]]),
     925.731827569374),
]


def volume(shape):
    properties = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(shape, properties, 1e-8, False, False)
    return properties.Mass()


def kite(edge, direction, normals, length):
    """The blend region of a convex straight edge between planes with these
    outward unit normals, radius 1, as a prism from one unit before the edge
    point along `direction` to `length` past it."""
    a, b = normals
    dot = sum(x * y for x, y in zip(a, b))
    centre = [e - (x + y) / (1 + dot) for e, x, y in zip(edge, a, b)]
    feet = [[c + x for c, x in zip(centre, n)] for n in (a, b)]
    toward = [e - c for e, c in zip(edge, centre)]
    scale = 1 / math.sqrt(sum(x * x for x in toward))
    middle = [c + x * scale for c, x in zip(centre, toward)]
    start = [-x for x in direction]
    moved = [gp_Pnt(*[p + s for p, s in zip(q, start)]) for q in (edge, feet[0], middle, feet[1])]
    wire = BRepBuilderAPI_MakeWire()
    wire.Add(BRepBuilderAPI_MakeEdge(moved[0], moved[1]).Edge())
    wire.Add(BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(moved[1], moved[2], moved[3]).Value()).Edge())
    wire.Add(BRepBuilderAPI_MakeEdge(moved[3], moved[0]).Edge())
    face = BRepBuilderAPI_MakeFace(wire.Wire()).Face()
    return BRepPrimAPI_MakePrism(face, gp_Vec(*[x * (length + 1) for x in direction])).Shape()


def cut(a, b):
    operation = BRepAlgoAPI_Cut(a, b)
    operation.Build()
    return operation.Shape()


def unequal_oracle():
    """Section 4.5 directly: the slope prism less the two kite prisms, each
    run past its end faces, so each is cut by the other's outer face."""
    polygon = BRepBuilderAPI_MakePolygon()
    for x, y in ((0, 0), (12, 0), (12, 5), (0, 10)):
        polygon.Add(gp_Pnt(x, y, 0))
    polygon.Close()
    prism = BRepPrimAPI_MakePrism(BRepBuilderAPI_MakeFace(polygon.Wire()).Face(),
                                  gp_Vec(0, 0, 10)).Shape()
    slope = (5 / 13, 12 / 13, 0)
    upright = kite((12, 5, 0), (0, 0, 1), ((1, 0, 0), slope), 12)
    top = kite((12, 5, 10), (-12 / 13, 5 / 13, 0), ((0, 0, 1), slope), 16)
    return volume(cut(cut(prism, upright), top))


class StripEndImports(unittest.TestCase):
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

    def test_ends_and_corners_write_raw_valid_solids(self):
        with tempfile.TemporaryDirectory(prefix='catallax-strip-ends-') as temporary:
            root = Path(temporary)
            for name, recipe, expected in FIXTURES:
                with self.subTest(name):
                    shape = self.read(root, name, recipe)
                    self.assertAlmostEqual(volume(shape), expected, delta=1e-6)

    def test_unequal_reach_matches_the_definition_built_by_booleans(self):
        recipe = rounded(SLOPE, 'body', '1', [[side('body', 'slope'), side('body', 'right')],
                                               [side('body', 'slope'), cap('body', 'end')]])
        with tempfile.TemporaryDirectory(prefix='catallax-strip-ends-') as temporary:
            shape = self.read(Path(temporary), 'unequal', recipe)
            self.assertAlmostEqual(volume(shape), unequal_oracle(), delta=1e-6)


if __name__ == '__main__':
    unittest.main()
