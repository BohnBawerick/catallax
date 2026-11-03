"""Torus rounds between coaxial surfaces of revolution, read without ShapeFix or sewing."""
import json
from math import atan2, cos, hypot, pi, sin, sqrt, tau
from pathlib import Path
import subprocess
import tempfile
import unittest

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID

from test_loopless_faces import BINARY, DRIVER, raw_solid, unique

UP = ['1', '0', '0', '0']


def pose(z):
    return {'rotation': UP, 'translation': ['0', '0', z]}


def face(step, role):
    return {'step': step, 'role': role}


def side(step, segment):
    return {'step': step, 'role': {'side': {'segment': segment}}}


def block(name, z, height):
    return [
        {'id': f'{name}-section', 'op': 'sketch', 'plane': {'rotation': UP, 'translation': ['0', '0', z]},
         'loops': [[{'kind': 'line', 'id': 's', 'to': ['20', '-20']},
                    {'kind': 'line', 'id': 'e', 'to': ['20', '20']},
                    {'kind': 'line', 'id': 'n', 'to': ['-20', '20']},
                    {'kind': 'line', 'id': 'w', 'to': ['-20', '-20']}]]},
        {'id': name, 'op': 'extrude', 'profile': f'{name}-section', 'distance': height},
    ]


def cone(base, top, height):
    return [{'id': 'c', 'op': 'cone', 'base_radius': base, 'top_radius': top, 'height': height,
             'pose': pose('0')}]


CHAMFERED = [
    {'id': 'c', 'op': 'cylinder', 'radius': '5', 'height': '10', 'pose': pose('0')},
    {'id': 'bevel', 'op': 'chamfer', 'body': 'c', 'distance': '1',
     'edges': [{'between': [face('c', 'end'), side('c', 'wall')]}]},
]
DOME = [{'id': 's', 'op': 'sphere', 'radius': '5', 'pose': pose('0')}] + block('b', '2', '10') + [
    {'id': 'j', 'op': 'boolean', 'kind': 'intersect', 'target': 's', 'tool': 'b'}]
SEAT = block('b', '0', '20') + [
    {'id': 's', 'op': 'sphere', 'radius': '5', 'pose': pose('23')},
    {'id': 'j', 'op': 'boolean', 'kind': 'union', 'target': 'b', 'tool': 's'}]
COUNTERSUNK = block('b', '0', '20') + [
    {'id': 'h', 'op': 'hole', 'target': 'b',
     'at': [{'id': 'p', 'point': ['0', '0', '20'], 'direction': ['0', '0', '-1']}],
     'diameter': '6', 'depth': 'through', 'tip': 'flat',
     'head': {'countersink': {'diameter': '12', 'angle': '90'}}}]


def hole(part):
    return {'step': 'h', 'role': {'hole': {'at': 'p', 'part': part}}}


def revolved(segments):
    return [
        {'id': 'section', 'op': 'sketch', 'plane': {'rotation': UP, 'translation': ['0', '0', '0']},
         'loops': [segments]},
        {'id': 'r', 'op': 'revolve', 'profile': 'section',
         'axis': {'point': ['0', '0', '0'], 'direction': ['0', '1', '0']}, 'angle': '360'},
    ]


GROOVE = revolved([
    {'kind': 'line', 'id': 'floor', 'to': ['6', '0']},
    {'kind': 'arc', 'id': 'groove', 'centre': ['10', '0'], 'to': ['10', '4'], 'turn': 'cw'},
    {'kind': 'line', 'id': 'rim', 'to': ['10', '6']},
    {'kind': 'line', 'id': 'top', 'to': ['0', '6']},
    {'kind': 'line', 'id': 'axis', 'to': ['0', '0']}])
ROD = revolved([
    {'kind': 'line', 'id': 'base', 'to': ['6', '0']},
    {'kind': 'line', 'id': 'rod', 'to': ['6', '8']},
    {'kind': 'arc', 'id': 'cap', 'centre': ['0', '0'], 'to': ['0', '10'], 'turn': 'ccw'},
    {'kind': 'line', 'id': 'axis', 'to': ['0', '0']}])


# The half-plane through the axis holds each profile as a line n.p = d, with
# unit n, or a circle. The region between an edge point and its ball turns
# about the axis; its volume is 2 pi times the boundary integral of rho^2/2 dz.
def moment(profile, p, q):
    if profile[0] == 'line':
        dr, dz = q[0] - p[0], q[1] - p[1]
        return dz / 2 * (p[0] ** 2 + p[0] * dr + dr ** 2 / 3)
    _, c, radius = profile
    start = atan2(p[1] - c[1], p[0] - c[0])
    turn = (atan2(q[1] - c[1], q[0] - c[0]) - start + pi) % tau - pi

    def f(t):
        return radius / 2 * (c[0] ** 2 * sin(t) + 2 * c[0] * radius * (t / 2 + sin(2 * t) / 4)
                             + radius ** 2 * (sin(t) - sin(t) ** 3 / 3))
    return f(start + turn) - f(start)


def foot(profile, c):
    if profile[0] == 'line':
        _, n, d = profile
        off = n[0] * c[0] + n[1] * c[1] - d
        return (c[0] - off * n[0], c[1] - off * n[1])
    _, m, radius = profile
    length = hypot(c[0] - m[0], c[1] - m[1])
    return (m[0] + (c[0] - m[0]) * radius / length, m[1] + (c[1] - m[1]) * radius / length)


def region(edge, profiles, centre, radius):
    a, b = foot(profiles[0], centre), foot(profiles[1], centre)
    total = (moment(profiles[0], edge, a) + moment(('circle', centre, radius), a, b)
             + moment(profiles[1], b, edge))
    return tau * abs(total)


def meet(a, d, b, e):
    det = a[0] * b[1] - a[1] * b[0]
    return ((d * b[1] - e * a[1]) / det, (a[0] * e - b[0] * d) / det)


def unit(x, y):
    length = hypot(x, y)
    return (x / length, y / length)


def line(n, d):
    return ('line', n, d)


FLAT = (0.0, 1.0)
WALL = (1.0, 0.0)
FRUSTUM = unit(3, 1), 15 / sqrt(10)
FLARE = unit(3, -1), 9 / sqrt(10)
BEVEL = unit(1, 1), 14 / sqrt(2)
SINK = unit(1, -1), -14 / sqrt(2)
TUBE = ('circle', (10.0, 0.0), 4.0)


def frustum(n, d, plane, ball, edge):
    # The plane z = plane, and the ball of radius 1 centred at height ball.
    return (edge, [line(FLAT, plane), line(n, d)], meet(n, d - 1.0, FLAT, ball), 1.0)


def bevel_top(r):
    n, d = BEVEL
    return (4.0, 10.0), [line(FLAT, 10.0), line(n, d)], meet(n, d - r, FLAT, 10.0 - r), r


def bevel_side(r):
    n, d = BEVEL
    return (5.0, 9.0), [line(WALL, 5.0), line(n, d)], meet(n, d - r, WALL, 5.0 - r), r


def mouth():
    n, d = SINK
    return (6.0, 20.0), [line(FLAT, 20.0), line(n, d)], meet(n, d + 1, FLAT, 19.0), 1.0


def throat():
    n, d = SINK
    return (3.0, 17.0), [line(n, d), line(WALL, 3.0)], meet(n, d + 1, WALL, 4.0), 1.0


SINK_VOLUME = 32000 - 216 * pi
BEVEL_VOLUME = 250 * pi - 14 * pi / 3
# name, steps, body, radius, edges, balls, volume before, sign of the change.
CASES = [
    ('frustum-top', cone('5', '3', '6'), 'c', '1', [[face('c', 'end'), side('c', 'wall')]],
     [frustum(*FRUSTUM, 6.0, 5.0, (3.0, 6.0))], 98 * pi, -1),
    ('frustum-base', cone('5', '3', '6'), 'c', '1', [[face('c', 'start'), side('c', 'wall')]],
     [frustum(*FRUSTUM, 0.0, 1.0, (5.0, 0.0))], 98 * pi, -1),
    ('flared-top', cone('3', '5', '6'), 'c', '1', [[face('c', 'end'), side('c', 'wall')]],
     [frustum(*FLARE, 6.0, 5.0, (5.0, 6.0))], 98 * pi, -1),
    ('chamfer-top', CHAMFERED, 'bevel', '1/2', [[face('c', 'end'), {'step': 'bevel'}]],
     [bevel_top(0.5)], BEVEL_VOLUME, -1),
    ('chamfer-side', CHAMFERED, 'bevel', '1/2', [[side('c', 'wall'), {'step': 'bevel'}]],
     [bevel_side(0.5)], BEVEL_VOLUME, -1),
    ('chamfer-both', CHAMFERED, 'bevel', '1/4',
     [[face('c', 'end'), {'step': 'bevel'}], [side('c', 'wall'), {'step': 'bevel'}]],
     [bevel_top(0.25), bevel_side(0.25)], BEVEL_VOLUME, -1),
    ('dome-rim', DOME, 'j', '1', [[face('b', 'start'), side('s', 'wall')]],
     [((sqrt(21), 2.0), [line(FLAT, 2.0), ('circle', (0.0, 0.0), 5.0)], (sqrt(7), 3.0), 1.0)],
     36 * pi, -1),
    ('sphere-seat', SEAT, 'j', '1', [[face('b', 'end'), side('s', 'wall')]],
     [((4.0, 20.0), [line(FLAT, 20.0), ('circle', (0.0, 23.0), 5.0)], (sqrt(32), 21.0), 1.0)],
     32000 + 448 * pi / 3, 1),
    ('countersink-both', COUNTERSUNK, 'h', '1',
     [[face('b', 'end'), hole('countersink')], [hole('countersink'), hole('wall')]],
     [mouth(), throat()], SINK_VOLUME, -1),
    ('groove', GROOVE, 'r', '1',
     [[side('r', 'floor'), side('r', 'groove')], [side('r', 'groove'), side('r', 'rim')]],
     [((6.0, 0.0), [line(FLAT, 0.0), TUBE], (10 - sqrt(24), 1.0), 1.0),
      ((10.0, 4.0), [TUBE, line(WALL, 10.0)], (9.0, sqrt(24)), 1.0)],
     600 * pi - 80 * pi ** 2 + 128 * pi / 3, -1),
    ('rod-cap', ROD, 'r', '1', [[side('r', 'rod'), side('r', 'cap')]],
     [((6.0, 8.0), [line(WALL, 6.0), ('circle', (0.0, 0.0), 10.0)], (5.0, sqrt(56)), 1.0)],
     288 * pi + 112 * pi / 3, -1),
]


# Revolved bodies whose torus bands lie beside a cone at tube angles other than
# quarter turns, or run across the outer equator. Their volume is the boundary
# moment of the profile, turned.
def profile_volume(pieces):
    return tau * sum(moment(*piece) for piece in pieces)


RIM = revolved([
    {'kind': 'line', 'id': 'bottom', 'to': ['8', '0']},
    {'kind': 'arc', 'id': 'rim', 'centre': ['8', '3'], 'to': ['8', '6'], 'turn': 'ccw'},
    {'kind': 'line', 'id': 'top', 'to': ['0', '6']},
    {'kind': 'line', 'id': 'axis', 'to': ['0', '0']}])
CONE_RIM = revolved([
    {'kind': 'line', 'id': 'bottom', 'to': ['8', '0']},
    {'kind': 'line', 'id': 'slope', 'to': ['7', '3']},
    {'kind': 'arc', 'id': 'rim', 'centre': ['4', '3'], 'to': ['4', '6'], 'turn': 'ccw'},
    {'kind': 'line', 'id': 'top', 'to': ['0', '6']},
    {'kind': 'line', 'id': 'axis', 'to': ['0', '0']}])
PLAIN = [
    ('rim', RIM, profile_volume([
        (line(FLAT, 0.0), (0.0, 0.0), (8.0, 0.0)),
        (('circle', (8.0, 3.0), 3.0), (8.0, 0.0), (11.0, 3.0)),
        (('circle', (8.0, 3.0), 3.0), (11.0, 3.0), (8.0, 6.0)),
        (line(FLAT, 6.0), (8.0, 6.0), (0.0, 6.0)),
        (line(WALL, 0.0), (0.0, 6.0), (0.0, 0.0))])),
    ('cone-rim', CONE_RIM, profile_volume([
        (line(FLAT, 0.0), (0.0, 0.0), (8.0, 0.0)),
        (line(WALL, 0.0), (8.0, 0.0), (7.0, 3.0)),
        (('circle', (4.0, 3.0), 3.0), (7.0, 3.0), (4.0, 6.0)),
        (line(FLAT, 6.0), (4.0, 6.0), (0.0, 6.0)),
        (line(WALL, 0.0), (0.0, 6.0), (0.0, 0.0))])),
]


class CoaxialFilletImports(unittest.TestCase):
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
        manifest = case['response']['manifest']
        self.assertTrue(manifest['exact'], manifest)
        self.assertEqual(manifest['marks'], [])
        shape = raw_solid(Path(case['path']) / 'answer.step')
        self.assertTrue(BRepCheck_Analyzer(shape).IsValid())
        self.assertEqual(len(unique(shape, TopAbs_SOLID)), 1)
        for each in unique(shape, TopAbs_FACE):
            self.assertTrue(BRepCheck_Analyzer(each).IsValid())
        properties = GProp_GProps()
        BRepGProp.VolumeProperties_s(shape, properties, 1e-8)
        return properties.Mass()

    def test_revolved_bands_beside_cones_and_across_the_equator(self):
        with tempfile.TemporaryDirectory(prefix='catallax-coaxial-band-') as temporary:
            for name, steps, want in PLAIN:
                with self.subTest(name):
                    recipe = {'format': 'catallax-recipe/3', 'name': name, 'steps': steps,
                              'result': 'r'}
                    found = self.read(Path(temporary), name, recipe)
                    self.assertAlmostEqual(found, want, delta=want * 1e-7)

    def test_coaxial_rounds(self):
        with tempfile.TemporaryDirectory(prefix='catallax-coaxial-fillet-') as temporary:
            root = Path(temporary)
            for name, steps, body, radius, edges, balls, before, sign in CASES:
                with self.subTest(name):
                    recipe = {'format': 'catallax-recipe/3', 'name': name,
                              'steps': steps + [{'id': 'round', 'op': 'fillet', 'body': body,
                                                 'radius': radius,
                                                 'edges': [{'between': e} for e in edges]}],
                              'result': 'round'}
                    found = self.read(root, name, recipe)
                    want = before + sign * sum(region(*ball) for ball in balls)
                    self.assertAlmostEqual(found, want, delta=want * 1e-7)


if __name__ == '__main__':
    unittest.main()
