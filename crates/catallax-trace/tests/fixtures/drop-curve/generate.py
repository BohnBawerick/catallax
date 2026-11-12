#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace curve drop gate.

usage: generate.py ORACLE [RAW_DIR]

ORACLE is the driver built from oracle.cpp against the pinned SMLib. Each
case runs in a fresh process under a memory and CPU cap and a timeout. RAW_DIR,
when given, receives every request, response and error stream; keep it outside
the repository.
"""
import importlib.util
import math
import pathlib
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent


def load(name, folder):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / folder / 'generate.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


surfaces = load('surfaces', 'surface-solve')
curves = load('curves', 'curve-global')
h, surface, curve = surfaces.h, surfaces.surface, curves.curve
KINDS = ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
CASES = []


def v(*xs):
    return ' '.join(h(x) for x in xs)


def basis(knots, p, u):
    """Cox-de Boor basis values at u (u inside the last span uses the span)."""
    n = len(knots) - p - 1
    span = max(i for i in range(p, n) if knots[i] <= u)
    values = [0.0] * n
    for i in range(n):
        values[i] = 1.0 if knots[i] <= u < knots[i + 1] or (i == span and knots[i] <= u <= knots[i + 1]) else 0.0
    for d in range(1, p + 1):
        for i in range(n):
            a = (u - knots[i]) / (knots[i + d] - knots[i]) * values[i] if knots[i + d] > knots[i] else 0.0
            b = ((knots[i + d + 1] - u) / (knots[i + d + 1] - knots[i + 1]) * values[i + 1]
                 if i + 1 < n and knots[i + d + 1] > knots[i + 1] else 0.0)
            values[i] = a + b
    return values


class Patch:
    """A non-rational B-spline surface with its control net, v fastest."""

    def __init__(self, du, dv, ku, kv, net):
        self.du, self.dv, self.ku, self.kv, self.net = du, dv, ku, kv, net
        self.text = 's ' + surface(du, dv, ku, kv, [p for row in net for p in row])

    def iso(self, axis, value):
        """The parameter line at u = value (axis 0) or v = value (axis 1) as a curve."""
        if axis == 0:
            b = basis(self.ku, self.du, value)
            points = [tuple(sum(b[i] * self.net[i][j][c] for i in range(len(self.net))) for c in range(3))
                      for j in range(len(self.net[0]))]
            return curve(self.dv, self.kv, points)
        b = basis(self.kv, self.dv, value)
        points = [tuple(sum(b[j] * self.net[i][j][c] for j in range(len(self.net[0]))) for c in range(3))
                  for i in range(len(self.net))]
        return curve(self.du, self.ku, points)


K = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3]
BUMP = Patch(3, 3, K, K, [[(i, j, 1.0 * math.sin(0.9 * i + 0.3) * math.cos(0.7 * j - 0.3) + 0.05 * i)
                           for j in range(6)] for i in range(6)])
HILL = Patch(3, 3, K, K, [[(i, j, 1.0 if 2 <= i <= 3 and 2 <= j <= 3 else 0.0) for j in range(6)] for i in range(6)])
# A bump with a C0 knot line at u = 2.
KC = [0, 0, 0, 0, 1, 2, 2, 2, 3, 3, 3, 3]
CREASE = Patch(3, 3, KC, K, [[(i * 3 / 7, j, (0.3 if i == 4 else 0.0) + 0.1 * math.sin(i + j)) for j in range(6)]
                             for i in range(8)])
BOX, PBOX = (0, 0, 3, 3), (-20, -20, 20, 20)


def plane(origin, x, y):
    return f'p 1 {v(*origin, *x, *y, 1, 1, 0, 0, 1, 1)}'


LEVEL = plane((0, 0, 0), (1, 0, 0), (0, 1, 0))
TILT = plane((0, 0, 0.35), (1, 0, 0.3), (0, 1, -0.2))


def drop(kind, s, domain, c, interval, tolerance=1e-3, keep_all=1, on_surface=1, op='drop'):
    assert kind in KINDS
    CASES.append((op, kind, f'{s} {v(*domain)} {c} {v(*interval)} {v(tolerance)} {keep_all} {on_surface}'))


# Parameter lines of the bump, whole, over part of the curve, over a part
# of the domain (the line leaves it), along v and u, at a knot, at sides,
# reversed, and at two tolerances; with and without keeping far drops.
for axis, value, kind in ((0, 1.3, 'regular'), (1, 0.7, 'regular'), (0, 2.0, 'regular'), (0, 0.0, 'boundary'),
                          (1, 3.0, 'boundary'), (0, 3.0, 'boundary')):
    c = BUMP.iso(axis, value)
    drop(kind, BUMP.text, BOX, c, (0, 3))
    drop(kind, BUMP.text, BOX, c, (0.4, 2.2))
    # A side line outside the smaller domain: SMLib finds no start and
    # tries its general curve drop, which is not ported.
    drop(kind, BUMP.text, (0.5, 0.5, 2.5, 2.5), c, (0, 3), op='drefuse' if value in (0.0, 3.0) else 'drop')
    drop(kind, BUMP.text, BOX, c, (0, 3), 1e-6)
    drop(kind, BUMP.text, BOX, c, (0, 3), 1e-3, 0)
    drop(kind, BUMP.text, BOX, c, (0, 3), 0.0)
# The hill's sides and middle lines on the hill and on the level plane
# through its base (the intersector's range answers).
for axis, value in ((0, 0.0), (1, 0.0), (0, 3.0), (0, 1.5)):
    c = HILL.iso(axis, value)
    drop('boundary', HILL.text, BOX, c, (0, 3))
    if value in (0.0, 3.0):
        drop('boundary', LEVEL, PBOX, c, (0, 3))
        drop('boundary', LEVEL, PBOX, c, (0.5, 2.5))
        drop('boundary', LEVEL, PBOX, c, (0, 3), 1e-6, 0)
# A creased surface: parameter lines across and along its C0 line.
for axis, value in ((1, 1.2), (0, 2.0), (0, 1.0)):
    c = CREASE.iso(axis, value)
    drop('touching', CREASE.text, BOX, c, (0, 3))
# Curves in planes: straight and curved, inside the domain, 1e-4 off with
# and without keeping far drops, on a tilted plane.
STRAIGHT = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(-1, -1, 0), (0, 0, 0), (1, 1, 0), (2, 2, 0)])
WAVY = curve(3, [0, 0, 0, 0, 1, 2, 2, 2, 2], [(-2, 0, 0), (-1, 2, 0), (0, -1, 0), (1, 1, 0), (2, 0, 0)])
for c in (STRAIGHT, WAVY):
    drop('regular', LEVEL, PBOX, c, (0, 1) if c is STRAIGHT else (0, 2))
    drop('regular', LEVEL, PBOX, c, (0.2, 0.7))
LIFTED = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(-1, -1, 1e-4), (0, 0, 1e-4), (1, 1, 1e-4), (2, 2, 1e-4)])
drop('near-singular', LEVEL, PBOX, LIFTED, (0, 1))
drop('near-singular', LEVEL, PBOX, LIFTED, (0, 1), 1e-5, 0)
drop('near-singular', LEVEL, PBOX, LIFTED, (0, 1), 1e-5, 1)
TILTED = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0, 0, 0.35), (1, 0.5, 0.55), (2, -0.5, 1.05), (3, 1, 1.05)])
drop('regular', TILT, PBOX, TILTED, (0, 1))
drop('regular', TILT, PBOX, TILTED, (0.3, 0.9))
# Refused: a curve 0.1 off the tilted plane at one control point, which
# SMLib's plane drop hands to the general curve drop.
OFF = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0, 0, 0.35), (1, 0.5, 0.55), (2, -0.5, 1.05), (3, 1, 0.95)])
drop('regular', TILT, PBOX, OFF, (0, 1), op='drefuse')
# Refused: a curve on the bump that is no parameter line (SMLib's general
# curve drop), and a curve not known to lie on the surface.
diag = curve(1, [0, 0, 1, 1], [(0.5, 0.5, 0), (2.5, 2.5, 0)])
drop('regular', BUMP.text, BOX, BUMP.iso(0, 1.3), (0, 3), on_surface=0, op='drefuse')


# Singular: the collapsed patch's parameter lines from its pole and its
# pole as a point curve. SMLib
# moves a start on the pole along the singular direction; the port refuses.
KP = [0, 0, 0, 0, 1, 1, 1, 1]
POLE = Patch(3, 3, KP, KP, [[(0, 0, 1) if i == 0 else (i, j - 1.5, 0.3 * i * j) for j in range(4)] for i in range(4)])
drop('singular', POLE.text, (0, 0, 1, 1), POLE.iso(1, 0.5), (0, 1), op='drefuse')
drop('singular', POLE.text, (0, 0, 1, 1), POLE.iso(0, 0.5), (0, 1))
drop('singular', POLE.text, (0, 0, 1, 1), POLE.iso(0, 0.0), (0, 1), op='drefuse')


def main():
    oracle = sys.argv[1]
    raw = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if raw:
        raw.mkdir(parents=True, exist_ok=True)
    lines = []
    for index, (op, kind, inputs) in enumerate(CASES):
        request = f'{op} {inputs}\n'
        run = subprocess.run(['prlimit', '--as=1073741824', '--cpu=20', 'timeout', '40', oracle],
                             input=request, capture_output=True, text=True)
        if raw:
            (raw / f'{index:04d}-{op}.txt').write_text(
                f'request: {request}exit: {run.returncode}\nstdout: {run.stdout}stderr: {run.stderr}')
        if run.returncode != 0:
            sys.exit(f'case {index} {op}: oracle exit {run.returncode}: {run.stderr.strip()}')
        response = run.stdout.split()
        assert response[0] == op, (index, response)
        lines.append(f'{op} {kind} | {inputs} | {" ".join(response[1:])}')
    header = ('# catallax-trace curve drop gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
