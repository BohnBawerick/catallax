#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace surface intersector gate.

usage: generate.py ORACLE [RAW_DIR]

ORACLE is the driver built from oracle.cpp against the pinned SMLib. Each
case runs in a fresh process under a memory and CPU cap and a timeout. RAW_DIR,
when given, receives every request, response and error stream; keep it outside
the repository. The cases follow the eleven surface cases of the M25 scope
study's coverage driver, with variations.
"""
import importlib.util
import math
import pathlib
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('surfaces', HERE.parent / 'surface-solve' / 'generate.py')
surfaces = importlib.util.module_from_spec(spec)
spec.loader.exec_module(surfaces)
h, surface = surfaces.h, surfaces.surface
KINDS = ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
CASES = []


def v(*xs):
    return ' '.join(h(x) for x in xs)


def plane(origin, x, y):
    return f'p 1 {v(*origin, *x, *y, 1, 1, 0, 0, 1, 1)}'


def bump(amp, shift, rational=False, hill=False):
    points = []
    for i in range(6):
        for j in range(6):
            z = amp * math.sin(0.9 * i + shift) * math.cos(0.7 * j - shift) + 0.05 * i
            if hill:
                z = 1.0 if 2 <= i <= 3 and 2 <= j <= 3 else 0.0
            points.append((i, j, z, 1.0 + 0.25 * ((i + j) % 3)))
    k = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3]
    return 's ' + surface(3, 3, k, k, points, rational)


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


def case(kind, a, da, b, db, tolerance=None):
    assert kind in KINDS
    t = tolerance if tolerance is not None else 0.0
    for op in ('search', 'planes'):
        CASES.append((op, kind, f'{a} {v(*da)} {b} {v(*db)} {int(tolerance is not None)} {v(t)}'))


BOX, UNIT, PBOX = (0, 0, 3, 3), (0, 0, 1, 1), (-20, -20, 20, 20)
BUMP, RBUMP, HILL = bump(1.0, 0.3), bump(1.0, 0.3, rational=True), bump(1.0, 0.3, hill=True)
OTHER = bump(-0.8, 1.1)
k = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3]
b = basis(k, 3, 1.5)
peak = (b[2] + b[3]) ** 2
level = lambda z: plane((0, 0, z), (1, 0, 0), (0, 1, 0))
TILT = plane((0, 0, 0.35), (1, 0, 0.3), (0, 1, -0.2))
c60, s60 = 0.5, math.sqrt(3.0) / 2
CYL = 's ' + surface(2, 1, [0, 0, 0, 1, 1, 1], [0, 0, 1, 1],
                     [(x, y, 2.0 * j, w) for x, y, w in ((c60, -s60, 1.0), (2.0, 0.0, 0.5), (c60, s60, 1.0))
                      for j in range(2)], rational=True)
bulge = [1.0, 2.0, 0.5, 1.5]
LOFT = 's ' + surface(2, 3, [0, 0, 0, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1],
                      [(i, bulge[j] if i == 1 else 0.0, j, 0.75 if i == 1 else 1.0) for i in range(3) for j in range(4)],
                      rational=True)
SLOPE = plane((2.5, 0, 0), (1, -1, 0), (1, 1, -2))

# The coverage driver's cases.
case('regular', BUMP, BOX, level(0.35), PBOX)
case('regular', BUMP, BOX, TILT, PBOX)
case('regular', RBUMP, BOX, level(0.35), PBOX)
case('regular', BUMP, BOX, level(0.35), PBOX, 1e-6)
case('regular', HILL, BOX, level(0.3), PBOX)
case('tangent', HILL, BOX, level(peak), PBOX)
case('near-singular', HILL, BOX, level(peak - 1e-6), PBOX)
case('tangent', CYL, UNIT, plane((1.0, 0, 0), (0, 1, 0), (0, 0, 1)), PBOX)
case('regular', CYL, UNIT, plane((0.8, 0, 0), (0, 1, 0), (0, 0, 1)), PBOX)
case('regular', BUMP, BOX, OTHER, BOX)
case('regular', LOFT, UNIT, SLOPE, PBOX)
case('regular', LOFT, UNIT, SLOPE, PBOX, 1e-5)
# Variations: planes missing the surfaces, touching them at a corner or
# edge, sub-domains, a smaller plane domain, and swapped surfaces.
case('regular', BUMP, BOX, level(5.0), PBOX)
case('touching', HILL, BOX, level(1.0 + 1e-9), PBOX)
case('touching', HILL, BOX, level(0.0), PBOX)
case('boundary', BUMP, (0.5, 0.5, 2.5, 2.5), level(0.35), PBOX)
case('boundary', BUMP, BOX, level(0.35), (0, 0, 1.5, 1.5))
case('boundary', BUMP, (2.0, 2.0, 3.0, 3.0), level(0.35), PBOX)
case('regular', level(0.35), PBOX, BUMP, BOX)
case('touching', CYL, UNIT, plane((1.0 + 1e-6, 0, 0), (0, 1, 0), (0, 0, 1)), PBOX)
case('near-singular', CYL, UNIT, plane((1.0 - 1e-6, 0, 0), (0, 1, 0), (0, 0, 1)), PBOX)
case('singular', 's ' + surfaces.POLE, UNIT, level(1.0), PBOX)
case('regular', 's ' + surfaces.POLE, UNIT, level(0.5), PBOX)
for z in (-0.6, 0.0, 0.7, 1.05):
    case('regular', BUMP, BOX, level(z), PBOX)
    case('regular', RBUMP, BOX, level(z), PBOX, 1e-5)


# Side planes: flat patches side by side with a gap (disjoint), sharing an
# edge (grazing, which needs the curve-curve intersection and refuses),
# overlapping, at an angle along the shared edge, and face to face.
def flat(x0, x1, z0=0.0, z1=0.0, y0=0.0, y1=1.0):
    return 's ' + surface(1, 1, [0, 0, 1, 1], [0, 0, 1, 1], [(x0, y0, z0), (x0, y1, z0), (x1, y0, z1), (x1, y1, z1)])


for gap, kind in ((5e-4, 'near-singular'), (2e-2, 'regular'), (0.0, 'touching'), (-5e-4, 'regular'),
                  (1e-9, 'touching')):
    case(kind, flat(0, 1), UNIT, flat(1 + gap, 2), UNIT)
case('touching', flat(0, 1), UNIT, flat(1, 2, 0.0, 0.5), UNIT)
case('regular', flat(0, 1), UNIT, flat(0.5, 1.5, 1e-3, 1e-3), UNIT)
case('boundary', flat(0, 1), UNIT, flat(1, 2, 0.0, 0.0, 0.5, 1.5), UNIT)
# A patch beside the flat one, turned 0.005 and 0.02 radians about z: SMLib
# compares side planes only when their normals are within 1e-4 of opposed.
for angle in (0.005, 0.02):
    c, sn = math.cos(angle), math.sin(angle)
    corners = [(1 + 5e-4 + x * c - y * sn, x * sn + y * c, 0.0) for x in (0, 1) for y in (0, 1)]
    case('near-singular', flat(0, 1), UNIT, 's ' + surface(1, 1, [0, 0, 1, 1], [0, 0, 1, 1], corners), UNIT)
# A requested domain beyond the bump's natural one.
case('boundary', BUMP, (-1, -1, 4, 4), level(0.35), PBOX)


# Intersection point values: one pair alone, and a second pair with the
# first as the previous point, at interior points, on domain sides, near a
# tangency and with the step reversed.
def pcase(kind, a, da, b, db, first, second=None, tolerance=None):
    assert kind in KINDS
    t = tolerance if tolerance is not None else 0.0
    two = second is not None
    second = second or ((0, 0), (0, 0))
    CASES.append(('point', kind, f'{a} {v(*da)} {b} {v(*db)} {int(tolerance is not None)} {v(t)} {int(two)} '
                                 f'{v(*first[0], *first[1], *second[0], *second[1])}'))


for a, da, b, db, pairs in (
        (BUMP, BOX, level(0.35), PBOX, (((1.0, 1.5), (1.0, 1.5)), ((1.02, 1.52), (1.02, 1.52)))),
        (BUMP, BOX, TILT, PBOX, (((2.2, 0.7), (2.2, 0.7)), ((2.18, 0.73), (2.18, 0.73)))),
        (RBUMP, BOX, level(0.35), PBOX, (((0.6, 2.1), (0.6, 2.1)), ((0.58, 2.14), (0.58, 2.14)))),
        (CYL, UNIT, plane((0.8, 0, 0), (0, 1, 0), (0, 0, 1)), PBOX, (((0.3, 0.5), (-0.5, 1.0)), ((0.3, 0.52), (-0.5, 1.04)))),
        (LOFT, UNIT, SLOPE, PBOX, (((0.4, 0.5), (0.3, 0.2)), ((0.42, 0.5), (0.31, 0.2)))),
        (BUMP, BOX, OTHER, BOX, (((1.5, 1.5), (1.4, 1.6)), ((1.52, 1.5), (1.42, 1.6))))):
    first, second = pairs
    pcase('regular', a, da, b, db, first)
    pcase('regular', a, da, b, db, first, second)
    pcase('regular', a, da, b, db, second, first)
# On the sides of the domains, where the tangent must point inward.
for uv in ((0.0, 1.5), (3.0, 1.2), (1.4, 0.0), (2.1, 3.0), (0.0, 0.0), (3.0, 3.0)):
    pcase('boundary', BUMP, BOX, level(0.35), PBOX, (uv, uv))
    pcase('boundary', BUMP, BOX, TILT, PBOX, (uv, uv), ((uv[0] * 0.98 + 0.03, uv[1] * 0.98 + 0.03), uv))
pcase('boundary', CYL, UNIT, plane((0.8, 0, 0), (0, 1, 0), (0, 0, 1)), PBOX, ((0.0, 0.5), (-0.5, 1.0)))
pcase('boundary', CYL, UNIT, plane((0.8, 0, 0), (0, 1, 0), (0, 0, 1)), PBOX, ((1.0, 1.0), (0.5, 2.0)))
# Normals within three degrees of parallel after a previous point: walls
# 1e-3 and 5e-4 inside the cylinder's ruling, where the point is refined by
# the normal solve.
for x in (0.999, 0.9995):
    w = plane((x, 0, 0), (0, 1, 0), (0, 0, 1))
    pcase('near-singular', CYL, UNIT, w, PBOX, ((0.48, 0.5), (0.02, 1.0)), ((0.47, 0.52), (0.03, 1.04)))
    pcase('near-singular', CYL, UNIT, w, PBOX, ((0.47, 0.52), (0.03, 1.04)))
pcase('regular', BUMP, BOX, level(0.35), PBOX, ((1.0, 1.5), (1.0, 1.5)), tolerance=1e-6)


def main():
    oracle = sys.argv[1]
    raw = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if raw:
        raw.mkdir(parents=True, exist_ok=True)
    lines = []
    for index, (op, kind, inputs) in enumerate(CASES):
        request = f'{op} {inputs}\n'
        run = subprocess.run(['prlimit', '--as=1073741824', '--cpu=10', 'timeout', '20', oracle],
                             input=request, capture_output=True, text=True)
        if raw:
            (raw / f'{index:04d}-{op}.txt').write_text(
                f'request: {request}exit: {run.returncode}\nstdout: {run.stdout}stderr: {run.stderr}')
        if run.returncode != 0:
            sys.exit(f'case {index} {op}: oracle exit {run.returncode}: {run.stderr.strip()}')
        response = run.stdout.split()
        assert response[0] == op, (index, response)
        lines.append(f'{op} {kind} | {inputs} | {" ".join(response[1:])}')
    header = ('# catallax-trace surface intersector gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
