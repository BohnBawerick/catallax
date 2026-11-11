#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace curve-surface gate.

usage: generate.py ORACLE [RAW_DIR]

ORACLE is the driver built from oracle.cpp against the pinned SMLib. Each
case runs in a fresh process under a memory and CPU cap and a timeout. RAW_DIR,
when given, receives every request, response and error stream; keep it outside
the repository.
"""
import math
import pathlib
import struct
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
KINDS = ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
CASES = []


def h(x):
    return format(struct.unpack('<Q', struct.pack('<d', float(x)))[0], '016x')


def surface(du, dv, ku, kv, points, rational=False):
    nu, nv = len(ku) - du - 1, len(kv) - dv - 1
    assert len(points) == nu * nv
    tokens = [h(k) for k in ku] + [h(k) for k in kv]
    for p in points:
        tokens += [h(p[0]), h(p[1]), h(p[2]), h(p[3] if len(p) > 3 else 1.0)]
    return f'{du} {dv} {nu} {nv} {int(rational)} ' + ' '.join(tokens)


def bump(rational):
    points = []
    for i in range(6):
        for j in range(6):
            z = math.sin(0.9 * i + 0.3) * math.cos(0.7 * j - 0.3) + 0.05 * i
            points.append((i, j, z, 1.0 + 0.25 * ((i + j) % 3)))
    k = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3]
    return surface(3, 3, k, k, points, rational)


BUMP = bump(False)
RBUMP = bump(True)
bulge = [1.0, 2.0, 0.5, 1.5]
LOFT = surface(2, 3, [0, 0, 0, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1],
               [(i, bulge[j] if i == 1 else 0.0, j, 0.75 if i == 1 else 1.0) for i in range(3) for j in range(4)],
               rational=True)
# A cubic patch whose u = 0 edge collapses to one point.
POLE = surface(3, 3, [0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1],
               [(0, 0, 1) if i == 0 else (i, j - 1.5, 0.3 * i * j) for i in range(4) for j in range(4)])
FLAT = surface(1, 1, [0, 0, 2, 2], [0, 0, 1, 1], [(0, 0, 0), (0, 1, 0), (2, 0, 0), (2, 1, 1)])
# A cylinder of radius 1 about z, closed in u with its seam at u = 0 and 4.
r = math.sqrt(0.5)
ring = [(1, 0, 1), (1, 1, r), (0, 1, 1), (-1, 1, r), (-1, 0, 1), (-1, -1, r), (0, -1, 1), (1, -1, r), (1, 0, 1)]
CYLINDER = surface(2, 1, [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4], [0, 0, 2, 2],
                   [(x, y, z, w) for x, y, w in ring for z in (0, 2)], rational=True)
MINIMIZE, MAXIMIZE, NORMALIZE = 3, 4, 5


def curve(degree, knots, points, rational=False):
    tokens = [h(k) for k in knots]
    for p in points:
        tokens += [h(p[0]), h(p[1]), h(p[2]), h(p[3] if len(p) > 3 else 1.0)]
    return f'{degree} {len(points)} {int(rational)} ' + ' '.join(tokens)


def line(p, q):
    return curve(1, [0, 0, 1, 1], [p, q])


def v(*xs):
    return ' '.join(h(x) for x in xs)


def intersect(kind, s, c, domain, interval, tolerance, t, uv):
    assert kind in KINDS
    CASES.append(('intersect', kind, f'{s} {c} {v(*domain)} {v(*interval)} {v(tolerance, t, *uv)}'))


def solve(kind, s, c, domain, interval, operation, t, uv):
    assert kind in KINDS
    CASES.append(('solve', kind, f'{s} {c} {v(*domain)} {v(*interval)} {operation} {v(1e-6, t, *uv)}'))


def sing(kind, s, uv, tolerance, point_only):
    assert kind in KINDS
    CASES.append(('sing', kind, f'{s} {v(*uv, tolerance)} {point_only}'))


UNIT = (0, 1)
BOX = (0, 0, 3, 3)
WAVE = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0.2, 0.3, -2), (1.0, 1.5, 4), (2.0, 1.0, -3), (2.8, 2.6, 2)])
for s in (BUMP, RBUMP):
    intersect('regular', s, line((1.2, 1.7, -3), (1.2, 1.7, 3)), BOX, UNIT, 1e-5, 0.5, (1.0, 1.5))
    intersect('regular', s, line((0.1, 0.2, -3), (2.9, 2.7, 3)), BOX, UNIT, 1e-5, 0.4, (1.0, 1.0))
    intersect('regular', s, WAVE, BOX, UNIT, 1e-5, 0.3, (0.8, 1.2))
    intersect('regular', s, WAVE, BOX, UNIT, 1e-5, 0.7, (2.2, 1.6))
    intersect('boundary', s, line((-1, 1.5, -3), (-1, 1.5, 3)), BOX, UNIT, 1e-5, 0.5, (0.0, 1.5))
    intersect('boundary', s, line((1.5, 1.5, 5), (1.5, 1.5, 9)), BOX, UNIT, 1e-5, 0.0, (1.5, 1.5))
    # A line along the surface normal direction missing it, and one skimming it.
    intersect('regular', s, line((0.5, 0.5, 4), (2.5, 2.5, 4)), BOX, UNIT, 1e-5, 0.5, (1.5, 1.5))
    for t, uv in ((0.5, (1.5, 1.5)), (0.1, (0.3, 2.0))):
        solve('regular', s, WAVE, BOX, UNIT, MINIMIZE, t, uv)
        solve('regular', s, WAVE, BOX, UNIT, NORMALIZE, t, uv)
        solve('regular', s, line((0.5, 0.5, 4), (2.5, 2.5, 4)), BOX, UNIT, MINIMIZE, t, uv)
    solve('regular', s, line((0.5, 0.5, 4), (2.5, 2.5, 4)), BOX, UNIT, MAXIMIZE, 0.5, (1.5, 1.5))
# Lines tangent to the cylinder, along it, just off it, and through it.
CBOX = (0, 0, 4, 2)
for offset, kind in ((0.0, 'tangent'), (1e-8, 'near-singular'), (-1e-8, 'touching'), (-0.3, 'regular')):
    tangent_line = line((1 + offset, -2, 1), (1 + offset, 2, 1))
    intersect(kind, CYLINDER, tangent_line, CBOX, UNIT, 1e-5, 0.5, (0.1, 1.0))
    intersect(kind, CYLINDER, tangent_line, CBOX, UNIT, 1e-5, 0.45, (3.9, 0.9))
    solve(kind, CYLINDER, tangent_line, CBOX, UNIT, MINIMIZE, 0.5, (0.1, 1.0))
intersect('singular', CYLINDER, line((1, 0, -1), (1, 0, 3)), CBOX, UNIT, 1e-5, 0.4, (0.0, 0.6))
intersect('touching', CYLINDER, line((0, -2, 0.5), (0, 2, 0.5)), CBOX, UNIT, 1e-5, 0.2, (2.9, 0.5))
# Lines into the pole of POLE.
PBOX = (0, 0, 1, 1)
for p, q, t, uv, kind in (((0, 0, 3), (0, 0, -1), 0.5, (0.05, 0.5), 'singular'),
                          ((-0.5, 0.3, 2), (0.5, -0.3, 0), 0.5, (0.02, 0.4), 'singular'),
                          ((0.5, -1, 1), (0.5, 1, 1), 0.5, (0.2, 0.5), 'near-singular'),
                          ((2, 0, 3), (2, 0, -3), 0.5, (0.6, 0.5), 'regular')):
    intersect(kind, POLE, line(p, q), PBOX, UNIT, 1e-5, t, uv)
    solve(kind, POLE, line(p, q), PBOX, UNIT, MINIMIZE, t, uv)
intersect('regular', LOFT, line((1, -1, 0.5), (1, 3, 0.5)), PBOX, UNIT, 1e-5, 0.6, (0.5, 0.5))
intersect('touching', FLAT, line((0, 0.5, 0), (2, 0.5, 0.5)), (0, 0, 2, 1), UNIT, 1e-5, 0.3, (0.6, 0.5))
# Coincidence: a ruling and a circle of the cylinder, a line half on it,
# copies 1e-6 off it, and a crossing line.
r = math.sqrt(0.5)
CIRCLE = curve(2, [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4],
               [(1, 0, 1, 1), (1, 1, 1, r), (0, 1, 1, 1), (-1, 1, 1, r), (-1, 0, 1, 1), (-1, -1, 1, r),
                (0, -1, 1, 1), (1, -1, 1, r), (1, 0, 1, 1)], rational=True)
ARC = curve(2, [0, 0, 0, 1, 1, 1], [(0, 1, 0.5, 1), (-1, 1, 0.5, r), (-1, 0, 0.5, 1)], rational=True)
OFF_ARC = curve(2, [0, 0, 0, 1, 1, 1], [(0, 1 + 1e-6, 0.5, 1), (-1 - 1e-6, 1 + 1e-6, 0.5, r),
                                        (-1 - 1e-6, 0, 0.5, 1)], rational=True)
HALF = curve(1, [0, 0, 0.5, 1, 1], [(-1, 0, 0), (-1, 0, 1), (-2, 0, 2)])


def coincide(kind, s, c, domain, interval, tolerance, partial):
    assert kind in KINDS
    CASES.append(('coincide', kind, f'{s} {c} {v(*domain)} {v(*interval)} {v(tolerance)} {partial}'))


for partial in (0, 1):
    coincide('singular', CYLINDER, line((-1, 0, 0), (-1, 0, 2)), CBOX, UNIT, 1e-5, partial)
    coincide('singular', CYLINDER, CIRCLE, CBOX, (0, 4), 1e-5, partial)
    coincide('singular', CYLINDER, ARC, CBOX, UNIT, 1e-5, partial)
    coincide('near-singular', CYLINDER, OFF_ARC, CBOX, UNIT, 1e-5, partial)
    coincide('regular', CYLINDER, OFF_ARC, CBOX, UNIT, 1e-7, partial)
    coincide('touching', CYLINDER, HALF, CBOX, UNIT, 1e-5, partial)
    coincide('regular', CYLINDER, line((-2, 0.1, 1), (2, 0.1, 1)), CBOX, UNIT, 1e-5, partial)
    coincide('regular', BUMP, WAVE, BOX, UNIT, 1e-5, partial)
# Global intersections over whole domains: cubics across and along the
# bumps, the loft and the cylinder, a cubic tangent to a bump, a curve into
# the pole, a point curve, and straight or coincident curves (which need the
# surface point drop and line intersection of later steps).


def global_cs(kind, s, c, domain, interval, tolerance=1e-5, op='global'):
    assert kind in KINDS
    CASES.append((op, kind, f'{s} {c} {v(*domain)} {v(*interval)} {v(tolerance)}'))


ARCH = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0.2, 0.3, 3), (1.0, 1.5, -2), (2.0, 1.0, -2), (2.8, 2.6, 3)])
LOOP = curve(3, [0, 0, 0, 0, 0.5, 1, 1, 1, 1], [(-1.5, -1.5, 0.2), (2.5, -1.5, 0.5), (2.5, 2.5, 1.0),
                                                 (-1.5, 2.5, 1.5), (-1.5, -1.5, 1.8)])
SPIRAL = curve(3, [0, 0, 0, 0, 1, 2, 3, 3, 3, 3], [(0, 2, 0.1), (2, 1, 0.5), (1, -2, 0.9), (-2, -1, 1.3),
                                                    (-1, 2, 1.7), (1.5, 0.5, 1.9)])
POINT = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(1.0, 1.0, 0.5)] * 4)
for s_, kind in ((BUMP, 'regular'), (RBUMP, 'regular')):
    global_cs(kind, s_, WAVE, BOX, UNIT)
    global_cs(kind, s_, ARCH, BOX, UNIT)
    global_cs(kind, s_, ARCH, BOX, (0.2, 0.8))
    global_cs('boundary', s_, ARCH, (0.5, 0.5, 2.5, 2.5), UNIT)
global_cs('regular', CYLINDER, LOOP, CBOX, UNIT)
global_cs('regular', CYLINDER, SPIRAL, CBOX, (0, 3))
global_cs('boundary', CYLINDER, SPIRAL, (0.5, 0, 3.5, 2), (0, 3))
global_cs('regular', LOFT, curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(1, -1, 0.2), (1, 1, 0.4), (1, 2, 0.6), (1, 4, 0.8)]),
          PBOX, UNIT)
# Through the pole: SMLib's answer comes from the singular point
# inversion's curve drop.
global_cs('singular', POLE, curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0, 0, 3), (0.1, 0, 1.5), (-0.1, 0, 0.5), (0, 0, -1)]),
          PBOX, UNIT)
global_cs('regular', BUMP, POINT, BOX, UNIT)
UNIT_PLANE = surface(1, 1, [0, 0, 1, 1], [0, 0, 1, 1], [(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0)])
global_cs('regular', UNIT_PLANE, curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0.25, 0.75, 0)] * 4), PBOX, UNIT)
CORNER = math.sin(0.3) * math.cos(-0.3)
global_cs('touching', BUMP, curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0.0, 0.0, CORNER)] * 4), BOX, UNIT)
CROSSER = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(-2, 0.3, 0.5), (-0.7, 0.2, 0.8), (0.7, 0.5, 1.2), (2, 0.3, 1.5)])
global_cs('regular', CYLINDER, CROSSER, CBOX, UNIT)
global_cs('boundary', CYLINDER, CROSSER, (0.5, 0, 3.5, 2), UNIT)
global_cs('regular', CYLINDER, CROSSER, CBOX, (0.3, 0.7))
# Near tangency: a slightly bent curve 0.999 from the cylinder axis meets it
# twice close together; curves crossing a knot line of the bump, its v = 0
# edge from outside, and its corner.
NEAR = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0.999, -2, 1), (0.999, -0.7, 1.001), (0.999, 0.7, 1.0), (0.999, 2, 1)])
global_cs('tangent', CYLINDER, NEAR, CBOX, UNIT)
global_cs('near-singular', CYLINDER,
          curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0.99999, -2, 1), (0.99999, -0.7, 1.0001), (0.99999, 0.7, 1.0),
                                              (0.99999, 2, 1)]), CBOX, UNIT)
for s_ in (BUMP, RBUMP):
    global_cs('boundary', s_, curve(3, [0, 0, 0, 0, 1, 1, 1, 1],
                                   [(1.0, 1.5, -3), (1.0, 1.52, -1), (1.0, 1.49, 1), (1.0, 1.5, 3)]), BOX, UNIT)
    global_cs('boundary', s_, curve(3, [0, 0, 0, 0, 1, 1, 1, 1],
                                   [(1.5, -1, 0.0), (1.5, -0.4, 0.3), (1.5, 0.3, 0.2), (1.5, 1, 0.5)]), BOX, UNIT)
    global_cs('boundary', s_, curve(3, [0, 0, 0, 0, 1, 1, 1, 1],
                                   [(-1, -1, 0.5), (-0.3, -0.4, 0.4), (0.2, 0.3, 0.2), (1, 1, 0.1)]), BOX, UNIT)
# A bent curve through the bump's v = 0 edge at u = 1.5, exactly (the edge
# point from SMLib's evaluation), and an arc 1e-4 off the cylinder's circle,
# meeting it twice 0.02 radians apart.
EDGE_Z = 0.5820453225519064
global_cs('touching', BUMP, curve(2, [0, 0, 0, 1, 1, 1], [(2.3, -1, EDGE_Z + 0.3), (2.7, 0, EDGE_Z),
                                                       (2.3, 1, EDGE_Z - 0.3)]), BOX, UNIT)
R, D = 0.9999, 0.0001
c45 = math.sqrt(0.5)
CLOSE_ARC = curve(2, [0, 0, 0, 1, 1, 1], [(D + R * c45, -R * c45, 1, 1), (D + R / c45, 0, 1, c45),
                                        (D + R * c45, R * c45, 1, 1)], rational=True)
# SMLib answers with a range, whose classification needs the surface point
# drop: the port refuses.
global_cs('tangent', CYLINDER, CLOSE_ARC, CBOX, UNIT, op='grefuse')
# A fixed pseudo-random sweep of cubic curves against the bumps, the loft
# and the cylinder. REFUSED lists the sweep cases the port refuses (each needs
# the surface point drop).
state = 4242
REFUSED = set()


def uniform(low, high):
    global state
    state = (1103515245 * state + 12345) % (1 << 31)
    return low + (high - low) * state / (1 << 31)


sweep = 0
for s_, box, lows, highs in ((BUMP, BOX, (-0.5, -0.5, -2), (3.5, 3.5, 3)), (RBUMP, BOX, (-0.5, -0.5, -2), (3.5, 3.5, 3)),
                             (LOFT, PBOX, (-0.5, -1, -0.5), (2.5, 3, 3.5)),
                             (CYLINDER, CBOX, (-2, -2, -0.5), (2, 2, 2.5))):
    for _ in range(12):
        points = [tuple(uniform(a, b) for a, b in zip(lows, highs)) for _ in range(4)]
        c = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], points)
        global_cs('regular', s_, c, box, UNIT, op='grefuse' if sweep in REFUSED else 'global')
        sweep += 1
# Singularity tests at and near the pole and elsewhere.
for uv, kind in (((0.0, 0.5), 'singular'), ((0.0, 0.0), 'singular'), ((1e-9, 0.5), 'near-singular'),
                 ((1e-4, 0.5), 'near-singular'), ((0.5, 0.5), 'regular'), ((1.0, 0.5), 'boundary')):
    for point_only in (0, 1):
        sing(kind, POLE, uv, 1e-5, point_only)
for s_, uv in ((BUMP, (1.5, 1.5)), (CYLINDER, (0.0, 1.0)), (FLAT, (0.0, 0.0)), (LOFT, (0.0, 0.5))):
    sing('regular', s_, uv, 1e-5, 0)


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
    header = ('# catallax-trace curve-surface gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
