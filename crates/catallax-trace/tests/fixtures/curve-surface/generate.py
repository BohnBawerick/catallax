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
# SMLib answers with a range from its coincidence check, whose curve drop
# meets two analytic curves (AnalyticCurveIntersect): the port refuses.
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
# Straight curves (SmSurface::GlobalLineIntersect): lines across the bumps,
# the cylinder, the loft and the pole patch, through and beside them, along
# a cylinder ruling and tangent to it, over sub-intervals and parameter
# ranges other than [0, 1], and a collinear cubic.
def collinear(p, q):
    return curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [tuple(p[k] + (q[k] - p[k]) * i / 3 for k in range(3)) for i in range(4)])


def straight(p, q, a=0.0, b=1.0):
    return curve(1, [a, a, b, b], [p, q])


for s_, box in ((BUMP, BOX), (RBUMP, BOX)):
    global_cs('regular', s_, line((0.1, 0.2, -3), (2.9, 2.7, 3)), box, UNIT)
    global_cs('regular', s_, line((1.2, 1.7, -3), (1.2, 1.7, 3)), box, UNIT)
    global_cs('regular', s_, line((-0.5, 0.3, 0.5), (3.5, 0.3, 0.5)), box, UNIT)
    global_cs('regular', s_, line((-0.5, 1.4, 0.2), (3.5, 1.6, 0.4)), box, UNIT)
    global_cs('boundary', s_, line((0.1, 0.2, -3), (2.9, 2.7, 3)), box, (0.0, 0.4))
    global_cs('boundary', s_, line((1.2, 1.7, -3), (1.2, 1.7, 3)), (0.5, 0.5, 2.5, 2.5), UNIT)
    global_cs('regular', s_, collinear((0.3, 2.8, -2), (2.6, 0.4, 3)), box, UNIT)
    global_cs('regular', s_, straight((0.1, 0.2, -3), (2.9, 2.7, 3), 2, 5), box, (2, 5))
    global_cs('boundary', s_, line((-1, 1.5, -3), (-1, 1.5, 3)), box, UNIT)
global_cs('regular', CYLINDER, line((-2, 0.3, 0.5), (2, 0.3, 1.5)), CBOX, UNIT)
global_cs('regular', CYLINDER, line((-2, -2, 0.2), (2, 2, 1.8)), CBOX, UNIT)
global_cs('tangent', CYLINDER, line((1, -2, 1), (1, 2, 1)), CBOX, UNIT)
# 1e-6 inside the ruling SMLib answers with ranges where neighbouring
# meetings join along the line.
global_cs('near-singular', CYLINDER, line((1 - 1e-6, -2, 1), (1 - 1e-6, 2, 1)), CBOX, UNIT)
global_cs('regular', CYLINDER, line((0, 0, -1), (0, 0, 3)), CBOX, UNIT)
global_cs('regular', LOFT, line((-2, 0.3, 0.5), (2, 0.3, 1.5)), PBOX, UNIT)
global_cs('regular', LOFT, line((0.5, -1, 1.5), (0.5, 3, 1.5)), PBOX, UNIT)
global_cs('singular', POLE, line((0, 0, 3), (0, 0, -1)), PBOX, UNIT)
global_cs('regular', POLE, line((0.5, -1, 0.1), (0.5, 1, 0.6)), PBOX, UNIT)


# Lines placed from the bump's own surface points: vertical through a point
# on the knot line u = 1 (found in two patches, merged), ending on the
# surface at (1.3, 1.7), and skimming 1e-3 under the surface's height at
# (1.5, 1.5) along v, meeting it twice inside one patch.
def basis_at(knots, p, u):
    n = len(knots) - p - 1
    values = [1.0 if knots[i] <= u < knots[i + 1] else 0.0 for i in range(n)]
    if u >= knots[-1]:
        values[n - 1] = 1.0
    for d in range(1, p + 1):
        for i in range(n):
            a = (u - knots[i]) / (knots[i + d] - knots[i]) * values[i] if knots[i + d] > knots[i] else 0.0
            b = ((knots[i + d + 1] - u) / (knots[i + d + 1] - knots[i + 1]) * values[i + 1]
                 if i + 1 < n and knots[i + d + 1] > knots[i + 1] else 0.0)
            values[i] = a + b
    return values


def bump_point(u, v):
    k = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3]
    bu, bv = basis_at(k, 3, u), basis_at(k, 3, v)
    pt = [0.0, 0.0, 0.0]
    for i in range(6):
        for j in range(6):
            z = math.sin(0.9 * i + 0.3) * math.cos(0.7 * j - 0.3) + 0.05 * i
            for c, x in enumerate((i, j, z)):
                pt[c] += bu[i] * bv[j] * x
    return pt


KNOT = bump_point(1.0, 1.5)
global_cs('boundary', BUMP, line((KNOT[0], KNOT[1], -3), (KNOT[0], KNOT[1], 3)), BOX, UNIT)
END = bump_point(1.3, 1.7)
global_cs('touching', BUMP, line((END[0] + 0.5, END[1] - 0.2, END[2] + 2), END), BOX, UNIT)
TOP = bump_point(1.5, 1.5)
global_cs('near-singular', BUMP, line((TOP[0], TOP[1] - 1.2, TOP[2] - 1e-3), (TOP[0], TOP[1] + 1.2, TOP[2] - 1e-3)), BOX, UNIT)
global_cs('near-singular', BUMP, line((TOP[0] - 1.2, TOP[1], TOP[2] - 1e-2), (TOP[0] + 1.2, TOP[1], TOP[2] - 1e-2)), BOX, UNIT)
# Singularity tests at and near the pole and elsewhere.
for uv, kind in (((0.0, 0.5), 'singular'), ((0.0, 0.0), 'singular'), ((1e-9, 0.5), 'near-singular'),
                 ((1e-4, 0.5), 'near-singular'), ((0.5, 0.5), 'regular'), ((1.0, 0.5), 'boundary')):
    for point_only in (0, 1):
        sing(kind, POLE, uv, 1e-5, point_only)
for s_, uv in ((BUMP, (1.5, 1.5)), (CYLINDER, (0.0, 1.0)), (FLAT, (0.0, 0.0)), (LOFT, (0.0, 0.5))):
    sing('regular', s_, uv, 1e-5, 0)


# Curves lying on the bump (coincidence: SMLib drops them with
# DropAndTrimCurve and answers with ranges): its sides and interior
# parameter lines, over whole and partial intervals and a sub-domain.
def iso_basis(knots, p, u):
    n = len(knots) - p - 1
    span = max(i for i in range(p, n) if knots[i] <= u)
    values = [0.0] * n
    for i in range(n):
        values[i] = 1.0 if knots[i] <= u < knots[i + 1] or (i == span and knots[i] <= u <= knots[i + 1]) else 0.0
    for d in range(1, p + 1):
        for i in range(n):
            a_ = (u - knots[i]) / (knots[i + d] - knots[i]) * values[i] if knots[i + d] > knots[i] else 0.0
            b_ = ((knots[i + d + 1] - u) / (knots[i + d + 1] - knots[i + 1]) * values[i + 1]
                  if i + 1 < n and knots[i + d + 1] > knots[i + 1] else 0.0)
            values[i] = a_ + b_
    return values


def bump_iso(axis, value):
    k = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3]
    net = [[(i, j, math.sin(0.9 * i + 0.3) * math.cos(0.7 * j - 0.3) + 0.05 * i) for j in range(6)] for i in range(6)]
    b = iso_basis(k, 3, value)
    if axis == 0:
        points = [tuple(sum(b[i] * net[i][j][c] for i in range(6)) for c in range(3)) for j in range(6)]
    else:
        points = [tuple(sum(b[j] * net[i][j][c] for j in range(6)) for c in range(3)) for i in range(6)]
    return curve(3, k, points)


for axis, value, kind in ((0, 0.0, 'boundary'), (1, 3.0, 'boundary'), (0, 1.3, 'regular'), (1, 2.0, 'regular')):
    global_cs(kind, BUMP, bump_iso(axis, value), BOX, (0, 3))
    global_cs(kind, BUMP, bump_iso(axis, value), BOX, (0.5, 2.5))
    global_cs(kind, BUMP, bump_iso(axis, value), BOX, (0, 3), 1e-3)

# Range answers in the solution combiner (sm_GlobalCurveIntersectCombineSolutions).
# The solver's ranges come from the coincidence checks along patch sides that
# are at most C0: here the v = 0 edge and the crease of folded bilinear
# sheets, and the v = 0 edge of the bump. Each curve lies on such a side over
# a short stretch between the whole-curve coincidence checker's samples (so
# that check finds nothing) and leaves the surface elsewhere, rising off it or
# running outside the domain. See README.md for what each group covers.
FOLD = surface(1, 1, [0, 0, 1, 2, 2], [0, 0, 1, 2, 2],
               [(x, y, z) for x, z in ((0, 0), (1, 0), (2, 1)) for y in (0, 1, 2)])
FOLDV = surface(1, 1, [0, 0, 1, 2, 2], [0, 0, 1, 2, 2],
                [(x, y, z) for x in (0, 1, 2) for y, z in ((0, 0), (1, 0), (2, 1))])
FDOM = (0, 0, 2, 2)


def polyline(points, knots):
    return curve(1, knots, points)


SHORT = [0, 0, 0.53, 0.58, 1, 1]
LIFT = [0, 0, 0.3, 0.33, 0.36, 0.39, 0.42, 1, 1]
STEPS = [0, 0, 0.53, 0.54, 0.55, 0.56, 0.57, 1, 1]
OUT, UP = (-0.5, 0), (0, 0.5)


def edge_run(inner, knots, leave, low=0.1, high=0.9):
    """A polyline from beside the v = 0 edge, along `inner` points on or near
    it, and away again; `leave` offsets the outer ends in y or z."""
    dy, dz = leave
    return polyline([(low, dy, dz)] + inner + [(high, dy, dz)], knots)


for leave in (OUT, UP):
    global_cs('boundary', FOLDV, edge_run([(0.45, 0, 0), (0.55, 0, 0)], SHORT, leave), FDOM, UNIT)
    global_cs('boundary', FOLDV, edge_run([(0.95, 0, 0), (1.05, 0, 0)], SHORT, leave, 0.5, 1.5), FDOM, UNIT)
    global_cs('boundary', FOLDV, edge_run([(0.4, 0, 0), (0.45, 0, 0), (0.5, 0, 1e-3), (0.55, 0, 0), (0.6, 0, 0)], LIFT,
                                          leave), FDOM, UNIT)
    for lift in (7e-6, 3e-5):
        global_cs('boundary', FOLDV, edge_run([(0.45, 0, 0), (0.48, 0, 0), (0.5, 0, lift), (0.52, 0, 0),
                                               (0.55, 0, 0)], STEPS, leave), FDOM, UNIT)
    global_cs('tangent', FOLDV, edge_run([(0.45, 0, 5e-6), (0.55, 0, 5e-6)], SHORT, leave), FDOM, UNIT)
    global_cs('tangent', FOLDV, edge_run([(0.45, 0, 0), (0.55, 0, 8e-6)], SHORT, leave), FDOM, UNIT)
    global_cs('boundary', FOLDV, edge_run([(0.45, 0, 0), (0.5, 0, 0), (0.55, 0, 0)], [0, 0, 0.53, 0.555, 0.58, 1, 1],
                                          leave), FDOM, UNIT)
    global_cs('touching', FOLDV, edge_run([(0.5, 0, 0), (0.5, 0, 0)], SHORT, leave), FDOM, UNIT)
global_cs('boundary', FOLDV, edge_run([(0.45, 0, 0), (0.55, 0, 0)], SHORT, OUT), FDOM, UNIT, 1e-3)
global_cs('boundary', FOLDV, edge_run([(0.4, 0, 0), (0.45, 0, 0), (0.5, -1e-3, 0), (0.55, 0, 0), (0.6, 0, 0)], LIFT,
                                      OUT), FDOM, UNIT)
global_cs('boundary', FOLDV, edge_run([(0.45, 0, 0), (0.45002, 0, 0), (0.55, 0, 0)], [0, 0, 0.53, 0.53002, 0.58, 1, 1],
                                      OUT), FDOM, UNIT)
global_cs('boundary', FOLDV, edge_run([(0.5, 0, 0), (0.5 + 1e-9, 0, 0)], [0, 0, 0.53, 0.53 + 1e-9, 1, 1], UP), FDOM,
          UNIT)
global_cs('touching', FOLDV, curve(2, [0, 0, 0, 1, 1, 1], [(0.2, -0.01, 0), (0.5, 0.01, 0), (0.8, -0.01, 0)]), FDOM,
          UNIT)
global_cs('regular', FOLDV, polyline([(0.2, 1, 0.5), (1.06, 1, 0), (1.16, 1, 0), (1.9, 1, 0.5)], SHORT), FDOM, UNIT)
global_cs('regular', FOLDV, polyline([(0.2, 1, 0.5), (0.95, 1, 0), (1.05, 1, 0), (1.9, 1, 0.5)], SHORT), FDOM, UNIT)
global_cs('regular', FOLD, polyline([(1, 0.2, 0.5), (1, 1.06, 0), (1, 1.16, 0), (1, 1.9, 0.5)], SHORT), FDOM, UNIT)
# On the unit plane the rebuild's boundary search meets the plane's straight
# edges, whose point search SMLib solves analytically (SmLine): the port
# refuses, and SMLib's answer stays in the raw evidence.
global_cs('boundary', UNIT_PLANE, edge_run([(0.45, 0, 0), (0.55, 0, 0)], SHORT, OUT), PBOX, UNIT, op='grefuse')

# The bump's v = 0 edge: stretches cut exactly from the edge curve (its
# blossoms), kept, lifted between pieces, raised 5e-6 or tilted to 8e-6 off
# it, joined to cubic legs that rise off the bump or run outside the domain.
EDGE = [(i, 0.0, math.sin(0.9 * i + 0.3) * math.cos(-0.3) + 0.05 * i) for i in range(6)]
EDGE_KNOTS = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3]


def edge_blossom(args):
    span = max(i for i in range(3, 6) if EDGE_KNOTS[i] <= min(args) and EDGE_KNOTS[i + 1] >= max(args))
    d = [EDGE[span - 3 + i] for i in range(4)]
    for r in range(1, 4):
        for i in range(3, r - 1, -1):
            k = span - 3 + i
            a = (args[r - 1] - EDGE_KNOTS[k]) / (EDGE_KNOTS[k + 4 - r] - EDGE_KNOTS[k])
            d[i] = tuple((1 - a) * p + a * q for p, q in zip(d[i - 1], d[i]))
    return d[3]


def edge_piece(u0, u1, lift=(0.0, 0.0, 0.0, 0.0)):
    points = [edge_blossom((u0, u0, u0)), edge_blossom((u0, u0, u1)), edge_blossom((u0, u1, u1)),
              edge_blossom((u1, u1, u1))]
    return [(x, y, z + dz) for (x, y, z), dz in zip(points, lift)]


def leg(p, q, leave):
    dy, dz = leave
    return [p] + [(p[0] + f * (q[0] - p[0]), dy, p[2] + f * (q[2] - p[2]) + dz) for f in (1 / 3, 2 / 3)] + [q]


def bump_run(pieces, breaks, leave):
    dy, dz = leave
    low = (0.4, dy, EDGE[0][2] + 2 * dz)
    high = (2.6, dy, 0.6 + 2 * dz)
    segments = [leg(low, pieces[0][0], leave)] + pieces + [leg(pieces[-1][3], high, leave)]
    points = list(segments[0])
    for segment in segments[1:]:
        points += segment[1:]
    knots = [breaks[0]] * 4 + [b for b in breaks[1:-1] for _ in range(3)] + [breaks[-1]] * 4
    return curve(3, knots, points)


BUMP_OUT, BUMP_UP = (-0.4, 0), (0, 0.4)
for leave in (BUMP_OUT, BUMP_UP):
    global_cs('boundary', BUMP, bump_run([edge_piece(1.2, 1.5)], [0, 0.53, 0.58, 1], leave), BOX, UNIT)
    global_cs('boundary', BUMP, bump_run([edge_piece(1.2, 1.35), edge_piece(1.35, 1.5)], [0, 0.53, 0.555, 0.58, 1],
                                         leave), BOX, UNIT)
    for lift in (1.2e-5, 2e-5):
        global_cs('boundary', BUMP, bump_run([edge_piece(1.2, 1.3), edge_piece(1.3, 1.4, (0, lift, lift, 0)),
                                              edge_piece(1.4, 1.5)], [0, 0.53, 0.54, 0.56, 0.57, 1], leave), BOX, UNIT)
    global_cs('tangent', BUMP, bump_run([edge_piece(1.2, 1.5, (5e-6,) * 4)], [0, 0.53, 0.58, 1], leave), BOX, UNIT)
    global_cs('tangent', BUMP, bump_run([edge_piece(1.2, 1.5, (0, 2.7e-6, 5.3e-6, 8e-6))], [0, 0.53, 0.58, 1], leave),
              BOX, UNIT)
global_cs('boundary', BUMP, bump_run([edge_piece(0.6, 1.0), edge_piece(1.0, 2.0), edge_piece(2.0, 2.4)],
                                     [0, 0.2, 0.3, 0.7, 0.8, 1], BUMP_OUT), BOX, UNIT)


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
