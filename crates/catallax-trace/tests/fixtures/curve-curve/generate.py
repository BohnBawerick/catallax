#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace curve-curve gate.

usage: generate.py ORACLE [RAW_DIR]

ORACLE is the driver built from oracle.cpp against the pinned SMLib. Each
case runs in a fresh process under a memory and CPU cap and a timeout. RAW_DIR,
when given, receives every request, response and error stream; keep it outside
the repository.
"""
import importlib.util
import math
import pathlib
import random
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('curves', HERE.parent / 'curve-global' / 'generate.py')
curves = importlib.util.module_from_spec(spec)
spec.loader.exec_module(curves)
h, v, curve = curves.h, curves.v, curves.curve
KINDS = ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
CASES = []
UNIT = (0, 1)


def cc(kind, a, ia, b, ib, tolerance=1e-5, flag=0, op='global'):
    assert kind in KINDS
    CASES.append((op, kind, f'{a} {v(*ia)} {b} {v(*ib)} {v(tolerance)} {flag}'))


def cubic(points, rational=False):
    return curve(3, [0, 0, 0, 0, 1, 1, 1, 1], points, rational)


def shifted(points, d):
    return [(p[0] + d[0], p[1] + d[1], p[2] + d[2]) + tuple(p[3:]) for p in points]


WAVE_P = [(0.2, 0.3, -2), (1.0, 1.5, 4), (2.0, 1.0, -3), (2.8, 2.6, 2)]
WAVE = cubic(WAVE_P)
SNAKE_P = [(0, 0, 0), (1, 2, 0), (2, -2, 0), (3, 0, 0)]
SNAKE = cubic(SNAKE_P)
AXIS = cubic([(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)])
LONG_AXIS = cubic([(-1, 0, 0), (0.5, 0, 0), (2, 0, 0), (3.5, 0, 0)])
PARABOLA = curve(2, [0, 0, 0, 1, 1, 1], [(-1, 1, 0), (0, -1, 0), (1, 1, 0)])
FLAT = cubic([(-2, 0, 0), (-0.5, 0, 0), (0.5, 0, 0), (2, 0, 0)])
LOOP = curves.LOOP
POINT = curves.POINT
SPIRAL = curves.SPIRAL
RWAVE = curves.RWAVE
CIRCLE = curves.CIRCLE
r = math.sqrt(0.5)
ARC = curve(2, [0, 0, 0, 1, 1, 1], [(1, 0, 0, 1), (1, 1, 0, r), (0, 1, 0, 1)], rational=True)
ARC2 = curve(2, [0, 0, 0, 1, 1, 1], [(2, 0.5, 0, 1), (2, 1.5, 0, r), (1, 1.5, 0, 1)], rational=True)
SEG1 = curve(1, [0, 0, 1, 1], [(0, 0, 0), (2, 1, 0)])
SEG2 = curve(1, [0, 0, 1, 1], [(0, 1, 0), (2, 0, 0)])
CUSP = cubic([(0, 0, 0), (2, 1, 0), (2, 1, 0), (0, 2, 0)])
TWO_SPAN = curve(3, [0, 0, 0, 0, 1, 2, 2, 2, 2], [(0, -1, 0), (0.5, 1, 0), (1.5, 1, 0), (2.5, -1, 0), (3, 0.5, 0)])
KINK = curve(3, [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
             [(0, 1, 0), (0.3, 0.5, 0), (0.6, 0.2, 0), (1, 0, 0), (1.4, 0.2, 0), (1.7, 0.5, 0), (2, 1, 0)])

# Crossings: a snake across an axis at three points (two of them its ends),
# a wave and a rational wave against curves through them, and two-span and
# kinked curves (a C0 knot) across the axis.
cc('regular', SNAKE, UNIT, AXIS, UNIT)
cc('regular', AXIS, UNIT, SNAKE, UNIT)
cc('regular', SNAKE, UNIT, LONG_AXIS, UNIT)
cc('regular', SNAKE, (0.2, 0.9), LONG_AXIS, UNIT)
cc('regular', SNAKE, UNIT, AXIS, (0.1, 0.6))
cc('regular', WAVE, UNIT, cubic([(0.2, 0.3, 1), (1.0, 1.5, 0), (2.0, 1.0, 0), (2.8, 2.6, -1)]), UNIT)
cc('regular', RWAVE, UNIT, cubic([(1.4, 0, 0), (1.4, 1, 0.5), (1.4, 2, 1), (1.4, 3, 1.5)]), UNIT)
cc('regular', TWO_SPAN, (0, 2), AXIS, UNIT)
cc('regular', KINK, (0, 2), cubic([(0, 0.1, 0), (1, 0.1, 0), (2, 0.1, 0), (3, 0.1, 0)]), UNIT)
cc('regular', SPIRAL, (0, 3), cubic([(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3)]), UNIT)
cc('regular', ARC, UNIT, AXIS, UNIT)
# Where the search drops a point onto the arc, SMLib hands it to the
# analytic circle (`DropPointFast`).
cc('regular', ARC, UNIT, cubic([(0, 0, 0), (0.5, 0.5, 0), (1, 1, 0), (1.5, 1.5, 0)]), UNIT)
cc('regular', LOOP, (0, 3), cubic([(-2, 0.2, 0), (-0.5, 0.2, 0), (0.5, 0.2, 0), (2, 0.2, 0)]), UNIT)
cc('regular', LOOP, (0, 3), cubic([(-2, 0, 0), (-0.5, 0, 0), (0.5, 0, 0), (2, 0, 0)]), UNIT)
# Misses: apart, 1e-4 apart and 1e-6 apart (inside the tolerance) at a
# crossing point.
cc('regular', SNAKE, UNIT, cubic(shifted([(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)], (0, 0, 1))), UNIT)
cc('near-singular', SNAKE, UNIT, cubic(shifted([(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)], (0, 0, 1e-4))), UNIT)
cc('near-singular', SNAKE, UNIT, cubic(shifted([(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)], (0, 0, 1e-6))), UNIT)
# Tangent meetings: the parabola on a line through its vertex, 1e-7 and
# 1e-4 below it, and two arcs touching; with a looser tolerance.
for d, kind in ((0.0, 'tangent'), (-1e-7, 'tangent'), (-1e-4, 'near-singular'), (1e-7, 'tangent')):
    for tolerance in (1e-5, 1e-3):
        cc(kind, PARABOLA, UNIT, cubic(shifted([(-2, 0, 0), (-0.5, 0, 0), (0.5, 0, 0), (2, 0, 0)], (0, d, 0))),
           UNIT, tolerance)
cc('tangent', ARC, UNIT, cubic([(1, -1, 0), (1, -0.3, 0), (1, 0.3, 0), (1, 1, 0)]), UNIT)
cc('tangent', ARC, UNIT, cubic([(1, -1, 0), (1, -0.3, 0), (1, 0.3, 0), (1, 1, 0)]), UNIT, 1e-3)
# Two arcs: SMLib meets them analytically (`AnalyticCurveIntersect`).
cc('tangent', ARC, UNIT, ARC2, UNIT, op='grefuse')
# Touching at ends: end to end, an end on the other's interior, an end
# 1e-7 off the other, ends turning back on each other.
cc('touching', AXIS, UNIT, cubic([(3, 0, 0), (3, 1, 0), (4, 1, 0), (4, 2, 0)]), UNIT)
cc('touching', AXIS, UNIT, cubic([(1.5, 0, 0), (1.5, 1, 0), (2.5, 1, 0), (2.5, 2, 0)]), UNIT)
cc('touching', AXIS, UNIT, cubic([(1.5, 1e-7, 0), (1.5, 1, 0), (2.5, 1, 0), (2.5, 2, 0)]), UNIT)
cc('touching', AXIS, UNIT, cubic([(3, 0, 0), (2, 0.5, 0), (1, 0.5, 0), (0, 1, 0)]), UNIT)
cc('touching', TWO_SPAN, (0, 1), TWO_SPAN, (1, 2))
# Coincidence: a copy of the wave, the reversed wave, the wave against its
# middle, two collinear lines overlapping by half, a copy 1e-7 away (within
# the tolerance) and 1e-3 away, the loop against its copy, and with the
# coincidence check skipped.
REV = cubic(list(reversed(WAVE_P)))
for kind, a, ia, b, ib in (
        ('boundary', WAVE, UNIT, cubic(WAVE_P), UNIT),
        ('boundary', WAVE, UNIT, REV, UNIT),
        ('boundary', WAVE, UNIT, cubic(WAVE_P), (0.25, 0.75)),
        ('boundary', WAVE, (0.25, 0.75), cubic(WAVE_P), UNIT),
        ('boundary', AXIS, UNIT, cubic([(1.5, 0, 0), (2.5, 0, 0), (3.5, 0, 0), (4.5, 0, 0)]), UNIT),
        ('boundary', AXIS, UNIT, cubic([(4.5, 0, 0), (3.5, 0, 0), (2.5, 0, 0), (1.5, 0, 0)]), UNIT),
        ('near-singular', WAVE, UNIT, cubic(shifted(WAVE_P, (0, 0, 1e-7))), UNIT),
        ('near-singular', WAVE, UNIT, cubic(shifted(WAVE_P, (0, 0, 1e-3))), UNIT),
        ('boundary', LOOP, (0, 3), curve(3, [0, 0, 0, 0, 1, 2, 3, 3, 3, 3],
                                         [(1, 0, 0), (1, 1, 0), (-1, 1, 0), (-1, -1, 0), (1, -1, 0), (1, 0, 0)]),
         (0, 3)),
        ('boundary', SPIRAL, (0, 3), curves.SPIRAL, (1, 2))):
    cc(kind, a, ia, b, ib)
    cc(kind, a, ia, b, ib, flag=1)
    cc(kind, a, ia, b, ib, op='coincide')
    cc(kind, a, ia, b, ib, flag=1, op='coincide')
# Partial coincidence the end-point check cannot settle: polylines and
# cubics that share a straight middle and part at both ends (a range from
# the span search, merged across spans), the shared part 1e-7 and 1e-3
# apart, reversed, and with one curve's interval ending inside it.
def poly(points):
    n = len(points)
    return curve(1, [0, 0] + list(range(1, n - 1)) + [n - 1, n - 1], points)


Y1 = poly([(0, 1, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0), (4, 1, 0)])
Y2 = poly([(0, -1, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0), (4, -1, 0)])
Y2R = poly([(4, -1, 0), (3, 0, 0), (2, 0, 0), (1, 0, 0), (0, -1, 0)])
for d, kind in ((0.0, 'boundary'), (1e-7, 'near-singular'), (1e-3, 'near-singular')):
    y2 = poly([(0, -1, d), (1, 0, d), (2, 0, d), (3, 0, d), (4, -1, d)])
    cc(kind, Y1, (0, 4), y2, (0, 4))
    cc(kind, Y1, (0, 4), y2, (0, 4), op='coincide')
cc('boundary', Y1, (0, 4), Y2R, (0, 4))
cc('boundary', Y1, (0.5, 2.5), Y2, (0, 4))
cc('boundary', Y1, (0, 4), Y2, (1.5, 3.5))
SHARE1 = curve(3, [0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4],
               [(0, 1, 0), (0.5, 0.5, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0), (3.5, 0.5, 0), (4, 1, 0)])
SHARE2 = curve(3, [0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4],
               [(0, -1, 0), (0.5, -0.5, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0), (3.5, -0.5, 0), (4, -1, 0)])
cc('boundary', SHARE1, (0, 4), SHARE2, (0, 4))
cc('boundary', SHARE1, (0, 4), SHARE2, (0, 4), op='coincide')
cc('boundary', SHARE1, (0, 4), SHARE2, (0, 4), 1e-3)
# Near-tangent stretches within the tolerance (ranges from the span search
# that are one tangent point): flat parabolas and a flat cubic on a line,
# and a line crossing another at a very small angle.
for a_, tolerance in ((1e-4, 1e-3), (1e-4, 1e-2), (1e-3, 1e-2), (1e-6, 1e-5)):
    flat = curve(2, [0, 0, 0, 1, 1, 1], [(-3, 9 * a_, 0), (0, -9 * a_, 0), (3, 9 * a_, 0)])
    cc('tangent', flat, UNIT, FLAT, UNIT, tolerance)
    cc('tangent', FLAT, UNIT, flat, UNIT, tolerance)
cc('tangent', cubic([(-2, 1e-3, 0), (-0.7, -1e-3, 0), (0.7, -1e-3, 0), (2, 1e-3, 0)]), UNIT, FLAT, UNIT, 1e-2)
cc('tangent', cubic([(-2, 2e-3, 0), (-0.7, 0, 0), (0.7, 0, 0), (2, 2e-3, 0)]), UNIT, FLAT, UNIT, 1e-2)
for slope, tolerance in ((1e-4, 1e-3), (1e-3, 1e-2), (1e-5, 1e-5)):
    tilted = cubic([(-2, -2 * slope, 0), (-0.5, -0.5 * slope, 0), (0.5, 0.5 * slope, 0), (2, 2 * slope, 0)])
    cc('near-singular', tilted, UNIT, FLAT, UNIT, tolerance)
# Coincident stretches the span search finds apart: a polyline on an axis
# polyline except for a bump of 0.3, 1e-3 and 1e-6 between them, leaving
# with or without coming back across, and leaving at a small angle.
AXIS5 = poly([(x, 0, 0) for x in range(6)])
for bump, kind in ((0.3, 'boundary'), (1e-3, 'boundary'), (1e-6, 'near-singular')):
    cc(kind, poly([(0, 0, 0), (1, 0, 0), (2, 0, 0), (2.5, bump, 0), (3, 0, 0), (4, 0, 0), (5, 0, 0)]), (0, 6),
       AXIS5, (0, 5))
    cc(kind, poly([(0, 1, 0), (1, 0, 0), (2, 0, 0), (2.5, bump, 0), (3, 0, 0), (4, 0, 0), (5, 1, 0)]), (0, 6),
       AXIS5, (0, 5))
cc('boundary', poly([(0, 1, 0), (1, 0, 0), (2, 0, 0), (3, 0.5, 0), (4, -0.5, 0)]), (0, 4), AXIS5, (0, 5))
cc('boundary', poly([(0, 1, 0), (1, 0, 0), (2, 0, 0), (4, 1e-3, 0), (5, -1, 0)]), (0, 4), AXIS5, (0, 5))
cc('boundary', poly([(0, 1, 0), (1, 0, 0), (2, 0, 0), (4, 1e-3, 0), (5, -1, 0)]), (0, 4), AXIS5, (0, 5), 1e-2)
# A parabola as a seven-span quadratic B-spline against a polyline on its
# tangent: within the tolerance over a stretch, a tangent point.
for a_, tolerance in ((1e-3, 1e-3), (1e-3, 1e-4), (1e-2, 1e-3), (1e-4, 1e-4)):
    knots = [-3, -3, -3, -2, -1, 0, 1, 2, 3, 3, 3]
    para = curve(2, knots, [((knots[i + 1] + knots[i + 2]) / 2, a_ * knots[i + 1] * knots[i + 2], 0)
                            for i in range(8)])
    line = poly([(x, 0, 0) for x in range(-3, 4)])
    cc('tangent', para, (-3, 3), line, (0, 6), tolerance)
    cc('tangent', line, (0, 6), para, (-3, 3), tolerance)
    shifted_line = poly([(x + 0.5, 0, 0) for x in range(-3, 4)])
    cc('tangent', para, (-3, 3), shifted_line, (0, 6), tolerance)
# A half parabola from its vertex along its tangent: a range that is one
# tangent point at its tight end.
for a_, tolerance in ((1e-3, 1e-3), (1e-2, 1e-3)):
    knots = [0, 0, 0, 1, 2, 3, 3, 3]
    half = curve(2, knots, [((knots[i + 1] + knots[i + 2]) / 2, a_ * knots[i + 1] * knots[i + 2], 0)
                            for i in range(5)])
    cc('tangent', half, (0, 3), poly([(x, 0, 0) for x in range(-3, 4)]), (0, 6), tolerance)
    cc('tangent', poly([(x, 0, 0) for x in range(-3, 4)]), (0, 6), half, (0, 3), tolerance)
# Curves sharing a stretch that covers neither: the end-point check finds
# the range and the search trims it off both intervals.
cc('boundary', poly([(0, 0, 0), (2, 0, 0), (3, 1, 0)]), (0, 2), poly([(-1, 1, 0), (0, 0, 0), (1, 0, 0), (1.5, 0, 0)]),
   (0, 3))
cc('boundary', poly([(-1, 1, 0), (0, 0, 0), (1, 0, 0), (1.5, 0, 0)]), (0, 3), poly([(0, 0, 0), (2, 0, 0), (3, 1, 0)]),
   (0, 2))
# Quadratic pieces whose ends meet with parallel tangents.
cc('touching', curve(2, [0, 0, 0, 1, 1, 1], [(0, 0, 0), (1, 1, 0), (2, 0, 0)]), UNIT,
   curve(2, [0, 0, 0, 1, 1, 1], [(2, 0, 0), (3, -1, 0), (3, -2, 0)]), UNIT)
cc('touching', curve(2, [0, 0, 0, 1, 1, 1], [(0, 0, 0), (1, 1, 0), (2, 0, 0)]), UNIT,
   curve(2, [0, 0, 0, 1, 1, 1], [(3, -2, 0), (3, -1, 0), (2, 0, 0)]), UNIT)
# Crossings at C0 knots of both curves and at a knot of one.
cc('touching', KINK, (0, 2), poly([(0, 0, 0), (1, 0, 0), (2, 0.5, 0)]), (0, 2))
cc('touching', KINK, (0, 2), cubic([(1, -1, 0), (1, -0.3, 0), (1, 0.3, 0), (1, 1, 0)]), UNIT)
cc('touching', poly([(0, 1, 0), (1, 0, 0), (2, 1, 0)]), (0, 2), poly([(0, -1, 0), (1, 0, 0), (2, -1, 0)]), (0, 2))
# Short pieces of low degree that end within the tolerance of each other.
cc('touching', curve(2, [0, 0, 0, 1, 1, 1], [(0, 0, 0), (0.5, 0.5, 0), (1, 0, 0)]), UNIT,
   curve(2, [0, 0, 0, 1, 1, 1], [(1, 1e-6, 0), (1.5, -0.5, 0), (2, 0, 0)]), UNIT)
cc('touching', curve(2, [0, 0, 0, 1, 1, 1], [(0, 0, 0), (0.5, 0.5, 0), (1, 0, 0)]), UNIT,
   curve(2, [0, 0, 0, 1, 1, 1], [(1, 0, 0), (1.5, 0.5, 0), (2, 0, 0)]), UNIT)
cc('touching', poly([(0, 0, 0), (1, 0, 0)]), UNIT, poly([(1, 0, 0), (2, 0, 0)]), UNIT)
# A point curve at an end of the other.
cc('singular', POINT, UNIT, cubic([(1, 1, 0.5), (2, 1, 0.5), (3, 2, 0.5), (4, 2, 0.5)]), UNIT)
# Singular curves: a point curve on and off the wave, a cusp crossed at
# and beside its cusp.
cc('singular', POINT, UNIT, cubic([(0, 0, 0.5), (1, 1, 0.5), (2, 2, 0.5), (3, 3, 0.5)]), UNIT)
cc('singular', POINT, UNIT, cubic([(0, 0, 0.6), (1, 1, 0.6), (2, 2, 0.6), (3, 3, 0.6)]), UNIT)
cc('singular', CUSP, UNIT, cubic([(1.5, 0, 0), (1.5, 0.7, 0), (1.5, 1.3, 0), (1.5, 2, 0)]), UNIT)
cc('singular', CUSP, UNIT, cubic([(1.0, 0, 0), (1.0, 0.7, 0), (1.0, 1.3, 0), (1.0, 2, 0)]), UNIT)
# SMLib meets analytic line pairs directly; arc pairs remain refused.
cc('regular', SEG1, UNIT, SEG2, UNIT)
ARC_Z1 = curve(2, [0, 0, 0, 1, 1, 1], [(1.5, 0, 1, 1), (1.5, 1, 1, r), (0.5, 1, 1, 1)], rational=True)
cc('regular', ARC_Z1, UNIT, CIRCLE, (0, 4), op='grefuse')
# A fixed pseudo-random sweep of cubic pairs in a box, at two tolerances.
rng = random.Random(20261112)
for i in range(40):
    # Mostly planar pairs, which meet; every fifth pair leaves the plane.
    lift = 0.5 if i % 5 == 0 else 0.0
    a = cubic([(rng.uniform(0, 3), rng.uniform(0, 3), rng.uniform(-lift, lift)) for _ in range(4)])
    b = cubic([(rng.uniform(0, 3), rng.uniform(0, 3), 0.0) for _ in range(4)])
    cc('regular', a, UNIT, b, UNIT, 1e-5 if i % 3 else 1e-2)



# SmLine::IntersectWithLine through the cache's analytic-curve dispatch.
def segment(start, end, natural=UNIT):
    return curve(1, [natural[0], natural[0], natural[1], natural[1]], [start, end])


def line_pair(kind, a, b, ia=UNIT, ib=UNIT, tolerance=1e-5, op='global'):
    for flag in (0, 1):
        cc(kind, a, ia, b, ib, tolerance, flag, op)
        cc(kind, b, ib, a, ia, tolerance, flag, op)


base = segment((0, 0, 0), (2, 0, 0))
for start, end in [((1, -1, 0), (1, 1, 0)),
                   ((0, -1, 0), (0, 1, 0)),
                   ((2, 0, 0), (3, 1, 0)),
                   ((2, 0, 0), (3, 0, 0)),
                   ((1, 0, 0), (3, 0, 0)),
                   ((3, 0, 0), (1, 0, 0)),
                   ((0, 0, 0), (2, 0, 0)),
                   ((2, 0, 0), (0, 0, 0)),
                   ((-1, 0, 0), (3, 0, 0)),
                   ((1, -1, 1), (1, 1, 1))]:
    line_pair('touching', base, segment(start, end))
for gap in (0., 1e-12, 1e-8, 1e-5, 1.0001e-5, 1e-3):
    line_pair('tangent', base, segment((1, gap, 0), (3, gap, 0)))
    line_pair('tangent', base, segment((1, gap, 0), (3, -gap, 0)))
for ia, ib in [((.25, .75), (.25, .75)), ((0., .25), UNIT),
               ((.5, 1.), (.5, 1.)), ((.5, .5), UNIT),
               (UNIT, (.5, .5))]:
    line_pair('singular' if ia[0] == ia[1] or ib[0] == ib[1] else 'boundary',
              base, segment((1, -1, 0), (1, 1, 0)), ia, ib)
for length in (1e-10, 1e-8, 1e-6):
    line_pair('near-singular', segment((0, 0, 0), (length, 0, 0)),
              segment((length / 2, -length, 0), (length / 2, length, 0)),
              tolerance=1e-12, op='known-miss' if length < 1e-6 else 'global')
for ia, ib in [((2., 5.), (-2., 3.)), ((-5., -2.), (4., 8.))]:
    line_pair('boundary', segment((0, 0, 0), (2, 1, 0), ia),
              segment((0, 1, 0), (2, 0, 0), ib), ia, ib)
line_rng = random.Random(2026111309)
for _ in range(40):
    center = [line_rng.uniform(-3, 3) for _ in range(3)]
    va = [line_rng.uniform(-2, 2) for _ in range(3)]
    vb = [line_rng.uniform(-2, 2) for _ in range(3)]
    a = segment([p-v for p,v in zip(center, va)], [p+v for p,v in zip(center, va)])
    b = segment([p-v for p,v in zip(center, vb)], [p+v for p,v in zip(center, vb)])
    line_pair('regular', a, b)


def main():
    oracle = sys.argv[1]
    raw = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if raw:
        raw.mkdir(parents=True, exist_ok=True)
    lines = []
    for index, (op, kind, inputs) in enumerate(CASES):
        oracle_op = 'global' if op == 'known-miss' else op
        request = f'{oracle_op} {inputs}\n'
        run = subprocess.run(['prlimit', '--as=1073741824', '--cpu=10', 'timeout', '20', oracle],
                             input=request, capture_output=True, text=True)
        if raw:
            (raw / f'{index:04d}-{op}.txt').write_text(
                f'request: {request}exit: {run.returncode}\nstdout: {run.stdout}stderr: {run.stderr}')
        if run.returncode != 0:
            sys.exit(f'case {index} {op}: oracle exit {run.returncode}: {run.stderr.strip()}')
        response = run.stdout.split()
        assert response[0] == oracle_op, (index, response)
        lines.append(f'{op} {kind} | {inputs} | {" ".join(response[1:])}')
    header = ('# catallax-trace curve-curve gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
