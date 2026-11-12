#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace global point-curve gate.

usage: generate.py ORACLE [RAW_DIR]

ORACLE is the driver built from oracle.cpp against the pinned SMLib. Each
case runs in a fresh process under a memory and CPU cap and a timeout. RAW_DIR,
when given, receives every request, response and error stream; keep it outside
the repository.
"""
import math
import pathlib
import random
import struct
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
KINDS = ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
CASES = []
INTERSECT, MINIMIZE, MAXIMIZE, NORMALIZE = 0, 3, 4, 5
SINGLE, ALL = 0, 1


def h(x):
    return format(struct.unpack('<Q', struct.pack('<d', float(x)))[0], '016x')


def v(*xs):
    return ' '.join(h(x) for x in xs)


def curve(degree, knots, points, rational=False):
    tokens = [h(k) for k in knots]
    for p in points:
        tokens += [h(p[0]), h(p[1]), h(p[2]), h(p[3] if len(p) > 3 else 1.0)]
    return f'{degree} {len(points)} {int(rational)} ' + ' '.join(tokens)


def solve(kind, c, interval, operation, test, tolerance=1e-5, requested=ALL, op='global'):
    assert kind in KINDS
    CASES.append((op, kind, f'{c} {v(*interval)} {operation} {v(*test, tolerance)} {requested}'))


def drop(kind, c, interval, test, inward=None, guess=None, tolerance=1e-5, operation=MINIMIZE, op='drop'):
    assert kind in KINDS
    flags = f'{int(inward is not None)} {v(*(inward or (0, 0, 0)))} {int(guess is not None)} {v(guess or 0.0)}'
    CASES.append((op, kind, f'{c} {v(*interval)} {operation} {v(*test, tolerance)} {flags}'))


UNIT = (0, 1)
WAVE = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0.2, 0.3, -2), (1.0, 1.5, 4), (2.0, 1.0, -3), (2.8, 2.6, 2)])
SPIRAL = curve(3, [0, 0, 0, 0, 1, 2, 3, 3, 3, 3], [(0, 2, 0.1), (2, 1, 0.5), (1, -2, 0.9), (-2, -1, 1.3),
                                                    (-1, 2, 1.7), (1.5, 0.5, 1.9)])
# A closed planar cubic whose ends meet with matching tangents.
LOOP = curve(3, [0, 0, 0, 0, 1, 2, 3, 3, 3, 3], [(1, 0, 0), (1, 1, 0), (-1, 1, 0), (-1, -1, 0), (1, -1, 0),
                                                  (1, 0, 0)])
# A closed cubic with a corner where its ends meet.
CORNER_LOOP = curve(3, [0, 0, 0, 0, 1, 2, 2, 2, 2], [(0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0), (0, 0, 0)])
POLYLINE = curve(1, [0, 0, 0.5, 1, 1], [(0, 0, 0), (1, 0, 0), (1, 1, 0)])
PARABOLA = curve(2, [0, 0, 0, 1, 1, 1], [(-1, 1, 0), (0, -1, 0), (1, 1, 0)])
STRAIGHT = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0, 0, 0), (1, 0.5, 0), (2, 1, 0), (3, 1.5, 0)])
POINT = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(1.0, 1.0, 0.5)] * 4)
SEGMENT = curve(1, [0, 0, 1, 1], [(0, 0, 0), (2, 1, 0)])
r = math.sqrt(0.5)
CIRCLE = curve(2, [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4],
               [(1, 0, 1, 1), (1, 1, 1, r), (0, 1, 1, 1), (-1, 1, 1, r), (-1, 0, 1, 1), (-1, -1, 1, r),
                (0, -1, 1, 1), (1, -1, 1, r), (1, 0, 1, 1)], rational=True)
# A rational cubic, which SMLib keeps as a spline.
RWAVE = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0.2, 0.3, -2, 1), (1.0, 1.5, 4, 2), (2.0, 1.0, -3, 0.5),
                                             (2.8, 2.6, 2, 1)], rational=True)

# Nearest, farthest, normal and meeting points from fixed points: on the
# curves, at and beyond their ends, at a corner, on the axis of the loops.
for c, interval, points in (
        (WAVE, UNIT, ((1.4, 1.2, 0.0), (0.2, 0.3, -2), (3.0, 3.0, 3.0), (-1, 0, 0), (1.5, 1.5, 10))),
        (SPIRAL, (0, 3), ((0, 0, 1), (0, 0, 5), (0, 2, 0.1), (1.5, 0.5, 1.9), (3, 3, 0))),
        (LOOP, (0, 3), ((0, 0, 0), (0, 0, 1), (1, 0, 0), (1.5, 0, 0), (0.3, -0.2, 0.4), (1.5, -0.3, 0),
                         (1.5, 0.3, 0), (1.2, 0.05, 0.2), (0.6, -0.1, 0))),
        (CORNER_LOOP, (0, 2), ((0, 0, 0), (-0.5, -0.5, 0), (1, 1, 0), (-1, 0.2, 0.3), (0.3, -0.1, 0),
                                (-0.1, 0.3, 0))),
        (POLYLINE, UNIT, ((1, 0, 0), (1.5, -0.5, 0), (0.5, 0.5, 0), (2, 0.5, 0), (0.9, 0.1, 0))),
        (PARABOLA, UNIT, ((0, 0, 0), (0, 2, 0), (0, -2, 0), (2, 2, 0))),
        (STRAIGHT, UNIT, ((1, 0.5, 0), (1, 2, 0), (-1, 0, 0), (4, 2, 1))),
        (RWAVE, UNIT, ((1.4, 1.2, 0.0), (3.0, 3.0, 3.0), (-1, 0, 0)))):
    for test in points:
        for operation in (MINIMIZE, MAXIMIZE, NORMALIZE, INTERSECT):
            kind = 'regular'
            if operation == INTERSECT:
                kind = 'touching'
            solve(kind, c, interval, operation, test)
        solve('regular', c, interval, NORMALIZE, test, requested=SINGLE)
# Points just off the seam of the smooth loop, along its tangent and across it.
for test in ((1, -5e-6, 0), (1, 5e-6, 0), (1 + 5e-6, -5e-6, 0)):
    for operation in (MINIMIZE, INTERSECT):
        solve('near-singular', LOOP, (0, 3), operation, test)
# Partial intervals, one ending at a corner.
for c, interval, test in ((WAVE, (0.2, 0.7), (3.0, 3.0, 3.0)), (WAVE, (0.2, 0.7), (1.4, 1.2, 0.0)),
                          (SPIRAL, (0.5, 2.5), (0, 0, 1)), (POLYLINE, (0, 0.5), (1.5, -0.5, 0)),
                          (POLYLINE, (0.25, 0.75), (2, 0.5, 0)), (LOOP, (0.5, 2.5), (1.5, 0, 0))):
    for operation in (MINIMIZE, MAXIMIZE, NORMALIZE):
        solve('boundary', c, interval, operation, test)
# Equidistant: the centre of the loop and of the circle, and points on a
# straight curve, where every parameter is an answer.
solve('singular', CIRCLE, (0, 4), MAXIMIZE, (0, 0, 1))
solve('singular', CIRCLE, (0, 4), NORMALIZE, (0, 0, 1))
solve('singular', CIRCLE, (0, 4), NORMALIZE, (0, 0, 3))
solve('singular', CIRCLE, (0, 4), NORMALIZE, (0.5, 0.2, 1))
solve('regular', CIRCLE, (0, 4), MAXIMIZE, (0.5, 0.2, 1))
solve('singular', SEGMENT, UNIT, NORMALIZE, (1, 0.5, 0))
solve('regular', SEGMENT, UNIT, NORMALIZE, (0, 1, 0))
solve('regular', SEGMENT, UNIT, MAXIMIZE, (0, 1, 0))
solve('near-singular', POINT, UNIT, MINIMIZE, (1, 1, 0.5))
solve('near-singular', POINT, UNIT, MINIMIZE, (1, 2, 0.5))
solve('near-singular', POINT, UNIT, NORMALIZE, (1, 2, 0.5))
solve('near-singular', POINT, UNIT, INTERSECT, (1, 1, 0.5 + 1e-7))
# Tolerance edges: a point 1e-6 and 1e-4 off the wave meets it within 1e-5
# only in the first case.
solve('tangent', WAVE, UNIT, INTERSECT, (1.4, 1.2, 0.0), 1e-4)
# Two-point lines SMLib solves with its analytic SmLine: the foot of the
# point clamped to the interval, refined on the spline when off it. Points
# on, beside, beyond either end and on the line's extension, sub-intervals
# that clamp the foot, meeting tolerances either side of the distance, a
# line over the parameter range [2, 5], a slanted line far from the origin,
# and a line shorter than the scaled zero (searched as a spline).
SLANT = curve(1, [2, 2, 5, 5], [(1e3, -2e3, 5e2), (1e3 + 3.7, -2e3 + 1.1, 5e2 - 2.9)])
TINY = curve(1, [0, 0, 1, 1], [(1, 1, 1), (1 + 1e-13, 1, 1)])
for test in ((1, 0.5, 0), (0, 1, 0), (3, 0, 0), (-1, -1, 0), (4, 2, 0), (1, 0.5, 1e-6), (1, 0.5, 1e-4)):
    solve('regular', SEGMENT, UNIT, MINIMIZE, test)
    solve('touching', SEGMENT, UNIT, INTERSECT, test)
    solve('boundary', SEGMENT, (0.2, 0.7), MINIMIZE, test)
    solve('boundary', SEGMENT, (0.2, 0.7), INTERSECT, test, 1e-3)
for test in ((1e3 + 1, -2e3 + 0.3, 5e2 - 0.8), (1e3, -2e3, 5e2), (1e3 + 9, -2e3, 5e2), (1e3 + 1.85, -2e3 + 0.55, 5e2 - 1.45)):
    solve('regular', SLANT, (2, 5), MINIMIZE, test)
    solve('touching', SLANT, (2, 5), INTERSECT, test)
    solve('boundary', SLANT, (3, 4), MINIMIZE, test)
solve('singular', TINY, UNIT, MINIMIZE, (1, 1, 1))
solve('singular', TINY, UNIT, INTERSECT, (1, 2, 1))
# A rational two-point line whose weights differ is not evenly parameterised,
# fails SMLib's sampled line test and is searched as a spline. SMLib solves
# an evenly weighted rational two-point line with the sampled test's SmLine,
# which the port refuses.
RSEGMENT = curve(1, [0, 0, 1, 1], [(0, 0, 0, 1), (2, 1, 0, 3)], rational=True)
REVEN = curve(1, [0, 0, 1, 1], [(0, 0, 0, 1), (2, 1, 0, 1)], rational=True)
solve('regular', RSEGMENT, UNIT, MINIMIZE, (1, 1, 0))
solve('touching', RSEGMENT, UNIT, INTERSECT, (1, 0.5, 0))
solve('regular', REVEN, UNIT, MINIMIZE, (1, 1, 0), op='grefuse')
# Circular arcs SMLib solves with its analytic SmCircle
# (SmBSplineCurve::CreateAnalyticCurve, SmEllipse::DropPointFast): the test
# point's angle about the circle's axis, clamped to the arc, converted to the
# spline parameter through the polar converter, then refined on the spline
# when the point is off the circle. A point on the axis answers both ends.


def arc(center, xaxis, yaxis, radius, start, sweep, spans, knots=None, scale=None):
    """A rational quadratic arc of `spans` equal conic spans with double
    interior knots; `knots` gives the span boundaries (default 0..spans) and
    `scale` reweights each span's middle and end (a Moebius change of
    parameter that keeps the circle)."""
    def at(degrees):
        a = math.radians(degrees)
        return tuple(c + radius * (math.cos(a) * x + math.sin(a) * y)
                     for c, x, y in zip(center, xaxis, yaxis))
    step = sweep / spans
    half = math.radians(step / 2)
    points = [at(start) + (1.0,)]
    for i in range(spans):
        a0 = start + i * step
        mid_dir = at(a0 + step / 2)
        corner = tuple(c + (m - c) / math.cos(half) for c, m in zip(center, mid_dir))
        points.append(corner + (math.cos(half),))
        points.append(at(a0 + step) + (1.0,))
    if scale is not None:
        points = [(p[0], p[1], p[2], p[3] * scale ** i) for i, p in enumerate(points)]
    bounds = knots or list(range(spans + 1))
    kv = [bounds[0]] * 3
    for b in bounds[1:-1]:
        kv += [b, b]
    kv += [bounds[-1]] * 3
    return curve(2, kv, points, rational=True)


def insert_once(knots, points, u):
    """Boehm insertion of `u` once into a rational quadratic, in homogeneous
    coordinates; returns the new knots and Euclidean points with weights."""
    p = 2
    k = max(i for i in range(len(knots) - 1) if knots[i] <= u < knots[i + 1])
    hom = [(x * w, y * w, z * w, w) for x, y, z, w in points]
    new = []
    for i in range(len(hom) + 1):
        if i <= k - p:
            new.append(hom[i])
        elif i >= k + 1:
            new.append(hom[i - 1])
        else:
            a = (u - knots[i]) / (knots[i + p] - knots[i])
            new.append(tuple(a * q + (1 - a) * r for q, r in zip(hom[i], hom[i - 1])))
    kv = knots[:k + 1] + [u] + knots[k + 1:]
    return kv, [(x / w, y / w, z / w, w) for x, y, z, w in new]


s2 = math.sqrt(0.5)
TILT_X = (s2, 0.0, s2)
TILT_Y = (0.0, 1.0, 0.0)
QUARTER = arc((1, 2, 3), TILT_X, TILT_Y, 2.0, 0.0, 90.0, 1)
HALF = arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 1.5, 30.0, 180.0, 2, knots=[2, 3.5, 5])
THREE_QUARTER = arc((-1, 0.5, 0), (0, 1, 0), (-1, 0, 0), 0.75, 0.0, 270.0, 3, knots=[0, 0.25, 0.5, 1])
SKEWED = arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 1.0, 0.0, 120.0, 2, scale=1.7)
# A one-span arc of 120 degrees with a knot inserted once: SMLib doubles it.
_k, _p = insert_once([0, 0, 0, 1, 1, 1],
                     [(1, 0, 0, 1), (1, math.sqrt(3), 0, 0.5), (-0.5, math.sqrt(3) / 2, 0, 1)], 0.4)
SINGLE_KNOT = curve(2, _k, _p, rational=True)
# A quarter ellipse is a rational quadratic but no arc: searched as a spline.
ELLIPSE = curve(2, [0, 0, 0, 1, 1, 1], [(2, 0, 0, 1), (2, 1, 0, r), (0, 1, 0, 1)], rational=True)

for test in ((2, 0, 1), (0.5, 0.5, 2), (-3, 0.2, 1), (0, -1, 1), (1, 1e-7, 1), (1, -1e-7, 1),
             (1.000001, 0, 1), (0.2, -0.9, 0.7)):
    solve('regular', CIRCLE, (0, 4), MINIMIZE, test)
    solve('touching', CIRCLE, (0, 4), INTERSECT, test)
for test in ((1, 1e-7, 1), (1, -1e-7, 1), (1, 0, 1)):
    solve('near-singular', CIRCLE, (0, 4), MINIMIZE, test, requested=SINGLE)
solve('singular', CIRCLE, (0, 4), MINIMIZE, (0, 0, 1))
solve('singular', CIRCLE, (0, 4), MINIMIZE, (0, 0, 1), requested=SINGLE)
solve('singular', CIRCLE, (0, 4), INTERSECT, (0, 0, 1), 2.0)
solve('tangent', CIRCLE, (0, 4), INTERSECT, (0, 1 + 1e-6, 1), 1e-5)
solve('tangent', CIRCLE, (0, 4), INTERSECT, (0, 1 + 1e-6, 1), 1e-7)
for interval in ((0.5, 1.5), (1.2, 3.9), (3, 4)):
    for test in ((2, 0, 1), (-1, -1, 1), (0, 0, 1)):
        solve('boundary', CIRCLE, interval, MINIMIZE, test)
        solve('boundary', CIRCLE, interval, INTERSECT, test, 0.5)
for c, interval, points in (
        (QUARTER, (0, 1), ((1 + 2 * s2, 2, 3 + 2 * s2), (1, 4, 3), (5, 3, 3), (1, 0, 3), (1 - 1, 2, 3 + 1),
                           (1, 2, 3), (1 + 3 * s2, 2.5, 3 + 3 * s2))),
        (HALF, (2, 5), ((0, 3, 0), (1, 1, 0.5), (0, -2, 0), (2, -0.5, 0), (0, 0, 2))),
        (THREE_QUARTER, (0, 1), ((-1, 1.5, 0), (0, 0, 0), (-1.5, -0.5, 0.3), (-0.5, 0.5, 0), (-1, 0.5, -1))),
        (SKEWED, (0, 2), ((2, 0.1, 0), (0, 2, 0), (-1, -1, 0), (0.3, 0.3, 0.1))),
        (SINGLE_KNOT, (0, 1), ((2, 1, 0), (-1, 1, 0), (0.5, -2, 0), (0, 0, 0)))):
    for test in points:
        solve('regular', c, interval, MINIMIZE, test)
        solve('touching', c, interval, INTERSECT, test, 1e-3)
        solve('regular', c, interval, MINIMIZE, test, requested=SINGLE)
solve('boundary', HALF, (2.5, 4), MINIMIZE, (0, 3, 0))
solve('boundary', HALF, (2.5, 4), MINIMIZE, (2, -0.5, 0))
solve('boundary', THREE_QUARTER, (0.1, 0.6), MINIMIZE, (-1.5, -0.5, 0.3))
for test in ((2, 0.5, 0), (0, 0, 0), (1.5, 1.5, 0)):
    solve('regular', ELLIPSE, (0, 1), MINIMIZE, test)
    solve('touching', ELLIPSE, (0, 1), INTERSECT, test)

# A fixed pseudo-random sweep.
rng = random.Random(9)
for c, interval, low, high in ((WAVE, UNIT, (-1, -1, -4), (4, 4, 5)), (SPIRAL, (0, 3), (-3, -3, -1), (3, 3, 3)),
                               (LOOP, (0, 3), (-2, -2, -1), (2, 2, 1)), (RWAVE, UNIT, (-1, -1, -4), (4, 4, 5)),
                               (PARABOLA, UNIT, (-2, -2, -1), (2, 2, 1))):
    for _ in range(8):
        test = tuple(rng.uniform(a, b) for a, b in zip(low, high))
        operation = rng.choice((MINIMIZE, MINIMIZE, MAXIMIZE, NORMALIZE))
        solve('regular', c, interval, operation, test)

# Point drops: from good and poor guesses, without a guess (the global
# solve), meeting within a tolerance, normal drops, and snaps to an end in
# and against the curve's inward direction.
for c, interval, test, guess in ((WAVE, UNIT, (1.4, 1.2, 0.0), 0.5), (WAVE, UNIT, (1.4, 1.2, 0.0), 0.0),
                                 (WAVE, UNIT, (3.0, 3.0, 3.0), 0.1), (SPIRAL, (0, 3), (0, 0, 1), 1.5),
                                 (LOOP, (0, 3), (1.5, 0, 0), 1.5), (POLYLINE, UNIT, (1.5, -0.5, 0), 0.2)):
    drop('regular', c, interval, test, guess=guess)
    drop('regular', c, interval, test)
    drop('regular', c, interval, test, guess=guess, operation=NORMALIZE)
    drop('regular', c, interval, test, operation=NORMALIZE)
    drop('touching', c, interval, test, guess=guess, operation=INTERSECT)
    drop('touching', c, interval, test, operation=INTERSECT, tolerance=0.5)
for inward, kind in (((0.8, 1.2, 6), 'boundary'), ((-0.8, -1.2, -6), 'boundary')):
    drop(kind, WAVE, UNIT, (0.2, 0.3, -2 + 1e-7), inward=inward)
    drop(kind, WAVE, UNIT, (0.2, 0.3, -2 + 1e-7), inward=inward, guess=0.01)
    drop(kind, WAVE, UNIT, (2.8, 2.6, 2 - 1e-7), inward=inward)
    drop(kind, WAVE, UNIT, (2.8, 2.6, 2 - 1e-7), inward=inward, guess=0.99)
    drop(kind, LOOP, (0, 3), (1, 1e-7, 0), inward=inward)
drop('near-singular', POINT, UNIT, (1, 2, 0.5))
drop('near-singular', POINT, UNIT, (1, 2, 0.5), guess=0.5)
# Without a guess the drop runs the global solve: the line's and the
# circle's analytic drops.
drop('regular', SEGMENT, UNIT, (1, 1, 0))
drop('regular', SEGMENT, UNIT, (3, 1, 0))
drop('touching', SEGMENT, UNIT, (1, 0.5, 2e-6), operation=INTERSECT)
drop('regular', CIRCLE, (0, 4), (2, 0.5, 1))
drop('touching', CIRCLE, (0, 4), (0, 1 + 1e-6, 1), operation=INTERSECT)
drop('regular', QUARTER, (0, 1), (5, 3, 3))
drop('regular', SEGMENT, UNIT, (1, 1, 0), guess=0.5)
drop('regular', CIRCLE, (0, 4), (2, 0.5, 1), guess=0.2)


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
    header = ('# catallax-trace global point-curve gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
