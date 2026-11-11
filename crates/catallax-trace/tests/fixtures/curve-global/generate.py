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
# Lines and circles SMLib solves with analytic curves: the port refuses.
for test in ((1, 0.5, 0), (0, 1, 0), (3, 0, 0)):
    solve('regular', SEGMENT, UNIT, MINIMIZE, test, op='grefuse')
    solve('regular', SEGMENT, UNIT, INTERSECT, test, op='grefuse')
for test in ((2, 0, 1), (0, 0, 1), (0.5, 0.5, 2)):
    solve('regular', CIRCLE, (0, 4), MINIMIZE, test, op='grefuse')
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
# The analytic line and circle need SMLib's analytic drops without a guess.
drop('regular', SEGMENT, UNIT, (1, 1, 0), op='drefuse')
drop('regular', CIRCLE, (0, 4), (2, 0.5, 1), op='drefuse')
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
