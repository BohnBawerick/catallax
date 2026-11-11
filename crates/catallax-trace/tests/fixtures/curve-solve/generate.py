#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace curve solver gate.

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
MINIMIZE, MAXIMIZE, NORMALIZE, INTERSECT = 3, 4, 5, 0


def h(x):
    return format(struct.unpack('<Q', struct.pack('<d', float(x)))[0], '016x')


def v(*xs):
    return ' '.join(h(x) for x in xs)


def curve(degree, knots, points, rational=False):
    tokens = [h(k) for k in knots]
    for p in points:
        tokens += [h(p[0]), h(p[1]), h(p[2]), h(p[3] if len(p) > 3 else 1.0)]
    return f'{degree} {len(points)} {int(rational)} ' + ' '.join(tokens)


def case(op, kind, *parts):
    assert kind in KINDS
    CASES.append((op, kind, ' '.join(str(p) for p in parts)))


bez = lambda pts: curve(3, [0, 0, 0, 0, 1, 1, 1, 1], pts)
WIGGLE = bez([(0.5, 0.4, -2), (1.5, 2, 3), (3, 3, -3), (4.5, 4.6, 2)])
ARCH = bez([(0, 0, 0), (1, 2, 0), (2, 2, 0), (3, 0, 0)])
SPLINE = curve(3, [0, 0, 0, 0, 1, 2, 3, 3, 3, 3], [(0, 0, 0), (1, 2, 0), (2, -1, 1), (3, 2, 0), (4, 0, 2), (5, 1, 0)])
LINE = curve(1, [0, 0, 1, 1], [(0, 0, 0), (4, 0, 0)])
LINE2 = curve(1, [0, 0, 1, 1], [(2, 0, 0), (6, 0, 0)])
CROSS = curve(1, [0, 0, 1, 1], [(1, -1, 0), (1, 1, 0)])
r = math.sqrt(0.5)
CIRCLE = curve(2, [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4],
               [(1, 0, 0, 1), (1, 1, 0, r), (0, 1, 0, 1), (-1, 1, 0, r), (-1, 0, 0, 1), (-1, -1, 0, r),
                (0, -1, 0, 1), (1, -1, 0, r), (1, 0, 0, 1)], rational=True)
PARABOLA = curve(2, [0, 0, 0, 1, 1, 1], [(-1, 1, 0), (0, -1, 0), (1, 1, 0)])
BIG = bez([(0, 0, 0), (400, 900, 0), (800, -500, 300), (1200, 100, 0)])
COLLAPSED = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0, 0, 0), (0, 0, 0), (2, 1, 0), (3, 3, 1)])
SQUARE = curve(1, [0, 0, 1, 2, 3, 4, 4], [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 0)])
UNIT = v(0, 1)

# Point drops.
for test, guess, kind in (((2, 3, 0), 0.5, 'regular'), ((5, 5, 0), 0.9, 'regular'), ((-1, 0, -3), 0.1, 'boundary'),
                          ((3.6, 3.9, -1.05), 0.6, 'touching')):
    case('point', kind, WIGGLE, UNIT, MINIMIZE, v(*test), h(guess))
case('point', 'regular', WIGGLE, UNIT, MAXIMIZE, v(2, 3, 0), h(0.5))
case('point', 'boundary', WIGGLE, v(0.2, 0.6), NORMALIZE, v(2, 3, 0), h(0.3))
case('point', 'regular', WIGGLE, UNIT, INTERSECT, v(2, 3, 0), h(0.5))
case('point', 'touching', ARCH, UNIT, MINIMIZE, v(1.5, 1.5, 0), h(0.3))
case('point', 'singular', CIRCLE, v(0, 4), MINIMIZE, v(0, 0, 0), h(1.3))
case('point', 'regular', CIRCLE, v(0, 4), MINIMIZE, v(3, 2, 1), h(0.2))
case('point', 'regular', SPLINE, v(0, 3), MINIMIZE, v(2.5, 3, 0), h(1.5))
# Points far inside a tight bend, where SMLib limits the second-derivative term.
for test, guess in (((-0.5, 0, 0), 0.3), ((-0.9, 0.1, 0), 0.05), ((0.2, -0.95, 0), 1.2)):
    case('point', 'regular', CIRCLE, v(0, 4), MINIMIZE, v(*test), h(guess))
case('point', 'regular', ARCH, UNIT, MINIMIZE, v(1.5, -3, 0), h(0.1))
case('point', 'regular', LINE, UNIT, MINIMIZE, v(1, 2, 0), h(0.9))
case('point', 'boundary', LINE, UNIT, MINIMIZE, v(7, 2, 0), h(0.5))
case('point', 'regular', BIG, UNIT, MINIMIZE, v(500, 200, 100), h(0.5))
case('point', 'singular', COLLAPSED, UNIT, MINIMIZE, v(-1, 0, 0), h(0.0))
case('point', 'regular', COLLAPSED, UNIT, MINIMIZE, v(2, 2, 0), h(0.5))

# Two curves: crossing, tangent, near-tangent, disjoint, end contact and
# coincident cases.
case('intersect', 'regular', ARCH, CROSS, UNIT, UNIT, h(1e-6), h(0.4), h(0.6))
case('intersect', 'regular', WIGGLE, SPLINE, UNIT, v(0, 3), h(1e-6), h(0.5), h(1.5))
case('intersect', 'tangent', PARABOLA, LINE, UNIT, UNIT, h(1e-6), h(0.45), h(0.2))
for offset, kind in ((-1e-8, 'touching'), (1e-8, 'near-singular'), (-1e-3, 'regular')):
    line = curve(1, [0, 0, 1, 1], [(-2, offset, 0), (2, offset, 0)])
    case('intersect', kind, PARABOLA, line, UNIT, UNIT, h(1e-6), h(0.45), h(0.55))
case('intersect', 'boundary', LINE, CROSS, UNIT, UNIT, h(1e-6), h(0.5), h(0.5))
case('intersect', 'boundary', ARCH, LINE, UNIT, UNIT, h(1e-6), h(0.05), h(0.1))
case('intersect', 'singular', LINE, LINE2, UNIT, UNIT, h(1e-6), h(0.7), h(0.2))
case('intersect', 'regular', ARCH, curve(1, [0, 0, 1, 1], [(0, 5, 0), (3, 5, 0)]), UNIT, UNIT, h(1e-6), h(0.5), h(0.5))
case('intersect', 'regular', CIRCLE, CROSS, v(0, 4), UNIT, h(1e-6), h(0.2), h(0.5))
# Starts far from the meeting, where SMLib limits the curvature terms of the
# Jacobian.
case('intersect', 'regular', ARCH, CROSS, UNIT, UNIT, h(1e-6), h(0.0), h(1.0))
case('intersect', 'regular', CIRCLE, CROSS, v(0, 4), UNIT, h(1e-6), h(1.9), h(0.0))
case('intersect', 'regular', CIRCLE, curve(1, [0, 0, 1, 1], [(-2, 0.5, 0), (2, 0.5, 0)]), v(0, 4), UNIT,
     h(1e-6), h(3.0), h(0.9))
case('intersect', 'regular', WIGGLE, SPLINE, UNIT, v(0, 3), h(1e-6), h(0.0), h(3.0))

# Closest and farthest pairs.
case('solve', 'regular', WIGGLE, SPLINE, UNIT, v(0, 3), MINIMIZE, h(1e-6), h(0.3), h(1.0))
case('solve', 'regular', ARCH, LINE, UNIT, UNIT, MAXIMIZE, h(1e-6), h(0.4), h(0.4))
case('solve', 'singular', LINE, LINE2, UNIT, UNIT, MINIMIZE, h(1e-6), h(0.7), h(0.2))
case('solve', 'tangent', PARABOLA, LINE, UNIT, UNIT, NORMALIZE, h(1e-6), h(0.4), h(0.3))

# Largest gaps between corresponding pieces.
shifted = bez([(0.5, 0.4, -2 + 1e-6), (1.5, 2, 3 + 1e-6), (3, 3, -3 + 1e-6), (4.5, 4.6, 2 + 1e-6)])
case('maxdist', 'touching', WIGGLE, shifted, UNIT, h(0), h(1), 5, 1, h(1e-3))
case('maxdist', 'touching', WIGGLE, shifted, UNIT, h(0), h(1), 0, 0, h(0))
case('maxdist', 'regular', ARCH, LINE, UNIT, h(0), h(0.75), 7, 0, h(0))
case('maxdist', 'regular', ARCH, LINE, UNIT, h(0.75), h(0), 7, 1, h(0.5))
case('maxdist', 'singular', LINE, LINE2, UNIT, h(0.5), h(0.5), 2, 0, h(0))

# Coincidence: overlapping lines, end contact, a curve against its copy on
# overlapping intervals, a near copy and a far copy, with and without boxes.
for boxes in (0, 1):
    case('coincide', 'singular', LINE, LINE2, UNIT, UNIT, h(1e-6), boxes)
    case('coincide', 'touching', LINE, curve(1, [0, 0, 1, 1], [(4, 0, 0), (5, 3, 0)]), UNIT, UNIT, h(1e-6), boxes)
    case('coincide', 'singular', WIGGLE, WIGGLE, v(0, 0.6), v(0.4, 1), h(1e-6), boxes)
    case('coincide', 'near-singular', WIGGLE, shifted, v(0.1, 0.9), v(0.2, 0.8), h(1e-5), boxes)
    case('coincide', 'regular', WIGGLE, shifted, v(0.1, 0.9), v(0.2, 0.8), h(1e-7), boxes)
    case('coincide', 'regular', ARCH, CROSS, UNIT, UNIT, h(1e-6), boxes)

# Closed and periodic tests.
case('closed', 'regular', CIRCLE, v(0, 4), h(0))
case('closed', 'boundary', CIRCLE, v(0, 2), h(0))
case('closed', 'touching', SQUARE, v(0, 4), h(0))
case('closed', 'regular', WIGGLE, UNIT, h(1e-3))


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
    header = ('# catallax-trace curve solver gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
