#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace curve-fit gate.

usage: generate.py ORACLE [RAW_DIR]

ORACLE is the driver built from oracle.cpp against the pinned SMLib. Each
case runs in a fresh process under a memory and CPU cap and a timeout. RAW_DIR,
when given, receives every request, response and error stream; keep it outside
the repository. Inputs are fixed: hand-picked tangent, touching, singular and
boundary cases plus seeded random regular cases.
"""
import math
import pathlib
import random
import struct
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent


def h(x):
    return format(struct.unpack('<Q', struct.pack('<d', float(x)))[0], '016x')


def v(*xs):
    return ' '.join(h(x) for x in xs)


CASES = []


def case(kind, par, dim, degree, points, vectors, ignore, more=None):
    assert kind in ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
    tokens = [par, dim, degree, len(points)] + [v(*p) for p in points] + [len(vectors)] + [v(*x) for x in vectors] + [ignore] + [0 if more is None else 1, 0 if more is None else len(more)] + [v(*x) for x in (more or [])]
    CASES.append(('interpolate', kind, ' '.join(str(t) for t in tokens)))


rng = random.Random(20261114)


def r3(scale=4.0):
    return [rng.uniform(-scale, scale) for _ in range(3)]


def path(count, step=1.0):
    pts, cur = [], [0.0, 0.0, 0.0]
    for _ in range(count):
        pts.append(list(cur))
        cur = [c + rng.uniform(.2, 1.5) * step * s for c, s in zip(cur, (1, rng.choice((-1, 1)), rng.choice((-1, 1))))]
    return pts


# Points with a vector each: every parameterisation, dimension, degree and size.
for par in (0, 1, 2):
    for dim in (2, 3):
        for degree in (2, 3):
            for count in (3, 4, 5, 8):
                pts = path(count)
                case('regular', par, dim, degree, pts, [r3(2.0) for _ in pts], 0)
                case('regular', par, dim, degree, pts, [r3(2.0) for _ in pts], 1)
# Vectors along the path, and straight data.
for degree in (2, 3):
    pts = [[float(i), 0.0, 0.0] for i in range(5)]
    case('tangent', 1, 3, degree, pts, [[1.0, 0.0, 0.0]] * 5, 0)
    case('tangent', 0, 3, degree, pts, [[2.0, 0.0, 0.0]] * 5, 1)
    case('tangent', 2, 3, degree, pts, [[1.0, 0.0, 0.0]] * 5, 0)
    pts = path(6)
    vec = [[b - a for a, b in zip(pts[max(i - 1, 0)], pts[min(i + 1, 5)])] for i in range(6)]
    case('tangent', 1, 3, degree, pts, vec, 0)
    case('tangent', 1, 3, degree, pts, vec, 1)
# Data scales and a long data set.
for scale in (1e-6, 1e-3, 1e3, 1e8):
    pts = [[c * scale for c in p] for p in path(5)]
    case('regular', 1, 3, 3, pts, [[c * scale for c in r3()] for _ in pts], 0)
case('regular', 2, 3, 3, path(30), [r3() for _ in range(30)], 0)
case('regular', 0, 3, 2, path(30), [r3() for _ in range(30)], 1)
# Repeated and nearly repeated points: chord lengths vanish.
base = path(5)
for par in (0, 1, 2):
    dup = [base[0], base[1], base[1], base[2], base[3]]
    case('touching', par, 3, 3, dup, [r3() for _ in dup], 0)
    dup = [base[0], base[0], base[2], base[3], base[4]]
    case('touching', par, 3, 2, dup, [r3() for _ in dup], 0)
    near = [base[0], [base[1][0] + 1e-13, base[1][1], base[1][2]], base[2], base[3], base[4]]
    case('near-singular', par, 3, 3, near, [r3() for _ in near], 0)
    case('singular', par, 3, 3, [base[1]] * 4, [r3() for _ in range(4)], 0)
    case('singular', par, 3, 2, [[0.0, 0.0, 0.0]] * 3, [r3() for _ in range(3)], 0)
case('near-singular', 1, 3, 3, [[0, 0, 0], [1e-7, 0, 0], [2e-7, 0, 0], [3e-7, 0, 0]], [r3() for _ in range(4)], 0)
case('near-singular', 1, 3, 3, [[0, 0, 0], [3e-13, 0, 0], [6e-13, 0, 0], [9e-13, 0, 0]], [r3() for _ in range(4)], 0)
# Zero vectors cannot be normalised.
pts = path(4)
case('singular', 1, 3, 3, pts, [[0.0, 0.0, 0.0]] + [r3() for _ in range(3)], 1)
case('singular', 1, 3, 3, pts, [r3(), r3(), [1e-13, 0, 0], r3()], 1)
case('boundary', 1, 3, 3, pts, [[0.0, 0.0, 0.0]] + [r3() for _ in range(3)], 0)
case('near-singular', 1, 3, 3, pts, [r3(), [1e-300, 0, 0], r3(), r3()], 1)
# Refused degrees and dimensions.
for degree in (0, 1, 4, 5, 17):
    case('singular', 1, 3, degree, pts, [r3() for _ in pts], 0)
for dim in (0, 1, 4):
    case('singular', 1, dim, 3, pts, [r3() for _ in pts], 0)
# The line shortcut.
for dim in (2, 3, 1, 4):
    case('boundary', 1, dim, 1, pts[:2], [], 0)
    case('boundary', 1, dim, 1, pts[:2], [r3(), r3()], 0)
case('singular', 1, 3, 1, [pts[0], pts[0]], [], 0)
case('boundary', 1, 3, 1, pts[:3], [r3() for _ in range(3)], 0)
# Derivatives at the ends only, higher derivatives at the ends, and none.
for count in (0, 1, 2):
    for degree in (1, 2, 3, 4, 5, 8):
        for par in (0, 1, 2):
            for npts in (2, 3, 5, 8):
                for nmore in ((None,) if count == 0 else (None, 0, 1, 2, 3, 4, 5)):
                    more = None if nmore is None else [r3() for _ in range(nmore)]
                    kind = 'regular'
                    if degree == 1 or npts <= 2 and degree > 2:
                        kind = 'boundary'
                    case(kind, par, 3, degree, path(npts), [r3() for _ in range(count)], 0, more)
        case('regular', 1, 2, degree, path(5), [r3() for _ in range(count)], 0, [r3(), r3()] if count else None)
case('regular', 1, 3, 3, path(2), [r3(), r3()], 0)
case('regular', 1, 3, 3, path(2), [], 0)
case('regular', 1, 3, 3, path(4), [r3() for _ in range(3)], 0)
case('singular', 1, 3, 3, path(4), [r3() for _ in range(5)], 0)
case('boundary', 1, 3, 3, [], [], 0)
case('boundary', 1, 3, 3, [[1.0, 2.0, 3.0]], [], 0)


# Fairing a blend's end tangents.
def fair(kind, p0, p1, t0, t1, higher, prevent):
    assert kind in ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
    CASES.append(('fair', kind, ' '.join([v(*p0), v(*p1), v(*t0), v(*t1), str(len(higher))] + [v(*x) for x in higher] + [str(prevent)])))


for prevent in (0, 1):
    for _ in range(12):
        fair('regular', r3(), r3(), r3(2.0), r3(2.0), [], prevent)
    for nhigher in (2, 4, 1, 3):
        fair('regular', r3(), r3(), r3(2.0), r3(2.0), [r3() for _ in range(nhigher)], prevent)
    # Parallel and opposite tangents, and tangents along the chord.
    fair('tangent', [0, 0, 0], [4, 0, 0], [1, 0, 0], [1, 0, 0], [], prevent)
    fair('tangent', [0, 0, 0], [4, 0, 0], [1, 0, 0], [-1, 0, 0], [], prevent)
    fair('tangent', [0, 0, 0], [4, 0, 0], [1, 1, 0], [1, -1, 0], [], prevent)
    fair('tangent', [0, 0, 0], [4, 1, 0], [3, 0.5, 0], [3, 0.5, 0], [], prevent)
    # An S-shaped blend takes the inflection branch.
    fair('regular', [0, 0, 0], [4, 0, 0], [1, 2, 0], [1, 2, 0], [], prevent)
    fair('regular', [0, 0, 0], [4, 0.5, 1], [0.2, 1, 0], [1, -3, 2], [], prevent)
    fair('regular', [0, 0, 0], [4, 1, 0], [1, 3, 0], [3, 1, 0], [], prevent)
    fair('boundary', [0, 0, 0], [4, 0, 0], [0.01, 3, 0], [0.01, 3, 0], [], prevent)
    # Coincident points, zero and tiny tangents.
    fair('touching', [1, 2, 3], [1, 2, 3], [1, 0, 0], [0, 1, 0], [], prevent)
    fair('singular', [0, 0, 0], [4, 0, 0], [0, 0, 0], [1, 0, 0], [], prevent)
    fair('singular', [0, 0, 0], [4, 0, 0], [1, 0, 0], [0, 0, 0], [], prevent)
    fair('near-singular', [0, 0, 0], [4, 0, 0], [1e-6, 1e-6, 0], [1, 0, 1], [], prevent)
    fair('near-singular', [0, 0, 0], [4, 0, 0], [1e-13, 0, 0], [1, 0, 1], [], prevent)
    fair('near-singular', [0, 0, 0], [1e-9, 0, 0], [1, 1, 0], [1, 0, 1], [], prevent)


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
                f'request: {request}exit: {run.returncode}\nstdout: {run.stderr and run.stdout}stderr: {run.stderr}')
        if run.returncode != 0:
            sys.exit(f'case {index} {op} {kind}: oracle exit {run.returncode}: {run.stderr.strip()}')
        response = run.stdout.split()
        assert response[0] == op, (index, response)
        if op == 'fair':
            lines.append(f'fair {kind} | {inputs} | {" ".join(response[1:])}')
            continue
        tokens = inputs.split()
        dim, degree, count = int(tokens[1]), int(tokens[2]), int(tokens[3])
        vector_count = int(tokens[4 + 3 * count])
        success = response[1] == 'i1000'
        # Results the checked curve cannot reproduce are pinned apart.
        pinned = op
        if success and dim == 2 and not (degree == 1 and count == 2):
            pinned = 'interpolate-planar'
        lines.append(f'{pinned} {kind} | {inputs} | {" ".join(response[1:])}')
    header = ('# catallax-trace curve-fit gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
