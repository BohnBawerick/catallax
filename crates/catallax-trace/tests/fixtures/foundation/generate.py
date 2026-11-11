#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace foundation gate.

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


def case(op, kind, *tokens):
    assert kind in ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
    CASES.append((op, kind, ' '.join(str(t) for t in tokens)))


rng = random.Random(20261111)


def r(scale=10.0):
    return rng.uniform(-scale, scale)


def r3(scale=10.0):
    return [r(scale), r(scale), r(scale)]


def unit(a):
    n = math.sqrt(sum(x * x for x in a))
    return [x / n for x in a]


# Vectors.
for _ in range(8):
    case('angle', 'regular', v(*r3(), *r3()))
case('angle', 'tangent', v(1, 2, 3, 2, 4, 6))
case('angle', 'tangent', v(1, 2, 3, -1, -2, -3.0000001))
case('angle', 'near-singular', v(1, 0, 0, 1, 1e-9, 0))
case('angle', 'near-singular', v(1, 0, 0, -1, 1e-9, 0))
case('angle', 'boundary', v(1, 0, 0, 0.9, math.sqrt(1 - 0.81), 0))
case('angle', 'regular', v(1, 0, 0, 0, 1, 0))
case('angle', 'singular', v(0, 0, 0, 1, 0, 0))
case('angle', 'singular', v(1e-13, 0, 0, 1, 0, 0))
for a in ([3, 4, 12], [1, 0, 0], [1e-12, 0, 0], [1e-13, 1e-13, 0], [0, 0, 0], r3(1e6)):
    case('unitize', 'singular' if max(map(abs, a)) < 1e-12 else 'regular', v(*a))
for b, deg in (([1, 0.005, 0], 0.5), ([1, 0.01, 0], 0.5), ([-1, 0.005, 0], 0.5), ([0, 0, 0], 0.5), ([1, 0, 0], 0.0)):
    case('parallel', 'tangent' if b != [0, 0, 0] else 'singular', v(1, 0, 0, *b, deg))
for b in ([0.005, 1, 0], [0.01, 1, 0], [0, 0, 0], [1, 1, 0]):
    case('perp', 'touching' if b[0] < 0.1 else 'regular', v(1, 0, 0, *b, 0.5))
case('ortho', 'regular', v(1, 2, 3), 1, v(0, 0, 1))
case('ortho', 'tangent', v(1, 2, 3), 1, v(2, 4, 6))
case('ortho', 'regular', v(1, 2, 3), 0, v(0, 0, 0))
case('ortho', 'boundary', v(0, 5, 0), 0, v(0, 0, 0))
case('ortho', 'singular', v(0, 0, 0), 0, v(0, 0, 0))
case('uderiv', 'regular', v(1, 2, 3, 0.5, -1, 2))
case('uderiv', 'singular', v(1e-13, 0, 0, 0.5, -1, 2))
for a, d in (([1, 2, 3], 7), ([1, 2, 3], 1e-13), ([1, 2, 3], 1e-11), ([0, 0, 0], 0), ([1, 0, 0], 0)):
    case('divassign3', 'singular' if abs(d) < 1e-11 else 'regular', v(*a, d))
for a, d in (([1, 2], 7), ([1, 2], 1e-13), ([0, 2], 1e-13), ([1, 2], 0)):
    case('div2', 'singular' if abs(d) < 1e-11 else 'regular', v(*a, d))
case('closer', 'regular', v(0, 0, 0, 1, 1, 1, 2))
case('closer', 'touching', v(0, 0, 0, 1, 1, 1, math.sqrt(3)))
case('closer', 'regular', v(0, 0, 0, 1, 1, 1, 1.5))
case('differs2', 'regular', v(1, 2, 1, 2))
case('differs2', 'touching', v(1, 2, 1 + 1e-12, 2))
case('differs2', 'regular', v(1, 2, 1.001, 2))

# Extents.
for lo, hi, x, tol, olo, ohi, u, kind in (
    (0, 1, 0.5, 1e-9, 0.25, 2, 0.3, 'regular'),
    (0, 1, 1, 1e-9, 1, 2, 1, 'touching'),
    (0, 1, 1 + 1e-10, 1e-9, 1 + 1e-10, 2, 0, 'touching'),
    (2, 2, 2, 1e-9, 2, 2, 0.5, 'singular'),
    (3, -1, -5, 0.0, -0.5, 0.5, 0.75, 'regular'),
    (-1234567.0, 1234567.0, 7, 1e-9, 0, 1, 0.5, 'boundary'),
    (0, 1234567.0, -1, 0.5, 2, 3, 0.25, 'boundary'),
    (0, 1, 0.4, 0.5, 0.5, 0.6, 1e-17, 'regular'),
):
    case('ext1', kind, v(lo, hi, x, tol, olo, ohi, u))
for lo, hi, p, tol, olo, ohi, uv, kind in (
    ([0, 0], [2, 1], [1, 0.5], 1e-9, [1, 0.5], [3, 3], [0.25, 0.75], 'regular'),
    ([0, 0], [2, 1], [2, 1], 1e-9, [2, 1], [3, 3], [1, 0], 'touching'),
    ([0, 0], [2, 1], [-1, 3], 0.0, [5, 5], [6, 6], [0, 1], 'regular'),
    ([0, 0], [2, 1], [2 + 1e-10, 0.5], 1e-9, [0.5, 0.25], [1, 0.5], [0.5, 0.5], 'touching'),
    ([1, 1], [1, 1], [1, 1], 1e-9, [1, 1], [1, 1], [0.5, 0.5], 'singular'),
    ([2, 1], [0, 0], [0, 0], 1e-12, [2, 1], [2, 3], [1, 1], 'boundary'),
):
    case('ext2', kind, v(*lo, *hi, *p, tol, *olo, *ohi, *uv))
for lo, hi, p, d, t, kind in (
    ([0, 0], [1, 1], [0.5, 0.5], [0.1, 0.2], 1, 'regular'),
    ([0, 0], [1, 1], [0.5, 0.5], [1, 0.2], 1, 'regular'),
    ([0, 0], [1, 1], [0.5, 0.5], [0.2, -1], 1, 'regular'),
    ([0, 0], [1, 1], [0.5, 0.5], [2, 3], -1, 'regular'),
    ([0, 0], [1, 1], [0, 0.5], [-1, 0], 1, 'touching'),
    ([0, 0], [1, 1], [1, 0.5], [1, 0.3], 1, 'touching'),
    ([0, 0], [1, 1], [0.5, 1], [0.1, -1], -1, 'touching'),
    ([0, 0], [1, 1], [-1, 0.5], [0, 1], 1, 'singular'),
    ([0, 0], [1, 1], [0.5, 0.5], [1, 1], 0.5 + 1e-11, 'near-singular'),
):
    case('clip', kind, v(*lo, *hi, *p, *d, t))
for raw, lo, hi, p, tol, olo, ohi, n, kind in (
    (0, [0, 0, 0], [1, 2, 3], [0.5, 1, 1.5], 1e-9, [0.5, 0.5, 0.5], [4, 4, 4], [0.25, 0.5, 0.75], 'regular'),
    (0, [0, 0, 0], [1, 1, 1], [1, 1, 1], 1e-9, [1, 1, 1], [2, 2, 2], [1, 1, 0], 'touching'),
    (0, [0, 0, 0], [1, 1, 1], [3, 0, 0], 0.0, [5, 5, 5], [6, 7, 8], [0, 0, 0], 'regular'),
    (0, [0, 0, 0], [1, 1, 1], [0, 0, 0], 0.0, [10, 0, 0], [11, 1, 1], [0.5, 0.5, 0.5], 'regular'),
    (0, [1, 1, 1], [1, 1, 1 + 1e-12], [1, 1, 1], 1e-9, [1, 1, 1], [1, 1, 1], [0.5, 0.5, 0.5], 'singular'),
    (1, [1, 0, 0], [0, 1, 1], [0.5, 0.5, 0.5], 1e-9, [0, 0, 0], [1, 1, 1], [0.5, 0.5, 0.5], 'boundary'),
    (0, [-1234567.0, 0, 0], [1, 1, 1], [0, 0.5, 0.5], 1e-3, [0, 0, 0], [1, 1, 1], [1, 0.5, 0.5], 'boundary'),
):
    case('ext3', kind, raw, v(*lo, *hi, *p, tol, *olo, *ohi, *n))

# SmGeomUtility.
for _ in range(4):
    case('triarea', 'regular', v(*r3(), *r3(), *r3()))
case('triarea', 'singular', v(0, 0, 0, 1, 1, 1, 2, 2, 2))
case('triarea', 'singular', v(1, 2, 3, 1, 2, 3, 1, 2, 3))
tri = [0, 0, 0, 2, 0, 0, 0, 2, 0]
for t, kind in (([0.5, 0.5, 1], 'regular'), ([-1, -1, 0], 'regular'), ([3, 3, 0], 'regular'),
                ([1, 1, 0], 'touching'), ([0, 0, 0], 'touching'), ([3, -1, 2], 'regular'),
                ([-1, 3, 0], 'regular')):
    case('tripoint', kind, v(*tri, *t))
case('tripoint', 'singular', v(0, 0, 0, 1, 1, 0, 2, 2, 0, 3, 0, 0))
case('tripoint', 'singular', v(0, 0, 0, 0, 0, 0, 1, 0, 0, 0.5, 1, 0))
for _ in range(4):
    case('tripoint', 'regular', v(*r3(), *r3(), *r3(), *r3()))
for p, d, t, kind in (([1, 2, 3], [1, 1, 0], [4, 0, 1], 'regular'), ([1, 2, 3], [0, 0, 0], [4, 0, 1], 'singular'),
                      ([1, 2, 3], [1e-13, 0, 0], [4, 0, 1], 'singular'), ([0, 0, 0], [1, 0, 0], [5, 0, 0], 'touching')):
    case('linepoint', kind, v(*p, *d, *t))
    case('lineclosest', kind, v(*p, *d, *t))
for s, e, t, kind in (([0, 0, 0], [2, 0, 0], [1, 1, 0], 'regular'), ([0, 0, 0], [2, 0, 0], [-1, 1, 0], 'regular'),
                      ([0, 0, 0], [2, 0, 0], [3, -1, 1], 'regular'), ([0, 0, 0], [2, 0, 0], [2, 0, 0], 'touching'),
                      ([0, 0, 0], [2, 0, 0], [5, 0, 0], 'tangent'), ([1, 1, 1], [1, 1, 1], [2, 2, 2], 'singular')):
    case('segpoint', kind, v(*s, *e, *t))
    case('segangle', kind, v(*s, *e, *t))
circle = [0, 0, 0, 0, 0, 1]
for p, d, radius, tol, kind in (
    ([-3, 0.5, 0], [1, 0, 0], 2, 1e-6, 'regular'),
    ([-3, 2, 0], [1, 0, 0], 2, 1e-6, 'tangent'),
    ([-3, 2 + 1e-7, 0], [1, 0, 0], 2, 1e-6, 'tangent'),
    ([-3, 2 - 1e-7, 0], [1, 0, 0], 2, 1e-6, 'tangent'),
    ([-3, 2 - 1e-13, 0], [2, 0, 0], 2, 1e-6, 'tangent'),
    ([-3, 2 - 4e-7, 0], [1, 0, 0], 2, 1e-6, 'near-singular'),
    ([-3, 2 + 1e-5, 0], [1, 0, 0], 2, 1e-6, 'regular'),
    ([-3, 0.5, 5], [1, 0.1, 3], 2, 1e-6, 'regular'),
    ([1, 1, 0], [0, 0, 1], 2, 1e-6, 'singular'),
    ([-3, 1.5, 0], [1, 0, 0], 2, 0.5, 'boundary'),
):
    case('linecircle', kind, v(*p, *d, *circle, radius, tol))
n2 = unit([1, 1, 0])
for p1, nn1, p2, nn2, kind in (
    ([0, 0, 0], [0, 0, 1], [1, 0, 0], [1, 0, 0], 'regular'),
    ([1, 2, 3], unit([1, 2, 2]), [-1, 0, 4], n2, 'regular'),
    ([0, 0, 0], [0, 0, 1], [5, 5, 0], [0, 0, 1], 'touching'),
    ([0, 0, 0], [0, 0, 1], [5, 5, 1], [0, 0, 1], 'singular'),
    ([0, 0, 0], [0, 0, 1], [5, 5, 0], [0, 0, -1], 'touching'),
    ([0, 0, 0], [0, 0, 1], [1, 0, 0], unit([1e-13, 0, 1]), 'near-singular'),
    ([0, 0, 0], [0, 0, 2], [1, 0, 0], [1, 0, 0], 'singular'),
):
    case('twoplanes', kind, v(*p1, *nn1, *p2, *nn2))
case('threeplanes', 'regular', v(1, 2, 3, 1, 0, 0, 4, 5, 6, 0, 1, 0, 7, 8, 9, 0, 0, 1))
case('threeplanes', 'regular', v(*r3(), *unit([1, 2, 3]), *r3(), *unit([-2, 1, 0.5]), *r3(), *unit([0.3, -0.7, 1])))
case('threeplanes', 'singular', v(1, 2, 3, 1, 0, 0, 4, 5, 6, 1, 0, 0, 7, 8, 9, 0, 0, 1))
case('threeplanes', 'singular', v(1, 2, 3, 2, 0, 0, 4, 5, 6, 0, 1, 0, 7, 8, 9, 0, 0, 1))
case('threeplanes', 'near-singular', v(0, 0, 0, 1, 0, 0, 0, 0, 0, *unit([1, 1e-13, 0]), 0, 0, 0, 0, 0, 1))
case('planedist', 'regular', v(1, 2, 3, *unit([1, 2, 2]), 5, -1, 0))
case('planedist', 'singular', v(1, 2, 3, 0, 0, 1.1, 5, -1, 0))
case('planedist', 'boundary', v(1, 2, 3, 0, 0, 1 + 1e-14, 5, -1, 3))
case('ptproj', 'regular', v(5, -1, 0, 1, 2, 3, 1, 2, 2))
case('ptproj', 'touching', v(1, 2, 3, 1, 2, 3, 0, 0, 1))
case('vecproj', 'regular', v(5, -1, 0, *unit([1, 2, 2])))
case('vecproj', 'tangent', v(0, 0, 3, 0, 0, 1))
case('lincomb', 'regular', v(1, 0, 0, 1, 1, 0, 3, 2, 5))
case('lincomb', 'regular', v(*r3(), *r3(), *r3()))
case('lincomb', 'singular', v(1, 2, 3, 2, 4, 6, 1, 0, 0))
case('lincomb', 'singular', v(0, 0, 0, 1, 0, 0, 1, 0, 0))
case('lincomb', 'near-singular', v(1, 0, 0, 1, 1e-9, 0, 0, 1, 0))
case('det3', 'regular', v(*r3(), *r3(), *r3()))
case('det3', 'singular', v(1, 2, 3, 2, 4, 6, 7, 8, 9))
helix = [[1, 0, 0], [0, 1, 0.5], [-1, 0, 0], [0, -1, 0]]
for count in range(4):
    case('frenet', 'regular', count, v(*sum(helix, [])))
case('frenet', 'singular', 3, v(1, 2, 3, 0, 0, 0, 1, 0, 0, 0, 1, 0))
case('frenet', 'singular', 3, v(1, 2, 3, 1, 0, 0, 0, 0, 0, 0, 1, 0))
case('frenet', 'singular', 3, v(1, 2, 3, 1, 0, 0, 2, 0, 0, 0, 1, 0))
case('frenet', 'regular', 3, v(*r3(), *r3(), *r3(), *r3()))
p, d1, s1, t1 = [1, 2, 3], [1, 0.5, 0], [0, 1, 0.25], [0.1, 0, 1]
for jet2, deg, has, kind in (
    ((p, d1, s1, t1), 0.5, 1, 'regular'),
    ((p, d1, s1, t1), 0.5, 0, 'regular'),
    ((p, [2, 1, 0], [0, 4, 1], [0.8, 0, 8]), 0.5, 1, 'regular'),
    ((p, [2, 1, 0], [0, 1, 0.25], t1), 0.5, 1, 'regular'),
    ((p, [1, 0.5, 0.009], s1, t1), 0.5, 1, 'touching'),
    ((p, [1, -0.5, 0], s1, t1), 0.5, 1, 'regular'),
    (([1, 2, 3.001], d1, s1, t1), 0.5, 1, 'regular'),
    ((p, [0, 0, 0], s1, t1), 0.5, 1, 'singular'),
    ((p, d1, [0, 1, 0.25 + 1e-6], t1), 0.5, 1, 'touching'),
    ((p, d1, s1, [0.1, 0, 1 + 1e-3]), 0.5, 1, 'touching'),
    ((p, [1 + 1e-6, 0.5, 0], s1, t1), 0.5, 1, 'touching'),
    ((p, [1 + 3e-6, 0.5, 0], s1, t1), 0.5, 1, 'touching'),
    ((p, d1, [1e-6, 1, 0.25], t1), 0.5, 1, 'touching'),
    ((p, d1, [3e-6, 1, 0.25], t1), 0.5, 1, 'touching'),
    ((p, d1, s1, [0.1 + 1e-6, 0, 1]), 0.5, 1, 'touching'),
    ((p, d1, s1, [0.1 + 3e-6, 0, 1]), 0.5, 1, 'touching'),
):
    case('continuity', kind, v(*p, *d1, *s1, *jet2[0], *jet2[1], *jet2[2], deg), has, v(*t1, *jet2[3]))
case('continuity', 'singular', v(*p, 0, 0, 0, *s1, *p, *d1, *s1, 0.5), 0, v(0, 0, 0, 0, 0, 0))


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
    header = ('# catallax-trace foundation gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
