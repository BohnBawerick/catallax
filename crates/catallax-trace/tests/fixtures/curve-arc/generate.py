#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace circle-arc and end-point gate.

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
NO_Z = sys.float_info.max
NO_W = -sys.float_info.max
NAN = float('nan')


def h(x):
    return format(struct.unpack('<Q', struct.pack('<d', float(x)))[0], '016x')


def v(*xs):
    return ' '.join(h(x) for x in xs)


CASES = []


def case(op, kind, *tokens):
    assert kind in ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
    CASES.append((op, kind, ' '.join(str(t) for t in tokens)))


rng = random.Random(20261113)


def r3(scale=4.0):
    return [rng.uniform(-scale, scale) for _ in range(3)]


def unit(a):
    n = math.sqrt(sum(x * x for x in a))
    return [x / n for x in a]


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def random_frame():
    x = unit(r3())
    t = unit(r3())
    z = unit(cross(x, t))
    y = cross(z, x)
    return r3(), x, y


FRAMES = [([0, 0, 0], [1, 0, 0], [0, 1, 0]), ([1, 2, 3], [0, 1, 0], [0, 0, 1]),
          ([-1, 0.5, 2], [1, 0, 0], [0, 0, 1]), random_frame(), random_frame(),
          ([0, 0, 0], [3, 0, 0], [0, 5, 0]),            # not unit
          ([0, 0, 0], [1, 0, 0], [0.001, 1, 0]),        # not perpendicular
          ([0, 0, 0], [1, 0, 0], [1e-13, 1, 0])]        # perpendicular within tolerance


def circle(kind, dim, frame, radius, start, end, quintic=0):
    op = 'circle' if (not quintic or (end - start) % 360.0 == 0 or abs(end - start) >= 360.0 - 0.01) else 'circle-quintic-partial'
    case(op if quintic else 'circle', kind, dim, v(*frame[0], *frame[1], *frame[2]), v(radius, start, end), quintic)


# Quadratic arcs: every span count and its thresholds.
sweeps = [10, 45, 90, 90.9999, 91, 91.0001, 135, 180, 180.9999, 181, 225, 270, 270.9999, 271, 300, 359, 359.98,
          359.9999, 360]
for frame in FRAMES:
    for sweep in (30, 90, 91, 200, 271, 360):
        for dim in (2, 3):
            circle('regular', dim, frame, rng.uniform(.5, 5), 0.0, sweep)
for sweep in sweeps:
    kind = 'boundary' if sweep in (91, 181, 271, 360, 90.9999, 180.9999, 270.9999, 359.98, 359.9999) else 'regular'
    circle(kind, 3, FRAMES[3], 2.0, 0.0, sweep)
    circle(kind, 3, FRAMES[3], 2.0, 37.5, 37.5 + sweep if 37.5 + sweep <= 360 else 360)
# Start and end angles: negative, wrapped, equal, at the limits.
for start, end in ((-90, 90), (-360, 0), (0, -90), (350, 10), (200, 100), (-300, 300), (-360, 360), (360, 360),
                   (90, 90), (0, 0), (-360, -360), (45, 45 + 1e-7), (45, 45 + 1e-6), (45, 45 + 2e-6),
                   (0, 359.99), (0, 359.995), (0, 360 - 1e-6), (0, 90 + 1e-6)):
    kind = 'touching' if start == end else 'boundary' if abs(start) in (360,) or abs(end) in (360,) else 'near-singular' if abs(end - start) < 1e-5 else 'regular'
    circle(kind, 3, FRAMES[3], 1.5, float(start), float(end))
# Refused inputs.
for radius in (0.0, -1.0, 1e-13, 0.9e-12, NAN):
    circle('singular', 3, FRAMES[0], radius, 0.0, 90.0)
circle('boundary', 3, FRAMES[0], 1e-12, 0.0, 90.0)
circle('near-singular', 3, FRAMES[0], 1.0000001e-12, 0.0, 90.0)
for dim in (1, 4, 0):
    circle('singular', dim, FRAMES[0], 1.0, 0.0, 90.0)
for start, end in ((-360.0001, 0.0), (0.0, 360.0001), (361.0, 0.0), (0.0, -361.0), (NAN, 0.0), (0.0, NAN)):
    circle('singular', 3, FRAMES[0], 1.0, start, end)
for radius in (1e-6, 1e6, 1e12):
    circle('regular', 3, FRAMES[3], radius, 0.0, 120.0)
# Quintic whole circles; partial arcs are not carried and are pinned apart.
for frame in FRAMES[:5]:
    for dim in (2, 3):
        circle('regular', dim, frame, rng.uniform(.5, 5), 0.0, 360.0, quintic=1)
for start, end in ((-360, 0), (90, 90), (45, 405 - 360 + 360), (0, 359.995), (-180, 180)):
    if end - start >= 360 - 0.01 or start == end:
        circle('boundary', 3, FRAMES[3], 2.0, float(start), float(end), quintic=1)
for radius in (0.0, 1e-13):
    circle('singular', 3, FRAMES[0], radius, 0.0, 360.0, quintic=1)
for start, end in ((0, 90), (0, 180), (-90, 90), (30, 300), (0, 359), (10, 20)):
    case('circle-quintic-partial', 'regular', 3, v(*FRAMES[3][0], *FRAMES[3][1], *FRAMES[3][2]), v(2.0, start, end), 1)
case('circle-quintic-partial', 'boundary', 3, v(*FRAMES[3][0], *FRAMES[3][1], *FRAMES[3][2]), v(2.0, 0.0, 359.98), 1)


def clamped(degree, interior):
    return [0.0] * (degree + 1) + list(interior) + [3.0] * (degree + 1)


def curve(dim, degree, interior, rational):
    knots = clamped(degree, interior)
    count = len(knots) - degree - 1
    vals = []
    for _ in range(count):
        x, y, z = r3()
        w = rng.uniform(.5, 2) if rational else None
        vals += [x * w, y * w, z * w if dim == 3 else NO_Z, w] if rational else [x, y, z if dim == 3 else NO_Z, NO_W]
    return f'{dim} {degree} {count} ' + v(*knots) + ' ' + v(*vals)


# End points.
for degree in (1, 2, 3):
    for rational in (False, True):
        for dim in (2, 3):
            for start in (1, 0):
                c = curve(dim, degree, [1.0] if degree > 1 else [], rational)
                p = r3()
                case('edit-end', 'regular', c, v(*p), start, 0, v(0.0))
                if rational:
                    case('edit-end', 'regular', c, v(*p), start, 1, v(rng.uniform(.5, 2)))
                    case('edit-end', 'boundary', c, v(*p), start, 1, v(1.0))
                else:
                    case('edit-end', 'singular', c, v(*p), start, 1, v(1.5))
                case('edit-end', 'touching', c, v(0.0, 0.0, 0.0), start, 0, v(0.0))
# Inputs the checked curve cannot hold.
for start in (1, 0):
    case('edit-end-invalid', 'singular', curve(3, 2, [1.0], True), v(1, 2, 3), start, 1, v(0.0))
    case('edit-end-invalid', 'singular', curve(3, 2, [1.0], False), v(NAN, 2, 3), start, 0, v(0.0))
    case('edit-end-invalid', 'singular', curve(3, 2, [1.0], False), v(1, float('inf'), 3), start, 0, v(0.0))


def main():
    oracle = sys.argv[1]
    raw = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if raw:
        raw.mkdir(parents=True, exist_ok=True)
    lines = []
    for index, (op, kind, inputs) in enumerate(CASES):
        request_op = {'circle-quintic-partial': 'circle', 'edit-end-invalid': 'edit-end'}.get(op, op)
        request = f'{request_op} {inputs}\n'
        run = subprocess.run(['prlimit', '--as=1073741824', '--cpu=10', 'timeout', '20', oracle],
                             input=request, capture_output=True, text=True)
        if raw:
            (raw / f'{index:04d}-{op}.txt').write_text(
                f'request: {request}exit: {run.returncode}\nstdout: {run.stdout}stderr: {run.stderr}')
        if run.returncode != 0:
            sys.exit(f'case {index} {op} {kind}: oracle exit {run.returncode}: {run.stderr.strip()}')
        response = run.stdout.split()
        assert response[0] == request_op, (index, response)
        lines.append(f'{op} {kind} | {inputs} | {" ".join(response[1:])}')
    header = ('# catallax-trace circle-arc and end-point gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
