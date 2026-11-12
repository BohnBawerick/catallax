#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace offset-surface gate.

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
NO_W = -sys.float_info.max


def h(x):
    return format(struct.unpack('<Q', struct.pack('<d', float(x)))[0], '016x')


def v(*xs):
    return ' '.join(h(x) for x in xs)


CASES = []


def case(op, kind, *tokens):
    assert kind in ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
    CASES.append((op, kind, ' '.join(str(t) for t in tokens)))


rng = random.Random(20261116)


def knots(degree, count, hi):
    interior = count - degree - 1
    return [0.0] * (degree + 1) + [hi * (i + 1) / (interior + 1) for i in range(interior)] + [hi] * (degree + 1)


def surface(du, dv, nu, nv, rational=False, poles=(), flat=False, scale=1.0):
    ku, kv = knots(du, nu, 1.0), knots(dv, nv, 2.0)
    vals = []
    for i in range(nu):
        for j in range(nv):
            x = i * 1.1 + rng.uniform(-.3, .3) * (0 if flat else 1)
            y = j * 0.9 + rng.uniform(-.3, .3) * (0 if flat else 1)
            z = 0.0 if flat else math.sin(i + 2 * j) + rng.uniform(-.4, .4)
            if i == 0 and 'u0' in poles:
                x, y, z = 0.5, 0.25, 0.7
            if i == nu - 1 and 'u1' in poles:
                x, y, z = 3.0, 1.0, -0.2
            if j == 0 and 'v0' in poles:
                x, y, z = 1.0, 0.0, 0.4
            w = rng.uniform(.6, 1.8) if rational else None
            vals += [x * scale * w, y * scale * w, z * scale * w, w] if rational else [x * scale, y * scale, z * scale, NO_W]
    return f'{du} {dv} {nu} {nv} {int(rational)} ' + v(*ku) + ' ' + v(*kv) + ' ' + v(*vals)


SURFACES = [
    surface(1, 1, 2, 2, flat=True), surface(1, 1, 3, 3), surface(2, 2, 3, 4), surface(3, 3, 5, 5),
    surface(2, 3, 4, 5), surface(3, 2, 6, 4), surface(2, 2, 4, 4, rational=True), surface(3, 3, 5, 6, rational=True),
    surface(2, 2, 3, 3, flat=True), surface(2, 2, 4, 4, poles=('u0',)), surface(3, 2, 5, 4, poles=('u0', 'u1')),
    surface(2, 3, 4, 5, poles=('v0',)), surface(2, 2, 4, 4, rational=True, poles=('u0',)),
    surface(3, 3, 5, 5, scale=1e-6), surface(3, 3, 5, 5, scale=1e6),
]
DISTANCES = (0.0, 0.5, -0.5, 2.0, 1e-9)
UVS = [(0.5, 1.0), (0.1, 0.2), (0.9, 1.7), (0.0, 0.0), (1.0, 2.0), (0.0, 1.0), (1.0, 0.5), (-0.5, 0.5), (0.5, 3.0),
       (1e-9, 1e-9), (0.5, 1.0 + 1e-12)]

for s in SURFACES:
    for d in DISTANCES:
        for uv in rng.sample(UVS, 3):
            kind = 'boundary' if uv in ((0.0, 0.0), (1.0, 2.0), (0.0, 1.0), (1.0, 0.5)) else 'regular'
            case('offset-point', kind, h(d), s, v(*uv))
            case('offset-normal', kind, h(d), s, v(*uv), rng.randint(0, 1), rng.randint(0, 1))

ORDERS_UPPER = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2), (2, 1), (1, 2), (3, 0), (0, 3), (2, 2), (3, 1), (3, 3), (4, 0)]
ORDERS_FULL = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2), (2, 1), (1, 2), (3, 0), (0, 3), (2, 2), (4, 0)]
for si, s in enumerate(SURFACES):
    for d in (0.5, -0.5, 2.0, 0.0, 1e-9):
        for uv in rng.sample(UVS, 4):
            for _ in range(2):
                upper = rng.randint(0, 1)
                hu, hv = rng.choice(ORDERS_UPPER if upper else ORDERS_FULL)
                simple = rng.randint(0, 1)
                kind = 'boundary' if uv in ((0.0, 0.0), (1.0, 2.0), (0.0, 1.0), (1.0, 0.5)) else \
                    'singular' if (hu + hv > 3 and not upper) or (max(hu, hv) > 3) else \
                    'tangent' if d == 0.0 else 'regular'
                case('offset-eval', kind, simple, h(d), s, v(*uv), simple and 0 or 1 if False else 0)
                CASES.pop()
                CASES.append(('offset-eval', kind, f'{h(d)} {s} {v(*uv)} {simple} {hu} {hv} {rng.randint(0, 1)} {rng.randint(0, 1)} {upper} {rng.randint(0, 1)}'))
# Poles, where the evaluation extrapolates from two steps off the surface.
for si in (9, 10, 11, 12):
    for d in (0.5, -0.5, 2.0):
        for uv in ((0.0, 0.5), (0.0, 1.0), (0.0, 1.9), (1.0, 1.0), (0.5, 0.0), (0.0, 0.0)):
            for simple in (0, 1):
                for hu, hv, upper in ((1, 1, 0), (2, 1, 0), (1, 2, 0), (3, 0, 0), (2, 2, 1), (3, 3, 1)):
                    CASES.append(('offset-eval', 'touching', f'{h(d)} {SURFACES[si]} {v(*uv)} {simple} {hu} {hv} 1 1 {upper} 1'))


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
            sys.exit(f'case {index} {op} {kind}: oracle exit {run.returncode}: {run.stderr.strip()}')
        response = run.stdout.split()
        assert response[0] == op, (index, response)
        lines.append(f'{op} {kind} | {inputs} | {" ".join(response[1:])}')
    header = ('# catallax-trace offset-surface gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
