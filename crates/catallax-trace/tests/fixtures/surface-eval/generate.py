#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace surface evaluation gate.

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


def evaluate(op, kind, s, uv, hu=0, hv=0, sides=(1, 1), upper=0, nonzero=0):
    assert kind in KINDS
    tail = ('' if op == 'point' else f' {sides[0]} {sides[1]}' if op == 'normal'
            else f' {hu} {hv} {sides[0]} {sides[1]} {upper} {nonzero}')
    CASES.append((op, kind, f'{s} {h(uv[0])} {h(uv[1])}{tail}'))


for s in (BUMP, RBUMP):
    for uv in ((0.37, 2.21), (1.5, 1.5), (2.9, 0.05)):
        evaluate('point', 'regular', s, uv)
        evaluate('eval', 'regular', s, uv, 2, 2)
        evaluate('eval', 'regular', s, uv, 2, 2, upper=1)
        evaluate('eval', 'regular', s, uv, 4, 4, upper=1)
    for sides in ((1, 1), (0, 1), (1, 0), (0, 0)):
        evaluate('eval', 'boundary', s, (1.0, 2.0), 2, 2, sides)
        evaluate('simple', 'boundary', s, (2.0, 1.0), 1, 1, sides)
    for uv in ((0, 0), (3, 3), (0, 3), (-0.5, 1.2), (3.5, 4.0)):
        evaluate('point', 'boundary', s, uv)
        evaluate('eval', 'boundary', s, uv, 1, 1)
    evaluate('eval', 'singular', s, (1.0, 1.0), 5, 5)
    evaluate('eval', 'singular', s, (1.0, 1.0), 2, 1, upper=1)
    evaluate('simple', 'regular', s, (0.7, 2.3), 3, 2)
for uv in ((0.25, 0.5), (0.5, 0.25), (0.75, 0.9), (0.0, 0.0), (1.0, 1.0)):
    evaluate('point', 'regular', LOFT, uv)
    evaluate('eval', 'regular', LOFT, uv, 2, 3)
    evaluate('eval', 'regular', LOFT, uv, 1, 1, nonzero=1)
for uv, kind in (((0.0, 0.5), 'singular'), ((0.0, 0.0), 'singular'), ((1e-9, 0.5), 'near-singular'),
                 ((0.5, 0.5), 'regular')):
    evaluate('point', kind, POLE, uv)
    evaluate('eval', kind, POLE, uv, 1, 1, nonzero=1)
    evaluate('eval', kind, POLE, uv, 2, 2, upper=1, nonzero=1)
    evaluate('eval', kind, POLE, uv, 2, 2, nonzero=0)
    evaluate('eval', kind, POLE, uv, 1, 0, nonzero=1)
# Normals: smooth points, knots from each side, the pole and its
# neighbourhood, where SMLib steps, uses second derivatives or samples.
for s in (BUMP, RBUMP, LOFT):
    for uv in ((0.37, 0.21), (1.0, 1.0), (0.0, 0.0)):
        for sides in ((1, 1), (0, 0)):
            evaluate('normal', 'regular', s, uv, sides=sides)
for uv, kind in (((0.0, 0.5), 'singular'), ((0.0, 0.0), 'singular'), ((0.0, 1.0), 'singular'),
                 ((1e-9, 0.5), 'near-singular'), ((1e-5, 0.3), 'near-singular')):
    for sides in ((1, 1), (0, 1), (1, 0)):
        evaluate('normal', kind, POLE, uv, sides=sides)
for uv in ((1.0, 0.5), (0.0, 1.0), (2.0, 0.0)):
    evaluate('eval', 'touching', FLAT, uv, 2, 2, nonzero=1)
    evaluate('point', 'touching', FLAT, uv)


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
    header = ('# catallax-trace surface evaluation gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
