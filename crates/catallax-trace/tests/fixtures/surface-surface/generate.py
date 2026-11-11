#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace local surface-surface gate.

usage: generate.py ORACLE [RAW_DIR]

ORACLE is the driver built from oracle.cpp against the pinned SMLib. Each
case runs in a fresh process under a memory and CPU cap and a timeout. RAW_DIR,
when given, receives every request, response and error stream; keep it outside
the repository. The B-spline surfaces are those of ../surface-solve/generate.py.
"""
import importlib.util
import pathlib
import random
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('surfaces', HERE.parent / 'surface-solve' / 'generate.py')
surfaces = importlib.util.module_from_spec(spec)
spec.loader.exec_module(surfaces)
h = surfaces.h
KINDS = ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
CASES = []
MINIMIZE, MAXIMIZE, NORMALIZE = 3, 4, 5


def v(*xs):
    return ' '.join(h(x) for x in xs)


def plane(kind, origin, x, y, scale=(1, 1), domain=(0, 0, 1, 1)):
    return f'p {kind} {v(*origin, *x, *y, *scale, *domain)}'


def horizontal(z):
    return plane(1, (0, 0, z), (1, 0, 0), (0, 1, 0))


def lsi(kind, a, da, b, db, ga, gb, tolerance=1e-5, normal=None, previous=None, boundary=1):
    assert kind in KINDS
    n = normal if normal is not None else (0, 0, 0)
    p = previous if previous is not None else (0, 0)
    CASES.append(('lsi', kind, f'{a} {v(*da)} {b} {v(*db)} {v(*ga, *gb, tolerance)} '
                               f'{int(normal is not None)} {v(*n)} {int(previous is not None)} {v(*p)} {boundary}'))


def lss(kind, a, da, b, db, ga, gb, operation):
    assert kind in KINDS
    CASES.append(('lss', kind, f'{a} {v(*da)} {b} {v(*db)} {v(*ga, *gb)} {operation}'))


BUMP, RBUMP, LOFT, POLE, CYLINDER = (f's {x}' for x in (surfaces.BUMP, surfaces.RBUMP, surfaces.LOFT,
                                                        surfaces.POLE, surfaces.CYLINDER))
BOX, UNIT, CBOX, PBOX = (0, 0, 3, 3), (0, 0, 1, 1), (0, 0, 4, 2), (-20, -20, 20, 20)
TILT = plane(1, (1.5, 1.5, 0.2), (1, 0, 0.3), (0, 1, -0.2))
WALL = plane(1, (1, 0, 0), (0, 1, 0), (0, 0, 1))

# Bumps against level and tilted planes from guesses near the meeting curve,
# with the plane moved each step or fixed, with and without side tests.
for s in (BUMP, RBUMP):
    for p, kind in ((horizontal(0.3), 'regular'), (horizontal(-0.2), 'regular'), (TILT, 'regular')):
        for g in ((1.0, 1.5), (2.0, 0.6), (0.4, 2.5), (1.5, 1.5)):
            lsi(kind, s, BOX, p, PBOX, g, g)
            lsi(kind, s, BOX, p, PBOX, g, g, normal=(1, 1, 0))
            lsi(kind, s, BOX, p, PBOX, g, g, boundary=0)
            lsi(kind, s, BOX, p, PBOX, g, g, previous=(g[0] + 0.01, g[1]))
        # Guesses on the bump's sides and the plane's.
        for g in ((0.0, 1.2), (3.0, 2.2), (1.4, 0.0), (2.1, 3.0)):
            lsi('boundary', s, BOX, p, PBOX, g, g)
            lsi('boundary', s, BOX, p, PBOX, g, g, previous=(g[0] + 0.5, g[1] + 0.5))
        lsi('boundary', s, (0.5, 0.5, 2.5, 2.5), p, PBOX, (0.5, 1.5), (0.5, 1.5))
# The loft and the cylinder against planes; the cylinder against a wall
# along its ruling (tangent), just inside it and just outside it.
for g in ((0.3, 0.4), (0.7, 0.2), (0.5, 0.9)):
    lsi('regular', LOFT, UNIT, horizontal(1.2), PBOX, g, (1.0, g[1] * 3))
    lsi('regular', LOFT, UNIT, TILT, PBOX, g, (1.0, g[1] * 3))
for x, kind in ((1.0, 'tangent'), (0.999, 'near-singular'), (1.001, 'touching'), (0.5, 'regular')):
    w = plane(1, (x, 0, 0), (0, 1, 0), (0, 0, 1))
    for g, gw in (((0.1, 1.0), (0.1, 1.0)), ((3.9, 0.5), (-0.1, 0.5)), ((0.3, 1.5), (0.4, 1.5))):
        lsi(kind, CYLINDER, CBOX, w, PBOX, g, gw)
        lsi(kind, CYLINDER, CBOX, w, PBOX, g, gw, normal=(0, 0, 1))
# Near tangency within loose tolerances, where the answer is refined by the
# normal solves: walls 1e-7 and 1e-5 off the cylinder's ruling.
for x, kind in ((1 + 1e-7, 'touching'), (1 - 1e-7, 'near-singular'), (1 + 1e-5, 'touching'),
                (1 - 1e-5, 'near-singular')):
    w = plane(1, (x, 0, 0), (0, 1, 0), (0, 0, 1))
    for tolerance in (1e-4, 1e-3):
        lsi(kind, CYLINDER, CBOX, w, PBOX, (0.02, 1.0), (0.02, 1.0), tolerance)
        lsi(kind, CYLINDER, CBOX, w, PBOX, (3.97, 0.6), (-0.03, 0.6), tolerance, normal=(0, 0, 1))
# Walls meeting the cylinder where its normal turns 2 and 4 degrees from
# theirs, either side of SMLib's 3-degree near-tangency test.
import math
for degrees in (2.0, 4.0):
    x = math.cos(math.radians(degrees))
    w = plane(1, (x, 0, 0), (0, 1, 0), (0, 0, 1))
    lsi('near-singular', CYLINDER, CBOX, w, PBOX, (0.06, 1.0), (0.07, 1.0), 1e-4)
    lsi('near-singular', CYLINDER, CBOX, w, PBOX, (0.06, 1.0), (0.07, 1.0), 1e-4, normal=(0, 0, 1))
# Interior guesses whose answer lies beyond a sub-domain side.
for z in (0.3, -0.2):
    for g in ((1.9, 1.5), (1.1, 1.2), (1.5, 1.95), (1.2, 1.05)):
        lsi('boundary', BUMP, (1, 1, 2, 2), horizontal(z), PBOX, g, g)
# The pole patch against planes through and near its pole.
for z, kind in ((1.0, 'singular'), (1.0 - 1e-9, 'near-singular'), (0.9, 'regular')):
    for g in ((0.05, 0.5), (0.0, 0.2), (0.3, 0.7)):
        lsi(kind, POLE, UNIT, horizontal(z), PBOX, g, (0.05, 0.0))
# Two splines: the cylinder and a bump, the bump and its rational copy.
for g, gb in (((0.4, 0.3), (0.95, 0.3)), ((0.1, 0.8), (0.9, 0.5))):
    lsi('regular', CYLINDER, CBOX, BUMP, BOX, g, gb)
lsi('touching', BUMP, BOX, RBUMP, BOX, (1.5, 1.5), (1.5, 1.5))
# Nearest, farthest and normal pairs.
for a, da, b, db, ga, gb in ((BUMP, BOX, horizontal(3.0), PBOX, (1.5, 1.5), (1.5, 1.5)),
                             (BUMP, BOX, horizontal(-3.0), PBOX, (0.5, 2.5), (0.5, 2.5)),
                             (CYLINDER, CBOX, BUMP, BOX, (0.4, 0.5), (1.0, 0.5)),
                             (LOFT, UNIT, TILT, PBOX, (0.5, 0.5), (1.0, 1.0)),
                             (BUMP, BOX, RBUMP, BOX, (1.0, 1.0), (2.0, 2.0)),
                             (POLE, UNIT, horizontal(2.0), PBOX, (0.0, 0.5), (0.0, 0.0))):
    for operation in (MINIMIZE, MAXIMIZE, NORMALIZE):
        kind = 'singular' if a is POLE else 'regular'
        lss(kind, a, da, b, db, ga, gb, operation)
# A seeded sweep of bump/plane guesses.
rng = random.Random(23)
for _ in range(30):
    g = (rng.uniform(0, 3), rng.uniform(0, 3))
    p = horizontal(rng.uniform(-0.5, 0.8))
    lsi('regular', BUMP, BOX, p, PBOX, g, (g[0] + rng.uniform(-0.2, 0.2), g[1] + rng.uniform(-0.2, 0.2)),
        normal=None if rng.random() < 0.7 else (0, 1, 0))


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
    header = ('# catallax-trace local surface-surface gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
