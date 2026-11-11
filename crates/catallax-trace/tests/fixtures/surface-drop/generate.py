#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace surface point drop gate.

usage: generate.py ORACLE [RAW_DIR]

ORACLE is the driver built from oracle.cpp against the pinned SMLib. Each
case runs in a fresh process under a memory and CPU cap and a timeout. RAW_DIR,
when given, receives every request, response and error stream; keep it outside
the repository. The surfaces are those of ../surface-solve/generate.py.
"""
import importlib.util
import pathlib
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('surfaces', HERE.parent / 'surface-solve' / 'generate.py')
surfaces = importlib.util.module_from_spec(spec)
spec.loader.exec_module(surfaces)
h = surfaces.h
KINDS = ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
CASES = []
INTERSECT, MINIMIZE, NORMALIZE = 0, 3, 5


def v(*xs):
    return ' '.join(h(x) for x in xs)


def plane(kind, origin, x, y, scale=(1, 1), domain=(0, 0, 1, 1)):
    return f'p {kind} {v(*origin, *x, *y, *scale, *domain)}'


def drop(kind, s, domain, point, guess=None, operation=MINIMIZE):
    assert kind in KINDS
    g = guess if guess is not None else (0, 0)
    CASES.append(('drop', kind, f'{s} {v(*domain, *point)} {int(guess is not None)} {v(*g)} {operation}'))


BUMP, RBUMP, LOFT, POLE, FLAT, CYLINDER = (f's {x}' for x in (surfaces.BUMP, surfaces.RBUMP, surfaces.LOFT,
                                                              surfaces.POLE, surfaces.FLAT, surfaces.CYLINDER))
INF = 1234567.0
PLANE = plane(1, (0, 0, 1), (1, 0, 0), (0, 1, 0))
SCALED = plane(0, (0.1, 0.7, -0.3), (0.36, 0.48, 0.8), (0.8, -0.6, 0), (0.3, 1.7), (-0.7, -1.1, 2.9, 3.3))
BOX, UNIT, CBOX = (0, 0, 3, 3), (0, 0, 1, 1), (0, 0, 4, 2)

for s, domain in ((BUMP, BOX), (RBUMP, BOX), (BUMP, (0.5, 0.5, 2.5, 2.5))):
    for point, kind in (((1.4, 1.6, 0.9), 'regular'), ((1.4, 1.6, 4.0), 'regular'), ((-1.0, 1.5, 0.0), 'boundary'),
                        ((3.5, 3.5, 2.0), 'boundary'), ((2.0, 1.0, 0.4), 'regular')):
        for guess in (None, (1.5, 1.5), (0.1, 2.9), (2.0, 1.0)):
            for operation in (MINIMIZE, NORMALIZE, INTERSECT):
                drop(kind, s, domain, point, guess, operation)
# The loft, the pole patch (its pole and points it turns towards), the
# bilinear patch, and the cylinder: its axis (every direction equally near),
# points across its seam and on it.
for s, domain, points in ((LOFT, UNIT, ((1, 0.5, 0.5), (1, 3, 0.5), (-1, -1, -1))),
                          (POLE, UNIT, ((0, 0, 1), (0.3, 0.5, 1.5), (0.5, 0, 0.725), (1e-9, 0, 1))),
                          (FLAT, (0, 0, 2, 1), ((1, 0.5, 2), (3, 2, 0), (0.5, 0.5, 0.125))),
                          (CYLINDER, CBOX, ((0, 0, 1), (2, 0.01, 1), (2, -0.01, 1), (1, 0, 1), (0.5, 0.5, 3)))):
        for point in points:
            kind = 'singular' if point in ((0, 0, 1), (1e-9, 0, 1)) else 'regular'
            for guess in (None, (0.05, 0.5), (3.9, 1.0)):
                for operation in (MINIMIZE, NORMALIZE, INTERSECT):
                    drop(kind, s, domain, point, guess, operation)
# Points meeting the surfaces exactly and 1e-6 off, against the
# intersection tolerance of 2e-5.
for s, domain, point, kind in ((BUMP, BOX, None, 'touching'), (PLANE, (-INF, -INF, INF, INF), (0.3, 0.4, 1.0), 'touching'),
                               (PLANE, (-INF, -INF, INF, INF), (0.3, 0.4, 1.0 + 1e-6), 'near-singular'),
                               (PLANE, (-INF, -INF, INF, INF), (0.3, 0.4, 1.0 + 1e-4), 'regular'),
                               (PLANE, (0, 0, 0.2, 0.2), (0.3, 0.4, 1.0), 'boundary'),
                               (SCALED, (-0.7, -1.1, 2.9, 3.3), (0.5, 0.2, 0.4), 'regular'),
                               (SCALED, (0.1, 0.2, 0.9, 0.7), (3, 3, 3), 'boundary')):
    if point is None:
        continue
    for guess in (None, (0.25, 0.35), (0.5, 0.5)):
        for operation in (MINIMIZE, NORMALIZE, INTERSECT):
            drop(kind, s, domain, point, guess, operation)
# The cylinder's seam: a point across it gives SMLib's two seam answers,
# picked by the guess side.
for guess in (None, (0.1, 1.0), (3.9, 1.0)):
    drop('boundary', CYLINDER, CBOX, (2, 0, 1), guess, MINIMIZE)
    drop('boundary', CYLINDER, CBOX, (1, 0, 0.5), guess, INTERSECT)
# A seeded sweep of points and guesses over sub-domains of the bumps and the
# loft, guesses on and near the sides.
import random
rng = random.Random(17)
for s, domain, low, high in ((BUMP, (0.5, 0.5, 2.5, 2.5), (-0.5, -0.5, -1.5), (3.5, 3.5, 2.0)),
                             (RBUMP, (1.0, 0.2, 2.0, 2.8), (-0.5, -0.5, -1.5), (3.5, 3.5, 2.0)),
                             (LOFT, (0.2, 0.1, 0.9, 0.8), (-0.5, -1.0, -0.5), (2.5, 3.0, 3.5))):
    for _ in range(25):
        point = tuple(rng.uniform(a, b) for a, b in zip(low, high))
        side = rng.choice((0, 1, 2))
        gu = (domain[0], domain[2], rng.uniform(domain[0], domain[2]))[side]
        gv = rng.uniform(domain[1], domain[3])
        drop('boundary' if side < 2 else 'regular', s, domain, point, (gu, gv),
             rng.choice((MINIMIZE, MINIMIZE, NORMALIZE)))
# A bump point 1e-7 above the surface, with a guess at its foot, and the
# tangent touch of a bump at a local top.
drop('touching', BUMP, BOX, (1.5, 1.5, 0.9), (1.5, 1.5), INTERSECT)
drop('tangent', BUMP, BOX, (1.5, 1.5, 0.9), None, NORMALIZE)


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
    header = ('# catallax-trace surface point drop gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
