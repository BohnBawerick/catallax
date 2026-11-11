#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace surface point solver gate.

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
# A cylinder of radius 1 about z, closed in u with its seam at u = 0 and 4.
r = math.sqrt(0.5)
ring = [(1, 0, 1), (1, 1, r), (0, 1, 1), (-1, 1, r), (-1, 0, 1), (-1, -1, r), (0, -1, 1), (1, -1, r), (1, 0, 1)]
CYLINDER = surface(2, 1, [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4], [0, 0, 2, 2],
                   [(x, y, z, w) for x, y, w in ring for z in (0, 2)], rational=True)
MINIMIZE, NORMALIZE, INTERSECT = 3, 5, 0


def solve(kind, s, domain, operation, test, guess, op='solve'):
    assert kind in KINDS
    tokens = ' '.join(h(x) for x in (*domain, )) + f' {operation} ' + ' '.join(h(x) for x in (*test, *guess))
    CASES.append((op, kind, f'{s} {tokens}'))


BOX = (0, 0, 3, 3)
for s in (BUMP, RBUMP):
    for test, guess in (((1.2, 1.7, 2.0), (1.5, 1.5)), ((2.2, 0.4, -1.5), (0.5, 2.5)), ((0.3, 2.8, 0.9), (2.9, 0.1))):
        for operation in (MINIMIZE, NORMALIZE):
            solve('regular', s, BOX, operation, test, guess)
    # Points beyond each side and corner: the answer stops on the boundary.
    for test, guess in (((-1.0, 1.5, 0.5), (1.5, 1.5)), ((4.0, 4.0, 1.0), (1.0, 1.0)), ((1.5, -2.0, 0.0), (1.5, 2.0)),
                        ((-0.5, -0.5, 3.0), (0.2, 0.2))):
        for operation in (MINIMIZE, NORMALIZE):
            solve('boundary', s, BOX, operation, test, guess)
    # A subdomain whose sides cut the surface.
    for test, guess in (((0.5, 1.5, 2.0), (1.5, 1.5)), ((1.5, 2.6, 0.0), (1.2, 1.2))):
        solve('boundary', s, (1, 1, 2, 2), MINIMIZE, test, guess)
    # A point on the surface: the gap closes.
    solve('touching', s, BOX, INTERSECT, (1.0, 1.0, 0.0), (0.8, 1.3))
    solve('touching', s, BOX, MINIMIZE, (0.0, 0.0, 0.0), (1.0, 1.0))
    solve('regular', s, BOX, INTERSECT, (1.0, 1.0, 5.0), (1.0, 1.0))
    # Starting on interior knots and on a corner.
    solve('regular', s, BOX, MINIMIZE, (2.0, 1.0, 1.5), (1.0, 2.0))
    solve('boundary', s, BOX, MINIMIZE, (2.0, 1.0, 1.5), (3.0, 0.0))
for test, guess in (((0.5, 1.0, 0.5), (0.5, 0.5)), ((1.5, 3.0, 2.0), (0.1, 0.9)), ((0.0, 0.0, 0.0), (0.9, 0.1))):
    for operation in (MINIMIZE, NORMALIZE):
        solve('regular', LOFT, (0, 0, 1, 1), operation, test, guess)
# Near the pole of POLE.
for test, guess, kind in (((0.5, 0.0, 2.0), (0.5, 0.5), 'regular'), ((0.1, 0.3, 1.5), (0.3, 0.6), 'near-singular'),
                          ((-0.5, 0.0, 1.2), (0.5, 0.5), 'singular'), ((0.0, 0.0, 1.0), (0.4, 0.4), 'singular')):
    solve(kind, POLE, (0, 0, 1, 1), MINIMIZE, test, guess)
# Iterates on the pole itself, where SMLib inverts the point near the
# singularity: at the pole, on its normal, and off it.
for guess in ((0.0, 0.2), (0.0, 0.5), (0.0, 0.8)):
    for test in ((0.3, 0.5, 1.5), (0.2, -0.4, 0.8), (-0.3, 0.0, 1.0), (0.0, 0.0, 1.5), (0.0, 0.0, 1.0 + 1e-9)):
        solve('singular', POLE, (0, 0, 1, 1), MINIMIZE, test, guess)
    solve('singular', POLE, (0, 0, 1, 1), NORMALIZE, (0.3, 0.5, 1.5), guess)
# From the pole: test points within the tight tolerance of the pole, and
# points the pole's cross tangent turns towards, which SMLib finishes by
# dropping the point on the iso-curve through the pole.
n = (0.4660849584677394, -0.25422815916422137, 0.8474271972140716)
for test in ((1e-13, 0.0, 1.0), (0.0, 1e-13, 1.0 - 1e-13)):
    solve('singular', POLE, (0, 0, 1, 1), MINIMIZE, test, (0.0, 0.5))
for test in ((0.5 * n[0], 0.5 * n[1], 1.0 + 0.5 * n[2]), (0.5, 0.0, 0.725), (0.4, 0.3, 0.8)):
    solve('singular', POLE, (0, 0, 1, 1), MINIMIZE, test, (0.0, 0.5))
    for guess in ((0.0, 0.1), (0.0, 0.9), (0.0, 0.0), (0.0, 1.0)):
        solve('singular', POLE, (0, 0, 1, 1), MINIMIZE, test, guess)
    solve('singular', POLE, (0, 0, 1, 1), NORMALIZE, test, (0.0, 0.5))
    solve('singular', POLE, (0, 0, 1, 1), INTERSECT, test, (0.0, 0.5))
# Points far beyond the pole along and across its cross tangents.
for test in ((-1.0, 0.2, 2.0), (2.0, 3.0, 0.0), (0.1, -3.0, 1.2)):
    solve('singular', POLE, (0, 0, 1, 1), MINIMIZE, test, (0.0, 0.5))
# Bilinear patch: steps may span the whole domain.
for test, guess in (((1.0, 0.5, 2.0), (0.0, 0.0)), ((3.0, 2.0, 0.0), (1.0, 0.5)), ((0.5, 0.5, 0.125), (2.0, 1.0))):
    solve('regular', FLAT, (0, 0, 2, 1), MINIMIZE, test, guess)
solve('touching', FLAT, (0, 0, 2, 1), INTERSECT, (1.0, 0.5, 0.25), (0.2, 0.9))
# Cylinder: the nearest point lies across the seam from the guess, and a
# point on the axis, where every direction is equally near.
for test, guess, kind in (((2.0, -0.3, 1.0), (0.4, 1.0), 'boundary'), ((2.0, 0.3, 1.0), (3.6, 1.0), 'boundary'),
                          ((3.0, 0.0, 1.0), (0.0, 1.0), 'boundary'), ((0.0, 0.0, 1.0), (1.3, 0.5), 'tangent'),
                          ((0.0, 0.0, -1.0), (2.5, 0.5), 'tangent'), ((0.5, 0.5, 1.0), (0.5, 1.5), 'regular')):
    solve(kind, CYLINDER, (0, 0, 4, 2), MINIMIZE, test, guess)
solve('regular', CYLINDER, (0, 0, 4, 2), NORMALIZE, (2.0, 2.0, 1.0), (0.5, 0.5))
solve('regular', CYLINDER, (0.5, 0, 3.5, 2), NORMALIZE, (2.0, -0.1, 1.0), (1.0, 0.5))

# Points deep inside the cylinder, beyond the centre of curvature of the
# guess, where SMLib limits the second-derivative terms of the Jacobian.
for test, guess in (((0.05, -0.02, 1.0), (1.0, 0.5)), ((0.05, -0.02, 1.0), (2.6, 1.5)), ((-0.08, 0.03, 0.5), (0.3, 1.0)),
                    ((0.02, 0.09, 1.7), (3.1, 0.2))):
    solve('regular', CYLINDER, (0, 0, 4, 2), MINIMIZE, test, guess)
    solve('regular', CYLINDER, (0, 0, 4, 2), NORMALIZE, test, guess)
# A fixed pseudo-random sweep of test points and guesses.
state = 12345


def uniform(low, high):
    global state
    state = (1103515245 * state + 12345) % (1 << 31)
    return low + (high - low) * state / (1 << 31)


for s_, box, lows, highs in ((BUMP, BOX, (-1, -1, -2), (4, 4, 3)), (RBUMP, BOX, (-1, -1, -2), (4, 4, 3)),
                             (LOFT, (0, 0, 1, 1), (-1, -1, -1), (3, 3, 4)),
                             (CYLINDER, (0, 0, 4, 2), (-2, -2, -0.5), (2, 2, 2.5)),
                             (POLE, (0, 0, 1, 1), (-0.5, -2, -0.5), (4, 2, 2))):
    for _ in range(16):
        test = tuple(uniform(a, b) for a, b in zip(lows, highs))
        guess = (uniform(box[0], box[2]), uniform(box[1], box[3]))
        operation = MINIMIZE if uniform(0, 1) < 0.7 else NORMALIZE
        solve('regular', s_, box, operation, test, guess)


def global_solve(kind, s, domain, operation, test, tolerance):
    assert kind in KINDS
    tokens = ' '.join(h(x) for x in domain) + f' {operation} ' + ' '.join(h(x) for x in (*test, tolerance))
    CASES.append(('global', kind, f'{s} {tokens}'))


# Global solves: points above and beside each surface, on it, beyond its
# sides, across the cylinder seam, on its axis, and at the pole.
for s_, box in ((BUMP, BOX), (RBUMP, BOX), (LOFT, (0, 0, 1, 1)), (FLAT, (0, 0, 2, 1))):
    for test, kind in (((1.2, 1.7, 2.0), 'regular'), ((2.2, 0.4, -1.5), 'regular'), ((-1.0, 1.5, 0.5), 'boundary'),
                       ((4.0, 4.0, 1.0), 'boundary'), ((0.5, 0.5, 0.1), 'regular')):
        for operation in (MINIMIZE, NORMALIZE):
            global_solve(kind, s_, box, operation, test, 1e-5)
    global_solve('touching', s_, box, INTERSECT, (1.0, 1.0, 0.0), 1e-5)
    global_solve('regular', s_, box, INTERSECT, (1.0, 1.0, 3.0), 1e-3)
global_solve('touching', BUMP, BOX, INTERSECT, (1.5, 1.5, 0.0), 1e-5)
for test, kind in (((2.0, -0.3, 1.0), 'boundary'), ((2.0, 0.3, 1.0), 'boundary'), ((3.0, 0.0, 1.0), 'boundary'),
                   ((0.0, 0.0, 1.0), 'tangent'), ((0.05, -0.02, 1.0), 'regular'), ((0.5, 0.5, 1.0), 'regular'),
                   ((1.0, 0.0, 2.5), 'boundary')):
    for operation in (MINIMIZE, NORMALIZE):
        global_solve(kind, CYLINDER, (0, 0, 4, 2), operation, test, 1e-5)
global_solve('touching', CYLINDER, (0, 0, 4, 2), INTERSECT, (1.0, 0.0, 1.0), 1e-5)
for test, kind in (((0.5, 0.0, 2.0), 'singular'), ((0.0, 0.0, 1.0), 'singular'), ((-0.5, 0.0, 1.2), 'singular'),
                   ((0.1, 0.3, 1.5), 'near-singular'), ((2.0, 0.0, 0.5), 'regular')):
    global_solve(kind, POLE, (0, 0, 1, 1), MINIMIZE, test, 1e-5)
# Points just off the cylinder seam, where the answer lands near it.
for offset in (-1e-6, 1e-6, -1e-9, 0.0):
    global_solve('boundary', CYLINDER, (0, 0, 4, 2), MINIMIZE, (2.0, offset, 1.0), 1e-5)
# A rational fan whose first row of control points lies on a circle about
# the test point: every point of that row ties.
FAN = surface(3, 1, [0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 1, 1],
              [(r * math.cos(a), r * math.sin(a), z, w) for a in (0, math.pi / 3, 2 * math.pi / 3, math.pi)
               for r, z, w in ((1.0, 0.0, 1.0), (2.0, 1.0, 1.5))], rational=True)
for test in ((0.0, 0.0, 0.0), (0.0, 0.0, -0.5), (0.0, 0.0, 0.3)):
    global_solve('tangent', FAN, (0, 0, 1, 1), MINIMIZE, test, 1e-5)
for s_, box, lows, highs in ((BUMP, BOX, (-1, -1, -2), (4, 4, 3)), (RBUMP, BOX, (-1, -1, -2), (4, 4, 3)),
                             (LOFT, (0, 0, 1, 1), (-1, -1, -1), (3, 3, 4)),
                             (CYLINDER, (0, 0, 4, 2), (-2, -2, -0.5), (2, 2, 2.5))):
    for _ in range(10):
        test = tuple(uniform(a, b) for a, b in zip(lows, highs))
        operation = MINIMIZE if uniform(0, 1) < 0.7 else NORMALIZE
        global_solve('regular', s_, box, operation, test, 1e-5)


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
    header = ('# catallax-trace surface point solver gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
