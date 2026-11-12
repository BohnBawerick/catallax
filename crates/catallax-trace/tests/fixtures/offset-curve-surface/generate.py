#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace offset curve-surface gate.

usage: generate.py ORACLE [RAW_DIR]

ORACLE is the driver built from oracle.cpp against the pinned SMLib. Each
case runs in a fresh process under a memory and CPU cap and a timeout. RAW_DIR,
when given, receives every request, response and error stream; keep it outside
the repository. Inputs are fixed: hand-picked tangent, touching, singular and
boundary cases plus seeded random regular cases.
"""
import pathlib
import random
import struct
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
NO_W = -sys.float_info.max
MINIMIZE, NORMALIZE = 3, 5


def h(x):
    return format(struct.unpack('<Q', struct.pack('<d', float(x)))[0], '016x')


def v(*xs):
    return ' '.join(h(x) for x in xs)


CASES = []


def case(op, kind, *tokens):
    assert kind in ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
    CASES.append((op, kind, ' '.join(str(t) for t in tokens)))


rng = random.Random(20261112)


def knots(degree, count):
    interior = count - degree - 1
    return [0.0] * (degree + 1) + [(i + 1) / (interior + 1) for i in range(interior)] + [1.0] * (degree + 1)


def net(du, dv, rows, rational=False):
    """A surface from its control rows `rows[i][j] = (x, y, z)`."""
    nu, nv = len(rows), len(rows[0])
    vals = []
    for row in rows:
        for x, y, z in row:
            if rational:
                w = rng.uniform(.8, 1.2)
                vals += [x * w, y * w, z * w, w]
            else:
                vals += [x, y, z, NO_W]
    return f'{du} {dv} {nu} {nv} {int(rational)} ' + v(*knots(du, nu)) + ' ' + v(*knots(dv, nv)) + ' ' + v(*vals)


def patch(degree, fn, bump=0.0, rational=False):
    """A bicubic (or other degree) patch of `fn(s, t)` over `[0,1]^2`, interior points perturbed by `bump`."""
    n = degree + 1
    rows = []
    for i in range(n):
        row = []
        for j in range(n):
            x, y, z = fn(i / (n - 1), j / (n - 1))
            if 0 < i < n - 1 and 0 < j < n - 1:
                x += rng.uniform(-bump, bump); y += rng.uniform(-bump, bump); z += rng.uniform(-bump, bump)
            row.append((x, y, z))
        rows.append(row)
    return net(degree, degree, rows, rational)


def floor(bump=0.0, rational=False):
    # normal +z
    return patch(3, lambda s, t: (3 * s, 2 * t, 0.0), bump, rational)


def wall(bump=0.0, rational=False):
    # normal +x
    return patch(3, lambda s, t: (0.0, 2 * s, 3 * t), bump, rational)


def arch(a):
    # (3s, 2t, 4a s(1-s)): highest at s = 1/2, height a, normal +z there
    return net(2, 1, [[(0.0, 0.0, 0.0), (0.0, 2.0, 0.0)],
                      [(1.5, 0.0, 2 * a), (1.5, 2.0, 2 * a)],
                      [(3.0, 0.0, 0.0), (3.0, 2.0, 0.0)]])


def fan():
    # a curved fan whose u = 0 side collapses to the pole (0, 1, 0)
    return net(1, 2, [[(0.0, 1.0, 0.0)] * 3,
                      [(3.0, 0.0, 0.0), (3.0, 1.0, 0.6), (3.0, 2.0, 0.0)]])


def fan_v():
    # the fan with u and v swapped: its v = 0 side collapses to the pole
    return net(2, 1, [[(0.0, 1.0, 0.0), (3.0, 0.0, 0.0)],
                      [(0.0, 1.0, 0.0), (3.0, 1.0, 0.6)],
                      [(0.0, 1.0, 0.0), (3.0, 2.0, 0.0)]])


UNIT = (0.0, 0.0, 1.0, 1.0)


def intersect(kind, da, a, db, b, axis, iso, clamp, domain, tol, t, uv):
    case('oc-intersect', kind, h(da), a, h(db), b, axis, h(iso), clamp, v(*domain), v(tol, t, *uv))


def solve(kind, da, a, db, b, axis, iso, clamp, domain, operation, t, uv):
    case('oc-solve', kind, h(da), a, h(db), b, axis, h(iso), clamp, v(*domain), operation, v(1e-6, t, *uv))


def sing(kind, d, s, uv, tol, point_only):
    case('oc-sing', kind, h(d), s, v(*uv, tol), point_only)


# Regular crossings: the iso-curve of one offset plane-like patch meets the
# other, as a fillet's clipped step asks it (floor at height r1, wall at x = r2).
PAIRS = [(floor(), wall()), (floor(0.06), wall(0.06)), (floor(0.08, True), wall(0.05, True))]
for f, w in PAIRS:
    for r1, r2 in ((0.4, 0.4), (0.3, 0.5), (-0.2, 0.6)):
        for tol in (1e-6, 1e-9):
            for iso in (0.0, 1.0, rng.uniform(0.1, 0.9)):
                # floor iso v = const runs along x and meets the wall offset at x = r2
                t_exact = r2 / 3.0
                for t, uv in ((t_exact, (iso, t_exact)), (rng.uniform(0, 1), (rng.uniform(0, 1), rng.uniform(0, 1)))):
                    intersect('regular', r1, f, r2, w, 1, iso, 1, UNIT, tol, t, uv)
                # wall iso u = const runs along z and meets the floor offset at z = r1
                t_exact = r1 / 3.0 if r1 > 0 else 0.1
                intersect('regular', r2, w, r1, f, 0, iso, 1, UNIT, tol, t_exact, (rng.uniform(0, 1), iso))
    # the curve parameter and surface parameters solved for together
    for operation in (MINIMIZE, NORMALIZE):
        for _ in range(3):
            solve('regular', 0.4, f, 0.4, w, 1, rng.uniform(0, 1), 1, UNIT, operation, rng.uniform(0, 1),
                  (rng.uniform(0, 1), rng.uniform(0, 1)))

# Tangent: the floor offset's iso-line at the arch offset's crest height grazes
# it at s = 1/2; touching just within tolerance; two crossings just above it.
A, D = 0.5, 0.2
for crest, kind in ((A + D, 'tangent'), (A + D - 1e-7, 'touching'), (A + D - 1e-3, 'regular'), (A + D + 1e-3, 'regular')):
    for iso in (0.5, 0.25):
        for t, uv in ((0.5, (0.5, iso)), (0.45, (0.55, iso)), (0.2, (0.3, iso))):
            for tol in (1e-6, 1e-3):
                intersect(kind, crest, floor(), D, arch(A), 1, iso, 1, UNIT, tol, t, uv)
        for operation in (MINIMIZE, NORMALIZE):
            solve(kind, crest, floor(), D, arch(A), 1, iso, 1, UNIT, operation, 0.45, (0.55, iso))

# Singular: the fan's offset collapses (and extrapolates its normal) at its
# pole; iso-lines of a floor offset pass near it.
for d in (0.3, -0.3):
    for s in (fan(), fan_v()):
        for uv in ((0.0, 0.5), (0.0, 0.0), (0.5, 0.0), (1e-9, 0.5), (0.5, 1e-9), (0.5, 0.5)):
            kind = 'near-singular' if 1e-9 in uv else 'singular' if 0.0 in uv else 'regular'
            for tol in (1e-6, 1e-2):
                sing(kind, d, s, uv, tol, 0)
            sing(kind, d, s, uv, 1e-6, 1)
        for height in (0.0, 0.25, -0.25):
            for iso in (0.5, 0.3):
                for t, uv in ((0.02, (0.0, 0.5)), (0.02, (0.5, 0.0)), (0.3, (0.1, 0.5)), (0.3, (0.5, 0.1))):
                    intersect('singular', height, floor(), d, s, 1, iso, 1, UNIT, 1e-6, t, uv)
            solve('singular', height, floor(), d, s, 1, 0.5, 1, UNIT, MINIMIZE, 0.02, (0.0, 0.5))
for d in (0.0, 0.5):
    for uv in ((0.5, 0.5), (0.0, 0.0), (1.0, 1.0)):
        sing('regular', d, floor(0.06), uv, 1e-6, 0)

# Boundary: iso-curves on and past the edges (SMLib's clipped step does not
# clamp its iso parameter; the offset surface clamps the evaluation), guesses
# outside the domain, sub-domains that hold the crossing and that miss it.
f, w = PAIRS[1]
for iso, clamp in ((0.0, 1), (1.0, 1), (1.2, 0), (-0.1, 0), (1.2, 1), (-0.1, 1), (1.0 + 2 ** -52, 0)):
    intersect('boundary', 0.4, f, 0.4, w, 1, iso, clamp, UNIT, 1e-6, 0.13, (0.5, 0.13))
    intersect('boundary', 0.4, w, 0.4, f, 0, iso, clamp, UNIT, 1e-6, 0.13, (0.13, 0.5))
for domain in ((0.0, 0.0, 0.5, 0.5), (0.5, 0.5, 1.0, 1.0), (0.2, 0.0, 0.6, 1.0), (0.0, 0.0, 0.1, 0.1)):
    for t, uv in ((0.13, (0.3, 0.13)), (1.5, (-0.5, 1.5)), (-1.0, (2.0, -1.0))):
        intersect('boundary', 0.4, f, 0.4, w, 1, 0.3, 1, domain, 1e-6, t, uv)
# Misses: the iso-line passes above the arch offset.
for t, uv in ((0.5, (0.5, 0.5)), (0.1, (0.9, 0.1))):
    intersect('regular', A + D + 0.5, floor(), D, arch(A), 1, 0.5, 1, UNIT, 1e-6, t, uv)

# The pole on the high side: a meeting at the pole moves to the side of the
# domain that holds it.
def fan_high():
    # the fan mirrored: its u = 1 side collapses to the pole (3, 1, 0)
    return net(1, 2, [[(0.0, 0.0, 0.0), (0.0, 1.0, 0.6), (0.0, 2.0, 0.0)],
                      [(3.0, 1.0, 0.0)] * 3])


# An offset of a pole is a curve, not a point, unless the distance is zero.
for d in (0.3, -0.3, 0.0):
    for uv in ((1.0, 0.5), (1.0 - 1e-9, 0.5)):
        sing('singular' if uv[0] == 1.0 else 'near-singular', d, fan_high(), uv, 1e-6, 0)
    for height in (0.0, 0.25, -0.25):
        for domain in (UNIT, (0.5, 0.0, 1.0, 1.0), (0.0, 0.0, 0.999, 1.0)):
            for t, uv in ((0.98, (1.0, 0.5)), (0.9, (0.9, 0.5)), (0.99, (0.99, 0.4))):
                intersect('singular', height, floor(), d, fan_high(), 1, 0.5, 1, domain, 1e-6, t, uv)


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
    header = ('# catallax-trace offset curve-surface gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
