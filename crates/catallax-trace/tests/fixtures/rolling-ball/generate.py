#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace rolling-ball gate.

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


rng = random.Random(20261117)


def knots(degree, count):
    interior = count - degree - 1
    return [0.0] * (degree + 1) + [(i + 1) / (interior + 1) for i in range(interior)] + [1.0] * (degree + 1)


def surface(du, dv, nu, nv, fn, bump=0.0, rational=False):
    """A patch of the map `fn(s, t)` over `[0,1]^2` with `bump`-sized perturbations."""
    ku, kv = knots(du, nu), knots(dv, nv)
    vals = []
    for i in range(nu):
        for j in range(nv):
            s, t = i / (nu - 1), j / (nv - 1)
            x, y, z = fn(s, t)
            if 0 < i < nu - 1 and 0 < j < nv - 1:
                x += rng.uniform(-bump, bump); y += rng.uniform(-bump, bump); z += rng.uniform(-bump, bump)
            w = rng.uniform(.8, 1.2) if rational else None
            vals += [x * w, y * w, z * w, w] if rational else [x, y, z, NO_W]
    return f'{du} {dv} {nu} {nv} {int(rational)} ' + v(*ku) + ' ' + v(*kv) + ' ' + v(*vals)


def floor(bump, rational=False, degree=3):
    return surface(degree, degree, degree + 1, degree + 1, lambda s, t: (3 * s, 2 * t, 0.0), bump, rational)


def wall(bump, rational=False, degree=3):
    return surface(degree, degree, degree + 1, degree + 1, lambda s, t: (0.0, 2 * s, 3 * t), bump, rational)


def ceiling(bump):
    # normal along -z: u runs along y, v along x
    return surface(3, 3, 4, 4, lambda s, t: (3 * t, 2 * s, 1.0), bump)


def head(r1, r2, o1, o2, tol, s1, s2):
    return f'{h(r1)} {h(r2)} {o1} {o2} {h(tol)} {s1} {s2}'


def domain():
    return v(0.0, 0.0, 1.0, 1.0)


PAIRS = [(floor(0.0), wall(0.0)), (floor(0.06), wall(0.06)), (floor(0.08, True), wall(0.08, degree=2))]
RADII = [(0.4, 0.4), (0.3, 0.5)]
ORIENT = [(1, 1), (0, 1), (0, 0)]


def consistent(r1, r2, y):
    """uv on each surface of the offset-surface crossing at height y (flat patches)."""
    return (r2 / 3.0, y / 2.0), (y / 2.0, r1 / 3.0)


for s1, s2 in PAIRS:
    for r1, r2 in RADII:
        for o1, o2 in ORIENT:
            tol = rng.choice((1e-6, 1e-9))
            for y in (0.3, 1.5):
                (u1, v1), (u2, v2) = consistent(r1, r2, y)
                # surface values at a point on the crossing, and off it
                case('rb-surface', 'regular', head(r1, r2, o1, o2, tol, s1, s2), v(u1, v1, u2, v2))
                case('rb-surface', 'regular', head(r1, r2, o1, o2, tol, s1, s2),
                     v(rng.uniform(0, 1), rng.uniform(0, 1), rng.uniform(0, 1), rng.uniform(0, 1)))
            if (o1, o2) == (1, 1):
                for y in (0.3, 1.5):
                    (u1, v1), (u2, v2) = consistent(r1, r2, y)
                    for direction in (1.0, -1.0):
                        case('rb-point', 'regular', head(r1, r2, o1, o2, tol, s1, s2), v(u1, v1, u2, v2),
                             h(direction), domain(), domain(), 0)
                    case('rb-point', 'regular', head(r1, r2, o1, o2, tol, s1, s2), v(u1, v1, u2, v2), h(1.0),
                         domain(), domain(), 1, v(0.0, 1.0, 0.0))
                    for step in (0.02, 0.2):
                        for sign in (1.0, -1.0):
                            case('rb-step', 'regular', head(r1, r2, o1, o2, tol, s1, s2),
                                 v(u1, v1, u2, v2), v(0.0, 0.5 * sign, 0.5 * sign, 0.0), v(r2, y, r1), v(0.0, sign, 0.0),
                                 h(step * sign), domain(), domain())
# Steps near the ends of the domain are clipped by it.
for s1, s2 in PAIRS[:2]:
    (u1, v1), (u2, v2) = consistent(0.4, 0.4, 1.95)
    for step in (0.2, 0.05, 0.5):
        case('rb-step', 'boundary', head(0.4, 0.4, 1, 1, 1e-6, s1, s2), v(u1, v1, u2, v2), v(0.0, 0.5, 0.5, 0.0),
             v(0.4, 1.95, 0.4), v(0.0, 1.0, 0.0), h(step), domain(), domain())
    (u1, v1), (u2, v2) = consistent(0.4, 0.4, 0.0)
    case('rb-step', 'boundary', head(0.4, 0.4, 1, 1, 1e-6, s1, s2), v(u1, v1, u2, v2), v(0.0, -0.5, -0.5, 0.0),
         v(0.4, 0.0, 0.4), v(0.0, -1.0, 0.0), h(0.1), domain(), domain())
# Clipped steps that end on the clipping boundary, where SMLib re-solves the
# step with LocalCurveIntersect on an offset surface: either surface clipping,
# on a u or a v side, forwards and backwards, flat, bumped and rational.
def box(u1, v1):
    return v(0.0, 0.0, u1, v1)


for s1, s2 in PAIRS:
    for r1, r2 in RADII:
        for y, sign in ((1.9, 1.0), (0.1, -1.0)):
            (u1, v1), (u2, v2) = consistent(r1, r2, y)
            point = v(r2, y, r1)
            tangent = v(0.0, sign, 0.0)
            for step in (0.2, 0.5):
                # surface 1 (the floor) clips on its v side
                case('rb-step', 'boundary', head(r1, r2, 1, 1, 1e-6, s1, s2), v(u1, v1, u2, v2),
                     v(0.0, 0.5 * sign, 0.5 * sign, 0.0), point, tangent, h(step * sign), domain(), domain())
                # surface 2 (the wall) clips first, on its u side
                low, high = (0.03, 0.97)
                d2 = box(high, 1.0) if sign > 0 else v(low, 0.0, 1.0, 1.0)
                case('rb-step', 'boundary', head(r1, r2, 1, 1, 1e-6, s1, s2), v(u1, v1, u2, v2),
                     v(0.0, 0.5 * sign, 0.5 * sign, 0.0), point, tangent, h(step * sign), domain(), d2)
                # the pair swapped: surface 1 (the wall) clips on its u side, surface 2 on its v side
                case('rb-step', 'boundary', head(r2, r1, 1, 1, 1e-6, s2, s1), v(u2, v2, u1, v1),
                     v(0.5 * sign, 0.0, 0.0, 0.5 * sign), point, tangent, h(step * sign), domain(), domain())
                case('rb-step', 'boundary', head(r2, r1, 1, 1, 1e-6, s2, s1), v(u2, v2, u1, v1),
                     v(0.5 * sign, 0.0, 0.0, 0.5 * sign), point, tangent, h(step * sign), domain(), box(1.0, 0.99))
        # a parameter tangent that misses the curve: the boundary solve finds no
        # meeting and the step goes on to the plane solve
        (u1, v1), (u2, v2) = consistent(r1, r2, 1.9)
        case('rb-step', 'boundary', head(r1, r2, 1, 1, 1e-6, s1, s2), v(u1, v1, u2, v2),
             v(0.4, 0.5, 0.5, 0.0), v(r2, 1.9, r1), v(0.0, 1.0, 0.0), h(0.3), domain(), domain())
# Domains past the natural ones, as SMLib's fillet allows. Where the other
# surface clips on its natural side the step is answered within the base
# surfaces; where the clipping side lies outside a base surface SMLib extends
# it (sm_ExtendSurface, not ported) and the port refuses.
(u1, v1), (u2, v2) = consistent(0.4, 0.4, 1.9)
case('rb-step', 'boundary', head(0.4, 0.4, 1, 1, 1e-6, *PAIRS[0]), v(u1, v1, u2, v2),
     v(0.0, 0.5, 0.5, 0.0), v(0.4, 1.9, 0.4), v(0.0, 1.0, 0.0), h(0.3), box(1.0, 1.05), domain())
for s1, s2 in PAIRS[:2]:
    for d1, d2 in ((box(1.0, 1.05), box(1.05, 1.0)), (box(1.0, 1.06), box(1.04, 1.0))):
        case('rb-step-extended', 'boundary', head(0.4, 0.4, 1, 1, 1e-6, s1, s2), v(u1, v1, u2, v2),
             v(0.0, 0.5, 0.5, 0.0), v(0.4, 1.9, 0.4), v(0.0, 1.0, 0.0), h(0.3), d1, d2)
# Parallel offset surfaces: the normals are parallel and the tangent comes from stepping.
for rc in (0.5, 0.45):
    for prev in (0, 1):
        for uv in ((0.3, 0.4, 0.5, 0.6), (0.5, 0.5, 0.5, 0.5), (0.1, 0.9, 0.9, 0.1)):
            args = [head(rc, 1.0 - rc, 1, 0, 1e-6, floor(0.0), ceiling(0.0)), v(*uv), h(1.0), domain(), domain(), prev]
            if prev:
                args.append(v(0.0, 1.0, 0.0))
            case('rb-point', 'tangent', *args)
        case('rb-surface', 'tangent', head(rc, 1.0 - rc, 1, 0, 1e-6, floor(0.0), ceiling(0.0)), v(0.3, 0.4, 0.5, 0.6))
# Degenerate: zero radii and equal surfaces.
for r1, r2 in ((0.0, 0.4), (0.4, 0.0), (0.0, 0.0), (1e-13, 0.4)):
    (u1, v1), (u2, v2) = consistent(0.4, 0.4, 0.9)
    case('rb-surface', 'singular', head(r1, r2, 1, 1, 1e-6, floor(0.0), wall(0.0)), v(u1, v1, u2, v2))
    case('rb-point', 'singular', head(r1, r2, 1, 1, 1e-6, floor(0.0), wall(0.0)), v(u1, v1, u2, v2), h(1.0),
         domain(), domain(), 0)
case('rb-surface', 'singular', head(0.4, 0.4, 1, 1, 1e-6, floor(0.0), floor(0.0)), v(0.5, 0.5, 0.5, 0.5))
case('rb-point', 'touching', head(0.4, 0.4, 1, 1, 1e-6, floor(0.0), floor(0.0)), v(0.5, 0.5, 0.5, 0.5), h(1.0),
     domain(), domain(), 0)


# An evaluator failure at the first guess (from review of this unit): constant
# patches have no normal and no nonsingular step-off, so SMLib's solver stops at
# the guess, "found" and unable to converge. The port refuses it.
def constant(x, y, z):
    return surface(1, 1, 2, 2, lambda s, t: (x, y, z))


for step in (0.01, -0.01):
    for pair in ((constant(0, 0, 0), constant(1, 0, 0)), (constant(0, 0, 0), constant(0, 0, 0)),
                 (floor(0.0), constant(0.4, 0.3, 0.2)), (constant(0.4, 0.3, 0.2), wall(0.0))):
        case('rb-step-failed', 'singular', head(0.4, 0.4, 1, 1, 1e-6, *pair), v(0.5, 0.5, 0.5, 0.5),
             v(1.0, 0.0, 1.0, 0.0), v(0.0, 0.0, 0.0), v(0.0, 1.0, 0.0), h(step), domain(), domain())


def main():
    oracle = sys.argv[1]
    raw = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if raw:
        raw.mkdir(parents=True, exist_ok=True)
    lines = []
    for index, (op, kind, inputs) in enumerate(CASES):
        # The evaluator-failure and extension steps are rb-step to SMLib and
        # pinned apart here.
        asked = 'rb-step' if op in ('rb-step-failed', 'rb-step-extended') else op
        request = f'{asked} {inputs}\n'
        run = subprocess.run(['prlimit', '--as=1073741824', '--cpu=10', 'timeout', '20', oracle],
                             input=request, capture_output=True, text=True)
        if raw:
            (raw / f'{index:04d}-{op}.txt').write_text(
                f'request: {request}exit: {run.returncode}\nstdout: {run.stdout}stderr: {run.stderr}')
        if run.returncode != 0:
            sys.exit(f'case {index} {op} {kind}: oracle exit {run.returncode}: {run.stderr.strip()}')
        response = run.stdout.split()
        assert response[0] == asked, (index, response)
        lines.append(f'{op} {kind} | {inputs} | {" ".join(response[1:])}')
    header = ('# catallax-trace rolling-ball gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
