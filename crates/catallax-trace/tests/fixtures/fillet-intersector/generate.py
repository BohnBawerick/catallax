#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace fillet-intersector gate.

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


rng = random.Random(20261118)


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



def consistent(r1, r2, y):
    """uv on each surface of the offset-surface crossing at height y (flat patches)."""
    return (r2 / 3.0, y / 2.0), (y / 2.0, r1 / 3.0)


def stretched_floor():
    # v runs as the square root of y, so the parameter speed changes fivefold
    return surface(3, 3, 4, 4, lambda s, t: (3 * s, 2 * t * t, 0.0))


def stretched_wall():
    return surface(3, 3, 4, 4, lambda s, t: (0.0, 2 * s * s, 3 * t))


def full(r1, r2, o1, o2, tol, s1, s2, d1=None, d2=None):
    return f'{head(r1, r2, o1, o2, tol, s1, s2)} {d1 or domain()} {d2 or domain()}'


def tail(approx, angle=math.radians(20.0), direction=1.0, through=0):
    return f'{h(approx)} {h(angle)} {h(direction)} {through}'


def consistent_stretched(r1, r2, y):
    s = math.sqrt(y / 2.0)
    return (r2 / 3.0, s), (s, r1 / 3.0)


PAIRS = [(floor(0.0), wall(0.0)), (floor(0.06), wall(0.06)), (floor(0.08, True), wall(0.08, degree=2)),
         (stretched_floor(), stretched_wall())]
RADII = [(0.4, 0.4), (0.3, 0.5)]
SUB = v(0.1, 0.1, 0.9, 0.9)


def point(pair_index, r1, r2, y):
    return consistent_stretched(r1, r2, y) if pair_index == 3 else consistent(r1, r2, y)


# Construction: the tolerance is a thousandth of the first offset surface's size.
for index, (s1, s2) in enumerate(PAIRS):
    for r1, r2 in RADII:
        for o1, o2 in ((1, 1), (0, 1)):
            case('fi-ctor', 'regular', full(r1, r2, o1, o2, 1e-6, s1, s2))
            case('fi-ctor', 'regular', full(r1, r2, o1, o2, 1e-6, s1, s2, SUB, SUB))
case('fi-ctor', 'singular', full(0.0, 0.4, 1, 1, 1e-6, floor(0.0), wall(0.0)))
case('fi-ctor', 'singular', full(0.4, 0.4, 1, 1, 1e-6, floor(0.0), floor(0.0)))

STEPS = (0.02, 0.2, 0.5)
for index, (s1, s2) in enumerate(PAIRS):
    for r1, r2 in RADII:
        for y in (0.3, 1.5):
            (u1, v1), (u2, v2) = point(index, r1, r2, y)
            uvs = v(u1, v1, u2, v2)
            base = full(r1, r2, 1, 1, 1e-6, s1, s2)
            for direction in (1.0, -1.0):
                case('fi-values', 'regular', base, tail(1e-3, direction=direction), uvs, 0)
            case('fi-values', 'regular', base, tail(1e-3), uvs, 1, v(0.0, 1.0, 0.0))
            # delta: a chord to a point a little further along the curve
            (du1, dv1), (du2, dv2) = point(index, r1, r2, y + 0.05)
            for param, delta_param in ((0.0, 0.05), (0.05, 0.0)):
                case('fi-delta', 'regular', base, tail(1e-3), uvs, v(du1, dv1, du2, dv2), h(param), h(delta_param))
            case('fi-delta', 'singular', base, tail(1e-3), uvs, uvs, h(0.0), h(0.0))
            # the plane normal to the curve at the point: the guess is on it, or off it
            origin = v(r2, y, r1)
            case('fi-law', 'regular', base, tail(1e-3), uvs, origin, v(0.0, 1.0, 0.0))
            case('fi-law', 'regular', base, tail(1e-3), v(u1 + 0.05, v1 - 0.03, u2 - 0.02, v2 + 0.04), origin, v(0.0, 1.0, 0.0))
            # next point and span, from the point
            for step in STEPS:
                for direction in (1.0, -1.0):
                    for through in (0, 1):
                        case('fi-next', 'regular', base, tail(1e-3, direction=direction, through=through), uvs, h(0.0),
                             h(step))
            # span between this point and one further along, with and without a middle point
            for far in (0.05, 0.3):
                (nu1, nv1), (nu2, nv2) = point(index, r1, r2, y + far)
                (mu1, mv1), (mu2, mv2) = point(index, r1, r2, y + far / 2)
                for approx in (1e-3, 1e-12):
                    for angle in (math.radians(20.0), 1e-4):
                        case('fi-span', 'regular', base, tail(approx, angle), uvs, h(0.0), v(nu1, nv1, nu2, nv2), h(far), 0)
                case('fi-span', 'regular', base, tail(1e-3), uvs, h(0.0), v(nu1, nv1, nu2, nv2), h(far), 1,
                     v(mu1, mv1, mu2, mv2))
                case('fi-span', 'regular', base, tail(1e-12), uvs, h(0.0), v(nu1, nv1, nu2, nv2), h(far), 1,
                     v(mu1 + 0.1, mv1, mu2, mv2))

# Steps near the ends of the domain are clipped by it.
for s1, s2 in PAIRS[:2]:
    for y, direction in ((1.95, 1.0), (0.05, -1.0)):
        (u1, v1), (u2, v2) = consistent(0.4, 0.4, y)
        base = full(0.4, 0.4, 1, 1, 1e-6, s1, s2)
        for step in (0.2, 0.05, 0.5):
            case('fi-next', 'boundary', base, tail(1e-3, direction=direction), v(u1, v1, u2, v2), h(0.0), h(step))
        case('fi-law', 'boundary', base, tail(1e-3), v(u1, v1, u2, v2), v(0.4, y + 0.3, 0.4), v(0.0, 1.0, 0.0))
case('fi-law', 'singular', full(0.4, 0.4, 1, 1, 1e-6, floor(0.0), wall(0.0)), tail(1e-3),
     v(*consistent(0.4, 0.4, 0.9)[0], *consistent(0.4, 0.4, 0.9)[1]), v(0.4, 0.9, 0.4), v(0.0, 0.0, 0.0))

# Parallel offset surfaces: the normals are parallel; the tangent comes from stepping.
for rc in (0.5, 0.45):
    base = full(rc, 1.0 - rc, 1, 0, 1e-6, floor(0.0), ceiling(0.0))
    for uv in ((0.3, 0.4, 0.5, 0.6), (0.5, 0.5, 0.5, 0.5)):
        case('fi-values', 'tangent', base, tail(1e-3), v(*uv), 0)
        case('fi-values', 'tangent', base, tail(1e-3), v(*uv), 1, v(0.0, 1.0, 0.0))
        case('fi-next', 'tangent', base, tail(1e-3), v(*uv), h(0.0), h(0.1))
        case('fi-span', 'tangent', base, tail(1e-3), v(*uv), h(0.0), v(0.35, 0.45, 0.55, 0.65), h(0.1), 0)
# Degenerate: zero radii and equal surfaces.
for r1, r2 in ((0.0, 0.4), (0.4, 0.0), (1e-13, 0.4)):
    (u1, v1), (u2, v2) = consistent(0.4, 0.4, 0.9)
    base = full(r1, r2, 1, 1, 1e-6, floor(0.0), wall(0.0))
    case('fi-values', 'singular', base, tail(1e-3), v(u1, v1, u2, v2), 0)
    case('fi-next', 'singular', base, tail(1e-3), v(u1, v1, u2, v2), h(0.0), h(0.1))
    case('fi-span', 'singular', base, tail(1e-3), v(u1, v1, u2, v2), h(0.0), v(u1, v1 + 0.05, u2 + 0.05, v2), h(0.1), 0)
same = full(0.4, 0.4, 1, 1, 1e-6, floor(0.0), floor(0.0))
case('fi-values', 'touching', same, tail(1e-3), v(0.5, 0.5, 0.5, 0.5), 0)
case('fi-next', 'touching', same, tail(1e-3), v(0.5, 0.5, 0.5, 0.5), h(0.0), h(0.1))
case('fi-span', 'touching', same, tail(1e-3), v(0.5, 0.5, 0.5, 0.5), h(0.0), v(0.5, 0.6, 0.5, 0.6), h(0.1), 0)


# The parameter speed changes more than threefold along a cubic: v runs as the
# cube root of y. A step across it changes a tangent's length, which a
# through-point clipped trace repairs from a short step.
def cubic_floor():
    return surface(3, 3, 4, 4, lambda s, t: (3 * s, 2.0 if round(3 * t) == 3 else 0.0, 0.0))


def cubic_wall():
    return surface(3, 3, 4, 4, lambda s, t: (0.0, 2.0 if round(3 * s) == 3 else 0.0, 3 * t))


for y in (0.05, 0.1, 0.3):
    t = (y / 2.0) ** (1.0 / 3.0)
    base = full(0.4, 0.4, 1, 1, 1e-6, cubic_floor(), cubic_wall())
    uvs = v(0.4 / 3.0, t, t, 0.4 / 3.0)
    for direction in (1.0, -1.0):
        for step in (0.05, 0.1, 0.2):
            for through in (0, 1):
                case('fi-next', 'regular', base, tail(1e-3, direction=direction, through=through), uvs, h(0.0), h(step))


# The same more-than-threefold change with a curved intersection. The floor's
# y runs as the cube of its parameter; the wall bends, x = 1 + v^2. The offset
# surfaces meet on a parabola, so the two chords a short step and a full step
# give are not collinear; `ComputeNextPoint`'s through-point repair evaluates
# its point values at the full step's UVs, as SMLib does, and these cases pin
# both the repaired current tangent and the next point where the directions
# differ. The start is the offset intersection at a wall parameter: the wall
# offsets by (1, -v)/sqrt(1 + v^2) and the floor's z offset plane cuts it at
# u = r1/3.
def bend_floor():
    # (3s, 2t^3, 0): degree 1 in s, degree 3 in t, four control rows in y.
    ku, kv = [0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    control = [(3 * s, y, 0.0) for s in (0.0, 1.0) for y in (0.0, 0.0, 0.0, 2.0)]
    return '1 3 2 4 0 ' + v(*ku) + ' ' + v(*kv) + ' ' + ' '.join(v(x, y, z, 1.0) for x, y, z in control)


def bend_wall():
    # (1 + v^2, 2v, 3u): degree 1 in u, degree 2 in v.
    ku, kv = [0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    control = [(x, y, 3 * u) for u in (0.0, 1.0) for x, y in ((1.0, 0.0), (1.0, 1.0), (2.0, 2.0))]
    return '1 2 2 3 0 ' + v(*ku) + ' ' + v(*kv) + ' ' + ' '.join(v(x, y, z, 1.0) for x, y, z in control)


def bend_uvs(bend, r1=0.4, r2=0.4):
    root = math.sqrt(1.0 + bend * bend)
    x, y = 1.0 + bend * bend - r2 / root, 2.0 * bend + r2 * bend / root
    return v(x / 3.0, (y / 2.0) ** (1.0 / 3.0), r1 / 3.0, bend)


bend_base = full(0.4, 0.4, 1, 1, 1e-6, bend_floor(), bend_wall())
for bend, step, direction, through in ((0.005, 0.01, 1.0, 1), (0.005, 0.01, 1.0, 0),
                                       (0.01, 0.2, -1.0, 1), (0.02, 0.3, -1.0, 1),
                                       (0.01, 0.02, 1.0, 1)):
    case('fi-next', 'regular', bend_base, tail(1e-3, direction=direction, through=through),
         bend_uvs(bend), h(0.0), h(step))


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
    header = ('# catallax-trace fillet-intersector gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
