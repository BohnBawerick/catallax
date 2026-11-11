#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace ruled-surface gate.

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
NO_Z = sys.float_info.max
NO_W = -sys.float_info.max


def h(x):
    return format(struct.unpack('<Q', struct.pack('<d', float(x)))[0], '016x')


def v(*xs):
    return ' '.join(h(x) for x in xs)


CASES = []


def case(op, kind, *tokens):
    assert kind in ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
    CASES.append((op, kind, ' '.join(str(t) for t in tokens)))


rng = random.Random(20261115)


def points_for(count, dim, rational):
    vals = []
    for _ in range(count):
        x, y, z = (rng.uniform(-4, 4) for _ in range(3))
        if rational:
            w = rng.uniform(.5, 2)
            vals += [x * w, y * w, z * w if dim == 3 else NO_Z, w]
        else:
            vals += [x, y, z if dim == 3 else NO_Z, NO_W]
    return vals


def curve(dim, degree, interior, rational, lo=0.0, hi=3.0):
    knots = [lo] * (degree + 1) + [lo + (hi - lo) * k / 3.0 for k in interior] + [hi] * (degree + 1)
    count = len(knots) - degree - 1
    return f'{dim} {degree} {count} ' + v(*knots) + ' ' + v(*points_for(count, dim, rational))


INTERIORS = {1: [[]], 2: [[], [1.0], [0.7, 2.2]], 3: [[], [1.5], [0.6, 1.4, 2.4]]}

# Ruled surfaces.
for direction in (0, 1):
    for da in (1, 2, 3):
        for db in (1, 2, 3):
            for rational_a, rational_b in ((False, False), (True, True), (False, True)):
                for dim_a, dim_b in ((3, 3), (2, 2), (2, 3)):
                    ia = rng.choice(INTERIORS[da])
                    ib = rng.choice(INTERIORS[db])
                    case('ruled', 'regular', direction, curve(dim_a, da, ia, rational_a), curve(dim_b, db, ib, rational_b))
    # Identical knots, so nothing to refine; same curve twice.
    for degree in (1, 2, 3):
        a = curve(3, degree, INTERIORS[degree][-1], False)
        case('ruled', 'tangent', direction, a, a)
        case('ruled', 'tangent', direction, a, curve(3, degree, INTERIORS[degree][-1], False))
    # Different parameter ranges are scaled to the unit interval.
    for lo, hi in ((1.0, 5.0), (-2.0, 7.5), (0.0, 1.0), (0.0, 3.0 + 1e-13), (1e-13, 3.0)):
        case('ruled', 'boundary' if lo or hi != 3.0 else 'regular', direction,
             curve(3, 2, [1.0], False), curve(3, 3, [1.5], False, lo, hi))
    case('ruled', 'near-singular', direction, curve(3, 2, [1.0], False), curve(3, 2, [1.0], False, 1e-9, 3.0))
    case('ruled', 'near-singular', direction, curve(3, 2, [1.0], False), curve(3, 2, [1.0 + 1e-10], False))
    case('ruled', 'touching', direction, curve(3, 2, [1.0], False), curve(3, 2, [1.0 + 1e-7], False))
    case('ruled', 'touching', direction, curve(3, 3, [1.0, 1.0], False), curve(3, 3, [1.0], False))
    # Degenerate curves: all control points coincide.
    degenerate = f'3 2 3 ' + v(0, 0, 0, 3, 3, 3) + ' ' + v(*([1.0, 2.0, 3.0, NO_W] * 3))
    case('ruled', 'singular', direction, degenerate, curve(3, 2, [], False))
    case('ruled', 'singular', direction, degenerate, degenerate)

# Reparameterising a surface.


def surface(du, dv, interior_u, interior_v, rational, lo_u=0.0, hi_u=3.0, lo_v=0.0, hi_v=2.0):
    ku = [lo_u] * (du + 1) + [lo_u + (hi_u - lo_u) * k / 3.0 for k in interior_u] + [hi_u] * (du + 1)
    kv = [lo_v] * (dv + 1) + [lo_v + (hi_v - lo_v) * k / 2.0 for k in interior_v] + [hi_v] * (dv + 1)
    nu, nv = len(ku) - du - 1, len(kv) - dv - 1
    return f'{du} {dv} {nu} {nv} {int(rational)} ' + v(*ku) + ' ' + v(*kv) + ' ' + v(*points_for(nu * nv, 3, rational))


for du, dv in ((1, 1), (2, 1), (1, 3), (3, 3), (2, 2)):
    for rational in (False, True):
        s = surface(du, dv, rng.choice(INTERIORS[du]), [0.5 * k for k in rng.choice(INTERIORS[dv] if dv < 3 else [[], [1.0]])][:0] or ([] if dv == 1 else [1.0]), rational)
        for dom in ((0.0, 1.0, 0.0, 1.0), (-2.0, 5.0, 10.0, 11.0), (0.0, 3.0, 0.0, 2.0), (1.0, 4.0, 0.0, 2.0),
                    (0.0, 3.0, 5.0, 9.0), (0.0, 3.0 + 1e-13, 0.0, 2.0), (3.0, 0.0, 2.0, 0.0), (1e-9, 3.0, 0.0, 2.0)):
            kind = 'boundary' if dom[:2] == (0.0, 3.0) or dom[2:] == (0.0, 2.0) else 'regular'
            case('reparam', kind, s, v(*dom))
        case('reparam', 'near-singular', s, v(0.5, 0.5 + 1e-9, 0.0, 1.0))
        case('reparam-single', 'singular', s, v(0.5, 0.5, 0.0, 1.0))
        case('reparam-single', 'singular', s, v(0.0, 1.0, 2.0, 2.0))
        case('reparam-single', 'singular', s, v(0.0, 3.0, 0.0, 0.0))


# Surfaces of revolution.
def generator(dim, degree, interior, rational, offset=(3.0, 0.0, 0.0), poles=False):
    knots = [0.0] * (degree + 1) + [3.0 * k / 3.0 for k in interior] + [3.0] * (degree + 1)
    count = len(knots) - degree - 1
    vals = []
    for i in range(count):
        x, y, z = (rng.uniform(-2, 2) + o for o, _ in zip(offset, range(3)))
        if poles and i in (0, count - 1):
            x, y = offset[0] * 0, offset[1] * 0
            x, y = 0.0, 0.0
        if rational:
            w = rng.uniform(.5, 2)
            vals += [x * w, y * w, z * w if dim == 3 else NO_Z, w]
        else:
            vals += [x, y, z if dim == 3 else NO_Z, NO_W]
    return f'{dim} {degree} {count} ' + v(*knots) + ' ' + v(*vals), count


AXES = [([0, 0, 0], [0, 0, 1]), ([1, -2, 0.5], [0, 0, 2]), ([0, 0, 0], [1, 1, 0.3]), ([0.5, 0.5, 0.5], [-1, 0.2, 0.7]),
        ([0, 0, 0], [1, 0, 0])]
for degree in (1, 2, 3):
    for rational in (False, True):
        for dim in (3, 2):
            for axis in AXES[:4]:
                for angle in (30.0, 90.0, 200.0, 360.0):
                    g_, _ = generator(dim, degree, rng.choice(INTERIORS[degree]), rational)
                    case('revolve', 'regular', g_, v(*axis[0], *axis[1], angle), rng.choice((0, 1)))
    for swap in (0, 1):
        g_, _ = generator(3, degree, INTERIORS[degree][-1], False)
        case('revolve', 'regular', g_, v(0, 0, 0, 0, 0, 1, 90.0), swap)
        case('revolve', 'regular', g_, v(0, 0, 0, 0, 0, 1, 360.0), swap)
    # Angles at the thresholds of the span count, and beyond the limits.
    for angle in (1e-6, 45.0, 90.9999, 91.0, 180.9999, 181.0, 270.9999, 271.0, 359.99, 359.995, 360.0, 0.0, -90.0, -360.0,
                  360.0001, 361.0, 720.0):
        kind = 'singular' if abs(angle) > 360 else 'boundary' if angle in (91.0, 181.0, 271.0, 360.0, 0.0, -360.0, 359.995) else 'regular'
        g_, _ = generator(3, degree, INTERIORS[degree][-1], False)
        case('revolve', kind, g_, v(0, 0, 0, 0, 0, 1, angle), 0)
    # Poles: end points on the axis, and a generator lying on the axis.
    g_, _ = generator(3, degree, INTERIORS[degree][-1], False, poles=True)
    case('revolve', 'touching', g_, v(0, 0, 0, 0, 0, 1, 360.0), 0)
    case('revolve', 'touching', g_, v(0, 0, 0, 0, 0, 1, 90.0), 1)
    on_axis = f'3 {degree} {degree + 1} ' + v(*([0.0] * (degree + 1) + [3.0] * (degree + 1))) + ' ' + v(*sum(([0.0, 0.0, float(i), NO_W] for i in range(degree + 1)), []))
    case('revolve', 'singular', on_axis, v(0, 0, 0, 0, 0, 1, 90.0), 0)
    # A point 1e-13 and 1e-11 from the axis.
    for off in (1e-13, 1e-11, 1e-9):
        c = f'3 {degree} {degree + 1} ' + v(*([0.0] * (degree + 1) + [3.0] * (degree + 1))) + ' ' + v(*sum(([off if i == 0 else 2.0, 0.0, float(i), NO_W] for i in range(degree + 1)), []))
        case('revolve', 'near-singular', c, v(0, 0, 0, 0, 0, 1, 90.0), 0)
    # Zero and tiny axes.
    g_, _ = generator(3, degree, INTERIORS[degree][-1], False)
    case('revolve', 'singular', g_, v(0, 0, 0, 0, 0, 0, 90.0), 0)
    case('revolve', 'near-singular', g_, v(0, 0, 0, 0, 0, 1e-6, 90.0), 0)
    case('revolve', 'singular', g_, v(0, 0, 0, 1e-13, 0, 0, 90.0), 0)
# Weights outside the range SMLib accepts, and equal weights.
for weight in (0.0005, 0.001, 100.0, 100.5, 2.0):
    g_, count = generator(3, 2, [1.0], False)
    toks = g_.split(' ')
    # rewrite as rational with a constant weight
    pts = []
    for _ in range(4):
        x, y, z = (rng.uniform(-2, 2) + o for o in (3.0, 0.0, 0.0))
        pts += [x * weight, y * weight, z * weight, weight]
    c = f'3 2 4 ' + v(0, 0, 0, 1.5, 3, 3, 3) + ' ' + v(*pts)
    case('revolve', 'boundary' if weight in (0.001, 100.0) else 'singular' if weight in (0.0005, 100.5) else 'tangent',
         c, v(0, 0, 0, 0, 0, 1, 90.0), 0)

# Unclamped generators: SMLib clamps them through N_IgesClampSrf.
for degree, knots in ((2, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]), (2, [0.0, 0.0, 1.0, 1.5, 2.0, 3.0, 3.0, 3.0]),
                      (3, [0.0, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 3.0, 3.0, 3.0]), (1, [0.0, 1.0, 2.0, 3.0, 4.0])):
    count = len(knots) - degree - 1
    c = f'3 {degree} {count} ' + v(*knots) + ' ' + v(*sum(([rng.uniform(1, 3), rng.uniform(-1, 1), rng.uniform(-1, 1), NO_W] for _ in range(count)), []))
    for swap in (0, 1):
        case('revolve-unclamped', 'singular', c, v(0, 0, 0, 0, 0, 1, 90.0), swap)


def main():
    oracle = sys.argv[1]
    raw = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if raw:
        raw.mkdir(parents=True, exist_ok=True)
    lines = []
    for index, (op, kind, inputs) in enumerate(CASES):
        request_op = {'reparam-single': 'reparam', 'revolve-unclamped': 'revolve'}.get(op, op)
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
        if op == 'revolve-unclamped' and response[1] != 'i1000':
            op = 'revolve'
        lines.append(f'{op} {kind} | {inputs} | {" ".join(response[1:])}')
    header = ('# catallax-trace ruled-surface gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
