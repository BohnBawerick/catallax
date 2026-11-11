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


def main():
    oracle = sys.argv[1]
    raw = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if raw:
        raw.mkdir(parents=True, exist_ok=True)
    lines = []
    for index, (op, kind, inputs) in enumerate(CASES):
        request_op = {'reparam-single': 'reparam'}.get(op, op)
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
    header = ('# catallax-trace ruled-surface gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
