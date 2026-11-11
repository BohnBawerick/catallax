#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace curve-construction gate.

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


rng = random.Random(20261111)


def r3(scale=10.0):
    return [rng.uniform(-scale, scale) for _ in range(3)]


def knots_for(degree, count, interior=None):
    """Clamped knot vector; `interior` lists the count - degree - 1 interior knots."""
    if interior is None:
        interior = sorted(round(rng.uniform(0.2, 2.8), 3) for _ in range(count - degree - 1))
        interior = sorted(set(interior))
        while len(interior) < count - degree - 1:
            interior.append(round(interior[-1] + 0.1, 3) if interior else 1.0)
    assert len(interior) == count - degree - 1
    return [0.0] * (degree + 1) + list(interior) + [3.0] * (degree + 1)


def curve(dim, degree, knots, rational, points=None):
    """`dim degree count knots... (x y z w)...`, stored values."""
    count = len(knots) - degree - 1
    if points is None:
        points = [r3(4) for _ in range(count)]
    values = []
    for x, y, z in points:
        if rational:
            w = rng.uniform(.5, 2)
            values += [x * w, y * w, z * w if dim == 3 else NO_Z, w]
        else:
            values += [x, y, z if dim == 3 else NO_Z, NO_W]
    return f'{dim} {degree} {count} ' + v(*knots) + ' ' + v(*values)


# Lines.
for dim in (2, 3):
    for _ in range(4):
        case('line', 'regular', dim, v(*r3(), *r3()))
    case('line', 'tangent', dim, v(0, 0, 0, 5, 0, 0))
    case('line', 'tangent', dim, v(1, 2, 3, 1, 2, 7))
    case('line', 'touching', dim, v(0, 0, 0, 1e-12, 0, 0))
    case('line', 'near-singular', dim, v(0, 0, 0, 0.9e-12, 0, 0))
    case('line', 'near-singular', dim, v(1, 1, 1, 1 + 1e-12, 1, 1))
    case('line', 'boundary', dim, v(1e15, 0, 0, 1e15 + 0.125, 0, 0))
    case('line', 'singular', dim, v(2, 3, 4, 2, 3, 4))
case('line', 'singular', 1, v(0, 0, 0, 1, 1, 1))
case('line', 'singular', 4, v(0, 0, 0, 1, 1, 1))
case('line', 'singular', 3, v(0, 0, 0, float('nan'), 1, 1))

# Degenerate and point curves.
for dim in (3,):
    for _ in range(3):
        case('point', 'regular', dim, v(*r3()))
    case('point', 'singular', dim, v(0, 0, 0))
case('point', 'singular', 1, v(0, 0, 0))
case('point', 'singular', 4, v(0, 0, 0))
for _ in range(2):
    case('point-planar', 'regular', 2, v(*r3()))

# Knot insertion.
for degree in (1, 2, 3, 4):
    for rational in (False, True):
        for dim in (2, 3):
            count = degree + 3
            kn = knots_for(degree, count)
            interior = kn[degree + 1:-degree - 1]
            for times in range(0, degree + 1):
                u = round(rng.uniform(0.05, 2.95), 4)
                if u in kn:
                    u += 0.0123
                case('insert', 'regular', curve(dim, degree, kn, rational), h(u), times)
            # An existing interior knot: tangent when the multiplicity reaches the degree.
            for knot in interior[:2]:
                mult = kn.count(knot)
                for times in range(1, degree - mult + 1):
                    kind = 'tangent' if times == degree - mult else 'touching'
                    case('insert', kind, curve(dim, degree, kn, rational), h(knot), times)
                case('insert', 'singular', curve(dim, degree, kn, rational), h(knot), degree - mult + 1)
            # A parameter a hair from an existing knot, which the fixup merges.
            if interior:
                for eps in (1e-13, -1e-13, 5e-13, 3e-12):
                    case('insert', 'near-singular', curve(dim, degree, kn, rational), h(interior[0] + eps), 1)
            # End parameters.
            case('insert', 'boundary', curve(dim, degree, kn, rational), h(0.0), 1)
            case('insert', 'boundary', curve(dim, degree, kn, rational), h(3.0), 1)
            case('insert', 'near-singular', curve(dim, degree, kn, rational), h(1e-13), 1)
            case('insert', 'near-singular', curve(dim, degree, kn, rational), h(3.0 - 1e-13), 1)
            case('insert', 'boundary', curve(dim, degree, kn, rational), h(-1.0), 1)
            # Too many insertions at a fresh parameter.
            case('insert', 'singular', curve(dim, degree, kn, rational), h(1.2345), degree + 1)
# A curve whose control points all coincide.
for degree in (2, 3):
    kn = knots_for(degree, degree + 2)
    case('insert', 'singular', curve(3, degree, kn, False, [[1., 2., 3.]] * (degree + 2)), h(1.5), 1)
    case('insert', 'singular', curve(3, degree, kn, True, [[1., 2., 3.]] * (degree + 2)), h(1.5), 2)
# Bezier curves.
for degree in (1, 2, 3, 4):
    case('insert', 'regular', curve(3, degree, [0.0] * (degree + 1) + [3.0] * (degree + 1), False), h(1.0), 1)
    case('insert', 'regular', curve(3, degree, [0.0] * (degree + 1) + [3.0] * (degree + 1), True), h(2.0), degree)

# Refinement.
for degree in (1, 2, 3, 4):
    for rational in (False, True):
        for dim in (2, 3):
            count = degree + 3
            kn = knots_for(degree, count)
            interior = kn[degree + 1:-degree - 1]
            for size in (1, 2, 4):
                xs = sorted(round(rng.uniform(0.05, 2.95), 4) for _ in range(size))
                case('refine', 'regular', curve(dim, degree, kn, rational), len(xs), v(*xs))
            if interior:
                case('refine', 'tangent', curve(dim, degree, kn, rational), 2,
                     v(interior[0], round(interior[0] + 0.5, 4)))
                case('refine', 'touching', curve(dim, degree, kn, rational), 2,
                     v(interior[0] - 1e-13, interior[0] + 1e-13))
            case('refine', 'regular', curve(dim, degree, kn, rational), 3, v(0.7, 0.7, 1.9))
            case('refine', 'boundary', curve(dim, degree, kn, rational), 2, v(0.0, 1.0))
            case('refine', 'boundary', curve(dim, degree, kn, rational), 2, v(1.0, 3.0))
            case('refine', 'boundary', curve(dim, degree, kn, rational), 1, v(-0.5))
            case('refine', 'near-singular', curve(dim, degree, kn, rational), 2, v(1e-13, 1.0))
            case('refine', 'singular', curve(dim, degree, kn, rational), 0)
            case('refine-unordered', 'singular', curve(dim, degree, kn, rational), 2, v(2.0, 1.0))
for degree in (2, 3):
    kn = knots_for(degree, degree + 2)
    case('refine', 'singular', curve(3, degree, kn, False), degree + 2, v(*([1.5] * (degree + 2))))

# Reparameterisation.
for degree in (1, 2, 3):
    for rational in (False, True):
        for dim in (2, 3):
            kn = knots_for(degree, degree + 3)
            for lo, hi in ((-1.0, 2.5), (10.0, 11.0), (0.0, 1.0), (-5.0, -4.0)):
                case('reparam', 'regular', curve(dim, degree, kn, rational), v(lo, hi))
            case('reparam', 'boundary', curve(dim, degree, kn, rational), v(0.0, 3.0))
            case('reparam', 'boundary', curve(dim, degree, kn, rational), v(1e-13, 3.0 + 1e-13))
            case('reparam', 'boundary', curve(dim, degree, kn, rational), v(0.0, 3.0 + 1e-9))
            case('reparam', 'boundary', curve(dim, degree, kn, rational), v(1e-9, 3.0))
            case('reparam', 'near-singular', curve(dim, degree, kn, rational), v(0.5, 0.5 + 1e-9))
            case('reparam', 'regular', curve(dim, degree, kn, rational), v(2.5, -1.0))
            case('reparam-nonincreasing', 'singular', curve(dim, degree, kn, rational), v(0.5, 0.5))
            case('reparam-nonincreasing', 'singular', curve(dim, degree, kn, rational), v(0.0, 0.0))


def main():
    oracle = sys.argv[1]
    raw = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if raw:
        raw.mkdir(parents=True, exist_ok=True)
    lines = []
    for index, (op, kind, inputs) in enumerate(CASES):
        request_op = {'point-planar': 'point', 'refine-unordered': 'refine',
                      'reparam-nonincreasing': 'reparam'}.get(op, op)
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
    header = ('# catallax-trace curve-construction gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
