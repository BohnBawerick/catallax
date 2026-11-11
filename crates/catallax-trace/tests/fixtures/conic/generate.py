#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace conic gate.

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
KINDS = ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
ELLIPSE, HYPERBOLA, PARABOLA, LINE = 0, 1, 2, 3
CASES = []


def h(x):
    return format(struct.unpack('<Q', struct.pack('<d', float(x)))[0], '016x')


def v(*xs):
    return ' '.join(h(x) for x in xs)


def case(op, kind, text):
    assert kind in KINDS
    CASES.append((op, kind, text))


def poly(op, kind, coefficients, zero=1e-12):
    case(op, kind, v(*coefficients, zero))


def expand(roots, lead=1.0):
    """Coefficients, constant first, of lead * prod(x - r)."""
    c = [lead]
    for r in roots:
        c = [0.0] + c
        for i in range(len(c) - 1):
            c[i] -= r * c[i + 1]
    return c


rng = random.Random(20261111)
# Cubics: three roots, one, a double and a triple root, a negligible leading
# coefficient, and random ones.
for roots in ((1, 2, 3), (-0.5, 0.25, 4), (0.1, 0.1, 2), (2, 2, 2), (0, 0, 0)):
    poly('cubic', 'singular' if len(set(roots)) < 3 else 'regular', expand(roots))
poly('cubic', 'regular', [1, 0, 1, 1])
poly('cubic', 'regular', [-2, 1, 0, 3])
poly('cubic', 'near-singular', [1, -3, 2, 1e-13])
poly('cubic', 'boundary', [0, 0, 0, 0])
for _ in range(8):
    poly('cubic', 'regular', [rng.uniform(-5, 5) for _ in range(4)])
# Quartics: four, two and no real roots, double roots, the r = 0 branch, a
# negligible leading coefficient, and random ones.
for roots in ((1, 2, 3, 4), (-1, -0.5, 0.5, 1), (0.3, 0.3, -2, 5), (1, 1, 1, 1), (0, 1, -1, 2)):
    poly('quartic', 'singular' if len(set(roots)) < 4 else 'regular', expand(roots))
poly('quartic', 'regular', [1, 0, 0, 0, 1])
poly('quartic', 'regular', [-1, 0, 0, 0, 1])
poly('quartic', 'regular', [0, 3, -2, 0, 1])
poly('quartic', 'regular', [4, 0, -5, 0, 1])
poly('quartic', 'near-singular', [1, -3, 2, 1, 1e-13])
poly('quartic', 'boundary', [0, 0, 0, 0, 0])
for _ in range(8):
    poly('quartic', 'regular', [rng.uniform(-5, 5) for _ in range(5)])


def conic(kind, a, b, *moves):
    return f'{kind} {v(a, b)} {len(moves)} ' + ' '.join(v(*m) for m in moves)


def pair(kind, one, two, tolerance=1e-4, op='conic'):
    case(op, kind, f'{one} {two} {v(tolerance)}')


# Dupin-like pairs: ellipses turned against each other (four meetings) and
# aligned (none, or tangent), ellipse and hyperbola, hyperbola branches,
# lines, a parabola, and moved origins.
E1 = conic(ELLIPSE, 2.0, 1.0)
pair('regular', E1, conic(ELLIPSE, 2.0, 1.0, (0.7, 0, 0)))
pair('regular', E1, conic(ELLIPSE, 1.5, 1.2, (0.3, 0, 0)))
pair('touching', E1, conic(ELLIPSE, 2.0, 0.5))
pair('tangent', E1, conic(ELLIPSE, 2.0, 1.0))
pair('regular', E1, conic(ELLIPSE, 1.0, 0.5))
pair('regular', E1, conic(HYPERBOLA, 1.0, 1.0))
pair('regular', E1, conic(HYPERBOLA, 1.0, 1.0, (math.pi, 0, 0)))
pair('regular', conic(HYPERBOLA, 1.0, 2.0), conic(HYPERBOLA, 1.0, 1.0, (0.4, 0, 0)))
pair('regular', conic(HYPERBOLA, 1.0, 2.0), conic(HYPERBOLA, 2.0, 1.0, (math.pi / 2, 0, 0)))
pair('regular', conic(LINE, 0.5, 0.0), E1)
pair('regular', conic(LINE, 0.5, 0.0, (math.pi / 2, 0, 0)), E1)
pair('touching', conic(LINE, 1.0, 0.0), E1)
pair('singular', conic(LINE, 0.5, 0.0), conic(LINE, 0.5, 0.0, (math.pi, 0, 0)))
pair('regular', conic(LINE, 0.5, 0.0), conic(LINE, 0.7, 0.0, (math.pi / 2, 0, 0)))
pair('regular', conic(PARABOLA, 0.5, 0.0), E1)
# Exceptions (`conicx`): the resolvent cubic's cube roots differ by one unit
# between the C library and the port's `libm`, and these meetings magnify that
# past the 16-unit gate; see README.md.
pair('regular', E1, conic(ELLIPSE, 1.0, 1.0, (0.0, 0.5, 0.25)), op='conicx')
pair('regular', conic(ELLIPSE, 2.0, 1.0, (0.2, 0.1, -0.3)), conic(HYPERBOLA, 0.8, 1.1, (-0.5, 0.2, 0.1)))
pair('near-singular', E1, conic(ELLIPSE, 2.0, 1.0, (1e-9, 0, 0)))
pair('boundary', E1, conic(ELLIPSE, 2.0, 1.0, (math.pi / 2, 0, 0)))
for n in range(8):
    one = conic(rng.choice((ELLIPSE, HYPERBOLA)), rng.uniform(0.5, 3), rng.uniform(0.5, 3),
                (rng.uniform(-3, 3), rng.uniform(-0.5, 0.5), rng.uniform(-0.5, 0.5)))
    two = conic(rng.choice((ELLIPSE, HYPERBOLA, LINE)), rng.uniform(0.5, 3), rng.uniform(0.5, 3),
                (rng.uniform(-3, 3), rng.uniform(-0.5, 0.5), rng.uniform(-0.5, 0.5)))
    pair('regular', one, two, op='conicx' if n == 1 else 'conic')


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
    header = ('# catallax-trace conic gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
