#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace local solver gate.

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


def v(*xs):
    return ' '.join(h(x) for x in xs)


def solve1(kind, family, params, guess, interval=None, periodic=0, algorithm=0,
           acceptable=1e-6, desired=-1.0):
    assert kind in KINDS
    lo, hi = interval if interval else (0.0, 0.0)
    params = list(params) + [0.0] * (4 - len(params))
    CASES.append(('solve1', kind, f'{family} {v(*params)} {h(guess)} {int(interval is not None)} '
                                  f'{v(lo, hi)} {periodic} {algorithm} {v(acceptable, desired)}'))


def solven(kind, family, params, guess, intervals=None, handler=1, max_hits=3, cutting=0,
           acceptable=1e-6):
    assert kind in KINDS
    params = list(params) + [0.0] * (4 - len(params))
    boxes = intervals or [(0.0, 0.0, 0)] * len(guess)
    box_tokens = ' '.join(f'{v(lo, hi)} {per}' for lo, hi, per in boxes)
    CASES.append(('solven', kind, f'{family} {v(*params)} {len(guess)} {v(*guess)} '
                                  f'{int(intervals is not None)} {box_tokens} {handler} {max_hits} '
                                  f'{cutting} {h(acceptable)}'))


# One variable: simple roots, extrema, double roots, near misses, oscillation,
# evaluation failures, bounds, periodic wrapping and the Brent algorithms.
for guess in (0.2, 1.6, 2.4, 3.7, 10.0):
    solve1('regular', 1, (1, 2, 3), guess)
solve1('regular', 1, (1, 2, 3), 2.2, interval=(1.5, 4.0))
solve1('singular', 1, (1, 2, 3), 2 - 1 / math.sqrt(3))
solve1('singular', 1, (1, 2, 3), 2 + 1 / math.sqrt(3), interval=(0.0, 4.0))
for guess in (0.0, 3.0, 1.0 + 1e-9):
    solve1('tangent', 2, (1, 0), guess)
solve1('tangent', 2, (1, 0), 0.0, acceptable=1e-12)
solve1('touching', 2, (1, 1e-10), 0.0)
solve1('touching', 2, (1, 1e-10), 2.5, interval=(0.0, 4.0))
solve1('regular', 2, (1, 1e-3), 0.0)
solve1('regular', 2, (1, -0.25), 0.0)
solve1('regular', 2, (1, -0.25), 1.0)
solve1('near-singular', 3, (0.3, 1e-12), 1.0)
solve1('near-singular', 3, (0.3, 1e-12), -0.7, interval=(-2.0, 2.0))
solve1('singular', 4, (2, 1.5), 1.0)
solve1('boundary', 4, (2, 1.5), 1.0, interval=(0.0, 1.4))
solve1('boundary', 1, (5, 6, 7), 1.0, interval=(0.0, 4.0))
solve1('boundary', 1, (1, 2, 3), 5.0, interval=(0.0, 4.0))
solve1('regular', 2, (0.5, -0.04), 0.95, interval=(0.0, 1.0), periodic=1)
solve1('regular', 2, (0.5, -0.04), 0.05, interval=(0.0, 1.0), periodic=1)
for algorithm in (1, 2, 3):
    solve1('regular', 5, (1, 3, 0.5), 2.1, interval=(0.0, 4.0), algorithm=algorithm)
    solve1('touching', 5, (1, 3, 0.0), 1.9, interval=(0.0, 4.0), algorithm=algorithm)
solve1('regular', 1, (1, 2, 3), 1.6, desired=1e-3)
solve1('regular', 1, (-1e-3, 1e-3, 5), 0.5, acceptable=1e-12)
solve1('near-singular', 1, (1, 1 + 1e-9, 3), 0.0)

# Several variables: crossing, tangent and missing line-circle cases, a
# finite-difference Jacobian, tangent and separate circles, a rank-deficient
# system, clamping, periodic wrapping, boundary handlers and step cutting.
solven('regular', 1, (1, 0.5, 0.2), (0.5, 0.5))
solven('regular', 1, (1, 0.5, 0.2), (-1.0, -0.2))
solven('tangent', 1, (1, 0, 1), (0.3, 0.8))
solven('touching', 1, (1, 0, 1 + 1e-9), (0.3, 0.8))
solven('regular', 1, (1, 0, 1.2), (0.3, 0.8))
solven('regular', 2, (1, 0.5, 0.2), (0.5, 0.5))
solven('tangent', 2, (1, 0, 1), (0.3, 0.8))
solven('regular', 3, (0, 1), (0.4, 0.7))
solven('tangent', 3, (0, 2), (1.1, 0.3))
solven('touching', 3, (0, 2 + 1e-8), (1.1, 0.3))
solven('regular', 3, (0, 2.5), (1.1, 0.3))
solven('singular', 4, (3, 6, 0), (0.1, 0.2, 0.3))
solven('singular', 4, (3, 7, 0), (0.1, 0.2, 0.3))
box = [(0.0, 0.5, 0), (0.0, 0.5, 0)]
solven('boundary', 1, (1, 0.5, 0.2), (0.4, 0.4), intervals=box)
solven('boundary', 1, (1, 0.5, 0.2), (0.4, 0.4), intervals=box, handler=0, max_hits=2)
solven('boundary', 1, (1, 0.5, 0.2), (0.4, 0.4), intervals=box, handler=2)
solven('regular', 1, (1, 0.5, 0.2), (0.4, 0.4), intervals=[(0.0, 1.0, 0), (0.0, 1.0, 0)])
solven('regular', 1, (1, 0.5, 0.2), (0.9, 0.1), intervals=[(-1.0, 1.0, 1), (-1.0, 1.0, 0)])
solven('regular', 1, (1, 0.5, 0.2), (0.4, 0.4), intervals=[(-1234567.0, 1234567.0, 0), (0.0, 1.0, 0)])
solven('regular', 1, (1, 0.5, 0.2), (3.0, -4.0), cutting=1)
solven('touching', 1, (1, 0, 1.2), (0.3, 0.8), cutting=1)
solven('regular', 3, (0, 1), (0.4, 0.7), acceptable=1e-14)


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
    header = ('# catallax-trace local solver gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
