#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace curve-on-surface gate.

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
# A bicubic with a triple interior u knot whose control net folds there.
KINK = surface(3, 3, [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2], [0, 0, 0, 0, 1, 1, 1, 1],
               [(i, j, abs(i - 3) * 0.6 + 0.1 * j * j) for i in range(7) for j in range(4)])
# A plane whose parametrisation kinks along u = 1 with the normal unchanged.
FLATKINK = surface(1, 1, [0, 0, 1, 2, 2], [0, 0, 1, 1],
                   [(0, 0, 0), (0, 1, 0), (1, 0.5, 0), (1, 1.5, 0), (2, 0, 0), (2, 1, 0)])
# Biquadratic waves with many knots, sampled at five v parameters.
WAVES = surface(2, 2, [0, 0, 0] + list(range(1, 8)) + [8, 8, 8], [0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 7, 7],
                [(i, j, math.sin(1.3 * i) * math.cos(0.9 * j)) for i in range(10) for j in range(9)])
# A long thin strip: subdivision stops at the smallest side ratio.
STRIP = surface(3, 1, [0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 1, 1],
                [(1000 * i, 0.001 * j, 50 * (i % 2)) for i in range(4) for j in range(2)])
# A tight double fold with control points that double back.
FOLD = surface(3, 2, [0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 0, 1, 1, 1],
               [((0, 3, 0, 3)[i], j, (0, 1, 2, 3)[i] * 0.2 + 0.05 * j) for i in range(4) for j in range(3)])
# A plane folded by half a degree along u = 1: SMLib raises the kink to G1.
FOLDED = surface(1, 1, [0, 0, 1, 2, 2], [0, 0, 1, 1],
                 [(0, 0, 0), (0, 1, 0), (1, 0.5, 0), (1, 1.5, 0),
                  (2, 0, math.tan(math.radians(0.5))), (2, 1, math.tan(math.radians(0.5)))])
# A net whose diagonal control points coincide.
PINCH = surface(3, 3, [0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1],
                [(0, 0, 0) if i == j else (i, j, 0.2 * (i - j) ** 2) for i in range(4) for j in range(4)])
# A net with control points 1e-7 apart, below the zone tolerance.
NEARDUP = surface(3, 3, [0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1],
                  [(1 + 1e-7, j + 1e-7, math.sin(1 + 0.5 * j) + 1e-7) if i == 2 else (i, j, math.sin(i + 0.5 * j))
                   for i in range(4) for j in range(4)])


NO_Z = sys.float_info.max
NO_W = -sys.float_info.max


def hermite(p1, d1, p2, d2):
    return '0 ' + ' '.join(h(x) for x in (*p1, *d1, *p2, *d2))


def spline(degree, knots, points, rational_weights=None):
    """A planar B-spline curve as stored: `dim degree count knots... (x y z w)...`."""
    tokens = [f'2 {degree} {len(points)}'] + [h(k) for k in knots]
    for i, (x, y) in enumerate(points):
        w = rational_weights[i] if rational_weights else None
        tokens += [h(x * w), h(y * w), h(NO_Z), h(w)] if w else [h(x), h(y), h(NO_Z), h(NO_W)]
    return '1 ' + ' '.join(tokens)


# Parameter-space curves, each inside the parameter domain of the surface it goes with.
BASES = [
    ('regular', BUMP, hermite((0.4, 0.5), (1.2, 0.6), (2.4, 2.1), (1.0, -0.4))),
    ('regular', BUMP, spline(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0.3, 0.3), (1.0, 2.5), (2.0, 0.4), (2.7, 2.6)])),
    ('regular', RBUMP, spline(2, [0, 0, 0, 1, 1, 1], [(0.2, 2.8), (1.5, 0.2), (2.8, 2.6)], [1.0, 0.7, 1.4])),
    ('regular', LOFT, hermite((0.1, 0.2), (0.8, 0.1), (0.9, 0.8), (0.2, 0.9))),
    ('regular', LOFT, spline(1, [0, 0, 1, 2, 2], [(0.05, 0.1), (0.9, 0.3), (0.2, 0.95)])),
    ('regular', FLAT, hermite((0.2, 0.1), (1.0, 0.0), (1.8, 0.9), (0.0, 1.0))),
    ('boundary', FLAT, hermite((0.0, 0.0), (2.0, 1.0), (2.0, 1.0), (0.0, 0.0))),
    ('regular', CYLINDER, hermite((0.3, 0.2), (1.4, 0.2), (3.5, 1.8), (0.3, 0.1))),
    ('singular', POLE, hermite((0.1, 0.0), (0.0, 0.0), (0.9, 0.0), (0.0, 0.0))),
    ('regular', POLE, spline(2, [0, 0, 0, 1, 1, 1], [(0.05, 0.1), (0.5, 0.9), (0.95, 0.2)])),
    ('regular', WAVES, hermite((0.5, 0.5), (3.0, 2.0), (6.5, 5.5), (-1.0, 2.0))),
]
TS = (0.0, 0.2, 0.5, 0.83, 1.0)


def call(oracle, request):
    run = subprocess.run(['prlimit', '--as=2147483648', '--cpu=20', 'timeout', '40', oracle],
                         input=request + '\n', capture_output=True, text=True)
    if run.returncode != 0:
        sys.exit(f'oracle exit {run.returncode} on {request[:60]}: {run.stderr.strip()}')
    return run.stdout.split()


def box(lo, hi):
    return ' '.join(h(x) for x in (*lo, *hi))


def main():
    oracle = sys.argv[1]
    raw = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if raw:
        raw.mkdir(parents=True, exist_ok=True)
    lines = []
    log = []

    def emit(op, kind, inputs, extra=''):
        request = f'{op} {inputs}'
        response = call(oracle, request)
        assert response[0] == op, response
        log.append(f'request: {request}\nstdout: {" ".join(response)}\n')
        lines.append(f'{op} {kind} | {inputs}{extra} | {" ".join(response[1:])}')
        return response[1:]

    for kind, surface, curve in BASES:
        for t in TS:
            for left in (0, 1):
                emit('cos-eval', 'boundary' if t in (0.0, 1.0) else kind, f'{surface} {curve} {h(t)} {left}')
        # the curve's own parameter range is read from the first call on it: use [0, 1]
        count = int(call(oracle, f'cos-node {surface} {curve} -1')[1][1:])
        for index in sorted({0, count // 2, count - 1}):
            node = call(oracle, f'cos-node {surface} {curve} {index}')
            box2 = [struct.unpack('<d', struct.pack('<Q', int(t[1:], 16)))[0] for t in node[2:8]]
            lo, hi = box2[:3], box2[3:]
            size = [b - a for a, b in zip(lo, hi)]
            mid = [(a + b) / 2 for a, b in zip(lo, hi)]
            node_tokens = ' '.join(token[1:] for token in node[2:10])
            for tol in (1e-3, 0.05):
                # a box that overlaps, and boxes a gap away in each direction and along the
                # diagonal: no gap (tangent), the tolerance (touching) and wider ones
                def around(c):
                    return [x - 0.02 for x in c], [x + 0.02 for x in c]

                cases = [(around(mid), kind)]
                for gap, gap_kind in ((0.0, 'tangent'), (tol, 'touching'), (0.01, kind), (0.2, kind), (3.0, kind)):
                    for axis in (0, 2):
                        for sign in (1, -1):
                            c = list(mid)
                            c[axis] += sign * (size[axis] / 2 + gap + 0.02)
                            cases.append((around(c), gap_kind))
                    c = [m + (s / 2 + gap + 0.02) for m, s in zip(mid, size)]
                    cases.append((around(c), gap_kind))
                for (blo, bhi), case_kind in cases:
                    inputs = f'{surface} {curve} {index} {box(blo, bhi)} {h(tol)}'
                    request = f'cos-branch {inputs}'
                    response = call(oracle, request)
                    assert response[0] == 'cos-branch', response
                    log.append(f'request: {request}\nstdout: {" ".join(response)}\n')
                    lines.append(f'cos-branch {case_kind} | {inputs} {node_tokens} | {" ".join(response[1:])}')
    if raw:
        for i, entry in enumerate(log):
            (raw / f'{i:05d}.txt').write_text(entry)
    header = ('# catallax-trace curve-on-surface gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
