#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace surface cache gate.

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
# Flat zigzag nets turning more than 540 degrees, so the chord test passes and
# the angle test falls back to sm_ComputeSurfConstants: a degree-7 net in u, and
# the same shape in other degrees, in v, in both directions, rational, with a
# small height, and across two patches.
def zigzag(du, dv=1, amplitude=3.0, rational=False, height=0.0, in_v=False):
    nu, nv = du + 1, dv + 1
    ku = [0] * (du + 1) + [1] * (du + 1)
    kv = [0] * (dv + 1) + [1] * (dv + 1)
    pts = []
    for i in range(nu):
        for j in range(nv):
            a = (0, amplitude)[i % 2]
            pts.append((a, j + 0.1 * i, height * ((i * 3 + j) % 4), 1.0 + 0.2 * ((i + j) % 3) if rational else 1.0))
    return surface(du, dv, ku, kv, pts, rational)


ZIGZAG = zigzag(7)
ZIGZAGS = [
    ('boundary', ZIGZAG), ('boundary', zigzag(5)), ('boundary', zigzag(6)), ('boundary', zigzag(8)),
    ('regular', zigzag(7, amplitude=1.5)), ('regular', zigzag(7, rational=True)),
    ('regular', zigzag(6, 2, rational=True)), ('near-singular', zigzag(7, height=0.01)),
    ('regular', zigzag(5, 3, amplitude=2.0, height=0.02)),
]
# The same in v: the surface with its two directions exchanged.
ZIGZAGS.append(('regular', surface(1, 7, [0, 0, 1, 1], [0] * 8 + [1] * 8,
                                   [(j + 0.1 * i, (0, 3)[j % 2], 0.0) for i in range(2) for j in range(8)])))
# Turning in both directions.
ZIGZAGS.append(('regular', surface(5, 5, [0] * 6 + [1] * 6, [0] * 6 + [1] * 6,
                                   [((0, 3)[i % 2] + 0.05 * j, (0, 2)[j % 2] + 0.05 * i, 0.0)
                                    for i in range(6) for j in range(6)])))
# Two patches whose nets each turn more than 540 degrees in u.
ZIGZAGS.append(('regular', surface(5, 1, [0] * 6 + [1] * 5 + [2] * 6, [0, 0, 1, 1],
                                   [((0, 3)[i % 2], j + 0.1 * i, 0.0) for i in range(11) for j in range(2)])))
# Control points that go forward and back along a nearly straight line: the net
# turns 180 degrees at every vertex, and the surface itself hardly turns, so the
# sampled angles are small where the net's are large.
def forward_back(du, dv=1, rational=False, in_v=False, tilt=0.0):
    nu, nv = du + 1, dv + 1
    ku = [0] * (du + 1) + [1] * (du + 1)
    kv = [0] * (dv + 1) + [1] * (dv + 1)
    pts = []
    for i in range(nu):
        for j in range(nv):
            x = i + (1.5 if i % 2 else 0.0)
            pts.append((x, j + tilt * i, 0.0, 1.0 + 0.15 * ((i + j) % 3) if rational else 1.0))
    return surface(du, dv, ku, kv, pts, rational)


ZIGZAGS += [
    ('regular', forward_back(5)), ('regular', forward_back(7)), ('regular', forward_back(6, 2)),
    ('regular', forward_back(7, rational=True)), ('regular', forward_back(5, tilt=0.05)),
    ('near-singular', forward_back(7, tilt=0.5)),
    ('regular', surface(1, 7, [0, 0, 1, 1], [0] * 8 + [1] * 8,
                        [(j + 0.1 * i, i * 0.0 + (j + (1.5 if j % 2 else 0.0)), 0.0) for i in range(2) for j in range(8)])),
]
for kind, z in ZIGZAGS:
    for op in ('cache', 'net', 'cont', 'bbox'):
        CASES.append((op, kind, z))
SURFACES = (('boundary', FOLDED), ('singular', PINCH), ('near-singular', NEARDUP), ('regular', BUMP), ('regular', RBUMP), ('regular', LOFT), ('singular', POLE), ('regular', FLAT),
            ('boundary', CYLINDER), ('boundary', KINK), ('boundary', FLATKINK), ('regular', WAVES),
            ('near-singular', STRIP), ('regular', FOLD))
for kind, s in SURFACES:
    for op in ('cache', 'net', 'cont', 'bbox', 'planes'):
        CASES.append((op, kind, s))
# Polar boxes of the normals: the leaves' boxes as GetPolarBox builds them
# on demand (preorder), and boxes over sub-domains without and with the
# whole surface's box as the guess. Two bilinear patches join: one flat
# (a single normal) and one twisted.
UNIT_FLAT = surface(1, 1, [0, 0, 1, 1], [0, 0, 1, 1], [(0, 0, 0), (0, 1, 0.5), (1, 0, 0), (1, 1, 0.5)])
TWIST = surface(1, 1, [0, 0, 1, 1], [0, 0, 1, 1], [(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0.4)])
# A single biquadratic patch (not sampled as analytic, unlike degree (2, 1)).
DOME = surface(2, 2, [0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1],
               [(i, j, 0.3 * (1 - (i - 1) ** 2) + 0.2 * (1 - (j - 1) ** 2)) for i in range(3) for j in range(3)])
RANGES = {id(BUMP): (3, 3), id(RBUMP): (3, 3), id(LOFT): (1, 1), id(POLE): (1, 1), id(FLAT): (2, 1),
          id(CYLINDER): (4, 2), id(KINK): (2, 1), id(FLATKINK): (2, 1), id(WAVES): (8, 7), id(STRIP): (1, 1),
          id(FOLD): (1, 1), id(FOLDED): (2, 1), id(PINCH): (1, 1), id(NEARDUP): (1, 1), id(UNIT_FLAT): (1, 1),
          id(TWIST): (1, 1), id(DOME): (1, 1)}


def hx(*xs):
    return ' '.join(h(x) for x in xs)


for kind, s in SURFACES + (('regular', UNIT_FLAT), ('regular', TWIST), ('regular', DOME)):
    CASES.append(('polar', kind, s))
    U, V = RANGES[id(s)]
    for a, b, c, d in ((0, 0, 1, 1), (0.25, 0.1, 0.6, 0.45), (0.5, 0.5, 1, 1), (0, 0.8, 0.2, 1)):
        CASES.append(('polarbox', kind, f'{s} {hx(a * U, b * V, c * U, d * V)}'))


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
    header = ('# catallax-trace surface cache gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
