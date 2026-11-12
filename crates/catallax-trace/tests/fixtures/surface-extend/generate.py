#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace surface-extend gate.

usage: generate.py ORACLE [RAW_DIR]

ORACLE is the driver built from oracle.cpp against the pinned SMLib. Each
case runs in a fresh process under a memory and CPU cap and a timeout. RAW_DIR,
when given, receives every request, response and error stream; keep it outside
the repository. Inputs are fixed: hand-built surfaces and seeded random ones.
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


rng = random.Random(20261120)


class Surface:
    """A tensor-product B-spline surface in homogeneous form `[x*w, y*w, z*w, w]` or `[x, y, z, NO_W]`."""

    def __init__(self, degrees, knots, net, rational):
        self.degrees = degrees
        self.knots = knots
        self.net = net  # net[i][j] = [x, y, z, w]
        self.rational = rational

    def tokens(self):
        nu, nv = len(self.net), len(self.net[0])
        flat = [c for row in self.net for point in row for c in point]
        return (f'{self.degrees[0]} {self.degrees[1]} {nu} {nv} {int(self.rational)} '
                + v(*self.knots[0]) + ' ' + v(*self.knots[1]) + ' ' + v(*flat))


def bezier_knots(p):
    return [0.0] * (p + 1) + [1.0] * (p + 1)


def random_net(pu, pv, rational, scale=1.0):
    net = []
    for i in range(pu + 1):
        row = []
        for j in range(pv + 1):
            x, y, z = i / pu * 3.0, j / pv * 2.0, scale * rng.uniform(-0.5, 0.5)
            if rational:
                w = rng.uniform(0.8, 1.25)
                row.append([x * w, y * w, z * w, w])
            else:
                row.append([x, y, z, NO_W])
        net.append(row)
    return net


def insert_knot(surface, axis, t):
    """Boehm's insertion of one copy of `t` in direction `axis`."""
    p = surface.degrees[axis]
    knots = surface.knots[axis]
    n = len(knots) - p - 2
    k = max(i for i in range(len(knots) - 1) if knots[i] <= t < knots[i + 1])
    new_knots = knots[:k + 1] + [t] + knots[k + 1:]
    new_axis_knots = list(surface.knots)
    new_axis_knots[axis] = new_knots
    # the lines of the net along the direction of insertion
    rows = [list(col) for col in zip(*surface.net)] if axis == 0 else surface.net
    out = []
    for line in rows:
        new_line = []
        for i in range(n + 2):
            if i <= k - p:
                new_line.append(line[i])
            elif i > k:
                new_line.append(line[i - 1])
            else:
                alpha = (t - knots[i]) / (knots[i + p] - knots[i])
                new_line.append([alpha * a + (1 - alpha) * b for a, b in zip(line[i], line[i - 1])])
        out.append(new_line)
    net = [list(col) for col in zip(*out)] if axis == 0 else out
    # the missing-weight marker must stay as it was
    if not surface.rational:
        for row in net:
            for point in row:
                point[3] = NO_W
    return Surface(surface.degrees, new_axis_knots, net, surface.rational)


def with_knots(pu, pv, rational, u_knots, v_knots, perturb=0.0):
    """A Bezier patch with the given interior knots inserted, then its interior net moved by `perturb`."""
    s = Surface([pu, pv], [bezier_knots(pu), bezier_knots(pv)], random_net(pu, pv, rational), rational)
    for t in u_knots:
        s = insert_knot(s, 0, t)
    for t in v_knots:
        s = insert_knot(s, 1, t)
    if perturb:
        nu, nv = len(s.net), len(s.net[0])
        for i in range(1, nu - 1):
            for j in range(1, nv - 1):
                for c in range(3):
                    s.net[i][j][c] += rng.uniform(-perturb, perturb) * (s.net[i][j][3] if rational else 1.0)
    return s


def request(surface, axis, t, times, tol):
    return f'{surface.tokens()} {axis} {h(t)} {times} {h(tol)}'


# Removable knots: inserted into a patch, so the surface does not move.
for pu, pv in ((1, 1), (2, 2), (3, 3), (1, 3), (3, 1), (2, 3)):
    for rational in (False, True):
        for u_knots, v_knots in (([0.4], [0.3]), ([0.4, 0.4], [0.3, 0.3]), ([0.25, 0.5, 0.75], [0.5]),
                                 ([0.4, 0.4, 0.4], [0.3, 0.6])):
            s = with_knots(pu, pv, rational, u_knots, v_knots)
            for axis, ks in ((0, u_knots), (1, v_knots)):
                for t in sorted(set(ks)):
                    for times in (1, 2, 5):
                        for tol in (1e-9, 1e-7):
                            case('srf-remove', 'regular', request(s, axis, t, times, tol))
# Near the limit: a perturbed net, removable only for a loose tolerance.
for pu, pv in ((2, 2), (3, 3)):
    for rational in (False, True):
        for perturb in (1e-9, 1e-6, 1e-3, 1e-1):
            s = with_knots(pu, pv, rational, [0.4, 0.4], [0.3, 0.3], perturb)
            for axis, t in ((0, 0.4), (1, 0.3)):
                for tol in (1e-7, 1e-5, 1e-3):
                    case('srf-remove', 'near-singular' if perturb in (1e-6, 1e-3) else 'regular',
                         request(s, axis, t, 2, tol))
# Knots that are not there, are at an end, or are outside the domain.
s = with_knots(3, 3, False, [0.4], [0.3])
for axis in (0, 1):
    for t in (0.0, 1.0, -0.1, 1.1, 0.123):
        case('srf-remove', 'boundary' if t in (0.0, 1.0) else 'singular', request(s, axis, t, 1, 1e-7))
case('srf-remove', 'singular', request(s, 0, 0.4, 0, 1e-7))
case('srf-remove', 'singular', request(s, 0, 0.4, 1, 0.0))
# A knot whose multiplicity is the degree, at the limit of the net.
for p in (1, 2, 3):
    s = with_knots(p, p, False, [0.5] * p, [0.5] * p)
    for axis in (0, 1):
        case('srf-remove', 'touching', request(s, axis, 0.5, p, 1e-7))
# Weights that a removal would take out of range.
s = with_knots(2, 2, True, [0.4, 0.4], [0.3, 0.3])
for row in s.net:
    for point in row:
        w = point[3]
        point[0] *= 150.0 / w; point[1] *= 150.0 / w; point[2] *= 150.0 / w; point[3] = 150.0
case('srf-remove', 'singular', request(s, 0, 0.4, 2, 1e-7))
# Random surfaces.
for _ in range(40):
    pu, pv = rng.randint(1, 3), rng.randint(1, 3)
    rational = rng.random() < 0.5
    u_knots = sorted(rng.choice([0.2, 0.35, 0.5, 0.65, 0.8]) for _ in range(rng.randint(1, 3)))
    v_knots = sorted(rng.choice([0.2, 0.35, 0.5, 0.65, 0.8]) for _ in range(rng.randint(1, 3)))
    s = with_knots(pu, pv, rational, u_knots, v_knots, rng.choice([0.0, 0.0, 1e-8, 1e-4]))
    axis = rng.randint(0, 1)
    t = rng.choice((u_knots, v_knots)[axis])
    case('srf-remove', 'regular', request(s, axis, t, rng.randint(1, 3), rng.choice([1e-9, 1e-7, 1e-4])))


# Extension through the unbounded evaluation: one extended surface carried through the calls.
def extend_request(surface, targets):
    flat = ' '.join(v(u, w) for u, w in targets)
    return f'{surface.tokens()} {len(targets)} {flat}'


def cylinder(radial_degree=2):
    """A rational quadratic circle along u (smoothly closed) times a line along v."""
    w = math.sqrt(0.5)
    ring = [(1, 0, 1), (1, 1, w), (0, 1, 1), (-1, 1, w), (-1, 0, 1), (-1, -1, w), (0, -1, 1), (1, -1, w), (1, 0, 1)]
    knots = [0.0] * 3 + [0.25, 0.25, 0.5, 0.5, 0.75, 0.75] + [1.0] * 3
    net = []
    for x, y, weight in ring:
        net.append([[x * weight, y * weight, z * weight, weight] for z in (0.0, 2.0)])
    return Surface([2, 1], [knots, [0.0, 0.0, 1.0, 1.0]], net, True)


def curved(pu, pv, rational, scale=0.4):
    return with_knots(pu, pv, rational, [], [], 0.0) if scale == 0 else Surface(
        [pu, pv], [bezier_knots(pu), bezier_knots(pv)], random_net(pu, pv, rational, scale), rational)


SEQUENCES = [
    [(0.5, 0.5)],
    [(1.2, 0.5)], [(-0.2, 0.5)], [(0.5, 1.2)], [(0.5, -0.2)], [(1.2, -0.2)], [(-0.4, 1.3)],
    [(1.0, 0.5)], [(0.0, 1.0)],
    [(2.5, 0.5)], [(-3.0, 0.5)], [(0.5, 4.0)], [(0.5, -2.0)], [(3.0, -2.0)],
    [(1.2, 0.5), (1.5, 0.5), (-0.3, 0.5), (0.5, 1.4), (0.5, -0.6), (0.2, 0.2)],
    [(1.2, 0.5), (3.7, 0.5), (0.5, 0.5), (-1.9, 1.1)],
    [(1e13, 0.5)], [(0.5, -1e13)],
]
for pu, pv in ((1, 1), (2, 2), (3, 3), (1, 3), (3, 2)):
    for rational in (False, True):
        for surface in (curved(pu, pv, rational),
                        with_knots(pu, pv, rational, [0.3, 0.7], [0.5], 0.0) if pu + pv > 2 else curved(pu, pv, rational)):
            for targets in SEQUENCES:
                kind = ('boundary' if targets[0] in ((1.0, 0.5), (0.0, 1.0))
                        else 'singular' if abs(targets[0][0]) > 1e12 or abs(targets[0][1]) > 1e12 else 'regular')
                case('srf-extend', kind, extend_request(surface, targets))
# Smoothly closed in U: the parameter wraps instead.
cyl = cylinder()
for targets in ([(1.3, 0.5)], [(-0.4, 0.5)], [(2.6, 0.5)], [(0.5, 1.3)], [(1.3, 1.3)], [(-0.4, -0.5)],
                [(1.3, 0.5), (0.5, 1.3), (-0.2, -0.2)], [(7.1, 0.5)]):
    case('srf-extend', 'regular', extend_request(cyl, targets))
# A flat surface with a collapsed edge: no tangent angle on it.
flat = curved(2, 2, False, 0.0)
for j in range(len(flat.net[0])):
    flat.net[0][j] = list(flat.net[0][0])
for targets in ([(-0.3, 0.5)], [(1.3, 0.5)], [(0.5, 1.3)]):
    case('srf-extend', 'singular', extend_request(flat, targets))
# A surface with a multiple interior knot, and one with excess multiplicity at the ends.
multi = with_knots(3, 3, False, [0.5, 0.5, 0.5], [0.5, 0.5], 0.0)
for targets in ([(1.3, 0.5)], [(-0.3, 0.5)], [(0.5, 1.3)], [(0.5, -0.3)]):
    case('srf-extend', 'touching', extend_request(multi, targets))
# Random surfaces, random targets.
for _ in range(60):
    pu, pv = rng.randint(1, 3), rng.randint(1, 3)
    rational = rng.random() < 0.5
    ks = lambda: sorted(rng.choice([0.25, 0.5, 0.75]) for _ in range(rng.randint(0, 2)))
    surface = with_knots(pu, pv, rational, ks(), ks(), rng.choice([0.0, 1e-3, 1e-1]))
    targets = [(rng.uniform(-1.5, 2.5), rng.uniform(-1.5, 2.5)) for _ in range(rng.randint(1, 4))]
    case('srf-extend', 'regular', extend_request(surface, targets))


def main():
    oracle = sys.argv[1]
    raw = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if raw:
        raw.mkdir(parents=True, exist_ok=True)
    lines = []
    for index, (op, kind, inputs) in enumerate(CASES):
        request_line = f'{op} {inputs}\n'
        run = subprocess.run(['prlimit', '--as=1073741824', '--cpu=10', 'timeout', '20', oracle],
                             input=request_line, capture_output=True, text=True)
        if raw:
            (raw / f'{index:04d}-{op}.txt').write_text(
                f'request: {request_line}exit: {run.returncode}\nstdout: {run.stdout}stderr: {run.stderr}')
        if run.returncode in (-11, -6, -7):
            # SMLib wrote or read past an array and died. The answer is undefined: the
            # case is pinned as one the port refuses, never compared.
            lines.append(f'{op}-overrun {kind} | {inputs} | i-1')
            continue
        if run.returncode != 0:
            sys.exit(f'case {index} {op} {kind}: oracle exit {run.returncode}: {run.stderr.strip()}')
        response = run.stdout.split()
        assert response[0] == op, (index, response)
        lines.append(f'{op} {kind} | {inputs} | {" ".join(response[1:])}')
    header = ('# catallax-trace surface-extend gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
