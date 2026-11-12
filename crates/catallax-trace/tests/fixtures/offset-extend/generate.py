#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace offset-extend gate.

usage: generate.py ORACLE [RAW_DIR]

ORACLE is the driver built from oracle.cpp against the pinned SMLib. Each
case runs in a fresh process under a memory and CPU cap and a timeout. RAW_DIR,
when given, receives every request, response and error stream; keep it outside
the repository. Inputs are fixed: hand-built surfaces and seeded random steps.
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


rng = random.Random(20261121)


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




def transposed(surface):
    """The same surface with U and V exchanged."""
    return Surface([surface.degrees[1], surface.degrees[0]], [surface.knots[1], surface.knots[0]],
                   [list(col) for col in zip(*surface.net)], surface.rational)


def curved(pu, pv, rational, scale=0.4):
    return Surface([pu, pv], [bezier_knots(pu), bezier_knots(pv)], random_net(pu, pv, rational, scale), rational)


def with_pole(pu, pv, rational, side):
    """A patch with one edge collapsed to a point: a pole."""
    s = curved(pu, pv, rational)
    nu, nv = len(s.net), len(s.net[0])
    if side == 'u0':
        s.net[0] = [list(s.net[0][0]) for _ in range(nv)]
    elif side == 'u1':
        s.net[nu - 1] = [list(s.net[nu - 1][0]) for _ in range(nv)]
    elif side == 'v0':
        for row in s.net:
            row[0] = list(s.net[0][0])
    return s


def step_point(u, w):
    return f'0 {v(u, w)}'


def step_normal(u, w, ul=1, vl=1):
    return f'1 {v(u, w)} {ul} {vl}'


def step_eval(u, w, simple, hu, hv, ul, vl, upper, nonzero):
    return f'2 {v(u, w)} {simple} {hu} {hv} {ul} {vl} {upper} {nonzero}'


def request(distance, allow, surface, steps):
    return f'{h(distance)} {allow} {surface.tokens()} {len(steps)} ' + ' '.join(steps)


SURFACES = [
    ('plain', curved(1, 1, False, 0.0)), ('curved-22', curved(2, 2, False)), ('curved-33', curved(3, 3, False)),
    ('curved-23r', curved(2, 3, True)),
    ('knotted', with_knots(3, 2, False, [0.4], [0.5], 0.0)),
    ('cylinder-u', cylinder()), ('cylinder-v', transposed(cylinder())),
    ('pole-u0', with_pole(2, 2, False, 'u0')), ('pole-v0', with_pole(2, 3, True, 'v0')),
]
TARGETS = [(1.3, 0.5), (-0.3, 0.5), (0.5, 1.4), (0.5, -0.4), (1.3, 1.3), (2.6, 0.5)]
ORDERS = [(0, 0), (1, 0), (1, 1), (2, 0), (2, 2), (3, 1)]
for name, surface in SURFACES:
    for d in (0.5, -0.5, 0.0):
        # one evaluation outside: a point, a normal, each order
        for target in TARGETS[:5]:
            case('ox-seq', 'regular', request(d, 1, surface, [step_point(*target)]))
        for target in ((1.3, 0.5), (-0.3, 1.2)):
            case('ox-seq', 'regular', request(d, 1, surface, [step_normal(*target, 1, 1)]))
            case('ox-seq', 'regular', request(d, 1, surface, [step_normal(*target, 0, 0)]))
        for target in ((1.3, 0.5), (0.5, -0.4)):
            for hu, hv in ORDERS:
                for simple in (0, 1):
                    upper = 1 if hu + hv <= 3 and rng.random() < 0.6 else 0
                    case('ox-seq', 'regular', request(
                        d, 1, surface, [step_eval(*target, simple, hu, hv, 1, 1, upper, 1)]))
        # the extension kept: outside, back inside, further outside, the other side
        case('ox-seq', 'regular', request(d, 1, surface, [
            step_point(0.5, 0.5), step_point(1.3, 0.5), step_point(0.5, 0.5), step_point(1.7, 0.5),
            step_normal(-0.3, 0.5), step_eval(0.5, 1.5, 0, 1, 1, 1, 1, 0, 1), step_point(0.2, 0.2)]))
        case('ox-seq', 'regular', request(d, 1, surface, [
            step_point(1.3, 0.5), step_point(1.3, 0.5), step_point(1.9, 0.5), step_point(0.5, -0.6), step_point(2.9, 3.1)]))
        # no extension allowed: the parameter is clamped
        for target in TARGETS[:2]:
            case('ox-seq', 'regular', request(d, 0, surface, [step_point(*target), step_normal(*target)]))
        # random sequences
        for _ in range(3):
            steps = []
            for _ in range(rng.randint(2, 4)):
                target = (rng.uniform(-1.2, 2.2), rng.uniform(-1.2, 2.2))
                pick = rng.randint(0, 2)
                if pick == 0:
                    steps.append(step_point(*target))
                elif pick == 1:
                    steps.append(step_normal(*target, rng.randint(0, 1), rng.randint(0, 1)))
                else:
                    hu, hv = rng.choice(ORDERS)
                    steps.append(step_eval(*target, rng.randint(0, 1), hu, hv, rng.randint(0, 1), rng.randint(0, 1),
                                           rng.randint(0, 1) if hu + hv <= 3 else 0, rng.randint(0, 1)))
            case('ox-seq', 'regular', request(d, 1, surface, steps))


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
        if run.returncode != 0:
            sys.exit(f'case {index} {op} {kind}: oracle exit {run.returncode}: {run.stderr.strip()}')
        response = run.stdout.split()
        assert response[0] == op, (index, response)
        lines.append(f'{op} {kind} | {inputs} | {" ".join(response[1:])}')
    header = ('# catallax-trace offset-extend gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
