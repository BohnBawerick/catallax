#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace curve-compatibility gate.

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


rng = random.Random(20261112)


def clamped(degree, interior):
    return [0.0] * (degree + 1) + list(interior) + [3.0] * (degree + 1)


def random_points(count, dim, rational):
    pts = []
    for _ in range(count):
        x, y, z = (rng.uniform(-4, 4) for _ in range(3))
        w = rng.uniform(.5, 2) if rational else None
        if rational:
            pts.append([x * w, y * w, z * w if dim == 3 else NO_Z, w])
        else:
            pts.append([x, y, z if dim == 3 else NO_Z, NO_W])
    return pts


class Curve:
    def __init__(self, dim, degree, knots, points):
        self.dim, self.degree, self.knots, self.points = dim, degree, list(knots), [list(p) for p in points]

    def token(self):
        flat = [c for p in self.points for c in p]
        return f'{self.dim} {self.degree} {len(self.points)} ' + v(*self.knots) + ' ' + v(*flat)

    def insert(self, u):
        """Boehm insertion of one knot, on the stored coordinates."""
        p, U, P = self.degree, self.knots, self.points
        k = max(i for i in range(len(U) - 1) if U[i] <= u < U[i + 1])
        newP = []
        for i in range(len(P) + 1):
            if i <= k - p:
                newP.append(P[i])
            elif i > k:
                newP.append(P[i - 1])
            else:
                a = (u - U[i]) / (U[i + p] - U[i])
                q = [a * P[i][c] + (1 - a) * P[i - 1][c] for c in range(4)]
                if P[i][2] == NO_Z:
                    q[2] = NO_Z
                if P[i][3] == NO_W:
                    q[3] = NO_W
                newP.append(q)
        return Curve(self.dim, p, U[:k + 1] + [u] + U[k + 1:], newP)


def make(dim, degree, interior, rational):
    knots = clamped(degree, interior)
    return Curve(dim, degree, knots, random_points(len(knots) - degree - 1, dim, rational))


# Degree elevation.
for degree in (1, 2, 3, 4):
    for rational in (False, True):
        for dim in (2, 3):
            for interior in ([], [1.0], [0.8, 2.1], [1.0, 1.0] if degree >= 2 else [1.0], [0.7, 0.7, 2.2] if degree >= 3 else [0.7, 2.2]):
                c = make(dim, degree, interior, rational)
                for t in (0, 1, 2, 3):
                    case('elevate', 'regular', c.token(), t)
# C0 interior knots (multiplicity equal to the degree) and a long elevation.
for degree in (2, 3):
    c = make(3, degree, [1.5] * degree, False)
    for t in (1, 2):
        case('elevate', 'tangent', c.token(), t)
    c = make(3, degree, [1.5] * degree, True)
    case('elevate', 'tangent', c.token(), 1)
c = make(3, 2, [1.0, 2.0], False)
case('elevate', 'boundary', c.token(), 15)
case('elevate', 'singular', c.token(), 16)
case('elevate', 'singular', make(3, 4, [1.0], False).token(), 14)
# Unclamped knots and a repeated end knot are refused by SMLib.
c = make(3, 2, [1.0, 2.0], False)
c.knots = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
case('elevate', 'singular', c.token(), 1)
c = make(3, 2, [1.0], False)
c.points = [[1., 2., 3., NO_W]] * len(c.points)
case('elevate', 'singular', c.token(), 2)

# Knot removal: a knot inserted into a curve is removable.
for degree in (2, 3, 4):
    for rational in (False, True):
        for dim in (2, 3):
            base = make(dim, degree, [1.2, 2.1], rational)
            for u in (1.6, 2.5, 0.4):
                once = base.insert(u)
                twice = once.insert(u)
                case('remove', 'regular', once.token(), h(u), 1, h(1e-6))
                case('remove', 'regular', twice.token(), h(u), 2, h(1e-6))
                case('remove', 'regular', twice.token(), h(u), 1, h(1e-6))
            # A genuine interior knot does not come out at a tight tolerance.
            case('remove', 'tangent', base.token(), h(1.2), 1, h(1e-9))
            case('remove', 'tangent', base.token(), h(2.1), 1, h(1e-12))
            # ... but does at a loose one.
            case('remove', 'touching', base.token(), h(1.2), 1, h(1e9))
            case('remove', 'touching', base.token(), h(1.2), 3, h(1e300))
            # A parameter that is not a knot, and the end knots.
            case('remove', 'boundary', base.token(), h(1.5), 1, h(1e-6))
            case('remove', 'boundary', base.token(), h(0.0), 1, h(1e-6))
            case('remove', 'boundary', base.token(), h(3.0), 1, h(1e-6))
            case('remove', 'near-singular', base.token(), h(1.2 + 1e-13), 1, h(1e9))
            case('remove', 'regular', base.token(), h(1.2), 0, h(1e9))
# Removal that C0 knots allow and weights forbid.
c = make(3, 3, [1.5, 1.5, 1.5], False)
case('remove', 'tangent', c.token(), h(1.5), 3, h(1e300))
rat = make(3, 3, [1.2], True)
rat.points[2][3] = 400.0
case('remove', 'singular', rat.insert(2.0).token(), h(2.0), 1, h(1e300))
rat = make(3, 3, [1.2], True)
for p in rat.points:
    p[3] = 0.0005
case('remove', 'singular', rat.token(), h(1.2), 1, h(1e300))

# Merging knot vectors.
def knotset(count, base, jitter=0.0):
    out = []
    for _ in range(count):
        ks = sorted(set(round(rng.choice(base) + (rng.uniform(-jitter, jitter) if jitter else 0), 12) for _ in range(rng.randint(1, 4))))
        degree = 2
        out.append([0.0] * (degree + 1) + ks + [3.0] * (degree + 1))
    return out


def knots_case(vectors, tol, kind):
    tokens = [str(len(vectors)), h(tol)]
    for vec in vectors:
        tokens += [str(len(vec)), v(*vec)]
    case('knots', kind, *tokens)


base = [0.5, 1.0, 1.5, 2.0, 2.5]
for _ in range(10):
    knots_case(knotset(2, base), 1e-8, 'regular')
for _ in range(6):
    knots_case(knotset(3, base), 1e-8, 'regular')
for _ in range(6):
    knots_case(knotset(3, base, 1e-9), 1e-8, 'near-singular')
for _ in range(4):
    knots_case(knotset(2, base, 0.02), 0.05, 'touching')
# Knots whose gap is below the tolerance.
knots_case([[0, 0, 0, 1.0, 1.0 + 1e-9, 3, 3, 3], [0, 0, 0, 1.0, 3, 3, 3]], 1e-8, 'near-singular')
knots_case([[0, 0, 0, 1.0, 1.0004, 3, 3, 3], [0, 0, 0, 1.0002, 3, 3, 3]], 1e-3, 'near-singular')
knots_case([[0, 0, 0, 1.0, 3, 3, 3], [0, 0, 0, 1.0, 3, 3, 3]], 1e-8, 'tangent')
knots_case([[0, 0, 0, 1.0, 1.0, 3, 3, 3], [0, 0, 0, 1.0, 3, 3, 3]], 1e-8, 'tangent')
knots_case([[0, 0, 0, 3, 3, 3], [0, 0, 0, 3, 3, 3]], 1e-8, 'boundary')
knots_case([[0, 0, 0, 1, 3, 3, 3]], 1e-8, 'boundary')
knots_case([[0, 0, 0, 1, 3, 3, 3], [0, 0, 0, 1, 3.0000001, 3.0000001, 3.0000001]], 1e-8, 'singular')
knots_case([[0, 0, 0, 1, 3, 3, 3], [0, 0, 0, 1, 2, 2, 2]], 1e-8, 'singular')
knots_case([[0, 0, 0, 1, 3, 3, 3], [1, 1, 1, 2, 3, 3, 3]], 1e-8, 'singular')

# Synchronising curves.
def sync_case(curves, tol, kind):
    # SMLib keeps each curve's dimension label though the data are made
    # three-dimensional; the gate pins those cases apart.
    op = 'sync-mixed-dimension' if len({c.dim for c in curves}) > 1 else 'sync'
    case(op, kind, len(curves), h(tol), ' '.join(c.token() for c in curves))


for _ in range(3):
    for rational in (False, True):
        for dim in (2, 3):
            a = make(dim, 2, [0.8, 2.1], rational)
            b = make(dim, 3, [1.4], rational)
            c = make(dim, 1, [], rational)
            sync_case([a, b], 1e-8, 'regular')
            sync_case([a, b, c], 1e-8, 'regular')
# Mixed rationality and dimension.
sync_case([make(2, 2, [1.0], False), make(3, 2, [2.0], True)], 1e-8, 'regular')
sync_case([make(3, 3, [1.0], True), make(2, 1, [], False)], 1e-8, 'regular')
# Different parameter ranges are scaled to the unit interval.
for degree in (2, 3):
    a = make(3, degree, [1.0], False)
    b = make(3, degree, [0.5], False)
    b.knots = [k / 3.0 * 7.0 + 1.0 for k in b.knots]
    sync_case([a, b], 1e-8, 'regular')
    c = make(3, degree, [1.0], False)
    c.knots = [k + 1e-9 for k in c.knots]
    sync_case([a, c], 1e-8, 'near-singular')
    c = make(3, degree, [1.0], False)
    c.knots = [k + 1e-7 for k in c.knots]
    sync_case([a, c], 1e-8, 'touching')
# Identical and nearly identical knots.
a = make(3, 3, [1.0, 2.0], False)
sync_case([a, make(3, 3, [1.0, 2.0], False)], 1e-8, 'tangent')
sync_case([a, make(3, 3, [1.0 + 1e-10, 2.0], False)], 1e-8, 'near-singular')
sync_case([a], 1e-8, 'boundary')
sync_case([a, a], 1e-8, 'tangent')
# A last span no longer than the tolerance: the knot is removed, or refused.
d = make(3, 2, [1.0, 3.0 - 1e-9], False)
sync_case([d, make(3, 2, [1.0], False)], 1e-8, 'singular')
d = make(3, 2, [1.0, 3.0 - 1e-9], False).insert(3.0 - 1e-9)
sync_case([d, make(3, 2, [1.0], False)], 1e-8, 'singular')
d = make(3, 2, [], False)
d.knots = [0.0, 0.0, 0.0, 3.0 - 1e-9, 3.0, 3.0, 3.0][:0] or d.knots
sync_case([make(3, 2, [3.0 - 1e-9], False), make(3, 3, [1.5], False)], 1e-8, 'singular')
# Bezier curves of different degree.
sync_case([make(3, 2, [], False), make(3, 4, [], False)], 1e-8, 'regular')
sync_case([make(3, 1, [], True), make(2, 3, [], False)], 1e-8, 'regular')
# A large tolerance merges distinct knots.
sync_case([make(3, 2, [1.0], False), make(3, 2, [1.05], False)], 0.1, 'touching')
sync_case([make(3, 2, [1.0, 1.04], False), make(3, 2, [1.02], False)], 0.1, 'touching')


def main():
    oracle = sys.argv[1]
    raw = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if raw:
        raw.mkdir(parents=True, exist_ok=True)
    lines = []
    for index, (op, kind, inputs) in enumerate(CASES):
        request_op = {'sync-mixed-dimension': 'sync'}.get(op, op)
        request = f'{request_op} {inputs}\n'
        run = subprocess.run(['prlimit', '--as=1073741824', '--cpu=10', 'timeout', '20', oracle],
                             input=request, capture_output=True, text=True)
        if raw:
            (raw / f'{index:04d}-{op}.txt').write_text(
                f'request: {request}exit: {run.returncode}\nstdout: {run.stdout}stderr: {run.stderr}')
        if run.returncode != 0:
            print(f'case {index} {op} {kind}: oracle exit {run.returncode}: {run.stderr.strip()}', file=sys.stderr)
            lines.append(f'{op} {kind} | {inputs} | crash{run.returncode}')
            continue
        response = run.stdout.split()
        assert response[0] == request_op, (index, response)
        lines.append(f'{op} {kind} | {inputs} | {" ".join(response[1:])}')
    header = ('# catallax-trace curve-compatibility gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
