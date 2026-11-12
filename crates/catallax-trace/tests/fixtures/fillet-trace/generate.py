#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace fillet-trace gate.

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
NO_W = -sys.float_info.max


def h(x):
    return format(struct.unpack('<Q', struct.pack('<d', float(x)))[0], '016x')


def v(*xs):
    return ' '.join(h(x) for x in xs)


CASES = []


def case(op, kind, *tokens):
    assert kind in ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
    CASES.append((op, kind, ' '.join(str(t) for t in tokens)))


rng = random.Random(20261119)


def knots(degree, count):
    interior = count - degree - 1
    return [0.0] * (degree + 1) + [(i + 1) / (interior + 1) for i in range(interior)] + [1.0] * (degree + 1)


def surface(du, dv, nu, nv, fn, bump=0.0, rational=False):
    """A patch of the map `fn(s, t)` over `[0,1]^2` with `bump`-sized perturbations."""
    ku, kv = knots(du, nu), knots(dv, nv)
    vals = []
    for i in range(nu):
        for j in range(nv):
            s, t = i / (nu - 1), j / (nv - 1)
            x, y, z = fn(s, t)
            if 0 < i < nu - 1 and 0 < j < nv - 1:
                x += rng.uniform(-bump, bump); y += rng.uniform(-bump, bump); z += rng.uniform(-bump, bump)
            w = rng.uniform(.8, 1.2) if rational else None
            vals += [x * w, y * w, z * w, w] if rational else [x, y, z, NO_W]
    return f'{du} {dv} {nu} {nv} {int(rational)} ' + v(*ku) + ' ' + v(*kv) + ' ' + v(*vals)


def floor(bump, rational=False, degree=3):
    return surface(degree, degree, degree + 1, degree + 1, lambda s, t: (3 * s, 2 * t, 0.0), bump, rational)


def wall(bump, rational=False, degree=3):
    return surface(degree, degree, degree + 1, degree + 1, lambda s, t: (0.0, 2 * s, 3 * t), bump, rational)


def ceiling(bump):
    # normal along -z: u runs along y, v along x
    return surface(3, 3, 4, 4, lambda s, t: (3 * t, 2 * s, 1.0), bump)


def head(r1, r2, o1, o2, tol, s1, s2):
    return f'{h(r1)} {h(r2)} {o1} {o2} {h(tol)} {s1} {s2}'


def domain():
    return v(0.0, 0.0, 1.0, 1.0)



def consistent(r1, r2, y):
    """uv on each surface of the offset-surface crossing at height y (flat patches)."""
    return (r2 / 3.0, y / 2.0), (y / 2.0, r1 / 3.0)


def stretched_floor():
    # v runs as the square root of y, so the parameter speed changes fivefold
    return surface(3, 3, 4, 4, lambda s, t: (3 * s, 2 * t * t, 0.0))


def stretched_wall():
    return surface(3, 3, 4, 4, lambda s, t: (0.0, 2 * s * s, 3 * t))


def full(r1, r2, o1, o2, tol, s1, s2, d1=None, d2=None):
    return f'{head(r1, r2, o1, o2, tol, s1, s2)} {d1 or domain()} {d2 or domain()}'


def tail(approx, angle=math.radians(20.0), direction=1.0, through=0):
    return f'{h(approx)} {h(angle)} {h(direction)} {through}'


def consistent_stretched(r1, r2, y):
    s = math.sqrt(y / 2.0)
    return (r2 / 3.0, s), (s, r1 / 3.0)



def cubic_floor():
    return surface(3, 3, 4, 4, lambda s, t: (3 * s, 2.0 if round(3 * t) == 3 else 0.0, 0.0))


def cubic_wall():
    return surface(3, 3, 4, 4, lambda s, t: (0.0, 2.0 if round(3 * s) == 3 else 0.0, 3 * t))


def consistent_cubic(r1, r2, y):
    t = (y / 2.0) ** (1.0 / 3.0)
    return (r2 / 3.0, t), (t, r1 / 3.0)


FAMILIES = [
    ('flat', floor(0.0), wall(0.0), consistent),
    ('stretched', stretched_floor(), stretched_wall(), consistent_stretched),
    ('cubic', cubic_floor(), cubic_wall(), consistent_cubic),
]
Y_DIR = v(0.0, 1.0, 0.0)


def request(extend_before=0, extend_after=0, start_dir=None, end_dir=None, approx=None, angle=None,
            extension=0.0, points=()):
    """The request tokens: flags, optional directions and tolerances, extension distance, through points."""
    out = [str(extend_before), str(extend_after)]
    out += ['1', start_dir] if start_dir else ['0']
    out += ['1', end_dir] if end_dir else ['0']
    out += ['1', h(approx)] if approx is not None else ['0']
    out += ['1', h(angle)] if angle is not None else ['0']
    out += [h(extension), str(len(points))]
    out += [v(*(p[0] + p[1])) for p in points]
    return ' '.join(out)


def on_curve(family, r1, r2, ys):
    return [family[3](r1, r2, y) for y in ys]


def case_for(kind, family, request_tokens, r1=0.4, r2=0.4, o1=1, o2=1, tol=1e-6):
    case('fi-trace', kind, full(r1, r2, o1, o2, tol, family[1], family[2]), request_tokens)


for family in FAMILIES:
    for r1, r2 in ((0.4, 0.4), (0.3, 0.5)):
        base = dict(r1=r1, r2=r2)
        # two points with a start direction: the trace runs from the first to the second
        for ys in ((0.3, 1.5), (0.1, 0.9)):
            for direction in (Y_DIR, v(0.0, -1.0, 0.0)):
                pts = on_curve(family, r1, r2, ys)
                case_for('regular', family, request(start_dir=direction, points=pts), **base)
                case_for('regular', family, request(start_dir=direction, points=pts[::-1]), **base)
        # several points, no directions, every pair of extension flags
        pts = on_curve(family, r1, r2, (0.2, 0.6, 1.0, 1.6))
        for before in (0, 1):
            for after in (0, 1):
                case_for('regular', family, request(before, after, points=pts), **base)
                case_for('regular', family, request(before, after, start_dir=Y_DIR, points=pts), **base)
        # one point
        for before in (0, 1):
            for after in (0, 1):
                case_for('regular', family, request(before, after, points=on_curve(family, r1, r2, (0.8,))), **base)
        # tolerances and an end direction
        pts = on_curve(family, r1, r2, (0.3, 1.2))
        case_for('regular', family, request(start_dir=Y_DIR, end_dir=Y_DIR, points=pts), **base)
        case_for('regular', family, request(start_dir=Y_DIR, approx=1e-4, points=pts), **base)
        case_for('regular', family, request(start_dir=Y_DIR, approx=1e-2, angle=math.radians(5.0), points=pts), **base)
        case_for('regular', family, request(1, 1, start_dir=Y_DIR, angle=1e-3, points=pts), **base)
        # a duplicated through point, and an extension distance
        dup = on_curve(family, r1, r2, (0.3, 0.8, 0.8, 1.4))
        case_for('regular', family, request(points=dup), **base)
        case_for('regular', family, request(0, 1, extension=0.5, points=pts[:1]), **base)
        case_for('regular', family, request(0, 1, extension=0.5, points=pts), **base)

# Starts near the ends of the domains.
for family in FAMILIES[:2]:
    for ys, direction in (((1.9, 1.98), Y_DIR), ((0.1, 0.02), v(0.0, -1.0, 0.0))):
        pts = on_curve(family, 0.4, 0.4, ys)
        case_for('boundary', family, request(start_dir=direction, points=pts))
        case_for('boundary', family, request(1, 1, points=pts[:1]))
    case_for('boundary', family, request(1, 1, points=on_curve(family, 0.4, 0.4, (1.995,))))

# Failures: no through points, a first point off the surfaces, zero radii, a zero direction.
flat = FAMILIES[0]
case_for('singular', flat, request(points=()))
(u1, v1), (u2, v2) = consistent(0.4, 0.4, 0.9)
case_for('singular', flat, request(points=[((u1, v1), (u2 + 0.2, v2))]))
case_for('singular', flat, request(points=[((u1, v1), (u2, v2)), ((u1, v1 + 0.2), (u2 + 0.2, v2))]))
for r1, r2 in ((0.0, 0.4), (0.4, 0.0), (1e-13, 0.4)):
    case_for('singular', flat, request(start_dir=Y_DIR, points=on_curve(flat, 0.4, 0.4, (0.3, 1.2))), r1=r1, r2=r2)
case_for('singular', flat, request(start_dir=v(0.0, 0.0, 0.0), points=on_curve(flat, 0.4, 0.4, (0.3, 1.2))))

# Parallel offset surfaces: a tangent point, and the trace's reversal from it.
for rc in (0.5, 0.45):
    pts = [((0.3, 0.4), (0.5, 0.6)), ((0.5, 0.5), (0.5, 0.5))]
    case('fi-trace', 'tangent', full(rc, 1.0 - rc, 1, 0, 1e-6, floor(0.0), ceiling(0.0)), request(points=pts[:1]))
    case('fi-trace', 'tangent', full(rc, 1.0 - rc, 1, 0, 1e-6, floor(0.0), ceiling(0.0)),
         request(start_dir=Y_DIR, points=pts[1:]))
    case('fi-trace', 'tangent', full(rc, 1.0 - rc, 1, 0, 1e-6, floor(0.0), ceiling(0.0)),
         request(1, 1, points=pts))
# Coincident surfaces.
case('fi-trace', 'touching', full(0.4, 0.4, 1, 1, 1e-6, floor(0.0), floor(0.0)),
     request(points=[((0.5, 0.5), (0.5, 0.5))]))
case('fi-trace', 'touching', full(0.4, 0.4, 1, 1, 1e-6, floor(0.0), floor(0.0)),
     request(1, 1, start_dir=Y_DIR, points=[((0.5, 0.5), (0.5, 0.5)), ((0.5, 0.6), (0.5, 0.6))]))


# Short spans inside the domains, which end on their end point before any step reaches a side.
for family in FAMILIES:
    for r1, r2 in ((0.4, 0.4), (0.3, 0.5)):
        for ys in ((0.6, 1.0), (0.8, 1.2), (0.5, 0.7), (0.9, 1.0)):
            pts = on_curve(family, r1, r2, ys)
            case_for('regular', family, request(start_dir=Y_DIR, points=pts), r1=r1, r2=r2)
            case_for('regular', family, request(start_dir=v(0.0, -1.0, 0.0), points=pts[::-1]), r1=r1, r2=r2)
            case_for('regular', family, request(start_dir=Y_DIR, approx=1e-4, points=pts), r1=r1, r2=r2)
            case_for('regular', family, request(start_dir=Y_DIR, angle=math.radians(2.0), points=pts), r1=r1, r2=r2)
        mid = on_curve(family, r1, r2, (0.6, 0.8, 1.0, 1.2))
        case_for('regular', family, request(start_dir=Y_DIR, points=mid), r1=r1, r2=r2)
        case_for('regular', family, request(start_dir=Y_DIR, end_dir=Y_DIR, points=mid), r1=r1, r2=r2)


def main():
    oracle = sys.argv[1]
    raw = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if raw:
        raw.mkdir(parents=True, exist_ok=True)
    lines = []
    for index, (op, kind, inputs) in enumerate(CASES):
        request_line = f'{op} {inputs}\n'
        run = subprocess.run(['prlimit', '--as=1073741824', '--cpu=20', 'timeout', '40', oracle],
                             input=request_line, capture_output=True, text=True)
        if raw:
            (raw / f'{index:04d}-{op}.txt').write_text(
                f'request: {request_line}exit: {run.returncode}\nstdout: {run.stdout}stderr: {run.stderr}')
        if run.returncode != 0:
            sys.exit(f'case {index} {op} {kind}: oracle exit {run.returncode}: {run.stderr.strip()}')
        response = run.stdout.split()
        assert response[0] == op, (index, response)
        lines.append(f'{op} {kind} | {inputs} | {" ".join(response[1:])}')
    header = ('# catallax-trace fillet-trace gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
