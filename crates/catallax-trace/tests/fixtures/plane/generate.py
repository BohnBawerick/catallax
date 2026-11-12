#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace plane gate.

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
INTERSECT, MINIMIZE, NORMALIZE = 0, 3, 5


def h(x):
    return format(struct.unpack('<Q', struct.pack('<d', float(x)))[0], '016x')


def v(*xs):
    return ' '.join(h(x) for x in xs)


def plane(kind, origin, x, y, scale=(1, 1), domain=(0, 0, 1, 1)):
    return f'{kind} {v(*origin, *x, *y, *scale, *domain)}'


def curve(degree, knots, points):
    tokens = [h(k) for k in knots]
    for p in points:
        tokens += [h(p[0]), h(p[1]), h(p[2]), h(1.0)]
    return f'{degree} {len(points)} 0 ' + ' '.join(tokens)


def case(op, kind, text):
    assert kind in KINDS
    CASES.append((op, kind, text))


INF = 1234567.0
FLAT = plane(1, (0, 0, 1), (1, 0, 0), (0, 1, 0))
TILT = plane(1, (1, 2, 0.5), (1, 1, 0), (-1, 1, 1))
BOUNDED = plane(0, (0.3, -0.2, 0.7), (0.6, 0.8, 0), (0, 0, 1), (2, 0.5), (-1, -2, 3, 4))
SKEW = plane(0, (0, 0, 0), (2, 0, 0), (1e-3, 3, 0), (1, 1), (0, 0, 2, 2))
NEAR_ORTHO = plane(0, (0, 0, 0), (1, 0, 0), (1e-13, 1, 0), (1, 1), (0, 0, 2, 2))
SCALED = plane(0, (0.1, 0.7, -0.3), (0.36, 0.48, 0.8), (0.8, -0.6, 0), (0.3, 1.7), (-0.7, -1.1, 2.9, 3.3))
PLANES = (FLAT, TILT, BOUNDED, SKEW, NEAR_ORTHO, SCALED)

for p in PLANES:
    case('net', 'regular', p)
    for uv, kind in (((0.25, 0.75), 'regular'), ((0, 0), 'boundary'), ((-2.5, 5), 'boundary'),
                     ((3, 4), 'boundary'), ((1e-9, -1e-9), 'near-singular')):
        case('point', kind, f'{p} {v(*uv)}')
        case('point', kind, f'{p} {v(uv[0] + 0.137, uv[1] - 0.291)}')
        case('normal', kind, f'{p} {v(*uv)}')
        case('flags', kind, f'{p} {v(*uv)}')
        for hu, hv, upper, simple in ((0, 0, 0, 0), (1, 1, 0, 0), (2, 2, 1, 0), (1, 0, 0, 1), (3, 2, 0, 0),
                                      (1, 1, 1, 1)):
            case('eval', kind, f'{p} {v(*uv)} {hu} {hv} {upper} {simple}')
# Point drops: on the plane, off it, beyond the domain, just off it, with
# meeting and target tolerances.
for p, domain in ((FLAT, (-INF, -INF, INF, INF)), (FLAT, (0, 0, 1, 1)), (TILT, (-INF, -INF, INF, INF)),
                  (BOUNDED, (-1, -2, 3, 4)), (BOUNDED, (0, 0, 1, 1)), (SKEW, (0, 0, 2, 2)),
                  (SCALED, (-0.7, -1.1, 2.9, 3.3)), (SCALED, (0.1, 0.2, 0.9, 0.7)),
                  (BOUNDED, (-10, -10, 10, 10))):
    for test, kind in (((0.4, 0.3, 1.0), 'touching'), ((0.4, 0.3, 2.5), 'regular'), ((5, -3, 0), 'boundary'),
                       ((0.5, 0.5, 1 + 1e-9), 'near-singular'), ((1.2, 0.7, 0.4), 'regular'),
                       ((0.0, 0.0, 1.0), 'singular')):
        for operation in (MINIMIZE, NORMALIZE, INTERSECT):
            case('drop', kind, f'{p} {operation} {v(*domain, *test, 1e-5)} 0 {v(0)}')
        case('drop', kind, f'{p} {MINIMIZE} {v(*domain, *test, 1e-5)} 1 {v(0.5)}')
        for operation in (MINIMIZE, INTERSECT, NORMALIZE):
            case('global', kind, f'{p} {operation} {v(*domain, *test, 1e-5)}')
# Curves against planes: across, tangent, touching at an end, a point curve,
# and curves through a bounded plane's domain edge.
CROSS = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0.2, 0.3, -1), (0.5, 0.9, 2), (0.8, 0.1, 0), (0.9, 0.8, 3)])
TOUCH = curve(2, [0, 0, 0, 1, 1, 1], [(0.1, 0.1, 2), (0.5, 0.5, 0), (0.9, 0.9, 2)])
END = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0.3, 0.3, 1), (0.4, 0.5, 1.5), (0.6, 0.5, 2), (0.7, 0.7, 3)])
POINT = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0.25, 0.75, 1.0)] * 4)
WAVE = curve(3, [0, 0, 0, 0, 1, 2, 2, 2, 2], [(-0.5, 0.2, 0), (0.2, 0.5, 2), (0.7, 0.1, 0), (1.2, 0.6, 2),
                                             (1.6, 0.4, 0)])
for c, kind in ((CROSS, 'regular'), (TOUCH, 'tangent'), (END, 'touching'), (POINT, 'singular'),
                (WAVE, 'regular')):
    for domain in ((-INF, -INF, INF, INF), (0, 0, 1, 1), (0.3, 0, 1, 0.45)):
        interval = (0, 2) if c is WAVE else (0, 1)
        case('cs', kind if domain[0] < 0 else 'boundary', f'{FLAT} {c} {v(*domain, *interval, 1e-5)}')
NEAR_TOUCH = curve(2, [0, 0, 0, 1, 1, 1], [(0.1, 0.1, 2), (0.5, 0.5, 1e-7), (0.9, 0.9, 2)])
case('cs', 'near-singular', f'{FLAT} {NEAR_TOUCH} {v(-INF, -INF, INF, INF, 0, 1, 1e-5)}')
TILT_CROSS = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0, 0, 0), (1, 3, 1), (2, 0, 3), (3, 3, 2)])
case('cs', 'regular', f'{TILT} {TILT_CROSS} {v(-INF, -INF, INF, INF, 0, 1, 1e-5)}')
case('cs', 'regular', f'{BOUNDED} {curve(3, [0, 0, 0, 0, 1, 1, 1, 1], [(1, -2, 0), (0, 1, 2), (1, 0, -1), (0, 1, 3)])} '
     f'{v(-1, -2, 3, 4, 0, 1, 1e-5)}')
# Straight curves (SmPlane::GlobalLineIntersect): crossing inside and
# outside the curve's interval and the domain, a collinear cubic, other
# parameter ranges, parallel lines away from the plane, beside its domain and
# beyond the tolerance, and lines lying in the plane or within the tolerance
# of it (range answers).
def straight(p, q, knots=(0, 1)):
    a, b = knots
    return curve(1, [a, a, b, b], [p, q])


def collinear(p, q):
    pts = [tuple(p[k] + (q[k] - p[k]) * i / 3 for k in range(3)) for i in range(4)]
    return curve(3, [0, 0, 0, 0, 1, 1, 1, 1], pts)


ALL = (-INF, -INF, INF, INF)
for c, interval, domain, kind in (
        (straight((0.2, 0.3, -1), (0.8, 0.6, 3)), (0, 1), ALL, 'regular'),
        (straight((0.2, 0.3, -1), (0.8, 0.6, 3)), (0, 1), (0, 0, 1, 1), 'boundary'),
        (straight((0.2, 0.3, -1), (0.8, 0.6, 3)), (0.6, 1), ALL, 'boundary'),
        (straight((1.5, 0.3, -1), (2.5, 0.6, 3)), (0, 1), (0, 0, 1, 1), 'boundary'),
        (straight((0.2, 0.3, 0.0), (0.8, 0.6, 1.0)), (0, 1), ALL, 'touching'),
        (straight((0.2, 0.3, -1), (0.8, 0.6, 3), (2, 5)), (2, 5), ALL, 'regular'),
        (collinear((0.1, 0.9, 3), (0.7, 0.2, -2)), (0, 1), ALL, 'regular'),
        (straight((0, 0, 2), (1, 1, 2)), (0, 1), ALL, 'tangent'),
        (straight((0, 0, 1 + 1e-4), (1, 1, 1 + 1e-4)), (0, 1), ALL, 'near-singular'),
        (straight((5, 5, 1), (6, 6, 1)), (0, 1), (0, 0, 1, 1), 'boundary')):
    case('cs', kind, f'{FLAT} {c} {v(*domain, *interval, 1e-5)}')
for p_, c in ((TILT, straight((0, 0, 0), (3, 1, 2))), (BOUNDED, straight((1, -2, 0), (0, 1, 3))),
              (SCALED, straight((0.5, 0.5, -1), (0.2, 0.9, 1)))):
    case('cs', 'regular', f'{p_} {c} {v(*ALL, 0, 1, 1e-5)}')
for c, domain, kind in ((straight((0.2, 0.2, 1), (0.8, 0.6, 1)), ALL, 'tangent'),
                        (straight((0.2, 0.2, 1 + 1e-6), (0.8, 0.6, 1 + 1e-6)), ALL, 'near-singular'),
                        (straight((-0.5, 0.2, 1), (1.5, 0.6, 1)), (0, 0, 1, 1), 'boundary'),
                        (straight((0.2, 0.3, 1 + 1e-9), (0.8, 0.6, 1 + 1e-9)), ALL, 'near-singular')):
    case('cs', kind, f'{FLAT} {c} {v(*domain, 0, 1, 1e-5)}')


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
    header = ('# catallax-trace plane gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
