#!/usr/bin/env python3
"""Regenerate cases.txt for the catallax-trace search gate.

usage: generate.py ORACLE [RAW_DIR]

ORACLE is the driver built from oracle.cpp against the pinned SMLib. Each
case runs in a fresh process under a memory and CPU cap and a timeout. RAW_DIR,
when given, receives every request, response and error stream; keep it outside
the repository. Inputs are fixed: hand-picked tangent, touching, singular and
boundary cases plus seeded random regular cases.

Object table, solution tokens and tree tokens are described in oracle.cpp.
"""
import math
import pathlib
import random
import struct
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent


def h(x):
    return format(struct.unpack('<Q', struct.pack('<d', float(x)))[0], '016x')


def v(*xs):
    return ' '.join(h(x) for x in xs)


CASES = []


def case(op, kind, *tokens):
    assert kind in ('regular', 'tangent', 'touching', 'singular', 'near-singular', 'boundary')
    CASES.append((op, kind, ' '.join(str(t) for t in tokens)))


rng = random.Random(20261112)

# Solution types.
NODES, SINGLE, RANGE, POLY, UNKNOWN = range(5)
# Objects.
LINE_A, LINE_B, PLANE_A, PLANE_B, VERTEX, EDGE, FACE, BASE = range(8)
# Solver operations (SmSolverOperationType order).
(INTERSECT, WIREFRAME, INTERSECTION_TEST, MINIMIZE, MAXIMIZE, NORMALIZE, AT_DISTANCE, FIND,
 RAYFIRE, DIRECTED3D_MIN, DIRECTED3D_MAX, SIGNED_DIRECTED_MIN, DIRECTED_MIN, DIRECTED_MAX,
 PROJECTED_MIN, PROJECTED_MAX, PROJECTED_INTERSECT, ROTATED_INTERSECT, PERSPECTIVE_INTERSECT,
 SURFACE_PROJECTED, ANGLE_MIN, SIGNED_ANGLE_MIN, PROJECTED_TANGENCY, SIGNED_PIVOT_MIN,
 TANGENT_THROUGH_POINT) = range(25)
# Operation categories, requests, sort keys.
CAT_MIN, CAT_MAX, CAT_SIGNED, CAT_OTHER = range(4)
SR_SINGLE, SR_ALL, SR_AMBIG, SR_NODES = range(4)
BY_VALUE, BY_FIRST, BY_FIRST_TWO, BY_VALUE_PARAM, BY_ALL = range(5)
# Auxiliary data types.
AUX_NONE, AUX_BEZIER_CURVE, AUX_CURVE, AUX_SURFACE, AUX_HEAD, AUX_SURFDATA, AUX_AUX = range(7)


def category(op):
    if op in (DIRECTED_MIN, PROJECTED_MIN, MINIMIZE, ANGLE_MIN, RAYFIRE):
        return CAT_MIN
    if op in (DIRECTED_MAX, PROJECTED_MAX, MAXIMIZE):
        return CAT_MAX
    if op in (SIGNED_DIRECTED_MIN, SIGNED_ANGLE_MIN, SIGNED_PIVOT_MIN, DIRECTED3D_MIN, DIRECTED3D_MAX):
        return CAT_SIGNED
    return CAT_OTHER


def S(kind, nvars, objs, value, start, end_value=0.0, end=(), nobjs=None):
    """One solution's tokens."""
    objs = list(objs)
    count = len(objs) if nobjs is None else nobjs
    objs += [-1] * (5 - len(objs))
    start = list(start) + [0.0] * (4 - len(start))
    end = list(end) + [0.0] * (4 - len(end))
    assert len(start) == 4 and len(end) == 4
    return ' '.join([str(kind), str(nvars), str(count)] + [str(o) for o in objs]
                    + [h(value)] + [h(x) for x in start] + [h(end_value)] + [h(x) for x in end])


def point(value, p, objs=(LINE_A, PLANE_A), nvars=3):
    return S(SINGLE, nvars, objs, value, p)


def span(v0, p0, v1, p1, objs=(LINE_A, PLANE_A), nvars=3):
    return S(RANGE, nvars, objs, v0, p0, v1, p1)


# Solution comparison: ==, classification of start and end, IsEqual, identity.
CP = (LINE_A, PLANE_A)
for a, b, kind in (
    (point(1e-9, [0.5, 1.0, 0.25]), point(1e-9, [0.5, 1.0, 0.25]), 'regular'),
    (point(1e-9, [0.5, 1.0, 0.25]), point(2e-9, [0.5, 1.0, 0.25]), 'regular'),
    (point(1e-9, [0.5, 1.0, 0.25]), point(1e-9, [0.5 + 1e-13, 1.0, 0.25]), 'touching'),
    (point(1e-9, [0.5, 1.0, 0.25]), point(1e-9, [0.5 + 1e-11, 1.0, 0.25]), 'touching'),
    (point(1e-9, [0.5, 1.0, 0.25]), point(1e-9, [0.5 + 5e-5, 1.0, 0.25]), 'touching'),
    (point(1e-9, [0.5, 1.0, 0.25]), point(1e-9, [0.5 + 2e-4, 1.0, 0.25]), 'touching'),
    (point(1e-9, [0.5, 1.0, 0.25]), point(1e-9, [0.5, 1.0 + 1.9e-4, 0.25]), 'touching'),
    (point(1e-9, [0.5, 1.0, 0.25]), point(1e-9, [0.5, 1.0, 0.25 + 3e-4]), 'touching'),
    (point(1e-9, [0.5, 1.0, 0.25]), point(1e-9 + 5e-7, [0.5, 1.0, 0.25]), 'boundary'),
    (point(1e-9, [0.5, 1.0, 0.25]), point(1e-9 + 2e-6, [0.5, 1.0, 0.25]), 'boundary'),
    (point(1e-9, [0.5, 1.0, 0.25]), point(1e-9, [0.5, 1.0, 0.25], (LINE_B, PLANE_A)), 'regular'),
    (point(1e-9, [0.5, 1.0, 0.25]), point(1e-9, [0.5, 1.0, 0.25], (LINE_A, PLANE_B)), 'regular'),
    (point(1e-9, [0.5, 1.0, 0.25]), point(1e-9, [0.5, 1.0], (LINE_A,), 2), 'regular'),
    (point(1e-9, [0.5, 1.0, 0.25]), span(1e-9, [0.5, 1.0, 0.25], 1e-9, [0.6, 1.0, 0.25]), 'regular'),
    (point(0.0, [0.5, 1.0, 0.25]), point(0.0, [0.5, 1.0, 0.25]), 'singular'),
    (point(1e-9, [0.0, 0.0, 0.0]), point(1e-9, [1e-13, 0.0, 0.0]), 'singular'),
    (point(1e-9, [1e6, 1.0, 0.25]), point(1e-9, [1e6 + 1e-7, 1.0, 0.25]), 'boundary'),
    (point(1e-9, [1e6, 1.0, 0.25]), point(1e-9, [1e6 + 1e-3, 1.0, 0.25]), 'boundary'),
):
    case('solcmp', kind, a, b)

# Range against point, range against range.
R = span(1e-9, [0.2, 0.5, 0.5], 2e-9, [0.6, 0.9, 0.7])
for b, kind in (
    (point(1e-9, [0.2, 0.5, 0.5]), 'touching'),
    (point(2e-9, [0.6, 0.9, 0.7]), 'touching'),
    (point(1.5e-9, [0.4, 0.7, 0.6]), 'regular'),
    (point(1e-9, [0.1, 0.4, 0.4]), 'regular'),
    (point(1e-9, [0.7, 1.0, 0.8]), 'regular'),
    (point(1e-9, [0.4, 0.7, 0.9]), 'near-singular'),
    (point(1e-9, [0.2, 0.7, 0.6]), 'near-singular'),
    (point(1e-9, [0.6, 0.5, 0.7]), 'near-singular'),
    (point(1e-3, [0.2, 0.5, 0.5]), 'boundary'),
    (point(2e-9, [0.6 + 1e-13, 0.9, 0.7]), 'touching'),
    (point(2e-9, [0.6 + 1e-11, 0.9, 0.7]), 'touching'),
    (span(1e-9, [0.3, 0.6, 0.55], 1e-9, [0.5, 0.8, 0.65]), 'regular'),
    (span(1e-9, [0.2, 0.5, 0.5], 1e-9, [0.4, 0.7, 0.6]), 'touching'),
    (span(1e-9, [0.4, 0.7, 0.6], 2e-9, [0.6, 0.9, 0.7]), 'touching'),
    (span(1e-9, [0.1, 0.4, 0.4], 1e-9, [0.3, 0.6, 0.55]), 'regular'),
    (span(1e-9, [0.5, 0.8, 0.65], 1e-9, [0.8, 1.0, 0.9]), 'regular'),
    (span(1e-9, [0.0, 0.3, 0.3], 1e-9, [0.1, 0.4, 0.4]), 'regular'),
    (span(1e-9, [0.7, 1.0, 0.8], 1e-9, [0.9, 1.2, 0.95]), 'regular'),
    (span(1e-9, [0.6, 0.9, 0.7], 1e-9, [0.8, 1.1, 0.8]), 'touching'),
    (span(1e-9, [0.1, 0.4, 0.4], 1e-9, [0.2, 0.5, 0.5]), 'touching'),
    (span(1e-9, [0.1, 0.4, 0.4], 1e-9, [0.8, 1.1, 0.8]), 'regular'),
    (span(1e-9, [0.2, 0.5, 0.5], 2e-9, [0.6, 0.9, 0.7]), 'tangent'),
    (span(1e-9, [0.2, 0.5, 0.5], 2e-3, [0.6, 0.9, 0.7]), 'tangent'),
    (span(1e-9, [0.3, 0.5, 0.5], 1e-9, [0.7, 0.9, 0.7]), 'regular'),
    (span(1e-9, [0.3, 0.6, 0.5], 1e-9, [0.7, 0.9, 0.7]), 'singular'),
    (span(1e-9, [0.4, 0.7, 0.6], 1e-9, [0.4, 0.7, 0.6]), 'singular'),
    (span(1e-9, [0.2, 0.5, 0.5], 2e-9, [0.6, 0.9, 0.7], (LINE_B, PLANE_A)), 'regular'),
    (S(NODES, 0, [], 0.0, []), 'regular'),
    (S(UNKNOWN, 3, CP, 1e-9, [0.4, 0.7, 0.6]), 'regular'),
):
    case('solcmp', kind, R, b)
Z = span(1e-9, [0.5, 0.5, 0.5], 1e-9, [0.5, 0.5, 0.5])
case('solcmp', 'singular', Z, point(1e-9, [0.5, 0.5, 0.5]))
case('solcmp', 'singular', Z, Z)
case('solcmp', 'singular', point(1e-9, [0.5, 0.5, 0.5]), Z)
case('solcmp', 'singular', S(NODES, 0, [], 0.0, []), S(NODES, 0, [], 0.0, []))
# Reverse-parametrised ranges and surface pairs.
RV = span(1e-9, [0.6, 0.9, 0.7], 2e-9, [0.2, 0.5, 0.5])
for b, kind in ((point(1e-9, [0.4, 0.7, 0.6]), 'regular'), (point(1e-9, [0.7, 1.0, 0.8]), 'regular'),
                (span(1e-9, [0.5, 0.8, 0.65], 1e-9, [0.1, 0.4, 0.4]), 'regular')):
    case('solcmp', kind, RV, b)
SS = (PLANE_A, PLANE_B)
case('solcmp', 'regular', point(1e-9, [0.5, 0.5, 1.5, 2.5], SS, 4), point(1e-9, [0.5, 0.5, 1.5, 2.5], SS, 4))
case('solcmp', 'touching', point(1e-9, [0.5, 0.5, 1.5, 2.5], SS, 4), point(1e-9, [0.5, 0.5 + 1e-3, 1.5, 2.5], SS, 4))
case('solcmp', 'touching', point(1e-9, [0.5, 0.5, 1.5, 2.5], SS, 4), point(1e-9, [0.5, 0.5, 1.5, 2.5 + 2e-4], SS, 4))
case('solcmp', 'regular', point(1e-9, [0.5, 0.5, 1.5, 2.5], SS, 4), point(1e-9, [0.5, 0.5, 1.5, 2.5], (PLANE_B, PLANE_A), 4))
case('solcmp', 'regular', point(1e-9, [0.5, 0.5, 1.5, 2.5], (EDGE, FACE), 3), point(1e-9, [0.5, 0.5, 1.5, 2.5], (EDGE, FACE), 3))
case('solcmp', 'singular', point(1e-9, [0.5], (VERTEX,), 1), point(1e-9, [0.5], (VERTEX,), 1))
case('solcmp', 'singular', point(1e-9, [0.5], (BASE,), 1), point(1e-9, [0.5 + 1e-3], (BASE,), 1))
case('solcmp', 'singular', S(SINGLE, 1, [], 1e-9, [0.5], nobjs=1), S(SINGLE, 1, [], 1e-9, [0.5], nobjs=1))
for obj in (BASE, VERTEX):
    for d in (0.0, 5e-13, 1.5e-12, 5e-11, 1.5e-10, 3e-10, 1e-8):
        case('solcmp', 'touching', point(1e-9, [0.5], (obj,), 1), point(1e-9, [0.5 + d], (obj,), 1))
        case('solcmp', 'touching', span(1e-9, [0.5], 1e-9, [0.7], (obj,), 1), span(1e-9, [0.5], 1e-9, [0.7 + d], (obj,), 1))
for _ in range(36):
    p1 = [rng.random() for _ in range(3)]
    p2 = [x + rng.choice((0, 0, 1e-13, 1e-11, 1e-9, 1e-5, 1e-3, 0.1)) for x in p1]
    v1, v2 = 10 ** rng.uniform(-12, -3), 10 ** rng.uniform(-12, -3)
    if rng.random() < 0.5:
        case('solcmp', 'regular', point(v1, p1), point(v2, p2))
    else:
        q1 = sorted(rng.random() for _ in range(2))
        case('solcmp', 'regular', span(v1, p1, v1, [q1[1]] + p1[1:]), point(v2, p2))

# Objects of a solution: index of an object, parameter object and interval.
CONFIGS = (((LINE_A, PLANE_A), 3), ((PLANE_A, LINE_A), 3), ((PLANE_A, PLANE_B), 4), ((LINE_A, LINE_B), 2),
           ((EDGE, FACE), 3), ((VERTEX, LINE_A), 1), ((BASE,), 1),
           ((LINE_B, PLANE_B, LINE_A), 4), ((PLANE_A, VERTEX, LINE_B), 3), ((VERTEX, PLANE_A), 2), ((), 1))
for objs, nvars in CONFIGS:
    for p in range(nvars + 2):
        target = objs[p % len(objs)] if objs else -1
        for flag in (1, 0):
            kind = 'singular' if p >= nvars or not objs or BASE in objs or VERTEX in objs else 'regular'
            case('solidx', kind, S(SINGLE, nvars, objs, 1e-9, [0.5] * 4), p, target, flag)

# AddSortedSolution.
ZERO_TOL = 1e-6


def adds(op, req, cat, tol, best, nvars, key, combine, initial, new, kind='regular'):
    case('addsorted', kind, op, req, cat, v(tol, best), nvars, key, int(combine), len(initial),
         *initial, len(new), *new)


for key in range(5):
    for combine in (True, False) if key in (BY_FIRST, BY_ALL) else (True,):
        # points in all orders, a duplicate, a near duplicate and a rejected big value
        pts = [point(1e-9, [0.5, 1.0, 0.25]), point(2e-9, [0.2, 0.3, 0.9]), point(1e-9, [0.8, 1.5, 0.4]),
               point(1e-9, [0.5, 1.0, 0.25]), point(5e-10, [0.5 + 1e-13, 1.0, 0.25]),
               point(1e20, [0.1, 0.1, 0.1]), point(3e-9, [0.2, 0.9, 0.1]), point(1e-9, [0.2 + 1e-7, 0.3, 0.9])]
        adds(INTERSECT, SR_ALL, CAT_OTHER, 1e-7, 0.0, 3, key, combine, [], pts)
        adds(INTERSECT, SR_SINGLE, CAT_OTHER, 1e-7, 0.0, 3, key, combine, [], pts[:3])
for key in (BY_FIRST, BY_ALL):
    for combine in (True, False):
        rngs = [span(1e-9, [0.1, 0.2, 0.3], 1e-9, [0.3, 0.4, 0.5]),     # first
                span(1e-9, [0.3, 0.4, 0.5], 1e-9, [0.5, 0.6, 0.7]),     # contiguous with first
                span(1e-9, [0.6, 0.6, 0.7], 1e-9, [0.8, 0.7, 0.8]),     # disjoint
                span(1e-9, [0.45, 0.55, 0.65], 1e-9, [0.7, 0.65, 0.75]),  # overlaps both
                point(1e-9, [0.35, 0.45, 0.55]),                        # inside
                point(1e-12, [0.15, 0.25, 0.35]),                       # inside, far better
                point(1e-9, [0.1, 0.2, 0.3]),                           # on the start
                span(1e-9, [0.0, 0.1, 0.2], 1e-9, [0.1, 0.2, 0.3])]     # touches start
        adds(INTERSECT, SR_ALL, CAT_OTHER, 1e-7, 0.0, 3, key, combine, [], rngs, 'touching')
        adds(MINIMIZE, SR_ALL, CAT_MIN, 1e-7, 1e20, 3, key, combine, [], rngs, 'touching')
# a point split of a range under INTERSECT and the same under another operation
for op in (INTERSECT, FIND, NORMALIZE, MINIMIZE):
    adds(op, SR_ALL, category(op), 1e-7, 1e20, 3, BY_ALL, True,
         [span(1e-3, [0.1, 0.2, 0.3], 2e-3, [0.9, 0.9, 0.9])], [point(1e-7, [0.5, 0.5, 0.5])], 'touching')
    adds(op, SR_ALL, category(op), 1e-7, 1e20, 3, BY_ALL, True,
         [span(1e-9, [0.1, 0.2, 0.3], 2e-9, [0.9, 0.9, 0.9])], [point(1e-9, [0.5, 0.5, 0.5])], 'touching')
    adds(op, SR_ALL, category(op), 1e-7, 1e20, 3, BY_ALL, True,
         [point(1e-7, [0.5, 0.5, 0.5])], [span(1e-3, [0.1, 0.2, 0.3], 2e-3, [0.9, 0.9, 0.9])], 'touching')
# minimise and maximise: best answer tracking and trimming
vals = [3e-3, 1e-3, 5e-4, 1e-4 + 1e-8, 1e-4, 2e-4, 0.5, 1e-4 + 5e-8]
mins = [point(x, [rng.random(), rng.random(), rng.random()]) for x in vals]
for op in (MINIMIZE, MAXIMIZE, DIRECTED_MIN, PROJECTED_MAX, RAYFIRE):
    for key in (BY_VALUE, BY_ALL):
        for req in (SR_ALL, SR_SINGLE):
            best = 1e20 if category(op) == CAT_MIN else (0.0 if category(op) == CAT_MAX else 0.0)
            adds(op, req, category(op), 1e-4, best, 3, key, True, [], mins)
# values at the keep/cull boundary of the best answer plus the tolerance
edge = [point(x, [0.1 * (i + 1), 0.3, 0.5]) for i, x in enumerate(
    (1e-4, 2e-4, 2e-4 + 1e-12, 2e-4 - 1e-12, 2e-4 + 1e-9, 1e-4 + 1e-12, 5e-5, 5e-5 + 1e-4, 5e-5 + 1e-4 + 1e-12))]
for op in (MINIMIZE, DIRECTED_MIN, RAYFIRE):
    adds(op, SR_ALL, CAT_MIN, 1e-4, 1e20, 3, BY_VALUE, True, [], edge, 'boundary')
for op in (MAXIMIZE, DIRECTED_MAX):
    adds(op, SR_ALL, CAT_MAX, 1e-4, 0.0, 3, BY_VALUE, True, [],
         [point(x, [0.1 * (i + 1), 0.3, 0.5]) for i, x in enumerate(
             (3e-4, 2e-4, 2e-4 - 1e-12, 2e-4 + 1e-12, 1e-4, 1e-4 - 1e-12, 1e-4 + 1e-12, 4e-4))], 'boundary')
for op in (SIGNED_DIRECTED_MIN, DIRECTED3D_MIN):
    adds(op, SR_ALL, CAT_SIGNED, 1e-4, 0.0, 3, BY_VALUE, True, [],
         [point(x, [0.1 * (i + 1), 0.3, 0.5]) for i, x in enumerate(
             (-1e-3, -1e-3 + 1e-4, -1e-3 + 1e-4 + 1e-12, -2e-3, 1e-4, -1e-3 - 1e-4))], 'boundary')
for op in (MINIMIZE, DIRECTED_MIN, PROJECTED_MIN):
    adds(op, SR_ALL, CAT_MIN, 0.25, 1e20, 3, BY_VALUE, True, [],
         [point(x, [0.1 * (i + 1), 0.3, 0.5]) for i, x in enumerate((0.5, 0.75, 0.75 + 2.0 ** -40, 0.625, 0.25))], 'boundary')
for op in (MAXIMIZE, DIRECTED_MAX, PROJECTED_MAX):
    adds(op, SR_ALL, CAT_MAX, 0.25, 0.0, 3, BY_VALUE, True, [],
         [point(x, [0.1 * (i + 1), 0.3, 0.5]) for i, x in enumerate((1.0, 0.75, 0.75 - 2.0 ** -40, 0.875, 2.0))], 'boundary')
for key in (BY_FIRST, BY_ALL, BY_FIRST_TWO):
    for dx in (1e-13, 1e-9, 1e-7, 1e-5):
        adds(INTERSECT, SR_ALL, CAT_OTHER, 1e-7, 0.0, 3, key, False,
             [span(1e-9, [0.5, 0.5, 0.5], 1e-9, [0.6, 0.6, 0.6])], [point(1e-9, [0.5 + dx, 0.3, 0.5])], 'touching')
        adds(INTERSECT, SR_ALL, CAT_OTHER, 1e-7, 0.0, 3, key, False,
             [point(1e-9, [0.5 + dx, 0.3, 0.5])], [span(1e-9, [0.5, 0.5, 0.5], 1e-9, [0.6, 0.6, 0.6])], 'touching')
for point_value, low, high in ((5e-8, 1e-9, 1e-10), (5e-8, 1e-10, 1e-9), (5e-7, 1e-9, 1e-9), (5e-9, 1e-9, 1e-9),
                               (2e-9, 1e-9, 2e-11), (2e-9, 2e-11, 1e-9)):
    adds(INTERSECT, SR_ALL, CAT_OTHER, 1e-7, 0.0, 3, BY_ALL, True,
         [span(low, [0.1, 0.2, 0.3], high, [0.9, 0.9, 0.9])], [point(point_value, [0.5, 0.5, 0.5])], 'touching')
    adds(INTERSECT, SR_ALL, CAT_OTHER, 1e-7, 0.0, 3, BY_ALL, True,
         [point(point_value, [0.5, 0.5, 0.5])], [span(low, [0.1, 0.2, 0.3], high, [0.9, 0.9, 0.9])], 'touching')
adds(MINIMIZE, SR_ALL, CAT_MIN, 0.0, 1e20, 3, BY_VALUE, True, [], mins, 'singular')
adds(MAXIMIZE, SR_ALL, CAT_MAX, 0.0, 0.0, 3, BY_VALUE, True, [], mins, 'singular')
adds(MINIMIZE, SR_ALL, CAT_MIN, 1e-4, 4e-8, 3, BY_VALUE, True, [], mins, 'boundary')
# signed categories
neg = [point(x, [rng.random(), rng.random(), rng.random()]) for x in (-1e-3, 2e-3, -4e-3, 0.0, -0.0, 1e-5, -2e-3)]
for op in (SIGNED_DIRECTED_MIN, DIRECTED3D_MIN, SIGNED_PIVOT_MIN):
    for best in (1e20, 0.0, -1e-6):
        for key in (BY_VALUE, BY_ALL):
            adds(op, SR_ALL, CAT_SIGNED, 1e-3, best, 3, key, True, [], neg)
# an existing array that is not sorted, with ranges
init = [point(1e-9, [0.9, 0.9, 0.9]), point(1e-9, [0.1, 0.1, 0.1]), span(1e-9, [0.3, 0.3, 0.3], 1e-9, [0.4, 0.4, 0.4])]
for key in range(5):
    adds(INTERSECT, SR_ALL, CAT_OTHER, 1e-7, 0.0, 3, key, True, init,
         [point(1e-9, [0.35, 0.35, 0.35]), point(1e-9, [0.05, 0.05, 0.05]), point(1e-9, [0.95, 0.95, 0.95])], 'touching')
# surface pairs
for key in range(5):
    adds(INTERSECT, SR_ALL, CAT_OTHER, 1e-7, 0.0, 4, key, True, [],
         [point(1e-9, [0.5, 0.5, 1.5, 2.5], SS, 4), point(1e-9, [0.5, 0.5, 1.5, 2.5 + 1e-13], SS, 4),
          point(1e-9, [0.2, 0.9, 1.0, 2.2], SS, 4), point(1e-9, [0.2, 0.1, 1.0, 2.2], SS, 4)])
# nodes-only solutions and unknown objects
adds(INTERSECT, SR_ALL, CAT_OTHER, 1e-7, 0.0, 0, BY_ALL, True, [],
     [S(NODES, 0, [], 0.0, []), S(NODES, 0, [], 0.0, [])], 'singular')
adds(INTERSECT, SR_ALL, CAT_OTHER, 1e-7, 0.0, 1, BY_FIRST, True, [],
     [point(1e-9, [0.5], (BASE,), 1), point(1e-9, [0.5 + 1e-3], (BASE,), 1)], 'singular')
# randomised sequences
for _ in range(28):
    op = rng.choice((INTERSECT, MINIMIZE, MAXIMIZE, AT_DISTANCE, FIND, RAYFIRE, DIRECTED_MIN, PROJECTED_MAX))
    key = rng.randrange(5)
    new = []
    for _ in range(rng.randrange(3, 7)):
        base = [round(rng.random(), 1) for _ in range(3)]
        value = 10 ** rng.uniform(-9, -3)
        if rng.random() < 0.3:
            new.append(span(value, base, value, [min(1.0, b + rng.choice((0.1, 0.2))) for b in base]))
        else:
            new.append(point(value, base))
    adds(op, rng.choice((SR_ALL, SR_ALL, SR_SINGLE)), category(op), 1e-4,
         1e20 if category(op) == CAT_MIN else 0.0, 3, key, rng.random() < 0.7, [], new)

# SmSolutionArray operations.
A0, A1, A2, A3 = (point(1e-9 * (i + 1), [0.1 * (i + 1), 0.2 * (i + 1), 0.3 * (i + 1)]) for i in range(4))


ARRAY_OPS = ('add', 'ins', 'rem', 'size', 'set', 'grow', 'swap', 'dedup', 'sorted', 'append', 'copy', 'insarr', 'reset')


def array(*ops):
    tokens = [len(ops)]
    for name, *args in ops:
        tokens.extend([ARRAY_OPS.index(name)] + args)
    case('array', 'regular', *tokens)


array(('add', A0), ('add', A1), ('add', A2), ('ins', 1, 1, A3), ('ins', 0, 2, A3), ('ins', 9, 2, A0),
      ('rem', 1, 2), ('rem', 0, 1), ('size', 2), ('size', 6), ('size', 0), ('add', A1))
array(('add', A0), ('add', A1), ('add', A2), ('swap', 0, 2), ('swap', 1, 1), ('set', 1, A3), ('grow', 5, A0),
      ('grow', 0, A2), ('size', 3), ('size', 8), ('ins', 3, 1, A1))
array(('add', A0), ('add', A0), ('add', A1), ('add', A1), ('add', A0), ('dedup',), ('add', A2), ('add', A2),
      ('dedup',))
dup = lambda a0, a1, a2: point(1e-9, [a0, a1, a2])
array(('add', dup(0.5, 0.5, 0.1)), ('add', dup(0.5, 0.6, 0.2)), ('add', dup(0.5, 0.7, 0.3)),
      ('add', dup(0.4, 0.7, 0.4)), ('add', dup(0.3, 0.7, 0.5)), ('add', dup(0.5, 0.5, 0.6)), ('dedup',))
for key in range(5):
    array(('sorted', key, A2), ('sorted', key, A0), ('sorted', key, A3), ('sorted', key, A1),
          ('sorted', key, A0), ('sorted', key, point(1e20, [0.1, 0.1, 0.1])),
          ('sorted', key, point(1e-9, [0.2 + 1e-9, 0.4, 0.6])), ('sorted', key, point(5e-10, [0.2 + 1e-13, 0.4, 0.6])),
          ('sorted', key, point(1e-9, [0.2, 0.4 + 1e-3, 0.6])), ('sorted', key, point(3e-10, [0.2, 0.9, 0.6])))
pb = lambda x, value=1e-9: point(value, [x], (BASE,), 1)
for key in (BY_VALUE, BY_FIRST, BY_ALL):
    array(('sorted', key, pb(0.5)), ('sorted', key, pb(0.5 + 5e-13)), ('sorted', key, pb(0.5 + 1.5e-12)),
          ('sorted', key, pb(0.5 + 1.5e-12)), ('sorted', key, pb(0.5, 5e-10)), ('sorted', key, pb(0.5 + 1e-11, 2e-9)))
    array(('sorted', key, point(1e-9, [0.2, 0.4, 0.6])), ('sorted', key, point(1e-9, [0.2 + 1e-13, 0.4, 0.6])),
          ('sorted', key, point(-1e-9, [0.2 + 2e-13, 0.4, 0.6])), ('sorted', key, point(2e-9, [0.2 + 3e-13, 0.4, 0.6])))
array(('append', 2, A0, A1), ('append', 1, A2), ('copy', 3, A3, A2, A1), ('insarr', 1, 2, A0, A0),
      ('append', 0), ('insarr', 0, 0), ('reset',), ('add', A2))
array(('size', 3), ('add', A0), ('size', 1), ('size', 4), ('rem', 0, 4), ('ins', 0, 3, A1), ('grow', 7, A2))
array(('ins', 0, 1, A0), ('ins', 5, 1, A1), ('rem', 4, 2), ('swap', 0, 3))


# Trees.
def tree_tokens(nodes):
    out = [len(nodes)]
    for box, aux, c1, c2 in nodes:
        out += [v(*box), aux, c1, c2]
    return out


def build_tree(box, depth, aux=AUX_NONE, leaf_aux=None, early=0.0, offset=0, nodes=None, overlap=0.05):
    """A binary box tree; returns the node list and the root index."""
    nodes = [] if nodes is None else nodes
    index = len(nodes) + offset
    nodes.append(None)
    if depth == 0 or rng.random() < early:
        nodes[index - offset] = (box, aux if leaf_aux is None else leaf_aux, -1, -1)
        return nodes, index
    lo, hi = box[:3], box[3:]
    axis = max(range(3), key=lambda i: hi[i] - lo[i])
    mid = (lo[axis] + hi[axis]) / 2
    margin = overlap * (hi[axis] - lo[axis])
    b1, b2 = list(box), list(box)
    b1[3 + axis] = mid + margin
    b2[axis] = mid - margin
    _, c1 = build_tree(tuple(b1), depth - 1, aux, leaf_aux, early, offset, nodes, overlap)
    _, c2 = build_tree(tuple(b2), depth - 1, aux, leaf_aux, early, offset, nodes, overlap)
    nodes[index - offset] = (box, aux, c1, c2)
    return nodes, index


def rbox(scale=4.0, size=2.0):
    lo = [rng.uniform(-scale, scale) for _ in range(3)]
    return tuple(lo + [x + rng.uniform(0.2, size) for x in lo])


def vectors3(special=None):
    if special is not None:
        return special
    out = []
    for _ in range(3):
        a = [rng.uniform(-1, 1) for _ in range(3)]
        n = math.sqrt(sum(x * x for x in a))
        out += [x / n for x in a]
    return out


def unitv(a):
    n = math.sqrt(sum(x * x for x in a))
    return [x / n for x in a]


# Branch culling: BranchMayContainAnswers for every operation, boxes regular and at tolerance.
def box_pair(gap):
    """Two boxes whose x gap is `gap` (negative overlaps)."""
    a = (0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
    b = (1.0 + gap, 0.2, 0.2, 2.0 + gap, 0.8, 0.9)
    return a, b


TOL = 1e-3
OPS_WITH_BOX_TESTS = [o for o in range(25) if o != SURFACE_PROJECTED]
for op in OPS_WITH_BOX_TESTS:
    best = 1e20 if category(op) == CAT_MIN else (0.0 if category(op) == CAT_MAX else 0.0)
    for gap, kind in ((-0.5, 'regular'), (0.0, 'tangent'), (TOL, 'touching'), (TOL * 0.999, 'boundary'),
                      (TOL * 1.001, 'boundary'), (0.5, 'regular'), (5.0, 'regular')):
        a, b = box_pair(gap)
        vec = vectors3()
        case('branch', kind, op, 2, v(TOL, best, 1.5), v(*vec), v(*a), v(*b), v(*a))
    for kind, best, atd in (('regular', 0.5, 0.7), ('regular', 4.0, 3.0), ('regular', 1e-6, 0.1), ('boundary', 1.0, 1.0)):
        for _ in range(3):
            a, b = rbox(), rbox()
            case('branch', kind, op, 2, v(TOL, best, atd), v(*vectors3()), v(*a), v(*b), v(*rbox()))
    a = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    case('branch', 'singular', op, 2, v(0.0, 0.0, 0.0), v(*vectors3()), v(*a), v(*a), v(*a))
    case('branch', 'singular', op, 2, v(0.0, 1.0, 1.0), v(*([0.0] * 9)), v(*a), v(*rbox()), v(*a))
    case('branch', 'singular', op, 2, v(TOL, 1.0, 1.0), v(*unitv([1, 2, 3]), *unitv([1, 2, 3]), *unitv([1, 2, 3])),
         v(*a), v(0.0, 0.0, 0.0, 2.0, 2.0, 2.0), v(*a))
for best, kind in ((1.0, 'boundary'), (0.99, 'boundary'), (1.01, 'boundary')):
    a, b = box_pair(1.5)
    case('branch', kind, MINIMIZE, 2, v(0.5, best, 0.0), v(*vectors3()), v(*a), v(*b), v(*a))
wa, wb = (0.0, 0.0, 0.0, 3.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0, 4.0, 0.0)
for best, kind in ((36.0, 'boundary'), (36.01, 'boundary'), (35.99, 'boundary')):
    case('branch', kind, MAXIMIZE, 2, v(1.0, best, 0.0), v(*vectors3()), v(*wa), v(*wb), v(*wa))
for atd, tol in ((6.0, 1.0), (7.0, 2.0), (6.0, 0.99), (4.0, 1.0), (4.0, 0.99), (3.0, -1.0), (3.0, -0.99)):
    case('branch', 'boundary', AT_DISTANCE, 2, v(tol, 0.0, atd), v(*vectors3()), v(*wa), v(*wb), v(*wa))
    a, b = box_pair(4.0)
    case('branch', 'boundary', AT_DISTANCE, 2, v(tol, 0.0, atd), v(*vectors3()), v(*a), v(*b), v(*a))
for op in (INTERSECT, INTERSECTION_TEST):
    for gap, kind in ((-0.5, 'regular'), (TOL, 'touching'), (TOL * 1.001, 'boundary'), (3.0, 'regular')):
        a, b = box_pair(gap)
        c = (0.5, 0.5, 0.5, 1.5, 1.5, 1.5)
        case('branch', kind, op, 3, v(TOL, 0.0, 0.0), v(*vectors3()), v(*a), v(*c), v(*b))
        case('branch', kind, op, 3, v(TOL, 0.0, 0.0), v(*vectors3()), v(*a), v(*b), v(*c))
        case('branch', kind, op, 1, v(TOL, 0.0, 0.0), v(*vectors3()), v(*a), v(*b), v(*c))
    for _ in range(6):
        case('branch', 'regular', op, 3, v(TOL, 0.0, 0.0), v(*vectors3()), v(*rbox(1.0, 3.0)), v(*rbox(1.0, 3.0)), v(*rbox(1.0, 3.0)))
# Surface/projected intersect: only the overlapping case that SMLib answers without sampling a curve.
for gap, kind in ((-0.5, 'regular'), (0.0, 'tangent'), (TOL * 0.5, 'touching')):
    a, b = box_pair(gap)
    case('branch', kind, SURFACE_PROJECTED, 2, v(TOL, 0.0, 0.0), v(*vectors3()), v(*a), v(*b), v(*a))
# Perspective and rotated projections with a sensible eye and axes.
eye = [0.5, 0.5, -10.0]
for op in (PERSPECTIVE_INTERSECT, ROTATED_INTERSECT):
    for gap in (-0.5, 0.2, 3.0):
        a, b = box_pair(gap)
        vec = ([0.0, 0.0, 0.0] + [0.0, 0.0, 1.0] + eye) if op == PERSPECTIVE_INTERSECT else (
            [0.0, 0.0, 0.0] + [0.0, 0.0, 1.0] + [1.0, 0.0, 0.0])
        case('branch', 'regular', op, 2, v(TOL, 0.0, 0.0), v(*vec), v(*a), v(*b), v(*a))

for op in (PERSPECTIVE_INTERSECT, ROTATED_INTERSECT, PROJECTED_INTERSECT, PROJECTED_MIN, PROJECTED_MAX):
    for tol in (0.0, 0.1, 1.0, 3.0):
        for _ in range(4):
            a, b = rbox(3, 2), rbox(3, 2)
            if op == PERSPECTIVE_INTERSECT:
                vec = [0.0, 0.0, 0.0] + unitv([0.1, 0.2, 1.0]) + [rng.uniform(-1, 1), rng.uniform(-1, 1), -12.0]
            elif op == ROTATED_INTERSECT:
                vec = [rng.uniform(-1, 1) for _ in range(3)] + unitv([0.0, 0.0, 1.0]) + unitv([1.0, 0.3, 0.0])
            else:
                vec = unitv([rng.uniform(-1, 1) for _ in range(3)]) + [0.0] * 6
            best = 1e20 if category(op) == CAT_MIN else 0.0
            case('branch', 'regular', op, 2, v(tol, best, 0.0), v(*vec), v(*a), v(*b), v(*a))


# Subdivision: AreReadyForLocalSolve, FindBestSubdivisionNode and Subdivide.
for op in OPS_WITH_BOX_TESTS + [SURFACE_PROJECTED]:
    if op in (DIRECTED3D_MAX,):
        pass
    for trial in range(4):
        trees = 2
        nodes = []
        _, root1 = build_tree(rbox(2, 3), rng.randrange(1, 4), AUX_NONE, None, 0.2, 0, nodes)
        _, root2 = build_tree(rbox(2, 3), rng.randrange(1, 4), AUX_NONE, None, 0.2, 0, nodes)
        if op == RAYFIRE:
            # SMLib dereferences a child of the second tree: make it a branch
            nodes = []
            _, root1 = build_tree(rbox(2, 3), rng.randrange(0, 3), AUX_NONE, None, 0.2, 0, nodes)
            _, root2 = build_tree(rbox(2, 3), 2, AUX_NONE, None, 0.0, 0, nodes)
        branch = [root1, root2, -1]
        kind = 'regular'
        case('tree', kind, op, trees, v(TOL, 1.0), v(*vectors3()), *tree_tokens(nodes), *branch)
# Branches deeper in a tree, three trees, one tree, leaves only and surface-patch markers.
for op in (INTERSECT, MINIMIZE, MAXIMIZE, FIND):
    nodes = []
    roots = []
    for _ in range(3):
        _, r = build_tree(rbox(2, 3), 3, AUX_NONE, None, 0.1, 0, nodes)
        roots.append(r)
    inner = [i for i, n in enumerate(nodes) if n[2] >= 0]
    for trees in (3, 2, 1):
        if trees > 2 and op in (MINIMIZE, MAXIMIZE):
            continue
        branch = [rng.choice(inner), rng.choice(inner), rng.choice(inner)]
        branch = branch[:trees] + [-1] * (3 - trees)
        case('tree', 'regular', op, trees, v(TOL, 1.0), v(*vectors3()), *tree_tokens(nodes), *branch)
leaf = [((0.0, 0.0, 0.0, 1.0, 1.0, 1.0), AUX_NONE, -1, -1), ((0.5, 0.5, 0.5, 2.0, 2.0, 2.0), AUX_NONE, -1, -1)]
for op in (INTERSECT, MINIMIZE, RAYFIRE):
    case('tree', 'singular', op, 2, v(TOL, 1.0), v(*vectors3()), *tree_tokens(leaf), 0, 1, -1)
case('tree', 'singular', INTERSECT, 1, v(TOL, 1.0), v(*vectors3()), *tree_tokens(leaf), 0, -1, -1)
# Equal box sizes: the later tree wins a tie; zero-size boxes still count.
ties = [((0.0, 0.0, 0.0, 1.0, 1.0, 1.0), AUX_NONE, 2, 3), ((3.0, 3.0, 3.0, 4.0, 4.0, 4.0), AUX_NONE, 4, 5),
        ((0.0, 0.0, 0.0, 0.5, 1.0, 1.0), AUX_NONE, -1, -1), ((0.5, 0.0, 0.0, 1.0, 1.0, 1.0), AUX_NONE, -1, -1),
        ((3.0, 3.0, 3.0, 3.5, 4.0, 4.0), AUX_NONE, -1, -1), ((3.5, 3.0, 3.0, 4.0, 4.0, 4.0), AUX_NONE, -1, -1)]
for op in (INTERSECT, MINIMIZE, MAXIMIZE, DIRECTED_MIN, DIRECTED3D_MIN):
    case('tree', 'tangent', op, 2, v(TOL, 1.0), v(*vectors3()), *tree_tokens(ties), 0, 1, -1)
flat = [((0.0, 0.0, 0.0, 0.0, 0.0, 0.0), AUX_NONE, 1, 2), ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0), AUX_NONE, -1, -1),
        ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0), AUX_NONE, -1, -1)]
for op in (INTERSECT, MINIMIZE):
    case('tree', 'singular', op, 2, v(TOL, 1.0), v(*vectors3()), *tree_tokens(flat), 0, 0, -1)
# Bezier-surface markers: ready with children that are not patches, and not ready when they are.
patch = [((0.0, 0.0, 0.0, 2.0, 2.0, 1.0), AUX_SURFACE, 1, 2),
         ((0.0, 0.0, 0.0, 1.0, 2.0, 1.0), AUX_AUX, 3, 4), ((1.0, 0.0, 0.0, 2.0, 2.0, 1.0), AUX_AUX, -1, -1),
         ((0.0, 0.0, 0.0, 1.0, 1.0, 1.0), AUX_AUX, -1, -1), ((0.0, 1.0, 0.0, 1.0, 2.0, 1.0), AUX_AUX, -1, -1),
         ((0.0, 0.0, 0.0, 2.0, 2.0, 1.0), AUX_SURFACE, 6, 7),
         ((0.0, 0.0, 0.0, 1.0, 2.0, 1.0), AUX_SURFACE, -1, -1), ((1.0, 0.0, 0.0, 2.0, 2.0, 1.0), AUX_SURFACE, -1, -1),
         ((0.0, 0.0, 0.0, 3.0, 3.0, 1.0), AUX_HEAD, 9, 10),
         ((0.0, 0.0, 0.0, 1.0, 3.0, 1.0), AUX_SURFACE, 11, 12), ((1.0, 0.0, 0.0, 3.0, 3.0, 1.0), AUX_NONE, -1, -1),
         ((0.0, 0.0, 0.0, 1.0, 1.0, 1.0), AUX_NONE, -1, -1), ((0.0, 1.0, 0.0, 1.0, 3.0, 1.0), AUX_NONE, -1, -1)]
for op in (INTERSECT, MINIMIZE):
    for branch in ((0, 5), (5, 0), (0, 0), (5, 5), (8, 0), (0, 8), (9, 5), (5, 9), (2, 6), (6, 2)):
        case('tree', 'regular', op, 2, v(TOL, 1.0), v(*vectors3()), *tree_tokens(patch), *branch, -1)

# SolveIt with SmGlobalSolver's own local solve: leaf node pairs, and sums over many trees.
SOLVE_OPS = (INTERSECT, INTERSECTION_TEST, MINIMIZE, MAXIMIZE, AT_DISTANCE, FIND, NORMALIZE, DIRECTED_MIN,
             DIRECTED_MAX, PROJECTED_MIN, PROJECTED_MAX, PROJECTED_INTERSECT, SIGNED_DIRECTED_MIN,
             SIGNED_PIVOT_MIN, DIRECTED3D_MIN, WIREFRAME, ANGLE_MIN, TANGENT_THROUGH_POINT)


def solve(op, req, tol, best, vec, trees, prior, nodes, tops, kind='regular', variables=3):
    tops = list(tops) + [-1] * (3 - len(tops))
    case('solve', kind, op, req, trees, variables, v(tol, best), v(*vec), len(prior), *prior,
         *tree_tokens(nodes), *tops)


def two_trees(box1, box2, depth1, depth2, early=0.0, overlap=0.05):
    nodes = []
    _, r1 = build_tree(box1, depth1, AUX_NONE, None, early, 0, nodes, overlap)
    _, r2 = build_tree(box2, depth2, AUX_NONE, None, early, 0, nodes, overlap)
    return nodes, (r1, r2)


for op in SOLVE_OPS:
    for req in (SR_ALL, SR_SINGLE, SR_NODES):
        if req != SR_ALL and op not in (INTERSECT, MINIMIZE, MAXIMIZE):
            continue
        nodes, tops = two_trees((0.0, 0.0, 0.0, 4.0, 1.0, 1.0), (1.0, -1.0, 0.2, 3.0, 2.0, 0.8), 2, 3, 0.1)
        best = 1e20 if category(op) == CAT_MIN else 0.0
        solve(op, req, 1e-3, best, vectors3(), 2, [], nodes, tops)
        if op in (INTERSECT, MINIMIZE, MAXIMIZE, AT_DISTANCE, DIRECTED_MIN, PROJECTED_INTERSECT, SIGNED_PIVOT_MIN):
            nodes, tops = two_trees(rbox(1, 3), rbox(1, 3), 3, 2, 0.2)
            solve(op, req, 1e-2, best, vectors3(), 2, [], nodes, tops)
# Touching, tangent, disjoint and degenerate trees.
for tol, kind in ((0.0, 'tangent'), (1e-9, 'touching'), (1e-3, 'touching')):
    for gap in (0.0, 1e-4, 1e-3, 0.1):
        nodes, tops = two_trees((0.0, 0.0, 0.0, 1.0, 1.0, 1.0), (1.0 + gap, 0.0, 0.0, 2.0 + gap, 1.0, 1.0), 2, 2, 0.1, 0.0)
        solve(INTERSECT, SR_ALL, tol, 0.0, vectors3(), 2, [], nodes, tops, kind)
nodes, tops = two_trees((0.0, 0.0, 0.0, 1.0, 1.0, 0.0), (0.0, 0.0, 0.0, 1.0, 1.0, 0.0), 2, 2, 0.0, 0.0)
solve(INTERSECT, SR_ALL, 0.0, 0.0, vectors3(), 2, [], nodes, tops, 'singular')
nodes, tops = two_trees((0.0, 0.0, 0.0, 1.0, 1.0, 1.0), (5.0, 5.0, 5.0, 6.0, 6.0, 6.0), 2, 2)
solve(INTERSECT, SR_ALL, 1e-3, 0.0, vectors3(), 2, [], nodes, tops, 'singular')
nodes, tops = two_trees((0.0, 0.0, 0.0, 1.0, 1.0, 1.0), (0.5, 0.5, 0.5, 0.5, 0.5, 0.5), 2, 0)
solve(INTERSECT, SR_ALL, 0.0, 0.0, vectors3(), 2, [], nodes, tops, 'singular')
# Three trees and one tree.
nodes = []
roots = []
for box in ((0.0, 0.0, 0.0, 2.0, 2.0, 2.0), (1.0, 1.0, 1.0, 3.0, 3.0, 3.0), (1.5, 0.0, 1.5, 2.5, 2.0, 2.5)):
    _, r = build_tree(box, 2, AUX_NONE, None, 0.1, 0, nodes)
    roots.append(r)
solve(INTERSECT, SR_ALL, 1e-3, 0.0, vectors3(), 3, [], nodes, roots)
solve(INTERSECTION_TEST, SR_ALL, 1e-3, 0.0, vectors3(), 3, [], nodes, roots)
solve(FIND, SR_ALL, 1e-3, 0.0, vectors3(), 3, [], nodes, roots)
solve(INTERSECT, SR_ALL, 1e-3, 0.0, vectors3(), 1, [], nodes, roots[:1])
solve(FIND, SR_ALL, 1e-3, 0.0, vectors3(), 1, [], nodes, roots[:1])
# Prior contents of the output array stay, and SolveIt resets the best answer for a maximisation.
prior = [S(SINGLE, 3, CP, 1e-9, [0.5, 0.5, 0.5])]
nodes, tops = two_trees((0.0, 0.0, 0.0, 2.0, 1.0, 1.0), (1.0, 0.0, 0.0, 3.0, 1.0, 1.0), 2, 2)
solve(INTERSECT, SR_ALL, 1e-3, 0.0, vectors3(), 2, prior, nodes, tops)
solve(MAXIMIZE, SR_ALL, 1e-3, 1e20, vectors3(), 2, [], nodes, tops)
solve(DIRECTED_MAX, SR_ALL, 1e-3, 1e20, vectors3(), 2, [], nodes, tops)
solve(MAXIMIZE, SR_ALL, 1e-3, 1.0, vectors3(), 2, [], nodes, tops)
# Ray fire and a 3D directed minimisation with swapped traversal order.
nodes, tops = two_trees((0.0, 0.0, 0.0, 1.0, 1.0, 1.0), (-1.0, 0.2, 0.2, 5.0, 0.8, 0.8), 2, 3)
solve(RAYFIRE, SR_ALL, 1e-3, 1e20, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 2, [], nodes, tops)
solve(RAYFIRE, SR_ALL, 1e-3, 1e20, [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 2, [], nodes, tops)
solve(DIRECTED3D_MIN, SR_ALL, 1e-3, 1e20, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 2, [], nodes, tops)
solve(DIRECTED3D_MIN, SR_ALL, 1e-3, 1e20, [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 2, [], nodes, tops)
solve(RAYFIRE, SR_ALL, 1e-3, 1e20, [0.0] * 9, 2, [], nodes, tops, 'singular')


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
    header = ('# catallax-trace search gate: operation kind | inputs | SMLib outputs.\n'
              '# Generated by generate.py from oracle.cpp; see README.md.\n')
    (HERE / 'cases.txt').write_text(header + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases')


if __name__ == '__main__':
    main()
