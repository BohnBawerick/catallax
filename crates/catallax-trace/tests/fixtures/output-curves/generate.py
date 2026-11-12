#!/usr/bin/env python3
"""Regenerate the output-curve fixtures with the pinned strict SMLib driver."""
import math
import pathlib
import random
import struct
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
RNG = random.Random(2026111107)
CASES = []
NO_Z = sys.float_info.max
NO_W = -sys.float_info.max
INF_PARAM = 1234567.0


def h(value):
    return struct.pack('>d', float(value)).hex()


def add(op, kind, tokens):
    CASES.append((op, kind, ' '.join(str(t) if isinstance(t, (int, str)) else h(t) for t in tokens)))


def unit_vector():
    while True:
        v = [RNG.uniform(-1, 1) for _ in range(3)]
        n = math.sqrt(sum(x * x for x in v))
        if n > .1:
            return [x / n for x in v]


def vec(scale=3.):
    return [RNG.uniform(-scale, scale) for _ in range(3)]


# ----------------------------------------------------------------- lines

def line_tail(natural, kind, reparam=None, der=None):
    """tol t der step_t part_a part_b reparam new_min new_max, for a line whose
    NURBS interval is `natural`."""
    lo, hi = natural
    if kind == 'boundary':
        t = RNG.choice([lo, hi, lo - 1., hi + 1.5, (lo + hi) / 2])
    else:
        t = RNG.uniform(lo, hi) if hi > lo else lo
    span = hi - lo
    a, b = sorted(RNG.uniform(lo, hi) for _ in range(2)) if span > 0 else (lo, hi)
    if kind == 'boundary':
        a, b = lo, hi
    step_t = RNG.uniform(-1, 2)
    do_reparam = 1 if reparam is None else reparam
    new = sorted([RNG.uniform(-5, 5), RNG.uniform(-5, 5)])
    if new[1] - new[0] < 1e-3:
        new[1] = new[0] + 1.
    return [RNG.choice([1e-12, 1e-10, 1e-6, 1e-3, 0.]), t,
            RNG.randint(0, 3) if der is None else der, step_t, a, b,
            do_reparam, new[0], new[1]]


def line_bounded(kind, point=None, vector=None, domain=None, scale=None, dim=None, inside=None):
    point = vec() if point is None else point
    vector = unit_vector() if vector is None else vector
    domain = sorted([RNG.uniform(-2, 2), RNG.uniform(-2, 2)]) if domain is None else domain
    if domain[1] - domain[0] < 1e-3 and kind not in ('singular', 'tangent'):
        domain[1] = domain[0] + 1.
    scale = RNG.uniform(.2, 5.) if scale is None else scale
    dim = RNG.choice([2, 3]) if dim is None else dim
    inside = RNG.randint(0, 1) if inside is None else inside
    natural = (domain[0], domain[1])
    if vector is not None and math.sqrt(sum(x * x for x in vector)) <= 1e-12:
        natural = (domain[0], domain[0])
    add('line', kind, [0] + point + vector + domain + [scale, dim, inside] + line_tail(natural, kind))


for _ in range(160):
    line_bounded('regular', vector=[x * RNG.choice([1., 1., 4., .3]) for x in unit_vector()])
for _ in range(60):
    line_bounded('boundary')
for inside in (0, 1):
    for domain in ([-INF_PARAM, INF_PARAM], [0., INF_PARAM], [-INF_PARAM, 0.]):
        for _ in range(4):
            line_bounded('boundary', domain=list(domain), inside=inside)
# Directions at and around SM_EFF_ZERO, and a vanishing domain.
for length in (0., 1e-13, 5e-13, 1e-12, 1.0000001e-12, 2e-12, 1e-9):
    for _ in range(4):
        d = [x * length for x in unit_vector()]
        line_bounded('singular' if length <= 1e-12 else 'near-singular', vector=d)
for length in (1e-13, 1e-12, 1e-11, 1e-9):
    for _ in range(3):
        line_bounded('tangent', scale=length)
for _ in range(10):
    p = vec()
    line_bounded('touching', point=p, domain=[0.7, 0.7], scale=RNG.uniform(.5, 3.))
for _ in range(12):
    line_bounded('regular', domain=[RNG.uniform(-1e3, 1e3), 0.], scale=RNG.choice([1e-3, 1e3]))

SEGMENT_REPARAM = [(0., 1.), (0., 1. + 1e-13), (1e-13, 1.), (-2., 3.), (5., 6.), (0.25, 0.75), (-1e3, 1e3),
                   (2., 2. + 5e-13), (0., 1e-12), (0., 2e-12), (3., 3.), (4., 4. + 1e-9)]


def line_segment(kind, start, end, dim, interval):
    has = 0 if interval is None else 1
    iv = interval or (0., 0.)
    natural = iv if has else (0., 1.)
    add('line', kind, [1] + start + end + [dim, has] + list(iv)
        + line_tail(natural, kind, reparam=0))


for _ in range(120):
    s, e = vec(), vec()
    interval = RNG.choice([None] + SEGMENT_REPARAM)
    line_segment('regular', s, e, RNG.choice([2, 3]), interval)
for _ in range(20):
    s = vec()
    line_segment('touching', s, list(s), RNG.choice([2, 3]), RNG.choice([None, (0., 1.)]))
for length in (1e-13, 9e-13, 1e-12, 1.1e-12, 1e-11, 1e-9):
    for _ in range(3):
        s = vec()
        d = [x * length for x in unit_vector()]
        line_segment('near-singular', s, [a + b for a, b in zip(s, d)], 3, None)
for interval in SEGMENT_REPARAM:
    for _ in range(3):
        line_segment('boundary', vec(), vec(), 3, interval)
for _ in range(10):
    line_segment('regular', [0., 0., 0.], [1e6, -1e6, 5e5], 3, RNG.choice([None, (-1., 1.)]))
for _ in range(10):
    line_segment('regular', [1e-6, 0., 0.], [2e-6, 1e-6, 0.], 2, RNG.choice([None, (-1., 1.)]))

# --------------------------------------------------------------- hermite


def hermite(kind, p1=None, d1=None, p2=None, d2=None, two_d=None, t=None, der=None, nonzero=None):
    p1 = vec() if p1 is None else p1
    p2 = vec() if p2 is None else p2
    d1 = vec(4.) if d1 is None else d1
    d2 = vec(4.) if d2 is None else d2
    two_d = RNG.randint(0, 1) if two_d is None else two_d
    if t is None:
        t = RNG.choice([0., 1., .5, RNG.uniform(0, 1), RNG.uniform(-.5, 1.5)]) if kind == 'boundary' \
            else RNG.uniform(0, 1)
    add('hermite', kind, [two_d] + p1 + d1 + p2 + d2
        + [t, RNG.randint(0, 3) if der is None else der, RNG.randint(0, 1) if nonzero is None else nonzero])


for _ in range(160):
    hermite('regular')
for _ in range(60):
    hermite('boundary')
# Zero end derivatives: the first derivative vanishes at an end.
for which in ('start', 'end', 'both'):
    for t in (0., 1e-9, .3, .5, .7, 1.):
        for der in (1, 2, 3):
            z = [0., 0., 0.]
            hermite('tangent', d1=z if which != 'end' else None, d2=z if which != 'start' else None,
                    t=t, der=der, nonzero=1)
# A first derivative that vanishes mid-interval: the replacement flips sign
# after the midpoint, and not at it.
for _ in range(12):
    p1, p2, d1 = vec(), vec(), vec(4.)
    d2 = [6 * (b - a) - c for a, b, c in zip(p1, p2, d1)]
    hermite('tangent', p1=p1, d1=d1, p2=p2, d2=d2, t=.5, der=RNG.randint(2, 3), nonzero=1)
    hermite('tangent', p1=p1, d1=d1, p2=p2, d2=d2, t=.5 + 1e-9, der=2, nonzero=1)
    hermite('tangent', p1=p1, d1=d1, p2=p2, d2=d2, t=.5 - 1e-9, der=2, nonzero=1)
# A start whose first derivative vanishes and whose second is a few scaled zeros.
for scale in (1.9e-13, 2.5e-13, 3e-13, 3.6e-13, 3.8e-13):
    for t in (0., 1.):
        delta = [x * scale for x in unit_vector()]
        hermite('near-singular', p1=[0., 0., 0.], d1=[0., 0., 0.], p2=delta, d2=[0., 0., 0.],
                t=t, der=2, nonzero=1, two_d=0)
# Closed and collapsed curves.
for _ in range(20):
    p = vec()
    hermite('touching', p1=p, p2=list(p))
for _ in range(10):
    p = vec()
    z = [0., 0., 0.]
    hermite('singular', p1=p, p2=list(p), d1=z, d2=z, t=RNG.choice([0., .3, 1.]))
for length in (1e-14, 1e-13, 9e-13, 1e-12, 1.2e-12, 2e-12, 1e-10):
    for t in (0., .2, .8, 1.):
        d = [x * length for x in unit_vector()]
        hermite('near-singular', d1=d, d2=[x * length for x in unit_vector()], t=t,
                der=RNG.randint(1, 3), nonzero=1)

# ------------------------------------------------------------------ clip


def clip(kind, lo, hi, point, vector, param):
    add('clip', kind, lo + hi + point + vector + [param])


for _ in range(80):
    lo = [RNG.uniform(-2, 0), RNG.uniform(-2, 0)]
    hi = [RNG.uniform(.5, 3), RNG.uniform(.5, 3)]
    point = [RNG.uniform(lo[0], hi[0]), RNG.uniform(lo[1], hi[1])]
    clip('regular', lo, hi, point, [RNG.uniform(-3, 3), RNG.uniform(-3, 3)], RNG.choice([1., -1., .5, 4., -6.]))
for _ in range(40):
    lo, hi = [0., 0.], [1., 1.]
    # A point on a side or a corner, with the vector pointing in, out or along.
    point = RNG.choice([[0., .3], [1., .6], [.4, 0.], [.2, 1.], [0., 0.], [1., 1.], [0., 1.], [1., 0.],
                        [1e-13, .5], [1. - 1e-13, .5]])
    clip('touching', lo, hi, point, RNG.choice([[1., 0.], [-1., 0.], [0., 1.], [0., -1.], [1., 1.], [-1., -1.],
                                                [1e-13, 1.], [-1e-13, -1.], [0., 0.]]),
         RNG.choice([1., -1., 2., 1e-3]))
for _ in range(30):
    lo, hi = [0., 0.], [1., 1.]
    point = [RNG.uniform(.1, .9), RNG.uniform(.1, .9)]
    # A vector that reaches past each side, through both, or only just reaches.
    clip('boundary', lo, hi, point, RNG.choice([[2., 0.], [-2., 0.], [0., 2.], [0., -2.], [2., 2.], [-2., 3.],
                                                [.5, 0.], [0., -.5]]), 1.)
for _ in range(20):
    lo, hi = [0., 0.], [1., 1.]
    # A start outside the box with a vector that cannot reach a side.
    point = RNG.choice([[-1., .5], [2., .5], [.5, -1.], [.5, 2.]])
    clip('singular', lo, hi, point, RNG.choice([[1e-14, 0.], [0., 1e-14], [1e-14, 1e-14], [0., 0.]]), 1.)
for _ in range(10):
    lo, hi = [0., 0.], [1., 1.]
    clip('near-singular', lo, hi, [.5, .5], [1., 1.], 1.00000001 * RNG.choice([1., 2., .5]))

# ----------------------------------------------------------- spline data


def clamped_knots(degree, count, lo=0., hi=None, interior=None):
    inner = count - degree - 1
    hi = float(inner + 1) if hi is None else hi
    if interior is None:
        interior = sorted(RNG.uniform(lo, hi) for _ in range(inner))
    return [lo] * (degree + 1) + list(interior) + [hi] * (degree + 1)


def spline_points(count, dim, rational, kind='regular'):
    points = [[RNG.uniform(-4, 4), RNG.uniform(-4, 4), RNG.uniform(-4, 4) if dim == 3 else NO_Z]
              for _ in range(count)]
    if kind == 'singular':
        points = [[2., 1., 3. if dim == 3 else NO_Z]] * count
    elif kind == 'tangent':
        points[1] = list(points[0])
    elif kind == 'touching':
        points[-1] = list(points[0])
    weights = [RNG.uniform(.3, 3.) for _ in range(count)] if rational else None
    return points, weights


def stored_tokens(degree, knots, points, weights):
    """Input tokens `degree count knots... (x y z w)...` for a stored curve."""
    tokens = [degree, len(points)] + list(knots)
    for i, point in enumerate(points):
        if weights is None:
            tokens += [point[0], point[1], point[2], NO_W]
        else:
            w = weights[i]
            tokens += [point[0] * w, point[1] * w, point[2] * w if point[2] != NO_Z else NO_Z, w]
    return tokens


def interior_knots(inner, degree, lo, hi, repeated):
    """`inner` interior knots in (lo, hi), with runs of up to `degree` equal knots."""
    runs, remaining = [], inner
    while remaining > 0:
        run = min(remaining, RNG.randint(1, degree) if repeated else 1)
        runs.append(run)
        remaining -= run
    step = (hi - lo) / (len(runs) + 1)
    values = [lo + step * (i + 1) + RNG.uniform(-.3, .3) * step for i in range(len(runs))]
    return [v for v, run in zip(values, runs) for _ in range(run)]


def curve_inputs(kind='regular', allow_rational=True, max_degree=4):
    degree = RNG.randint(1, max_degree)
    count = degree + 1 + RNG.choice([0, 0, 1, 2, 3, 5])
    rational = allow_rational and RNG.random() < .4
    points, weights = spline_points(count, RNG.choice([2, 3]), rational, kind)
    lo = RNG.choice([0., 0., -1., 2.])
    hi = lo + RNG.uniform(1., 4.)
    knots = clamped_knots(degree, count, lo, hi,
                          interior_knots(count - degree - 1, degree, lo, hi, RNG.random() < .3))
    return degree, knots, points, weights


# ---------------------------------------------------------------- reparam


def reparam(kind, degree, knots, points, weights, new):
    add('reparam', kind, stored_tokens(degree, knots, points, weights) + list(new))


for _ in range(80):
    degree, knots, points, weights = curve_inputs()
    reparam('regular', degree, knots, points, weights, sorted([RNG.uniform(-5, 5), RNG.uniform(6, 9)]))
for _ in range(20):
    degree, knots, points, weights = curve_inputs()
    # Equal to, and within the scaled zero of, the current range: no change.
    lo, hi = knots[0], knots[-1]
    reparam('boundary', degree, knots, points, weights,
            RNG.choice([(lo, hi), (lo + 1e-13, hi), (lo, hi + 1e-13), (lo + 1e-11, hi), (lo, hi + 2e-12),
                        (lo - 1e-13, hi - 1e-13)]))
for kind in ('singular', 'tangent', 'touching'):
    for _ in range(8):
        degree, knots, points, weights = curve_inputs(kind)
        reparam(kind, degree, knots, points, weights, [RNG.uniform(-3, 0), RNG.uniform(1, 3)])
for _ in range(10):
    degree, knots, points, weights = curve_inputs()
    reparam('near-singular', degree, knots, points, weights, [3., 3. + RNG.choice([1e-9, 1e-6, 1e-3])])
for _ in range(10):
    degree, knots, points, weights = curve_inputs()
    reparam('regular', degree, knots, points, weights, [-1e6, 1e6])

# ------------------------------------------------------------------ trim


def trim(kind, degree, knots, points, weights, a, b, mode=0):
    add('trim', kind, stored_tokens(degree, knots, points, weights) + [a, b, mode])


def distinct(knots):
    out = []
    for k in knots:
        if not out or k != out[-1]:
            out.append(k)
    return out


for _ in range(220):
    degree, knots, points, weights = curve_inputs(max_degree=4)
    lo, hi = knots[0], knots[-1]
    a, b = sorted(RNG.uniform(lo, hi) for _ in range(2))
    if b - a < .05:
        b = min(hi, a + .3)
    trim('regular', degree, knots, points, weights, a, b, 1 if RNG.random() < .1 else 0)
def near_knot(knot):
    """A knot, or a value at a small offset or at a fraction of the trim's
    snapping tolerance (100 scaled zeros) from it."""
    pick = RNG.random()
    if pick < .3:
        return knot
    if pick < .5:
        return knot + RNG.choice([1e-13, -1e-13, 1e-11, 5e-11])
    return knot + RNG.choice([.5, .95, 1.05, 2.]) * 100e-12 * (1 + abs(knot)) * RNG.choice([1, -1])


for _ in range(160):
    # Ends at or near existing knots.
    degree, knots, points, weights = curve_inputs(max_degree=4)
    uniq = distinct(knots)
    lo, hi = uniq[0], uniq[-1]
    a = near_knot(RNG.choice(uniq[:-1]))
    b = near_knot(RNG.choice(uniq[1:]))
    if b <= a:
        a, b = lo, hi
    trim('touching', degree, knots, points, weights, a, b, 1 if RNG.random() < .15 else 0)
for _ in range(40):
    # Cuts that cover the whole curve or only barely trim it.
    degree, knots, points, weights = curve_inputs(max_degree=4)
    lo, hi = knots[0], knots[-1]
    a, b = RNG.choice([(lo, hi), (lo - 1., hi + 1.), (lo - 1e-13, hi), (lo, hi + 1e-13), (lo + 1e-13, hi - 1e-13),
                       (lo + 1e-11, hi), (lo, hi - 1e-11), (lo + 1e-9, hi), (lo, hi - 1e-9), (lo + 1e-7, hi)])
    trim('boundary', degree, knots, points, weights, a, b)
for _ in range(60):
    # Slivers: an end just past a knot, so the cut leaves a very short span.
    degree, knots, points, weights = curve_inputs(max_degree=4)
    uniq = distinct(knots)
    if len(uniq) < 3:
        continue
    i = RNG.randint(1, len(uniq) - 2)
    sliver = RNG.choice([1e-9, 1e-8, 1e-7, 1e-6, 4e-6, 1e-5, 1e-4])
    if RNG.random() < .5:
        a, b = uniq[0], uniq[i] + sliver
    else:
        a, b = uniq[i] - sliver, uniq[-1]
    trim('near-singular', degree, knots, points, weights, a, b)
for _ in range(40):
    # A single span trimmed to a short piece, which SMLib leaves alone.
    degree = RNG.randint(1, 4)
    points, weights = spline_points(degree + 1, RNG.choice([2, 3]), RNG.random() < .4)
    lo = RNG.choice([0., 1.])
    hi = lo + RNG.uniform(1., 3.)
    a = RNG.uniform(lo + .1, hi - .1)
    width = RNG.choice([1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-3])
    trim('near-singular', degree, clamped_knots(degree, degree + 1, lo, hi, []), points, weights, a, a + width)
for kind in ('singular', 'tangent'):
    for _ in range(15):
        degree, knots, points, weights = curve_inputs(kind, max_degree=4)
        lo, hi = knots[0], knots[-1]
        a, b = sorted(RNG.uniform(lo, hi) for _ in range(2))
        if b - a < .05:
            b = min(hi, a + .3)
        trim(kind, degree, knots, points, weights, a, b)
for _ in range(25):
    # Intervals outside the curve or of no length: refused.
    degree, knots, points, weights = curve_inputs(max_degree=4)
    lo, hi = knots[0], knots[-1]
    mid = (lo + hi) / 2
    a, b = RNG.choice([(lo - 1., mid), (mid, hi + 1.), (mid, mid), (mid, mid + 1e-13), (lo - 1e-9, mid),
                       (mid, hi + 1e-9), (lo - 1e-12, hi - .1)])
    trim('singular', degree, knots, points, weights, a, b, RNG.randint(0, 1))

# ----------------------------------------------------------------- canon


def canon(kind, dim, degree, points, knots, mults, weights=None, trim_interval=None):
    tokens = [dim, degree, len(points), 0 if weights is None else 1, 0 if trim_interval is None else 1]
    if trim_interval is not None:
        tokens += list(trim_interval)
    tokens.append(len(knots))
    for point in points:
        tokens += point
    if weights is not None:
        tokens += weights
    for k, m in zip(knots, mults):
        tokens += [k, m]
    add('canon', kind, tokens)


def knot_pattern(degree, count, interior_runs, end_extra=(0, 0), lo=0., hi=None):
    """Unique knots and multiplicities for `count` points of `degree`, with the
    given interior runs and extra end multiplicity, padding with a final
    interior run so the knots total count + degree + 1."""
    total = count + degree + 1
    first, last = degree + 1 + end_extra[0], degree + 1 + end_extra[1]
    runs = list(interior_runs)
    used = first + last + sum(runs)
    if used > total:
        return None
    if used < total:
        runs.append(total - used)
    hi = float(len(runs) + 1) if hi is None else hi
    uniq = [lo] + [lo + (hi - lo) * (i + 1) / (len(runs) + 1) for i in range(len(runs))] + [hi]
    return uniq, [first] + runs + [last]


def canon_points(count, dim):
    return [[RNG.uniform(-4, 4), RNG.uniform(-4, 4), RNG.uniform(-4, 4) if dim == 3 else NO_Z]
            for _ in range(count)]


for _ in range(160):
    degree = RNG.randint(1, 5)
    count = degree + 1 + RNG.choice([0, 1, 2, 3, 6])
    pattern = None
    while pattern is None:
        runs = [RNG.randint(1, degree) for _ in range(RNG.randint(0, max(0, count - degree - 1)))]
        pattern = knot_pattern(degree, count, runs, lo=RNG.choice([0., -2.]), hi=RNG.choice([1., 5.]))
    dim = RNG.choice([1, 2, 3, 3])
    pdim = 2 if dim == 2 else 3
    weights = None
    if RNG.random() < .4:
        weights = [RNG.uniform(.2, 5.) for _ in range(count)]
    canon('regular', dim, degree, canon_points(count, pdim), pattern[0], pattern[1], weights)
# Planar points asked for in three dimensions, and the reverse.
for _ in range(20):
    degree = RNG.randint(1, 4)
    count = degree + 1 + RNG.randint(0, 3)
    pattern = knot_pattern(degree, count, [])
    canon('boundary', 3, degree, canon_points(count, 2), pattern[0], pattern[1])
    canon('boundary', 2, degree, canon_points(count, 3), pattern[0], pattern[1])
# Interior multiplicities above the degree, and end multiplicities above degree + 1.
for _ in range(80):
    degree = RNG.randint(1, 4)
    extra_first, extra_last = RNG.choice([(0, 0), (1, 0), (0, 1), (2, 1), (1, 1)])
    run = degree + RNG.randint(1, 2)
    count = degree + 1 + extra_first + extra_last + run + RNG.randint(0, 3)
    pattern = knot_pattern(degree, count, [run], (extra_first, extra_last))
    if pattern is None:
        continue
    weights = [RNG.uniform(.5, 2.) for _ in range(count)] if RNG.random() < .3 else None
    canon('touching', 3, degree, canon_points(count, 3), pattern[0], pattern[1], weights)
# Unclamped knots, with the interval trimmed to the clamped range or inside it.
for _ in range(60):
    degree = RNG.randint(1, 4)
    count = degree + 1 + RNG.randint(1, 5)
    total = count + degree + 1
    uniq = [float(i) for i in range(total)]
    base = canon_points(count, 3)
    bounds = RNG.choice([None, (uniq[degree], uniq[total - degree - 1]),
                         (uniq[degree] + .3, uniq[total - degree - 1] - .2)])
    canon('regular', 3, degree, base, uniq, [1] * total, None, bounds)
# Junk end knots: a first knot below, or a last knot above, the bounds that
# repeat degree times.
for _ in range(40):
    degree = RNG.randint(1, 4)
    junk_first, junk_last = RNG.choice([(1, 0), (0, 1), (1, 1)])
    interior = [RNG.randint(1, degree) for _ in range(RNG.randint(1, 3))]
    first = [1, degree] if junk_first else [degree + 1]
    last = [degree, 1] if junk_last else [degree + 1]
    mults = first + interior + last
    uniq = ([-.5] if junk_first else []) + [0.] + [float(i + 1) for i in range(len(interior))]
    uniq += [float(len(interior) + 1)] + ([len(interior) + 1.5] if junk_last else [])
    count = sum(mults) - degree - 1
    weights = [RNG.uniform(.5, 2.) for _ in range(count)] if RNG.random() < .3 else None
    canon('touching', 3, degree, canon_points(count, 3), uniq, mults, weights)
# Trim intervals outside the knot range and reversed or empty weights/degrees: refused.
for _ in range(30):
    degree = RNG.randint(1, 4)
    count = degree + 1 + RNG.randint(0, 3)
    pattern = knot_pattern(degree, count, [])
    lo, hi = pattern[0][0], pattern[0][-1]
    canon('singular', 3, degree, canon_points(count, 3), pattern[0], pattern[1], None,
          RNG.choice([(lo - 1., hi), (lo, hi + 1.), (lo - 1e-9, hi), (hi, hi), (lo, lo)]))
for degree in (17, 18, 20, 32, 33, 0):
    count = degree + 2
    pattern = knot_pattern(degree, count, [])
    canon('regular' if degree <= 17 else 'singular', 3, max(degree, 0), canon_points(count, 3), pattern[0], pattern[1])
for weight in (1e-4, 9e-4, 1e-3, 100., 100.5, 200., 1e-9):
    for _ in range(2):
        degree, count = 2, 5
        pattern = knot_pattern(degree, count, [])
        weights = [RNG.uniform(1., 2.) for _ in range(count)]
        weights[RNG.randint(0, count - 1)] = weight
        canon('boundary' if 1e-3 <= weight <= 100. else 'singular', 3, degree, canon_points(count, 3),
              pattern[0], pattern[1], weights)
for dim in (0, 4, 7):
    degree, count = 2, 4
    pattern = knot_pattern(degree, count, [])
    canon('singular', dim, degree, canon_points(count, 3), pattern[0], pattern[1])
for _ in range(10):
    # All weights equal: a non-rational curve with unscaled points.
    degree, count = 3, 6
    pattern = knot_pattern(degree, count, [1, 1])
    w = RNG.uniform(.5, 4.)
    canon('tangent', 3, degree, canon_points(count, 3), pattern[0], pattern[1], [w] * count)
for _ in range(10):
    degree, count = RNG.randint(1, 3), 5
    pattern = knot_pattern(degree, count, [])
    pts = canon_points(count, 3)
    pts[1] = list(pts[0])
    pts[-1] = list(pts[0])
    canon('touching', 3, degree, pts, pattern[0], pattern[1])
# Knot spacing at the scaled zero, which the constructor merges.
for gap in (1e-14, 5e-13, 1e-12, 2e-12, 1e-11, 1e-9):
    for end in (0., 1., 1000.):
        degree, count = 2, 5
        pattern = knot_pattern(degree, count, [1, 1], lo=end, hi=end + 3.)
        uniq = list(pattern[0])
        uniq[2] = uniq[1] + gap
        canon('near-singular', 3, degree, canon_points(count, 3), uniq, pattern[1])


# ------------------------------------------------------------------ main


# ------------------------------------------------------------ short spans
# Cut curves that keep short segments (no longer than the zone tolerance) for
# SMLib's N_CrvRemoveDegenSegs. A separate generator keeps the cases above
# unchanged.

SHORT = random.Random(2026111220)
ZONE = 1e-5


def short_curve(degree, spans, runs, rational, dim, scale=4.):
    """A clamped curve of `spans` spans whose runs `(first span, span count,
    size)` have their control points within `size` of one point (zero for
    coincident points)."""
    count = degree + spans
    points = [[SHORT.uniform(-scale, scale) for _ in range(2)] + [SHORT.uniform(-scale, scale) if dim == 3 else NO_Z]
              for _ in range(count)]
    for first, length, size in runs:
        base = [SHORT.uniform(-scale, scale) for _ in range(3)]
        for index in range(first, first + length + degree):
            offset = [SHORT.uniform(-1, 1) * size for _ in range(3)]
            points[index] = [base[0] + offset[0], base[1] + offset[1], base[2] + offset[2] if dim == 3 else NO_Z]
    lo = SHORT.choice([0., -1., 2.])
    hi = lo + SHORT.uniform(1., 4.)
    inner = sorted(SHORT.uniform(lo, hi) for _ in range(spans - 1))
    knots = [lo] * (degree + 1) + inner + [hi] * (degree + 1)
    weights = [SHORT.uniform(.3, 3.) for _ in range(count)] if rational else None
    return knots, points, weights


def sliver_interval(knots, ends):
    """A trim leaving a sliver of a span at one or both ends."""
    uniq = distinct(knots)
    sliver = lambda: SHORT.choice([1e-9, 1e-8, 1e-7, 5e-7, 1e-6, 2e-6])
    a, b = uniq[0], uniq[-1]
    if ends in ('start', 'both'):
        a = uniq[SHORT.randint(1, max(1, len(uniq) // 2 - 1))] - sliver()
    if ends in ('end', 'both'):
        b = uniq[SHORT.randint(len(uniq) // 2 + 1, len(uniq) - 2)] + sliver()
    return a, b


def short_trim(kind, degree, knots, points, weights, ends):
    a, b = sliver_interval(knots, ends)
    trim(kind, degree, knots, points, weights, a, b)


# Interior runs of one to three short spans, coincident (singular), near
# the tolerance (touching) or well within it (near-singular), with slivers
# at the ends.
for kind, sizes in (('singular', [0.]), ('touching', [1e-6, 2e-6, 3e-6, 4e-6]),
                    ('near-singular', [1e-9, 1e-8, 1e-7])):
    for _ in range(40):
        degree = SHORT.randint(1, 4)
        spans = SHORT.randint(6, 10)
        runs = []
        first = SHORT.randint(2, 3)
        while first + 3 < spans - 2 and len(runs) < 2:
            length = SHORT.randint(1, 3)
            runs.append((first, length, SHORT.choice(sizes)))
            first += length + SHORT.randint(2, 3)
        knots, points, weights = short_curve(degree, spans, runs, SHORT.random() < .4, SHORT.choice([2, 3]))
        short_trim(kind, degree, knots, points, weights, SHORT.choice(['start', 'end', 'both']))

# Short runs at the curve's own ends, kept by a sliver at the other end, so
# the run is reparameterised at the start or the end.
for _ in range(30):
    degree = SHORT.randint(1, 4)
    spans = SHORT.randint(5, 8)
    at_start = SHORT.random() < .5
    length = SHORT.randint(1, 2)
    first = 0 if at_start else spans - length
    knots, points, weights = short_curve(degree, spans, [(first, length, SHORT.choice([0., 1e-7, 2e-6]))],
                                         SHORT.random() < .4, SHORT.choice([2, 3]))
    short_trim('boundary', degree, knots, points, weights, 'end' if at_start else 'start')

# Short spans shaped as lines and circular arcs, which the arc length
# measures exactly: collinear points, and rational quadratic Bezier arcs.
for _ in range(15):
    degree = SHORT.randint(1, 3)
    spans = SHORT.randint(6, 9)
    knots, points, weights = short_curve(degree, spans, [], False, 3)
    first = SHORT.randint(2, spans - 4)
    base, direction = vec(), unit_vector()
    for k, index in enumerate(range(first, first + 1 + degree)):
        step = k * SHORT.choice([1e-6, 2e-6]) / max(degree, 1)
        points[index] = [base[i] + step * direction[i] for i in range(3)]
    short_trim('tangent', degree, knots, points, weights, SHORT.choice(['start', 'end']))
for _ in range(15):
    spans = SHORT.randint(5, 8)
    lo = 0.
    breaks = sorted(SHORT.uniform(lo + .2, 3.8) for _ in range(spans - 1))
    knots = [lo] * 3 + [b for b in breaks for _ in range(2)] + [4.] * 3
    count = len(knots) - 3
    points = [[SHORT.uniform(-4, 4) for _ in range(3)] for _ in range(count)]
    weights = [SHORT.uniform(.5, 2.) for _ in range(count)]
    span = SHORT.randint(2, spans - 3)
    radius = SHORT.choice([1e-6, 2e-6, 3e-6])
    sweep = SHORT.choice([math.pi / 3, math.pi / 2, 2 * math.pi / 3])
    center, x_axis = vec(), unit_vector()
    y_axis = [x_axis[1] - x_axis[2], x_axis[2] - x_axis[0], x_axis[0] - x_axis[1]]
    n = math.sqrt(sum(c * c for c in y_axis))
    y_axis = [c / n for c in y_axis]
    dot = sum(a * b for a, b in zip(x_axis, y_axis))
    y_axis = [y - dot * x for x, y in zip(x_axis, y_axis)]
    n = math.sqrt(sum(c * c for c in y_axis))
    y_axis = [c / n for c in y_axis]
    at = lambda angle, r: [center[i] + r * (math.cos(angle) * x_axis[i] + math.sin(angle) * y_axis[i]) for i in range(3)]
    i0 = 2 * span
    points[i0] = at(0., radius)
    points[i0 + 1] = at(sweep / 2, radius / math.cos(sweep / 2))
    points[i0 + 2] = at(sweep, radius)
    weights[i0], weights[i0 + 1], weights[i0 + 2] = 1., math.cos(sweep / 2), 1.
    short_trim('tangent', 2, knots, points, weights, SHORT.choice(['start', 'end']))

# Whole cut curves within the tolerance: one merged short segment, which
# SMLib leaves alone.
for _ in range(10):
    degree = SHORT.randint(1, 3)
    spans = SHORT.randint(4, 6)
    knots, points, weights = short_curve(degree, spans, [], SHORT.random() < .4, 3, scale=SHORT.choice([2e-7, 5e-7]))
    short_trim('singular', degree, knots, points, weights, SHORT.choice(['start', 'end']))
# Whole cut curves of several runs of short segments, longer together than
# the tolerance: SMLib removes their interior knots with N_CrvRemoveKnots,
# which the port refuses.
for _ in range(12):
    degree = SHORT.randint(1, 3)
    spans = SHORT.randint(5, 7)
    knots, points, weights = short_curve(degree, spans, [], SHORT.random() < .4, 3)
    base, direction = vec(), unit_vector()
    step = SHORT.choice([2e-6, 3e-6, 4e-6])
    for index in range(len(points)):
        points[index] = [base[i] + index * step * direction[i] + SHORT.uniform(-1, 1) * 1e-7 for i in range(3)]
    short_trim('touching', degree, knots, points, weights, SHORT.choice(['start', 'end']))


def main():
    oracle = sys.argv[1]
    raw = pathlib.Path(sys.argv[2])
    raw.mkdir(parents=True, exist_ok=True)

    def ask(request):
        run = subprocess.run(['prlimit', '--as=1073741824', '--cpu=10', 'timeout', '20', oracle],
                             input=request, capture_output=True, text=True)
        return run

    consts = ask('consts\n')
    assert consts.stdout.split()[1] == 'd' + h(1e-5), 'zone tolerance is not 1e-5'
    lines = []
    refused = short = 0
    for index, (op, kind, inputs) in enumerate(CASES):
        request = f'{op} {inputs}\n'
        run = ask(request)
        again = ask(request)
        (raw / f'{index:04d}-{op}.txt').write_text(
            f'request: {request}exit: {run.returncode}\nstdout: {run.stdout}stderr: {run.stderr}')
        if run.returncode or again.returncode:
            sys.exit(f'case {index} {op}: oracle exit {run.returncode}: {run.stderr}')
        if run.stdout != again.stdout:
            sys.exit(f'case {index} {op}: oracle output is not deterministic')
        response = run.stdout.split()
        assert response[0] == op
        outputs = response[1:]
        op_name = op
        if op == 'trim' and len(outputs) > 3 and outputs[0] == 'i1000' and outputs[3] == 'i1':
            # SMLib entered its short-span branch (N_CrvRemoveDegenSegs).
            whole = outputs.pop(4)
            if whole == 'i1':
                # One run of short segments covers the cut curve, so SMLib
                # calls N_CrvRemoveKnots, which is not ported: the port
                # refuses, and only the status, interval and flag are kept.
                outputs = outputs[:4]
                op_name = 'trim-whole'
                refused += 1
            else:
                op_name = 'trim-short'
                short += 1
        lines.append(f'{op_name} {kind} | {inputs} | {" ".join(outputs)}')
    (HERE / 'cases.txt').write_text(
        '# Generated by generate.py from the pinned strict SMLib; see README.md.\n' + '\n'.join(lines) + '\n')
    print(f'{len(lines)} cases, {short} short-span trims, {refused} whole-curve refusals')


if __name__ == '__main__':
    main()
