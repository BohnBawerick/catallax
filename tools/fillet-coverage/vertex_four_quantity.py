#!/usr/bin/env python3
"""Bound the volume of vertex-4-all-unequal-slopes independently of CAD output.

Requires mpmath 1.3.0. Run with one output JSON path. All surface fluxes have
closed-form antiderivatives; no quadrature or kernel-derived quantities enter.
"""

import argparse
import json
from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal, localcontext
from fractions import Fraction
from pathlib import Path

from mpmath import iv


def q(a, b=1):
    return iv.mpf(a) / b


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def scale(a, s):
    return [x * s for x in a]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def norm(a):
    return iv.sqrt(dot(a, a))


def unit(a):
    return scale(a, 1 / norm(a))


def angle(a, b):
    return iv.atan2(norm(cross(a, b)), dot(a, b))


def area(points, normal):
    twice = [q(0), q(0), q(0)]
    for a, b in zip(points, points[1:] + points[:1]):
        twice = add(twice, cross(a, b))
    return abs(dot(twice, normal)) / 2


def base_foot(centre, normal, axis, radius):
    foot = add(centre, scale(normal, radius))
    return add(foot, scale(axis, -foot[2] / axis[2]))


def incoming_flux(centre, a, b, axis, r):
    cosine = dot(a, b)
    sine = norm(cross(a, b))
    alpha = iv.atan2(sine, cosine)
    w = scale(sub(b, scale(a, cosine)), 1 / sine)
    integral_cos = sine
    integral_sin = 1 - cosine
    integral_cos2 = alpha / 2 + sine * cosine / 2
    integral_sin2 = alpha / 2 - sine * cosine / 2
    integral_sincos = sine * sine / 2
    u, v = dot(centre, a), dot(centre, w)
    linear = centre[2] * (u * integral_cos + v * integral_sin + r * alpha)
    quadratic = r * (u * a[2] * integral_cos2 + (u * w[2] + v * a[2]) * integral_sincos + v * w[2] * integral_sin2)
    radial = r * r * (a[2] * integral_cos + w[2] * integral_sin)
    return r * (linear + quadratic + radial) / (-axis[2])


def bridge_flux(centre, a, b, length, r):
    cosine = dot(a, b)
    sine = norm(cross(a, b))
    alpha = iv.atan2(sine, cosine)
    w = scale(sub(b, scale(a, cosine)), 1 / sine)
    return r * length * (dot(centre, a) * sine + dot(centre, w) * (1 - cosine) + r * alpha)


def sphere_flux(centre, normals, r):
    a, b, c = normals
    omega = 2 * iv.atan2(dot(a, cross(b, c)), 1 + dot(a, b) + dot(b, c) + dot(c, a))
    moment = [q(0), q(0), q(0)]
    for p, z in [(a, b), (b, c), (c, a)]:
        moment = add(moment, scale(unit(cross(p, z)), angle(p, z) / 2))
    return r * r * dot(centre, moment) + r * r * r * omega


def exact_endpoint(bits):
    sign, mantissa, exponent, count = bits
    assert count >= 0
    value = Fraction(-mantissa if sign else mantissa)
    return value * 2**exponent if exponent >= 0 else value / Fraction(2 ** (-exponent))


def decimal_endpoint(value, rounding):
    with localcontext() as context:
        context.prec = 75
        context.rounding = rounding
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def interval_record(value):
    lo, hi = map(exact_endpoint, value._mpi_)
    assert lo <= hi
    return {
        "lower": decimal_endpoint(lo, ROUND_FLOOR),
        "upper": decimal_endpoint(hi, ROUND_CEILING),
        "width_upper": decimal_endpoint(hi - lo, ROUND_CEILING),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    iv.dps = 70
    h, r, k = q(10), q(1, 2), q(1, 2)
    sqrt2 = iv.sqrt(2)
    sqrt_d = iv.sqrt(1 + k * k)
    a = [1 / sqrt2, q(0), 1 / sqrt2]
    b = [q(0), 1 / sqrt2, 1 / sqrt2]
    d = [q(0), -1 / sqrt_d, k / sqrt_d]
    t = r * (sqrt_d - k * sqrt2) / (1 + k)
    z = h - r * (sqrt2 + sqrt_d) / (1 + k)
    right, left = [t, t, z], [-t, t, z]
    u_ab = unit([q(1), q(1), q(-1)])
    u_ad = unit([q(1), -k, q(-1)])
    foot_a = add(right, scale(a, r))
    foot_b_right, foot_b_left = add(right, scale(b, r)), add(left, scale(b, r))
    foot_d_right, foot_d_left = add(right, scale(d, r)), add(left, scale(d, r))
    base_a_ab, base_a_ad = base_foot(right, a, u_ab, r), base_foot(right, a, u_ad, r)
    base_b_right = base_foot(right, b, u_ab, r)
    base_b_left = [-base_b_right[0], base_b_right[1], base_b_right[2]]
    base_d_right = base_foot(right, d, u_ad, r)
    base_d_left = [-base_d_right[0], base_d_right[1], base_d_right[2]]
    parts = {
        "input_A_and_C": 2 * (h / sqrt2) * area([foot_a, base_a_ab, base_a_ad], a),
        "input_B": (h / sqrt2) * area([foot_b_right, foot_b_left, base_b_left, base_b_right], b),
        "input_D": (k * h / sqrt_d) * area([foot_d_right, foot_d_left, base_d_left, base_d_right], d),
        "AB_and_CB_cylinders": 2 * incoming_flux(right, a, b, u_ab, r),
        "AD_and_CD_cylinders": 2 * incoming_flux(right, a, d, u_ad, r),
        "BD_bridge_cylinder": bridge_flux(left, b, d, 2 * t, r),
        "ABD_and_CBD_spheres": 2 * sphere_flux(right, [a, b, d], r),
    }
    result = {
        "method": "Divergence theorem with closed-form planar, cylindrical and spherical flux integrals; outward 70-decimal interval arithmetic, no quadrature or kernel calibration.",
        "parameters": {"height": "10", "radius": "1/2", "lower_y_slope": "1/2"},
        "right_sphere_centre": [interval_record(x) for x in right],
        "left_sphere_centre": [interval_record(x) for x in left],
        "flux_mm3": {name: interval_record(value) for name, value in parts.items()},
        "volume_mm3": interval_record(sum(parts.values()) / 3),
        "base_flux": "zero, because the retained base lies in z=0",
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
