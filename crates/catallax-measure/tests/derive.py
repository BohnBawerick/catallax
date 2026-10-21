#!/usr/bin/env python3
"""Closed-form checks using only Python integer/Fraction arithmetic, no kernel.

Run: python3 crates/catallax-measure/tests/derive.py
Polynomial coefficients below omit a common pi factor where noted.
"""
from fractions import Fraction as F
from math import isqrt


def integral(coefficients, lo, hi):
    return sum((F(c, i + 1) * (hi ** (i + 1) - lo ** (i + 1))
                for i, c in enumerate(coefficients)), F(0))


def frustum(base, top, height):
    slope = F(top - base, height)
    section = [base * base, 2 * base * slope, slope * slope]
    volume_over_pi = integral(section, F(0), F(height))
    moment_over_pi = integral([0] + section, F(0), F(height))
    return volume_over_pi, moment_over_pi / volume_over_pi


assert frustum(2, 1, 2) == (F(14, 3), F(11, 14))
assert frustum(3, 0, 4) == (F(12), F(1))
# Slant heights sqrt((R-r)^2+h^2); lateral area/pi = (R+r)*slant.
assert (2 - 1) ** 2 + 2 ** 2 == 5
assert (3 - 0) ** 2 + 4 ** 2 == 25
assert 3 * isqrt(25) + 3 ** 2 == 24


def cap(radius, floor):
    section = [radius ** 2, 0, -1]
    volume = integral(section, F(floor), F(radius))
    moment = integral([0] + section, F(floor), F(radius))
    # Curved sphere area/pi = 2R(R-floor), plus closing disc radius squared.
    area = 2 * radius * (radius - floor) + radius ** 2 - floor ** 2
    return volume, area, moment / volume


assert cap(2, 0) == (F(16, 3), 12, F(3, 4))
assert cap(5, 3) == (F(52, 3), 36, F(48, 13))
R = F(3)
octant_volume = integral([R * R / 4, 0, F(-1, 4)], F(0), R)
octant_moment = integral([0, R * R / 4, 0, F(-1, 4)], F(0), R)
assert octant_volume == F(9, 2)
assert octant_moment / octant_volume == F(9, 8)
assert 4 * R * R / 8 + 3 * R * R / 4 == F(45, 4)

# Half-disc first moment: t=r^2-y^2 gives integral_0^r 2y*sqrt(r^2-y^2)dy=2r^3/3.
r, h = F(2), F(3)
half_volume_over_pi = r * r * h / 2
half_y_moment = F(2, 3) * r ** 3 * h
assert half_volume_over_pi == 6
assert half_y_moment / half_volume_over_pi == F(8, 3)  # divide by pi
assert r * h + r * r == 10  # curved side and two half discs, coefficient of pi
assert 2 * r * h == 12  # rectangle

# Steinmetz at y: square section of side 2*sqrt(1-y^2).
assert integral([4, 0, -4], F(-1), F(1)) == F(16, 3)
# One cylinder patch area: integral_-pi/2^pi/2 2*cos(u) du = 4.
assert 4 * 4 == 16

for size, volume, area in [((F(2), F(3), F(5)), F(30), F(62)),
                           ((F(1, 3), F(7), F(1, 1000)), F(7, 3000), F(7022, 1500))]:
    a, b, c = size
    assert a * b * c == volume
    assert 2 * (a * b + b * c + c * a) == area
assert 10 ** 3 - 2 ** 3 == 992
assert 6 * 10 ** 2 + 6 * 2 ** 2 == 624
assert 8 + 8 - 1 == 15
assert 24 + 24 - 6 == 42
assert (F(8) * 1 - F(3, 2)) / 7 == F(13, 14)
assert F(8) * 1 + F(8) * 2 - F(3, 2) == F(15) * F(3, 2)
# Hole area coefficient = lateral wall - two removed discs.
assert 2 * 1 * 4 - 2 * 1 ** 2 == 6
assert 3 ** 2 * 10 == 90
assert 2 * 3 * 10 - 2 * 3 ** 2 == 42
assert 2 * 3 * 10 + 2 * 3 ** 2 == 78
assert 5 ** 2 * 10 == 250
assert 2 * 5 * 10 - 2 * 5 ** 2 == 50
# Sphere/cylinder intersection, cylinder side + two spherical caps.
assert 3 ** 2 * 8 + 2 * cap(5, 4)[0] == F(244, 3)
assert 2 * 3 * 8 + 2 * 2 * 5 * 1 == 68
assert 2 * 2 ** 2 * 10 == 80
assert F(16, 3) * 2 ** 3 == F(128, 3)
# Pappus: torus volume = (pi*r^2)*(2*pi*R), area=(2*pi*r)*(2*pi*R).
assert 1 ** 2 * 2 * 3 == 6
assert 2 * 1 * 2 * 3 == 12


def atan_bounds(inverse, terms):
    partial = sum((F((-1) ** k, (2 * k + 1) * inverse ** (2 * k + 1))
                   for k in range(terms)), F(0))
    next_term = F((-1) ** terms, (2 * terms + 1) * inverse ** (2 * terms + 1))
    return min(partial, partial + next_term), max(partial, partial + next_term)


# Alternating-series remainders give a rigorous independent pi interval.
a, b = atan_bounds(5, 180), atan_bounds(239, 60)
pi_lo, pi_hi = 16 * a[0] - 4 * b[1], 16 * a[1] - 4 * b[0]
assert F(314159265358979323846264338327950288419716939937510, 10 ** 50) < pi_lo
assert pi_hi < F(314159265358979323846264338327950288419716939937511, 10 ** 50)
scale = 10 ** 100
sqrt5_lo = F(isqrt(5 * scale * scale), scale)
sqrt5_hi = sqrt5_lo + F(1, scale)
area_lo = (5 + 3 * sqrt5_lo) * pi_lo
area_hi = (5 + 3 * sqrt5_hi) * pi_hi
assert F(3678240746107114537178, 10 ** 20) < area_lo
assert area_hi < F(3678240746107114537179, 10 ** 20)
print("exact derivations passed: blocks, frusta, apex, cap, hemisphere, octant, half cylinders, Steinmetz, holes, void, torus, all nine boolean rows")
print("independent rational bounds passed: Machin pi, sqrt(5), frustum area")
