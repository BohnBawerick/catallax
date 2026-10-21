#!/usr/bin/env python3
"""Independent fixture derivations. Uses no Rust code or intersection implementation."""
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from math import isqrt


@dataclass(frozen=True)
class Radical:
    """a + b sqrt(d), with rational coefficients; enough for the listed radicals."""
    a: Q
    b: Q
    d: int

    def __add__(self, other):
        if not isinstance(other, Radical):
            other = Radical(Q(other), Q(0), self.d)
        assert self.d == other.d
        return Radical(self.a + other.a, self.b + other.b, self.d)

    __radd__ = __add__

    def __neg__(self):
        return Radical(-self.a, -self.b, self.d)

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if not isinstance(other, Radical):
            other = Radical(Q(other), Q(0), self.d)
        assert self.d == other.d
        return Radical(self.a * other.a + self.b * other.b * self.d,
                       self.a * other.b + self.b * other.a, self.d)

    __rmul__ = __mul__

    def square(self):
        return self * self

    def rational(self):
        assert self.b == 0
        return self.a


def root(d):
    return Radical(Q(0), Q(1), d)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return tuple(a[(i + 1) % 3] * b[(i + 2) % 3]
                 - a[(i + 2) % 3] * b[(i + 1) % 3] for i in range(3))


def exact_checks():
    # PC2: 3y+4z=0, x²+y²=25. The major vertex is (0,-5,15/4).
    normal = (Q(0), Q(3), Q(4))
    major = (Q(0), Q(-12, 25), Q(9, 25))
    vertex = (Q(0), Q(-5), Q(15, 4))
    assert dot(normal, vertex) == dot(normal, major) == 0
    assert dot(vertex, vertex) == Q(625, 16)
    assert vertex[0] ** 2 + vertex[1] ** 2 == 25
    assert cross(normal, tuple(x / Q(625, 16) for x in vertex)) == (Q(4, 5), 0, 0)
    # PC3 irrational witnesses; CP3 rational witnesses and h².
    assert ((2 * root(6)).square() + 1).rational() == 25
    assert 4**2 + 3**2 == (4 - 8)**2 + 3**2 == 25
    assert Q(25 - 4**2, 8**2) == Q(9, 64)
    # SS3 radical plane x=3, circle radius²=16, reference +z.
    assert 3**2 + 4**2 == (3 - 6)**2 + 4**2 == 25
    assert cross((6, 0, 0), (0, 1, 0)) == (0, 0, 6)
    # CC4 node and witnesses: x²+y²=4, (y-1)²+z²=1.
    assert root(3).square().rational() + 1 == 4
    assert (1 - 1)**2 + 1**2 == 1
    assert 0**2 + 2**2 == 4 and (2 - 1)**2 + 0**2 == 1
    assert cross((0, 4, 0), (0, 2, 0)) == (0, 0, 0)
    # CC6 oblique angle. The normals are orthogonal and major_sq are 10,10/9.
    minus = (-15, 0, 5)
    plus = (15, 0, 45)
    assert dot(minus, plus) == 0
    assert Q(2) / (1 - Q(4, 5)) == 10
    assert Q(2) / (1 + Q(4, 5)) == Q(10, 9)
    # A major vertex with squared length 10 and direction plus has x²=1,z²=9.
    assert Q(plus[0]**2 * 10, dot(plus, plus)) == 1
    assert Q(plus[2]**2 * 10, dot(plus, plus)) == 9
    # CS1 sphere centred at z=1, radius 6; circle centres are 1 +/- sqrt(11).
    for sign in (-1, 1):
        z = 1 + sign * root(11)
        assert ((z - 1).square() + 25).rational() == 36
    # Viviani z=1 hits, and CS6 witness (1,0,sqrt(3)).
    y = Q(1, 2) * root(3)
    assert (y.square() + Q(1, 4)).rational() == 1
    assert (y.square() + Q(9, 4) + 1).rational() == 4
    assert (root(3).square() + 1).rational() == 4
    # Case 27 corrected inner near misses; exact, even at 10^-30.
    eps = Q(1, 10**30)
    assert 4 - (3 - (1 + eps))**2 == 4 * eps - eps**2 > 0
    assert 4 - (3 - (1 - eps))**2 == -4 * eps - eps**2 < 0
    # Near PC/PS tangent radius² is 10eps-eps² and is not a rational square.
    h2 = 10 * eps - eps**2
    assert isqrt(h2.numerator)**2 != h2.numerator or isqrt(h2.denominator)**2 != h2.denominator
    # Case 30 elimination: x²+y²=25, y²+2(x-4)²=8.
    x = 8 - root(15)
    y2 = 25 - x.square()
    assert (y2 + 2 * (x - 4).square()).rational() == 8
    assert (x.square() - 16 * x + 49).rational() == 0
    assert 25 - 12 * 5 + 35 == 0
    print("Exact Fraction/radical checks passed: cases 4,6,7,12,15,19,22-25,27,30,38.")


def simpson(f, n):
    h = Decimal(2) / n
    total = f(Decimal(-1)) + f(Decimal(1))
    for i in range(1, n):
        total += (4 if i % 2 else 2) * f(-1 + h * i)
    return total * h / 3


def length_checks():
    # Rational substitution v=tan(t/4) removes trigonometric functions for L24.
    # For L25 first set sin(theta/2)=sqrt(3/8) sin(u), then v=tan(u/2).
    # These are independent smooth integrands on [-1,1], with no endpoint poles.
    # Decimal Simpson is a numerical cross-check, not a certified length oracle.
    with localcontext() as ctx:
        ctx.prec = 50
        def f24(v):
            return 4 * (1 + 6 * v*v + v**4).sqrt() / (1 + v*v)**2
        def f25(v):
            s = 2*v / (1 + v*v)
            c = (1 - v*v) / (1 + v*v)
            return 4 / (1 + v*v) * (Decimal('1.5')*c*c / (1 - Decimal('.375')*s*s) + 3*s*s).sqrt()
        for name, f, expected in [
            ('L24', f24, Decimal('7.64039557805542403580952416434288658')),
            ('L25', f25, Decimal('9.53883361768658472990916036073734440')),
        ]:
            coarse, fine = simpson(f, 8192), simpson(f, 16384)
            assert abs(coarse - fine) < Decimal('1e-12')
            assert abs(fine - expected) < Decimal('1e-12')
            print(f'{name}: {fine}, coarse/fine difference {abs(coarse-fine)}')


if __name__ == '__main__':
    exact_checks()
    length_checks()
