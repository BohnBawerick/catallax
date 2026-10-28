//! Exact rational helpers for recognised data.

use catallax_numbers::{Integer, Rational, Sign};

use crate::Vec3;

pub(crate) fn sub(a: &Vec3, b: &Vec3) -> Vec3 {
    std::array::from_fn(|i| &a[i] - &b[i])
}

pub(crate) fn dot(a: &Vec3, b: &Vec3) -> Rational {
    (0..3).fold(Rational::zero(), |s, i| s + &a[i] * &b[i])
}

pub(crate) fn cross(a: &Vec3, b: &Vec3) -> Vec3 {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

pub(crate) fn is_zero(a: &Vec3) -> bool {
    a.iter().all(Rational::is_zero)
}

/// Floor of the square root of a non-negative integer.
pub(crate) fn isqrt(n: &Integer) -> Integer {
    if n.sign() != Sign::Positive {
        return Integer::zero();
    }
    let two = Integer::from_i64(2);
    let mut x = Integer::one().mul_pow2(n.bit_length().div_ceil(2));
    loop {
        let (q, _) = n.div_rem_floor(&x).expect("nonzero");
        let (y, _) = (&x + &q).div_rem_floor(&two).expect("two");
        if y >= x {
            return x;
        }
        x = y;
    }
}

/// The square root when it is rational.
pub(crate) fn rational_sqrt(q: &Rational) -> Option<Rational> {
    if q.sign() == Sign::Negative {
        return None;
    }
    let (n, d) = (isqrt(q.numerator()), isqrt(q.denominator()));
    (&n * &n == *q.numerator() && &d * &d == *q.denominator())
        .then(|| Rational::new(n, d).expect("nonzero denominator"))
}

/// The direction scaled to length one, when its length is rational.
pub(crate) fn unit(v: &Vec3) -> Option<Vec3> {
    let len = rational_sqrt(&dot(v, v))?;
    Some(std::array::from_fn(|i| &v[i] / &len))
}
