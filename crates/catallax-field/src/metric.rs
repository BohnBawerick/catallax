//! Analytic Euclidean distances. Boolean composition lives in the plan, not here.

use catallax_numbers::{Ball, Rational};

pub(crate) const PRECISION: u32 = 192;

pub(crate) trait Scalar: Clone {
    fn integer(n: i64) -> Self;
    fn add(&self, b: &Self) -> Self;
    fn sub(&self, b: &Self) -> Self;
    fn mul(&self, b: &Self) -> Self;
    fn div(&self, b: &Self) -> Self;
    fn neg(&self) -> Self;
    fn sqrt(&self) -> Self;
    fn min(&self, b: &Self) -> Self;
    fn max(&self, b: &Self) -> Self;
    fn square(&self) -> Self {
        self.mul(self)
    }
}

impl Scalar for f64 {
    fn integer(n: i64) -> Self {
        n as f64
    }
    fn add(&self, b: &Self) -> Self {
        self + b
    }
    fn sub(&self, b: &Self) -> Self {
        self - b
    }
    fn mul(&self, b: &Self) -> Self {
        self * b
    }
    fn div(&self, b: &Self) -> Self {
        self / b
    }
    fn neg(&self) -> Self {
        -self
    }
    fn sqrt(&self) -> Self {
        software_sqrt(*self)
    }
    fn min(&self, b: &Self) -> Self {
        if self.is_nan() || b.is_nan() {
            f64::NAN
        } else {
            f64::min(*self, *b)
        }
    }
    fn max(&self, b: &Self) -> Self {
        if self.is_nan() || b.is_nan() {
            f64::NAN
        } else {
            f64::max(*self, *b)
        }
    }
}

// Fixed Newton iteration, with a bit-defined seed. Approximate, deterministic and
// display-only. No libm call, fused operation, architecture-specific SIMD or JIT.
fn software_sqrt(x: f64) -> f64 {
    if x == 0.0 || !x.is_finite() {
        return x;
    }
    if x < f64::MIN_POSITIVE {
        return software_sqrt(x * 18014398509481984.0) / 134217728.0;
    }
    let mut y = f64::from_bits((x.to_bits() >> 1) + 0x1ff8000000000000);
    for _ in 0..8 {
        y = 0.5 * (y + x / y);
    }
    y
}

fn enclosure(lo: Rational, hi: Rational) -> Ball {
    Ball::from_rational(&lo, PRECISION).hull(&Ball::from_rational(&hi, PRECISION))
}

impl Scalar for Ball {
    fn integer(n: i64) -> Self {
        Ball::from_rational(&Rational::from_i64(n), PRECISION)
    }
    fn add(&self, b: &Self) -> Self {
        self + b
    }
    fn sub(&self, b: &Self) -> Self {
        self - b
    }
    fn mul(&self, b: &Self) -> Self {
        self * b
    }
    fn div(&self, b: &Self) -> Self {
        self.checked_div(b)
            .expect("metric denominator is a positive exact constant")
    }
    fn neg(&self) -> Self {
        -self
    }
    fn sqrt(&self) -> Self {
        Ball::sqrt(self).expect("sum of squares has a non-negative upper bound")
    }
    fn square(&self) -> Self {
        Ball::square(self)
    }
    fn min(&self, b: &Self) -> Self {
        if self.hi() <= b.lo() {
            return self.clone();
        }
        if b.hi() <= self.lo() {
            return b.clone();
        }
        enclosure(self.lo().min(b.lo()), self.hi().min(b.hi()))
    }
    fn max(&self, b: &Self) -> Self {
        if self.lo() >= b.hi() {
            return self.clone();
        }
        if b.lo() >= self.hi() {
            return b.clone();
        }
        enclosure(self.lo().max(b.lo()), self.hi().max(b.hi()))
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Kind {
    Block,
    Round,
    Sphere,
    Torus,
}

/// Parameters: box sizes; or base radius, top radius, height, exact radius
/// difference and exact squared side length; or sphere radius; or torus radii.
pub(crate) fn distance<T: Scalar>(kind: Kind, a: &[T; 5], p: &[T; 3]) -> T {
    let zero = T::integer(0);
    match kind {
        Kind::Block => {
            let q: [T; 3] = std::array::from_fn(|i| p[i].neg().max(&p[i].sub(&a[i])));
            let outside = q[0]
                .max(&zero)
                .square()
                .add(&q[1].max(&zero).square())
                .add(&q[2].max(&zero).square())
                .sqrt();
            outside.add(&q[0].max(&q[1]).max(&q[2]).min(&zero))
        }
        Kind::Sphere => p[0]
            .square()
            .add(&p[1].square())
            .add(&p[2].square())
            .sqrt()
            .sub(&a[0]),
        Kind::Torus => {
            let radial = p[0].square().add(&p[1].square()).sqrt().sub(&a[0]);
            radial.square().add(&p[2].square()).sqrt().sub(&a[1])
        }
        Kind::Round => {
            let rho = p[0].square().add(&p[1].square()).sqrt();
            let t = rho
                .sub(&a[0])
                .mul(&a[3])
                .add(&p[2].mul(&a[2]))
                .div(&a[4])
                .max(&zero)
                .min(&T::integer(1));
            let side = rho
                .sub(&a[0].add(&t.mul(&a[3])))
                .square()
                .add(&p[2].sub(&t.mul(&a[2])).square())
                .sqrt();
            let bottom = rho
                .sub(&a[0])
                .max(&zero)
                .square()
                .add(&p[2].square())
                .sqrt();
            let top = rho
                .sub(&a[1])
                .max(&zero)
                .square()
                .add(&p[2].sub(&a[2]).square())
                .sqrt();
            // Unsigned distance. The caller supplies the exact or tentative sign.
            side.min(&bottom).min(&top)
        }
    }
}
