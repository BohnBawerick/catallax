//! Allocation-free outward-rounded double intervals.

#![allow(clippy::float_cmp)] // Exact endpoint identities, not approximate equality.

use crate::{Integer, Rational, StrictSign};

/// Closed double interval. Arithmetic never stores NaN.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Interval {
    lo: f64,
    hi: f64,
}
impl Interval {
    /// Constructs an interval; rejects NaN and reversed ends.
    pub fn new(lo: f64, hi: f64) -> Option<Self> {
        if lo.is_nan() || hi.is_nan() || lo > hi {
            None
        } else {
            Some(Self {
                lo: if lo == 0.0 { 0.0 } else { lo },
                hi: if hi == 0.0 { 0.0 } else { hi },
            })
        }
    }
    /// Single point, or `None` for NaN.
    pub fn point(v: f64) -> Option<Self> {
        Self::new(v, v)
    }
    /// Whole real line.
    pub fn entire() -> Self {
        Self {
            lo: f64::NEG_INFINITY,
            hi: f64::INFINITY,
        }
    }
    /// Tight neighboring-double enclosure of a rational.
    pub fn from_rational(r: &Rational) -> Self {
        let v = r.to_f64_nearest();
        if v == f64::INFINITY {
            return Self {
                lo: f64::MAX,
                hi: v,
            };
        }
        if v == f64::NEG_INFINITY {
            return Self {
                lo: v,
                hi: -f64::MAX,
            };
        }
        let q = Rational::from_f64_exact(v).unwrap();
        Self::new(
            if q > *r { v.next_down() } else { v },
            if q < *r { v.next_up() } else { v },
        )
        .unwrap()
    }
    /// Tight neighboring-double enclosure of an integer.
    pub fn from_integer(n: &Integer) -> Self {
        Self::from_rational(&Rational::from_integer(n.clone()))
    }
    /// Lower end.
    pub fn lo(&self) -> f64 {
        self.lo
    }
    /// Upper end.
    pub fn hi(&self) -> f64 {
        self.hi
    }
    /// Membership, false for NaN.
    pub fn contains(&self, v: f64) -> bool {
        self.lo <= v && v <= self.hi
    }
    /// Whether zero lies inside.
    pub fn contains_zero(&self) -> bool {
        self.contains(0.0)
    }
    /// Whether both ends are finite.
    pub fn is_finite(&self) -> bool {
        self.lo.is_finite() && self.hi.is_finite()
    }
    /// Width rounded upwards.
    pub fn width(&self) -> f64 {
        if self.lo == self.hi {
            0.0
        } else {
            up(self.hi - self.lo)
        }
    }
    /// Finite representative near the midpoint.
    pub fn mid(&self) -> f64 {
        if self.lo == self.hi {
            return self.lo;
        }
        match (self.lo.is_infinite(), self.hi.is_infinite()) {
            (true, true) => 0.0,
            (true, false) => -f64::MAX,
            (false, true) => f64::MAX,
            _ => self.lo / 2.0 + self.hi / 2.0,
        }
    }
    /// Strict sign only when zero is excluded.
    pub fn sign(&self) -> Option<StrictSign> {
        if self.lo > 0.0 {
            Some(StrictSign::Positive)
        } else if self.hi < 0.0 {
            Some(StrictSign::Negative)
        } else {
            None
        }
    }
    /// Enclosure of the square.
    pub fn square(&self) -> Self {
        let a = self.lo.abs();
        let b = self.hi.abs();
        let lo = if self.contains_zero() {
            0.0
        } else {
            down(a.min(b) * a.min(b)).max(0.0)
        };
        Self::new(lo, up(a.max(b) * a.max(b))).unwrap()
    }
    /// Square root of the non-negative part; `None` if entirely negative.
    pub fn sqrt(&self) -> Option<Self> {
        if self.hi < 0.0 {
            return None;
        }
        let lo = self.lo.max(0.0).sqrt();
        let hi = self.hi.sqrt();
        let lower = if exact_square(lo, self.lo.max(0.0)) {
            lo
        } else {
            down(lo).max(0.0)
        };
        let upper = if exact_square(hi, self.hi) {
            hi
        } else {
            up(hi)
        };
        Self::new(lower, upper)
    }
    /// Smallest interval containing both.
    pub fn hull(&self, other: &Self) -> Self {
        Self::new(self.lo.min(other.lo), self.hi.max(other.hi)).unwrap()
    }
    fn add(&self, rhs: &Self) -> Self {
        if self.lo == 0.0 && self.hi == 0.0 {
            return *rhs;
        }
        if rhs.lo == 0.0 && rhs.hi == 0.0 {
            return *self;
        }
        Self::new(down(self.lo + rhs.lo), up(self.hi + rhs.hi)).unwrap_or(Self::entire())
    }
    fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }
    fn mul(&self, rhs: &Self) -> Self {
        if (self.lo == 0.0 && self.hi == 0.0) || (rhs.lo == 0.0 && rhs.hi == 0.0) {
            return Self::point(0.0).unwrap();
        }
        let values = [
            self.lo * rhs.lo,
            self.lo * rhs.hi,
            self.hi * rhs.lo,
            self.hi * rhs.hi,
        ];
        if values.iter().any(|x| x.is_nan()) {
            return Self::entire();
        }
        let lo = values.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Self::new(down(lo), up(hi)).unwrap()
    }
    fn div(&self, rhs: &Self) -> Self {
        if rhs.contains_zero() {
            return Self::entire();
        }
        self.mul(&Self::new(down(1.0 / rhs.hi), up(1.0 / rhs.lo)).unwrap_or(Self::entire()))
    }
    fn neg(&self) -> Self {
        Self::new(-self.hi, -self.lo).unwrap()
    }
}
fn down(v: f64) -> f64 {
    if v == f64::NEG_INFINITY {
        v
    } else {
        v.next_down()
    }
}
fn up(v: f64) -> f64 {
    if v == f64::INFINITY { v } else { v.next_up() }
}
// Compare the exact integer mantissas, not a rounded product.
fn exact_square(v: f64, x: f64) -> bool {
    if !v.is_finite() {
        return x == f64::INFINITY;
    }
    if v == 0.0 {
        return x == 0.0;
    }
    fn parts(v: f64) -> (u128, i32) {
        let bits = v.to_bits();
        let e = ((bits >> 52) & 2047) as i32;
        let m = u128::from(bits & ((1u64 << 52) - 1)) | if e == 0 { 0 } else { 1u128 << 52 };
        (m, if e == 0 { -1074 } else { e - 1075 })
    }
    fn normal(m: u128, e: i32) -> (u128, i32) {
        let z = m.trailing_zeros();
        (m >> z, e + z as i32)
    }
    let (m, e) = parts(v);
    let (n, f) = parts(x);
    normal(m * m, 2 * e) == normal(n, f)
}
macro_rules! interval_ops { ($($tr:ident::$m:ident),+) => { $(
    impl std::ops::$tr for Interval { type Output = Self; fn $m(self, rhs: Self) -> Self { Interval::$m(&self, &rhs) } }
)+ }; }
interval_ops!(Add::add, Sub::sub, Mul::mul, Div::div);
impl std::ops::Neg for Interval {
    type Output = Self;
    fn neg(self) -> Self {
        Interval::neg(&self)
    }
}
