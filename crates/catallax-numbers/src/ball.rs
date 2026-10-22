//! Binary midpoint-radius arithmetic and elementary functions.
//!
//! Series are evaluated with outward binary rounding at guard precision.
//! Every truncation includes an explicit Taylor or geometric tail bound.

use crate::{Decider, Integer, Interval, Rational, Refusal, Sign, StrictSign};

mod binary;
use binary::Binary;

/// Binary midpoint-radius enclosure with a working precision in bits.
#[derive(Clone, Debug)]
pub struct Ball {
    mid: Binary,
    rad: Binary,
    precision: u32,
}
impl Ball {
    /// Contains a rational, rounding its midpoint to the requested significant bits.
    pub fn from_rational(r: &Rational, precision: u32) -> Self {
        let p = precision.max(1);
        if r.is_integer() {
            return Self::from_integer(r.numerator(), p);
        }
        let numerator = Binary::from_integer(r.numerator().clone());
        let denominator = Binary::from_integer(r.denominator().clone());
        let mid = numerator.div_rounded(&denominator, p, false);
        // |r - mid| = |n - mid*d| / d. Keep that residual dyadic;
        // normalizing a new rational would only add redundant GCDs.
        let error = (&numerator - &(&mid * &denominator)).into_abs();
        let rad = error.div_rounded(&denominator, p, true);
        Self {
            mid,
            rad,
            precision: p,
        }
    }
    /// Contains an integer at the requested precision.
    pub fn from_integer(n: &Integer, precision: u32) -> Self {
        Self::centered(&Binary::from_integer(n.clone()), &Binary::zero(), precision)
    }
    /// Working precision.
    pub fn precision(&self) -> u32 {
        self.precision
    }
    /// Reduces precision with outward rounding, or only changes the working precision.
    pub fn with_precision(&self, precision: u32) -> Self {
        let p = precision.max(1);
        if p >= self.precision {
            let mut b = self.clone();
            b.precision = p;
            b
        } else {
            Self::centered(&self.mid, &self.rad, p)
        }
    }
    /// Exact binary midpoint.
    pub fn mid(&self) -> Rational {
        self.mid.to_rational()
    }
    /// Exact binary non-negative radius.
    pub fn rad(&self) -> Rational {
        self.rad.to_rational()
    }
    /// Exact lower bound.
    pub fn lo(&self) -> Rational {
        self.lower().to_rational()
    }
    /// Exact upper bound.
    pub fn hi(&self) -> Rational {
        self.upper().to_rational()
    }
    fn lower(&self) -> Binary {
        &self.mid - &self.rad
    }
    fn upper(&self) -> Binary {
        &self.mid + &self.rad
    }
    /// Rational membership.
    pub fn contains(&self, r: &Rational) -> bool {
        self.lo() <= *r && *r <= self.hi()
    }
    /// Whether zero lies inside.
    pub fn contains_zero(&self) -> bool {
        self.mid.cmp_abs(&self.rad) != std::cmp::Ordering::Greater
    }
    /// Whether this ball is exactly the singleton zero, without copying its fields.
    pub fn is_exact_zero(&self) -> bool {
        self.mid.is_zero() && self.rad.is_zero()
    }
    /// Strict sign when zero is excluded.
    pub fn sign(&self) -> Option<StrictSign> {
        if self.contains_zero() {
            None
        } else if self.mid.sign() == Sign::Negative {
            Some(StrictSign::Negative)
        } else {
            Some(StrictSign::Positive)
        }
    }
    /// Double enclosure of the entire ball.
    pub fn to_interval(&self) -> Interval {
        Interval::new(
            Interval::from_rational(&self.lo()).lo(),
            Interval::from_rational(&self.hi()).hi(),
        )
        .unwrap()
    }
    /// Smallest binary ball containing both; uses the larger working precision.
    pub fn hull(&self, other: &Self) -> Self {
        let lo = self.lower().min(other.lower());
        let hi = self.upper().max(other.upper());
        Self {
            mid: (&lo + &hi).mul_pow2(-1),
            rad: (&hi - &lo).mul_pow2(-1),
            precision: self.precision.max(other.precision),
        }
    }
    /// Enclosure of the square.
    pub fn square(&self) -> Self {
        let lo = self.lower();
        let hi = self.upper();
        let a = &lo * &lo;
        let b = &hi * &hi;
        Self::bounds_binary(
            if self.contains_zero() {
                Binary::zero()
            } else {
                a.clone().min(b.clone())
            },
            a.max(b),
            self.precision,
        )
    }
    /// Square root of the non-negative part; `None` only when entirely negative.
    pub fn sqrt(&self) -> Option<Self> {
        let hi = self.hi();
        if hi.sign() == Sign::Negative {
            return None;
        }
        let lo = self.lo().max(Rational::zero());
        let a = sqrt_bounds(&lo, self.precision);
        let b = sqrt_bounds(&hi, self.precision);
        Some(Self::bounds(a.0, b.1, self.precision))
    }
    /// Reciprocal, or `None` if zero is inside.
    pub fn recip(&self) -> Option<Self> {
        if self.contains_zero() {
            return None;
        }
        let denominator = &(&self.mid * &self.mid) - &(&self.rad * &self.rad);
        let mid = self.mid.div_rounded(&denominator, self.precision, false);
        let error = (&self.mid - &(&mid * &denominator)).into_abs();
        let rad = (&self.rad + &error).div_rounded(&denominator, self.precision, true);
        Some(Self {
            mid,
            rad,
            precision: self.precision,
        })
    }
    /// Quotient, or `None` if the divisor contains zero.
    pub fn checked_div(&self, other: &Self) -> Option<Self> {
        Some(self * &other.recip()?)
    }
    /// Non-negative integer power.
    pub fn pow_u32(&self, mut exp: u32) -> Self {
        let mut b = self.clone();
        let mut r = Self::from_rational(&Rational::one(), self.precision);
        while exp != 0 {
            if exp & 1 != 0 {
                r = &r * &b;
            }
            exp >>= 1;
            if exp != 0 {
                b = b.square();
            }
        }
        r
    }
    /// Pi, from Machin's identity `16 atan(1/5) - 4 atan(1/239)`.
    pub fn pi(precision: u32) -> Self {
        pi(precision, None).expect("unbounded pi enclosure")
    }
    /// Pi enclosure, abandoning the series when the caller's deadline passes.
    pub fn pi_checked(precision: u32, dec: &Decider) -> Result<Self, Refusal> {
        pi(precision, Some(dec))
    }
    /// Exponential enclosure.
    pub fn exp(&self) -> Self {
        let p = self.precision;
        let a = exp_point(&self.lo(), p);
        let b = exp_point(&self.hi(), p);
        Self::bounds(a.lo(), b.hi(), p)
    }
    /// Natural logarithm, requiring a strictly positive lower bound.
    pub fn ln(&self) -> Option<Self> {
        if self.lo().sign() != Sign::Positive {
            return None;
        }
        let a = ln_point(&self.lo(), self.precision);
        let b = ln_point(&self.hi(), self.precision);
        Some(Self::bounds(a.lo(), b.hi(), self.precision))
    }
    /// Sine enclosure with argument reduction by a proven pi ball.
    pub fn sin(&self) -> Self {
        self.trig(false)
    }
    /// Cosine enclosure with argument reduction by a proven pi ball.
    pub fn cos(&self) -> Self {
        self.trig(true)
    }
    /// Arctangent enclosure.
    pub fn atan(&self) -> Self {
        atan(self, None).expect("unbounded arctangent enclosure")
    }
    /// Arctangent enclosure, abandoning the series when the caller's deadline passes.
    pub fn atan_checked(&self, dec: &Decider) -> Result<Self, Refusal> {
        atan(self, Some(dec))
    }
    /// Arcsine, requiring the entire ball inside `[-1, 1]`.
    pub fn asin(&self) -> Option<Self> {
        if self.lo() < Rational::from_i64(-1) || self.hi() > Rational::one() {
            return None;
        }
        let a = asin_point(&self.lo(), self.precision);
        let b = asin_point(&self.hi(), self.precision);
        Some(Self::bounds(a.lo(), b.hi(), self.precision))
    }
    /// Arccosine, requiring the entire ball inside `[-1, 1]`.
    pub fn acos(&self) -> Option<Self> {
        let p = self.precision;
        Some(
            (&(Self::pi(p.saturating_add(16)) * scalar_ratio(1, 2, p.saturating_add(16)))
                - &self.with_precision(p.saturating_add(16)).asin()?)
                .with_precision(p),
        )
    }
    /// Principal angle, or `None` when both coordinates contain zero.
    pub fn atan2(y: &Self, x: &Self) -> Option<Self> {
        if x.contains_zero() && y.contains_zero() {
            return None;
        }
        let p = x.precision.max(y.precision);
        let g = p.saturating_add(16);
        let x = x.with_precision(g);
        let y = y.with_precision(g);
        let pi = Self::pi(g);
        let result = if !y.contains_zero() && x.contains_zero() {
            let half = &pi * &scalar_ratio(1, 2, g);
            let half = if y.sign() == Some(StrictSign::Negative) {
                -half
            } else {
                half
            };
            &half - &x.checked_div(&y).unwrap().atan()
        } else if x.sign() == Some(StrictSign::Positive) {
            y.checked_div(&x).unwrap().atan()
        } else if y.lo().sign() == Sign::Negative && y.hi().sign() != Sign::Negative {
            Self::bounds(-pi.hi(), pi.hi(), g)
        } else {
            let angle = y.checked_div(&x).unwrap().atan();
            if y.hi().sign() == Sign::Negative {
                &angle - &pi
            } else {
                &angle + &pi
            }
        };
        Some(result.with_precision(p))
    }
    pub(crate) fn bounds(lo: Rational, hi: Rational, precision: u32) -> Self {
        assert!(lo <= hi);
        let p = precision.max(1);
        let center = (&lo + &hi).mul_pow2(-1);
        let mid = rounded(&center, p, false);
        // If mid lies outside [lo, hi], the farther endpoint still wins.
        let rad = rounded(&(&mid - &lo).max(&hi - &mid), p, true);
        Self {
            mid: Binary::from_rational(mid),
            rad: Binary::from_rational(rad),
            precision: p,
        }
    }
    fn bounds_binary(lo: Binary, hi: Binary, precision: u32) -> Self {
        debug_assert!(lo <= hi);
        Self::centered(
            &(&lo + &hi).mul_pow2(-1),
            &(&hi - &lo).mul_pow2(-1),
            precision,
        )
    }
    fn centered(center: &Binary, radius: &Binary, precision: u32) -> Self {
        let p = precision.max(1);
        let mid = center.rounded(p, false);
        // max(|mid-lo|, |hi-mid|) = radius + |mid-center| exactly.
        // An unchanged midpoint needs no subtraction to certify zero error.
        let rad = if &mid == center {
            radius.rounded(p, true)
        } else {
            let error = (center - &mid).into_abs();
            (radius + &error).rounded(p, true)
        };
        Self {
            mid,
            rad,
            precision: p,
        }
    }
    pub(crate) fn magnitude(&self) -> Rational {
        (&self.mid.abs() + &self.rad).to_rational()
    }
    fn inflate(&self, error: &Rational) -> Self {
        Self::bounds(&self.lo() - error, &self.hi() + error, self.precision)
    }
    fn add(&self, rhs: &Self) -> Self {
        Self::centered(
            &(&self.mid + &rhs.mid),
            &(&self.rad + &rhs.rad),
            self.precision.max(rhs.precision),
        )
    }
    fn sub(&self, rhs: &Self) -> Self {
        Self::centered(
            &(&self.mid - &rhs.mid),
            &(&self.rad + &rhs.rad),
            self.precision.max(rhs.precision),
        )
    }
    fn mul(&self, rhs: &Self) -> Self {
        let precision = self.precision.max(rhs.precision);
        if self.rad.is_zero() && rhs.rad.is_zero() {
            return Self::centered(&(&self.mid * &rhs.mid), &Binary::zero(), precision);
        }
        if self.rad.is_zero() {
            return Self::centered(
                &(&self.mid * &rhs.mid),
                &self.mid.abs_product(&rhs.rad),
                precision,
            );
        }
        if rhs.rad.is_zero() {
            return Self::centered(
                &(&self.mid * &rhs.mid),
                &rhs.mid.abs_product(&self.rad),
                precision,
            );
        }
        let left_stable = self.mid.cmp_abs(&self.rad) != std::cmp::Ordering::Less;
        let right_stable = rhs.mid.cmp_abs(&rhs.rad) != std::cmp::Ordering::Less;
        match (left_stable, right_stable) {
            (true, true) => {
                // The exact extremal products have center mn +/- rs and
                // radius |m|s + |n|r. Form these before the same final rounding,
                // without expanding four much wider dyadic endpoints.
                let product = &self.mid * &rhs.mid;
                let cross = &self.rad * &rhs.rad;
                let center = if self.mid.sign() == rhs.mid.sign() {
                    &product + &cross
                } else {
                    &product - &cross
                };
                let radius = &self.mid.abs_product(&rhs.rad) + &rhs.mid.abs_product(&self.rad);
                Self::centered(&center, &radius, precision)
            }
            (true, false) | (false, true) => {
                let (stable, crossing) = if left_stable {
                    (self, rhs)
                } else {
                    (rhs, self)
                };
                let far = if stable.mid.sign() == Sign::Negative {
                    stable.lower()
                } else {
                    stable.upper()
                };
                Self::centered(
                    &(&far * &crossing.mid),
                    &far.abs_product(&crossing.rad),
                    precision,
                )
            }
            (false, false) => {
                let a = self.lower();
                let b = self.upper();
                let c = rhs.lower();
                let d = rhs.upper();
                Self::bounds_binary((&a * &d).min(&b * &c), (&a * &c).max(&b * &d), precision)
            }
        }
    }
    fn neg(&self) -> Self {
        Self {
            mid: -&self.mid,
            rad: self.rad.clone(),
            precision: self.precision,
        }
    }
    fn trig(&self, cosine: bool) -> Self {
        let p = self.precision;
        let extra = if self.mid.is_zero() {
            0
        } else {
            self.mid.binary_exponent().max(0).min(i64::from(u32::MAX)) as u32
        };
        let g = p.saturating_add(extra).saturating_add(32);
        let pi = Self::pi(g);
        let period = &pi * &scalar(2, g);
        let turns = (&self.mid() / &period.mid()).mul_pow2(0);
        let n = (&turns + &Rational::from_ratio(1, 2).unwrap()).floor();
        let reduced = &self.with_precision(g) - &(&period * &Self::from_integer(&n, g));
        if reduced.magnitude() > Rational::from_i64(4)
            || reduced.rad() > Rational::from_ratio(1, 4).unwrap()
        {
            return Self::bounds(Rational::from_i64(-1), Rational::one(), p);
        }
        let square = -reduced.square();
        let mut term = if cosine {
            scalar(1, g)
        } else {
            reduced.clone()
        };
        let mut sum = term.clone();
        let mut k = 1u64;
        let eps = epsilon(g);
        loop {
            let a = if cosine { 2 * k - 1 } else { 2 * k };
            let divisor =
                Self::from_rational(&Rational::from_integer(Integer::from_u64(a * (a + 1))), g);
            term = (&term * &square).checked_div(&divisor).unwrap();
            let error = term.magnitude();
            if a > 4 && error <= eps {
                return sum.inflate(&error).with_precision(p);
            }
            sum = &sum + &term;
            k += 1;
        }
    }
}
fn rounded(r: &Rational, p: u32, upwards: bool) -> Rational {
    if r.is_zero() {
        return Rational::zero();
    }
    let numerator = &r.numerator().0;
    let denominator = &r.denominator().0;
    let power = denominator.trailing_zeros().unwrap();
    if denominator.bits() == power + 1 {
        let discard = numerator.bits().saturating_sub(u64::from(p));
        if discard == 0 || numerator.trailing_zeros().unwrap() >= discard {
            return r.clone();
        }
        // Signed shifts round down. The second expression implements
        // floor(scaled + 1/2), including negative ties toward +infinity.
        let q = if upwards {
            (numerator >> discard) + 1u32
        } else {
            ((numerator >> (discard - 1)) + 1u32) >> 1u32
        };
        return Rational::from_integer(Integer(q)).mul_pow2(discard as i64 - power as i64);
    }
    let step = r.binary_exponent() - i64::from(p) + 1;
    let scaled = r.mul_pow2(-step);
    let q = if upwards {
        scaled.ceil()
    } else {
        (&scaled + &Rational::from_ratio(1, 2).unwrap()).floor()
    };
    Rational::from_integer(q).mul_pow2(step)
}
pub(crate) fn sqrt_bounds(r: &Rational, p: u32) -> (Rational, Rational) {
    if r.is_zero() {
        return (Rational::zero(), Rational::zero());
    }
    assert!(r.sign() == Sign::Positive);
    let step = r.binary_exponent().div_euclid(2) - i64::from(p.max(1));
    let t = r.mul_pow2(-2 * step);
    let n = t.floor();
    let root = Integer(n.0.sqrt());
    let lo = Rational::from_integer(root.clone()).mul_pow2(step);
    let hi = if &lo * &lo == *r {
        lo.clone()
    } else {
        Rational::from_integer(&root + &Integer::one()).mul_pow2(step)
    };
    (lo, hi)
}
fn scalar(v: i64, p: u32) -> Ball {
    Ball::from_rational(&Rational::from_i64(v), p)
}
fn scalar_ratio(n: i64, d: i64, p: u32) -> Ball {
    Ball::from_rational(&Rational::from_ratio(n, d).unwrap(), p)
}
fn epsilon(p: u32) -> Rational {
    Rational::one().mul_pow2(-i64::from(p))
}
fn pi(precision: u32, dec: Option<&Decider>) -> Result<Ball, Refusal> {
    let p = precision.max(1);
    let g = p.saturating_add(32);
    let a = atan_series(&Rational::from_ratio(1, 5).unwrap(), g, dec)?;
    let b = atan_series(&Rational::from_ratio(1, 239).unwrap(), g, dec)?;
    Ok((&(&a * &scalar(16, g)) - &(&b * &scalar(4, g))).with_precision(p))
}
fn atan(value: &Ball, dec: Option<&Decider>) -> Result<Ball, Refusal> {
    let a = atan_point(&value.lo(), value.precision, dec)?;
    let b = atan_point(&value.hi(), value.precision, dec)?;
    Ok(Ball::bounds(a.lo(), b.hi(), value.precision))
}
fn atan_series(x: &Rational, p: u32, dec: Option<&Decider>) -> Result<Ball, Refusal> {
    let b = Ball::from_rational(x, p);
    let square = -b.square();
    let mut power = b.clone();
    let mut sum = b;
    let mut k = 1i64;
    loop {
        if let Some(dec) = dec {
            dec.check_deadline()?;
        }
        power = &power * &square;
        let term = power.checked_div(&scalar(2 * k + 1, p)).unwrap();
        let error = term.magnitude();
        if error <= epsilon(p) {
            return Ok(sum.inflate(&error));
        }
        sum = &sum + &term;
        k += 1;
    }
}
fn atan_point(x: &Rational, p: u32, dec: Option<&Decider>) -> Result<Ball, Refusal> {
    if let Some(dec) = dec {
        dec.check_deadline()?;
    }
    if x.sign() == Sign::Negative {
        return Ok(-atan_point(&(-x), p, dec)?);
    }
    let g = p.saturating_add(24);
    if *x > Rational::one() {
        return Ok((&(&pi(g, dec)? * &scalar_ratio(1, 2, g))
            - &atan_point(&x.recip().unwrap(), g, dec)?)
            .with_precision(p));
    }
    // Half-angle identity makes the alternating series geometrically convergent.
    let b = Ball::from_rational(x, g);
    let divisor = &scalar(1, g) + &(&scalar(1, g) + &b.square()).sqrt().unwrap();
    let t = b.checked_div(&divisor).unwrap();
    let a = atan_series(&t.lo(), g, dec)?;
    let c = atan_series(&t.hi(), g, dec)?;
    Ok((&Ball::bounds(a.lo(), c.hi(), g) * &scalar(2, g)).with_precision(p))
}
fn exp_point(x: &Rational, p: u32) -> Ball {
    if x.sign() == Sign::Negative {
        return exp_point(&(-x), p.saturating_add(16))
            .recip()
            .unwrap()
            .with_precision(p);
    }
    let k = if x.is_zero() {
        0
    } else {
        (x.binary_exponent() + 2).max(0) as u32
    };
    let g = p.saturating_add(k).saturating_add(32);
    let b = Ball::from_rational(&x.mul_pow2(-i64::from(k)), g);
    let mut sum = scalar(1, g);
    let mut term = sum.clone();
    let mut n = 1i64;
    loop {
        term = (&term * &b).checked_div(&scalar(n, g)).unwrap();
        let error = term.magnitude().mul_pow2(1);
        if error <= epsilon(g) {
            sum = sum.inflate(&error);
            break;
        }
        sum = &sum + &term;
        n += 1;
    }
    for _ in 0..k {
        sum = sum.square();
    }
    sum.with_precision(p)
}
fn log_unit(x: &Rational, p: u32) -> Ball {
    let t = (x - &Rational::one()) / (x + &Rational::one());
    let b = Ball::from_rational(&t, p);
    let square = b.square();
    let mut power = b.clone();
    let mut sum = b;
    let mut n = 1i64;
    loop {
        power = &power * &square;
        let term = power.checked_div(&scalar(2 * n + 1, p)).unwrap();
        // |t| <= 1/3, so the remaining odd-power tail is <= 2*next term.
        let error = term.magnitude().mul_pow2(1);
        if error <= epsilon(p) {
            return &sum.inflate(&error) * &scalar(2, p);
        }
        sum = &sum + &term;
        n += 1;
    }
}
fn ln_point(x: &Rational, p: u32) -> Ball {
    let e = x.binary_exponent();
    let extra = 64 - e.unsigned_abs().leading_zeros();
    let g = p.saturating_add(extra).saturating_add(32);
    let unit = log_unit(&x.mul_pow2(-e), g);
    if e == 0 {
        return unit.with_precision(p);
    }
    (&unit + &(&log_unit(&Rational::from_i64(2), g) * &scalar(e, g))).with_precision(p)
}
fn asin_point(x: &Rational, p: u32) -> Ball {
    if x.sign() == Sign::Negative {
        return -asin_point(&(-x), p);
    }
    let g = p.saturating_add(24);
    if *x == Rational::one() {
        return (&Ball::pi(g) * &scalar_ratio(1, 2, g)).with_precision(p);
    }
    let b = Ball::from_rational(x, g);
    let denominator = Ball::from_rational(&(&Rational::one() - &(x * x)), g)
        .sqrt()
        .unwrap();
    Ball::atan2(&b, &denominator).unwrap().with_precision(p)
}
binops!(Ball; Add::add, Sub::sub, Mul::mul);
negop!(Ball);
