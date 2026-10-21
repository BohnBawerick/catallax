//! Exact rationals in normal form.

use crate::{Integer, ParseErrorKind, ParseNumberError, Sign};
use num_bigint::BigInt;
use num_integer::Integer as _;
use num_traits::{One, Signed, ToPrimitive, Zero};
use std::{cmp::Ordering, fmt, str::FromStr};

/// Arbitrary-precision exact rational with positive coprime denominator.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Rational {
    pub(crate) n: Integer,
    pub(crate) d: Integer,
}
impl Default for Rational {
    fn default() -> Self {
        Self::zero()
    }
}
impl Ord for Rational {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.d == other.d {
            return self.n.cmp(&other.n);
        }
        let sign = self.sign().cmp(&other.sign());
        if sign != Ordering::Equal {
            return sign;
        }
        (&self.n.0 * &other.d.0).cmp(&(&other.n.0 * &self.d.0))
    }
}
impl PartialOrd for Rational {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Rational {
    /// Zero.
    pub fn zero() -> Self {
        Self::from_i64(0)
    }
    /// One.
    pub fn one() -> Self {
        Self::from_i64(1)
    }
    /// The rational equal to `n`.
    pub fn from_integer(n: Integer) -> Self {
        Self {
            n,
            d: Integer::one(),
        }
    }
    /// The rational equal to `v`.
    pub fn from_i64(v: i64) -> Self {
        Self::from_integer(Integer::from_i64(v))
    }
    /// Normalizes a fraction; `None` when the denominator is zero.
    pub fn new(numerator: Integer, denominator: Integer) -> Option<Self> {
        if denominator.is_zero() {
            return None;
        }
        // Binary enclosures repeatedly normalize dyadics. Their only possible
        // common factor is a power of two; do not run a general Euclidean GCD.
        let twos = denominator.0.trailing_zeros().unwrap();
        if denominator.0.bits() == twos + 1 {
            let cancel = numerator.0.trailing_zeros().unwrap_or(twos).min(twos);
            let n = numerator.0 >> cancel;
            return Some(Self {
                n: Integer(if denominator.0.is_negative() { -n } else { n }),
                d: Integer(denominator.0.abs() >> cancel),
            });
        }
        let mut g = numerator.0.gcd(&denominator.0);
        if denominator.0.is_negative() {
            g = -g;
        }
        Some(Self {
            n: Integer(numerator.0 / &g),
            d: Integer(denominator.0 / g),
        })
    }
    /// `numerator / denominator`; `None` when the denominator is zero.
    pub fn from_ratio(numerator: i64, denominator: i64) -> Option<Self> {
        Self::new(numerator.into(), denominator.into())
    }
    /// Exact decimal or fraction parse using the piece's grammar.
    pub fn parse(s: &str) -> Result<Self, ParseNumberError> {
        let err = |kind| ParseNumberError { kind };
        if s.is_empty() {
            return Err(err(ParseErrorKind::Empty));
        }
        let negative = s.starts_with('-');
        let body = s.strip_prefix(['+', '-']).unwrap_or(s);
        let digits = |s: &str| !s.is_empty() && s.bytes().all(|c| c.is_ascii_digit());
        if let Some((n, d)) = body.split_once('/') {
            if !digits(n) || !digits(d) {
                return Err(err(ParseErrorKind::InvalidSyntax));
            }
            let den = Integer::parse(d)?;
            if den.is_zero() {
                return Err(err(ParseErrorKind::ZeroDenominator));
            }
            let n = Integer::parse(n)?;
            return Ok(Self::new(if negative { -n } else { n }, den).expect("nonzero denominator"));
        }
        let mut parts = body.split(['e', 'E']);
        let mantissa = parts.next().unwrap();
        let exponent = parts.next();
        if parts.next().is_some() {
            return Err(err(ParseErrorKind::InvalidSyntax));
        }
        let mut point = mantissa.split('.');
        let whole = point.next().unwrap();
        let fraction = point.next();
        if point.next().is_some() || !digits(whole) || fraction.is_some_and(|f| !digits(f)) {
            return Err(err(ParseErrorKind::InvalidSyntax));
        }
        let exp = if let Some(e) = exponent {
            let ds = e.strip_prefix(['+', '-']).unwrap_or(e);
            if !digits(ds) {
                return Err(err(ParseErrorKind::InvalidSyntax));
            }
            let mut v = 0i64;
            for c in ds.bytes() {
                v = (v * 10 + i64::from(c - b'0')).min(100001);
            }
            if v > 100000 {
                return Err(err(ParseErrorKind::ExponentTooLarge));
            }
            if e.starts_with('-') { -v } else { v }
        } else {
            0
        };
        let fraction = fraction.unwrap_or("");
        let n = Integer::parse(&format!("{whole}{fraction}"))?;
        let r = Self::from_integer(if negative { -n } else { n });
        let scale = exp - i64::try_from(fraction.len()).expect("decimal length fits i64");
        let ten = Self::from_integer(
            Integer::from_i64(10).pow(
                scale
                    .unsigned_abs()
                    .try_into()
                    .expect("decimal length fits u32"),
            ),
        );
        Ok(if scale >= 0 { &r * &ten } else { &r / &ten })
    }
    /// Numerator in lowest terms.
    pub fn numerator(&self) -> &Integer {
        &self.n
    }
    /// Positive denominator in lowest terms.
    pub fn denominator(&self) -> &Integer {
        &self.d
    }
    /// Exact sign.
    pub fn sign(&self) -> Sign {
        self.n.sign()
    }
    /// Whether the value is zero.
    pub fn is_zero(&self) -> bool {
        self.n.is_zero()
    }
    /// Whether the denominator is one.
    pub fn is_integer(&self) -> bool {
        self.d.0.is_one()
    }
    /// Absolute value.
    pub fn abs(&self) -> Self {
        Self {
            n: self.n.abs(),
            d: self.d.clone(),
        }
    }
    /// Reciprocal; `None` for zero.
    pub fn recip(&self) -> Option<Self> {
        if self.n.is_zero() {
            None
        } else {
            Some(Self {
                n: if self.n.0.is_negative() {
                    -&self.d
                } else {
                    self.d.clone()
                },
                d: self.n.abs(),
            })
        }
    }
    /// Checked quotient.
    pub fn checked_div(&self, other: &Self) -> Option<Self> {
        Some(self * &other.recip()?)
    }
    /// Integer power; `None` for zero to a negative power.
    pub fn pow(&self, exp: i32) -> Option<Self> {
        let r = if exp < 0 { self.recip()? } else { self.clone() };
        Some(Self {
            n: r.n.pow(exp.unsigned_abs()),
            d: r.d.pow(exp.unsigned_abs()),
        })
    }
    /// Largest integer not above the value.
    pub fn floor(&self) -> Integer {
        Integer(self.n.0.div_floor(&self.d.0))
    }
    /// Smallest integer not below the value.
    pub fn ceil(&self) -> Integer {
        Integer(self.n.0.div_ceil(&self.d.0))
    }
    /// `self * 2^k`.
    pub fn mul_pow2(&self, k: i64) -> Self {
        if self.is_zero() || k == 0 {
            return self.clone();
        }
        if k > 0 {
            let cancel = self.d.0.trailing_zeros().unwrap().min(k as u64);
            Self {
                n: Integer(&self.n.0 << (k as u64 - cancel)),
                d: Integer(&self.d.0 >> cancel),
            }
        } else {
            let cancel = self.n.0.trailing_zeros().unwrap().min(k.unsigned_abs());
            Self {
                n: Integer(&self.n.0 >> cancel),
                d: Integer(&self.d.0 << (k.unsigned_abs() - cancel)),
            }
        }
    }
    /// Exact finite double value; `None` for NaN or infinity.
    pub fn from_f64_exact(v: f64) -> Option<Self> {
        if !v.is_finite() {
            return None;
        }
        let bits = v.to_bits();
        let e = ((bits >> 52) & 2047) as i64;
        let m = (bits & ((1u64 << 52) - 1)) | if e == 0 { 0 } else { 1u64 << 52 };
        let r = Self::from_integer(Integer::from_u64(m)).mul_pow2(if e == 0 {
            -1074
        } else {
            e - 1075
        });
        Some(if bits >> 63 == 0 { r } else { -r })
    }
    /// Nearest double, ties to even, with overflow to infinity.
    pub fn to_f64_nearest(&self) -> f64 {
        if self.is_zero() {
            return 0.0;
        }
        let r = self.abs();
        let e = r.binary_exponent();
        let sign = if self.sign() == Sign::Negative {
            1u64 << 63
        } else {
            0
        };
        if e > 1023 {
            return f64::from_bits(sign | (2047u64 << 52));
        }
        if e < -1075 {
            return f64::from_bits(sign);
        }
        let shift = if e < -1022 { 1074 } else { 52 - e };
        let scaled = r.mul_pow2(shift);
        let (mut q, rem) = scaled.n.0.div_rem(&scaled.d.0);
        let twice = rem << 1u32;
        if twice > scaled.d.0 || (twice == scaled.d.0 && q.is_odd()) {
            q += 1;
        }
        let m = q.to_u64().expect("53-bit rounded mantissa");
        if e < -1022 {
            return f64::from_bits(sign | m);
        }
        let (m, e) = if m == 1u64 << 53 {
            (m >> 1, e + 1)
        } else {
            (m, e)
        };
        if e > 1023 {
            f64::from_bits(sign | (2047u64 << 52))
        } else {
            f64::from_bits(sign | (((e + 1023) as u64) << 52) | (m - (1u64 << 52)))
        }
    }
    /// Exact terminating decimal, or `None` for other denominators.
    pub fn to_decimal_exact(&self) -> Option<String> {
        let mut d = self.d.0.clone();
        let mut twos = 0u32;
        let mut fives = 0u32;
        while d.is_even() {
            d >>= 1u32;
            twos += 1;
        }
        while (&d % 5u32).is_zero() {
            d /= 5u32;
            fives += 1;
        }
        if !d.is_one() {
            return None;
        }
        let places = twos.max(fives);
        let n =
            &self.n.0 * BigInt::from(2).pow(places - twos) * BigInt::from(5).pow(places - fives);
        Some(format_plain(&n, -(i64::from(places))))
    }
    /// Tight directed decimal bounds with the requested significant digits.
    pub fn to_decimal_bounds(&self, significant_digits: u32) -> (String, String) {
        if self.is_zero() {
            return ("0".into(), "0".into());
        }
        let r = self.abs();
        let mut e = r.n.0.to_str_radix(10).len() as i64 - r.d.0.to_str_radix(10).len() as i64;
        let ten = |k: i64| {
            Self::from_integer(Integer::from_i64(10).pow(k.unsigned_abs().try_into().unwrap()))
        };
        let power = if e >= 0 {
            ten(e)
        } else {
            ten(e).recip().unwrap()
        };
        if r < power {
            e -= 1;
        }
        let step = e - i64::from(significant_digits.max(1)) + 1;
        let scale = ten(step);
        let scaled = if step >= 0 {
            self / &scale
        } else {
            self * &scale
        };
        (
            format_significant(&scaled.floor().0, step),
            format_significant(&scaled.ceil().0, step),
        )
    }
    pub(crate) fn binary_exponent(&self) -> i64 {
        assert!(!self.is_zero());
        let mut e = self.n.bit_length() as i64 - self.d.bit_length() as i64;
        let less = if e >= 0 {
            self.n.0.abs() < (&self.d.0 << e as u64)
        } else {
            (self.n.0.abs() << e.unsigned_abs()) < self.d.0
        };
        if less {
            e -= 1;
        }
        e
    }
    pub(crate) fn add(&self, rhs: &Self) -> Self {
        self.add_sub(rhs, false)
    }
    fn add_sub(&self, rhs: &Self, subtract: bool) -> Self {
        if rhs.is_zero() {
            return self.clone();
        }
        if self.is_zero() {
            return if subtract { -rhs } else { rhs.clone() };
        }
        let a = self.d.0.trailing_zeros().unwrap();
        let b = rhs.d.0.trailing_zeros().unwrap();
        if self.d.0.bits() == a + 1 && rhs.d.0.bits() == b + 1 {
            let n = match a.cmp(&b) {
                Ordering::Less => {
                    let left = &self.n.0 << (b - a);
                    if subtract {
                        left - &rhs.n.0
                    } else {
                        left + &rhs.n.0
                    }
                }
                Ordering::Equal => {
                    if subtract {
                        &self.n.0 - &rhs.n.0
                    } else {
                        &self.n.0 + &rhs.n.0
                    }
                }
                Ordering::Greater => {
                    let right = &rhs.n.0 << (a - b);
                    if subtract {
                        &self.n.0 - right
                    } else {
                        &self.n.0 + right
                    }
                }
            };
            let exponent = a.max(b);
            let cancel = n.trailing_zeros().unwrap_or(exponent).min(exponent);
            return Self {
                n: Integer(n >> cancel),
                d: Integer(BigInt::one() << (exponent - cancel)),
            };
        }
        let g = self.d.0.gcd(&rhs.d.0);
        let a = &self.d.0 / &g;
        let b = &rhs.d.0 / &g;
        let left = &self.n.0 * &b;
        let right = &rhs.n.0 * &a;
        // The reduced numerator is coprime to both denominator quotients.
        // Only a factor of the shared denominator GCD can remain.
        let n = if subtract { left - right } else { left + right };
        let cancel = n.gcd(&g);
        Self {
            n: Integer(n / &cancel),
            d: Integer(a * (&rhs.d.0 / cancel)),
        }
    }
    pub(crate) fn sub(&self, rhs: &Self) -> Self {
        self.add_sub(rhs, true)
    }
    pub(crate) fn mul(&self, rhs: &Self) -> Self {
        let a = self.d.0.trailing_zeros().unwrap();
        let b = rhs.d.0.trailing_zeros().unwrap();
        if self.d.0.bits() == a + 1 && rhs.d.0.bits() == b + 1 {
            let left = self.n.0.trailing_zeros().unwrap_or(b).min(b);
            let right = rhs.n.0.trailing_zeros().unwrap_or(a).min(a);
            // Normalized dyadics can cancel on at most one side.
            let n = if left > 0 {
                (&self.n.0 >> left) * &rhs.n.0
            } else if right > 0 {
                &self.n.0 * (&rhs.n.0 >> right)
            } else {
                &self.n.0 * &rhs.n.0
            };
            return Self {
                n: Integer(n),
                d: Integer(BigInt::one() << (a + b - left - right)),
            };
        }
        let (a, d) = cancel(&self.n.0, &rhs.d.0);
        let (b, c) = cancel(&rhs.n.0, &self.d.0);
        Self {
            n: Integer(a * b),
            d: Integer(c * d),
        }
    }
    pub(crate) fn div(&self, rhs: &Self) -> Self {
        self.checked_div(rhs).expect("rational division by zero")
    }
    pub(crate) fn neg(&self) -> Self {
        Self {
            n: -&self.n,
            d: self.d.clone(),
        }
    }
}
// Denominators are positive. Dyadic cancellation only needs trailing bits.
fn cancel(n: &BigInt, d: &BigInt) -> (BigInt, BigInt) {
    let twos = d.trailing_zeros().unwrap();
    if d.bits() == twos + 1 {
        let shift = n.trailing_zeros().unwrap_or(twos).min(twos);
        (n >> shift, d >> shift)
    } else {
        let g = n.gcd(d);
        (n / &g, d / g)
    }
}
fn format_plain(n: &BigInt, step: i64) -> String {
    if n.is_zero() {
        return "0".into();
    }
    let sign = if n.is_negative() { "-" } else { "" };
    let mut digits = n.abs().to_str_radix(10);
    if step >= 0 {
        digits.push_str(&"0".repeat(step as usize));
    } else {
        let pos = digits.len() as i64 + step;
        if pos > 0 {
            digits.insert(pos as usize, '.');
        } else {
            digits = format!("0.{}{digits}", "0".repeat((-pos) as usize));
        }
        while digits.ends_with('0') {
            digits.pop();
        }
        if digits.ends_with('.') {
            digits.pop();
        }
    }
    format!("{sign}{digits}")
}
fn format_significant(n: &BigInt, step: i64) -> String {
    if n.is_zero() {
        return "0".into();
    }
    let mut digits = n.abs().to_str_radix(10);
    let e = digits.len() as i64 - 1 + step;
    if (-6..=20).contains(&e) {
        return format_plain(n, step);
    }
    while digits.len() > 1 && digits.ends_with('0') {
        digits.pop();
    }
    if digits.len() > 1 {
        digits.insert(1, '.');
    }
    format!("{}{digits}e{e:+}", if n.is_negative() { "-" } else { "" })
}
binops!(Rational; Add::add, Sub::sub, Mul::mul, Div::div);
negop!(Rational);
impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_integer() {
            self.n.fmt(f)
        } else {
            write!(f, "{}/{}", self.n, self.d)
        }
    }
}
impl FromStr for Rational {
    type Err = ParseNumberError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}
impl From<i64> for Rational {
    fn from(v: i64) -> Self {
        Self::from_i64(v)
    }
}
impl From<i32> for Rational {
    fn from(v: i32) -> Self {
        Self::from_i64(i64::from(v))
    }
}
impl From<Integer> for Rational {
    fn from(n: Integer) -> Self {
        Self::from_integer(n)
    }
}
