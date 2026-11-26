//! Exact rationals in normal form.

use crate::gcd::{gcd, word};
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
        if let (Some((a, b)), Some((c, d))) = (self.small(), other.small()) {
            let wide = i128::from;
            return (wide(a) * wide(d)).cmp(&(wide(c) * wide(b)));
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
            let negative = denominator.0.is_negative();
            let d = if negative {
                -denominator.0
            } else {
                denominator.0
            };
            let cancel = numerator.0.trailing_zeros().unwrap_or(twos).min(twos);
            let n = numerator.0 >> cancel;
            return Some(Self {
                n: Integer(if negative { -n } else { n }),
                d: Integer(d >> cancel),
            });
        }
        let mut g = gcd(&numerator.0, &denominator.0);
        if g.is_one() {
            return Some(if denominator.0.is_negative() {
                Self {
                    n: Integer(-numerator.0),
                    d: Integer(-denominator.0),
                }
            } else {
                Self {
                    n: numerator,
                    d: denominator,
                }
            });
        }
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
        let power = Self {
            n: r.n.pow(exp.unsigned_abs()),
            d: r.d.pow(exp.unsigned_abs()),
        };
        // Counted as the products a repeated multiplication takes, one fewer
        // than the exponent, each priced as the largest of them: the power by
        // the base.
        crate::operation::count_arithmetic(
            u64::from(exp.unsigned_abs().saturating_sub(1))
                .saturating_mul(power.limbs().saturating_mul(self.limbs())),
        );
        Some(power)
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
        if m == 0 {
            return Some(Self::zero());
        }
        // Normalize the 53-bit significand before allocating big integers.
        // Every finite binary64 denominator is a power of two.
        let shift = if e == 0 { -1074 } else { e - 1075 };
        let cancel = i64::from(m.trailing_zeros()).min((-shift).max(0));
        let mantissa = (m >> cancel) as i64;
        let numerator = BigInt::from(if bits >> 63 == 0 { mantissa } else { -mantissa });
        let shift = shift + cancel;
        Some(if shift >= 0 {
            Self::from_integer(Integer(numerator << shift as u64))
        } else {
            Self {
                n: Integer(numerator),
                d: Integer(BigInt::one() << shift.unsigned_abs()),
            }
        })
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
        let _whole = self.decimal_conversion();
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
        let _whole = self.decimal_conversion();
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
    // Counts the writing of this value in decimal as one operation. That
    // converts the numerator and the denominator to base ten and divides one
    // by the other, each quadratic in the value's own limbs with schoolbook
    // arithmetic (Knuth, volume 2, section 4.3.1 algorithm D and section 4.4),
    // so it is priced as the value times itself.
    fn decimal_conversion(&self) -> crate::operation::Uncounted {
        crate::operation::Uncounted::whole(self.limbs().saturating_mul(self.limbs()))
    }
    /// The limbs of the numerator and of the denominator, added: the size
    /// the decider prices this value's arithmetic by (`decider.rs`, "Logical
    /// arithmetic costs").
    pub fn limbs(&self) -> u64 {
        self.n
            .bit_length()
            .div_ceil(64)
            .max(1)
            .saturating_add(self.d.bit_length().div_ceil(64).max(1))
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
        crate::operation::count_arithmetic(self.limbs().saturating_add(rhs.limbs()));
        if rhs.is_zero() {
            return self.clone();
        }
        if self.is_zero() {
            return if subtract { -rhs } else { rhs.clone() };
        }
        if let (Some((a, b)), Some((c, d))) = (self.small(), rhs.small()) {
            return Self::sum_of_small(a, b, if subtract { -c } else { c }, d);
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
        if self.d == rhs.d {
            let numerator = if subtract {
                &self.n - &rhs.n
            } else {
                &self.n + &rhs.n
            };
            return Self::new(numerator, self.d.clone()).expect("positive denominator");
        }
        let g = gcd(&self.d.0, &rhs.d.0);
        if g.is_one() {
            // Coprime denominators: the cross sum is already in lowest terms.
            let left = &self.n.0 * &rhs.d.0;
            let right = &rhs.n.0 * &self.d.0;
            return Self {
                n: Integer(if subtract { left - right } else { left + right }),
                d: Integer(&self.d.0 * &rhs.d.0),
            };
        }
        let a = &self.d.0 / &g;
        let b = &rhs.d.0 / &g;
        let left = &self.n.0 * &b;
        let right = &rhs.n.0 * &a;
        // The reduced numerator is coprime to both denominator quotients.
        // Only a factor of the shared denominator GCD can remain.
        let n = if subtract { left - right } else { left + right };
        let cancel = gcd(&n, &g);
        if cancel.is_one() {
            return Self {
                n: Integer(n),
                d: Integer(a * &rhs.d.0),
            };
        }
        Self {
            n: Integer(n / &cancel),
            d: Integer(a * (&rhs.d.0 / cancel)),
        }
    }
    pub(crate) fn sub(&self, rhs: &Self) -> Self {
        self.add_sub(rhs, true)
    }
    // Parts below 2^62 in magnitude keep every product and sum of two
    // operands below 2^125, so machine arithmetic is exact.
    fn small(&self) -> Option<(i64, i64)> {
        const LIMIT: u64 = 1 << 62;
        let n = self.n.0.to_i64()?;
        let d = self.d.0.to_i64()?;
        (n.unsigned_abs() < LIMIT && d.unsigned_abs() < LIMIT).then_some((n, d))
    }
    // `a/b + c/d` in lowest terms, for two small operands in lowest terms.
    //
    // Knuth's method for fractions (The Art of Computer Programming, volume
    // 2, section 4.5.1): with g the divisor of the two denominators, only a
    // factor of g can be cancelled from the sum over the common denominator.
    // Every divisor is then taken of single words, and none of the
    // double-word cross products, which is what the reduction cost before.
    fn sum_of_small(a: i64, b: i64, c: i64, d: i64) -> Self {
        let wide = i128::from;
        if (b as u64).is_power_of_two() && (d as u64).is_power_of_two() {
            let (b, d) = (b.trailing_zeros(), d.trailing_zeros());
            let exponent = b.max(d);
            let n = (wide(a) << (exponent - b)) + (wide(c) << (exponent - d));
            let cancel = n.trailing_zeros().min(exponent);
            return Self {
                n: Integer(BigInt::from(n >> cancel)),
                d: Integer(BigInt::from(1u64 << (exponent - cancel))),
            };
        }
        let g = word(b.unsigned_abs(), d.unsigned_abs()) as i64;
        let (n, d) = if g == 1 {
            (wide(a) * wide(d) + wide(c) * wide(b), wide(b) * wide(d))
        } else {
            let (b, d_over_g) = (b / g, d / g);
            let n = wide(a) * wide(d_over_g) + wide(c) * wide(b);
            // Double-word division is a library call; divide in one word
            // where the numerator fits one.
            let (cancel, n) = match i64::try_from(n) {
                Ok(n) => {
                    let cancel = word((n % g).unsigned_abs(), g.unsigned_abs()) as i64;
                    (cancel, wide(n / cancel))
                }
                Err(_) => {
                    let rest = (n % wide(g)).unsigned_abs() as u64;
                    let cancel = word(rest, g.unsigned_abs()) as i64;
                    (cancel, if cancel == 1 { n } else { n / wide(cancel) })
                }
            };
            (n, wide(b) * wide(d / cancel))
        };
        if n == 0 {
            return Self::zero();
        }
        Self {
            n: Integer(BigInt::from(n)),
            d: Integer(BigInt::from(d)),
        }
    }
    // `a/b * c/d` in lowest terms, for two small non-zero operands in lowest
    // terms: each numerator is cancelled against the other denominator
    // first, and the products of what is left share no factor (the same
    // section of Knuth).
    fn product_of_small(a: i64, b: i64, c: i64, d: i64) -> Self {
        let wide = i128::from;
        let (a, d) = cancel_small(a, d);
        let (c, b) = cancel_small(c, b);
        Self {
            n: Integer(BigInt::from(wide(a) * wide(c))),
            d: Integer(BigInt::from(wide(b) * wide(d))),
        }
    }
    pub(crate) fn mul(&self, rhs: &Self) -> Self {
        crate::operation::count_arithmetic(self.limbs().saturating_mul(rhs.limbs()));
        if self.is_zero() || rhs.is_zero() {
            return Self::zero();
        }
        if self.n.0.is_one() && self.d.0.is_one() {
            return rhs.clone();
        }
        if rhs.n.0.is_one() && rhs.d.0.is_one() {
            return self.clone();
        }
        if let (Some((a, b)), Some((c, d))) = (self.small(), rhs.small()) {
            return Self::product_of_small(a, b, c, d);
        }
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
fn cancel_small(n: i64, d: i64) -> (i64, i64) {
    if (d as u64).is_power_of_two() {
        let shift = n.trailing_zeros().min(d.trailing_zeros());
        return (n >> shift, d >> shift);
    }
    let g = word(n.unsigned_abs(), d as u64) as i64;
    if g == 1 { (n, d) } else { (n / g, d / g) }
}

// Denominators are positive. Dyadic cancellation only needs trailing bits.
fn cancel(n: &BigInt, d: &BigInt) -> (BigInt, BigInt) {
    let twos = d.trailing_zeros().unwrap();
    if d.bits() == twos + 1 {
        let shift = n.trailing_zeros().unwrap_or(twos).min(twos);
        (n >> shift, d >> shift)
    } else {
        let g = gcd(n, d);
        if g.is_one() {
            (n.clone(), d.clone())
        } else {
            (n / &g, d / g)
        }
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
