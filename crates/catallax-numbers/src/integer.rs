//! Arbitrary-precision integers, exact signs and number parse errors.

use num_bigint::BigInt;
use num_integer::Integer as _;
use num_traits::{Signed, ToPrimitive, Zero};
use std::{fmt, str::FromStr};

/// Exact sign of an exact value.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Sign {
    /// Below zero.
    Negative,
    /// Exactly zero.
    Zero,
    /// Above zero.
    Positive,
}

/// What made a number string unacceptable.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ParseErrorKind {
    /// The string was empty.
    Empty,
    /// The string does not match the grammar.
    InvalidSyntax,
    /// A fraction had a zero denominator.
    ZeroDenominator,
    /// The decimal exponent's magnitude exceeds 100,000.
    ExponentTooLarge,
}

/// A rejected exact-number string.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct ParseNumberError {
    /// Why the string was rejected.
    pub kind: ParseErrorKind,
}
impl fmt::Display for ParseNumberError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self.kind {
            ParseErrorKind::Empty => "empty number",
            ParseErrorKind::InvalidSyntax => "invalid number syntax",
            ParseErrorKind::ZeroDenominator => "zero denominator",
            ParseErrorKind::ExponentTooLarge => "decimal exponent exceeds 100000",
        })
    }
}
impl std::error::Error for ParseNumberError {}

/// Arbitrary-precision signed integer.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub struct Integer(pub(crate) BigInt);
impl Integer {
    /// Zero.
    pub fn zero() -> Self {
        Self(BigInt::from(0))
    }
    /// One.
    pub fn one() -> Self {
        Self(BigInt::from(1))
    }
    /// The integer equal to `v`.
    pub fn from_i64(v: i64) -> Self {
        Self(BigInt::from(v))
    }
    /// The integer equal to `v`.
    pub fn from_u64(v: u64) -> Self {
        Self(BigInt::from(v))
    }
    /// Parses `sign? digits`, with no surrounding whitespace.
    pub fn parse(s: &str) -> Result<Self, ParseNumberError> {
        let digits = s.strip_prefix(['+', '-']).unwrap_or(s);
        if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
            return Err(ParseNumberError {
                kind: if s.is_empty() {
                    ParseErrorKind::Empty
                } else {
                    ParseErrorKind::InvalidSyntax
                },
            });
        }
        Ok(Self(
            BigInt::parse_bytes(s.as_bytes(), 10).expect("validated integer"),
        ))
    }
    /// The value as `i64` when it fits.
    pub fn to_i64(&self) -> Option<i64> {
        self.0.to_i64()
    }
    /// Exact sign.
    pub fn sign(&self) -> Sign {
        if self.0.is_zero() {
            Sign::Zero
        } else if self.0.is_negative() {
            Sign::Negative
        } else {
            Sign::Positive
        }
    }
    /// Whether the value is zero.
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    /// Absolute value.
    pub fn abs(&self) -> Self {
        Self(self.0.abs())
    }
    /// Number of bits of the absolute value; 0 for zero.
    pub fn bit_length(&self) -> u64 {
        self.0.bits()
    }
    /// `self` raised to `exp`; `0^0` is 1.
    pub fn pow(&self, exp: u32) -> Self {
        Self(self.0.pow(exp))
    }
    /// Non-negative greatest common divisor; `gcd(0, 0)` is 0.
    pub fn gcd(&self, other: &Self) -> Self {
        Self(self.0.gcd(&other.0))
    }
    /// Floor quotient and remainder; the remainder has the divisor's sign.
    /// Returns `None` when the divisor is zero.
    pub fn div_rem_floor(&self, d: &Self) -> Option<(Self, Self)> {
        if d.is_zero() {
            return None;
        }
        let (q, r) = self.0.div_mod_floor(&d.0);
        Some((Self(q), Self(r)))
    }
    /// `self * 2^k`.
    pub fn mul_pow2(&self, k: u64) -> Self {
        Self(&self.0 << k)
    }
    pub(crate) fn add(&self, rhs: &Self) -> Self {
        Self(&self.0 + &rhs.0)
    }
    pub(crate) fn sub(&self, rhs: &Self) -> Self {
        Self(&self.0 - &rhs.0)
    }
    pub(crate) fn mul(&self, rhs: &Self) -> Self {
        Self(&self.0 * &rhs.0)
    }
    pub(crate) fn neg(&self) -> Self {
        Self(-&self.0)
    }
}
binops!(Integer; Add::add, Sub::sub, Mul::mul);
negop!(Integer);
impl fmt::Display for Integer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}
impl FromStr for Integer {
    type Err = ParseNumberError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}
impl From<i64> for Integer {
    fn from(v: i64) -> Self {
        Self::from_i64(v)
    }
}
impl From<u64> for Integer {
    fn from(v: u64) -> Self {
        Self::from_u64(v)
    }
}
impl From<i32> for Integer {
    fn from(v: i32) -> Self {
        Self::from_i64(i64::from(v))
    }
}
