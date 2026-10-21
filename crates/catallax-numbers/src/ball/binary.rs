//! Exact signed dyadics. Nonzero mantissas are odd; zero has exponent zero.

use crate::{Integer, Rational, Sign};
use std::{
    cmp::Ordering,
    ops::{Add, Mul, Neg, Sub},
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Binary {
    mantissa: Integer,
    exponent: i64,
}
impl Binary {
    fn new(mut mantissa: Integer, mut exponent: i64) -> Self {
        if mantissa.is_zero() {
            exponent = 0;
        } else {
            let twos = mantissa.0.trailing_zeros().unwrap();
            if twos != 0 {
                mantissa.0 >>= twos;
                exponent += twos as i64;
            }
        }
        Self { mantissa, exponent }
    }

    pub(super) fn from_rational(value: Rational) -> Self {
        let power = value.d.bit_length() - 1;
        debug_assert_eq!(value.d.0.trailing_zeros(), Some(power));
        Self::new(value.n, -(power as i64))
    }

    pub(super) fn to_rational(&self) -> Rational {
        if self.exponent >= 0 {
            Rational::from_integer(self.mantissa.mul_pow2(self.exponent as u64))
        } else {
            Rational {
                n: self.mantissa.clone(),
                d: Integer::one().mul_pow2(self.exponent.unsigned_abs()),
            }
        }
    }

    pub(super) fn zero() -> Self {
        Self {
            mantissa: Integer::zero(),
            exponent: 0,
        }
    }

    pub(super) fn is_zero(&self) -> bool {
        self.mantissa.is_zero()
    }

    pub(super) fn sign(&self) -> Sign {
        self.mantissa.sign()
    }

    pub(super) fn abs(&self) -> Self {
        Self {
            mantissa: self.mantissa.abs(),
            exponent: self.exponent,
        }
    }

    pub(super) fn binary_exponent(&self) -> i64 {
        assert!(!self.is_zero());
        self.exponent + self.mantissa.bit_length() as i64 - 1
    }

    pub(super) fn mul_pow2(mut self, power: i64) -> Self {
        if !self.is_zero() {
            self.exponent += power;
        }
        self
    }

    pub(super) fn rounded(&self, precision: u32, upwards: bool) -> Self {
        let discard = self
            .mantissa
            .bit_length()
            .saturating_sub(u64::from(precision));
        if discard == 0 {
            return self.clone();
        }
        // The odd mantissa has a nonzero discarded part. Signed shifts floor;
        // nearest ties use the same positive-infinity convention as Ball.
        let mantissa = if upwards {
            (&self.mantissa.0 >> discard) + 1u32
        } else {
            ((&self.mantissa.0 >> (discard - 1)) + 1u32) >> 1u32
        };
        Self::new(Integer(mantissa), self.exponent + discard as i64)
    }

    pub(super) fn div_rounded(&self, denominator: &Self, precision: u32, upwards: bool) -> Self {
        debug_assert_eq!(denominator.sign(), Sign::Positive);
        debug_assert!(precision > 0);
        if self.is_zero() {
            return Self::zero();
        }
        if denominator.mantissa.bit_length() == 1 {
            return self
                .rounded(precision, upwards)
                .mul_pow2(-denominator.exponent);
        }
        let numerator = self.mantissa.0.magnitude();
        let divisor = denominator.mantissa.0.magnitude();
        let mut power = numerator.bits() as i64 - divisor.bits() as i64;
        let below = match power.cmp(&0) {
            Ordering::Equal => numerator < divisor,
            Ordering::Greater => numerator < &(divisor << power as u64),
            Ordering::Less => &(numerator << (-power) as u64) < divisor,
        };
        if below {
            power -= 1;
        }
        let shift = i64::from(precision) - power - 1;
        let divide = |numerator: &Integer, denominator: &Integer| {
            let (mut quotient, remainder) = numerator.div_rem_floor(denominator).unwrap();
            if !remainder.is_zero() && (upwards || remainder.mul_pow2(1) >= *denominator) {
                quotient.0 += 1u32;
            }
            quotient
        };
        let mantissa = match shift.cmp(&0) {
            Ordering::Greater => {
                divide(&self.mantissa.mul_pow2(shift as u64), &denominator.mantissa)
            }
            Ordering::Less => divide(
                &self.mantissa,
                &denominator.mantissa.mul_pow2((-shift) as u64),
            ),
            Ordering::Equal => divide(&self.mantissa, &denominator.mantissa),
        };
        Self::new(mantissa, self.exponent - denominator.exponent - shift)
    }

    pub(super) fn cmp_abs(&self, other: &Self) -> Ordering {
        if self.is_zero() || other.is_zero() {
            return self
                .mantissa
                .0
                .magnitude()
                .cmp(other.mantissa.0.magnitude());
        }
        let magnitude = self.binary_exponent().cmp(&other.binary_exponent());
        if magnitude != Ordering::Equal {
            return magnitude;
        }
        match self.exponent.cmp(&other.exponent) {
            Ordering::Equal => self
                .mantissa
                .0
                .magnitude()
                .cmp(other.mantissa.0.magnitude()),
            Ordering::Less => self
                .mantissa
                .0
                .magnitude()
                .cmp(&(other.mantissa.0.magnitude() << (other.exponent - self.exponent) as u64)),
            Ordering::Greater => (self.mantissa.0.magnitude()
                << (self.exponent - other.exponent) as u64)
                .cmp(other.mantissa.0.magnitude()),
        }
    }

    fn add_sub(&self, other: &Self, subtract: bool) -> Self {
        if other.is_zero() {
            return self.clone();
        }
        if self.is_zero() {
            return if subtract { -other } else { other.clone() };
        }
        let exponent = self.exponent.min(other.exponent);
        let mantissa = match self.exponent.cmp(&other.exponent) {
            Ordering::Equal => {
                if subtract {
                    &self.mantissa.0 - &other.mantissa.0
                } else {
                    &self.mantissa.0 + &other.mantissa.0
                }
            }
            Ordering::Greater => {
                let mut value = &self.mantissa.0 << (self.exponent - exponent) as u64;
                if subtract {
                    value -= &other.mantissa.0;
                } else {
                    value += &other.mantissa.0;
                }
                value
            }
            Ordering::Less => {
                let shifted = &other.mantissa.0 << (other.exponent - exponent) as u64;
                if subtract {
                    &self.mantissa.0 - shifted
                } else {
                    shifted + &self.mantissa.0
                }
            }
        };
        Self::new(Integer(mantissa), exponent)
    }
}
impl Ord for Binary {
    fn cmp(&self, other: &Self) -> Ordering {
        let sign = self.sign().cmp(&other.sign());
        if sign != Ordering::Equal {
            return sign;
        }
        let magnitude = self.cmp_abs(other);
        if self.sign() == Sign::Negative {
            magnitude.reverse()
        } else {
            magnitude
        }
    }
}
impl PartialOrd for Binary {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Add for &Binary {
    type Output = Binary;
    fn add(self, other: Self) -> Binary {
        self.add_sub(other, false)
    }
}
impl Sub for &Binary {
    type Output = Binary;
    fn sub(self, other: Self) -> Binary {
        self.add_sub(other, true)
    }
}
impl Mul for &Binary {
    type Output = Binary;
    fn mul(self, other: Self) -> Binary {
        if self.is_zero() || other.is_zero() {
            Binary::zero()
        } else {
            Binary {
                mantissa: &self.mantissa * &other.mantissa,
                exponent: self.exponent + other.exponent,
            }
        }
    }
}
impl Neg for &Binary {
    type Output = Binary;
    fn neg(self) -> Binary {
        Binary {
            mantissa: -&self.mantissa,
            exponent: self.exponent,
        }
    }
}
