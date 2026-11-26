//! Eight-point Gauss integration with certified roots, weights and remainder.

use super::{Ball, Checked, Coefficients, Program, Rational, Refusal, Stop, Work, turning};
use catallax_numbers::{Integer, Sign};

mod tensor;
pub(crate) use tensor::NativeArea;

pub(super) const ORDER: usize = 16;
const SIZE: usize = ORDER + 1;
// 128 P_8, in descending powers. Its four positive roots and their negatives
// are all eight roots, certified below by disjoint sign-changing intervals.
const POLYNOMIAL: [i64; 9] = [6435, 0, -12012, 0, 6930, 0, -1260, 0, 35];
const DERIVATIVE: [i64; 8] = [51480, 0, -72072, 0, 27720, 0, -2520, 0];

pub(crate) struct Rule {
    points: Vec<(Ball, Ball)>,
    native: std::sync::OnceLock<Vec<(catallax_numbers::Interval, catallax_numbers::Interval)>>,
}

fn certificate() -> Stop {
    Stop::Refusal(Refusal::new(
        "measure.quadrature-certificate",
        "Legendre root intervals do not certify eight distinct roots",
    ))
}

fn limbs(numerator: &Integer, denominator: &Integer) -> u64 {
    numerator.bit_length().div_ceil(64).max(1) + denominator.bit_length().div_ceil(64).max(1)
}

fn charge_exact(work: &mut Work, cost: u64) -> Checked<()> {
    let unit = u64::from(work.precision.div_ceil(64)).pow(2);
    // Keep the original price of each exact rational operation.
    for _ in 0..cost.div_ceil(unit).max(1) {
        work.charge()?;
    }
    Ok(())
}

fn polynomial(x: &Rational, work: &mut Work) -> Checked<Sign> {
    const SMALL: u64 = 2; // Two one-limb integers, including zero's numerator.
    let x_num = x.numerator();
    let x_den = x.denominator();
    let x_limbs = limbs(x_num, x_den);
    charge_exact(work, SMALL * x_limbs)?;
    charge_exact(work, SMALL + SMALL)?;
    charge_exact(work, SMALL * x_limbs)?;
    let first = &Rational::from_i64(POLYNOMIAL[0]) * x;
    charge_exact(work, limbs(first.numerator(), first.denominator()) + SMALL)?;

    // Brackets and their midpoints have only 2 and 5 in the denominator.
    // The leading coefficient has one factor of 5, and the next coefficient
    // is zero. After the first product its numerator is coprime to x_den.
    // The remaining Horner denominators include x_den, so adding an integer
    // cannot introduce a common factor. There is nothing left to reduce.
    let mut numerator = first.numerator().clone();
    let mut denominator = first.denominator().clone();
    for coefficient in POLYNOMIAL.into_iter().skip(2) {
        charge_exact(work, limbs(&numerator, &denominator) * x_limbs)?;
        numerator = &numerator * x_num;
        denominator = &denominator * x_den;
        charge_exact(work, limbs(&numerator, &denominator) + SMALL)?;
        if coefficient != 0 {
            numerator = numerator + Integer::from_i64(coefficient) * &denominator;
        }
    }
    Ok(numerator.sign())
}

impl Rule {
    fn native_points(&self) -> &[(catallax_numbers::Interval, catallax_numbers::Interval)] {
        use catallax_numbers::Interval;
        self.native.get_or_init(|| {
            let interval = |a: &Ball| {
                Interval::new(
                    Interval::from_rational(&a.lo()).lo(),
                    Interval::from_rational(&a.hi()).hi(),
                )
                .unwrap()
            };
            self.points
                .iter()
                .map(|(x, w)| (interval(x), interval(w)))
                .collect()
        })
    }

    pub(crate) fn new(work: &mut Work) -> Checked<Self> {
        let threshold = Rational::one().mul_pow2(-i64::from(work.precision) - 16);
        let mut points = Vec::with_capacity(8);
        let mut previous = Rational::zero();
        for (lower, upper) in super::gauss_nodes::NODE_INTERVALS {
            let mut low = Rational::parse(lower).expect("generated rational endpoint");
            let mut high = Rational::parse(upper).expect("generated rational endpoint");
            if low < previous || high > Rational::one() || low >= high {
                return Err(certificate());
            }
            previous = high.clone();
            let lower_sign = polynomial(&low, work)?;
            let upper_sign = polynomial(&high, work)?;
            if lower_sign == Sign::Zero || upper_sign == Sign::Zero || lower_sign == upper_sign {
                return Err(certificate());
            }
            while &high - &low > threshold {
                let middle = (&low + &high).mul_pow2(-1);
                let sign = polynomial(&middle, work)?;
                if sign == Sign::Zero {
                    low = middle.clone();
                    high = middle;
                    break;
                }
                if sign == lower_sign {
                    low = middle;
                } else {
                    high = middle;
                }
            }
            let point = work.interval(&low, &high)?;
            let mut derivative = work.rational(&Rational::zero())?;
            for coefficient in DERIVATIVE {
                let product = work.mul(&derivative, &point)?;
                let coefficient = work.rational(&Rational::from_i64(coefficient))?;
                derivative = work.add(&product, &coefficient)?;
            }
            let one = work.rational(&Rational::one())?;
            let square = work.mul(&point, &point)?;
            let gap = work.sub(&one, &square)?;
            let derivative_squared = work.mul(&derivative, &derivative)?;
            let denominator = work.mul(&gap, &derivative_squared)?;
            // Q = 128 P_8, so 2 / ((1-x^2) P'_8(x)^2) has numerator 32768.
            let numerator = work.rational(&Rational::from_i64(32768))?;
            let weight = work.div(&numerator, &denominator)?;
            if weight.lo() <= Rational::zero() {
                return Err(Stop::Split);
            }
            let zero = work.rational(&Rational::zero())?;
            points.push((work.sub(&zero, &point)?, weight.clone()));
            points.push((point, weight));
        }
        Ok(Self {
            points,
            native: std::sync::OnceLock::new(),
        })
    }
}

impl Program {
    /// Integrates f(cos(pi x), sin(pi x)). The Gauss remainder is the enclosed
    /// coefficient f^(16)/16! times the squared monic polynomial's integral.
    /// On a panel of width h that factor is h^17 / (17 binom(16,8)^2).
    /// See https://dlmf.nist.gov/3.5.E19.
    pub(crate) fn turning_gauss_panel(
        &self,
        a: &Rational,
        b: &Rational,
        rule: &Rule,
        work: &mut Work,
    ) -> Checked<Vec<Ball>> {
        if a >= b {
            return Err(certificate());
        }
        let midpoint = (a + b).mul_pow2(-1);
        let width = b - a;
        let center = work.rational(&midpoint)?;
        let half = work.rational(&width.mul_pow2(-1))?;
        let zero = work.rational(&Rational::zero())?;
        let mut sums = vec![zero; self.outputs.len()];
        for (point, weight) in &rule.points {
            let offset = work.mul(&half, point)?;
            let x = work.add(&center, &offset)?;
            let input = turning::<1>(&x, work)?;
            let values = self.evaluate(&input, Coefficients::All, work)?;
            for (sum, value) in sums.iter_mut().zip(values) {
                let term = work.mul(&value[0], weight)?;
                *sum = work.add(sum, &term)?;
            }
        }
        let range = work.interval(a, b)?;
        let input = turning::<SIZE>(&range, work)?;
        let remainder = self.evaluate(&input, Coefficients::GaussRemainder, work)?;
        let coefficient =
            width.pow(17).expect("positive exponent") / Rational::from_i64(17 * 12870 * 12870);
        let coefficient = work.rational(&coefficient)?;
        sums.iter()
            .zip(remainder)
            .map(|(sum, remainder)| {
                let quadrature = work.mul(&half, sum)?;
                let error = work.mul(&coefficient, &remainder[ORDER])?;
                work.add(&quadrature, &error)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integrate::Builder;
    use catallax_numbers::{Decider, with_operation_budget};

    #[test]
    fn integer_root_sign_keeps_the_rational_budget_at_generated_brackets() {
        let dec = Decider::default();
        let policy = *dec.policy();
        for (lo, hi) in super::super::gauss_nodes::NODE_INTERVALS {
            let lo = Rational::parse(lo).unwrap();
            let hi = Rational::parse(hi).unwrap();
            for x in [&lo, &hi, &(&lo + &hi).mul_pow2(-1)] {
                let mut old = Work::new(&policy, 128, policy.max_work, &dec).unwrap();
                let mut new = Work::new(&policy, 128, policy.max_work, &dec).unwrap();
                let mut value = Rational::zero();
                for c in POLYNOMIAL {
                    let x_limbs = limbs(x.numerator(), x.denominator());
                    charge_exact(
                        &mut old,
                        limbs(value.numerator(), value.denominator()) * x_limbs,
                    )
                    .unwrap();
                    value = &value * x;
                    charge_exact(&mut old, limbs(value.numerator(), value.denominator()) + 2)
                        .unwrap();
                    value = value + Rational::from_i64(c);
                }
                assert_eq!(polynomial(x, &mut new).unwrap(), value.sign());
                assert_eq!(new.remaining(), old.remaining());
            }
        }
    }

    #[test]
    fn certified_rule_integrates_monomials_and_its_first_remainder() {
        let ((), used) = with_operation_budget(1_000_000, || {
            let dec = Decider::default();
            let policy = *dec.policy();
            let mut work = Work::new(&policy, 96, policy.max_work, &dec).unwrap();
            let rule = Rule::new(&mut work).unwrap();
            for degree in 0..=16 {
                let mut sum = work.rational(&Rational::zero()).unwrap();
                for (point, weight) in &rule.points {
                    let mut power = work.rational(&Rational::one()).unwrap();
                    for _ in 0..degree {
                        power = work.mul(&power, point).unwrap();
                    }
                    let term = work.mul(&power, weight).unwrap();
                    sum = work.add(&sum, &term).unwrap();
                }
                if degree == 16 {
                    let gamma = Rational::from_i64(2).pow(17).unwrap()
                        / Rational::from_i64(17 * 12870 * 12870);
                    let gamma = work.rational(&gamma).unwrap();
                    sum = work.add(&sum, &gamma).unwrap();
                }
                let exact = if degree % 2 == 0 {
                    Rational::from_ratio(2, degree + 1).unwrap()
                } else {
                    Rational::zero()
                };
                assert!(sum.lo() <= exact && exact <= sum.hi(), "degree {degree}");
            }
        });
        assert!(!used.exhausted);
    }

    #[test]
    fn turning_rule_includes_the_signed_remainder() {
        let ((), used) = with_operation_budget(1_000_000, || {
            let dec = Decider::default();
            let policy = *dec.policy();
            let mut work = Work::new(&policy, 96, policy.max_work, &dec).unwrap();
            let rule = Rule::new(&mut work).unwrap();
            let builder = Builder::default();
            let cosine = builder.input(0);
            let program = builder.finish(&[cosine.square()]);
            let low = Rational::from_ratio(-1, 2).unwrap();
            let high = Rational::from_ratio(1, 2).unwrap();
            let result = program
                .turning_gauss_panel(&low, &high, &rule, &mut work)
                .unwrap();
            let exact = Rational::from_ratio(1, 2).unwrap();
            assert!(result[0].lo() <= exact && exact <= result[0].hi());
            assert!(result[0].hi() - result[0].lo() < Rational::parse("1e-6").unwrap());
        });
        assert!(!used.exhausted);
    }
}
