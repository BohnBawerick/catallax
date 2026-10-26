//! Parameter-coupled degree-seven UV fits, with outward Taylor bounds.

use catallax_numbers::{Ball, Decider, Rational, Refusal};

use crate::{fit, number};
use number::Number;

/// Taylor coefficients, not derivatives. Constant and polynomial jets retain
/// their degree so affine coordinate changes do not multiply known zero tails.
#[derive(Clone)]
pub(crate) struct Jet<const N: usize> {
    pub coefficients: [Ball; N],
    degree: usize,
}

fn scalar(value: i64, precision: u32) -> Ball {
    Ball::from_rational(&Rational::from_i64(value), precision)
}

impl<const N: usize> Jet<N> {
    pub fn constant(value: Ball) -> Self {
        let precision = value.precision();
        Self {
            coefficients: std::array::from_fn(|i| {
                if i == 0 {
                    value.clone()
                } else {
                    scalar(0, precision)
                }
            }),
            degree: 0,
        }
    }

    pub fn variable(value: Ball) -> Self {
        let precision = value.precision();
        let mut result = Self::constant(value);
        if N > 1 {
            result.coefficients[1] = scalar(1, precision);
            result.degree = 1;
        }
        result
    }

    pub fn from_coefficients(coefficients: [Ball; N]) -> Self {
        Self {
            coefficients,
            degree: N - 1,
        }
    }

    pub fn add(&self, other: &Self) -> Self {
        let mut result = Self::constant(scalar(0, self.coefficients[0].precision()));
        result.degree = self.degree.max(other.degree);
        for i in 0..=result.degree {
            result.coefficients[i] = &self.coefficients[i] + &other.coefficients[i];
        }
        result
    }

    pub fn sub(&self, other: &Self) -> Self {
        let mut result = Self::constant(scalar(0, self.coefficients[0].precision()));
        result.degree = self.degree.max(other.degree);
        for i in 0..=result.degree {
            result.coefficients[i] = &self.coefficients[i] - &other.coefficients[i];
        }
        result
    }

    pub fn scale(&self, factor: &Ball) -> Self {
        let mut result = Self::constant(scalar(0, self.coefficients[0].precision()));
        if factor.is_exact_zero() {
            return result;
        }
        result.degree = self.degree;
        for i in 0..=self.degree {
            result.coefficients[i] = &self.coefficients[i] * factor;
        }
        result
    }

    pub fn mul(&self, other: &Self) -> Self {
        let mut result = Self::constant(scalar(0, self.coefficients[0].precision()));
        result.degree = (self.degree + other.degree).min(N - 1);
        for n in 0..=result.degree {
            for k in n.saturating_sub(other.degree)..=n.min(self.degree) {
                result.coefficients[n] =
                    &result.coefficients[n] + &(&self.coefficients[k] * &other.coefficients[n - k]);
            }
        }
        result
    }

    pub fn square(&self) -> Self {
        let precision = self.coefficients[0].precision();
        let mut result = Self::constant(self.coefficients[0].square());
        result.degree = (self.degree * 2).min(N - 1);
        for n in 1..=result.degree {
            for k in n.saturating_sub(self.degree)..=self.degree.min(n / 2) {
                let other = n - k;
                let term = if k == other {
                    self.coefficients[k].square()
                } else {
                    &(&self.coefficients[k] * &self.coefficients[other]) * &scalar(2, precision)
                };
                result.coefficients[n] = &result.coefficients[n] + &term;
            }
        }
        result
    }

    pub fn sqrt(&self) -> Result<Self, Refusal> {
        if self.coefficients[0].lo() <= Rational::zero() {
            return Err(number::budget("UV radial jet is not separated from a pole"));
        }
        let root = self.coefficients[0]
            .sqrt()
            .expect("positive UV radial square");
        let inverse = (&root * &scalar(2, root.precision()))
            .recip()
            .ok_or_else(|| number::budget("UV radial derivative is not separated"))?;
        let mut result = Self::constant(root);
        if self.degree == 0 {
            return Ok(result);
        }
        result.degree = N - 1;
        for n in 1..N {
            let mut residual = self.coefficients[n].clone();
            for k in 1..n {
                residual = residual - &result.coefficients[k] * &result.coefficients[n - k];
            }
            result.coefficients[n] = &residual * &inverse;
        }
        Ok(result)
    }

    /// Integrate d atan2(y,x) = (x dy - y dx)/(x²+y²). The caller selects
    /// the continuous angular lift for coefficient zero; derivatives do not
    /// depend on that lift and need no transcendental series for remainder work.
    pub fn angle(y: &Self, x: &Self, constant: Ball) -> Result<Self, Refusal> {
        let mut result = Self::constant(constant);
        if N == 1 || (x.degree == 0 && y.degree == 0) {
            return Ok(result);
        }
        let precision = result.coefficients[0].precision();
        let mut denominator: [Ball; N] = std::array::from_fn(|_| scalar(0, precision));
        let mut quotient: [Ball; N] = std::array::from_fn(|_| scalar(0, precision));
        let inverse = (x.coefficients[0].square() + y.coefficients[0].square())
            .recip()
            .ok_or_else(|| number::budget("UV angular jet is not separated from its axis"))?;
        for n in 0..N - 1 {
            let mut numerator = scalar(0, precision);
            for k in 0..=n {
                if n != 0 && k <= x.degree && n - k <= x.degree {
                    denominator[n] =
                        &denominator[n] + &(&x.coefficients[k] * &x.coefficients[n - k]);
                }
                if n != 0 && k <= y.degree && n - k <= y.degree {
                    denominator[n] =
                        &denominator[n] + &(&y.coefficients[k] * &y.coefficients[n - k]);
                }
                let order = n - k + 1;
                let factor = scalar(order as i64, precision);
                if k <= x.degree && order <= y.degree {
                    numerator =
                        numerator + &(&x.coefficients[k] * &y.coefficients[order]) * &factor;
                }
                if k <= y.degree && order <= x.degree {
                    numerator =
                        numerator - &(&y.coefficients[k] * &x.coefficients[order]) * &factor;
                }
            }
            for k in 1..=n {
                numerator = numerator - &denominator[k] * &quotient[n - k];
            }
            quotient[n] = &numerator * &inverse;
            result.coefficients[n + 1] = &quotient[n]
                * &Ball::from_rational(
                    &Rational::from_ratio(1, (n + 1) as i64).unwrap(),
                    precision,
                );
        }
        result.degree = N - 1;
        Ok(result)
    }
}

/// A local curve has an actual STEP parameter interval. Unlike the existing
/// space-curve fit, UV refinement must never reparameterize its Bezier spans.
pub(crate) trait Parametric {
    fn range(&self) -> [Rational; 2];
    /// `false` requests interval subdivision before evaluating derivatives;
    /// arithmetic refusals still propagate without being swallowed or retried.
    fn regular(&self, _parameter: &Ball, _dec: &Decider) -> Result<bool, Refusal> {
        Ok(true)
    }
    fn coefficients<const N: usize>(
        &self,
        parameter: &Ball,
        constants: bool,
        dec: &Decider,
    ) -> Result<[[Ball; 2]; N], Refusal>;
    fn value(
        &self,
        parameter: &Rational,
        precision: u32,
        dec: &Decider,
    ) -> Result<[Ball; 2], Refusal> {
        Ok(
            self.coefficients::<1>(&Ball::from_rational(parameter, precision), true, dec)?[0]
                .clone(),
        )
    }
}

pub(crate) struct Fit {
    pub points: Vec<[Number; 2]>,
    pub knots: Vec<Rational>,
    pub multiplicities: Vec<usize>,
    pub error: Rational,
}

struct Span {
    range: [Rational; 2],
    points: Vec<[Number; 2]>,
    ideals: Vec<[[Rational; 2]; 2]>,
    remainder: Rational,
}

impl Span {
    fn error(&self) -> Rational {
        self.remainder.clone()
            + self
                .points
                .iter()
                .zip(&self.ideals)
                .map(|(p, bounds)| {
                    (0..2).fold(Rational::zero(), |sum, i| {
                        sum + (&p[i].value - &bounds[i][0])
                            .abs()
                            .max((&p[i].value - &bounds[i][1]).abs())
                    })
                })
                .max()
                .expect("nonempty UV controls")
    }
}

fn span<P: Parametric>(
    piece: &P,
    range: [Rational; 2],
    precision: u32,
    digits: u32,
    domain: &Ball,
    dec: &Decider,
) -> Result<Span, Refusal> {
    let [lo, hi] = &range;
    let width = hi - lo;
    let center = Ball::from_rational(&(lo + hi).mul_pow2(-1), precision);
    let coefficients = piece.coefficients::<{ fit::DEGREE + 1 }>(&center, true, dec)?;
    let remainder_coefficient = piece.coefficients::<{ fit::DEGREE + 2 }>(domain, false, dec)?;
    let remainder = remainder_coefficient[fit::DEGREE + 1]
        .iter()
        .fold(Rational::zero(), |sum, value| {
            sum + number::magnitude(value)
        })
        * width.mul_pow2(-1).pow((fit::DEGREE + 1) as i32).unwrap();
    let mut points = Vec::with_capacity(fit::DEGREE + 1);
    let mut ideals = Vec::with_capacity(fit::DEGREE + 1);
    fit::taylor_controls(&coefficients, &width, digits, |point, ideal| {
        points.push(point);
        ideals.push(ideal);
    });
    for (index, parameter) in [(0, lo), (fit::DEGREE, hi)] {
        points[index] = piece
            .value(parameter, precision, dec)?
            .each_ref()
            .map(|value| number::ball_number(value, digits));
    }
    Ok(Span {
        range,
        points,
        ideals,
        remainder,
    })
}

pub(crate) struct Budget {
    pub precision: u32,
    pieces: usize,
    cost: u64,
    remaining: u64,
}

impl Budget {
    pub fn new(pieces: usize, digits: u32, dec: &Decider) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        assert!(pieces != 0, "empty UV trace");
        // Two coordinate jets, two endpoints and Taylor-to-Bezier controls.
        // The request-wide arithmetic meter remains active as well.
        const OPERATIONS: u64 = 32_768;
        let count = u64::try_from(pieces).unwrap();
        let initial = count
            .checked_mul(OPERATIONS)
            .ok_or_else(|| number::budget("too many initial UV spans"))?;
        let precision = number::work_precision(
            dec.policy(),
            digits.max(1).saturating_mul(4).saturating_add(32),
            initial,
        )
        .ok_or_else(|| number::budget("no work available for an initial UV fit"))?;
        let cost = u64::from(precision.div_ceil(64)).pow(2) * OPERATIONS;
        Ok(Self {
            precision,
            pieces,
            cost,
            remaining: dec.policy().max_work - cost * count,
        })
    }
}

pub(crate) fn fit<P: Parametric>(
    pieces: &[P],
    digits: u32,
    tolerance: &Rational,
    budget: Budget,
    dec: &Decider,
) -> Result<Fit, Refusal> {
    dec.check_deadline()?;
    assert_eq!(pieces.len(), budget.pieces, "changed UV reservation");
    let Budget {
        precision,
        cost,
        mut remaining,
        ..
    } = budget;
    let mut spans = Vec::new();
    for piece in pieces {
        let mut pending = vec![piece.range()];
        while let Some(range) = pending.pop() {
            dec.check_deadline()?;
            let domain = Ball::from_rational(&range[0], precision)
                .hull(&Ball::from_rational(&range[1], precision));
            let mut current = None;
            let refine = if piece.regular(&domain, dec)? {
                let span = span(piece, range.clone(), precision, digits, &domain, dec)?;
                let error = span.error();
                let useful =
                    &span.remainder * &Rational::from_i64(4) > error.clone().max(tolerance.clone());
                current = Some(span);
                error > *tolerance && useful
            } else {
                true
            };
            if refine && remaining >= cost.saturating_mul(2) {
                let middle = number::rational(&(&range[0] + &range[1]).mul_pow2(-1), digits).value;
                if middle > range[0] && middle < range[1] {
                    remaining -= cost * 2;
                    pending.push([middle.clone(), range[1].clone()]);
                    pending.push([range[0].clone(), middle]);
                    continue;
                }
            }
            spans.push(current.ok_or_else(|| {
                number::budget("UV chart domain cannot be separated from its axis")
            })?);
        }
    }
    let mut points: Vec<[Number; 2]> = Vec::with_capacity(spans.len() * fit::DEGREE + 1);
    let mut knots = Vec::with_capacity(spans.len() + 1);
    let mut error = Rational::zero();
    for mut span in spans {
        if let Some(previous) = points.last() {
            span.points[0] = previous.clone();
            assert_eq!(
                knots.last(),
                Some(&span.range[0]),
                "disjoint UV parameter pieces"
            );
        } else {
            knots.push(span.range[0].clone());
        }
        error = error.max(span.error());
        let skip = usize::from(!points.is_empty());
        points.extend(span.points.into_iter().skip(skip));
        let [_, hi] = span.range;
        knots.push(hi);
    }
    let multiplicities = (0..knots.len())
        .map(|i| {
            if i == 0 || i + 1 == knots.len() {
                fit::DEGREE + 1
            } else {
                fit::DEGREE
            }
        })
        .collect();
    Ok(Fit {
        points,
        knots,
        multiplicities,
        error,
    })
}

/// A locally continuous atan2 branch. Periodic callers choose its integer lift.
pub(crate) fn angle(y: &Ball, x: &Ball, pi: &Ball, dec: &Decider) -> Result<Ball, Refusal> {
    dec.check_deadline()?;
    if x.lo() > Rational::zero() {
        return y.checked_div(x).unwrap().atan_checked(dec);
    }
    if x.hi() < Rational::zero() {
        return Ok(pi - &y.checked_div(&(-x)).unwrap().atan_checked(dec)?);
    }
    if y.lo() > Rational::zero() {
        return Ok(pi * &scalar_ratio(1, 2, pi.precision())
            - x.checked_div(y).unwrap().atan_checked(dec)?);
    }
    if y.hi() < Rational::zero() {
        return Ok(-(pi * &scalar_ratio(1, 2, pi.precision()))
            - x.checked_div(y).unwrap().atan_checked(dec)?);
    }
    Err(number::budget("UV angle is not separated from its axis"))
}

pub(crate) fn positive_angle(
    y: &Ball,
    x: &Ball,
    pi: &Ball,
    dec: &Decider,
) -> Result<Ball, Refusal> {
    dec.check_deadline()?;
    if x.lo() > Rational::zero() {
        if y.lo() >= Rational::zero() {
            return y.checked_div(x).unwrap().atan_checked(dec);
        }
        if y.hi() < Rational::zero() {
            return Ok(
                pi * &scalar(2, pi.precision()) + y.checked_div(x).unwrap().atan_checked(dec)?
            );
        }
    }
    if x.hi() < Rational::zero() {
        return Ok(pi - &y.checked_div(&(-x)).unwrap().atan_checked(dec)?);
    }
    if y.lo() > Rational::zero() {
        return Ok(pi * &scalar_ratio(1, 2, pi.precision())
            - x.checked_div(y).unwrap().atan_checked(dec)?);
    }
    if y.hi() < Rational::zero() {
        return Ok(pi * &scalar_ratio(3, 2, pi.precision())
            + x.checked_div(&(-y)).unwrap().atan_checked(dec)?);
    }
    Err(number::budget(
        "UV angular lift is not separated from its cut",
    ))
}

fn scalar_ratio(numerator: i64, denominator: i64, precision: u32) -> Ball {
    Ball::from_rational(
        &Rational::from_ratio(numerator, denominator).unwrap(),
        precision,
    )
}
