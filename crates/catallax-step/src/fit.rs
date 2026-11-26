//! Piecewise Taylor polynomials written as degree-seven Bezier spans.
//! Each span carries a uniform Taylor remainder and a convex-hull bound for
//! decimal rounding, enclosure width, and endpoint replacement.
//!
//! The exact rational arithmetic written here is counted as it runs, by the
//! number types themselves, at the decider's limb rate in
//! `OperationWork::uncharged`, and nothing is charged
//! (`catallax-numbers/src/operation.rs`, "The count the number types keep").
//! The fixed control weights are made once for the process and are not
//! counted.

use catallax_intersect::{CurveTrace, TracePiece};
use catallax_numbers::{Ball, Decider, Rational, Refusal, Sign};

use crate::StepOptions;
use crate::number::{self, Point};
use std::collections::BTreeMap;

pub(crate) const DEGREE: usize = 7;

pub(crate) struct Fit {
    pub(crate) points: Vec<Point>,
    pub(crate) degree: usize,
    pub(crate) knots: Vec<Rational>,
    pub(crate) multiplicities: Vec<usize>,
    pub(crate) error: Rational,
    pub(crate) bounds: [[Rational; 2]; 3],
}

struct Span {
    points: Vec<Point>,
    ideals: Vec<[[Rational; 2]; 3]>,
    remainder: Rational,
    bounds: [[Rational; 2]; 3],
}

impl Span {
    fn error(&self) -> Rational {
        let controls = self
            .points
            .iter()
            .zip(&self.ideals)
            .map(|(p, b)| {
                (0..3).fold(Rational::zero(), |s, i| {
                    let value = &p.xyz[i].value;
                    let distance = (value - &b[i][0]).abs().max((value - &b[i][1]).abs());
                    &s + &distance
                })
            })
            .max()
            .unwrap();
        &self.remainder + &controls
    }
}

pub(crate) fn fit(
    trace: &CurveTrace,
    start: &Point,
    end: &Point,
    reverse: bool,
    options: &StepOptions,
    dec: &Decider,
    order: usize,
) -> Result<Fit, Refusal> {
    fit_spans(
        trace,
        start,
        end,
        reverse,
        (options.significant_digits, &options.curve_tolerance, order),
        dec,
        |_, _, _| {},
    )
}

pub(crate) struct TraceSpan {
    pub piece: usize,
    pub from: Rational,
    pub to: Rational,
}

pub(crate) fn fit_with_ranges(
    trace: &CurveTrace,
    start: &Point,
    end: &Point,
    reverse: bool,
    options: &StepOptions,
    dec: &Decider,
    order: usize,
) -> Result<(Fit, Vec<TraceSpan>), Refusal> {
    let mut ranges = Vec::new();
    // A chart projects the written space curve back onto its carrier. Reserve
    // the remaining requested error for that projection and its pcurve fit.
    let tolerance = &options.curve_tolerance / &Rational::from_i64(4);
    let fit = fit_spans(
        trace,
        start,
        end,
        reverse,
        (options.significant_digits, &tolerance, order),
        dec,
        |piece, from, to| {
            ranges.push(TraceSpan {
                piece,
                from: from.clone(),
                to: to.clone(),
            });
        },
    )?;
    if reverse {
        ranges.reverse();
        for range in &mut ranges {
            std::mem::swap(&mut range.from, &mut range.to);
        }
    }
    Ok((fit, ranges))
}

fn fit_spans(
    trace: &CurveTrace,
    start: &Point,
    end: &Point,
    reverse: bool,
    (digits, tolerance, order): (u32, &Rational, usize),
    dec: &Decider,
    mut accepted: impl FnMut(usize, &Rational, &Rational),
) -> Result<Fit, Refusal> {
    dec.check_deadline()?;
    let policy = dec.policy();
    assert!(!trace.pieces().is_empty(), "empty curve trace");
    assert!((1..=DEGREE).contains(&order));
    // TracePiece's internal arithmetic is excluded, as allowed by piece 6's
    // convention. Each span constructs five input balls and one hull.
    let pieces = u64::try_from(trace.pieces().len()).unwrap();
    let initial_operations = pieces
        .checked_mul(6)
        .ok_or_else(|| number::budget("too many initial spline spans"))?;
    let p = number::work_precision(
        policy,
        digits.max(1).saturating_mul(4).saturating_add(32),
        initial_operations,
    )
    .ok_or_else(|| number::budget("no work available for an initial spline fit"))?;
    let cost = u64::from(p.div_ceil(64)).pow(2) * 6;
    // Reserve every initial span before spending anything on refinement.
    let mut remaining = policy.max_work - cost * pieces;
    let mut spans = Vec::new();
    for (piece_index, piece) in trace.pieces().iter().enumerate() {
        let mut endpoints = BTreeMap::new();
        let (lo, hi) = piece.range();
        let whole = Ball::from_rational(&lo, p).hull(&Ball::from_rational(&hi, p));
        if whole.lo() < lo || whole.hi() > hi {
            return Err(number::budget(
                "trace range cannot be enclosed at working precision",
            ));
        }
        let initial = span(
            piece,
            &lo,
            &hi,
            (p, digits, order),
            &whole,
            dec,
            &mut endpoints,
        )?;
        // Its work was reserved above.
        let mut pending = vec![(lo, hi, initial)];
        while let Some((lo, hi, current)) = pending.pop() {
            dec.check_deadline()?;
            let error = current.error();
            // Once decimal/enclosure error dominates, subdivision cannot help.
            let useful = (&current.remainder * &Rational::from_i64(4))
                > error.clone().max(tolerance.clone());
            if &error > tolerance && useful && remaining >= cost.saturating_mul(2) {
                remaining -= cost * 2;
                let mid = (&lo + &hi).mul_pow2(-1);
                let left = span(
                    piece,
                    &lo,
                    &mid,
                    (p, digits, order),
                    &whole,
                    dec,
                    &mut endpoints,
                )?;
                let right = span(
                    piece,
                    &mid,
                    &hi,
                    (p, digits, order),
                    &whole,
                    dec,
                    &mut endpoints,
                )?;
                if left.error().max(right.error()) < error {
                    pending.push((mid.clone(), hi, right));
                    pending.push((lo, mid, left));
                } else {
                    accepted(piece_index, &lo, &hi);
                    spans.push(current);
                }
            } else {
                accepted(piece_index, &lo, &hi);
                spans.push(current);
            }
        }
    }
    let (first, last) = if reverse { (end, start) } else { (start, end) };
    spans.first_mut().unwrap().points[0] = first.clone();
    spans.last_mut().unwrap().points[DEGREE] = last.clone();
    let mut points: Vec<Point> = Vec::with_capacity(spans.len() * DEGREE + 1);
    let mut error = Rational::zero();
    let mut bounds = spans[0].bounds.clone();
    for mut span in spans {
        if let Some(previous) = points.last() {
            span.points[0] = previous.clone();
        }
        error = error.max(span.error());
        for (b, other) in bounds.iter_mut().zip(&span.bounds) {
            b[0] = b[0].clone().min(other[0].clone());
            b[1] = b[1].clone().max(other[1].clone());
        }
        let skip = usize::from(!points.is_empty());
        points.extend(span.points.into_iter().skip(skip));
    }
    if reverse {
        points.reverse();
    }
    let count = (points.len() - 1) / DEGREE;
    let knots = knots(count, digits);
    let multiplicities = (0..=count)
        .map(|i| {
            if i == 0 || i == count {
                DEGREE + 1
            } else {
                DEGREE
            }
        })
        .collect();
    Ok(Fit {
        points,
        degree: DEGREE,
        knots,
        multiplicities,
        error,
        bounds,
    })
}

fn span(
    piece: &TracePiece,
    lo: &Rational,
    hi: &Rational,
    (precision, digits, order): (u32, u32, usize),
    whole: &Ball,
    dec: &Decider,
    endpoints: &mut BTreeMap<Rational, [Ball; 3]>,
) -> Result<Span, Refusal> {
    let mid = (lo + hi).mul_pow2(-1);
    let width = hi - lo;
    let lower = whole.lo();
    let upper = whole.hi();
    let inside = |ball: &Ball| ball.lo() >= lower && ball.hi() <= upper;
    let raw_center = Ball::from_rational(&mid, precision);
    let center = if inside(&raw_center) {
        &raw_center
    } else {
        whole
    };
    // Outward rounding at low precision may leave a piece's range. Fall
    // back to the certified whole range rather than evaluate beyond it.
    let candidate = (lo != &lower || hi != &upper)
        .then(|| Ball::from_rational(lo, precision).hull(&Ball::from_rational(hi, precision)));
    let domain = candidate
        .as_ref()
        .filter(|ball| inside(ball))
        .unwrap_or(whole);
    let mut coefficients = piece.taylor(center, order, dec)?;
    let remainder_coefficient = piece.taylor(domain, order + 1, dec)?.pop().unwrap();
    let magnitude = remainder_coefficient
        .iter()
        .fold(Rational::zero(), |s, b| &s + &number::magnitude(b));
    let remainder = &magnitude * &power(&width.mul_pow2(-1), order + 1);
    // Elevate a lower-order polynomial into the existing degree-seven format.
    coefficients.resize_with(DEGREE + 1, || {
        std::array::from_fn(|_| Ball::from_rational(&Rational::zero(), precision))
    });
    let mut ideals = Vec::with_capacity(DEGREE + 1);
    let mut points = Vec::with_capacity(DEGREE + 1);
    taylor_controls(&coefficients, &width, digits, |xyz, control| {
        points.push(Point { xyz });
        ideals.push(control);
    });
    for (index, t) in [(0, lo), (DEGREE, hi)] {
        let raw = Ball::from_rational(t, precision);
        let parameter = if inside(&raw) { &raw } else { whole };
        let position = if order < DEGREE {
            if let Some(position) = endpoints.get(t) {
                position.clone()
            } else {
                let position = piece.evaluate(parameter, dec)?;
                endpoints.insert(t.clone(), position.clone());
                position
            }
        } else {
            piece.evaluate(parameter, dec)?
        };
        points[index] = Point {
            xyz: position.each_ref().map(|b| number::ball_number(b, digits)),
        };
    }
    let bounds = piece
        .evaluate(domain, dec)?
        .each_ref()
        .map(|b| [b.lo(), b.hi()]);
    Ok(Span {
        points,
        ideals,
        remainder,
        bounds,
    })
}

/// The knots of `count` spans, from zero. A repeated Bezier knot isolates
/// each span's parameterization. Each knot is the next exactly representable
/// decimal past the one before, so low digit counts cannot collapse distinct
/// knots and change the fitted curve. Each costs one sum and one writing in
/// decimal.
fn knots(count: usize, digits: u32) -> Vec<Rational> {
    let mut knots = Vec::with_capacity(count + 1);
    let mut knot = Rational::zero();
    knots.push(knot.clone());
    for _ in 0..count {
        let next = (&knot + &Rational::one())
            .to_decimal_bounds(digits.max(1))
            .1;
        knot = number::formatted_value(&next);
        knots.push(knot.clone());
    }
    knots
}

/// Convert an enclosed Taylor polynomial on a centered interval to rounded
/// Bezier controls, retaining the unrounded enclosure for its error proof.
pub(crate) fn taylor_controls<const N: usize>(
    coefficients: &[[Ball; N]],
    width: &Rational,
    digits: u32,
    mut push: impl FnMut([number::Number; N], [[Rational; 2]; N]),
) {
    use std::sync::LazyLock;

    // Reflection gives w[n-j,k] = (-1)^k w[j,k]. Evaluate the even and odd
    // terms once for each pair of controls, retaining exact interval endpoints.
    // Made once for the process, so not counted (`uncounted_memo`).
    static WEIGHTS: LazyLock<[[Rational; DEGREE + 1]; DEGREE / 2 + 1]> = LazyLock::new(|| {
        catallax_numbers::uncounted_memo(|| {
            std::array::from_fn(|j| {
                std::array::from_fn(|k| {
                    (0..=k.min(j)).fold(Rational::zero(), |sum, l| {
                        sum + Rational::from_ratio(
                            binomial(k, l) * binomial(j, l),
                            binomial(DEGREE, l),
                        )
                        .unwrap()
                            * Rational::from_ratio(-1, 2)
                                .unwrap()
                                .pow((k - l) as i32)
                                .unwrap()
                    })
                })
            })
        })
    });
    assert_eq!(coefficients.len(), DEGREE + 1);
    let zeros = || -> [[[Rational; 2]; N]; DEGREE + 1] {
        std::array::from_fn(|_| std::array::from_fn(|_| [Rational::zero(), Rational::zero()]))
    };
    // Scale each coefficient enclosure once, before applying the control weights.
    // All operations are rational, so reassociation introduces no rounding.
    // Each power of the width is one product more than the last.
    let mut scaled = zeros();
    let mut power = Rational::one();
    for (k, (row, coefficient)) in scaled.iter_mut().zip(coefficients).enumerate() {
        if k != 0 {
            power = &power * width;
        }
        for (cell, ball) in row.iter_mut().zip(coefficient) {
            if !ball.is_exact_zero() {
                *cell = [&ball.lo() * &power, &ball.hi() * &power];
                if power.sign() == Sign::Negative {
                    cell.swap(0, 1);
                }
            }
        }
    }
    let mut controls = zeros();
    for (j, weights) in WEIGHTS.iter().enumerate() {
        let mut parts: [[[Rational; 2]; N]; 2] =
            std::array::from_fn(|_| std::array::from_fn(|_| [Rational::zero(), Rational::zero()]));
        for (k, (coefficient, weight)) in scaled.iter().zip(weights).enumerate() {
            let [lo, hi] = match weight.sign() {
                Sign::Zero => continue,
                Sign::Positive => [0, 1],
                Sign::Negative => [1, 0],
            };
            for (interval, coefficient) in parts[k % 2].iter_mut().zip(coefficient) {
                if !coefficient[0].is_zero() || !coefficient[1].is_zero() {
                    interval[0] = &interval[0] + &(&coefficient[lo] * weight);
                    interval[1] = &interval[1] + &(&coefficient[hi] * weight);
                }
            }
        }
        for i in 0..N {
            let (even, odd) = (&parts[0][i], &parts[1][i]);
            controls[j][i] = [&even[0] + &odd[0], &even[1] + &odd[1]];
            if DEGREE - j != j {
                controls[DEGREE - j][i] = [&even[0] - &odd[1], &even[1] - &odd[0]];
            }
        }
    }
    for control in controls {
        push(
            control
                .each_ref()
                .map(|b| number::rounded(&b[0], &b[1], digits)),
            control,
        );
    }
}

/// `base` to the power `exponent`, as `exponent` products.
fn power(base: &Rational, exponent: usize) -> Rational {
    let mut power = Rational::one();
    for _ in 0..exponent {
        power = &power * base;
    }
    power
}

pub(crate) fn bezier_power<const N: usize>(
    mut differences: [[Rational; N]; DEGREE + 1],
) -> [[Rational; N]; DEGREE + 1] {
    let mut powers = std::array::from_fn(|_| std::array::from_fn(|_| Rational::zero()));
    for (k, power) in powers.iter_mut().enumerate() {
        let factor = Rational::from_i64(binomial(DEGREE, k));
        for (value, difference) in power.iter_mut().zip(&differences[0]) {
            *value = difference * &factor;
        }
        for j in 0..DEGREE - k {
            differences[j] =
                std::array::from_fn(|axis| &differences[j + 1][axis] - &differences[j][axis]);
        }
    }
    powers
}

fn binomial(n: usize, k: usize) -> i64 {
    (0..k).fold(1, |v, i| v * (n - i) as i64 / (i + 1) as i64)
}

#[cfg(test)]
#[path = "fit/tests.rs"]
mod tests;
