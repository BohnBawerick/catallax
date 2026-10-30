//! Bounded joint fits at a pinched regular spherical chart corner.

use crate::{fit, number, uv};
use catallax_numbers::{Decider, Rational, Refusal};

type Point = [Rational; 2];
type Controls = [Point; fit::DEGREE + 1];
const GRAPH_DEGREE: usize = 4;
const RESIDUAL_DEGREE: usize = GRAPH_DEGREE * fit::DEGREE;

pub(crate) struct Repair {
    pub start: bool,
    pub points: Vec<[number::Number; 2]>,
    pub knot: Rational,
    pub displacement: Rational,
}

impl Repair {
    pub fn apply(mut self, curve: &mut uv::Fit) {
        let index = if self.start { 1 } else { curve.knots.len() - 1 };
        let range = if self.start {
            0..fit::DEGREE + 1
        } else {
            self.points.reverse();
            curve.points.len() - fit::DEGREE - 1..curve.points.len()
        };
        curve.points.splice(range, self.points);
        curve.knots.insert(index, self.knot);
        curve.multiplicities.insert(index, fit::DEGREE);
        curve.error = &curve.error + &self.displacement;
    }
}

struct Frame {
    origin: Point,
    axis: usize,
    positive: bool,
    max_depth: u32,
}

impl Frame {
    fn point(&self, point: &Point) -> Point {
        let x = &point[self.axis] - &self.origin[self.axis];
        let y = &point[1 - self.axis] - &self.origin[1 - self.axis];
        [
            if self.positive { x } else { -x },
            if self.positive == (self.axis == 0) {
                y
            } else {
                -y
            },
        ]
    }

    fn powers(&self, controls: &Controls) -> [[Rational; 2]; fit::DEGREE + 1] {
        fit::bezier_power(controls.each_ref().map(|point| self.point(point)))
    }
}

fn controls(curve: &uv::Fit, start: bool, offset: usize) -> Controls {
    std::array::from_fn(|i| {
        let index = if start {
            offset + i
        } else {
            curve.points.len() - 1 - offset - i
        };
        curve.points[index].each_ref().map(|n| n.value.clone())
    })
}

fn split(mut points: Controls, at: &Rational) -> (Controls, Controls) {
    let mut left: Controls = std::array::from_fn(|_| std::array::from_fn(|_| Rational::zero()));
    let mut right = left.clone();
    left[0] = points[0].clone();
    right[fit::DEGREE] = points[fit::DEGREE].clone();
    for level in 1..=fit::DEGREE {
        for i in 0..=fit::DEGREE - level {
            points[i] = std::array::from_fn(|axis| {
                let delta = &points[i + 1][axis] - &points[i][axis];
                &points[i][axis] + &(&delta * at)
            });
        }
        left[level] = points[0].clone();
        right[fit::DEGREE - level] = points[fit::DEGREE - level].clone();
    }
    (left, right)
}

fn compose<const N: usize>(coefficients: &[Rational], variable: &[Rational]) -> [Rational; N] {
    let mut value: [Rational; N] = std::array::from_fn(|_| Rational::zero());
    let mut length = 1;
    for coefficient in coefficients.iter().rev() {
        let mut next: [Rational; N] = std::array::from_fn(|_| Rational::zero());
        for (i, a) in value[..length].iter().enumerate() {
            if a.is_zero() {
                continue;
            }
            for (j, b) in variable.iter().take(N - i).enumerate() {
                if !b.is_zero() {
                    next[i + j] = &next[i + j] + &(a * b);
                }
            }
        }
        next[0] = &next[0] + coefficient;
        value = next;
        length = (length + variable.len() - 1).min(N);
    }
    value
}

fn graph(powers: &[[Rational; 2]; fit::DEGREE + 1]) -> Option<[Rational; GRAPH_DEGREE + 1]> {
    let x = powers.each_ref().map(|p| p[0].clone());
    let y = powers.each_ref().map(|p| p[1].clone());
    if x[1].is_zero() {
        return None;
    }
    let mut inverse: [Rational; GRAPH_DEGREE + 1] = std::array::from_fn(|_| Rational::zero());
    inverse[1] = &Rational::one() / &x[1];
    for order in 2..=GRAPH_DEGREE {
        inverse[order] =
            -compose::<{ GRAPH_DEGREE + 1 }>(&x, &inverse)[order].clone() / x[1].clone();
    }
    Some(compose(&y, &inverse))
}

fn residual(
    controls: &Controls,
    frame: &Frame,
    separator: &[Rational; GRAPH_DEGREE + 1],
) -> [Rational; RESIDUAL_DEGREE + 1] {
    let powers = frame.powers(controls);
    let x = powers.each_ref().map(|p| p[0].clone());
    let mut result = compose::<{ RESIDUAL_DEGREE + 1 }>(separator, &x).map(|v| -v);
    for (value, power) in result.iter_mut().zip(powers) {
        *value = &*value + &power[1];
    }
    result
}

fn bernstein_nonnegative(powers: &[Rational; RESIDUAL_DEGREE + 1], upper: bool) -> bool {
    let mut positive = false;
    for i in 0..=RESIDUAL_DEGREE {
        let mut coefficient = Rational::zero();
        let mut factor = Rational::one();
        for (j, power) in powers.iter().take(i + 1).enumerate() {
            coefficient = coefficient + power * &factor;
            if j != i {
                factor = factor
                    * Rational::from_ratio((i - j) as i64, (RESIDUAL_DEGREE - j) as i64).unwrap();
            }
        }
        if !upper {
            coefficient = -coefficient;
        }
        if coefficient < Rational::zero() {
            return false;
        }
        positive |= !coefficient.is_zero();
    }
    positive
}

fn separates(
    controls: &Controls,
    frame: &Frame,
    separator: &[Rational; GRAPH_DEGREE + 1],
    upper: bool,
    limit: Option<&Rational>,
    depth: u32,
    dec: &Decider,
) -> Result<bool, Refusal> {
    dec.check_deadline()?;
    if limit.is_some_and(|limit| controls.iter().all(|p| frame.point(p)[0] > *limit)) {
        return Ok(true);
    }
    let powers = residual(controls, frame, separator);
    let end = powers
        .iter()
        .fold(Rational::zero(), |sum, value| &sum + value);
    // Bernstein basis functions are strictly positive inside a span. A zero
    // at any subdivision endpoint is permitted only at the common vertex.
    for (point, value) in [(&controls[0], &powers[0]), (&controls[fit::DEGREE], &end)] {
        if value.is_zero() {
            if *point != frame.origin {
                return Ok(false);
            }
        } else if (value > &Rational::zero()) != upper {
            return Ok(false);
        }
    }
    if bernstein_nonnegative(&powers, upper) {
        return Ok(true);
    }
    // An unresolved extra contact cannot consume unbounded exact arithmetic.
    if powers.iter().all(Rational::is_zero) || depth >= frame.max_depth {
        return Ok(false);
    }
    let (left, right) = split(controls.clone(), &Rational::from_ratio(1, 2).unwrap());
    Ok(
        separates(&left, frame, separator, upper, limit, depth + 1, dec)?
            && separates(&right, frame, separator, upper, limit, depth + 1, dec)?,
    )
}

fn corner(
    curve: (&uv::Fit, bool),
    original: &Controls,
    frame: &Frame,
    upper: bool,
    amplitude: &Rational,
    digits: u32,
    dec: &Decider,
) -> Result<Option<Repair>, Refusal> {
    let (curve, start) = curve;
    let index = if start { 0 } else { curve.knots.len() - 2 };
    let [lo, hi] = [&curve.knots[index], &curve.knots[index + 1]];
    let mut fraction = Rational::from_ratio(1, 2).unwrap();
    for _ in 0..dec.policy().max_precision_bits {
        dec.check_deadline()?;
        let knot = number::rational(
            &(if start { lo.clone() } else { hi.clone() }
                + (hi - lo)
                    * if start {
                        fraction.clone()
                    } else {
                        -fraction.clone()
                    }),
            digits,
        )
        .value;
        if knot <= *lo || knot >= *hi {
            return Ok(None);
        }
        let at = if start { &knot - lo } else { hi - &knot } / (hi - lo);
        let (mut left, mut right) = split(original.clone(), &at);
        let inside = left.iter().all(|point| {
            point
                .iter()
                .zip(&frame.origin)
                .fold(Rational::zero(), |sum, (a, b)| sum + (a - b).abs())
                <= *amplitude
        });
        if !inside {
            fraction = fraction.mul_pow2(-1);
            continue;
        }
        let mut displacement = Rational::zero();
        let mut points = Vec::with_capacity(2 * fit::DEGREE + 1);
        for (part, ramp) in [(&mut left, true), (&mut right, false)] {
            for (i, point) in part.iter_mut().enumerate() {
                let mut desired = point.clone();
                let fraction = Rational::from_ratio(
                    if ramp {
                        i as i64
                    } else {
                        (fit::DEGREE - i) as i64
                    },
                    fit::DEGREE as i64,
                )
                .unwrap();
                let shift = amplitude * &fraction;
                let positive_y = upper == (frame.positive == (frame.axis == 0));
                desired[1 - frame.axis] = if positive_y {
                    &desired[1 - frame.axis] + &shift
                } else {
                    &desired[1 - frame.axis] - &shift
                };
                let written = desired.each_ref().map(|v| number::rational(v, digits));
                let error = written
                    .iter()
                    .zip(&*point)
                    .fold(Rational::zero(), |sum, (a, b)| sum + (&a.value - b).abs());
                displacement = displacement.max(error);
                *point = written.each_ref().map(|n| n.value.clone());
                if ramp || i != 0 {
                    points.push(written);
                }
            }
        }
        return Ok(Some(Repair {
            start,
            points,
            knot,
            displacement,
        }));
    }
    Ok(None)
}

/// Inputs are the outgoing lower branch and incoming upper branch of a
/// surface-CCW corner. The coordinate change preserves that orientation.
pub(crate) fn pair(
    lower: (&uv::Fit, bool),
    upper: (&uv::Fit, bool),
    amplitude: &Rational,
    digits: u32,
    dec: &Decider,
) -> Result<Option<[Repair; 2]>, Refusal> {
    dec.check_deadline()?;
    if amplitude <= &Rational::zero() {
        return Ok(None);
    }
    let curves = [lower, upper];
    let original = curves.map(|(curve, start)| controls(curve, start, 0));
    if original[0][0] != original[1][0] {
        return Ok(None);
    }
    let origin = original[0][0].clone();
    let magnitude = |axis: usize| {
        original.iter().fold(Rational::zero(), |sum, c| {
            sum + (&c[1][axis] - &origin[axis]).abs()
        })
    };
    let axis = usize::from(magnitude(1) > magnitude(0));
    let direction = original.iter().fold(Rational::zero(), |sum, c| {
        sum + (&c[fit::DEGREE][axis] - &origin[axis])
    });
    if direction.is_zero() {
        return Ok(None);
    }
    let frame = Frame {
        origin,
        axis,
        positive: direction > Rational::zero(),
        max_depth: dec
            .policy()
            .max_precision_bits
            .min(digits.max(1).saturating_mul(4).saturating_add(32)),
    };
    let ends = original
        .each_ref()
        .map(|c| frame.point(&c[fit::DEGREE])[0].clone());
    for ((curve, start), (points, end)) in curves.iter().zip(original.iter().zip(&ends)) {
        if end <= &Rational::zero()
            || points
                .windows(2)
                .any(|p| frame.point(&p[0])[0] > frame.point(&p[1])[0])
            || (fit::DEGREE..curve.points.len()).any(|i| {
                let index = if *start {
                    i
                } else {
                    curve.points.len() - 1 - i
                };
                frame.point(&curve.points[index].each_ref().map(|n| n.value.clone()))[0] < *end
            })
        {
            return Ok(None);
        }
    }
    let Some(a) = graph(&frame.powers(&original[0])) else {
        return Ok(None);
    };
    let Some(b) = graph(&frame.powers(&original[1])) else {
        return Ok(None);
    };
    let separator = std::array::from_fn(|i| (&a[i] + &b[i]).mul_pow2(-1));
    let Some(a) = corner(lower, &original[0], &frame, false, amplitude, digits, dec)? else {
        return Ok(None);
    };
    let Some(b) = corner(upper, &original[1], &frame, true, amplitude, digits, dec)? else {
        return Ok(None);
    };
    let repairs = [a, b];
    for (index, repair) in repairs.iter().enumerate() {
        let upper = index != 0;
        for offset in [0, fit::DEGREE] {
            let points = std::array::from_fn(|i| {
                repair.points[offset + i]
                    .each_ref()
                    .map(|n| n.value.clone())
            });
            if !separates(&points, &frame, &separator, upper, None, 0, dec)? {
                return Ok(None);
            }
        }
        let (curve, start) = curves[index];
        for offset in (fit::DEGREE..curve.points.len() - 1).step_by(fit::DEGREE) {
            if !separates(
                &controls(curve, start, offset),
                &frame,
                &separator,
                upper,
                Some(&ends[1 - index]),
                0,
                dec,
            )? {
                return Ok(None);
            }
        }
    }
    Ok(Some(repairs))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn graph_separator_allows_only_the_named_corner() {
        // Degree-seven Bernstein controls for t^4 and t^2(t-1/2)^2.
        // The second graph touches the separator again at x=1/2.
        let cases = [
            (
                true,
                [
                    (0, 1),
                    (0, 1),
                    (0, 1),
                    (0, 1),
                    (1, 35),
                    (1, 7),
                    (3, 7),
                    (1, 1),
                ],
            ),
            (
                false,
                [
                    (0, 1),
                    (0, 1),
                    (1, 84),
                    (1, 140),
                    (-1, 70),
                    (-1, 42),
                    (1, 28),
                    (1, 4),
                ],
            ),
            // t^2(t^2-1/2)^2 also touches at the irrational x=1/sqrt(2).
            (
                false,
                [
                    (0, 1),
                    (0, 1),
                    (1, 84),
                    (1, 28),
                    (3, 70),
                    (-1, 42),
                    (-3, 28),
                    (1, 4),
                ],
            ),
        ];
        let frame = Frame {
            origin: [Rational::zero(), Rational::zero()],
            axis: 0,
            positive: true,
            max_depth: 32,
        };
        let separator = std::array::from_fn(|_| Rational::zero());
        for (expected, y) in cases {
            let points = std::array::from_fn(|i| {
                [
                    Rational::from_ratio(i as i64, fit::DEGREE as i64).unwrap(),
                    Rational::from_ratio(y[i].0, y[i].1).unwrap(),
                ]
            });
            assert_eq!(
                separates(
                    &points,
                    &frame,
                    &separator,
                    true,
                    None,
                    0,
                    &Decider::default()
                )
                .unwrap(),
                expected,
            );
        }
    }
}
