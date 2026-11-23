//! Enclosures of a proved cap cut, without replacing it by its STEP fit.

use super::*;
use crate::tessellation::{RationalSpan, rational_span_range};
use catallax_numbers::Interval;

/// A proved plane cut, parametrised by the carrier's second parameter.
/// Its private coefficients retain the monotonicity proof of `plane_cut`.
#[derive(Clone, Debug)]
pub struct FilletCutTrace {
    fillet: SplineFillet,
    cut: FilletCut,
    coefficients: Vec<[Vec<Rational>; 3]>,
}

/// Binary64 enclosures of the same cap equation and carrier used by the trace.
pub struct DisplayCut {
    spans: Vec<DisplaySpan>,
    before: Sign,
}

/// A point enclosure and its parameter enclosure on the cut's own carrier.
pub struct DisplayCutSample {
    /// Cartesian coordinate bounds in millimetres.
    pub position: [Interval; 3],
    /// Axial and across-fillet parameters, in carrier coordinates.
    pub parameters: [Interval; 2],
}

struct DisplaySpan {
    g: [Vec<Interval>; 3],
    speed: Interval,
    carrier: RationalSpan,
}

fn evaluate(mut row: Vec<Ball>, t: &Ball) -> Ball {
    for level in 1..row.len() {
        for i in 0..row.len() - level {
            row[i] = &row[i] + &(t * &(&row[i + 1] - &row[i]));
        }
    }
    row[0].clone()
}

fn difference(row: &[Interval]) -> Vec<Interval> {
    let degree = Interval::point(row.len().saturating_sub(1) as f64).unwrap();
    row.windows(2)
        .map(|pair| degree * (pair[1] - pair[0]))
        .collect()
}

fn hull(values: impl IntoIterator<Item = Interval>) -> Interval {
    values
        .into_iter()
        .reduce(|a, b| a.hull(&b))
        .unwrap_or_else(|| Interval::point(0.0).unwrap())
}

// Woźny and Chudy (2020), Algorithm 1.1. Reverse below one half
// so the denominator in the parameter ratio is at least one half.
fn display_bernstein(values: &[Interval], at: f64, dec: &mut Decider) -> Result<Interval, Refusal> {
    if at <= 0.0 {
        return Ok(values[0]);
    }
    if at >= 1.0 {
        return Ok(*values.last().unwrap());
    }
    let one = Interval::point(1.0).unwrap();
    let parameter = Interval::point(at).unwrap();
    let reversed = at < 0.5;
    let ratio = if reversed {
        parameter / (one - parameter)
    } else {
        (one - parameter) / parameter
    };
    let degree = values.len() - 1;
    let control = |i| values[if reversed { degree - i } else { i }];
    let mut weight = one;
    let mut value = control(0);
    dec.charge_work(3 * degree as u64)?;
    for i in 1..=degree {
        let order = Interval::point(i as f64).unwrap();
        let remaining = Interval::point((degree - i + 1) as f64).unwrap();
        weight = one / (one + (order / remaining) * ratio / weight);
        value = (one - weight) * value + weight * control(i);
    }
    Ok(value)
}

impl SplineFillet {
    /// Proves a plane cut and prepares its exact equation for point queries.
    pub fn cut_trace(
        &self,
        normal: [Rational; 3],
        offset: Rational,
        dec: &mut Decider,
    ) -> Result<Option<FilletCutTrace>, Refusal> {
        let Some(cut) = self.plane_cut(normal, offset, dec)? else {
            return Ok(None);
        };
        let (normal, offset) = surface::local_plane(&self.data.frame, &cut.normal, &cut.offset)
            .expect("plane_cut checked the rational frame");
        let coefficients = (cut.spans[0]..=cut.spans[1])
            .map(|span| cut::coefficients(&self.data.fit.rows[span], &normal, &offset, dec))
            .collect::<Result<_, _>>()?;
        Ok(Some(FilletCutTrace {
            fillet: self.clone(),
            cut,
            coefficients,
        }))
    }
}

impl FilletCutTrace {
    /// The checked cut, including its exact rail endpoints.
    pub fn cut(&self) -> &FilletCut {
        &self.cut
    }

    /// The fillet whose carrier supplies the cut's parameters.
    pub fn fillet(&self) -> &SplineFillet {
        &self.fillet
    }

    /// The exact endpoint on rail `side`.
    pub fn endpoint(&self, side: usize) -> Point3 {
        self.fillet
            .rail_point(side, &self.cut.ends[side])
            .expect("a cut ends at rail nodes")
    }

    /// Prepares outward-rounded coefficients for repeated display sampling.
    pub fn display_chart(&self, dec: &mut Decider) -> Result<DisplayCut, Refusal> {
        let [first, last] = self.cut.spans;
        let spans = rational_span_range(self.fillet.carrier(), first..=last, dec)?
            .into_iter()
            .zip(&self.coefficients)
            .map(|(carrier, rows)| {
                let g = rows
                    .each_ref()
                    .map(|row| row.iter().map(Interval::from_rational).collect::<Vec<_>>());
                let speed = hull(g.iter().flat_map(|row| difference(row)));
                dec.charge_work(g.iter().map(|row| row.len() - 1).sum::<usize>() as u64)?;
                Ok(DisplaySpan { g, speed, carrier })
            })
            .collect::<Result<Vec<_>, Refusal>>()?;
        Ok(DisplayCut {
            spans,
            before: self.cut.before,
        })
    }

    /// Encloses every cut point for `theta` in `[0,1]`. Bisection only
    /// tightens a box; an undecided sign keeps it, and never discards a root.
    pub fn enclose(&self, theta: &Ball, dec: &Decider) -> Result<[Ball; 3], Refusal> {
        let bits = theta.precision().min(dec.policy().max_precision_bits);
        let zero = Rational::zero();
        let one = Rational::one();
        if bits == 0 || theta.lo() < zero || theta.hi() > one {
            return Err(Refusal::new(
                SPLINE_FILLET_CUT,
                "cut parameter is outside [0,1] or has no precision",
            ));
        }
        if theta.lo() == theta.hi() && (theta.lo() == zero || theta.lo() == one) {
            let endpoint = self.endpoint(usize::from(theta.lo() == one));
            let [x, y, z] =
                [&endpoint.x, &endpoint.y, &endpoint.z].map(|v| dec.enclose_at(v, bits));
            return Ok([x?, y?, z?]);
        }
        let fit = &self.fillet.data.fit;
        let [first, last] = self.cut.spans;
        let ball = |q: &Rational| Ball::from_rational(q, bits);
        let theta = theta.with_precision(bits);
        let side = |tau: &Rational| -> Result<Option<Sign>, Refusal> {
            dec.check_deadline()?;
            let span = (first..=last)
                .find(|&i| tau <= &fit.nodes[i + 1].tau)
                .unwrap_or(last);
            let start = &fit.nodes[span].tau;
            let width = &fit.nodes[span + 1].tau - start;
            let s = ball(&((tau - start) / width));
            let rows = self.coefficients[span - first]
                .each_ref()
                .map(|row| evaluate(row.iter().map(&ball).collect(), &s));
            let value = evaluate(rows.to_vec(), &theta);
            Ok(if value.lo() > zero {
                Some(Sign::Positive)
            } else if value.hi() < zero {
                Some(Sign::Negative)
            } else {
                None
            })
        };
        let range = [
            fit.nodes[first].tau.clone(),
            fit.nodes[last + 1].tau.clone(),
        ];
        // The lower bracket stays before every root; the upper stays after
        // every root, even when theta is an interval with many root values.
        let mut lower = range.clone();
        let mut upper = range;
        for _ in 0..24 {
            let middle = (&lower[0] + &lower[1]).mul_pow2(-1);
            if side(&middle)? == Some(self.cut.before) {
                lower[0] = middle;
            } else {
                lower[1] = middle;
            }
            let middle = (&upper[0] + &upper[1]).mul_pow2(-1);
            if side(&middle)?.is_some_and(|sign| sign != self.cut.before) {
                upper[1] = middle;
            } else {
                upper[0] = middle;
            }
        }
        let tau = ball(&lower[0]).hull(&ball(&upper[1]));
        self.fillet.carrier().evaluate(&tau, &theta, dec)
    }
}

impl DisplayCut {
    /// Bounds the spatial acceleration, in millimetres per squared theta.
    /// The carrier has shared first jets across its internal knots.
    pub fn acceleration_bound(&self, dec: &mut Decider) -> Result<Rational, Refusal> {
        let twice = Interval::point(2.0).unwrap();
        let magnitude =
            |bound: Interval| Interval::point(bound.lo().abs().max(bound.hi().abs())).unwrap();
        let mut maximum = 0.0_f64;
        for span in &self.spans {
            let g = &span.g;
            let degree = g[0].len() - 1;
            if span.speed.contains_zero() {
                return Err(Refusal::new(
                    SPLINE_FILLET_CUT,
                    "binary64 could not separate the plane cut's monotone derivative",
                ));
            }
            let gs = g.each_ref().map(|row| difference(row));
            let gt = std::array::from_fn::<_, 2, _>(|i| {
                g[i + 1]
                    .iter()
                    .zip(&g[i])
                    .map(|(&b, &a)| twice * (b - a))
                    .collect::<Vec<_>>()
            });
            dec.charge_work(
                (3 * degree + 2 * (degree + 1) + 3 * (degree - 1) + 2 * degree) as u64,
            )?;
            let gss = hull(gs.iter().flat_map(|row| difference(row)));
            let gst = hull(gt.iter().flat_map(|row| difference(row)));
            let gtt = hull(gt[1].iter().zip(&gt[0]).map(|(&b, &a)| b - a));
            // Differentiate the monotone Bernstein equation g(s(theta), theta) = 0.
            // The derivatives and the carrier bounds both use the local span coordinates.
            let rate = -hull(gt.into_iter().flatten()) / span.speed;
            let rate_sq = rate.square();
            let acceleration = -(gtt + twice * gst * rate + gss * rate_sq) / span.speed;
            let [ss, st, tt] = span.carrier.affine_hessian();
            let [speed, _] = span.carrier.first_derivative_bounds();
            dec.charge_work(1 + 3 + 4)?;
            let bound =
                ss * rate_sq + twice * st * magnitude(rate) + tt + speed * magnitude(acceleration);
            maximum = maximum.max(bound.hi());
        }
        Rational::from_f64_exact(maximum).ok_or_else(|| {
            Refusal::new(SPLINE_FILLET_CUT, "plane cut acceleration exceeds binary64")
        })
    }

    /// Encloses the cut point by interval Newton inside the Bernstein spans.
    /// A sample that does not converge in binary64 uses the original trace.
    pub fn enclose(
        &self,
        theta: Interval,
        output_width: f64,
        dec: &mut Decider,
    ) -> Result<Option<DisplayCutSample>, Refusal> {
        let blend = |a: Interval, b: Interval| a + (b - a) * theta;
        let side = |g: Interval, sign: Sign| match sign {
            Sign::Positive => g.lo() > 0.0,
            Sign::Negative => g.hi() < 0.0,
            Sign::Zero => false,
        };
        let after = if self.before == Sign::Positive {
            Sign::Negative
        } else {
            Sign::Positive
        };
        let mut found: Option<DisplayCutSample> = None;
        for span in &self.spans {
            let [a, b, c] = &span.g;
            let column = |i: usize| blend(blend(a[i], b[i]), blend(b[i], c[i]));
            let end = a.len() - 1;
            dec.charge_work(6)?;
            if side(column(0), after) || side(column(end), self.before) {
                continue;
            }
            if span.speed.lo() <= 0.0 && span.speed.hi() >= 0.0 {
                return Ok(None);
            }
            let values: Vec<_> = (0..=end).map(column).collect();
            dec.charge_work((3 * values.len()) as u64)?;
            let mut root = Interval::new(0.0, 1.0).unwrap();
            let speed = span.carrier.first_derivative_bounds()[0].hi();
            for _ in 0..f64::MANTISSA_DIGITS {
                // The Euclidean speed bounds the L1 output width after
                // multiplication by sqrt(3), bounded here by two.
                if root.width() * speed <= output_width * 0.5 {
                    break;
                }
                let middle = Interval::point(root.mid()).unwrap();
                let value = display_bernstein(&values, middle.mid(), dec)?;
                let next = middle - value / span.speed;
                let Some(next) = Interval::new(root.lo().max(next.lo()), root.hi().min(next.hi()))
                else {
                    return Ok(None);
                };
                if next.width() >= root.width() {
                    return Ok(None);
                }
                root = next;
            }
            let position = span.carrier.evaluate_box([root, theta], dec)?;
            let [low, high] = &span.carrier.domain()[0];
            let origin = Interval::from_rational(low);
            dec.charge_work(1)?;
            let parameter = origin + Interval::from_rational(&(high - low)) * root;
            let (position, parameter) = found.map_or((position, parameter), |cover| {
                (
                    std::array::from_fn(|i| cover.position[i].hull(&position[i])),
                    cover.parameters[0].hull(&parameter),
                )
            });
            found = Some(DisplayCutSample {
                position,
                parameters: [parameter, theta],
            });
        }
        Ok(found)
    }
}

#[cfg(test)]
mod display_tests {
    use super::*;

    #[test]
    fn scalar_bernstein_encloses_the_cubic_from_either_end() {
        let controls = [0.0, 0.0, 0.0, 1.0].map(|value| Interval::point(value).unwrap());
        let mut dec = Decider::default();
        for at in [0.0_f64, 0.25, 0.5, 0.75, 1.0] {
            let value = display_bernstein(&controls, at, &mut dec).unwrap();
            assert!(value.contains(at.powi(3)));
            assert!(value.width() < f64::EPSILON.sqrt());
        }
    }
}
