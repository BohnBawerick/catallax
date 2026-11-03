//! A box round a stretch of one spline span, from its rational Bezier controls.
//!
//! Interval arithmetic on a span's polynomial in the knot parameter loses the
//! cancellation between its large coefficients. On a degree twelve rail it read a
//! stretch of a cylinder as thousands of millimetres tall, so a trace beside a
//! query ran out of cells before any box cleared it. Controls cut down to the
//! stretch hold it by the convex hull property, so the box shrinks with the curve.

use catallax_numbers::{Ball, Decider, Rational, Real, Refusal};
use catallax_shapes::BSplineCurve;
use std::sync::{Arc, Mutex};

type Homogeneous = [Ball; 4];
type Controls = Arc<Vec<Homogeneous>>;

#[derive(Debug)]
pub(crate) struct SpanHull {
    curve: BSplineCurve,
    span: usize,
    range: [Real; 2],
    controls: Mutex<Option<(u32, Controls)>>,
}

// A mix of two controls.
fn between(a: &Homogeneous, c: &Homogeneous, at: &Ball, rest: &Ball) -> Homogeneous {
    std::array::from_fn(|k| &(&a[k] * rest) + &(&c[k] * at))
}

// The controls of the stretch `[low, high]` of a Bezier segment. Both ends may be
// balls, so the answer holds for every stretch between them.
fn cut(mut controls: Vec<Homogeneous>, low: &Ball, high: &Ball) -> Option<Vec<Homogeneous>> {
    let degree = controls.len() - 1;
    let one = Ball::from_rational(&Rational::one(), high.precision());
    let rest = &one - high;
    for round in 1..=degree {
        for j in (round..=degree).rev() {
            controls[j] = between(&controls[j - 1], &controls[j], high, &rest);
        }
    }
    let share = low.checked_div(high)?;
    let rest = &one - &share;
    for round in 1..=degree {
        for j in 0..=degree - round {
            controls[j] = between(&controls[j], &controls[j + 1], &share, &rest);
        }
    }
    Some(controls)
}

impl SpanHull {
    pub(crate) fn new(curve: &BSplineCurve, span: usize, range: [Real; 2]) -> Self {
        Self {
            curve: curve.clone(),
            span,
            range,
            controls: Mutex::new(None),
        }
    }

    // The span's own Bezier controls, as the blossom at its two knots.
    fn span_controls(&self, bits: u32, dec: &Decider) -> Result<Vec<Homogeneous>, Refusal> {
        let (curve, span) = (&self.curve, self.span);
        let degree = curve.degree();
        let knots = curve.knots();
        let mut row = Vec::with_capacity(degree + 1);
        for i in span - degree..=span {
            dec.check_deadline()?;
            let (point, weight) = (curve.control(i), curve.weight(i));
            let scaled = [&point.x, &point.y, &point.z].map(|value| value * weight);
            let [x, y, z] = scaled.each_ref().map(|value| dec.enclose_at(value, bits));
            row.push([x?, y?, z?, dec.enclose_at(weight, bits)?]);
        }
        let mut controls = Vec::with_capacity(degree + 1);
        for index in 0..=degree {
            dec.check_deadline()?;
            let mut values = row.clone();
            for round in 1..=degree {
                let at = &knots[span + usize::from(round > degree - index)];
                for j in (round..=degree).rev() {
                    let from = span - degree + j;
                    let low = &knots[from];
                    // A nonzero span lies inside each of these supports.
                    let share = (at - low) / (&knots[from + degree + 1 - round] - low);
                    let rest = Rational::one() - share.clone();
                    let [share, rest] = [&share, &rest].map(|v| Ball::from_rational(v, bits));
                    values[j] = between(&values[j - 1], &values[j], &share, &rest);
                }
            }
            controls.push(values.pop().expect("a span has controls"));
        }
        Ok(controls)
    }

    // The piece's controls at a precision, cached for its latest one.
    fn controls(&self, bits: u32, dec: &Decider) -> Result<Option<Controls>, Refusal> {
        let mut cache = self.controls.lock().expect("no panic holds this lock");
        if let Some((held, controls)) = cache.as_ref()
            && *held == bits
        {
            return Ok(Some(controls.clone()));
        }
        let knots = self.curve.knots();
        let start = Real::from_rational(knots[self.span].clone());
        let width = Real::from_rational(&knots[self.span + 1] - &knots[self.span]);
        let [low, high] = self
            .range
            .each_ref()
            .map(|end| dec.enclose_at(&(&(end - &start) / &width), bits));
        let Some(controls) = cut(self.span_controls(bits, dec)?, &low?, &high?) else {
            return Ok(None);
        };
        let controls = Arc::new(controls);
        *cache = Some((bits, controls.clone()));
        Ok(Some(controls))
    }

    /// A box round the piece over `t`, or `None` if the controls' weights do not
    /// separate from zero at this precision.
    pub(crate) fn enclose(&self, t: &Ball, dec: &Decider) -> Result<Option<[Ball; 3]>, Refusal> {
        dec.check_deadline()?;
        let bits = t.precision();
        let Some(controls) = self.controls(bits, dec)? else {
            return Ok(None);
        };
        let [low, high] = [t.lo(), t.hi()].map(|end| Ball::from_rational(&end, bits));
        let Some(stretch) = cut(controls.to_vec(), &low, &high) else {
            return Ok(None);
        };
        let mut found: Option<[Ball; 3]> = None;
        for [x, y, z, weight] in &stretch {
            if weight.lo() <= Rational::zero() {
                return Ok(None);
            }
            let [Some(x), Some(y), Some(z)] = [x, y, z].map(|value| value.checked_div(weight))
            else {
                return Ok(None);
            };
            found = Some(match found {
                None => [x, y, z],
                Some([a, c, e]) => [a.hull(&x), c.hull(&y), e.hull(&z)],
            });
        }
        Ok(found)
    }
}
