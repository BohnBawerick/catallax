//! Certified readouts for a pair of finite straight edges, from exact endpoints.

use catallax_numbers::{Answer, Answered, Decider, MarkLog, Rational, Real, Refusal, Sign};
use catallax_shapes::{Curve, RPoint3};
use catallax_topology::{EdgeId, Model};

use crate::{Accuracy, Enclosure};

/// The scalar requested from two finite straight edges.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum EdgePairRequest {
    /// Minimum distance over both closed segments, in model units.
    MinimumDistance,
    /// Distance between both half-arc-length points, in model units.
    MidpointDistance,
    /// Unoriented acute carrier angle in degrees. Only proven parallel (0) and
    /// perpendicular (90) carriers are answered.
    DirectionAngleDegrees,
}

pub(crate) fn measure(
    model: &Model,
    a: EdgeId,
    b: EdgeId,
    request: EdgePairRequest,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<Enclosure> {
    accuracy.validate();
    dec.check_deadline()?;
    let p = segment(model, a, "first")?;
    let q = segment(model, b, "second")?;
    match request {
        EdgePairRequest::MinimumDistance => distance(minimum_square(&p, &q), accuracy, dec),
        EdgePairRequest::MidpointDistance => {
            let m = midpoint(model, a, accuracy, dec)?;
            let n = midpoint(model, b, accuracy, dec)?;
            distance(m.minus(&n).norm_squared(), accuracy, dec)
        }
        EdgePairRequest::DirectionAngleDegrees => {
            let (u, v) = (p[1].minus(&p[0]), q[1].minus(&q[0]));
            let degrees = if u.cross(&v).is_zero() {
                0
            } else if u.dot(&v).is_zero() {
                90
            } else {
                return Err(Refusal::new(
                    "measure.angle-domain",
                    format!("edges {a:?} and {b:?} are neither parallel nor perpendicular"),
                ));
            };
            Ok(Answer::exact(Enclosure::exact(Rational::from_i64(degrees))))
        }
    }
}

/// Exact rational ends of a nondegenerate finite line edge.
fn segment(model: &Model, edge: EdgeId, role: &str) -> Result<[RPoint3; 2], Refusal> {
    let e = model.edge(edge);
    let (Curve::Line(_), Some(start), Some(end)) = (model.geom().curve(e.curve), e.start, e.end)
    else {
        return Err(Refusal::new(
            "measure.unsupported-curve",
            format!("{role} edge {edge:?} is not a finite straight edge"),
        ));
    };
    let ends = [start, end].map(|v| model.geom().position(model.vertex(v).point).as_rational());
    let [Some(start), Some(end)] = ends else {
        return Err(Refusal::new(
            "measure.irrational-endpoint",
            format!("{role} edge {edge:?} has an endpoint without exact rational coordinates"),
        ));
    };
    if start == end {
        return Err(Refusal::new(
            "measure.degenerate-edge",
            format!("{role} edge {edge:?} has zero length"),
        ));
    }
    Ok([start, end])
}

/// The native half-arc-length point. Rational line ends make it exact.
fn midpoint(
    model: &Model,
    edge: EdgeId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Result<RPoint3, Refusal> {
    let answer = crate::edge::midpoint(model, edge, accuracy, dec)?;
    match answer.value.map(|c| c.is_exact().then_some(c.lo)) {
        [Some(x), Some(y), Some(z)] => Ok(RPoint3::new(x, y, z)),
        _ => Err(Refusal::new(
            "measure.irrational-endpoint",
            format!("edge {edge:?} has no exact rational midpoint"),
        )),
    }
}

/// Squared minimum of `|p(s) - q(t)|` over `s, t` in `[0, 1]`. The objective is a
/// convex quadratic, so its minimum is an interior critical point or the clamped
/// minimum along one side of the parameter square.
fn minimum_square(p: &[RPoint3; 2], q: &[RPoint3; 2]) -> Rational {
    let u = p[1].minus(&p[0]);
    let v = q[1].minus(&q[0]);
    let w = p[0].minus(&q[0]);
    let (a, b, c) = (u.dot(&u), u.dot(&v), v.dot(&v));
    let (d, e) = (u.dot(&w), v.dot(&w));
    let unit = |x: Rational| x.max(Rational::zero()).min(Rational::one());
    let at = |s: &Rational, t: &Rational| w.plus(&u.scaled(s)).minus(&v.scaled(t)).norm_squared();
    let mut candidates = Vec::with_capacity(5);
    for side in [Rational::zero(), Rational::one()] {
        let t = unit(&(&(&b * &side) + &e) / &c);
        candidates.push(at(&side, &t));
        let s = unit(&(&(&b * &side) - &d) / &a);
        candidates.push(at(&s, &side));
    }
    let det = &(&a * &c) - &(&b * &b);
    if det.sign() == Sign::Positive {
        let s = &(&(&b * &e) - &(&c * &d)) / &det;
        let t = &(&(&a * &e) - &(&b * &d)) / &det;
        if s == unit(s.clone()) && t == unit(t.clone()) {
            candidates.push(at(&s, &t));
        }
    }
    candidates.into_iter().min().unwrap()
}

fn distance(square: Rational, accuracy: &Accuracy, dec: &mut Decider) -> Answered<Enclosure> {
    let enclosure = root(square, dec)?;
    let mut marks = MarkLog::new();
    accuracy.mark(&mut marks, &[&enclosure], "edge pair distance");
    Ok(Answer::with_marks(enclosure, marks))
}

/// The exact root of a nonnegative rational when it is rational, otherwise an
/// outward enclosure of nonzero width.
pub(crate) fn root(square: Rational, dec: &mut Decider) -> Result<Enclosure, Refusal> {
    if square.is_zero() {
        return Ok(Enclosure::exact(square));
    }
    let root = Real::from(square).sqrt();
    if let Some(value) = dec.rational_value(&root)? {
        return Ok(Enclosure::exact(value));
    }
    let bits = dec.policy().max_precision_bits.min(256);
    Ok(Enclosure::from_ball(&dec.enclose_budgeted_at(&root, bits)?))
}

#[cfg(test)]
mod tests;
