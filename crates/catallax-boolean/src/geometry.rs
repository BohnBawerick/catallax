//! Certified helpers shared by the stages and point queries.

use catallax_intersect::{Cyclic, curve_is_closed, cyclic_order, on_curve, order_on_curve};
use catallax_numbers::{
    Decider, Fact, FactKind, FactOrigin, Outcome, Rational, Real, Refusal, Sign, StrictOrder,
    StrictSign,
};
use catallax_shapes::{CurveId, GeomKey, PointId, RVector3, Vector3};
use catallax_topology::{EdgeId, Model, Sense};

// Geometric search counts are independent of arithmetic working precision.
// Exhausting these searches refuses; it never turns an unresolved sign into
// a geometric answer.
pub(crate) const DIRECTION_ATTEMPTS: u32 = 64;
pub(crate) const PROBE_DEPTHS: u32 = 256;

pub(crate) fn refusal(code: &'static str, detail: impl Into<String>) -> Refusal {
    Refusal::new(code, detail)
}

pub(crate) fn sign(x: &Real, dec: &mut Decider, code: &'static str) -> Result<Sign, Refusal> {
    match dec.sign(x) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(refusal(code, n.detail)),
        Outcome::Refused(r) => Err(r),
    }
}

pub(crate) fn coincident<T>(outcome: Outcome<T>, code: &'static str) -> Result<bool, Refusal> {
    match outcome {
        Outcome::Proven(_) => Ok(false),
        Outcome::ProvenCoincident(_) => Ok(true),
        Outcome::NotExact(n) => Err(refusal(code, n.detail)),
        Outcome::Refused(r) => Err(r),
    }
}

pub(crate) fn reverse(sense: Sense) -> Sense {
    match sense {
        Sense::Same => Sense::Reversed,
        Sense::Reversed => Sense::Same,
    }
}

pub(crate) fn product(a: Sense, b: Sense) -> Sense {
    if a == b { Sense::Same } else { Sense::Reversed }
}

pub(crate) fn oriented(v: &Vector3, sense: Sense) -> Vector3 {
    if sense == Sense::Same {
        v.clone()
    } else {
        v.scaled(&Real::from_i64(-1))
    }
}

pub(crate) fn natural_ends(model: &Model, edge: EdgeId) -> (Option<PointId>, Option<PointId>) {
    let edge = model.edge(edge);
    let a = edge.start.map(|v| model.vertex(v).point);
    let b = edge.end.map(|v| model.vertex(v).point);
    if edge.sense == Sense::Same {
        (a, b)
    } else {
        (b, a)
    }
}

pub(crate) fn on_arc(
    model: &Model,
    curve: CurveId,
    from: Option<PointId>,
    to: Option<PointId>,
    point: PointId,
    dec: &mut Decider,
    code: &'static str,
) -> Result<bool, Refusal> {
    if !coincident(on_curve(model.geom(), curve, point, dec), code)? {
        return Ok(false);
    }
    let (Some(a), Some(b)) = (from, to) else {
        assert!(from.is_none() && to.is_none());
        return Ok(true);
    };
    if coincident(model.geom().same_point(a, point, dec), code)?
        || coincident(model.geom().same_point(b, point, dec), code)?
    {
        return Ok(true);
    }
    if coincident(model.geom().same_point(a, b, dec), code)? {
        return Ok(true);
    }
    if !curve_is_closed(model.geom(), curve, dec)? {
        let after = match order_on_curve(model.geom(), curve, a, point, dec) {
            Outcome::Proven(StrictOrder::Less) => true,
            Outcome::Proven(_) | Outcome::ProvenCoincident(_) => false,
            Outcome::NotExact(n) => return Err(refusal(code, n.detail)),
            Outcome::Refused(r) => return Err(r),
        };
        let before = match order_on_curve(model.geom(), curve, point, b, dec) {
            Outcome::Proven(StrictOrder::Less) => true,
            Outcome::Proven(_) | Outcome::ProvenCoincident(_) => false,
            Outcome::NotExact(n) => return Err(refusal(code, n.detail)),
            Outcome::Refused(r) => return Err(r),
        };
        Ok(after && before)
    } else {
        match cyclic_order(model.geom(), curve, [a, point, b], dec) {
            Outcome::Proven(Cyclic::Forward) => Ok(true),
            Outcome::Proven(Cyclic::Backward) => Ok(false),
            Outcome::ProvenCoincident(_) => Ok(true),
            Outcome::NotExact(n) => Err(refusal(code, n.detail)),
            Outcome::Refused(r) => Err(r),
        }
    }
}

pub(crate) fn record_incident(model: &mut Model, a: GeomKey, b: GeomKey) {
    if !model.geom().facts().holds(FactKind::Incident, a, b) {
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a,
            b,
            at: None,
            origin: FactOrigin::Construction("boolean.sew"),
        });
    }
}

pub(crate) fn record_same_point(model: &mut Model, a: PointId, b: PointId) {
    if a != b {
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Coincident,
            a: GeomKey::Point(a),
            b: GeomKey::Point(b),
            at: None,
            origin: FactOrigin::Proof(catallax_numbers::ZeroProof::Algebraic),
        });
    }
}

pub(crate) fn canonical_point(
    model: &mut Model,
    points: &mut Vec<PointId>,
    point: PointId,
    dec: &mut Decider,
    code: &'static str,
) -> Result<PointId, Refusal> {
    for &existing in points.iter() {
        if coincident(model.geom().same_point(existing, point, dec), code)? {
            record_same_point(model, existing, point);
            return Ok(existing);
        }
    }
    points.push(point);
    Ok(point)
}

pub(crate) fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
pub(crate) fn unit(axis: usize) -> Vector3 {
    RVector3::new(
        r(i64::from(axis == 0)),
        r(i64::from(axis == 1)),
        r(i64::from(axis == 2)),
    )
    .to_real()
}

pub(crate) fn perpendicular(
    vector: &Vector3,
    dec: &mut Decider,
    code: &'static str,
) -> Result<Vector3, Refusal> {
    for axis in 0..3 {
        let direction = vector.cross(&unit(axis));
        if sign(&direction.norm_squared(), dec, code)? != Sign::Zero {
            return Ok(direction);
        }
    }
    Err(refusal(code, "zero direction has no transverse axis"))
}

/// A diagonal enumeration of two signed Calkin-Wilf sequences. Its rational
/// pairs are dense, unlike a single polynomial family of ray directions.
pub(crate) fn probe_parameters(index: u32) -> (Rational, Rational) {
    fn fraction(index: u64) -> Rational {
        if index == 0 {
            return Rational::zero();
        }
        let positive = index.div_ceil(2);
        let mut numerator = 1i64;
        let mut denominator = 1i64;
        let highest = 63 - positive.leading_zeros();
        for bit in (0..highest).rev() {
            if positive & (1u64 << bit) == 0 {
                denominator += numerator;
            } else {
                numerator += denominator;
            }
        }
        if index.is_multiple_of(2) {
            numerator = -numerator;
        }
        Rational::from_ratio(numerator, denominator).expect("positive Calkin-Wilf denominator")
    }
    let index = u64::from(index);
    let diagonal = ((index * 8 + 1).isqrt() - 1) / 2;
    let left = index - diagonal * (diagonal + 1) / 2;
    (fraction(left), fraction(diagonal - left))
}
