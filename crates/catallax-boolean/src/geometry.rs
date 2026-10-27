//! Certified helpers shared by the stages and point queries.

use catallax_intersect::{Cyclic, curve_is_closed, cyclic_order, on_curve, order_on_curve};
use catallax_numbers::{
    Decider, Fact, FactKind, FactOrigin, Outcome, Rational, Real, RealPoly, Refusal, Sign,
    StrictOrder, StrictSign,
};
use catallax_shapes::{Curve, CurveId, GeomKey, PointId, RVector3, Vector3};
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

/// Whether an existing exact relation proves this equation zero at the source.
/// This says nothing about a different selected component or the rest of a path.
pub(crate) fn equation_is_source_contact(
    model: &Model,
    source: PointId,
    source_curve: CurveId,
    barrier: CurveId,
    equation_index: usize,
    equation: &RealPoly,
    source_equations: &[RealPoly],
) -> bool {
    let facts = model.geom().facts();
    if facts.holds(
        FactKind::Coincident,
        GeomKey::Curve(source_curve),
        GeomKey::Curve(barrier),
    ) || facts.holds(
        FactKind::Incident,
        GeomKey::Point(source),
        GeomKey::Curve(barrier),
    ) || source_equations.iter().any(|e| e.same_definition(equation))
    {
        return true;
    }
    if let Curve::Intersection(branch) = model.geom().curve(barrier) {
        let surface = GeomKey::Surface(branch.surfaces[equation_index]);
        if facts.holds(FactKind::Incident, GeomKey::Point(source), surface)
            || facts.holds(FactKind::Incident, GeomKey::Curve(source_curve), surface)
        {
            return true;
        }
    }
    let opposite = -equation;
    source_equations
        .iter()
        .any(|e| e.same_definition(&opposite))
}

pub(crate) fn canonical_point(
    model: &mut Model,
    points: &mut Vec<PointId>,
    point: PointId,
    dec: &mut Decider,
    code: &'static str,
) -> Result<PointId, Refusal> {
    for &existing in points.iter() {
        let identity = if model.geom().facts().holds(
            FactKind::Coincident,
            GeomKey::Point(existing),
            GeomKey::Point(point),
        ) {
            None
        } else {
            crate::point_identity::same_cylinder_point(model.geom(), existing, point, dec)?
        };
        let equal = match identity {
            Some(equal) => equal,
            None => coincident(model.geom().same_point(existing, point, dec), code)?,
        };
        if equal {
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

#[cfg(test)]
mod contact_fact_tests {
    use super::*;
    use catallax_shapes::{Circle, Line, Plane, Point3, PointDef, Surface};

    fn fixture() -> (Model, PointId, CurveId, CurveId) {
        let mut model = Model::new();
        let mut dec = Decider::default();
        let source = model.geom_mut().add_curve(Curve::Line(
            Line::new(Point3::origin(), unit(0), &mut dec).unwrap(),
        ));
        let barrier = model.geom_mut().add_curve(Curve::Circle(
            Circle::new(Point3::origin(), unit(2), Real::one(), unit(0), &mut dec).unwrap(),
        ));
        let point = model.geom_mut().add_point(PointDef::Explicit(Point3::new(
            Real::one(),
            Real::zero(),
            Real::zero(),
        )));
        (model, point, source, barrier)
    }
    #[test]
    fn direct_point_curve_fact_retains_analytic_contact() {
        let (mut model, point, source, barrier) = fixture();
        let source_equations = model.geom().curve(source).equations(model.geom()).unwrap();
        let equations = model.geom().curve(barrier).equations(model.geom()).unwrap();
        assert!(!equation_is_source_contact(
            &model,
            point,
            source,
            barrier,
            1,
            &equations[1],
            &source_equations
        ));
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Point(point),
            b: GeomKey::Curve(barrier),
            at: None,
            origin: FactOrigin::Proof(catallax_numbers::ZeroProof::ExactArithmetic),
        });
        assert!(equation_is_source_contact(
            &model,
            point,
            source,
            barrier,
            1,
            &equations[1],
            &source_equations
        ));
    }
    #[test]
    fn a_shared_plane_does_not_prove_the_other_circle_equation() {
        let (mut model, _, source, barrier) = fixture();
        let point = model
            .geom_mut()
            .add_point(PointDef::Explicit(Point3::origin()));
        let plane = model.geom_mut().add_surface(Surface::Plane(
            Plane::new(Point3::origin(), unit(2), unit(0), &mut Decider::default()).unwrap(),
        ));
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Point(point),
            b: GeomKey::Surface(plane),
            at: None,
            origin: FactOrigin::Proof(catallax_numbers::ZeroProof::ExactArithmetic),
        });
        let source_equations = model.geom().curve(source).equations(model.geom()).unwrap();
        let equations = model.geom().curve(barrier).equations(model.geom()).unwrap();
        assert!(!equation_is_source_contact(
            &model,
            point,
            source,
            barrier,
            1,
            &equations[1],
            &source_equations
        ));
    }
    #[test]
    fn opposite_source_equation_retains_the_same_zero_contact() {
        let (model, point, source, barrier) = fixture();
        let equations = model.geom().curve(barrier).equations(model.geom()).unwrap();
        let opposite = -&equations[0];
        assert!(equation_is_source_contact(
            &model,
            point,
            source,
            barrier,
            0,
            &equations[0],
            &[opposite]
        ));
    }
    #[test]
    fn coincident_curve_fact_retains_contact_under_a_different_definition() {
        let (mut model, point, source, _) = fixture();
        let barrier = model.geom_mut().add_curve(Curve::Line(
            Line::new(
                Point3::new(Real::from_i64(2), Real::zero(), Real::zero()),
                unit(0).scaled(&Real::from_i64(2)),
                &mut Decider::default(),
            )
            .unwrap(),
        ));
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Coincident,
            a: GeomKey::Curve(source),
            b: GeomKey::Curve(barrier),
            at: None,
            origin: FactOrigin::Proof(catallax_numbers::ZeroProof::ExactArithmetic),
        });
        let source_equations = model.geom().curve(source).equations(model.geom()).unwrap();
        let equations = model.geom().curve(barrier).equations(model.geom()).unwrap();
        for (index, equation) in equations.iter().enumerate() {
            assert!(equation_is_source_contact(
                &model,
                point,
                source,
                barrier,
                index,
                equation,
                &source_equations
            ));
        }
    }
}
