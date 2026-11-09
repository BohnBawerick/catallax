//! Exact normal-plane sections of finite polynomial extrusion carriers.

use super::{Contact, CurveForm, Plan};
use catallax_numbers::{Decider, Outcome, Real, Refusal, Sign, StrictSign, ZeroProof};
use catallax_shapes::{BSplineCurve, BezierExtrusion, Plane};

fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(p) if p.proof != ZeroProof::SeparationBound => Ok(Sign::Zero),
        Outcome::ProvenCoincident(_) => Err(Refusal::new(
            "numbers.geometric-gap-unproven",
            "a finite Bezier section needs an exact incidence decision",
        )),
        Outcome::NotExact(error) => Err(Refusal::new(error.code, error.detail)),
        Outcome::Refused(error) => Err(error),
    }
}

pub(super) fn plane(
    plane: &Plane,
    support: &BezierExtrusion,
    swap: bool,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let direction = support.direction();
    let mut sides = Vec::with_capacity(support.controls().len() * 2);
    for point in support
        .controls()
        .iter()
        .flat_map(|p| [p.clone(), p.plus(direction)])
    {
        sides.push(sign(&point.minus(&plane.origin).dot(&plane.normal), dec)?);
    }
    // Every point has a positive Bernstein weight at some corner. Strict
    // corner values and weakly same-sided controls exclude the whole patch.
    if let Some(side) = sides.iter().copied().find(|s| *s != Sign::Zero)
        && sides.iter().all(|s| *s == Sign::Zero || *s == side)
        && [0, 1, sides.len() - 2, sides.len() - 1]
            .into_iter()
            .all(|i| sides[i] == side)
    {
        return Ok(Plan::disjoint());
    }
    if sign(&plane.normal.cross(direction).norm_squared(), dec)? != Sign::Zero {
        return Err(Refusal::new(
            "intersect.bezier-plane-oblique",
            "Bezier extrusion sections currently require a plane normal to the generators",
        ));
    }
    let denominator = plane.normal.dot(direction);
    let parameter = &plane
        .origin
        .minus(&support.controls()[0])
        .dot(&plane.normal)
        / &denominator;
    if sign(&parameter, dec)? == Sign::Negative
        || sign(&(&Real::one() - &parameter), dec)? == Sign::Negative
    {
        return Ok(Plan::disjoint());
    }
    let offset = direction.scaled(&parameter);
    let mut controls: Vec<_> = support.controls().iter().map(|p| p.plus(&offset)).collect();
    if (sign(&denominator, dec)? == Sign::Negative) != swap {
        controls.reverse();
    }
    let lower = support.lower_edge();
    let curve = BSplineCurve::new(
        lower.degree(),
        lower.knots().to_vec(),
        controls,
        vec![Real::one(); lower.control_count()],
        dec,
    )?;
    Ok(Plan::form(CurveForm::BSpline(curve), Contact::Crossing))
}
