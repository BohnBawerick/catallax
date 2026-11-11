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
        Outcome::NotExact(error) => Err(Refusal::from(error)),
        Outcome::Refused(error) => Err(error),
    }
}

enum Hull {
    Separated,
    Contact,
    Unproved,
}

fn hull(rows: &[Vec<Real>; 2], dec: &mut Decider) -> Result<Hull, Refusal> {
    let sides: [Vec<_>; 2] = [
        rows[0]
            .iter()
            .map(|value| sign(value, dec))
            .collect::<Result<_, _>>()?,
        rows[1]
            .iter()
            .map(|value| sign(value, dec))
            .collect::<Result<_, _>>()?,
    ];
    let last = rows[0].len() - 1;
    let corners = [sides[0][0], sides[0][last], sides[1][0], sides[1][last]];
    let side = corners[0];
    // These are actual patch points. A zero or opposite corner prevents a
    // complete separation proof, so genuine sections do not exhaust the budget.
    if side == Sign::Zero || corners.iter().any(|corner| *corner != side) {
        return Ok(Hull::Contact);
    }
    Ok(
        if sides
            .iter()
            .flatten()
            .all(|s| *s == Sign::Zero || *s == side)
        {
            Hull::Separated
        } else {
            Hull::Unproved
        },
    )
}

fn split_row(values: &[Real]) -> [Vec<Real>; 2] {
    let mut row = values.to_vec();
    let mut left = vec![row[0].clone()];
    let mut right = vec![row.last().unwrap().clone()];
    for count in (1..row.len()).rev() {
        for i in 0..count {
            row[i] = (&row[i] + &row[i + 1]) / Real::from_i64(2);
        }
        left.push(row[0].clone());
        right.push(row[count - 1].clone());
    }
    right.reverse();
    [left, right]
}

fn subdivided_separation(
    plane: &Plane,
    support: &BezierExtrusion,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let lower: Vec<_> = support
        .controls()
        .iter()
        .map(|point| point.minus(&plane.origin).dot(&plane.normal))
        .collect();
    let shift = support.direction().dot(&plane.normal);
    let upper = lower.iter().map(|value| value + &shift).collect();
    let mut pending = vec![([lower, upper], 0)];
    let mut cells = 0;
    while let Some((rows, depth)) = pending.pop() {
        dec.check_deadline()?;
        cells += 1;
        if cells > 4096 {
            return Ok(false);
        }
        match hull(&rows, dec)? {
            Hull::Separated => continue,
            Hull::Contact => return Ok(false),
            Hull::Unproved => {
                if depth >= 48 {
                    return Ok(false);
                }
                let [a, b] = split_row(&rows[0]);
                let [c, d] = split_row(&rows[1]);
                pending.push(([b, d], depth + 1));
                pending.push(([a, c], depth + 1));
            }
        }
    }
    Ok(true)
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
        if subdivided_separation(plane, support, dec)? {
            return Ok(Plan::disjoint());
        }
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
