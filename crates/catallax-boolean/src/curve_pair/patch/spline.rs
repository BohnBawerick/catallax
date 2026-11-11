//! Complete zeros of an analytic cutting equation on finite spline spans.

use super::*;

fn nonnegative(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(catallax_numbers::StrictSign::Positive) => Ok(true),
        Outcome::Proven(catallax_numbers::StrictSign::Negative) => Ok(false),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => Ok(true),
        Outcome::ProvenCoincident(_) => Err(unresolved(
            "a spline meeting needs exact knot-span membership",
        )),
        Outcome::NotExact(error) => Err(Refusal::from(error)),
        Outcome::Refused(error) => Err(error),
    }
}

fn substitute(equation: &RealPoly, coordinates: &[RealPoly; 4]) -> RealPoly {
    let terms = equation.terms();
    let degree = terms.iter().map(|(m, _)| m.degree()).max().unwrap_or(0);
    let mut sum = RealPoly::zero(1);
    for (monomial, coefficient) in terms {
        let mut term = RealPoly::constant(1, coefficient);
        for (axis, &power) in monomial.exponents().iter().enumerate() {
            for _ in 0..power {
                term = &term * &coordinates[axis];
            }
        }
        for _ in 0..degree - monomial.degree() {
            term = &term * &coordinates[3];
        }
        sum = sum + term;
    }
    sum
}

// None means every span satisfies the implicit cutting equation. No surface
// incidence fact is published here: a cone may also have an excluded nappe.
pub(super) fn hits(
    curve: &BSplineCurve,
    equation: &RealPoly,
    dec: &mut Decider,
) -> Result<Option<Vec<Point3>>, Refusal> {
    let mut contained = 0;
    let mut spans = 0;
    let mut points = Vec::new();
    for (span, [lo, hi]) in curve.spans() {
        dec.check_deadline()?;
        spans += 1;
        let coordinates = curve.query_span_polynomials(span, dec)?;
        // Positive weights clear the denominator without losing a root.
        let equation = substitute(equation, &coordinates);
        let rational = equation.rational().ok_or_else(|| {
            unresolved("finite spline/cutter meetings need rational homogeneous cut coefficients")
        })?;
        if rational.is_zero() {
            contained += 1;
            continue;
        }
        for parameter in equation.real_roots(dec)? {
            if nonnegative(&(&parameter - &Real::from_rational(lo.clone())), dec)?
                && nonnegative(&(&Real::from_rational(hi.clone()) - &parameter), dec)?
            {
                points.push(curve.point_at(&parameter, dec)?);
            }
        }
    }
    if contained == spans {
        Ok(None)
    } else if contained != 0 {
        Err(unresolved(
            "a spline lies on the cutting equation over only some spans; a partial-overlap solver is missing",
        ))
    } else {
        Ok(Some(points))
    }
}
