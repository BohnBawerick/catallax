//! Complete zeros of an analytic cutting equation on finite spline spans.

use super::*;
use catallax_numbers::Rational;
use std::collections::{BTreeMap, btree_map::Entry};

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

// Polynomial operators do not reserve work for their eager rational arithmetic.
fn exact(a: &Real, b: &Real, multiply: bool, dec: &mut Decider) -> Result<Real, Refusal> {
    if let (Some(a), Some(b)) = (a.rational_ref(), b.rational_ref()) {
        dec.reserve_exact_arithmetic(a, b, multiply)?;
    } else {
        dec.charge_work(1)?;
    }
    Ok(if multiply { a * b } else { a + b })
}

fn accumulate(
    coefficients: &mut BTreeMap<Vec<u32>, Real>,
    exponents: Vec<u32>,
    coefficient: Real,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if coefficient.rational_ref().is_some_and(Rational::is_zero) {
        return Ok(());
    }
    match coefficients.entry(exponents) {
        Entry::Vacant(entry) => {
            entry.insert(coefficient);
        }
        Entry::Occupied(mut entry) => {
            let sum = exact(entry.get(), &coefficient, false, dec)?;
            if sum.rational_ref().is_some_and(Rational::is_zero) {
                entry.remove();
            } else {
                entry.insert(sum);
            }
        }
    }
    Ok(())
}

fn product(a: &RealPoly, b: &RealPoly, dec: &mut Decider) -> Result<RealPoly, Refusal> {
    let mut coefficients = BTreeMap::new();
    let right = b.terms();
    for (m, x) in a.terms() {
        for (n, y) in &right {
            let exponents = m
                .exponents()
                .iter()
                .zip(n.exponents())
                .map(|(a, b)| a.checked_add(*b).expect("polynomial exponent overflow"))
                .collect();
            let coefficient = exact(&x, y, true, dec)?;
            accumulate(&mut coefficients, exponents, coefficient, dec)?;
        }
    }
    Ok(RealPoly::from_terms(
        a.nvars(),
        coefficients.into_iter().collect(),
    ))
}

fn sum(a: &RealPoly, b: &RealPoly, dec: &mut Decider) -> Result<RealPoly, Refusal> {
    let mut coefficients = a
        .terms()
        .into_iter()
        .map(|(m, x)| (m.exponents().to_vec(), x))
        .collect();
    for (m, x) in b.terms() {
        accumulate(&mut coefficients, m.exponents().to_vec(), x, dec)?;
    }
    Ok(RealPoly::from_terms(
        a.nvars(),
        coefficients.into_iter().collect(),
    ))
}

fn substitute(
    equation: &RealPoly,
    coordinates: &[RealPoly; 4],
    dec: &mut Decider,
) -> Result<RealPoly, Refusal> {
    let terms = equation.terms();
    let degree = terms.iter().map(|(m, _)| m.degree()).max().unwrap_or(0);
    let mut result = RealPoly::zero(1);
    for (monomial, coefficient) in terms {
        let mut term = RealPoly::constant(1, coefficient);
        for (axis, &power) in monomial.exponents().iter().enumerate() {
            for _ in 0..power {
                term = product(&term, &coordinates[axis], dec)?;
            }
        }
        for _ in 0..degree - monomial.degree() {
            term = product(&term, &coordinates[3], dec)?;
        }
        result = sum(&result, &term, dec)?;
    }
    Ok(result)
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
        let equation = substitute(equation, &coordinates, dec)?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    fn cut(scale: i64, contained: bool) -> (RealPoly, [RealPoly; 4]) {
        let factor = Real::from_rational(Rational::one().mul_pow2(scale));
        let dense = RealPoly::from_terms(
            1,
            (0..=6)
                .map(|power| (vec![power], Real::from_i64(i64::from(power) + 1)))
                .collect(),
        )
        .scale(&factor);
        let coordinates = [
            dense.clone(),
            dense.clone(),
            &(&dense * &dense) + &RealPoly::constant(1, Real::from_i64(i64::from(!contained))),
            RealPoly::constant(1, Real::one()),
        ];
        let x = RealPoly::var(3, 0);
        let z = RealPoly::var(3, 2);
        (&(&x * &x) - &z, coordinates)
    }

    #[test]
    fn zero_cut_expansion_cannot_bypass_a_zero_work_allowance() {
        let (cutter, coordinates) = cut(0, true);
        let (answer, work) = with_operation_budget(0, || {
            substitute(&cutter, &coordinates, &mut Decider::default())
        });
        assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
        assert!(
            work.exhausted,
            "coefficient expansion escaped its allowance"
        );
    }

    #[test]
    fn dense_cut_expansion_reserves_coefficient_work_and_keeps_exact_answers() {
        for contained in [false, true] {
            let mut costs = Vec::new();
            for scale in [0, 256] {
                let (cutter, coordinates) = cut(scale, contained);
                let (answer, work) = with_operation_budget(20_000_000, || {
                    substitute(&cutter, &coordinates, &mut Decider::default())
                });
                let expected = RealPoly::constant(1, Real::from_i64(-i64::from(!contained)));
                assert_eq!(answer.unwrap().rational(), expected.rational());
                assert!(!work.exhausted);
                assert!(work.used > 1);
                for limit in [1, work.used - 1] {
                    let (answer, stopped) = with_operation_budget(limit, || {
                        substitute(&cutter, &coordinates, &mut Decider::default())
                    });
                    assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
                    assert!(stopped.exhausted);
                    assert!(stopped.used <= limit);
                }
                println!("contained={contained} scale={scale}: {} units", work.used);
                costs.push(work.used);
            }
            assert!(
                costs[1] > costs[0],
                "large coefficients need more reserved work"
            );
        }
    }
}
