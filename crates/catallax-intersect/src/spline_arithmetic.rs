//! Caller-budgeted polynomial arithmetic for finite spline queries.

use catallax_numbers::{Decider, Outcome, Real, RealPoly, Refusal, Sign, StrictSign, ZeroProof};
use std::collections::{BTreeMap, btree_map::Entry};

pub(super) fn arithmetic(
    a: &Real,
    b: &Real,
    multiply: bool,
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    if let (Some(a), Some(b)) = (a.rational_ref(), b.rational_ref()) {
        dec.reserve_exact_arithmetic(a, b, multiply)?;
    } else {
        dec.charge_work(1)?;
    }
    Ok(if multiply { a * b } else { a + b })
}

pub(super) fn difference(a: &Real, b: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
    arithmetic(a, &(-b), false, dec)
}

pub(super) fn quotient(a: &Real, b: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
    if let (Some(a), Some(b)) = (a.rational_ref(), b.rational_ref()) {
        dec.reserve_exact_arithmetic(a, b, true)?;
    } else {
        dec.charge_work(1)?;
    }
    Ok(a / b)
}

pub(super) fn sign(value: &Real, dec: &mut Decider, code: &'static str) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => {
            Ok(Sign::Zero)
        }
        Outcome::ProvenCoincident(_) => Err(Refusal::new(
            code,
            "finite spline meetings need exact scalar equality",
        )),
        Outcome::NotExact(error) => Err(Refusal::from(error)),
        Outcome::Refused(error) => Err(error),
    }
}

pub(super) fn cut(
    coordinates: &[RealPoly; 4],
    factors: &[Real; 4],
    dec: &mut Decider,
) -> Result<RealPoly, Refusal> {
    let mut terms = BTreeMap::<Vec<u32>, Real>::new();
    for (coordinate, factor) in coordinates.iter().zip(factors) {
        for (monomial, coefficient) in coordinate.terms() {
            let term = arithmetic(&coefficient, factor, true, dec)?;
            match terms.entry(monomial.exponents().to_vec()) {
                Entry::Vacant(entry) => {
                    entry.insert(term);
                }
                Entry::Occupied(mut entry) => {
                    let sum = arithmetic(entry.get(), &term, false, dec)?;
                    entry.insert(sum);
                }
            }
        }
    }
    Ok(RealPoly::from_terms(1, terms.into_iter().collect()))
}

pub(super) fn value(
    poly: &RealPoly,
    at: &Real,
    derivative: bool,
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    derivative_value(poly, at, usize::from(derivative), dec)
}

pub(super) fn derivative_value(
    poly: &RealPoly,
    at: &Real,
    order: usize,
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    let mut sum = Real::zero();
    for (monomial, coefficient) in poly.terms() {
        let mut power = monomial.exponents()[0];
        if u64::from(power) < order as u64 {
            continue;
        }
        let mut term = coefficient;
        for _ in 0..order {
            term = arithmetic(&term, &Real::from_i64(i64::from(power)), true, dec)?;
            power -= 1;
        }
        for _ in 0..power {
            term = arithmetic(&term, at, true, dec)?;
        }
        sum = arithmetic(&sum, &term, false, dec)?;
    }
    Ok(sum)
}
