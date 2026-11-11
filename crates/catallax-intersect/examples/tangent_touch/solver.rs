//! Opt-in R2 experiment. No production caller uses these reports.

use catallax_numbers::{
    Decider, MPoly, Outcome, Rational, Real, RealPoly, Refusal, Sign, StrictSign,
};
use catallax_shapes::{Plane, RationalPatch};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    IsolatedTouch,
    TwoBranchesCrossing,
    DegenerateContact,
}

impl Kind {
    pub fn local_outcome(self) -> Result<(), Refusal> {
        match self {
            Self::IsolatedTouch => Ok(()),
            Self::TwoBranchesCrossing => Err(refusal("prototype.tangent-self-touch")),
            Self::DegenerateContact => Err(refusal("prototype.tangent-degenerate-contact")),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Contact {
    pub parameters: [Real; 2],
    pub kind: Kind,
}

pub fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(n.into()),
        Outcome::Refused(r) => Err(r),
    }
}

fn refusal(reason: &'static str) -> Refusal {
    Refusal::new(
        reason,
        "R2 contact enumeration is not a complete Boolean certificate",
    )
}

// Reserve every coefficient operation with the existing arithmetic meter.
// MPoly's public resultant does not reserve its own arithmetic.
fn combine(a: &MPoly, b: &MPoly, multiply: bool, dec: &mut Decider) -> Result<MPoly, Refusal> {
    let mut terms = BTreeMap::<Vec<u32>, Rational>::new();
    if multiply {
        for (x, c) in a.terms() {
            for (y, d) in b.terms() {
                dec.reserve_exact_arithmetic(&c, &d, true)?;
                let value = &c * &d;
                let powers = x
                    .exponents()
                    .iter()
                    .zip(y.exponents())
                    .map(|(i, j)| i + j)
                    .collect();
                let previous = terms.entry(powers).or_insert_with(Rational::zero);
                dec.reserve_exact_arithmetic(previous, &value, false)?;
                *previous = &*previous + &value;
            }
        }
    } else {
        for (m, c) in a.terms().into_iter().chain(b.terms()) {
            let previous = terms
                .entry(m.exponents().to_vec())
                .or_insert_with(Rational::zero);
            dec.reserve_exact_arithmetic(previous, &c, false)?;
            *previous = &*previous + &c;
        }
    }
    Ok(MPoly::from_terms(2, terms.into_iter().collect()))
}

fn scaled(poly: &MPoly, scalar: &Rational, dec: &mut Decider) -> Result<MPoly, Refusal> {
    let mut terms = Vec::new();
    for (monomial, coefficient) in poly.terms() {
        dec.reserve_exact_arithmetic(&coefficient, scalar, true)?;
        terms.push((monomial.exponents().to_vec(), &coefficient * scalar));
    }
    Ok(MPoly::from_terms(2, terms))
}

fn derivative(poly: &MPoly, axis: usize, dec: &mut Decider) -> Result<MPoly, Refusal> {
    let mut terms = Vec::new();
    for (monomial, coefficient) in poly.terms() {
        let mut powers = monomial.exponents().to_vec();
        if powers[axis] == 0 {
            continue;
        }
        let factor = Rational::from_i64(i64::from(powers[axis]));
        dec.reserve_exact_arithmetic(&coefficient, &factor, true)?;
        powers[axis] -= 1;
        terms.push((powers, &coefficient * &factor));
    }
    Ok(MPoly::from_terms(2, terms))
}

fn coefficients(poly: &MPoly, axis: usize) -> Vec<MPoly> {
    let terms = poly.terms();
    let degree = terms
        .iter()
        .map(|(m, _)| m.exponents()[axis])
        .max()
        .unwrap_or(0);
    (0..=degree)
        .map(|power| {
            MPoly::from_terms(
                2,
                terms
                    .iter()
                    .filter_map(|(m, c)| {
                        let mut e = m.exponents().to_vec();
                        if e[axis] != power {
                            return None;
                        }
                        e[axis] = 0;
                        Some((e, c.clone()))
                    })
                    .collect(),
            )
        })
        .collect()
}

fn eliminant(a: &MPoly, b: &MPoly, axis: usize, dec: &mut Decider) -> Result<MPoly, Refusal> {
    if a.is_zero() || b.is_zero() {
        return Ok(MPoly::zero(2));
    }
    let a = coefficients(a, axis);
    let b = coefficients(b, axis);
    let (m, n) = (a.len() - 1, b.len() - 1);
    // Two equations independent of this axis give no projection certificate.
    if m + n == 0 {
        return Ok(MPoly::zero(2));
    }
    let size = m + n;
    if size > 10 {
        return Err(refusal("prototype.tangent-elimination-size"));
    }
    let mut matrix = vec![vec![MPoly::zero(2); size]; size];
    for i in 0..n {
        for j in 0..=m {
            matrix[i][i + j] = a[m - j].clone();
        }
    }
    for i in 0..m {
        for j in 0..=n {
            matrix[n + i][i + j] = b[n - j].clone();
        }
    }
    // Subset determinant avoids unmetered polynomial division in Bareiss.
    let mut sums = vec![MPoly::zero(2); 1 << size];
    sums[0] = MPoly::constant(2, Rational::one());
    for mask in 0usize..(1 << size) - 1 {
        dec.charge_work(1)?;
        if sums[mask].is_zero() {
            continue;
        }
        let row = mask.count_ones() as usize;
        for (col, coefficient) in matrix[row].iter().enumerate() {
            if mask & (1 << col) != 0 || coefficient.is_zero() {
                continue;
            }
            let mut term = combine(&sums[mask], coefficient, true, dec)?;
            if (mask >> (col + 1)).count_ones() % 2 == 1 {
                term = scaled(&term, &Rational::from_i64(-1), dec)?;
            }
            let next = mask | (1 << col);
            sums[next] = combine(&sums[next], &term, false, dec)?;
        }
    }
    Ok(sums.pop().expect("nonempty determinant"))
}

fn roots(poly: &MPoly, axis: usize, dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
    let poly = poly.to_upoly(axis).expect("eliminated coordinate");
    RealPoly::from(MPoly::from_terms(
        1,
        poly.coefficients()
            .iter()
            .enumerate()
            .map(|(i, c)| (vec![i as u32], c.clone()))
            .collect(),
    ))
    .real_roots(dec)
}

fn in_range(value: &Real, bounds: [&Rational; 2], dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(
        sign(&(value - &Real::from_rational(bounds[0].clone())), dec)? != Sign::Negative
            && sign(&(&Real::from_rational(bounds[1].clone()) - value), dec)? != Sign::Negative,
    )
}

fn plane_equation(h: &[RealPoly; 4], plane: &Plane, dec: &mut Decider) -> Result<MPoly, Refusal> {
    let h = h
        .iter()
        .map(|coordinate| {
            coordinate
                .rational()
                .ok_or_else(|| refusal("prototype.tangent-nonrational"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z];
    let origin = [&plane.origin.x, &plane.origin.y, &plane.origin.z];
    let mut equation = MPoly::zero(2);
    for i in 0..3 {
        let origin = origin[i]
            .rational_ref()
            .ok_or_else(|| refusal("prototype.tangent-nonrational"))?;
        let normal = normal[i]
            .rational_ref()
            .ok_or_else(|| refusal("prototype.tangent-nonrational"))?;
        let offset = scaled(h[3], origin, dec)?;
        let negative = scaled(&offset, &Rational::from_i64(-1), dec)?;
        let residual = combine(h[i], &negative, false, dec)?;
        let term = scaled(&residual, normal, dec)?;
        equation = combine(&equation, &term, false, dec)?;
    }
    Ok(equation)
}

pub fn contacts(
    patch: &RationalPatch,
    plane: &Plane,
    dec: &mut Decider,
) -> Result<Vec<Contact>, Refusal> {
    if sign(&plane.normal.norm_squared(), dec)? != Sign::Positive {
        return Err(refusal("prototype.tangent-plane-normal"));
    }
    let carrier = patch.carrier();
    let knots = carrier.knots();
    for (axis, axis_knots) in knots.iter().enumerate() {
        let mut previous = None;
        let mut multiplicity = 0;
        let domain = carrier.domain()[axis];
        for knot in *axis_knots {
            dec.charge_work(1)?;
            multiplicity = if previous == Some(knot) {
                multiplicity + 1
            } else {
                1
            };
            previous = Some(knot);
            if knot > domain[0] && knot < domain[1] && multiplicity >= carrier.degrees()[axis] {
                return Err(refusal("prototype.tangent-nonsmooth-knot"));
            }
        }
    }
    let mut found: Vec<Contact> = Vec::new();
    for span in carrier.span_indices() {
        dec.check_deadline()?;
        let mut side = None;
        let mut separated = true;
        for [x, y, z, w] in carrier.span_bezier_controls(span, dec)? {
            let value = &(&x - &(&plane.origin.x * &w)) * &plane.normal.x
                + &(&y - &(&plane.origin.y * &w)) * &plane.normal.y
                + &(&z - &(&plane.origin.z * &w)) * &plane.normal.z;
            let current = sign(&value, dec)?;
            if current == Sign::Zero || side.is_some_and(|previous| previous != current) {
                separated = false;
                break;
            }
            side = Some(current);
        }
        if separated {
            continue;
        }
        let h = carrier.span_polynomials(span, dec)?;
        let f = plane_equation(&h, plane, dec)?;
        if f.is_zero() {
            return Err(refusal("prototype.tangent-flush-span"));
        }
        if f.total_degree() == Some(0) {
            continue;
        }
        let df = [derivative(&f, 0, dec)?, derivative(&f, 1, dec)?];
        if df.iter().any(|d| d.total_degree() == Some(0)) {
            continue;
        }
        let projections = [
            eliminant(&df[0], &df[1], 1, dec)?,
            eliminant(&df[0], &df[1], 0, dec)?,
        ];
        if projections.iter().any(MPoly::is_zero) {
            return Err(refusal("prototype.tangent-stationary-component"));
        }
        let us = roots(&projections[0], 0, dec)?;
        let vs = roots(&projections[1], 1, dec)?;
        let bounds = [0, 1].map(|i| [&knots[i][span[i]], &knots[i][span[i] + 1]]);
        let equations: [RealPoly; 3] =
            [f.clone().into(), df[0].clone().into(), df[1].clone().into()];
        for u in &us {
            if !in_range(u, bounds[0], dec)? {
                continue;
            }
            for v in &vs {
                if !in_range(v, bounds[1], dec)? {
                    continue;
                }
                let at = [u.clone(), v.clone()];
                let mut valid = true;
                for equation in &equations {
                    if sign(&equation.eval_real(&at), dec)? != Sign::Zero {
                        valid = false;
                        break;
                    }
                }
                if !valid {
                    continue;
                }
                // A one-sided Hessian cannot certify a joined or trimmed neighbourhood.
                for i in 0..2 {
                    for bound in bounds[i] {
                        if sign(&(&at[i] - &Real::from_rational(bound.clone())), dec)? == Sign::Zero
                        {
                            return Err(refusal("prototype.tangent-boundary-contact"));
                        }
                    }
                }
                let a = RealPoly::from(derivative(&df[0], 0, dec)?).eval_real(&at);
                let b = RealPoly::from(derivative(&df[0], 1, dec)?).eval_real(&at);
                let c = RealPoly::from(derivative(&df[1], 1, dec)?).eval_real(&at);
                let kind = match sign(&(&a * &c - &b * &b), dec)? {
                    Sign::Positive => Kind::IsolatedTouch,
                    Sign::Negative => Kind::TwoBranchesCrossing,
                    Sign::Zero => Kind::DegenerateContact,
                };
                // Shared span edges were refused above, so spans cannot duplicate points.
                found.push(Contact {
                    parameters: at,
                    kind,
                });
            }
        }
    }
    Ok(found)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derivative_reserves_before_multiplying_coefficients() {
        let u = MPoly::var(2, 0);
        let poly = &u * &u;
        let (answer, work) = catallax_numbers::with_operation_budget(0, || {
            derivative(&poly, 0, &mut Decider::default())
        });
        assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
        assert!(work.exhausted);
        let (answer, work) = catallax_numbers::with_operation_budget(100, || {
            derivative(&poly, 0, &mut Decider::default())
        });
        assert_eq!(answer.unwrap(), poly.partial(0));
        assert!(work.used > 0);
        assert!(!work.exhausted);
    }

    #[test]
    fn determinant_negation_reserves_before_scaling_coefficients() {
        let poly = MPoly::var(2, 0);
        let minus_one = Rational::from_i64(-1);
        let (answer, work) = catallax_numbers::with_operation_budget(0, || {
            scaled(&poly, &minus_one, &mut Decider::default())
        });
        assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
        assert!(work.exhausted);
        let (answer, work) = catallax_numbers::with_operation_budget(100, || {
            scaled(&poly, &minus_one, &mut Decider::default())
        });
        assert_eq!(answer.unwrap(), -poly);
        assert!(work.used > 0);
        assert!(!work.exhausted);
    }

    #[test]
    fn plane_equation_reserves_its_rational_arithmetic() {
        let h = [
            RealPoly::var(2, 0),
            RealPoly::var(2, 1),
            RealPoly::var(2, 0).pow_u32(2),
            RealPoly::constant(2, Real::one()),
        ];
        let plane = Plane {
            origin: catallax_shapes::Point3::new(Real::one(), Real::one(), Real::one()),
            normal: catallax_shapes::Vector3::new(Real::one(), Real::one(), Real::one()),
            ref_dir: catallax_shapes::Vector3::new(Real::one(), Real::zero(), Real::zero()),
        };
        let (answer, work) = catallax_numbers::with_operation_budget(0, || {
            plane_equation(&h, &plane, &mut Decider::default())
        });
        assert!(
            answer.is_err(),
            "plane equation arithmetic ran without an allowance"
        );
        assert!(work.exhausted);
    }

    #[test]
    fn metered_determinants_match_the_existing_fraction_free_resultant() {
        let u = MPoly::var(2, 0);
        let v = MPoly::var(2, 1);
        let one = MPoly::constant(2, Rational::one());
        for (a, b) in [
            (&(&u * &u) + &(&v * &v) - one.clone(), &u - &v),
            (&u * &v + one.clone(), &u * &v + u.clone()),
            (&(&u - &one) * &(&v - &one), &(&u - &one) * &(&v + &one)),
        ] {
            for axis in 0..2 {
                let mut dec = Decider::default();
                assert_eq!(
                    eliminant(&a, &b, axis, &mut dec).unwrap(),
                    a.resultant(&b, axis)
                );
            }
        }
    }
}
