//! Exact Bernstein-to-power conversion without a global parameter round trip.

use super::{Dense, Value, mul, quotient, sub};
use crate::BSplineSurface;
use catallax_numbers::{Decider, Rational, Refusal};

pub(super) fn polynomials(
    source: &BSplineSurface,
    span: [usize; 2],
    dec: &mut Decider,
) -> Result<[Dense; 4], Refusal> {
    let controls = source.span_bezier_controls(span, dec)?;
    let mut coefficients: [Vec<Rational>; 4] =
        std::array::from_fn(|_| Vec::with_capacity(controls.len()));
    for row in controls {
        for (axis, value) in row.into_iter().enumerate() {
            coefficients[axis].push(value.as_rational().ok_or_else(|| {
                Refusal::new(
                    "shapes.rational-atlas-polynomial-coefficients",
                    "polynomial-bound admission needs exact rational homogeneous span coefficients",
                )
            })?);
        }
    }
    let degrees = source.degrees();
    let mut local = coefficients.map(|c| Dense { deg: degrees, c });
    let mut work = Vec::with_capacity(degrees[0].max(degrees[1]) + 1);
    for (axis, &degree) in degrees.iter().enumerate() {
        let lines = degrees[1 - axis] + 1;
        for polynomial in &mut local {
            for line in 0..lines {
                dec.check_deadline()?;
                let index = |k: usize| {
                    if axis == 0 {
                        k * (degrees[1] + 1) + line
                    } else {
                        line * (degrees[1] + 1) + k
                    }
                };
                work.clear();
                work.extend((0..=degree).map(|k| polynomial.c[index(k)].clone()));
                let mut factor = Rational::one();
                // Power coefficient k is C(degree,k) times forward difference k.
                for k in 0..=degree {
                    polynomial.c[index(k)] = mul(&work[0], &factor, dec)?;
                    let active = degree - k + 1;
                    if active == 1 {
                        break;
                    }
                    let mut constant = true;
                    for value in &work[1..active] {
                        if !<Rational as Value>::same(&work[0], value, dec)? {
                            constant = false;
                            break;
                        }
                    }
                    if constant {
                        for remaining in k + 1..=degree {
                            polynomial.c[index(remaining)] = Rational::zero();
                        }
                        break;
                    }
                    for i in 0..active - 1 {
                        work[i] = sub(&work[i + 1], &work[i], dec)?;
                    }
                    let numerator = mul(&factor, &Rational::from_i64((degree - k) as i64), dec)?;
                    factor = quotient(&numerator, &Rational::from_i64((k + 1) as i64), dec)?;
                }
            }
        }
    }
    Ok(local.map(Dense::trimmed))
}
