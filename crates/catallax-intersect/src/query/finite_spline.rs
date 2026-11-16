//! Whole finite span identities and exact off-curve witnesses.

use catallax_numbers::{Decider, Rational, Real, RealPoly, Refusal};
use catallax_shapes::BSplineCurve;

struct Span {
    domain: [Rational; 2],
    coordinates: [RealPoly; 4],
}

fn spans(curve: &BSplineCurve, dec: &mut Decider) -> Result<Vec<Span>, Refusal> {
    curve
        .spans()
        .map(|(index, domain)| {
            Ok(Span {
                domain: domain.map(Clone::clone),
                coordinates: curve.query_span_polynomials(index, dec)?,
            })
        })
        .collect()
}

fn affine_identity(
    a: &[Span],
    b: &[Span],
    offset: &Rational,
    scale: &Rational,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let parameter = RealPoly::var(1, 0).scale(&Real::from_rational(scale.clone()))
        + RealPoly::constant(1, Real::from_rational(offset.clone()));
    for second in b {
        let mut domain = second.domain.each_ref().map(|t| &(t - offset) / scale);
        if domain[0] > domain[1] {
            domain.swap(0, 1);
        }
        let coordinates = second
            .coordinates
            .each_ref()
            .map(|p| p.substitute(0, &parameter));
        for first in a {
            if std::cmp::max(&first.domain[0], &domain[0])
                >= std::cmp::min(&first.domain[1], &domain[1])
            {
                continue;
            }
            for axis in 0..3 {
                if !(&first.coordinates[axis] * &coordinates[3]
                    - &coordinates[axis] * &first.coordinates[3])
                    .is_zero(dec)?
                {
                    return Ok(false);
                }
            }
        }
    }
    // The affine map covers both closed domains. Span interiors partition
    // them, and admitted C0 continuity supplies all shared knot endpoints.
    Ok(true)
}

pub(super) fn same(
    a: &BSplineCurve,
    b: &BSplineCurve,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    let first = spans(a, dec)?;
    let second = spans(b, dec)?;
    let [a_lo, a_hi] = a.domain();
    let [b_lo, b_hi] = b.domain();
    let scale = (b_hi - b_lo) / (a_hi - a_lo);
    for forward in [true, false] {
        let scale = if forward { scale.clone() } else { -&scale };
        let offset = if forward { b_lo } else { b_hi } - &(a_lo * &scale);
        if affine_identity(&first, &second, &offset, &scale, dec)? {
            return Ok(Some(true));
        }
    }
    for (source, target) in [(a, b), (b, a)] {
        let mut parameters = vec![source.domain()[0].clone(), source.domain()[1].clone()];
        parameters.extend(source.spans().map(|(_, [lo, hi])| (lo + hi).mul_pow2(-1)));
        for parameter in parameters {
            let point = source.point_at(&Real::from_rational(parameter), dec)?;
            if target.parameter_of(&point, dec)?.is_none() {
                return Ok(Some(false));
            }
        }
    }
    // Matching witness points never establishes equality of complete curves.
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;
    use catallax_shapes::Point3;
    use std::time::{Duration, Instant};

    #[test]
    fn many_span_affine_identity_charges_visits_and_checks_the_deadline() {
        let count = 128_i64;
        let mut dec = Decider::default();
        let controls: Vec<_> = (0..=2 * count)
            .map(|k| {
                let i = k / 2;
                Point3::new(
                    Rational::from_ratio(k, 2).unwrap().into(),
                    Real::from_i64(i * i + if k % 2 == 1 { i } else { 0 }),
                    Real::zero(),
                )
            })
            .collect();
        let knots: Vec<_> = (0..=count)
            .flat_map(|i| {
                std::iter::repeat_n(
                    Rational::from_i64(i),
                    if i == 0 || i == count { 3 } else { 2 },
                )
            })
            .collect();
        let first = BSplineCurve::new(
            2,
            knots.clone(),
            controls.clone(),
            vec![Real::one(); controls.len()],
            &mut dec,
        )
        .unwrap();
        let first = spans(&first, &mut dec).unwrap();
        for reversed in [false, true] {
            let mut controls = controls.clone();
            if reversed {
                controls.reverse();
            }
            let second = BSplineCurve::new(
                2,
                knots
                    .iter()
                    .map(|t| t * &Rational::from_i64(3) + Rational::from_i64(7))
                    .collect(),
                controls.clone(),
                vec![Real::one(); controls.len()],
                &mut dec,
            )
            .unwrap();
            let second = spans(&second, &mut dec).unwrap();
            let scale = Rational::from_i64(if reversed { -3 } else { 3 });
            let offset = Rational::from_i64(if reversed { 7 + 3 * count } else { 7 });
            let (result, work) = with_operation_budget(100_000, || {
                affine_identity(&first, &second, &offset, &scale, &mut Decider::default())
            });
            assert!(result.unwrap());
            assert!(!work.exhausted);
            assert!(
                work.used >= count as u64,
                "span visits must be charged: {work:?}"
            );
            let (result, work) = with_operation_budget(16, || {
                affine_identity(&first, &second, &offset, &scale, &mut Decider::default())
            });
            assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
            assert!(work.exhausted);
            let (result, _) = with_operation_budget(100_000, || {
                let mut dec = Decider::default();
                dec.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
                affine_identity(&first, &second, &offset, &scale, &mut dec)
            });
            assert_eq!(result.unwrap_err().code, "numbers.deadline");
        }
    }
}
