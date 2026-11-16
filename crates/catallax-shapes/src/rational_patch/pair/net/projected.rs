//! Coupled rational Bernstein certificates for a chart's existing projection.
//!
//! The projected Jacobian is N/W². Bounds are taken after the exact coordinate
//! combinations and quotient-rule cancellations, not on independent XYZ jets.
//! Positive diagonal bounds and ab - c²/4 > 0 prove a positive symmetric part.
//! Its determinant/trace bound, divided by max(W)², retains the same global
//! strong-monotonicity invariant as interval admission. No sample is a proof.

use super::admission::{bernstein, quotient, restricted_range, restriction};
use super::*;

pub(crate) struct ProjectedJacobian {
    domain: Domain,
    widths: [Rational; 2],
    projection_scale: Rational,
    weight: Grid,
    diagonal: [Grid; 2],
    off_sum: Grid,
}

fn direction(vector: &Vector3) -> Result<[&Rational; 3], Refusal> {
    let missing = || {
        Refusal::new(
            "shapes.rational-patch-unproven",
            "the coupled Bernstein certificate needs a rational centre projection",
        )
    };
    Ok([
        vector.x.rational_ref().ok_or_else(missing)?,
        vector.y.rational_ref().ok_or_else(missing)?,
        vector.z.rational_ref().ok_or_else(missing)?,
    ])
}

// Both projection rows use one positive scale, so the symmetric Jacobian's
// positivity is unchanged. Convert its margin back before publishing it.
fn integral_projection(
    projection: &[Vector3; 2],
    dec: &mut Decider,
) -> Result<([Vector3; 2], Rational), Refusal> {
    let mut scale = Rational::one();
    for vector in projection {
        for component in direction(vector)? {
            let denominator = Rational::from_integer(component.denominator().clone());
            // For positive integers C,D: lcm(C,D) = numerator(C/D) * D.
            // Division and multiplication retain their existing work charges.
            let ratio = quotient(&scale, &denominator, dec)?;
            let cofactor = Rational::from_integer(ratio.numerator().clone());
            scale = mul(&cofactor, &denominator, dec)?;
        }
    }
    let row = |vector: &Vector3, dec: &mut Decider| -> Result<Vector3, Refusal> {
        let components = direction(vector)?;
        Ok(Vector3::new(
            mul(components[0], &scale, dec)?.into(),
            mul(components[1], &scale, dec)?.into(),
            mul(components[2], &scale, dec)?.into(),
        ))
    };
    Ok((
        [row(&projection[0], dec)?, row(&projection[1], dec)?],
        scale,
    ))
}

fn combine(terms: &[(&Dense, &Rational)], dec: &mut Decider) -> Result<Dense, Refusal> {
    let mut degree = [0usize; 2];
    for (poly, factor) in terms {
        if !factor.is_zero() {
            degree = std::array::from_fn(|axis| degree[axis].max(poly.deg[axis]));
        }
    }
    let mut result = Dense::zero(degree);
    for (poly, factor) in terms {
        if factor.is_zero() {
            continue;
        }
        for i in 0..=poly.deg[0] {
            for j in 0..=poly.deg[1] {
                let coefficient = poly.at(i, j);
                if coefficient.is_zero() {
                    continue;
                }
                let term = mul(coefficient, factor, dec)?;
                let value = add(result.at(i, j), &term, dec)?;
                *result.at_mut(i, j) = value;
            }
        }
    }
    Ok(result.trimmed())
}

impl ProjectedJacobian {
    pub(crate) fn new(
        homogeneous: &[MPoly; 4],
        projection: &[Vector3; 2],
        domain: Domain,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let (projection, projection_scale) = integral_projection(projection, dec)?;
        let local = [
            Dense::local(&homogeneous[0], &domain, dec)?,
            Dense::local(&homogeneous[1], &domain, dec)?,
            Dense::local(&homogeneous[2], &domain, dec)?,
            Dense::local(&homogeneous[3], &domain, dec)?,
        ];
        let weight = bernstein(&local[3], dec)?;
        if !weight.all_positive() {
            return Err(Refusal::new(
                "shapes.rational-atlas-polynomial-weights",
                "the exact span weight has no strictly positive Bernstein hull",
            ));
        }
        let mut project = |vector: &Vector3| {
            let a = direction(vector)?;
            combine(
                &[(&local[0], a[0]), (&local[1], a[1]), (&local[2], a[2])],
                dec,
            )
        };
        let projected = [project(&projection[0])?, project(&projection[1])?];
        let dw = [local[3].derivative(0, dec)?, local[3].derivative(1, dec)?];
        let derivative = |coordinate: &Dense, axis: usize, dec: &mut Decider| {
            let dn = coordinate.derivative(axis, dec)?;
            let left = local[3].mul(&dn, dec)?;
            let right = coordinate.mul(&dw[axis], dec)?;
            left.sub(&right, dec)
        };
        let n11 = derivative(&projected[0], 0, dec)?;
        let n22 = derivative(&projected[1], 1, dec)?;
        let n12 = derivative(&projected[0], 1, dec)?;
        let n21 = derivative(&projected[1], 0, dec)?;
        let one = Rational::one();
        let off = combine(&[(&n12, &one), (&n21, &one)], dec)?;
        let diagonal = [bernstein(&n11, dec)?, bernstein(&n22, dec)?];
        let off_sum = bernstein(&off, dec)?;
        let widths = [
            sub(&domain[0][1], &domain[0][0], dec)?,
            sub(&domain[1][1], &domain[1][0], dec)?,
        ];
        Ok(Self {
            domain,
            widths,
            projection_scale,
            weight,
            diagonal,
            off_sum,
        })
    }

    pub(crate) fn margin(&self, cell: &Domain, dec: &mut Decider) -> Result<Rational, Refusal> {
        let cuts = restriction(&self.domain, &self.widths, cell, dec)?;
        let weight = restricted_range(&self.weight, &cuts, dec)?;
        let a = restricted_range(&self.diagonal[0], &cuts, dec)?;
        let b = restricted_range(&self.diagonal[1], &cuts, dec)?;
        if weight[0].sign() != Sign::Positive
            || a[0].sign() != Sign::Positive
            || b[0].sign() != Sign::Positive
        {
            return Ok(Rational::zero());
        }
        let off = restricted_range(&self.off_sum, &cuts, dec)?;
        let ends = off.each_ref().map(|value| value.abs());
        let upper = if less(&ends[0], &ends[1], dec)? {
            &ends[1]
        } else {
            &ends[0]
        };
        let squared = mul(upper, upper, dec)?;
        let quarter = quotient(&squared, &Rational::from_i64(4), dec)?;
        let along = mul(&a[0], &b[0], dec)?;
        let determinant = sub(&along, &quarter, dec)?;
        if determinant.sign() != Sign::Positive {
            return Ok(Rational::zero());
        }
        let trace = add(&a[1], &b[1], dec)?;
        let weight_squared = mul(&weight[1], &weight[1], dec)?;
        let denominator = mul(&trace, &weight_squared, dec)?;
        let denominator = mul(&denominator, &self.projection_scale, dec)?;
        quotient(&determinant, &denominator, dec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BSplineSurface, Point3, RationalPatch};
    use catallax_numbers::with_operation_budget;

    fn q(n: i64, d: i64) -> Rational {
        Rational::from_ratio(n, d).unwrap()
    }

    #[test]
    fn weighted_bending_chart_keeps_translated_geometry_and_nonunit_parameters() {
        let (result, work) = with_operation_budget(20_000_000, || -> Result<(), Refusal> {
            let mut dec = Decider::default();
            let mut points = Vec::new();
            let mut weights = Vec::new();
            for u in 0..3 {
                for v in 0..3 {
                    let z: i64 = 300 + [0, 1, 0][u] + if v == 2 { 1 } else { 0 };
                    points.push(Point3::new(
                        (100 + u as i64).into(),
                        Real::from(q(-400 + v as i64, 2)),
                        z.into(),
                    ));
                    weights.push(Real::from(if u == 1 { q(1, 2) } else { q(1, 1) }));
                }
            }
            let carrier = BSplineSurface::new(
                [2, 2],
                [
                    [2, 2, 2, 6, 6, 6].map(Rational::from_i64).to_vec(),
                    [-3, -3, -3, 5, 5, 5].map(Rational::from_i64).to_vec(),
                ],
                [3, 3],
                points,
                weights,
                &mut dec,
            )?;
            let patch = RationalPatch::new(carrier, &mut dec)?;
            let point = patch.point_at(&[Real::from(3), Real::from(3)], &mut dec)?;
            assert_eq!(point.x.as_rational(), Some(q(1305, 13)));
            assert_eq!(point.y.as_rational(), Some(q(-797, 4)));
            assert_eq!(point.z.as_rational(), Some(q(62565, 208)));
            assert_eq!(patch.embedding_lower_bound().sign(), Sign::Positive);
            Ok(())
        });
        assert!(!work.exhausted);
        result.unwrap();
    }

    #[test]
    fn homogeneous_projection_normalization_keeps_a_valid_geometric_margin() {
        let (result, work) = with_operation_budget(20_000_000, || -> Result<(), Refusal> {
            let mut dec = Decider::default();
            let points = (0..3)
                .flat_map(|u| {
                    (0..3).map(move |v| {
                        Point3::new(
                            q(3 * u, 2).into(),
                            q(5 * v, 2).into(),
                            (&q(if u == 1 { -1 } else { 1 }, 4)
                                + &q(if v == 1 { -1 } else { 1 }, 4))
                                .into(),
                        )
                    })
                })
                .collect();
            let knots = [0, 0, 0, 1, 1, 1].map(Rational::from_i64).to_vec();
            let carrier = BSplineSurface::new(
                [2, 2],
                [knots.clone(), knots],
                [3, 3],
                points,
                vec![Real::one(); 9],
                &mut dec,
            )?;
            let patch = RationalPatch::new(carrier, &mut dec)?;
            // X=3u,Y=5v,Z=(u-1/2)^2+(v-1/2)^2. Its centre-dual
            // projection is (X/3,Y/5), whose Jacobian is exactly I.
            assert_eq!(patch.embedding_lower_bound().sign(), Sign::Positive);
            assert!(patch.embedding_lower_bound() <= &Rational::one());
            Ok(())
        });
        assert!(!work.exhausted);
        result.unwrap();
    }

    #[test]
    fn a_positive_centre_frame_does_not_admit_an_interior_fold() {
        let (result, work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let points = [q(1, 4), q(1, 12), q(11, 12), q(3, 4)]
                .into_iter()
                .flat_map(|x| {
                    [0, 1].map(move |v| Point3::new(Real::from(x.clone()), v.into(), Real::zero()))
                })
                .collect();
            // X = u - 2(u-1/2)^3, Y = v. X_u(1/2)=1, but X_u(0)=-1/2.
            let carrier = BSplineSurface::new(
                [3, 1],
                [
                    [0, 0, 0, 0, 1, 1, 1, 1].map(Rational::from_i64).to_vec(),
                    [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
                ],
                [4, 2],
                points,
                vec![Real::one(); 8],
                &mut dec,
            )?;
            RationalPatch::new(carrier, &mut dec)
        });
        assert!(!work.exhausted);
        assert_eq!(result.unwrap_err().code, "shapes.rational-patch-unproven");
    }
}
