//! Complete finite preimages through two affine world coordinates.

use super::*;

pub(super) enum Preimage {
    Unsupported,
    Absent,
    Point([Real; 2]),
}

struct Frame {
    axes: [usize; 2],
    coordinates: [[Real; 3]; 2],
    weight: Real,
    determinant: Real,
}

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

fn difference(a: &Real, b: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
    arithmetic(a, &(-b), false, dec)
}

fn quotient(a: &Real, b: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
    if let (Some(a), Some(b)) = (a.rational_ref(), b.rational_ref()) {
        dec.reserve_exact_arithmetic(a, b, true)?;
    } else {
        dec.charge_work(1)?;
    }
    Ok(a / b)
}

// Constant, u and v coefficients. Unrecognised terms keep the general solver.
fn affine(poly: &RealPoly) -> Option<[Real; 3]> {
    let mut coefficients = std::array::from_fn(|_| Real::zero());
    for (monomial, coefficient) in poly.terms() {
        let slot = match monomial.exponents() {
            [0, 0] => 0,
            [1, 0] => 1,
            [0, 1] => 2,
            _ => return None,
        };
        coefficients[slot] = coefficient;
    }
    Some(coefficients)
}

fn frame(span: &Span, dec: &mut Decider) -> Result<Option<Frame>, Refusal> {
    let Some([weight, u, v]) = affine(&span.coordinates[3]) else {
        return Ok(None);
    };
    if !u.rational_ref().is_some_and(Rational::is_zero)
        || !v.rational_ref().is_some_and(Rational::is_zero)
        || !weight.rational_ref().is_some_and(|w| !w.is_zero())
    {
        return Ok(None);
    }
    let coordinates = span.coordinates[..3].iter().map(affine).collect::<Vec<_>>();
    for axes in [[0, 1], [0, 2], [1, 2]] {
        let (Some(a), Some(b)) = (&coordinates[axes[0]], &coordinates[axes[1]]) else {
            continue;
        };
        if [&a[1], &a[2], &b[1], &b[2]]
            .iter()
            .any(|value| value.rational_ref().is_none())
        {
            continue;
        }
        let ad = arithmetic(&a[1], &b[2], true, dec)?;
        let bc = arithmetic(&a[2], &b[1], true, dec)?;
        let determinant = difference(&ad, &bc, dec)?;
        if determinant.rational_ref().unwrap().is_zero() {
            continue;
        }
        return Ok(Some(Frame {
            axes,
            coordinates: [a.clone(), b.clone()],
            weight,
            determinant,
        }));
    }
    Ok(None)
}

fn candidate(frame: &Frame, point: &Point3, dec: &mut Decider) -> Result<[Real; 2], Refusal> {
    let position = [&point.x, &point.y, &point.z];
    let mut rhs = Vec::new();
    for (axis, coordinate) in frame.axes.iter().zip(&frame.coordinates) {
        let scaled = arithmetic(&frame.weight, position[*axis], true, dec)?;
        rhs.push(difference(&scaled, &coordinate[0], dec)?);
    }
    let [a, b] = &frame.coordinates;
    let du = arithmetic(&b[2], &rhs[0], true, dec)?;
    let bv = arithmetic(&a[2], &rhs[1], true, dec)?;
    let av = arithmetic(&a[1], &rhs[1], true, dec)?;
    let cu = arithmetic(&b[1], &rhs[0], true, dec)?;
    Ok([
        quotient(&difference(&du, &bv, dec)?, &frame.determinant, dec)?,
        quotient(&difference(&av, &cu, dec)?, &frame.determinant, dec)?,
    ])
}

fn value(poly: &RealPoly, at: &[Real; 2], dec: &mut Decider) -> Result<Real, Refusal> {
    let mut result = Real::zero();
    for (monomial, coefficient) in poly.terms() {
        let mut term = coefficient;
        for (axis, &power) in monomial.exponents().iter().enumerate() {
            for _ in 0..power {
                term = arithmetic(&term, &at[axis], true, dec)?;
            }
        }
        result = arithmetic(&result, &term, false, dec)?;
    }
    Ok(result)
}

impl RationalPatch {
    pub(crate) fn common_affine_projection(
        patches: &[&Self],
        dec: &mut Decider,
    ) -> Result<Option<[Vector3; 2]>, Refusal> {
        let mut reference: Option<([usize; 2], [[Real; 3]; 2], bool)> = None;
        for patch in patches {
            for span in &patch.data.spans {
                let Some(frame) = frame(span, dec)? else {
                    return Ok(None);
                };
                let mut coordinates = [
                    [Real::zero(), Real::zero(), Real::zero()],
                    [Real::zero(), Real::zero(), Real::zero()],
                ];
                for (row, source) in coordinates.iter_mut().zip(&frame.coordinates) {
                    for (target, coefficient) in row.iter_mut().zip(source) {
                        *target = quotient(coefficient, &frame.weight, dec)?;
                    }
                }
                if let Some((axes, expected, _)) = &reference {
                    if *axes != frame.axes {
                        return Ok(None);
                    }
                    for (actual, expected) in
                        coordinates.iter().flatten().zip(expected.iter().flatten())
                    {
                        if exact_sign(&difference(actual, expected, dec)?, dec)? != Sign::Zero {
                            return Ok(None);
                        }
                    }
                } else {
                    reference = Some((
                        frame.axes,
                        coordinates,
                        frame.determinant.rational_ref().unwrap().sign() == Sign::Positive,
                    ));
                }
            }
        }
        Ok(reference.map(|(axes, _, positive)| {
            let unit = |axis, sense| {
                let [x, y, z] =
                    std::array::from_fn(|i| Real::from_i64(if i == axis { sense } else { 0 }));
                Vector3::new(x, y, z)
            };
            [
                unit(axes[0], 1),
                unit(axes[1], if positive { 1 } else { -1 }),
            ]
        }))
    }
}

pub(super) fn parameters(
    patch: &RationalPatch,
    point: &Point3,
    dec: &mut Decider,
) -> Result<Preimage, Refusal> {
    let mut frames = Vec::new();
    for span in &patch.data.spans {
        let Some(frame) = frame(span, dec)? else {
            return Ok(Preimage::Unsupported);
        };
        frames.push(frame);
    }
    math::validate_point(point, dec)?;
    let mut found: Option<[Real; 2]> = None;
    for (span, frame) in patch.data.spans.iter().zip(frames) {
        let at = candidate(&frame, point, dec)?;
        let mut inside = true;
        for (parameter, bounds) in at.iter().zip(&span.domain) {
            let lower = difference(parameter, &Real::from_rational(bounds[0].clone()), dec)?;
            let upper = difference(&Real::from_rational(bounds[1].clone()), parameter, dec)?;
            inside &= exact_sign(&lower, dec)? != Sign::Negative
                && exact_sign(&upper, dec)? != Sign::Negative;
        }
        if !inside {
            continue;
        }
        // The nonsingular affine pair proves two coordinate equations and
        // uniqueness. The remaining original homogeneous equation decides
        // incidence; a projection-only preimage is never accepted.
        let remaining = (0..3).find(|axis| !frame.axes.contains(axis)).unwrap();
        let actual = value(&span.coordinates[remaining], &at, dec)?;
        let position = [&point.x, &point.y, &point.z];
        let expected = arithmetic(&frame.weight, position[remaining], true, dec)?;
        if exact_sign(&difference(&actual, &expected, dec)?, dec)? != Sign::Zero {
            continue;
        }
        if let Some(previous) = &found {
            for axis in 0..2 {
                if exact_sign(&difference(&at[axis], &previous[axis], dec)?, dec)? != Sign::Zero {
                    return Err(unproven(
                        "a checked embedding acquired distinct finite preimages",
                    ));
                }
            }
        } else {
            found = Some(at);
        }
    }
    Ok(match found {
        Some(at) => Preimage::Point(at),
        None => Preimage::Absent,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    fn q(n: i64, d: i64) -> Real {
        Rational::from_ratio(n, d).unwrap().into()
    }

    fn graph(coupled: bool, curved_x: bool, dec: &mut Decider) -> RationalPatch {
        let mut controls = Vec::new();
        for i in 0..3 {
            for j in 0..4 {
                let u = q(i, 2);
                let v = q(j, 3);
                let (x, z) = if coupled {
                    (
                        &(&Real::from_i64(2) * &u) + &v,
                        &u + &(&Real::from_i64(3) * &v),
                    )
                } else {
                    (u, v)
                };
                let x = if curved_x {
                    x + q(i64::from(i == 2), 64)
                } else {
                    x
                };
                let y = q(i64::from(i == 2) + i64::from(j == 3), 64);
                controls.push(Point3::new(x, y, z));
            }
        }
        let knots = |values: &[i64]| values.iter().map(|n| Rational::from_i64(*n)).collect();
        let source = BSplineSurface::new(
            [2, 3],
            [knots(&[0, 0, 0, 1, 1, 1]), knots(&[0, 0, 0, 0, 1, 1, 1, 1])],
            [3, 4],
            controls,
            vec![Real::from_i64(if coupled { 2 } else { 1 }); 12],
            dec,
        )
        .unwrap();
        RationalPatch::new(source, dec).unwrap()
    }

    #[test]
    fn common_graph_projection_keeps_orientation_and_rejects_mismatched_maps() {
        let ((patch, other), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let patch = graph(false, false, &mut dec);
            let other = graph(true, false, &mut dec);
            let rotation = crate::Rotation::from_quaternion(
                Rational::zero(),
                Rational::one(),
                Rational::zero(),
                Rational::zero(),
            )
            .unwrap();
            let pose = crate::Pose::new(
                rotation,
                Vector3::new(Real::zero(), Real::zero(), Real::zero()),
                &mut dec,
            )
            .unwrap();
            let rotated = patch.transformed(&pose);
            for source in [&patch, &other, &rotated] {
                let [a, b] = RationalPatch::common_affine_projection(&[source, source], &mut dec)
                    .unwrap()
                    .unwrap();
                let [u, v] = source
                    .derivatives_at(&[q(1, 3), q(1, 4)], &mut dec)
                    .unwrap();
                assert_eq!(
                    exact_sign(&a.cross(&b).dot(&u.cross(&v)), &mut dec).unwrap(),
                    Sign::Positive
                );
            }
            assert!(
                RationalPatch::common_affine_projection(&[&patch, &other], &mut dec)
                    .unwrap()
                    .is_none()
            );
            let curved = graph(false, true, &mut dec);
            assert!(
                RationalPatch::common_affine_projection(&[&curved], &mut dec)
                    .unwrap()
                    .is_none()
            );
            (patch, other)
        });
        assert!(!work.exhausted, "{} units", work.used);
        let (answer, stopped) = with_operation_budget(0, || {
            RationalPatch::common_affine_projection(&[&patch, &other], &mut Decider::default())
        });
        assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
        assert!(stopped.exhausted);
    }

    #[test]
    fn exact_graph_preimages_keep_algebraic_coordinates_and_reject_off_carrier_points() {
        for coupled in [false, true] {
            let ((), work) = with_operation_budget(20_000_000, || {
                let mut dec = Decider::default();
                let patch = graph(coupled, false, &mut dec);
                let variable = RealPoly::var(1, 0);
                let equation = &variable.pow_u32(3) - &RealPoly::constant(1, Real::from_i64(2));
                let v = equation.real_roots(&mut dec).unwrap().remove(0) / Real::from_i64(2);
                let u = q(1, 3);
                let (x, z) = if coupled {
                    (
                        &(&Real::from_i64(2) * &u) + &v,
                        &u + &(&Real::from_i64(3) * &v),
                    )
                } else {
                    (u.clone(), v.clone())
                };
                let on = Point3::new(x.clone(), q(13, 2304), z.clone());
                let at = patch.parameters_at(&on, &mut dec).unwrap().unwrap();
                assert_eq!(exact_sign(&(&at[0] - &u), &mut dec).unwrap(), Sign::Zero);
                assert_eq!(exact_sign(&(&at[1] - &v), &mut dec).unwrap(), Sign::Zero);
                let off = Point3::new(x, q(49, 2304), z);
                assert!(patch.parameters_at(&off, &mut dec).unwrap().is_none());
            });
            assert!(!work.exhausted, "coupled={coupled}: {} units", work.used);
        }
    }

    #[test]
    fn finite_domain_and_zero_budget_remain_enforced() {
        let mut dec = Decider::default();
        let patch = graph(false, false, &mut dec);
        let outside = Point3::new(q(2, 1), q(4, 64), Real::zero());
        assert!(patch.parameters_at(&outside, &mut dec).unwrap().is_none());
        let on = Point3::new(q(1, 2), q(1, 256), Real::zero());
        let (answer, work) = with_operation_budget(0, || patch.parameters_at(&on, &mut dec));
        assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
        assert!(work.exhausted);
    }

    #[test]
    fn non_affine_coordinates_keep_the_existing_preimage_solver() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let patch = graph(false, true, &mut dec);
            let at = [q(1, 3), q(1, 2)];
            let point = patch.point_at(&at, &mut dec).unwrap();
            assert!(matches!(
                parameters(&patch, &point, &mut dec),
                Ok(Preimage::Unsupported)
            ));
            let actual = patch.parameters_at(&point, &mut dec).unwrap().unwrap();
            for axis in 0..2 {
                assert_eq!(
                    exact_sign(&(&actual[axis] - &at[axis]), &mut dec).unwrap(),
                    Sign::Zero
                );
            }
        });
        assert!(!work.exhausted, "{} units", work.used);
    }

    #[test]
    fn different_plane_supports_can_still_share_points() {
        let ((geom, spline, plane), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let patch = graph(false, false, &mut dec);
            let point = patch.point_at(&[q(1, 2), q(1, 2)], &mut dec).unwrap();
            let mut geom = crate::GeomStore::new();
            let spline = geom.add_surface(crate::Surface::BSpline(patch.into()));
            let plane = geom.add_surface(crate::Surface::Plane(crate::Plane {
                origin: Point3::new(Real::zero(), Real::zero(), q(1, 2)),
                normal: Vector3::new(Real::zero(), Real::zero(), Real::one()),
                ref_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
            }));
            for (a, b) in [(spline, plane), (plane, spline)] {
                assert!(matches!(
                    geom.same_surface(a, b, &mut dec),
                    Outcome::Proven(_)
                ));
            }
            let point = geom.add_point(crate::PointDef::Explicit(point));
            assert!(matches!(
                geom.side(plane, point, &mut dec),
                Outcome::ProvenCoincident(_)
            ));
            (geom, spline, plane)
        });
        assert!(!work.exhausted, "{} units", work.used);
        let (answer, stopped) = with_operation_budget(0, || {
            geom.same_surface(spline, plane, &mut Decider::default())
        });
        assert!(
            matches!(answer, Outcome::Refused(error) if error.code == "numbers.operation-budget")
        );
        assert!(stopped.exhausted);
    }

    #[test]
    fn a_zero_plane_polynomial_keeps_the_unproved_equality_refusal() {
        let mut dec = Decider::default();
        let knots = [0, 0, 1, 1].map(Rational::from_i64).to_vec();
        let source = BSplineSurface::new(
            [1, 1],
            [knots.clone(), knots],
            [2, 2],
            vec![
                Point3::new(q(0, 1), q(0, 1), q(0, 1)),
                Point3::new(q(0, 1), q(0, 1), q(1, 1)),
                Point3::new(q(1, 1), q(0, 1), q(0, 1)),
                Point3::new(q(1, 1), q(0, 1), q(1, 1)),
            ],
            vec![Real::one(); 4],
            &mut dec,
        )
        .unwrap();
        let patch = RationalPatch::new(source, &mut dec).unwrap();
        let plane = crate::Plane {
            origin: Point3::new(Real::zero(), Real::zero(), Real::zero()),
            normal: Vector3::new(Real::zero(), Real::one(), Real::zero()),
            ref_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
        };
        assert!(!patch.proves_off_plane(&plane, &mut dec).unwrap());
        let mut geom = crate::GeomStore::new();
        let a = geom.add_surface(crate::Surface::BSpline(patch.into()));
        let b = geom.add_surface(crate::Surface::Plane(plane));
        assert!(
            matches!(geom.same_surface(a, b, &mut dec), Outcome::Refused(error) if error.code == "shapes.parametric-query")
        );
    }
}
