//! Exact outgoing normal acceleration of a finite native spline.

use crate::spline_arithmetic::{arithmetic, derivative_value, difference, quotient};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{BSplineCurve, Curve, CurveId, GeomStore, Point3, Vector3};

fn missing(detail: &str) -> Refusal {
    Refusal::new("intersect.spline-germ", detail)
}

fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    crate::spline_arithmetic::sign(value, dec, "intersect.spline-germ")
}

fn components(vector: &Vector3) -> [&Real; 3] {
    [&vector.x, &vector.y, &vector.z]
}

fn dot(a: &Vector3, b: &Vector3, dec: &mut Decider) -> Result<Real, Refusal> {
    let mut sum = Real::zero();
    for (a, b) in components(a).into_iter().zip(components(b)) {
        let term = arithmetic(a, b, true, dec)?;
        sum = arithmetic(&sum, &term, false, dec)?;
    }
    Ok(sum)
}

fn cross(a: &Vector3, b: &Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
    let a = components(a);
    let b = components(b);
    let mut result = Vec::new();
    for (i, j) in [(1, 2), (2, 0), (0, 1)] {
        let first = arithmetic(a[i], b[j], true, dec)?;
        let second = arithmetic(a[j], b[i], true, dec)?;
        result.push(difference(&first, &second, dec)?);
    }
    let [x, y, z] = result.try_into().expect("three cross-product components");
    Ok(Vector3::new(x, y, z))
}

fn scaled(vector: &Vector3, scale: &Real, dec: &mut Decider) -> Result<Vector3, Refusal> {
    Ok(Vector3::new(
        arithmetic(&vector.x, scale, true, dec)?,
        arithmetic(&vector.y, scale, true, dec)?,
        arithmetic(&vector.z, scale, true, dec)?,
    ))
}

fn minus(a: &Vector3, b: &Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
    Ok(Vector3::new(
        difference(&a.x, &b.x, dec)?,
        difference(&a.y, &b.y, dec)?,
        difference(&a.z, &b.z, dec)?,
    ))
}

fn same_position(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    for (a, b) in [(&a.x, &b.x), (&a.y, &b.y), (&a.z, &b.z)] {
        if sign(&difference(a, b, dec)?, dec)? != Sign::Zero {
            return Ok(false);
        }
    }
    Ok(true)
}

fn exposed_endpoint(
    curve: &BSplineCurve,
    point: &Point3,
    dec: &mut Decider,
) -> Result<Option<Real>, Refusal> {
    if curve.degree() == 0 {
        return Ok(None);
    }
    let controls = curve.controls().collect::<Vec<_>>();
    for end in [0, controls.len() - 1] {
        if !same_position(point, controls[end], dec)? {
            continue;
        }
        let opposite = if end == 0 { controls.len() - 1 } else { 0 };
        let seam = same_position(point, controls[opposite], dec)?;
        let endpoint = [&point.x, &point.y, &point.z];
        for (axis, endpoint) in endpoint.into_iter().enumerate() {
            let mut side = None;
            let mut exposed = true;
            for (index, control) in controls.iter().enumerate() {
                if index == end || (seam && index == opposite) {
                    continue;
                }
                let coordinate = [&control.x, &control.y, &control.z][axis];
                let direction = sign(&difference(coordinate, endpoint, dec)?, dec)?;
                if direction == Sign::Zero || side.is_some_and(|side| side != direction) {
                    exposed = false;
                    break;
                }
                side = Some(direction);
            }
            if exposed && side.is_some() {
                // Positive weights preserve rational B-spline convex combinations.
                // At every interior parameter, positive total basis weight lies
                // outside the matched endpoint controls. Strict exposure of all
                // those controls excludes every interior point preimage.
                let parameter = if end == 0 { 0 } else { 1 };
                return Ok(Some(Real::from_rational(curve.domain()[parameter].clone())));
            }
        }
    }
    Ok(None)
}

fn jet(
    curve: &BSplineCurve,
    span: usize,
    parameter: &Real,
    dec: &mut Decider,
) -> Result<[Vector3; 2], Refusal> {
    let coordinates = curve.query_span_polynomials(span, dec)?;
    let mut homogeneous = [Vec::new(), Vec::new(), Vec::new()];
    for (order, values) in homogeneous.iter_mut().enumerate() {
        for coordinate in &coordinates {
            values.push(derivative_value(coordinate, parameter, order, dec)?);
        }
    }
    let [position, first, second] = &homogeneous;
    let mut velocity = Vec::new();
    let mut acceleration = Vec::new();
    for axis in 0..3 {
        // H = W C, H' = W' C + W C',
        // H'' = W'' C + 2 W' C' + W C''. Positive weights prove W nonzero.
        let p = quotient(&position[axis], &position[3], dec)?;
        let correction = arithmetic(&p, &first[3], true, dec)?;
        let v = quotient(
            &difference(&first[axis], &correction, dec)?,
            &position[3],
            dec,
        )?;
        let correction = arithmetic(&p, &second[3], true, dec)?;
        let numerator = difference(&second[axis], &correction, dec)?;
        let twice = arithmetic(&Real::from_i64(2), &v, true, dec)?;
        let correction = arithmetic(&twice, &first[3], true, dec)?;
        let a = quotient(
            &difference(&numerator, &correction, dec)?,
            &position[3],
            dec,
        )?;
        velocity.push(v);
        acceleration.push(a);
    }
    let vector = |values: Vec<Real>| {
        let [x, y, z] = values.try_into().expect("three Cartesian derivatives");
        Vector3::new(x, y, z)
    };
    Ok([vector(velocity), vector(acceleration)])
}

/// Exact normal acceleration with the supplied nonzero outgoing tangent as
/// velocity. A unique finite preimage and a regular knot side are required.
/// Closed seams consider both clamped ends. An unproved branch selector or
/// regular reparameterization refuses with intersect.spline-germ.
/// Other curve kinds return None and retain their existing query route.
pub fn native_spline_normal_acceleration(
    geom: &GeomStore,
    curve: CurveId,
    position: &Point3,
    tangent: &Vector3,
    dec: &mut Decider,
) -> Result<Option<Vector3>, Refusal> {
    let Curve::BSpline(curve) = geom.curve(curve) else {
        return Ok(None);
    };
    let tangent_speed = dot(tangent, tangent, dec)?;
    if sign(&tangent_speed, dec)? != Sign::Positive {
        return Err(missing(
            "normal acceleration needs a nonzero outgoing tangent",
        ));
    }
    let at = match exposed_endpoint(curve, position, dec)? {
        Some(at) => at,
        None => curve
            .unique_parameter_of(position, dec)?
            .ok_or_else(|| {
                missing("normal acceleration needs an exact point on the finite spline")
            })?
            .value()
            .clone(),
    };
    let mut parameters = vec![at.clone()];
    let [lo, hi] = curve.domain();
    let at_lo = sign(
        &difference(&at, &Real::from_rational(lo.clone()), dec)?,
        dec,
    )? == Sign::Zero;
    let at_hi = sign(
        &difference(&Real::from_rational(hi.clone()), &at, dec)?,
        dec,
    )? == Sign::Zero;
    if at_lo || at_hi {
        // Clamping makes these controls actual curve endpoints. The unique
        // preimage query has already excluded unrelated interior preimages.
        let other = if at_lo {
            curve.controls().last().expect("finite spline endpoint")
        } else {
            curve.controls().next().expect("finite spline endpoint")
        };
        if same_position(position, other, dec)? {
            parameters.push(Real::from_rational(if at_lo {
                hi.clone()
            } else {
                lo.clone()
            }));
        }
    }
    let mut found: Option<Vector3> = None;
    for parameter in parameters {
        for (span, [lo, hi]) in curve.spans() {
            let lower = sign(
                &difference(&parameter, &Real::from_rational(lo.clone()), dec)?,
                dec,
            )?;
            let upper = sign(
                &difference(&Real::from_rational(hi.clone()), &parameter, dec)?,
                dec,
            )?;
            if lower == Sign::Negative || upper == Sign::Negative {
                continue;
            }
            let [first, second] = jet(curve, span, &parameter, dec)?;
            let speed = dot(&first, &first, dec)?;
            if sign(&speed, dec)? != Sign::Positive {
                return Err(missing(
                    "a zero knot-side speed needs a regular reparameterization proof",
                ));
            }
            let transverse = cross(&first, tangent, dec)?;
            if sign(&dot(&transverse, &transverse, dec)?, dec)? != Sign::Zero {
                continue;
            }
            let direction = sign(&dot(&first, tangent, dec)?, dec)?;
            if direction == Sign::Zero {
                return Err(missing(
                    "the outgoing tangent needs an exact knot-side orientation",
                ));
            }
            if (lower == Sign::Zero && direction == Sign::Negative)
                || (upper == Sign::Zero && direction == Sign::Positive)
            {
                continue;
            }
            let projection = quotient(&dot(&first, &second, dec)?, &speed, dec)?;
            let parallel = scaled(&first, &projection, dec)?;
            let normal = minus(&second, &parallel, dec)?;
            let scale = quotient(&tangent_speed, &speed, dec)?;
            let acceleration = scaled(&normal, &scale, dec)?;
            if let Some(previous) = &found {
                let delta = minus(previous, &acceleration, dec)?;
                for component in components(&delta) {
                    if sign(component, dec)? != Sign::Zero {
                        return Err(missing(
                            "the tangent selects unequal knot-side accelerations; an explicit branch selector is missing",
                        ));
                    }
                }
            } else {
                found = Some(acceleration);
            }
        }
    }
    found.map(Some).ok_or_else(|| {
        missing("the tangent does not select a regular outgoing side of the finite spline")
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Rational, with_operation_budget};

    fn q(n: i64, d: i64) -> Real {
        Rational::from_ratio(n, d).unwrap().into()
    }

    fn point(x: i64, y: i64) -> Point3 {
        Point3::new(Real::from_i64(x), Real::from_i64(y), Real::zero())
    }

    fn native(
        controls: Vec<Point3>,
        weights: Vec<Real>,
        knots: &[i64],
        degree: usize,
        dec: &mut Decider,
    ) -> BSplineCurve {
        BSplineCurve::new(
            degree,
            knots.iter().copied().map(Rational::from_i64).collect(),
            controls,
            weights,
            dec,
        )
        .unwrap()
    }

    fn parabola(dec: &mut Decider) -> BSplineCurve {
        native(
            vec![
                point(0, 0),
                Point3::new(q(1, 2), Real::zero(), Real::zero()),
                point(1, 1),
            ],
            vec![Real::one(); 3],
            &[0, 0, 0, 1, 1, 1],
            2,
            dec,
        )
    }

    fn equal(a: &Vector3, b: &Vector3, dec: &mut Decider) {
        let delta = minus(a, b, dec).unwrap();
        for value in components(&delta) {
            assert_eq!(sign(value, dec).unwrap(), Sign::Zero);
        }
    }

    #[test]
    fn polynomial_normal_acceleration_uses_the_requested_speed_and_direction() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let curve = parabola(&mut dec);
            let position = curve.point_at(&q(1, 2), &mut dec).unwrap();
            let mut geom = GeomStore::new();
            let curve = geom.add_curve(Curve::BSpline(curve));
            for speed in [1, -1, 2] {
                let tangent =
                    Vector3::new(Real::from_i64(speed), Real::from_i64(speed), Real::zero());
                let acceleration =
                    native_spline_normal_acceleration(&geom, curve, &position, &tangent, &mut dec)
                        .unwrap()
                        .unwrap();
                let squared = speed * speed;
                equal(
                    &acceleration,
                    &Vector3::new(
                        Real::from_i64(-squared),
                        Real::from_i64(squared),
                        Real::zero(),
                    ),
                    &mut dec,
                );
                assert_eq!(
                    sign(&dot(&acceleration, &tangent, &mut dec).unwrap(), &mut dec).unwrap(),
                    Sign::Zero
                );
            }
        });
        assert!(!work.exhausted, "{} units", work.used);
    }

    #[test]
    fn rational_circle_acceleration_includes_weight_derivatives() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let curve = native(
                vec![
                    point(1, 0),
                    Point3::new(Real::one(), q(1, 2), Real::zero()),
                    Point3::new(q(3, 5), q(4, 5), Real::zero()),
                ],
                vec![Real::one(), Real::one(), q(5, 4)],
                &[0, 0, 0, 1, 1, 1],
                2,
                &mut dec,
            );
            let mut geom = GeomStore::new();
            let id = geom.add_curve(Curve::BSpline(curve.clone()));
            for (parameter, sense) in [(q(0, 1), 2), (q(1, 2), 2), (q(1, 1), -2)] {
                let position = curve.point_at(&parameter, &mut dec).unwrap();
                let tangent = Vector3::new(
                    -&position.y * Real::from_i64(sense),
                    &position.x * &Real::from_i64(sense),
                    Real::zero(),
                );
                let acceleration =
                    native_spline_normal_acceleration(&geom, id, &position, &tangent, &mut dec)
                        .unwrap()
                        .unwrap();
                equal(
                    &acceleration,
                    &Vector3::new(
                        &position.x * &Real::from_i64(-4),
                        &position.y * &Real::from_i64(-4),
                        Real::zero(),
                    ),
                    &mut dec,
                );
            }
        });
        assert!(!work.exhausted, "{} units", work.used);
    }

    #[test]
    fn knot_sides_keep_their_different_second_derivatives() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let curve = native(
                vec![
                    point(0, 0),
                    Point3::new(q(1, 2), Real::zero(), Real::zero()),
                    point(1, 0),
                    Point3::new(q(3, 2), Real::zero(), Real::zero()),
                    point(2, 1),
                ],
                vec![Real::one(); 5],
                &[0, 0, 0, 1, 1, 2, 2, 2],
                2,
                &mut dec,
            );
            let mut geom = GeomStore::new();
            let curve = geom.add_curve(Curve::BSpline(curve));
            for (direction, expected) in [(1, 2), (-1, 0)] {
                let tangent = Vector3::new(Real::from_i64(direction), Real::zero(), Real::zero());
                let acceleration = native_spline_normal_acceleration(
                    &geom,
                    curve,
                    &point(1, 0),
                    &tangent,
                    &mut dec,
                )
                .unwrap()
                .unwrap();
                equal(
                    &acceleration,
                    &Vector3::new(Real::zero(), Real::from_i64(expected), Real::zero()),
                    &mut dec,
                );
            }
        });
        assert!(!work.exhausted, "{} units", work.used);
    }

    #[test]
    fn closed_seam_selects_the_finite_side_from_the_outgoing_tangent() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let origin = point(0, 0);
            let curve = native(
                vec![
                    origin.clone(),
                    point(1, 0),
                    point(1, 1),
                    point(0, 1),
                    origin,
                ],
                vec![Real::one(); 5],
                &[0, 0, 1, 2, 3, 4, 4],
                1,
                &mut dec,
            );
            let mut geom = GeomStore::new();
            let curve = geom.add_curve(Curve::BSpline(curve));
            for tangent in [
                Vector3::new(Real::one(), Real::zero(), Real::zero()),
                Vector3::new(Real::zero(), Real::one(), Real::zero()),
            ] {
                let acceleration = native_spline_normal_acceleration(
                    &geom,
                    curve,
                    &point(0, 0),
                    &tangent,
                    &mut dec,
                )
                .unwrap()
                .unwrap();
                equal(
                    &acceleration,
                    &Vector3::new(Real::zero(), Real::zero(), Real::zero()),
                    &mut dec,
                );
            }
        });
        assert!(!work.exhausted, "{} units", work.used);
    }

    #[test]
    fn unproved_directions_and_zero_speed_refuse_and_zero_budget_stops() {
        let ((geom, curve, position), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let curve = parabola(&mut dec);
            let position = curve.point_at(&q(1, 2), &mut dec).unwrap();
            let mut geom = GeomStore::new();
            let curve = geom.add_curve(Curve::BSpline(curve));
            let tangent = Vector3::new(Real::zero(), Real::zero(), Real::one());
            assert_eq!(
                native_spline_normal_acceleration(&geom, curve, &position, &tangent, &mut dec)
                    .unwrap_err()
                    .code,
                "intersect.spline-germ"
            );
            let singular = native(
                vec![point(0, 0), point(0, 0), point(1, 0)],
                vec![Real::one(); 3],
                &[0, 0, 0, 1, 1, 1],
                2,
                &mut dec,
            );
            let singular = geom.add_curve(Curve::BSpline(singular));
            let tangent = Vector3::new(Real::one(), Real::zero(), Real::zero());
            assert_eq!(
                native_spline_normal_acceleration(
                    &geom,
                    singular,
                    &point(0, 0),
                    &tangent,
                    &mut dec
                )
                .unwrap_err()
                .code,
                "intersect.spline-germ"
            );
            (geom, curve, position)
        });
        assert!(!work.exhausted, "{} units", work.used);
        let (answer, stopped) = with_operation_budget(0, || {
            native_spline_normal_acceleration(
                &geom,
                curve,
                &position,
                &Vector3::new(Real::one(), Real::one(), Real::zero()),
                &mut Decider::default(),
            )
        });
        assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
        assert!(stopped.exhausted);
    }

    #[test]
    fn exposed_endpoints_keep_positive_algebraic_weights_and_both_closed_sides() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let curve = native(
                vec![
                    point(0, 0),
                    Point3::new(q(1, 2), Real::zero(), Real::zero()),
                    point(1, 1),
                ],
                vec![Real::one(), Real::from_i64(2).sqrt(), Real::one()],
                &[0, 0, 0, 1, 1, 1],
                2,
                &mut dec,
            );
            assert_eq!(
                exposed_endpoint(&curve, &point(0, 0), &mut dec)
                    .unwrap()
                    .unwrap()
                    .rational_ref(),
                Some(&Rational::zero())
            );
            let mut geom = GeomStore::new();
            let id = geom.add_curve(Curve::BSpline(curve));
            let acceleration = native_spline_normal_acceleration(
                &geom,
                id,
                &point(0, 0),
                &Vector3::new(Real::one(), Real::zero(), Real::zero()),
                &mut dec,
            )
            .unwrap()
            .unwrap();
            equal(
                &acceleration,
                &Vector3::new(Real::zero(), Real::one(), Real::zero()),
                &mut dec,
            );

            let origin = point(0, 0);
            let curve = native(
                vec![origin.clone(), point(1, 1), point(2, 1), origin.clone()],
                vec![Real::one(); 4],
                &[0, 0, 1, 2, 3, 3],
                1,
                &mut dec,
            );
            assert!(
                exposed_endpoint(&curve, &origin, &mut dec)
                    .unwrap()
                    .is_some()
            );
            let id = geom.add_curve(Curve::BSpline(curve));
            for (x, y) in [(1, 1), (2, 1)] {
                let acceleration = native_spline_normal_acceleration(
                    &geom,
                    id,
                    &origin,
                    &Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::zero()),
                    &mut dec,
                )
                .unwrap()
                .unwrap();
                equal(
                    &acceleration,
                    &Vector3::new(Real::zero(), Real::zero(), Real::zero()),
                    &mut dec,
                );
            }
        });
        assert!(!work.exhausted, "{} units", work.used);
    }

    #[test]
    fn endpoint_hulls_do_not_hide_an_interior_preimage_or_a_constant_span() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let origin = point(0, 0);
            let folded = native(
                vec![
                    origin.clone(),
                    Point3::new(q(-1, 6), Real::zero(), Real::zero()),
                    Point3::new(q(1, 6), Real::zero(), Real::zero()),
                    origin.clone(),
                ],
                vec![Real::one(); 4],
                &[0, 0, 0, 0, 1, 1, 1, 1],
                3,
                &mut dec,
            );
            let constant = native(
                vec![origin.clone(), origin.clone(), point(1, 0)],
                vec![Real::one(); 3],
                &[0, 0, 1, 2, 2],
                1,
                &mut dec,
            );
            let mut geom = GeomStore::new();
            for (curve, expected) in [
                (folded, "shapes.spline-ambiguous-parameter"),
                (constant, "shapes.spline-inverse"),
            ] {
                assert!(
                    exposed_endpoint(&curve, &origin, &mut dec)
                        .unwrap()
                        .is_none()
                );
                let id = geom.add_curve(Curve::BSpline(curve));
                assert_eq!(
                    native_spline_normal_acceleration(
                        &geom,
                        id,
                        &origin,
                        &Vector3::new(Real::one(), Real::zero(), Real::zero()),
                        &mut dec
                    )
                    .unwrap_err()
                    .code,
                    expected
                );
            }
        });
        assert!(!work.exhausted, "{} units", work.used);
    }
}
