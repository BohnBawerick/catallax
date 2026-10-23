#![allow(missing_docs)]

use catallax_features::fit_elliptic_offset;
use catallax_numbers::{Ball, Decider, Rational, Real, with_operation_budget};
use catallax_shapes::{EllipticCylinder, EllipticOffsetDefinition, Point3, Vector3};

fn r(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn definition(distance: &str) -> EllipticOffsetDefinition {
    let mut dec = Decider::default();
    let original = EllipticCylinder::new(
        Point3::origin(),
        Vector3::new(Real::zero(), Real::zero(), Real::one()),
        Vector3::new(Real::one(), Real::zero(), Real::zero()),
        Real::from_i64(4),
        Real::one(),
        &mut dec,
    )
    .unwrap();
    EllipticOffsetDefinition::new(original, r(distance).into(), &mut dec).unwrap()
}
fn rational_source_definition(dec: &mut Decider) -> EllipticOffsetDefinition {
    let original = EllipticCylinder::new(
        Point3::origin(),
        Vector3::new(3.into(), 0.into(), 4.into()),
        Vector3::new(0.into(), 1.into(), 0.into()),
        4.into(),
        r("64/25").into(),
        dec,
    )
    .unwrap();
    EllipticOffsetDefinition::new(original, r("-1/4").into(), dec).unwrap()
}
fn ball(value: &Rational) -> Ball {
    Ball::from_rational(value, 128)
}

#[test]
fn closed_cubic_offset_has_shared_tangents_and_bounded_interior_error() {
    let definition = definition("-1/4");
    let mut dec = Decider::default();
    let target = r("0.0001");
    let fit = fit_elliptic_offset(
        &definition,
        [r("1/4").into(), r("15/4").into()],
        &target,
        &mut dec,
    )
    .unwrap();
    assert!(fit.position_error.upper <= target);
    let measured = catallax_field::measurement::measure_elliptic_offset_sweep(
        &fit.surface,
        &definition,
        [&fit.axial[0], &fit.axial[1]],
        &catallax_field::measurement::MeasurementOptions {
            width: target.clone(),
            max_cells: 1,
            path_segments: 1,
        },
        &mut dec,
    )
    .unwrap();
    assert!(measured.converged);
    assert_eq!(measured.upper, fit.position_error.upper);
    assert_eq!(measured.spans, fit.spans.len());
    let zero = ball(&r("0"));
    let first = fit.surface.jet(&zero, &zero, &dec).unwrap();
    let last = fit.surface.jet(&ball(&r("4")), &zero, &dec).unwrap();
    for axis in 0..3 {
        assert!((&first.value[axis] - &last.value[axis]).contains_zero());
        assert!((&first.first[0][axis] - &last.first[0][axis]).contains_zero());
        assert!(&first.first[0][axis].hi() - &first.first[0][axis].lo() < r("1e-20"));
    }
    let chart = definition.chart();
    for (i, span) in fit.spans.iter().enumerate() {
        // At a shared internal knot, both one-sided derivative enclosures must
        // be the same tight value; a kink would widen the hull.
        let knot = Rational::from_i64(span.quarter as i64) + span.parameter[0].clone();
        let joint = fit.surface.jet(&ball(&knot), &zero, &dec).unwrap();
        for derivative in &joint.first[0] {
            assert!(&derivative.hi() - &derivative.lo() < r("1e-20"));
        }
        if i % 5 != 0 {
            continue;
        }
        let middle = (&span.parameter[0] + &span.parameter[1]).mul_pow2(-1);
        let (exact, _) = chart
            .point_and_tangent(span.quarter, &middle, &r("1/4").into(), &mut dec)
            .unwrap();
        let u = Rational::from_i64(span.quarter as i64) + middle;
        let candidate = fit.surface.jet(&ball(&u), &zero, &dec).unwrap();
        let mut actual_upper = Rational::zero();
        for (actual, expected) in candidate.value.iter().zip([exact.x, exact.y, exact.z]) {
            let error = actual - &dec.enclose_at(&expected, 128).unwrap();
            actual_upper = actual_upper + error.lo().abs().max(error.hi().abs());
        }
        assert!(actual_upper <= span.error.upper);
    }
}

#[test]
fn unsupported_normal_branch_and_exhausted_work_do_not_return_a_fit() {
    let mut dec = Decider::default();
    let result = fit_elliptic_offset(
        &definition("-3/4"),
        [Real::zero(), Real::one()],
        &r("0.0001"),
        &mut dec,
    );
    assert_eq!(result.unwrap_err().code, "features.offset-fit");
    let definition = definition("-1/4");
    let (result, work) = with_operation_budget(0, || {
        fit_elliptic_offset(
            &definition,
            [Real::zero(), Real::one()],
            &r("0.0001"),
            &mut dec,
        )
    });
    assert!(result.is_err());
    assert!(work.exhausted);
}

#[test]
fn independently_measured_sweep_rejects_wrong_side_detached_column_and_forged_seam() {
    use catallax_field::measurement::{MeasurementOptions, measure_elliptic_offset_sweep};
    use catallax_shapes::BSplineSurface;
    let source = definition("-1/4");
    let mut dec = Decider::default();
    let fit =
        fit_elliptic_offset(&source, [Real::zero(), Real::one()], &r("0.0001"), &mut dec).unwrap();
    let options = MeasurementOptions {
        width: r("0.0001"),
        ..Default::default()
    };
    let result = measure_elliptic_offset_sweep(
        &fit.surface,
        &definition("1/4"),
        [&fit.axial[0], &fit.axial[1]],
        &options,
        &mut dec,
    );
    assert_eq!(result.unwrap_err().code, "field.offset-spline");
    for index in [1, fit.surface.controls().len() - 2] {
        let mut controls = fit.surface.controls().to_vec();
        controls[index].x = &controls[index].x + &Real::from_rational(r("0.00001"));
        let forged = BSplineSurface::new(
            fit.surface.degrees(),
            fit.surface.knots().map(<[Rational]>::to_vec),
            fit.surface.control_count(),
            controls,
            fit.surface.weights().to_vec(),
            &mut dec,
        )
        .unwrap();
        let result = measure_elliptic_offset_sweep(
            &forged,
            &source,
            [&fit.axial[0], &fit.axial[1]],
            &options,
            &mut dec,
        );
        assert_eq!(result.unwrap_err().code, "field.offset-spline");
    }
    let shifted_knots = fit
        .surface
        .knots()
        .map(|knots| knots.iter().map(|x| x + &Rational::one()).collect());
    let cropped = BSplineSurface::new(
        fit.surface.degrees(),
        shifted_knots,
        fit.surface.control_count(),
        fit.surface.controls().to_vec(),
        fit.surface.weights().to_vec(),
        &mut dec,
    )
    .unwrap();
    let result = measure_elliptic_offset_sweep(
        &cropped,
        &source,
        [&fit.axial[0], &fit.axial[1]],
        &options,
        &mut dec,
    );
    assert_eq!(result.unwrap_err().code, "field.offset-spline");
}

#[test]
fn fitted_algebraic_boundary_recovers_exact_non_dyadic_parameter_without_elimination() {
    use catallax_shapes::{BSplineCurve, SplineBoundary};
    let mut dec = Decider::default();
    let definition = rational_source_definition(&mut dec);
    let (result, work) = with_operation_budget(100_000_000, || {
        let fit = fit_elliptic_offset(
            &definition,
            [r("9/4").into(), r("23/4").into()],
            &r("0.0001"),
            &mut dec,
        )?;
        let edge = BSplineCurve::boundary(&fit.surface, SplineBoundary::VMin);
        let expected = Real::from_rational(r("1/3"));
        let point = edge.point_at(&expected, &mut dec)?;
        let inverse = edge
            .parameter_of(&point, &mut dec)?
            .expect("point on its finite spline");
        assert!(dec.compare(inverse.value(), &expected).is_coincident());
        let tangent = inverse.first(&dec)?;
        let parameter = r("1/3");
        let zero = r("0");
        let chart = fit
            .surface
            .jet_box([&parameter, &parameter], [&zero, &zero], 128, &dec)?;
        for (component, expected) in [tangent.x, tangent.y, tangent.z]
            .iter()
            .zip(&chart.first[0])
        {
            let actual = dec.enclose(component, &r("1e-12"))?;
            assert!(actual.lo() <= expected.hi() && actual.hi() >= expected.lo());
        }
        Ok::<_, catallax_numbers::Refusal>(())
    });
    result.unwrap();
    assert!(!work.exhausted);
}

#[test]
fn checked_finite_spline_support_has_exact_incidence_outward_side_and_rigid_transport() {
    use catallax_numbers::{Outcome, Sign, StrictSign};
    use catallax_shapes::{Pose, Rotation, SplineSurface, certify_elliptic_offset_sweep};
    let mut dec = Decider::default();
    let definition = rational_source_definition(&mut dec);
    let (result, work) = with_operation_budget(100_000_000, || {
        let fit = fit_elliptic_offset(
            &definition,
            [r("9/4").into(), r("23/4").into()],
            &r("0.0001"),
            &mut dec,
        )?;
        let proof = certify_elliptic_offset_sweep(
            &fit.surface,
            &definition,
            [&fit.axial[0], &fit.axial[1]],
            &mut dec,
        )?;
        let support = SplineSurface::from_offset(proof);
        let u = Real::from_rational(r("1/3"));
        let v = Real::from_rational(r("1/3"));
        let base = support.lower_edge().point_at(&u, &mut dec)?;
        let point = base.plus(&support.axis().scaled(&(&v * support.height())));
        let inverse = support
            .parameters_at(&point, &mut dec)
            .expect("finite incidence proof")
            .expect("finite represented point");
        assert!(dec.compare(inverse.u(), &u).is_coincident());
        assert!(dec.compare(inverse.v(), &v).is_coincident());
        let rebuilt = inverse.point(&dec)?;
        for (a, b) in [&rebuilt.x, &rebuilt.y, &rebuilt.z]
            .into_iter()
            .zip([&point.x, &point.y, &point.z])
        {
            assert!(dec.compare(a, b).is_coincident());
        }
        let normal = inverse.normal(&dec)?;
        assert!(
            dec.compare(&normal.dot(support.axis()), &Real::zero())
                .is_coincident()
        );
        let centre = support
            .frame()
            .0
            .plus(&support.axis().scaled(&Real::from_i64(4)));
        assert!(matches!(
            dec.sign(&normal.dot(&point.minus(&centre))),
            Outcome::Proven(StrictSign::Positive)
        ));
        assert_eq!(
            support
                .cross_section_side(&centre, &mut dec)
                .expect("interior side proof"),
            Sign::Negative
        );
        let outside = centre.plus(&Vector3::new(0.into(), 3.into(), 0.into()));
        assert_eq!(
            support
                .cross_section_side(&outside, &mut dec)
                .expect("exterior side proof"),
            Sign::Positive
        );
        assert_eq!(
            support
                .cross_section_side(&point, &mut dec)
                .expect("boundary side proof"),
            Sign::Zero
        );
        let near_outside = point.plus(&normal.scaled(&r("0.000001").into()));
        assert_eq!(
            support
                .cross_section_side(&near_outside, &mut dec)
                .expect("near exterior side proof"),
            Sign::Positive
        );
        let beyond_cap = point.plus(
            &support
                .axis()
                .scaled(&(support.height() * &Real::from_i64(2))),
        );
        assert!(
            support
                .parameters_at(&beyond_cap, &mut dec)
                .expect("finite cap exclusion")
                .is_none()
        );
        let pose = Pose::new(
            Rotation::from_quaternion(r("1"), r("1"), r("1"), r("1")).unwrap(),
            Vector3::new(17.into(), (-4).into(), 9.into()),
            &mut dec,
        )?;
        let moved = support.transformed(&pose);
        let moved_point = pose.apply_point(&point);
        let inverse = moved
            .parameters_at(&moved_point, &mut dec)
            .expect("transported incidence proof")
            .expect("rigidly transported point");
        assert!(dec.compare(inverse.u(), &u).is_coincident());
        assert!(dec.compare(inverse.v(), &v).is_coincident());
        let moved_normal = inverse.normal(&dec)?;
        let expected = pose.apply_vector(&normal);
        for (a, b) in [&moved_normal.x, &moved_normal.y, &moved_normal.z]
            .into_iter()
            .zip([&expected.x, &expected.y, &expected.z])
        {
            assert!(dec.compare(a, b).is_coincident());
        }
        assert_eq!(
            moved
                .cross_section_side(&pose.apply_point(&centre), &mut dec)
                .expect("transported interior side proof"),
            Sign::Negative
        );
        Ok::<_, catallax_numbers::Refusal>(())
    });
    result.unwrap();
    assert!(!work.exhausted);
}

#[test]
fn finite_interval_inverse_covers_a_knot_and_refuses_to_choose_the_periodic_seam() {
    use catallax_shapes::{SplineSupport, SplineSurface, certify_elliptic_offset_sweep};
    let mut dec = Decider::default();
    let definition = rational_source_definition(&mut dec);
    let (result, work) = with_operation_budget(100_000_000, || {
        let fit = fit_elliptic_offset(
            &definition,
            [r("9/4").into(), r("23/4").into()],
            &r("0.0001"),
            &mut dec,
        )?;
        let knot = &fit.spans[1];
        let parameter = Rational::from_i64(knot.quarter as i64) + knot.parameter[0].clone();
        let axial = Real::from_rational(r("13/4"));
        let (point, tangent) = definition.chart().point_and_tangent(
            knot.quarter,
            &knot.parameter[0],
            &axial,
            &mut dec,
        )?;
        let proof = certify_elliptic_offset_sweep(
            &fit.surface,
            &definition,
            [&fit.axial[0], &fit.axial[1]],
            &mut dec,
        )?;
        let support = SplineSurface::from_offset(proof);
        let carrier: SplineSupport = support.clone().into();
        let point_box =
            [&point.x, &point.y, &point.z].map(|value| dec.enclose_at(value, 128).unwrap());
        let parameters = carrier
            .parameters_of(&point_box, &dec)?
            .expect("one connected finite knot branch");
        assert!(parameters[0].contains(&parameter));
        assert!(parameters[1].contains(&r("2/7")));
        assert!(&parameters[0].hi() - &parameters[0].lo() < r("1e-12"));
        let evaluated = carrier
            .carrier()
            .evaluate(&parameters[0], &parameters[1], &dec)?;
        let derivatives = carrier
            .carrier()
            .derivatives(&parameters[0], &parameters[1], &dec)?;
        let axial_tangent = support.axis().scaled(support.height());
        for (actual, expected) in evaluated.iter().zip([&point.x, &point.y, &point.z]) {
            let expected = dec.enclose_at(expected, 160)?;
            assert!(actual.lo() <= expected.lo() && actual.hi() >= expected.hi());
        }
        for (actual, expected) in derivatives[0]
            .iter()
            .zip([&tangent.x, &tangent.y, &tangent.z])
            .chain(derivatives[1].iter().zip([
                &axial_tangent.x,
                &axial_tangent.y,
                &axial_tangent.z,
            ]))
        {
            let expected = dec.enclose_at(expected, 160)?;
            assert!(actual.lo() <= expected.lo() && actual.hi() >= expected.hi());
        }
        let (seam, _) = definition
            .chart()
            .point_and_tangent(0, &r("0"), &axial, &mut dec)?;
        let seam_box = [&seam.x, &seam.y, &seam.z].map(|value| dec.enclose_at(value, 128).unwrap());
        assert!(carrier.parameters_of(&seam_box, &dec)?.is_none());
        let beyond = point.plus(&support.axis().scaled(&Real::from_i64(10)));
        assert_eq!(
            carrier.side_at(&beyond, &mut dec).unwrap_err().code,
            "shapes.bspline-domain"
        );
        Ok::<_, catallax_numbers::Refusal>(())
    });
    result.unwrap();
    assert!(!work.exhausted);
}
