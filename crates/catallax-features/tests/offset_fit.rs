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
    let original = EllipticCylinder::new(
        Point3::origin(),
        Vector3::new(Real::from_i64(3), Real::zero(), Real::from_i64(4)),
        Vector3::new(Real::zero(), Real::one(), Real::zero()),
        Real::from_i64(4),
        r("64/25").into(),
        &mut dec,
    )
    .unwrap();
    let definition = EllipticOffsetDefinition::new(original, r("-1/4").into(), &mut dec).unwrap();
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
