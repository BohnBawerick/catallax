use super::checked::pay_trace_work;
use super::*;
use catallax_numbers::current_operation_work;
use catallax_trace::{surface_curves::IntersectionCurve, surface_intersector::Intersector};

fn sampled_residual(
    fits: &[IntersectionCurve],
    sources: &[BSplineSurface; 2],
    dec: &Decider,
) -> Rational {
    let mut maximum = r(0);
    let exact = |value: f64| Ball::from_rational(&Rational::from_f64_exact(value).unwrap(), BITS);
    for fit in fits {
        let uv = fit
            .uv_curves
            .each_ref()
            .map(|curve| curve.as_ref().unwrap());
        for span in fit
            .curve
            .knots()
            .windows(2)
            .filter(|span| span[0] < span[1])
        {
            for i in 0..=4 {
                let t = span[0] + (span[1] - span[0]) * f64::from(i) / 4.;
                let point = fit.curve.evaluate_point(t).unwrap();
                let point = [exact(point.x), exact(point.y), exact(point.z)];
                for k in 0..2 {
                    let at = uv[k].evaluate_point(t).unwrap();
                    let on = sources[k]
                        .evaluate(&exact(at.x), &exact(at.y), dec)
                        .unwrap();
                    let d: [Ball; 3] = std::array::from_fn(|j| &point[j] - &on[j]);
                    let residual = (d[0].square() + d[1].square() + d[2].square())
                        .sqrt()
                        .unwrap()
                        .hi();
                    maximum = maximum.max(residual);
                }
            }
        }
    }
    maximum
}

#[test]
fn checked_flush_refreshes_coupled_endpoint_derivatives() {
    let ((a, b, trace_work, trace_bill), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let surfaces = [
            Surface::Cone(Cone {
                apex: Point3::new(s(0), s(-11), s(0)),
                axis: Vector3::new(s(0), s(1), s(0)),
                slope: s(1),
                ref_dir: Vector3::new(s(1), s(0), s(0)),
            }),
            Surface::Cylinder(Cylinder {
                origin: Point3::new(s(10), s(12), s(0)),
                axis: Vector3::new(s(0), s(-1), s(0)),
                radius_sq: s(1),
                ref_dir: Vector3::new(s(1), s(0), s(0)),
            }),
        ];
        let windows = numerical_native_windows([&surfaces[0], &surfaces[1]], &mut dec)
            .unwrap()
            .unwrap();
        let sources = windows.sources();
        let a =
            crate::pair::trace_candidates::surface(&sources[0], &Isometry::identity(), &mut dec)
                .unwrap();
        let b =
            crate::pair::trace_candidates::surface(&sources[1], &Isometry::identity(), &mut dec)
                .unwrap();
        let mut intersector = Intersector::new(&a, a.natural_domain(), &b, b.natural_domain());
        intersector.approx_tolerance = 2.5e-8;
        let before = intersector.do_intersection().unwrap();
        let old_residual = sampled_residual(&before, sources, &dec);
        assert!(old_residual > Rational::from_ratio(1, 1_000_000).unwrap());
        intersector.check_output_spans = true;
        let start = current_operation_work().used;
        let (after, trace_work) = intersector.do_intersection_with_work(None, pay_trace_work);
        let trace_bill = current_operation_work().used - start;
        let after = after.unwrap();
        assert!(!trace_work.cancelled);
        assert_eq!(trace_work.accepted, trace_work.attempted);
        assert!(trace_work.search_branches > 0);
        assert!(trace_work.local_iterations > 0);
        assert!(trace_work.surface_evaluations > 0);
        assert!(trace_work.curve_evaluations > 0);
        assert!(trace_work.source_items > 0);
        assert!(trace_work.trace_visits > 0);
        assert_eq!(before.len(), after.len());
        let unmetered = intersector.do_intersection().unwrap();
        for (a, b) in after.iter().zip(&unmetered) {
            assert_eq!(a.curve.knots(), b.curve.knots());
            assert_eq!(a.curve.control_points(), b.curve.control_points());
            for (a, b) in a.uv_curves.iter().zip(&b.uv_curves) {
                let (a, b) = (a.as_ref().unwrap(), b.as_ref().unwrap());
                assert_eq!(a.knots(), b.knots());
                assert_eq!(a.control_points(), b.control_points());
            }
        }
        eprintln!("metered native fit work: {trace_work:?}; request bill={trace_bill}");
        let new_residual = sampled_residual(&after, sources, &dec);
        eprintln!(
            "coupled residual before={} after={}",
            old_residual.to_f64_nearest(),
            new_residual.to_f64_nearest()
        );
        assert!(new_residual <= Rational::from_ratio(1, 10_000_000).unwrap());
        eprintln!(
            "checked native fit sampled residual: {} -> {} mm",
            old_residual.to_f64_nearest(),
            new_residual.to_f64_nearest()
        );
        (a, b, trace_work, trace_bill)
    });
    assert!(!work.exhausted, "{work:?}");
    let mut intersector = Intersector::new(&a, a.natural_domain(), &b, b.natural_domain());
    intersector.approx_tolerance = 2.5e-8;
    intersector.check_output_spans = true;
    for limit in [0, trace_bill / 2, trace_bill - 1, trace_bill] {
        let ((result, report), budget) = with_operation_budget(limit, || {
            intersector.do_intersection_with_work(None, pay_trace_work)
        });
        if limit < trace_bill {
            assert!(result.is_err(), "cancelled trace kept partial curves");
            assert!(report.cancelled);
            assert!(budget.exhausted);
        } else {
            assert!(!result.unwrap().is_empty());
            assert_eq!(report, trace_work);
            assert!(!budget.exhausted);
        }
    }
}
