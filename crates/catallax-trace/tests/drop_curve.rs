#![allow(missing_docs)]
//! Differential gate for SMLib's curve drop onto a surface
//! (`SmSurface::DropAndTrimCurve` with `SmSurfaceDropCurve`) against saved
//! outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::extent::Extent2;
use catallax_trace::global_solver::SolverRefusal;
use catallax_trace::spline::{BSplineCurve, BSplineSurface};
use catallax_trace::surface_drop_curve::{
    DroppedCurves, UvCurve, drop_and_trim_curve, drop_curve_onto,
};
use catallax_trace::vector::Vector2;
use common::{Input, Output, run_gate, status_code};

fn domain(input: &mut Input) -> Extent2 {
    let (a, b, c, d) = (input.d(), input.d(), input.d(), input.d());
    Extent2::new(Vector2::new(a, b), Vector2::new(c, d))
}

fn surface(input: &mut Input) -> BSplineSurface {
    if input.word() == "s" {
        return input.surface();
    }
    let kind = input.i();
    let (origin, x, y) = (input.v3(), input.v3(), input.v3());
    let scale = input.v2();
    let bounds = domain(input);
    if kind == 1 {
        BSplineSurface::canonical_plane(origin, x, y).expect("plane")
    } else {
        BSplineSurface::plane(origin, x, y, scale, bounds).expect("plane")
    }
}

/// A curve as SMLib stores it: degree, counts, knots, homogeneous points.
fn raw(out: &mut Output, curve: &BSplineCurve) {
    out.i(curve.degree() as i64);
    out.i(curve.control_points().len() as i64);
    out.i(curve.knots().len() as i64);
    for &knot in curve.knots() {
        out.d(knot);
    }
    for point in curve.control_points() {
        for value in point.coordinates {
            out.d(value);
        }
    }
}

fn dropped_curves(out: &mut Output, dropped: &DroppedCurves) {
    out.i(common::SUCCESS);
    out.i(dropped.curves.len() as i64);
    for curve in &dropped.curves {
        match curve {
            UvCurve::Spline(spline) => raw(out, spline),
            UvCurve::Line(line) => raw(out, &line.make_nurb().expect("line nurb")),
        }
    }
    out.d(dropped.max_drop);
    out.d(dropped.max_deviation);
}

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    let s = surface(input);
    let d = domain(input);
    let curve = input.curve();
    let interval = input.interval();
    let tolerance = input.d();
    let keep_all = input.i() != 0;
    let on_surface = input.i() != 0;
    let (answer, refuse) = if op == "general" || op == "grefuse" {
        (
            drop_curve_onto(&s, d, &curve, interval, tolerance, keep_all),
            "grefuse",
        )
    } else {
        (
            drop_and_trim_curve(&s, d, &curve, interval, tolerance, keep_all, on_surface),
            "drefuse",
        )
    };
    match answer {
        Ok(dropped) if op != refuse => dropped_curves(out, &dropped),
        Err(SolverRefusal::NotPorted(_)) if op == refuse => out.i(-1),
        Err(SolverRefusal::Failed(fault)) if op != refuse => out.i(status_code(fault)),
        // Mismatched expectations show as a failing case.
        Ok(_) => out.i(-2),
        Err(_) => out.i(-3),
    }
}

#[test]
fn curve_drop_matches_smlib() {
    let report = run_gate(
        include_str!("fixtures/drop-curve/cases.txt"),
        |_| 16,
        evaluate,
    );
    report.assert_passed(&["touching", "singular", "near-singular", "boundary"]);
    common::assert_same_with_workers(
        include_str!("fixtures/drop-curve/cases.txt"),
        |_| 16,
        evaluate,
    );
}

mod review {
    use catallax_trace::extent::{Extent1, Extent2};
    use catallax_trace::spline::{BSplineCurve, BSplineSurface, ControlPoint};
    use catallax_trace::surface_drop_curve::drop_and_trim_curve;
    use catallax_trace::surface_eval::evaluate_point;
    use catallax_trace::vector::{Vector2, Vector3};

    /// A curved planar curve whose start tangent and knot samples lie on a
    /// parameter line must not be dropped as that straight line: every
    /// parameter curve returned keeps the curve within the tolerance.
    /// A parameter line of a surface stored with unit weights (a rational
    /// flag, polynomial geometry) is still proven and dropped as that line.
    #[test]
    fn a_unit_weight_parameter_line_is_proven() {
        let mut points = Vec::new();
        for (x, z) in [(-1.0, 1.0), (0.0, -1.0), (1.0, 1.0)] {
            for y in [-2.0, 2.0] {
                points.push(ControlPoint {
                    coordinates: [x, y, z, 1.0],
                });
            }
        }
        let surface = BSplineSurface::new(
            [2, 1],
            [
                vec![-1.0, -1.0, -1.0, 1.0, 1.0, 1.0],
                vec![-2.0, -2.0, 2.0, 2.0],
            ],
            [3, 2],
            points,
        )
        .expect("surface");
        let line = BSplineCurve::new(
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                ControlPoint {
                    coordinates: [0.5, 0.0, 0.25, 1.0],
                },
                ControlPoint {
                    coordinates: [0.5, 1.0, 0.25, 1.0],
                },
            ],
        )
        .expect("line");
        let domain = Extent2::new(Vector2::new(-1.0, -2.0), Vector2::new(1.0, 2.0));
        let dropped = drop_and_trim_curve(
            &surface,
            domain,
            &line,
            Extent1::new(0.0, 1.0),
            1.0e-5,
            true,
            true,
        )
        .expect("the line is a parameter line of the surface");
        assert_eq!(dropped.curves.len(), 1);
    }

    #[test]
    fn a_curved_curve_is_not_dropped_as_a_parameter_line() {
        let plane = BSplineSurface::plane(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector2::new(1.0, 1.0),
            Extent2::new(Vector2::new(-2.0, -2.0), Vector2::new(2.0, 2.0)),
        )
        .expect("plane");
        let points = [
            (0.0, 0.0),
            (0.25, 0.0),
            (0.5, 1.0),
            (0.75, -1.5),
            (1.0, 0.0),
        ]
        .map(|(x, y)| ControlPoint {
            coordinates: [x, y, 0.0, 1.0],
        })
        .to_vec();
        let curve = BSplineCurve::new(
            4,
            vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            points,
        )
        .expect("curve");
        let domain = Extent2::new(Vector2::new(-2.0, -2.0), Vector2::new(2.0, 2.0));
        let tolerance = 1.0e-6;
        for keep_all in [false, true] {
            let dropped = drop_and_trim_curve(
                &plane,
                domain,
                &curve,
                Extent1::new(0.0, 1.0),
                tolerance,
                keep_all,
                true,
            )
            .expect("drop");
            assert!(!dropped.curves.is_empty(), "the curve lies on the plane");
            for uv_curve in &dropped.curves {
                let span = uv_curve.natural_interval();
                for i in 0..=16 {
                    let t = span.evaluate(f64::from(i) / 16.0);
                    let uv = uv_curve.evaluate_point(t).expect("uv");
                    let on_plane =
                        evaluate_point(&plane, Vector2::new(uv.x, uv.y)).expect("plane point");
                    let target = curve.evaluate_point(t).expect("curve point");
                    assert!(
                        on_plane.distance(target) <= tolerance,
                        "t = {t}: {on_plane:?} is {} from {target:?}",
                        on_plane.distance(target)
                    );
                }
            }
        }
    }
}
