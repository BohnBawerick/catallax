#![allow(missing_docs)]
//! Differential gate for SMLib's local curve-surface solvers
//! (`SmSurface::LocalCurveIntersect`, `LocalCurveSolve`) and its singularity
//! tests (`IsSingularity`, `GetSingularities`) against saved outputs of the
//! pinned C++ SMLib.

mod common;

use catallax_trace::curve_surface::{
    CurveSurfacePoint, Singular, get_singularities, is_singularity, local_curve_intersect,
    local_curve_solve, simple_coincidence_check,
};
use catallax_trace::curve_surface_global::global_curve_intersect;
use catallax_trace::extent::{Extent1, Extent2};
use catallax_trace::global_solver::{SolverOperation, SolverRefusal};
use catallax_trace::scalar::ZONE_TOL_3D;
use catallax_trace::solution::SolutionType;
use common::{Input, Output, run_gate, status_code};

fn found(out: &mut Output, result: Result<Option<CurveSurfacePoint>, SolverRefusal>) {
    match result {
        Ok(point) => {
            out.i(common::SUCCESS);
            out.flag(point.is_some());
            if let Some(point) = point {
                out.d(point.t);
                out.v2(point.uv);
                out.d(point.distance);
            }
        }
        Err(SolverRefusal::Failed(fault)) => out.i(status_code(fault)),
        Err(SolverRefusal::NotPorted(operation)) => panic!("not ported: {operation:?}"),
    }
}

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    let surface = input.surface();
    if op == "sing" {
        let uv = input.v2();
        let tolerance = input.d();
        let point_only = input.i() != 0;
        let direction = is_singularity(&surface, uv, tolerance, ZONE_TOL_3D, point_only)
            .expect("singularity test");
        out.flag(direction.is_some());
        out.i(match direction {
            Some(Singular::U) => 0,
            Some(Singular::V) => 1,
            Some(Singular::Both) => 3,
            None => 4,
        });
        let (flags, poles) = get_singularities(&surface, ZONE_TOL_3D);
        out.i(i64::from(flags));
        for pole in poles {
            match pole {
                Some(point) => out.v3(point),
                None => out.i(0),
            }
        }
        return;
    }
    let curve = input.curve();
    let (low, high) = (input.v2(), input.v2());
    let domain = Extent2::new(low, high);
    let (t0, t1) = (input.d(), input.d());
    let interval = Extent1::new(t0, t1);
    match op {
        "intersect" => {
            let (tolerance, t, uv) = (input.d(), input.d(), input.v2());
            found(
                out,
                local_curve_intersect(
                    &surface,
                    domain,
                    &curve,
                    interval,
                    tolerance,
                    (t, uv),
                    ZONE_TOL_3D,
                ),
            );
        }
        "global" | "grefuse" => {
            let tolerance = input.d();
            match global_curve_intersect(&surface, domain, &curve, interval, tolerance, None) {
                Ok(solutions) if op == "global" => {
                    out.i(common::SUCCESS);
                    out.i(solutions.len() as i64);
                    for solution in solutions.as_slice() {
                        out.i(solution.kind as i64);
                        for k in 0..3 {
                            out.d(solution.start.parameters[k]);
                        }
                        out.d(solution.start.value);
                        if solution.kind == SolutionType::RangeOfValues {
                            for k in 0..3 {
                                out.d(solution.end.parameters[k]);
                            }
                            out.d(solution.end.value);
                        }
                    }
                }
                Err(SolverRefusal::NotPorted(_)) if op == "grefuse" => out.i(-1),
                Err(SolverRefusal::Failed(fault)) if op == "global" => out.i(status_code(fault)),
                // Mismatched expectations show as a failing case.
                Ok(_) => out.i(-2),
                Err(_) => out.i(-3),
            }
        }
        "coincide" => {
            let tolerance = input.d();
            let partial = input.i() != 0;
            match simple_coincidence_check(
                &surface, domain, &curve, interval, tolerance, partial, None,
            ) {
                Ok(range) => {
                    out.i(common::SUCCESS);
                    out.flag(range.is_some());
                    for end in range.iter().flatten() {
                        out.d(end.t);
                        out.v2(end.uv);
                        out.d(end.distance);
                    }
                }
                Err(SolverRefusal::Failed(fault)) => out.i(status_code(fault)),
                Err(SolverRefusal::NotPorted(operation)) => panic!("not ported: {operation:?}"),
            }
        }
        "solve" => {
            let operation =
                SolverOperation::from_code(input.i() as usize).expect("SMLib operation code");
            let (_tolerance, t, uv) = (input.d(), input.d(), input.v2());
            found(
                out,
                local_curve_solve(&surface, domain, &curve, interval, operation, (t, uv)),
            );
        }
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn curve_surface_solvers_match_smlib() {
    let report = run_gate(
        include_str!("fixtures/curve-surface/cases.txt"),
        |_| 16,
        evaluate,
    );
    report.assert_passed(&["tangent", "touching", "singular", "near-singular"]);
}

#[test]
fn curve_surface_candidates_are_bit_identical_with_one_two_and_four_workers() {
    common::assert_same_with_workers(
        include_str!("fixtures/curve-surface/cases.txt"),
        |_| 16,
        evaluate,
    );
}

/// Regressions from the step 6 review: typed failures in the point-curve
/// fallback, the solution layout it returns, and answers outside the
/// requested surface domain.
mod review {
    use catallax_trace::curve_surface_global::{
        CURVE_OBJECT, SURFACE_OBJECT, global_curve_intersect,
    };
    use catallax_trace::extent::{Extent1, Extent2};
    use catallax_trace::global_solver::SolverRefusal;
    use catallax_trace::spline::{BSplineCurve, BSplineSurface, ControlPoint, NO_WEIGHT};
    use catallax_trace::vector::Vector2;

    fn point(x: f64, y: f64, z: f64, w: Option<f64>) -> ControlPoint {
        ControlPoint {
            coordinates: match w {
                Some(w) => [x * w, y * w, z * w, w],
                None => [x, y, z, NO_WEIGHT],
            },
        }
    }

    fn bilinear(points: [ControlPoint; 4]) -> BSplineSurface {
        let knots = vec![0.0, 0.0, 1.0, 1.0];
        BSplineSurface::new([1, 1], [knots.clone(), knots], [2, 2], points.to_vec())
            .expect("bilinear")
    }

    fn unit_plane() -> BSplineSurface {
        bilinear([
            point(0.0, 0.0, 0.0, None),
            point(0.0, 1.0, 0.0, None),
            point(1.0, 0.0, 0.0, None),
            point(1.0, 1.0, 0.0, None),
        ])
    }

    fn constant_cubic(x: f64, y: f64, z: f64) -> BSplineCurve {
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        BSplineCurve::new(3, knots, vec![point(x, y, z, None); 4]).expect("cubic")
    }

    fn unit() -> Extent2 {
        Extent2::new(Vector2::new(0.0, 0.0), Vector2::new(1.0, 1.0))
    }

    /// A degenerate curve on a surface whose rational net collapses to the
    /// origin: the surface point search fails, and the failure must reach
    /// the caller instead of an empty answer. (SMLib ignores the search
    /// status and answers `(0, 0, 0)`.)
    #[test]
    fn degenerate_curve_fallback_keeps_typed_failures() {
        let zero = bilinear([point(0.0, 0.0, 0.0, Some(1.0)); 4]);
        let curve = constant_cubic(0.0, 0.0, 0.0);
        let result =
            global_curve_intersect(&zero, unit(), &curve, Extent1::new(0.0, 1.0), 1e-6, None);
        assert!(
            matches!(result, Err(SolverRefusal::Failed(_))),
            "{result:?}"
        );
    }

    /// The point-curve fallback answers in the curve-surface layout: three
    /// variables `(t, u, v)` and both objects.
    #[test]
    fn degenerate_curve_answer_has_three_variables() {
        let curve = constant_cubic(0.25, 0.75, 0.0);
        let found = global_curve_intersect(
            &unit_plane(),
            unit(),
            &curve,
            Extent1::new(0.0, 1.0),
            1e-5,
            None,
        )
        .expect("answer");
        assert_eq!(found.len(), 1);
        let solution = found.as_slice()[0];
        assert_eq!(solution.variables, 3);
        assert_eq!(solution.object_count, 2);
        assert_eq!(solution.objects[0], Some(CURVE_OBJECT));
        assert_eq!(solution.objects[1], Some(SURFACE_OBJECT));
        assert_eq!(solution.start.parameters[..3], [0.0, 0.25, 0.75]);
    }

    /// A curve that stays outside the requested half of the unit plane but
    /// starts on the other half: the end-point fallback must not answer
    /// outside the requested domain. (SMLib answers `(0, 0.75, 0.5)`.)
    #[test]
    fn end_point_fallback_stays_in_the_requested_domain() {
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let curve = BSplineCurve::new(
            2,
            knots,
            vec![
                point(0.75, 0.5, 0.0, None),
                point(0.75, 0.5, 0.5, None),
                point(0.85, 0.5, 1.0, None),
            ],
        )
        .expect("quadratic");
        let half = Extent2::new(Vector2::new(0.0, 0.0), Vector2::new(0.5, 1.0));
        let found = global_curve_intersect(
            &unit_plane(),
            half,
            &curve,
            Extent1::new(0.0, 1.0),
            1e-5,
            None,
        )
        .expect("search");
        for solution in found.as_slice() {
            let uv = Vector2::new(solution.start.parameters[1], solution.start.parameters[2]);
            assert!(half.contains_point(uv, 0.0), "{solution:?}");
        }
    }
}

mod unlimited_line {
    use catallax_trace::extent::Extent2;
    use catallax_trace::global_solver::SolverRefusal;
    use catallax_trace::line_surface::global_line_intersect;
    use catallax_trace::solution::SolutionType;
    use catallax_trace::spline::{BSplineSurface, ControlPoint};
    use catallax_trace::vector::{Vector2, Vector3};

    /// A line without an interval lying in the surface over `t` in
    /// `[-2, 2]`: SMLib probes its coincidence on the unit segment only, so
    /// ranges built there would omit the overlap outside `[0, 1]`. The port
    /// must not answer with ranges confined to that segment.
    #[test]
    fn an_unlimited_coincident_line_is_not_cut_to_the_unit_segment() {
        // S(u, v) = (u, v, u^2) over [-1, 1] x [-2, 2].
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
        let domain = Extent2::new(Vector2::new(-1.0, -2.0), Vector2::new(1.0, 2.0));
        let found = global_line_intersect(
            &surface,
            domain,
            Vector3::new(0.5, 0.0, 0.25),
            Vector3::new(0.0, 1.0, 0.0),
            None,
            1.0e-5,
        );
        match found {
            Err(SolverRefusal::NotPorted(_)) => {}
            Ok(solutions) => {
                for solution in solutions.as_slice() {
                    if solution.kind == SolutionType::RangeOfValues {
                        let (a, b) = (solution.start.parameters[0], solution.end.parameters[0]);
                        assert!(
                            a.min(b) < 0.0 && a.max(b) > 1.0,
                            "range [{a}, {b}] is cut to the unit segment"
                        );
                    }
                }
            }
            Err(other) => panic!("unexpected failure {other:?}"),
        }
    }
}
