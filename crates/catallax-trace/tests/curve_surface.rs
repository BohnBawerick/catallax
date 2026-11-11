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
