#![allow(missing_docs)]
//! Differential gate for SMLib's `SmPlane` (its bilinear net, analytic
//! evaluation, normal, flags and `DropPointFast`) and for the global point
//! and curve solvers on a plane, against saved outputs of the pinned C++
//! SMLib.

mod common;

use catallax_trace::curve_surface::is_singularity;
use catallax_trace::curve_surface_global::global_curve_intersect;
use catallax_trace::extent::{Extent1, Extent2};
use catallax_trace::global_solver::{SolverOperation, SolverRefusal};
use catallax_trace::plane::drop_point_fast;
use catallax_trace::scalar::ZONE_TOL_3D;
use catallax_trace::solution::SolutionArray;
use catallax_trace::spline::BSplineSurface;
use catallax_trace::surface_eval::{evaluate, evaluate_normal, evaluate_point, evaluate_simple};
use catallax_trace::surface_global::global_point_solve;
use catallax_trace::surface_solve::{is_closed, is_on_seam};
use catallax_trace::vector::Vector2;
use common::{Input, Output, run_gate, status_code};

fn domain(input: &mut Input) -> Extent2 {
    let (a, b, c, d) = (input.d(), input.d(), input.d(), input.d());
    Extent2::new(Vector2::new(a, b), Vector2::new(c, d))
}

fn plane(input: &mut Input) -> BSplineSurface {
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

fn solutions(out: &mut Output, found: &SolutionArray, variables: usize) {
    out.i(found.len() as i64);
    for solution in found.as_slice() {
        out.i(solution.kind as i64);
        for k in 0..variables {
            out.d(solution.start.parameters[k]);
        }
        out.d(solution.start.value);
    }
}

fn refused(out: &mut Output, refusal: SolverRefusal) {
    match refusal {
        SolverRefusal::Failed(fault) => out.i(status_code(fault)),
        SolverRefusal::NotPorted(what) => panic!("not ported: {what:?}"),
    }
}

fn run(op: &str, input: &mut Input, out: &mut Output) {
    let surface = plane(input);
    match op {
        "net" => {
            out.i(common::SUCCESS);
            let d = surface.natural_domain();
            out.d(d.min.x);
            out.d(d.min.y);
            out.d(d.max.x);
            out.d(d.max.y);
            for point in surface.control_points() {
                out.i(common::SUCCESS);
                out.v3(point.euclidean());
            }
        }
        "point" => {
            let uv = input.v2();
            if let Some(point) = out.status(evaluate_point(&surface, uv)) {
                out.v3(point);
            }
        }
        "eval" => {
            let uv = input.v2();
            let (hu, hv) = (input.i() as usize, input.i() as usize);
            let (upper, simple) = (input.i() != 0, input.i() != 0);
            let result = if simple {
                evaluate_simple(&surface, uv, hu, hv, [true, true], upper)
            } else {
                evaluate(&surface, uv, hu, hv, [true, true], upper, true, true)
            };
            if let Some(cells) = out.status(result) {
                for k in 0..=hu {
                    for l in 0..=hv {
                        if !upper || k + l <= hu {
                            out.v3(cells[k * (hv + 1) + l]);
                        }
                    }
                }
            }
        }
        "normal" => {
            let uv = input.v2();
            out.i(common::SUCCESS);
            out.v3(evaluate_normal(&surface, uv, [true, true]));
        }
        "flags" => {
            let uv = input.v2();
            out.i(common::SUCCESS);
            let singular = is_singularity(&surface, uv, 1e-5, ZONE_TOL_3D, false).expect("plane");
            out.flag(singular.is_some());
            // SMLib's SM_SP_NEITHER.
            out.i(4);
            out.flag(is_closed(&surface, surface.natural_domain(), 0, None));
            let (_, on_seam, _) = is_on_seam(&surface, uv, ZONE_TOL_3D);
            out.flag(on_seam);
        }
        "drop" => {
            let operation = SolverOperation::from_code(input.i() as usize).expect("operation");
            let bounds = domain(input);
            let test = input.v3();
            let tolerance = input.d();
            let has_target = input.i() != 0;
            let target = input.d();
            match drop_point_fast(
                &surface,
                bounds,
                operation,
                test,
                tolerance,
                has_target.then_some(target),
            ) {
                Ok(found) => {
                    out.i(common::SUCCESS);
                    out.i(i64::from(found.is_some()));
                    if let Some((uv, distance)) = found {
                        out.i(1);
                        out.v2(uv);
                        out.d(distance);
                    }
                }
                Err(fault) => out.i(status_code(fault)),
            }
        }
        "global" => {
            let operation = SolverOperation::from_code(input.i() as usize).expect("operation");
            let bounds = domain(input);
            let test = input.v3();
            let tolerance = input.d();
            match global_point_solve(&surface, bounds, operation, test, tolerance, None) {
                Ok(found) => {
                    out.i(common::SUCCESS);
                    solutions(out, &found, 2);
                }
                Err(refusal) => refused(out, refusal),
            }
        }
        "cs" => {
            let curve = input.curve();
            let bounds = domain(input);
            let interval = Extent1::new(input.d(), input.d());
            let tolerance = input.d();
            match global_curve_intersect(&surface, bounds, &curve, interval, tolerance, None) {
                Ok(found) => {
                    out.i(common::SUCCESS);
                    solutions(out, &found, 3);
                }
                Err(refusal) => refused(out, refusal),
            }
        }
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn plane_matches_smlib() {
    let report = run_gate(include_str!("fixtures/plane/cases.txt"), |_| 16, run);
    report.assert_passed(&[
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}
