#![allow(missing_docs)]
//! Differential gate for SMLib's local surface-surface solvers
//! (`SmSurface::LocalSurfaceIntersect` and `LocalSurfaceSolve`) on B-spline
//! surfaces and planes, against saved outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::extent::Extent2;
use catallax_trace::global_solver::{SolverOperation, SolverRefusal};
use catallax_trace::spline::BSplineSurface;
use catallax_trace::surface_surface::{SurfacePair, local_surface_intersect, local_surface_solve};
use catallax_trace::vector::Vector2;
use common::{Input, Output, run_gate, status_code};

fn domain(input: &mut Input) -> Extent2 {
    let (a, b, c, d) = (input.d(), input.d(), input.d(), input.d());
    Extent2::new(Vector2::new(a, b), Vector2::new(c, d))
}

/// `s` and a surface as `Input::surface` reads it, or `p` and a plane as
/// the plane gate reads it.
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

fn answer(out: &mut Output, result: Result<Option<SurfacePair>, SolverRefusal>) {
    match result {
        Ok(found) => {
            out.i(common::SUCCESS);
            out.flag(found.is_some());
            if let Some(pair) = found {
                out.v2(pair.uv);
                out.v2(pair.other_uv);
                out.d(pair.value);
            }
        }
        Err(SolverRefusal::Failed(fault)) => out.i(status_code(fault)),
        Err(SolverRefusal::NotPorted(what)) => panic!("not ported: {what:?}"),
    }
}

fn run(op: &str, input: &mut Input, out: &mut Output) {
    let a = surface(input);
    let da = domain(input);
    let b = surface(input);
    let db = domain(input);
    let (guess, other_guess) = (input.v2(), input.v2());
    match op {
        "lsi" => {
            let tolerance = input.d();
            let has_normal = input.i() != 0;
            let normal = input.v3();
            let has_previous = input.i() != 0;
            let previous = input.v2();
            let boundary = input.i() != 0;
            answer(
                out,
                local_surface_intersect(
                    &a,
                    da,
                    &b,
                    db,
                    tolerance,
                    guess,
                    other_guess,
                    has_normal.then_some(normal),
                    has_previous.then_some(previous),
                    boundary,
                ),
            );
        }
        "lss" => {
            let operation = SolverOperation::from_code(input.i() as usize).expect("operation");
            answer(
                out,
                local_surface_solve(&a, da, &b, db, operation, guess, other_guess),
            );
        }
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn local_surface_solvers_match_smlib() {
    let report = run_gate(
        include_str!("fixtures/surface-surface/cases.txt"),
        |_| 16,
        run,
    );
    report.assert_passed(&[
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}
