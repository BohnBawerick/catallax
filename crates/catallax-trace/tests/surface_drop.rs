#![allow(missing_docs)]
//! Differential gate for SMLib's surface point drop (`SmSurface::DropPoint`
//! with `sm_GoodDropResult`) on B-spline surfaces and planes, against saved
//! outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::extent::Extent2;
use catallax_trace::global_solver::{SolverOperation, SolverRefusal};
use catallax_trace::spline::BSplineSurface;
use catallax_trace::surface_drop::drop_point;
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

fn run(op: &str, input: &mut Input, out: &mut Output) {
    assert_eq!(op, "drop");
    let surface = surface(input);
    let bounds = domain(input);
    let point = input.v3();
    let has_guess = input.i() != 0;
    let guess = input.v2();
    let operation = SolverOperation::from_code(input.i() as usize).expect("operation");
    match drop_point(
        &surface,
        bounds,
        point,
        has_guess.then_some(guess),
        operation,
    ) {
        Ok(drop) => {
            out.i(common::SUCCESS);
            out.flag(drop.success);
            if drop.success {
                out.v2(drop.uv);
                out.d(drop.gap);
                out.flag(drop.multivalued);
            }
        }
        Err(SolverRefusal::Failed(fault)) => out.i(status_code(fault)),
        Err(SolverRefusal::NotPorted(what)) => panic!("not ported: {what:?}"),
    }
}

#[test]
fn surface_point_drop_matches_smlib() {
    let report = run_gate(include_str!("fixtures/surface-drop/cases.txt"), |_| 16, run);
    report.assert_passed(&[
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}
