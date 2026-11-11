#![allow(missing_docs)]
//! Differential gate for SMLib's surface-surface intersector
//! (`SmAdvSurfaceIntersector`) stage by stage, against saved outputs of the
//! pinned C++ SMLib read through a probe subclass.

mod common;

use catallax_trace::extent::Extent2;
use catallax_trace::global_solver::SolverRefusal;
use catallax_trace::spline::BSplineSurface;
use catallax_trace::surface_intersector::{BoundaryPlaneTest, Intersector};
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

fn extent2(out: &mut Output, e: Extent2) {
    out.v2(e.min);
    out.v2(e.max);
}

fn run(op: &str, input: &mut Input, out: &mut Output) {
    let a = surface(input);
    let da = domain(input);
    let b = surface(input);
    let db = domain(input);
    let has_tolerance = input.i() != 0;
    let tolerance = input.d();
    let mut intersector = Intersector::new(&a, da, &b, db);
    if has_tolerance {
        intersector.approx_tolerance = tolerance;
    }
    assert!(op == "search" || op == "planes", "unknown operation {op}");
    out.i(common::SUCCESS);
    out.d(intersector.approx_tolerance);
    for i in 0..2 {
        out.flag(intersector.closed_u[i]);
        out.flag(intersector.closed_v[i]);
        extent2(out, intersector.domains[i]);
    }
    let result = intersector.trees().and_then(|trees| {
        Ok((
            trees.root_boxes_meet(intersector.approx_tolerance)?,
            intersector.search(&trees)?,
        ))
    });
    match result {
        Ok((false, _)) => out.i(-1),
        Ok((true, found)) => {
            out.i(common::SUCCESS);
            out.i(i64::from(found.is_some()));
            if let Some(found) = found {
                let p = found.bbox;
                for value in [
                    p.min.x,
                    p.min.y,
                    p.min.z,
                    p.max.x,
                    p.max.y,
                    p.max.z,
                    found.count,
                ] {
                    out.d(value);
                }
                extent2(out, found.domains[0]);
                extent2(out, found.domains[1]);
                if op == "planes" {
                    match intersector.boundary_plane_test(found.bbox) {
                        Ok(BoundaryPlaneTest::Unknown) => {
                            out.i(common::SUCCESS);
                            out.i(0);
                        }
                        Ok(BoundaryPlaneTest::Disjoint) => {
                            out.i(common::SUCCESS);
                            out.i(1);
                        }
                        // SMLib grazes here; the port refuses.
                        Err(SolverRefusal::NotPorted(_)) => {
                            out.i(common::SUCCESS);
                            out.i(-1);
                        }
                        Err(SolverRefusal::Failed(fault)) => out.i(status_code(fault)),
                    }
                }
            }
        }
        Err(SolverRefusal::Failed(fault)) => out.i(status_code(fault)),
        Err(SolverRefusal::NotPorted(what)) => panic!("not ported: {what:?}"),
    }
}

#[test]
fn intersector_search_matches_smlib() {
    let report = run_gate(include_str!("fixtures/intersector/cases.txt"), |_| 16, run);
    report.assert_passed(&[
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}
