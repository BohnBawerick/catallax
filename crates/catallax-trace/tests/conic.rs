#![allow(missing_docs)]
//! Differential gate for SMLib's cubic and quartic solvers
//! (`SmPolynomial`) and the planar conics of the Dupin indicatrix
//! (`SmConic`) against saved outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::conic::{Conic, ConicKind};
use catallax_trace::polynomial::{solve_cubic, solve_quartic};
use catallax_trace::vector::Vector2;
use common::{Input, Output, run_gate};

fn conic(input: &mut Input) -> Conic {
    let kind = match input.i() {
        0 => ConicKind::Ellipse,
        1 => ConicKind::Hyperbola,
        2 => ConicKind::Parabola,
        _ => ConicKind::Line,
    };
    let (a, b) = (input.d(), input.d());
    let mut conic = Conic::new(kind, a, b);
    for _ in 0..input.i() {
        let (angle, x, y) = (input.d(), input.d(), input.d());
        conic.transform_2d(angle, Vector2::new(x, y));
    }
    conic
}

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    match op {
        "cubic" => {
            let c = [input.d(), input.d(), input.d(), input.d()];
            let (count, roots) = solve_cubic(c, input.d());
            out.i(common::SUCCESS);
            out.i(count as i64);
            for &root in &roots[..count] {
                out.d(root);
            }
        }
        "quartic" => {
            let c = [input.d(), input.d(), input.d(), input.d(), input.d()];
            let (count, roots) = solve_quartic(c, input.d());
            out.i(common::SUCCESS);
            out.i(count as i64);
            for &root in &roots[..count] {
                out.d(root);
            }
        }
        "conic" => {
            let one = conic(input);
            let two = conic(input);
            let tolerance = input.d();
            for t in [-0.75, -0.25, 0.3, 0.9, 1.5, 2.6] {
                let p = one.evaluate_point(t);
                out.v2(p);
                let back = one.invert_point(p);
                out.flag(back.is_some());
                if let Some(back) = back {
                    out.d(back);
                }
            }
            for c in two.implicit_equation() {
                out.d(c);
            }
            let found = one.intersect_conic(&two, tolerance);
            out.i(common::SUCCESS);
            out.i(found.len() as i64);
            for (t, t2, distance) in found {
                out.d(t);
                out.d(t2);
                out.flag(distance <= tolerance);
            }
        }
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn conic_matches_smlib() {
    let report = run_gate(include_str!("fixtures/conic/cases.txt"), |_| 16, evaluate);
    report.assert_passed(&[
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}
