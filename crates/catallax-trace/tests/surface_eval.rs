#![allow(missing_docs)]
//! Differential gate for B-spline surface evaluation
//! (`SmBSplineSurface::Evaluate`, `EvaluatePoint`, `EvaluateSimple`) against
//! saved outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::surface_eval::{evaluate, evaluate_normal, evaluate_point, evaluate_simple};
use catallax_trace::vector::Vector2;
use common::{Input, Output, run_gate};

fn run(op: &str, input: &mut Input, out: &mut Output) {
    let surface = input.surface();
    let uv = Vector2::new(input.d(), input.d());
    if op == "point" {
        if let Some(point) = out.status(evaluate_point(&surface, uv)) {
            out.v3(point);
        }
        return;
    }
    if op == "normal" {
        let sides = [input.i() != 0, input.i() != 0];
        out.status(Ok(()));
        out.v3(evaluate_normal(&surface, uv, sides));
        return;
    }
    let (hu, hv) = (input.i() as usize, input.i() as usize);
    let sides = [input.i() != 0, input.i() != 0];
    let (upper, nonzero) = (input.i() != 0, input.i() != 0);
    let result = match op {
        "eval" => evaluate(&surface, uv, hu, hv, sides, upper, nonzero, true),
        "simple" => evaluate_simple(&surface, uv, hu, hv, sides, upper),
        other => panic!("unknown operation {other}"),
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

#[test]
fn surface_evaluation_matches_smlib() {
    let report = run_gate(include_str!("fixtures/surface-eval/cases.txt"), |_| 4, run);
    report.assert_passed(&["touching", "singular", "boundary"]);
}
