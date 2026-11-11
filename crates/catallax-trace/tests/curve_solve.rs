#![allow(missing_docs)]
//! Differential gate for SMLib's local curve solvers (`LocalPointSolve`,
//! `LocalCurveIntersect`, `LocalCurveSolve`, `CurveMaxDistanceBetween`,
//! `SimpleCoincidenceChecker`, `IsClosed`, `IsPeriodic`) against saved
//! outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::curve_solve::{
    Coincidence, CurvePair, curve_max_distance_between, is_closed, is_periodic,
    local_curve_intersect, local_curve_solve, local_point_solve, simple_coincidence_check,
};
use catallax_trace::global_solver::{SolverOperation, SolverRefusal};
use common::{Input, Output, run_gate, status_code};

fn refusal<T>(out: &mut Output, result: Result<T, SolverRefusal>) -> Option<T> {
    match result {
        Ok(value) => {
            out.i(common::SUCCESS);
            Some(value)
        }
        Err(SolverRefusal::Failed(fault)) => {
            out.i(status_code(fault));
            None
        }
        Err(SolverRefusal::NotPorted(operation)) => panic!("not ported: {operation:?}"),
    }
}

fn operation(code: i64) -> SolverOperation {
    SolverOperation::from_code(code as usize).expect("SMLib operation code")
}

fn pair(out: &mut Output, pair: &CurvePair) {
    out.d(pair.parameters[0]);
    out.d(pair.parameters[1]);
    out.d(pair.value);
}

fn found_pair(out: &mut Output, found: Option<CurvePair>) {
    out.flag(found.is_some());
    if let Some(found) = found {
        pair(out, &found);
    }
}

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    match op {
        "point" => {
            let curve = input.curve();
            let interval = input.interval();
            let operation = operation(input.i());
            let (test, guess) = (input.v3(), input.d());
            let result = local_point_solve(&curve, interval, operation, test, guess);
            if let Some(found) = refusal(out, result) {
                out.flag(found.is_some());
                if let Some(hit) = found {
                    out.d(hit.parameter);
                    out.d(hit.distance);
                }
            }
        }
        "intersect" => {
            let (a, b) = (input.curve(), input.curve());
            let (ia, ib) = (input.interval(), input.interval());
            let (tolerance, ga, gb) = (input.d(), input.d(), input.d());
            let result = local_curve_intersect(&a, ia, &b, ib, tolerance, [ga, gb]);
            if let Some(found) = refusal(out, result) {
                found_pair(out, found);
            }
        }
        "solve" => {
            let (a, b) = (input.curve(), input.curve());
            let (ia, ib) = (input.interval(), input.interval());
            let operation = operation(input.i());
            let (tolerance, ga, gb) = (input.d(), input.d(), input.d());
            let result = local_curve_solve(&a, ia, &b, ib, operation, tolerance, [ga, gb]);
            if let Some(found) = refusal(out, result) {
                found_pair(out, found);
            }
        }
        "maxdist" => {
            let (a, b) = (input.curve(), input.curve());
            let interval = input.interval();
            let (start, end) = (input.d(), input.d());
            let samples = input.i() as usize;
            let has_needed = input.i() != 0;
            let needed = input.d();
            let result = curve_max_distance_between(
                &a,
                interval,
                &b,
                start,
                end,
                samples,
                has_needed.then_some(needed),
            );
            if let Some(found) = refusal(out, result) {
                out.d(found.distance);
                out.d(found.this_parameter);
                out.d(found.other_parameter);
                out.i(found.incorrect_direction.map_or(2, i64::from));
            }
        }
        "coincide" => {
            let (a, b) = (input.curve(), input.curve());
            let (ia, ib) = (input.interval(), input.interval());
            let tolerance = input.d();
            let boxes = input.i() != 0;
            let box_a = boxes.then(|| a.bounding_boxes(ia).expect("box").1);
            let box_b = boxes.then(|| b.bounding_boxes(ib).expect("box").1);
            let result =
                simple_coincidence_check(&a, ia, &b, ib, tolerance, box_a.as_ref(), box_b.as_ref());
            if let Some(found) = refusal(out, result) {
                let single = |out: &mut Output, p: &CurvePair| {
                    out.i(1);
                    pair(out, p);
                };
                match found {
                    Coincidence::None => out.i(0),
                    Coincidence::Point(p) => {
                        out.i(1);
                        single(out, &p);
                    }
                    Coincidence::Points([p, q]) => {
                        out.i(2);
                        single(out, &p);
                        single(out, &q);
                    }
                    Coincidence::Range([p, q]) => {
                        out.i(1);
                        out.i(2);
                        pair(out, &p);
                        pair(out, &q);
                    }
                }
            }
        }
        "closed" => {
            let curve = input.curve();
            let interval = input.interval();
            let tolerance = input.d();
            out.flag(is_closed(&curve, interval, tolerance));
            out.flag(is_periodic(&curve, interval));
        }
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn curve_solvers_match_smlib() {
    let report = run_gate(
        include_str!("fixtures/curve-solve/cases.txt"),
        |_| 16,
        evaluate,
    );
    report.assert_passed(&["tangent", "touching", "singular"]);
}
