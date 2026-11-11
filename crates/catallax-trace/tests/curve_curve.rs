#![allow(missing_docs)]
//! Differential gate for SMLib's global curve-curve intersector
//! (`SmCurve::GlobalCurveIntersect` with `SmCCIGlobalSolver`) and its
//! whole-curve coincidence check (`SmCurve::GlobalCoincidenceChecker`)
//! against saved outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::curve_curve::{global_coincidence_check, global_curve_intersect};
use catallax_trace::global_solver::SolverRefusal;
use catallax_trace::solution::{SolutionArray, SolutionType};
use common::{Input, Output, run_gate, status_code};

fn solutions(out: &mut Output, solutions: &SolutionArray) {
    out.i(solutions.len() as i64);
    for solution in solutions.as_slice() {
        out.i(solution.kind as i64);
        out.d(solution.start.parameters[0]);
        out.d(solution.start.parameters[1]);
        out.d(solution.start.value);
        if solution.kind == SolutionType::RangeOfValues {
            out.d(solution.end.parameters[0]);
            out.d(solution.end.parameters[1]);
            out.d(solution.end.value);
        }
    }
}

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    let a = input.curve();
    let ia = input.interval();
    let b = input.curve();
    let ib = input.interval();
    let tolerance = input.d();
    let flag = input.i() != 0;
    match op {
        "global" | "grefuse" => match global_curve_intersect(&a, ia, &b, ib, tolerance, flag) {
            Ok(found) if op == "global" => {
                out.i(common::SUCCESS);
                solutions(out, &found);
            }
            Err(SolverRefusal::NotPorted(_)) if op == "grefuse" => out.i(-1),
            Err(SolverRefusal::Failed(fault)) if op == "global" => out.i(status_code(fault)),
            // Mismatched expectations show as a failing case.
            Ok(_) => out.i(-2),
            Err(_) => out.i(-3),
        },
        "coincide" => match global_coincidence_check(&a, ia, &b, ib, tolerance, flag) {
            Ok(check) => {
                out.i(common::SUCCESS);
                out.flag(check.more_work);
                solutions(out, &check.solutions);
            }
            Err(SolverRefusal::Failed(fault)) => out.i(status_code(fault)),
            Err(SolverRefusal::NotPorted(_)) => out.i(-3),
        },
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn global_curve_curve_intersector_matches_smlib() {
    let report = run_gate(
        include_str!("fixtures/curve-curve/cases.txt"),
        |_| 16,
        evaluate,
    );
    report.assert_passed(&[
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}
