#![allow(missing_docs)]
//! Differential gate for SMLib's global point-curve solver
//! (`SmCurve::GlobalPointSolve` with `SmCPGlobalSolver`) and the curve point
//! drop (`SmCurve::DropPoint`) against saved outputs of the pinned C++
//! SMLib.

mod common;

use catallax_trace::curve_global::{drop_point, global_point_solve};
use catallax_trace::global_solver::{SolutionRequested, SolverOperation, SolverRefusal};
use catallax_trace::solution::SolutionType;
use common::{Input, Output, run_gate, status_code};

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    let curve = input.curve();
    let interval = input.interval();
    let operation = SolverOperation::from_code(input.i() as usize).expect("operation");
    let test = input.v3();
    let tolerance = input.d();
    match op {
        "global" | "grefuse" => {
            let requested = SolutionRequested::from_code(input.i() as usize).expect("request");
            match global_point_solve(&curve, interval, operation, test, tolerance, requested) {
                Ok(solutions) if op == "global" => {
                    out.i(common::SUCCESS);
                    out.i(solutions.len() as i64);
                    for solution in solutions.as_slice() {
                        out.i(solution.kind as i64);
                        out.d(solution.start.parameters[0]);
                        out.d(solution.start.value);
                        if solution.kind == SolutionType::RangeOfValues {
                            out.d(solution.end.parameters[0]);
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
        "drop" | "drefuse" => {
            let has_inward = input.i() != 0;
            let inward = input.v3();
            let has_guess = input.i() != 0;
            let guess = input.d();
            match drop_point(
                &curve,
                interval,
                test,
                has_inward.then_some(inward),
                tolerance,
                has_guess.then_some(guess),
                operation,
            ) {
                Ok(dropped) if op == "drop" => {
                    out.i(common::SUCCESS);
                    out.flag(dropped.success);
                    out.d(dropped.parameter);
                    out.d(dropped.distance);
                }
                Err(SolverRefusal::NotPorted(_)) if op == "drefuse" => out.i(-1),
                Err(SolverRefusal::Failed(fault)) if op == "drop" => out.i(status_code(fault)),
                Ok(_) => out.i(-2),
                Err(_) => out.i(-3),
            }
        }
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn global_point_curve_solver_matches_smlib() {
    let report = run_gate(
        include_str!("fixtures/curve-global/cases.txt"),
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
