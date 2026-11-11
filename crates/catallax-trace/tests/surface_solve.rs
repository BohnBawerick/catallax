#![allow(missing_docs)]
//! Differential gate for SMLib's local and global surface point solvers
//! (`SmSurface::LocalPointSolve_1`, `GlobalPointSolve_1`) against saved
//! outputs of the pinned C++ SMLib. Local cases the port refuses are listed
//! as `refuse` and must refuse.

mod common;

use catallax_trace::extent::Extent2;
use catallax_trace::global_solver::{SolverOperation, SolverRefusal};
use catallax_trace::surface_global::global_point_solve;
use catallax_trace::surface_solve::local_point_solve;
use common::{Input, Output, run_gate, status_code};

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    let surface = input.surface();
    let (low, high) = (input.v2(), input.v2());
    let operation = SolverOperation::from_code(input.i() as usize).expect("SMLib operation code");
    let test = input.v3();
    if op == "global" {
        let tolerance = input.d();
        match global_point_solve(
            &surface,
            Extent2::new(low, high),
            operation,
            test,
            tolerance,
            None,
        ) {
            Ok(solutions) => {
                out.i(common::SUCCESS);
                out.i(solutions.len() as i64);
                for solution in solutions.as_slice() {
                    out.d(solution.start.parameters[0]);
                    out.d(solution.start.parameters[1]);
                    out.d(solution.start.value);
                }
            }
            Err(SolverRefusal::Failed(fault)) => out.i(status_code(fault)),
            Err(SolverRefusal::NotPorted(_)) => out.i(-1),
        }
        return;
    }
    let guess = input.v2();
    let result = local_point_solve(
        &surface,
        Extent2::new(low, high),
        operation,
        test,
        guess,
        None,
    );
    match (op, result) {
        ("solve", Ok(found)) => {
            out.i(common::SUCCESS);
            out.flag(found.found);
            out.v2(found.uv);
            out.d(found.value);
            out.i(found.reason as i64);
            out.i(found.iterations as i64);
        }
        ("solve", Err(SolverRefusal::Failed(fault))) => out.i(status_code(fault)),
        ("refuse", Err(SolverRefusal::NotPorted(_))) => out.i(-1),
        (op, result) => panic!("{op}: unexpected {result:?}"),
    }
}

#[test]
fn surface_point_solver_matches_smlib() {
    let report = run_gate(
        include_str!("fixtures/surface-solve/cases.txt"),
        |_| 16,
        evaluate,
    );
    report.assert_passed(&["tangent", "touching", "singular", "boundary"]);
}
