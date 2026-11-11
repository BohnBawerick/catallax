#![allow(missing_docs)]
//! Differential gate for SMLib's local solvers (`SmLocalSolve1d`,
//! `SmLocalSolveNd`) against saved outputs of the pinned C++ SMLib. The test
//! functions repeat `fixtures/local-solve/oracle.cpp` operation for operation.

mod common;

use catallax_trace::Fault;
use catallax_trace::extent::Extent1;
use catallax_trace::local_solve::{
    BoundaryHandler, EvalFunction, EvalNFunctions, ExtentN, LocalSolve1d, LocalSolveNd,
    SolverAlgorithm, compute_jacobian,
};
use catallax_trace::matrix::Matrix;
use common::{Input, Output, run_gate};

struct Family1 {
    kind: i64,
    p: [f64; 4],
}

impl EvalFunction for Family1 {
    fn evaluate(
        &mut self,
        t: f64,
        f: &mut f64,
        fp: &mut f64,
        found: &mut bool,
        _signal_errors: bool,
    ) -> Result<(), Fault> {
        let p = self.p;
        *found = false;
        match self.kind {
            1 => {
                *f = (t - p[0]) * (t - p[1]) * (t - p[2]);
                *fp = (t - p[1]) * (t - p[2]) + (t - p[0]) * (t - p[2]) + (t - p[0]) * (t - p[1]);
            }
            2 => {
                *f = (t - p[0]) * (t - p[0]) + p[1];
                *fp = 2.0 * (t - p[0]);
            }
            3 => {
                let d = t - p[0];
                *f = if d >= 0.0 { d.sqrt() } else { -(-d).sqrt() };
                *fp = 0.5 / (d.abs() + p[1]).sqrt();
            }
            4 => {
                if t > p[1] {
                    return Err(Fault::Degenerate);
                }
                *f = t - p[0];
                *fp = 1.0;
            }
            5 => {
                *f = (t - p[0]) * (t - p[0]) * (t - p[1]) * (t - p[1]) - p[2];
                *fp = 2.0 * (t - p[0]) * (t - p[1]) * (t - p[1])
                    + 2.0 * (t - p[0]) * (t - p[0]) * (t - p[1]);
            }
            _ => return Err(Fault::Degenerate),
        }
        if *f == 0.0 {
            *found = true;
        }
        Ok(())
    }
}

struct FamilyN {
    kind: i64,
    p: [f64; 4],
}

impl EvalNFunctions for FamilyN {
    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        jacobian: Option<&mut Matrix>,
        found: &mut bool,
    ) -> Result<(), Fault> {
        let p = self.p;
        *found = false;
        match self.kind {
            1 | 2 => {
                f[0] = x[0] * x[0] + x[1] * x[1] - p[0];
                f[1] = x[1] - p[1] * x[0] - p[2];
                if let Some(jacobian) = jacobian {
                    if self.kind == 1 {
                        jacobian[0][0] = 2.0 * x[0];
                        jacobian[0][1] = 2.0 * x[1];
                        jacobian[1][0] = -p[1];
                        jacobian[1][1] = 1.0;
                    } else {
                        let values = f.to_vec();
                        let mut xx = x.to_vec();
                        compute_jacobian(self, &mut xx, &values, jacobian, None)?;
                    }
                }
            }
            3 => {
                f[0] = (x[0] - p[0]) * (x[0] - p[0]) + x[1] * x[1] - 1.0;
                f[1] = (x[0] - p[1]) * (x[0] - p[1]) + x[1] * x[1] - 1.0;
                if let Some(jacobian) = jacobian {
                    jacobian[0][0] = 2.0 * (x[0] - p[0]);
                    jacobian[0][1] = 2.0 * x[1];
                    jacobian[1][0] = 2.0 * (x[0] - p[1]);
                    jacobian[1][1] = 2.0 * x[1];
                }
            }
            4 => {
                f[0] = x[0] + x[1] + x[2] - p[0];
                f[1] = 2.0 * x[0] + 2.0 * x[1] + 2.0 * x[2] - p[1];
                f[2] = x[0] - x[2] - p[2];
                if let Some(jacobian) = jacobian {
                    let rows = [[1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [1.0, 0.0, -1.0]];
                    for (i, row) in rows.iter().enumerate() {
                        jacobian[i].copy_from_slice(row);
                    }
                }
            }
            _ => return Err(Fault::Degenerate),
        }
        Ok(())
    }
}

fn params(input: &mut Input) -> [f64; 4] {
    std::array::from_fn(|_| input.d())
}

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    match op {
        "solve1" => {
            let family = Family1 {
                kind: input.i(),
                p: params(input),
            };
            let guess = input.d();
            let has_interval = input.i() != 0;
            let (lo, hi) = (input.d(), input.d());
            let periodic = input.i() != 0;
            let algorithm = [
                SolverAlgorithm::Newton,
                SolverAlgorithm::Bisection,
                SolverAlgorithm::BrentMaximize,
                SolverAlgorithm::BrentMinimize,
            ][input.i() as usize];
            let (acceptable, desired) = (input.d(), input.d());
            let mut solver = if has_interval {
                LocalSolve1d::with_interval(family, Extent1::new(lo, hi), periodic)
            } else {
                LocalSolve1d::new(family)
            };
            solver.set_solver_algorithm(algorithm);
            if desired > 0.0 {
                solver.set_desired_accuracy(desired);
            }
            if let Some(root) = out.status(solver.solve_it(guess, acceptable)) {
                out.flag(root.found);
                out.d(root.t);
                out.d(solver.found_accuracy());
                out.i(solver.termination_reason() as i64);
            }
        }
        "solven" => {
            let family = FamilyN {
                kind: input.i(),
                p: params(input),
            };
            let n = input.i() as usize;
            let guess: Vec<f64> = (0..n).map(|_| input.d()).collect();
            let has_intervals = input.i() != 0;
            let mut intervals = Vec::new();
            let mut periodic = Vec::new();
            for _ in 0..n {
                intervals.push(Extent1::new(input.d(), input.d()));
                periodic.push(input.i() != 0);
            }
            let handler = [
                BoundaryHandler::TotalBoundaryHits,
                BoundaryHandler::HitsAllBoundaries,
                BoundaryHandler::IgnoreBoundaryHits,
            ][input.i() as usize];
            let max_hits = input.i() as u64;
            let cutting = input.i() != 0;
            let acceptable = input.d();
            let mut solver = if has_intervals {
                LocalSolveNd::new(family, Some(ExtentN { intervals }), Some(periodic))
            } else {
                LocalSolveNd::new(family, None, None)
            };
            solver.set_boundary_handler(handler, max_hits);
            solver.set_step_cutting_enabled(cutting);
            if let Some(root) = out.status(solver.solve_it(&guess, acceptable)) {
                out.flag(root.found);
                for value in root.x {
                    out.d(value);
                }
                out.i(solver.termination_reason() as i64);
            }
        }
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn local_solvers_match_smlib() {
    let report = run_gate(
        include_str!("fixtures/local-solve/cases.txt"),
        |_| 4,
        evaluate,
    );
    report.assert_passed(&["tangent", "touching", "singular"]);
}
