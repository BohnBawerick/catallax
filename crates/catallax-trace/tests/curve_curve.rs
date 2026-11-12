#![allow(missing_docs)]
//! Differential gate for SMLib's global curve-curve intersector
//! (`SmCurve::GlobalCurveIntersect` with `SmCCIGlobalSolver`) and its
//! whole-curve coincidence check (`SmCurve::GlobalCoincidenceChecker`)
//! against saved outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::Fault;
use catallax_trace::curve_curve::{global_coincidence_check, global_curve_intersect};
use catallax_trace::extent::Extent1;
use catallax_trace::global_solver::SolverRefusal;
use catallax_trace::solution::{SolutionArray, SolutionType};
use catallax_trace::spline::{BSplineCurve, ControlPoint, NO_WEIGHT};
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

// Review witnesses for #454.

/// A degree-1 curve whose two controls sit at the same point, over [0, 1].
fn point_curve(point: [f64; 3]) -> BSplineCurve {
    BSplineCurve::new(
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        [point, point]
            .map(|[x, y, z]| ControlPoint {
                coordinates: [x, y, z, NO_WEIGHT],
            })
            .to_vec(),
    )
    .unwrap()
}

#[test]
fn a_degenerate_curve_keeps_smlibs_empty_coincidence_answer() {
    let point = point_curve([0.0, 0.0, 0.0]);
    let line = BSplineCurve::new(
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        [(-1.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
            .map(|(x, y, z)| ControlPoint {
                coordinates: [x, y, z, NO_WEIGHT],
            })
            .to_vec(),
    )
    .unwrap();
    let interval = Extent1::new(0.0, 1.0);
    // The limitation documented on `CoincidenceCheck::more_work`, pinned as
    // the curve-curve gate pins SMLib's `global singular` answer: a
    // degenerate curve ends the coincidence check without a search, even
    // when the point curve meets the other curve at the origin.
    let check = global_coincidence_check(&point, interval, &line, interval, 1.0e-6, false).unwrap();
    assert!(!check.more_work);
    assert!(check.solutions.is_empty());
    let found = global_curve_intersect(&point, interval, &line, interval, 1.0e-6, false).unwrap();
    assert!(found.is_empty());
}

#[test]
fn coincidence_refuses_a_subdivision_that_makes_no_progress() {
    // 2^52 + 1 and 2^52 + 2 round midpoint evaluation back to the ends.
    let a = 4_503_599_627_370_497.0;
    let b = a + 1.0;
    let curve = BSplineCurve::new(
        2,
        vec![a, a, a, b, b, b],
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 0.0)]
            .map(|(x, y, z)| ControlPoint {
                coordinates: [x, y, z, NO_WEIGHT],
            })
            .to_vec(),
    )
    .unwrap();
    let interval = Extent1::new(a, b);
    assert!(matches!(
        global_coincidence_check(&curve, interval, &curve.clone(), interval, 1.0e-6, false),
        Err(SolverRefusal::Failed(Fault::InvalidInput))
    ));
}
