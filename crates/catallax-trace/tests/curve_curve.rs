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
        "global" | "grefuse" | "known-miss" => {
            match global_curve_intersect(&a, ia, &b, ib, tolerance, flag) {
                Ok(found) if op == "global" || op == "known-miss" => {
                    out.i(common::SUCCESS);
                    solutions(out, &found);
                }
                Err(SolverRefusal::NotPorted(_)) if op == "grefuse" => out.i(-1),
                Err(SolverRefusal::Failed(fault)) if op == "global" => out.i(status_code(fault)),
                // Mismatched expectations show as a failing case.
                Ok(_) => out.i(-2),
                Err(_) => out.i(-3),
            }
        }
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
    let fixture = include_str!("fixtures/curve-curve/cases.txt");
    let ordinary = fixture
        .lines()
        .filter(|line| !line.starts_with("known-miss "))
        .collect::<Vec<_>>()
        .join("\n");
    let report = run_gate(&ordinary, |_| 16, evaluate);
    assert_eq!(report.cases, 460);
    report.assert_passed(&[
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
    common::assert_same_with_workers(
        include_str!("fixtures/curve-curve/cases.txt"),
        |_| 16,
        evaluate,
    );
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

fn analytic_segment(start: [f64; 3], finish: [f64; 3]) -> BSplineCurve {
    BSplineCurve::new(
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        [start, finish]
            .map(|p| ControlPoint {
                coordinates: [p[0], p[1], p[2], NO_WEIGHT],
            })
            .to_vec(),
    )
    .unwrap()
}

#[test]
fn analytic_line_pair_crossing() {
    let a = analytic_segment([0.0, 0.0, 0.0], [2.0, 1.0, 0.0]);
    let b = analytic_segment([0.0, 1.0, 0.0], [2.0, 0.0, 0.0]);
    let found = global_curve_intersect(
        &a,
        Extent1::new(0.0, 1.0),
        &b,
        Extent1::new(0.0, 1.0),
        1e-5,
        false,
    )
    .expect("the analytic line pair is ported");
    assert_eq!(found.len(), 1);
    let point = found.as_slice()[0];
    assert_eq!(point.kind, SolutionType::SingleValue);
    assert_eq!(point.variables, 2);
    assert_eq!(point.object_count, 2);
    for parameter in &point.start.parameters[..2] {
        assert!(common::ulps(*parameter, 0.5) <= 2);
    }
    assert_eq!(point.start.value, 0.0);
}

#[test]
fn analytic_line_range_keeps_reverse_parameters() {
    let a = analytic_segment([0.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
    let b = analytic_segment([3.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
    let found = global_curve_intersect(
        &a,
        Extent1::new(0.0, 1.0),
        &b,
        Extent1::new(0.0, 1.0),
        1e-5,
        false,
    )
    .unwrap();
    assert_eq!(found.len(), 1);
    let range = found.as_slice()[0];
    assert_eq!(range.kind, SolutionType::RangeOfValues);
    assert_eq!(&range.start.parameters[..2], &[0.5, 1.0]);
    assert_eq!(&range.end.parameters[..2], &[1.0, 0.5]);
    assert_eq!(range.start.value, 0.0);
    assert_eq!(range.end.value, 0.0);
}

#[test]
fn smlibs_small_crossings_are_known_misses() {
    let fixture = include_str!("fixtures/curve-curve/cases.txt");
    let misses = fixture
        .lines()
        .filter(|line| line.starts_with("known-miss "))
        .collect::<Vec<_>>()
        .join("\n");
    let report = run_gate(&misses, |_| 16, evaluate);
    assert_eq!(
        report.cases, 8,
        "known misses must stay outside the correctness count"
    );
    assert!(report.failures.is_empty(), "{:?}", report.failures);
    for length in [1e-10, 1e-8] {
        let a = analytic_segment([0.0, 0.0, 0.0], [length, 0.0, 0.0]);
        let b = analytic_segment([length / 2.0, -length, 0.0], [length / 2.0, length, 0.0]);
        // Both affine curves meet at the same dyadic midpoint, despite the
        // candidate helper's absolute cross-product threshold returning none.
        assert_eq!(
            a.evaluate_point(0.5).unwrap(),
            b.evaluate_point(0.5).unwrap()
        );
        for (first, second) in [(&a, &b), (&b, &a)] {
            for skip in [false, true] {
                let candidates = global_curve_intersect(
                    first,
                    Extent1::new(0.0, 1.0),
                    second,
                    Extent1::new(0.0, 1.0),
                    1e-12,
                    skip,
                )
                .unwrap();
                assert!(candidates.is_empty(), "the known SMLib miss changed");
            }
        }
    }
}
