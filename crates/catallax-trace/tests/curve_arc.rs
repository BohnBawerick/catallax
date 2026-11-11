#![allow(missing_docs)]
//! Differential gate for `SmBSplineCurve::CreateCircleSegment` (the
//! quadratic arcs and the quintic circle) and `EditEndPoint` against saved
//! outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::curve_arc::CircleParameterization;
use catallax_trace::curve_nlib::{
    closest_point_on_curve, curve_derivatives, extract_curve_segment,
};
use catallax_trace::spline::{BSplineCurve, ControlPoint};
use common::{Input, Output, run_gate};

fn raw_curve(input: &mut Input) -> BSplineCurve {
    let _dimension = input.i();
    let degree = input.i() as usize;
    let count = input.i() as usize;
    let knots = (0..count + degree + 1).map(|_| input.d()).collect();
    let points = (0..count)
        .map(|_| ControlPoint {
            coordinates: [input.d(), input.d(), input.d(), input.d()],
        })
        .collect();
    BSplineCurve::new(degree, knots, points).expect("valid test curve")
}

fn emit(out: &mut Output, curve: &BSplineCurve) {
    out.i(curve.dimension() as i64);
    out.i(curve.degree() as i64);
    out.i(curve.control_points().len() as i64);
    out.i(curve.knots().len() as i64);
    for &knot in curve.knots() {
        out.d(knot);
    }
    for point in curve.control_points() {
        for value in point.coordinates {
            out.d(value);
        }
    }
}

fn run(op: &str, input: &mut Input, out: &mut Output) {
    match op {
        "circle" => {
            let dimension = input.i() as usize;
            let (origin, x, y) = (input.v3(), input.v3(), input.v3());
            let (radius, start, end) = (input.d(), input.d(), input.d());
            let parameterization = if input.i() == 0 {
                CircleParameterization::Quadratic
            } else {
                CircleParameterization::Quintic
            };
            let circle = BSplineCurve::circle_segment(
                dimension,
                origin,
                x,
                y,
                radius,
                start,
                end,
                parameterization,
            );
            if let Some(curve) = out.status(circle) {
                emit(out, &curve);
            }
        }
        "edit-end" | "edit-end-invalid" => {
            let curve = raw_curve(input);
            let point = input.v3();
            let start = input.i() != 0;
            let weighted = input.i() != 0;
            let weight = input.d();
            let edited = curve.edit_end_point(point, start, weighted.then_some(weight));
            if let Some(curve) = out.status(edited) {
                emit(out, &curve);
            }
        }
        "derivs" => {
            let curve = raw_curve(input);
            let (u, left, order) = (input.d(), input.i() != 0, input.i() as usize);
            if let Some(values) = out.status(curve_derivatives(&curve, u, left, order)) {
                for value in values {
                    out.v3(value);
                }
            }
        }
        "closest" => {
            let curve = raw_curve(input);
            let point = input.v3();
            let (start, distance, cosine) = (input.d(), input.d(), input.d());
            if let Some((u, found)) = out.status(closest_point_on_curve(
                &curve, point, start, distance, cosine,
            )) {
                out.d(u);
                out.v3(found);
            }
        }
        "extract" => {
            let curve = raw_curve(input);
            let (left, right) = (input.d(), input.d());
            if let Some(segment) = out.status(extract_curve_segment(&curve, left, right)) {
                emit(out, &segment);
            }
        }
        other => panic!("unknown operation {other}"),
    }
}

/// Binary64 outputs that differ from SMLib in their last bits: `libm`'s sine
/// and cosine against the C library's.
const INEXACT_VALUES: usize = 23;

const FIXTURE: &str = include_str!("fixtures/curve-arc/cases.txt");

/// Operations whose SMLib results the port does not reproduce. Each is pinned
/// individually below, never compared and never dropped.
const PINNED: [&str; 1] = ["edit-end-invalid"];

fn is_pinned(line: &str) -> bool {
    PINNED.iter().any(|op| {
        line.strip_prefix(op)
            .is_some_and(|rest| rest.starts_with(' '))
    })
}

fn lines(pinned: bool) -> String {
    FIXTURE
        .lines()
        .filter(|line| line.starts_with('#') || is_pinned(line) == pinned)
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn curve_arc_matches_smlib() {
    // Sines and cosines come from `libm`; SMLib's come from the C library.
    let report = run_gate(&lines(false), |_| 4, run);
    assert!(report.cases >= 270, "the gate lost cases: {}", report.cases);
    // Every value that is not bit-identical comes from a sine or cosine; the
    // count is pinned so that no other difference can join them unseen.
    assert_eq!(
        report.inexact.len(),
        INEXACT_VALUES,
        "{:#?}",
        report.inexact
    );
    assert!(report.max_ulps <= 3);
    report.assert_passed(&[
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}

/// One output token: an integer or binary64 value.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Tok {
    I(i64),
    D(f64),
}

/// The port's and SMLib's outputs for each pinned line of `op`.
fn pinned(op: &str, expected: usize) -> Vec<(Vec<Tok>, Vec<Tok>, String)> {
    let cases: Vec<_> = FIXTURE
        .lines()
        .filter(|line| {
            line.strip_prefix(op)
                .is_some_and(|rest| rest.starts_with(' '))
        })
        .map(|line| {
            let mut parts = line.splitn(3, '|');
            let _ = parts.next();
            let (inputs, outputs) = (parts.next().unwrap(), parts.next().unwrap());
            let report = run_gate(line, |_| u64::MAX, run);
            let port = report
                .candidate_bits
                .chunks(2)
                .map(|pair| match pair[0] {
                    0 => Tok::D(f64::from_bits(pair[1])),
                    _ => Tok::I(pair[1] as i64),
                })
                .collect();
            let oracle = outputs
                .split_whitespace()
                .map(|token| match token.split_at(1) {
                    ("d", hex) => Tok::D(f64::from_bits(u64::from_str_radix(hex, 16).unwrap())),
                    ("i", decimal) => Tok::I(decimal.parse().unwrap()),
                    _ => panic!("bad output token {token}"),
                })
                .collect();
            (port, oracle, inputs.to_owned())
        })
        .collect();
    assert_eq!(cases.len(), expected, "pinned {op} cases");
    cases
}

/// SMLib writes a zero weight or a non-finite coordinate into an end control
/// point and succeeds; a checked `BSplineCurve` refuses both.
#[test]
fn end_points_a_checked_curve_cannot_hold_are_refused() {
    for (port, oracle, inputs) in pinned("edit-end-invalid", 6) {
        assert_eq!(port, [Tok::I(1007)], "{inputs}");
        assert_eq!(oracle[0], Tok::I(1000), "{inputs}");
    }
}
