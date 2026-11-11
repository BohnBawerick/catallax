#![allow(missing_docs)]
//! Differential gate for `SmBSplineCurve::CreateInterpolatingCurve` through
//! points with a vector at each (`N_FitCrvFirstDeriv`) against saved outputs
//! of the pinned C++ SMLib.

mod common;

use catallax_trace::curve_fit::{Parameterization, fair_blend_curve_derivatives};
use catallax_trace::spline::{BSplineCurve, NO_Z};
use common::{Input, Output, run_gate};

fn run(op: &str, input: &mut Input, out: &mut Output) {
    if op == "fair" {
        let points = [input.v3(), input.v3()];
        let tangents = [input.v3(), input.v3()];
        let higher: Vec<_> = (0..input.i()).map(|_| input.v3()).collect();
        let prevent = input.i() != 0;
        if let Some(faired) = out.status(fair_blend_curve_derivatives(
            points, tangents, &higher, prevent,
        )) {
            for tangent in faired.tangents {
                out.v3(tangent);
            }
            out.i(faired.higher.len() as i64);
            for vector in faired.higher {
                out.v3(vector);
            }
        }
        return;
    }
    assert!(matches!(op, "interpolate" | "interpolate-planar"));
    let parameterization = match input.i() {
        0 => Parameterization::Uniform,
        1 => Parameterization::ChordLength,
        _ => Parameterization::Centripetal,
    };
    let dimension = input.i() as usize;
    let degree = input.i() as usize;
    let points: Vec<_> = (0..input.i()).map(|_| input.v3()).collect();
    let vectors: Vec<_> = (0..input.i()).map(|_| input.v3()).collect();
    let ignore = input.i() != 0;
    let has_more = input.i() != 0;
    let more: Vec<_> = (0..input.i()).map(|_| input.v3()).collect();
    let result = BSplineCurve::interpolating_curve(
        parameterization,
        dimension,
        degree,
        &points,
        &vectors,
        has_more.then_some(&more[..]),
        ignore,
    );
    if let Some(result) = out.status(result) {
        let curve = &result.curve;
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
        out.i(result.vectors.len() as i64);
        for vector in &result.vectors {
            out.v3(*vector);
        }
    }
}

const FIXTURE: &str = include_str!("fixtures/curve-fit/cases.txt");
const PINNED: [&str; 1] = ["interpolate-planar"];

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
fn curve_fit_matches_smlib() {
    // The fit is plain arithmetic and a banded solve: no budget.
    let report = run_gate(&lines(false), |_| 0, run);
    assert!(
        report.cases >= 1200,
        "the gate lost cases: {}",
        report.cases
    );
    report.assert_passed(&[
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Tok {
    I(i64),
    D(f64),
}

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

/// SMLib builds a two-dimensional curve from three-dimensional data and
/// labels it two-dimensional, keeping the points' z. The port stores the
/// planar marker. Every other value, and the vectors SMLib leaves, are
/// bit-identical.
#[test]
fn planar_curves_store_the_planar_marker() {
    for (port, oracle, inputs) in pinned("interpolate-planar", PLANAR_CASES) {
        assert_eq!(port.len(), oracle.len(), "{inputs}");
        assert_eq!(oracle[0], Tok::I(1000));
        let (Tok::I(points), Tok::I(knots)) = (oracle[3], oracle[4]) else {
            unreachable!()
        };
        let first_point = 5 + knots as usize;
        let z_slots: Vec<usize> = (0..points as usize)
            .map(|i| first_point + 4 * i + 2)
            .collect();
        for (index, (a, b)) in port.iter().zip(&oracle).enumerate() {
            if z_slots.contains(&index) {
                assert_eq!(*a, Tok::D(NO_Z), "{inputs}");
            } else if let (Tok::D(a), Tok::D(b)) = (a, b) {
                assert_eq!(a.to_bits(), b.to_bits(), "{inputs} output {index}");
            } else {
                assert_eq!(a, b, "{inputs} output {index}");
            }
        }
    }
}

const PLANAR_CASES: usize = 60;

/// A degree above NLib's limit of 17 reaches no fixed-size buffer: the
/// zero-vector and one-vector fits refuse it, with enough points that the
/// count test alone would let it through.
#[test]
fn degrees_above_the_nlib_limit_are_refused() {
    use catallax_trace::Fault;
    use catallax_trace::vector::Vector3;
    let points: Vec<Vector3> = (0..40)
        .map(|i| Vector3::new(f64::from(i), f64::from(i * i % 7), 0.5 * f64::from(i)))
        .collect();
    let start = [Vector3::new(1.0, 0.0, 0.0)];
    for degree in [18, 19, 31, 32] {
        for vectors in [&[][..], &start[..]] {
            let result = BSplineCurve::interpolating_curve(
                Parameterization::Uniform,
                3,
                degree,
                &points,
                vectors,
                None,
                false,
            );
            assert_eq!(result.err(), Some(Fault::Degenerate), "degree {degree}");
        }
    }
    // The limit itself still works.
    let result = BSplineCurve::interpolating_curve(
        Parameterization::Uniform,
        3,
        17,
        &points,
        &[],
        None,
        false,
    );
    assert!(result.is_ok());
}
