#![allow(missing_docs)]
//! Differential gate for the curve construction and knot-editing routines
//! (`SmBSplineCurve::CreateLineSegment`, `CreateDegenerateCurve`,
//! `InsertOneKnot`, `InsertKnots` and `SmCurve::EditParameterization`)
//! against saved outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::extent::Extent1;
use catallax_trace::spline::{BSplineCurve, ControlPoint, NO_Z};
use common::{Input, Output, run_gate};

/// `dim degree count knots... (x y z w)...`, every value as stored.
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

/// The curve as SMLib's driver prints it: dimension, degree, counts, every
/// knot and every control point as stored.
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
        "line" => {
            let dimension = input.i() as usize;
            let (start, end) = (input.v3(), input.v3());
            if let Some(curve) = out.status(BSplineCurve::line_segment(dimension, start, end)) {
                emit(out, &curve);
            }
        }
        "point" | "point-planar" => {
            let dimension = input.i() as usize;
            let point = input.v3();
            if let Some(curve) = out.status(BSplineCurve::degenerate(dimension, point)) {
                emit(out, &curve);
            }
        }
        "insert" => {
            let curve = raw_curve(input);
            let parameter = input.d();
            let times = input.i() as usize;
            if let Some(curve) = out.status(curve.insert_knot(parameter, times)) {
                emit(out, &curve);
            }
        }
        "refine" | "refine-unordered" => {
            let curve = raw_curve(input);
            let count = input.i() as usize;
            let knots: Vec<f64> = (0..count).map(|_| input.d()).collect();
            if let Some(curve) = out.status(curve.refine(&knots)) {
                emit(out, &curve);
            }
        }
        "reparam" | "reparam-nonincreasing" => {
            let curve = raw_curve(input);
            let interval = Extent1::new(input.d(), input.d());
            if let Some(curve) = out.status(curve.reparameterized(interval)) {
                emit(out, &curve);
            }
        }
        other => panic!("unknown operation {other}"),
    }
}

const FIXTURE: &str = include_str!("fixtures/curve-build/cases.txt");

/// Operations whose SMLib results the checked curve storage cannot hold. Each
/// is pinned individually below, never compared and never dropped.
const PINNED: [&str; 3] = ["point-planar", "refine-unordered", "reparam-nonincreasing"];

fn is_pinned(line: &str) -> bool {
    PINNED.iter().any(|op| {
        line.strip_prefix(op)
            .is_some_and(|rest| rest.starts_with(' '))
    })
}

/// The fixture lines on one side of the pinned split, comments kept.
fn lines(pinned: bool) -> String {
    FIXTURE
        .lines()
        .filter(|line| line.starts_with('#') || is_pinned(line) == pinned)
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn curve_build_matches_smlib() {
    let report = run_gate(&lines(false), |_| 0, run);
    assert!(report.cases >= 600, "the gate lost cases: {}", report.cases);
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

/// One pinned fixture line: its operation, port outputs and SMLib outputs.
struct Pinned {
    op: String,
    port: Vec<Tok>,
    oracle: Vec<Tok>,
    inputs: String,
}

fn pinned_cases() -> Vec<Pinned> {
    FIXTURE
        .lines()
        .filter(|line| is_pinned(line))
        .map(|line| {
            let mut parts = line.splitn(3, '|');
            let op = parts.next().unwrap().split_whitespace().next().unwrap();
            let (inputs, outputs) = (parts.next().unwrap(), parts.next().unwrap());
            // The gate's candidate record lists each port output as a tag
            // (0 binary64, 1 integer) and its bits.
            let report = run_gate(line, |_| u64::MAX, run);
            let port = report
                .candidate_bits
                .chunks(2)
                .map(|pair| match pair[0] {
                    0 => Tok::D(f64::from_bits(pair[1])),
                    _ => Tok::I(pair[1] as i64),
                })
                .collect();
            Pinned {
                op: op.to_owned(),
                port,
                oracle: parse_tokens(outputs),
                inputs: inputs.to_owned(),
            }
        })
        .collect()
}

fn parse_tokens(text: &str) -> Vec<Tok> {
    text.split_whitespace()
        .map(|token| match token.split_at(1) {
            ("d", hex) => Tok::D(f64::from_bits(u64::from_str_radix(hex, 16).unwrap())),
            ("i", decimal) => Tok::I(decimal.parse().unwrap()),
            _ => panic!("bad output token {token}"),
        })
        .collect()
}

fn doubles(tokens: &[Tok]) -> Vec<f64> {
    tokens
        .iter()
        .filter_map(|token| match token {
            Tok::D(value) => Some(*value),
            Tok::I(_) => None,
        })
        .collect()
}

fn cases_of(op: &str, expected: usize) -> Vec<Pinned> {
    let cases: Vec<_> = pinned_cases()
        .into_iter()
        .filter(|case| case.op == op)
        .collect();
    assert_eq!(cases.len(), expected, "pinned {op} cases");
    cases
}

/// SMLib's `SmCurve::EditParameterization` accepts a single-valued interval
/// and returns a curve whose active knot interval is empty. A checked
/// `BSplineCurve` has a positive-length interval, so the port refuses.
#[test]
fn reparameterisation_to_a_single_value_is_refused() {
    for case in cases_of("reparam-nonincreasing", 24) {
        assert_eq!(case.port, [Tok::I(1007)], "{}", case.inputs);
        assert_eq!(case.oracle[0], Tok::I(1000));
        let (Tok::I(degree), Tok::I(count)) = (case.oracle[2], case.oracle[3]) else {
            unreachable!()
        };
        let knots = doubles(&case.oracle);
        let (degree, count) = (degree as usize, count as usize);
        assert_eq!(
            knots[degree].to_bits(),
            knots[count].to_bits(),
            "{}",
            case.inputs
        );
    }
}

/// `N_CrvRefine` documents an ascending knot list. Given a descending pair it
/// returns a curve with decreasing knots; the port refuses it.
#[test]
fn refinement_by_unordered_knots_is_refused() {
    for case in cases_of("refine-unordered", 16) {
        assert_eq!(case.port, [Tok::I(1007)], "{}", case.inputs);
        assert_eq!(case.oracle[0], Tok::I(1000));
        let Tok::I(total) = case.oracle[4] else {
            unreachable!()
        };
        let knots = &doubles(&case.oracle)[..total as usize];
        assert!(
            knots.windows(2).any(|pair| pair[0] > pair[1]),
            "{}",
            case.inputs
        );
    }
}

/// `CreateDegenerateCurve` builds a two-dimensional curve through
/// `CreateCanonical`, which keeps the point's stored z. The port stores the
/// two-dimensional marker. Every other output is bit-identical.
#[test]
fn planar_point_curves_store_the_planar_marker() {
    for case in cases_of("point-planar", 2) {
        assert_eq!(case.port.len(), case.oracle.len());
        // After the status and four integers: four knots, then two control
        // points of (x, y, z, w) whose third values are the z slots.
        let z_slots = [5 + 4 + 2, 5 + 4 + 6];
        for (index, (a, b)) in case.port.iter().zip(&case.oracle).enumerate() {
            if z_slots.contains(&index) {
                assert_eq!(*a, Tok::D(NO_Z));
                assert_ne!(*b, Tok::D(NO_Z));
            } else {
                match (a, b) {
                    (Tok::D(a), Tok::D(b)) => {
                        assert_eq!(a.to_bits(), b.to_bits(), "output {index}")
                    }
                    _ => assert_eq!(a, b, "output {index}"),
                }
            }
        }
    }
}
