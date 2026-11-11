#![allow(missing_docs)]
//! Differential gate for degree elevation, knot removal, knot-vector merging
//! and `SmBSplineCurve::SyncronizeKnotsOfCurves` against saved outputs of
//! the pinned C++ SMLib.

mod common;

use catallax_trace::curve_compat::{compatible_knot_vectors, synchronize_knots};
use catallax_trace::spline::{BSplineCurve, ControlPoint};
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

/// The curve as SMLib's driver prints it.
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
        "elevate" => {
            let curve = raw_curve(input);
            let times = input.i() as usize;
            if let Some(curve) = out.status(curve.elevate_degree(times)) {
                emit(out, &curve);
            }
        }
        "remove" => {
            let curve = raw_curve(input);
            let parameter = input.d();
            let times = input.i() as usize;
            let tolerance = input.d();
            if let Some((curve, removed)) =
                out.status(curve.remove_knot(parameter, times, tolerance))
            {
                out.i(removed as i64);
                emit(out, &curve);
            }
        }
        "knots" => {
            let count = input.i() as usize;
            let tolerance = input.d();
            let vectors: Vec<Vec<f64>> = (0..count)
                .map(|_| {
                    let length = input.i() as usize;
                    (0..length).map(|_| input.d()).collect()
                })
                .collect();
            if let Some(merged) = out.status(compatible_knot_vectors(&vectors, tolerance)) {
                for vector in merged {
                    out.i(vector.adjusted.len() as i64);
                    for knot in vector.adjusted {
                        out.d(knot);
                    }
                    out.i(vector.additions.len() as i64);
                    for knot in vector.additions {
                        out.d(knot);
                    }
                }
            }
        }
        "sync" | "sync-mixed-dimension" => {
            let count = input.i() as usize;
            let tolerance = input.d();
            let curves: Vec<BSplineCurve> = (0..count).map(|_| raw_curve(input)).collect();
            if let Some(curves) = out.status(synchronize_knots(&curves, tolerance)) {
                for curve in &curves {
                    emit(out, curve);
                }
            }
        }
        other => panic!("unknown operation {other}"),
    }
}

const FIXTURE: &str = include_str!("fixtures/curve-compat/cases.txt");

/// The operation whose SMLib dimension labels the checked curves cannot keep.
const PINNED: &str = "sync-mixed-dimension";

fn is_pinned(line: &str) -> bool {
    line.strip_prefix(PINNED)
        .is_some_and(|rest| rest.starts_with(' '))
}

fn lines(pinned: bool) -> String {
    FIXTURE
        .lines()
        .filter(|line| line.starts_with('#') || is_pinned(line) == pinned)
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn curve_compat_matches_smlib() {
    let report = run_gate(&lines(false), |_| 0, run);
    assert!(report.cases >= 620, "the gate lost cases: {}", report.cases);
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

fn parse_tokens(text: &str) -> Vec<Tok> {
    text.split_whitespace()
        .map(|token| match token.split_at(1) {
            ("d", hex) => Tok::D(f64::from_bits(u64::from_str_radix(hex, 16).unwrap())),
            ("i", decimal) => Tok::I(decimal.parse().unwrap()),
            _ => panic!("bad output token {token}"),
        })
        .collect()
}

/// SMLib's `SyncronizeKnotsOfCurves` makes a planar curve three-dimensional
/// when it is synchronised with a spatial one (z becomes zero) but keeps the
/// curve's dimension label of 2. A `BSplineCurve` has no separate label: its
/// dimension is whether its points carry z. Every value is bit-identical;
/// only the label of each planar input differs, and the port reports 3.
#[test]
fn a_planar_curve_synchronised_with_a_spatial_one_reports_three_dimensions() {
    let cases: Vec<&str> = FIXTURE.lines().filter(|line| is_pinned(line)).collect();
    assert_eq!(cases.len(), 3);
    for line in cases {
        let mut parts = line.splitn(3, '|');
        let _ = parts.next();
        let (inputs, outputs) = (parts.next().unwrap(), parts.next().unwrap());
        let report = run_gate(line, |_| u64::MAX, run);
        let port: Vec<Tok> = report
            .candidate_bits
            .chunks(2)
            .map(|pair| match pair[0] {
                0 => Tok::D(f64::from_bits(pair[1])),
                _ => Tok::I(pair[1] as i64),
            })
            .collect();
        let oracle = parse_tokens(outputs);
        assert_eq!(port.len(), oracle.len(), "{inputs}");
        // The input dimensions, in order.
        let mut tokens = inputs.split_whitespace();
        let count: usize = tokens.next().unwrap().parse().unwrap();
        let mut rest: Vec<&str> = tokens.skip(1).collect();
        let mut input_dimensions = Vec::new();
        for _ in 0..count {
            let (dimension, degree, points): (usize, usize, usize) = (
                rest[0].parse().unwrap(),
                rest[1].parse().unwrap(),
                rest[2].parse().unwrap(),
            );
            input_dimensions.push(dimension);
            rest = rest.split_off(3 + points + degree + 1 + 4 * points);
        }
        assert!(input_dimensions.contains(&2) && input_dimensions.contains(&3));
        assert_eq!(port[0], Tok::I(1000));
        assert_eq!(oracle[0], Tok::I(1000));
        let mut at = 1;
        for dimension in input_dimensions {
            let Tok::I(degree) = oracle[at + 1] else {
                unreachable!()
            };
            let Tok::I(points) = oracle[at + 2] else {
                unreachable!()
            };
            let Tok::I(knots) = oracle[at + 3] else {
                unreachable!()
            };
            let _ = (degree, points);
            assert_eq!(oracle[at], Tok::I(dimension as i64));
            assert_eq!(port[at], Tok::I(3));
            let length = 4 + knots as usize + 4 * points as usize;
            for index in at + 1..at + length {
                match (port[index], oracle[index]) {
                    (Tok::D(a), Tok::D(b)) => {
                        assert_eq!(a.to_bits(), b.to_bits(), "output {index}")
                    }
                    (a, b) => assert_eq!(a, b, "output {index}"),
                }
            }
            at += length;
        }
        assert_eq!(at, port.len());
    }
}

/// Counts beyond the signed range are refused (elevation) or clamped to the
/// knot's multiplicity (removal) without panicking.
#[test]
fn oversized_counts_do_not_panic() {
    use catallax_trace::Fault;
    let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0];
    let points = (0..5)
        .map(|i| ControlPoint {
            coordinates: [f64::from(i), f64::from(i * i), 0.0, -f64::MAX],
        })
        .collect();
    let curve = BSplineCurve::new(2, knots, points).unwrap();
    for times in [usize::MAX, usize::MAX - 1, i64::MAX as usize, 16] {
        assert_eq!(
            curve.elevate_degree(times),
            Err(Fault::Degenerate),
            "{times}"
        );
    }
    let (_, removed) = curve.remove_knot(1.0, usize::MAX, 1e9).unwrap();
    assert!(removed <= 2);
}

/// Runs `work` on a thread and fails the test if it does not finish, so that
/// a regression to SMLib's endless loops fails instead of hanging.
fn finishes<T: Send + 'static>(work: impl FnOnce() -> T + Send + 'static) -> T {
    let (sender, receiver) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = sender.send(work());
    });
    receiver
        .recv_timeout(std::time::Duration::from_secs(20))
        .expect("the call did not finish")
}

/// SMLib's knot merge makes no progress for a tolerance that is zero,
/// negative or NaN, or for non-finite knots; the port refuses them.
#[test]
fn the_knot_merge_refuses_what_cannot_progress() {
    use catallax_trace::Fault;
    for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY, -0.0] {
        let result =
            finishes(move || compatible_knot_vectors(&[vec![0.0, 0.0, 1.0, 1.0]], tolerance));
        assert_eq!(result, Err(Fault::InvalidInput), "tolerance {tolerance}");
    }
    let result = finishes(|| compatible_knot_vectors(&[vec![0.0, f64::NAN, 1.0, 1.0]], 1e-8));
    assert_eq!(result, Err(Fault::InvalidInput));
    let result = finishes(|| compatible_knot_vectors(&[], 1e-8));
    assert_eq!(result, Err(Fault::InvalidInput));
}

/// A knot removal whose basis function has a support too short to step
/// across at its magnitude ends with SMLib's error flag instead of looping.
#[test]
fn a_tiny_support_does_not_stall_knot_removal() {
    let eps = f64::EPSILON;
    let knots = vec![
        0.0,
        0.0,
        0.0,
        1.0,
        1.0 + eps,
        1.0 + 2.0 * eps,
        1.0 + 3.0 * eps,
        2.0,
        2.0,
        2.0,
    ];
    let points = (0..7)
        .map(|i| ControlPoint {
            coordinates: [f64::from(i), f64::from(i * i), 0.0, -f64::MAX],
        })
        .collect();
    let curve = BSplineCurve::new(2, knots, points).unwrap();
    for knot in [1.0 + eps, 1.0 + 2.0 * eps, 1.0] {
        let curve = curve.clone();
        let _ = finishes(move || curve.remove_knot(knot, 1, 1e300));
    }
}
