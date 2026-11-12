#![allow(missing_docs)]
//! Differential gate for NLib's surface knot removal and surface extension
//! against saved outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::spline::{BSplineSurface, ControlPoint};
use catallax_trace::surface_extend::{Closed, evaluate_unbounded, remove_knot_conditional};
use catallax_trace::vector::Vector2;
use common::{Input, Output, run_gate};

fn raw_surface(input: &mut Input) -> BSplineSurface {
    let degrees = [input.i() as usize, input.i() as usize];
    let counts = [input.i() as usize, input.i() as usize];
    let _rational = input.i();
    let knots = [0, 1].map(|axis| {
        (0..counts[axis] + degrees[axis] + 1)
            .map(|_| input.d())
            .collect()
    });
    let points = (0..counts[0] * counts[1])
        .map(|_| ControlPoint {
            coordinates: [input.d(), input.d(), input.d(), input.d()],
        })
        .collect();
    BSplineSurface::new(degrees, knots, counts, points).expect("valid test surface")
}

/// The surface as the oracle prints it: degrees, net and knot counts, knots, net.
fn emit_surface(out: &mut Output, surface: &BSplineSurface) {
    let (degrees, counts, knots) = (surface.degrees(), surface.counts(), surface.knots());
    out.i(degrees[0] as i64);
    out.i(degrees[1] as i64);
    out.i(counts[0] as i64);
    out.i(counts[1] as i64);
    out.i(knots[0].len() as i64);
    out.i(knots[1].len() as i64);
    for axis_knots in knots {
        for &knot in axis_knots {
            out.d(knot);
        }
    }
    for point in surface.control_points() {
        for coordinate in point.coordinates {
            out.d(coordinate);
        }
    }
}

fn run(op: &str, input: &mut Input, out: &mut Output) {
    match op {
        "srf-remove" | "srf-remove-overrun" => {
            let surface = raw_surface(input);
            let axis = input.i() as usize;
            let (t, times, tolerance) = (input.d(), input.i() as usize, input.d());
            match remove_knot_conditional(&surface, axis, t, times, tolerance) {
                Ok(done) => {
                    out.i(i64::from(done.flagged));
                    out.i(done.removed as i64);
                    emit_surface(out, &done.surface);
                }
                // Where the port refuses, SMLib's surface is as it was.
                Err(_) => {
                    out.i(1);
                    out.i(0);
                    emit_surface(out, &surface);
                }
            }
        }
        "srf-extend" => {
            let surface = raw_surface(input);
            let count = input.i() as usize;
            let mut extension = None;
            for _ in 0..count {
                let uv = Vector2::new(input.d(), input.d());
                match evaluate_unbounded(&surface, &mut extension, uv) {
                    Ok(point) => {
                        out.i(0);
                        out.v3(point);
                    }
                    Err(_) => out.i(1),
                }
                match &extension {
                    None => out.i(0),
                    Some(extended) => {
                        out.i(1);
                        out.i(match extended.closed {
                            Closed::No => 0,
                            Closed::U => 1,
                            Closed::V => 2,
                            Closed::UV => 3,
                        });
                        emit_surface(out, &extended.surface().expect("a valid extension"));
                    }
                }
            }
        }
        other => panic!("unknown operation {other}"),
    }
}

const FIXTURE: &str = include_str!("fixtures/surface-extend/cases.txt");

fn is_case_of(line: &str, name: &str) -> bool {
    line.strip_prefix(name)
        .is_some_and(|rest| rest.starts_with(' '))
}

fn compared() -> String {
    FIXTURE
        .lines()
        .filter(|line| !is_case_of(line, "srf-remove-overrun"))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn surface_extend_matches_smlib() {
    let report = run_gate(&compared(), |_| 0, run);
    assert!(report.cases >= 800, "the gate lost cases: {}", report.cases);
    report.assert_passed(&[
        "regular",
        "near-singular",
        "singular",
        "boundary",
        "touching",
    ]);
}

/// SMLib's removal sizes its scratch rows from the larger degree and reads and
/// writes past them when the span reaches further, as it does for the
/// surfaces of degree 1 and 2 here; it dies, so there is no answer to compare.
/// The port refuses.
#[test]
fn removals_that_overrun_smlib_scratch_are_refused() {
    let cases: Vec<&str> = FIXTURE
        .lines()
        .filter(|line| is_case_of(line, "srf-remove-overrun"))
        .collect();
    assert_eq!(cases.len(), OVERRUN_CASES);
    for line in cases {
        assert!(line.ends_with("| i-1"), "{line}");
        // `run` emits `1, 0, surface...` for a refusal; the first two tokens say so.
        let report = run_gate(line, |_| u64::MAX, run);
        assert_eq!(&report.candidate_bits[..4], &[1, 1, 1, 0], "{line}");
    }
}

const OVERRUN_CASES: usize = 25;
