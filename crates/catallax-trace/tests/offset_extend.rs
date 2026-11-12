#![allow(missing_docs)]
//! Differential gate for `SmOffsetSurface` with extension allowed: the
//! extension of the base surface, its keeping between evaluations, the
//! periodic wrapping and the stop at a singular side, against saved outputs of
//! the pinned C++ SMLib.

mod common;

use catallax_trace::spline::{BSplineSurface, ControlPoint};
use catallax_trace::surface_offset::OffsetSurface;
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

/// The extension as the oracle prints it: degrees, net and knot counts, knots, net.
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

fn run(_op: &str, input: &mut Input, out: &mut Output) {
    let distance = input.d();
    let allow = input.i() != 0;
    let surface = raw_surface(input);
    let count = input.i() as usize;
    let offset = OffsetSurface::new(&surface, distance).allow_extension(allow);
    out.flag(offset.is_periodic(0));
    out.flag(offset.is_periodic(1));
    out.i(i64::from(offset.singularities()));
    for _ in 0..count {
        let kind = input.i();
        let uv = Vector2::new(input.d(), input.d());
        match kind {
            0 => {
                if let Some(point) = out.status(offset.evaluate_point(uv)) {
                    out.v3(point);
                }
            }
            1 => {
                let sides = [input.i() != 0, input.i() != 0];
                out.i(1000);
                out.v3(offset.evaluate_normal(uv, sides));
            }
            _ => {
                let simple = input.i() != 0;
                let (hu, hv) = (input.i() as usize, input.i() as usize);
                let sides = [input.i() != 0, input.i() != 0];
                let upper = input.i() != 0;
                let non_zero = input.i() != 0;
                let cells = if simple {
                    offset.evaluate_simple(uv, hu, hv, sides, upper)
                } else {
                    offset.evaluate(uv, hu, hv, sides, upper, non_zero)
                };
                if let Some(cells) = out.status(cells) {
                    for cell in cells {
                        out.v3(cell);
                    }
                }
            }
        }
        let domain = offset.natural_domain();
        out.v2(domain.min);
        out.v2(domain.max);
        out.flag(offset.extended_surface().is_some());
    }
    if let Some(extended) = offset.extended_surface() {
        emit_surface(out, &extended);
    }
}

const FIXTURE: &str = include_str!("fixtures/offset-extend/cases.txt");

#[test]
fn offset_extension_matches_smlib() {
    let report = run_gate(FIXTURE, |_| 0, run);
    assert!(
        report.cases >= 1000,
        "the gate lost cases: {}",
        report.cases
    );
    report.assert_passed(&["regular"]);
}
