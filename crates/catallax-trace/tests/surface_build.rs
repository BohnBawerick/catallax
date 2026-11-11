#![allow(missing_docs)]
//! Differential gate for `SmBSplineSurface::CreateRuledSurface` and
//! `Reparameterize` against saved outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::extent::Extent2;
use catallax_trace::spline::{BSplineCurve, BSplineSurface, ControlPoint};
use catallax_trace::surface_build::SurfaceParameter;
use catallax_trace::vector::Vector2;
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

/// `du dv nu nv rational knots_u... knots_v... (x y z w)...`, as stored.
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

fn emit(out: &mut Output, surface: &BSplineSurface) {
    let (degrees, counts, knots) = (surface.degrees(), surface.counts(), surface.knots());
    out.i(degrees[0] as i64);
    out.i(degrees[1] as i64);
    out.i(counts[0] as i64);
    out.i(counts[1] as i64);
    out.i(knots[0].len() as i64);
    out.i(knots[1].len() as i64);
    for &knot in knots[0].iter().chain(knots[1]) {
        out.d(knot);
    }
    for point in surface.control_points() {
        for value in point.coordinates {
            out.d(value);
        }
    }
}

fn run(op: &str, input: &mut Input, out: &mut Output) {
    match op {
        "ruled" => {
            let direction = if input.i() == 0 {
                SurfaceParameter::U
            } else {
                SurfaceParameter::V
            };
            let (first, second) = (raw_curve(input), raw_curve(input));
            if let Some(surface) = out.status(BSplineSurface::ruled(&first, &second, direction)) {
                emit(out, &surface);
            }
        }
        "revolve" => {
            let generator = raw_curve(input);
            let (origin, axis) = (input.v3(), input.v3());
            let angle = input.d();
            let swap = input.i() != 0;
            let surface =
                BSplineSurface::surface_of_revolution(&generator, origin, axis, angle, swap);
            if let Some(surface) = out.status(surface) {
                emit(out, &surface);
            }
        }
        "reparam" | "reparam-single" => {
            let surface = raw_surface(input);
            let (a, b, c, d) = (input.d(), input.d(), input.d(), input.d());
            let domain = Extent2::new(Vector2::new(a, c), Vector2::new(b, d));
            if let Some(surface) = out.status(surface.reparameterized(domain)) {
                emit(out, &surface);
            }
        }
        other => panic!("unknown operation {other}"),
    }
}

/// Binary64 outputs that differ from SMLib in their last bits.
const INEXACT_VALUES: usize = 210;

const FIXTURE: &str = include_str!("fixtures/surface-build/cases.txt");
const PINNED: &str = "reparam-single";

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
fn surface_build_matches_smlib() {
    // The circles of a revolution take their sines and cosines from `libm`,
    // SMLib from the C library; only `revolve` may differ in the last bits.
    let report = run_gate(&lines(false), |op| if op == "revolve" { 4 } else { 0 }, run);
    assert!(report.cases >= 570, "the gate lost cases: {}", report.cases);
    assert_eq!(report.inexact.len(), INEXACT_VALUES);
    assert!(report.inexact.iter().all(|line| line.contains(" revolve ")));
    assert!(report.max_ulps <= 2);
    report.assert_passed(&[
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}

/// SMLib reparameterises a surface to a single-valued range and returns a
/// surface whose active knot interval is empty; a checked `BSplineSurface`
/// has a positive-length domain, so the port refuses.
#[test]
fn a_single_valued_parameter_range_is_refused() {
    let cases: Vec<&str> = FIXTURE.lines().filter(|line| is_pinned(line)).collect();
    assert_eq!(cases.len(), 30);
    for line in cases {
        let outputs = line.rsplit('|').next().unwrap();
        assert!(outputs.trim_start().starts_with("i1000"), "{line}");
        let report = run_gate(line, |_| u64::MAX, run);
        assert_eq!(report.candidate_bits, [1, 1007], "{line}");
    }
}
