#![allow(missing_docs)]
//! Differential gate for the constant-radius rolling-ball equations
//! (`SmConstantRadiusFS` with `SmFilletSolver::ComputeSurfaceValues`,
//! `ComputePointValues` and `StepSolve`) against saved outputs of the pinned
//! C++ SMLib.

mod common;

use catallax_trace::extent::Extent2;
use catallax_trace::rolling_ball::{RollingBall, StepRefusal};
use catallax_trace::spline::{BSplineSurface, ControlPoint};
use catallax_trace::surface_trace::{PointKind, TsectPoint};
use catallax_trace::vector::{Vector2, Vector3};
use common::{Input, Output, run_gate, status_code};

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

fn domain(input: &mut Input) -> Extent2 {
    let (a, b, c, d) = (input.d(), input.d(), input.d(), input.d());
    Extent2::new(Vector2::new(a, b), Vector2::new(c, d))
}

fn kind_number(kind: PointKind) -> i64 {
    match kind {
        PointKind::Undefined => 0,
        PointKind::Unknown => 1,
        PointKind::Crossing => 2,
        PointKind::TangentPoint => 3,
        PointKind::TangentCurve => 4,
        PointKind::Singularity => 5,
        PointKind::Coincidence => 6,
        PointKind::Touching => 7,
    }
}

fn emit_point(out: &mut Output, point: &TsectPoint) {
    out.i(kind_number(point.kind));
    for surface in 0..2 {
        out.v3(point.position(surface));
        out.v3(point.du(surface));
        out.v3(point.dv(surface));
        out.v3(point.normal(surface));
        out.v2(point.uv[surface]);
        out.v2(point.uv_deriv[surface]);
    }
    out.v3(point.curve_point);
    out.v3(point.curve_tangent);
    out.d(point.tangent_plane_angle);
}

/// `SmTsectPnt()` as SMLib leaves it: the tangents are undefined until a
/// solve sets them, where the port's default point holds zeros.
fn blank_point() -> TsectPoint {
    TsectPoint {
        uv_deriv: [Vector2::default(); 2],
        curve_tangent: Vector3::default(),
        ..TsectPoint::default()
    }
}

fn uv_pair(input: &mut Input) -> [Vector2; 2] {
    [input.v2(), input.v2()]
}

fn refusal_status(refusal: StepRefusal) -> i64 {
    match refusal {
        StepRefusal::Failed(fault) => status_code(fault),
        StepRefusal::ClippedStep => -1,
        StepRefusal::NotPorted(_) => -2,
    }
}

fn run(op: &str, input: &mut Input, out: &mut Output) {
    let (r1, r2) = (input.d(), input.d());
    let (o1, o2) = (input.i() != 0, input.i() != 0);
    let tolerance = input.d();
    let (s1, s2) = (raw_surface(input), raw_surface(input));
    let ball = RollingBall::new([&s1, &s2], [r1, r2], [o1, o2], tolerance);
    match op {
        "rb-surface" => {
            let uv = uv_pair(input);
            let mut point = blank_point();
            if out
                .status(ball.compute_surface_values(uv, &mut point))
                .is_some()
            {
                emit_point(out, &point);
            }
        }
        "rb-point" => {
            let uv = uv_pair(input);
            let direction = input.d();
            let domains = [domain(input), domain(input)];
            let previous = (input.i() != 0).then(|| {
                let mut previous = blank_point();
                previous.curve_tangent = input.v3();
                previous
            });
            let mut point = blank_point();
            let result =
                ball.compute_point_values(uv, domains, direction, &mut point, previous.as_ref());
            match result {
                Ok(()) => {
                    out.i(1000);
                    emit_point(out, &point);
                }
                Err(refusal) => out.i(refusal_status(refusal)),
            }
        }
        "rb-step" | "rb-step-clipped" | "rb-step-failed" => {
            let mut current = blank_point();
            current.uv = uv_pair(input);
            current.uv_deriv = uv_pair(input);
            current.curve_point = input.v3();
            current.curve_tangent = input.v3();
            let step = input.d();
            let (d1, d2) = (domain(input), domain(input));
            match ball.step_solve(&current, step, d1, d2) {
                Ok(outcome) => {
                    out.i(1000);
                    out.flag(outcome.found);
                    out.flag(outcome.clipped);
                    out.flag(outcome.boundary_hit);
                    out.d(outcome.step);
                    if outcome.found {
                        out.v2(outcome.next.uv[0]);
                        out.v2(outcome.next.uv[1]);
                    }
                }
                Err(refusal) => out.i(refusal_status(refusal)),
            }
        }
        other => panic!("unknown operation {other}"),
    }
}

const FIXTURE: &str = include_str!("fixtures/rolling-ball/cases.txt");
const PINNED: [&str; 2] = ["rb-step-clipped", "rb-step-failed"];

fn is_case_of(line: &str, name: &str) -> bool {
    line.strip_prefix(name)
        .is_some_and(|rest| rest.starts_with(' '))
}

fn is_pinned(line: &str) -> bool {
    PINNED.iter().any(|name| is_case_of(line, name))
}

fn pinned_cases(name: &str) -> Vec<&'static str> {
    FIXTURE
        .lines()
        .filter(|line| is_case_of(line, name))
        .collect()
}

fn lines(pinned: bool) -> String {
    FIXTURE
        .lines()
        .filter(|line| line.starts_with('#') || is_pinned(line) == pinned)
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn rolling_ball_matches_smlib() {
    let report = run_gate(&lines(false), |_| 0, run);
    assert!(report.cases >= 150, "the gate lost cases: {}", report.cases);
    report.assert_passed(&["tangent", "touching", "singular", "boundary"]);
}

/// A step that a parameter domain clips, and that is not a boundary hit, goes
/// on in SMLib to `LocalCurveIntersect` on an offset surface. The port does
/// not carry that branch: it refuses with `StepRefusal::ClippedStep` where
/// SMLib answers. SMLib's answer is found with the clipped flag set.
#[test]
fn clipped_steps_are_refused() {
    let cases = pinned_cases("rb-step-clipped");
    assert_eq!(cases.len(), CLIPPED_CASES);
    for line in cases {
        let outputs = line.rsplit('|').next().unwrap();
        assert!(outputs.trim_start().starts_with("i1000 i1 i1 i0"), "{line}");
        // `run_gate` records tag and bits per token: tag 1 for an integer.
        let report = run_gate(line, |_| u64::MAX, run);
        assert_eq!(report.candidate_bits, [1, (-1i64) as u64], "{line}");
    }
}

const CLIPPED_CASES: usize = 4;

/// An offset evaluator that fails at the first guess leaves SMLib's solver
/// "found" and unable to converge, so SMLib's `StepSolve` answers found, with
/// the guessed parameters. The port refuses it, as the surface-surface wrappers
/// do, rather than report a failed guess as a good point.
#[test]
fn an_evaluator_failure_at_the_guess_is_refused() {
    let cases = pinned_cases("rb-step-failed");
    assert_eq!(cases.len(), FAILED_CASES);
    for line in cases {
        let outputs = line.rsplit('|').next().unwrap();
        assert!(outputs.trim_start().starts_with("i1000 i1 i0 i0"), "{line}");
        let report = run_gate(line, |_| u64::MAX, run);
        assert_eq!(report.candidate_bits, [1, 1001u64], "{line}");
    }
}

const FAILED_CASES: usize = 8;
