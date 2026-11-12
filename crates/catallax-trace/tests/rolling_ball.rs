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
        StepRefusal::NotPorted(_) => -2,
        StepRefusal::SurfaceExtension => -3,
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
        "rb-step" | "rb-step-failed" | "rb-step-extended" => {
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
const PINNED: [&str; 2] = ["rb-step-failed", "rb-step-extended"];

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

/// A step that a parameter domain clips, and that is not a boundary hit,
/// ends on the clipping boundary: SMLib re-solves it with
/// `LocalCurveIntersect` on an offset surface. The gate above compares these
/// steps bit for bit; this keeps enough of them in it.
#[test]
fn the_gate_keeps_its_clipped_steps() {
    let clipped = FIXTURE
        .lines()
        .filter(|line| is_case_of(line, "rb-step"))
        .filter(|line| {
            let outputs = line.rsplit('|').next().unwrap().trim_start();
            outputs.starts_with("i1000 i1 i1 i0")
        })
        .count();
    assert!(clipped >= CLIPPED_CASES, "clipped steps: {clipped}");
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

/// A clipped step whose clipping side lies past the base surface's domain
/// (SMLib's fillet lets the domains exceed the natural ones): SMLib extends
/// the base surface there (`sm_ExtendSurface`) and re-solves the step on the
/// extension. The port does not extend a surface: it refuses with
/// `StepRefusal::SurfaceExtension` rather than clamp. SMLib's answer is found
/// with the clipped flag set.
#[test]
fn a_clipped_step_past_the_domain_is_refused() {
    let cases = pinned_cases("rb-step-extended");
    assert_eq!(cases.len(), EXTENDED_CASES);
    for line in cases {
        let outputs = line.rsplit('|').next().unwrap();
        assert!(outputs.trim_start().starts_with("i1000 i1 i1 i0"), "{line}");
        let report = run_gate(line, |_| u64::MAX, run);
        assert_eq!(report.candidate_bits, [1, (-3i64) as u64], "{line}");
    }
}

const EXTENDED_CASES: usize = 4;

/// A plane and a quadratic sheet meeting on `B(s, 0)` at `A(s - s^2, s)`,
/// both of zero radius. The caller's domain for the plane stops at
/// `u = 0.249, v = 0.45`, but the iso-curve at the `u` boundary meets the
/// sheet at `v = 0.468...`; SMLib's natural interval admits that meeting.
/// The port restricts the iso-curve to the caller's domain and refuses the
/// point outside it instead of reporting it as the step's answer.
#[test]
fn a_clipped_boundary_meeting_stays_in_the_callers_domain() {
    fn surface(degrees: [usize; 2], knots: [Vec<f64>; 2], points: Vec<[f64; 4]>) -> BSplineSurface {
        let counts = [points.len() / (degrees[1] + 1), degrees[1] + 1];
        let points = points
            .into_iter()
            .map(|coordinates| ControlPoint { coordinates })
            .collect();
        BSplineSurface::new(degrees, knots, counts, points).expect("valid test surface")
    }
    let plane = surface(
        [1, 1],
        [vec![0.0, 0.0, 1.0, 1.0], vec![0.0, 0.0, 1.0, 1.0]],
        vec![
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0, 1.0],
        ],
    );
    let sheet = surface(
        [2, 1],
        [vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], vec![0.0, 0.0, 1.0, 1.0]],
        vec![
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
            [0.5, 0.5, 0.0, 1.0],
            [0.5, 0.5, 1.0, 1.0],
            [0.0, 1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0, 1.0],
        ],
    );
    let ball = RollingBall::new([&plane, &sheet], [0.0, 0.0], [true, true], 1e-9);
    let mut current = TsectPoint {
        uv: [Vector2::new(0.24, 0.4), Vector2::new(0.4, 0.0)],
        uv_deriv: [Vector2::new(0.2, 1.0), Vector2::new(1.0, 0.0)],
        curve_point: Vector3::new(0.24, 0.4, 0.0),
        curve_tangent: Vector3::new(0.2, 1.0, 0.0),
        ..TsectPoint::default()
    };
    current.curve_tangent.unitize().expect("nonzero tangent");
    current.kind = PointKind::Crossing;
    let plane_domain = Extent2::new(Vector2::new(0.0, 0.0), Vector2::new(0.249, 0.45));
    let sheet_domain = Extent2::new(Vector2::new(0.0, 0.0), Vector2::new(1.0, 1.0));
    let outcome = ball
        .step_solve(&current, 0.06, plane_domain, sheet_domain)
        .expect("a step");
    assert!(
        plane_domain.contains_point(outcome.next.uv[0], 0.0),
        "the step answered outside the caller's plane domain: {:?}",
        outcome.next.uv[0]
    );
    assert!(
        sheet_domain.contains_point(outcome.next.uv[1], 0.0),
        "the step answered outside the caller's sheet domain: {:?}",
        outcome.next.uv[1]
    );
}
