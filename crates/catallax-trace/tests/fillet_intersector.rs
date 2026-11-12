#![allow(missing_docs)]
//! Differential gate for the stepping half of `SmFilletIntersector`
//! (`ComputePointValues`, `ComputePointValuesDelta`, `TestSpanAccuracy`,
//! `ComputeNextPoint`, `EvaluateLawPoint` and the constructor's tolerances)
//! against saved outputs of the pinned C++ SMLib.

mod common;

use catallax_trace::extent::Extent2;
use catallax_trace::fillet_intersector::FilletIntersector;
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

/// `SmTsectPnt()` as SMLib leaves it: every vector undefined, where the
/// port's default point holds zeros.
fn blank_point() -> TsectPoint {
    let undefined = Vector3::default();
    TsectPoint {
        uv: [Vector2::default(); 2],
        uv_deriv: [Vector2::default(); 2],
        surface: [[undefined; 9]; 2],
        curve_point: undefined,
        curve_tangent: undefined,
        ..TsectPoint::default()
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
    out.d(point.curve_parameter);
}

fn uv_pair(input: &mut Input) -> [Vector2; 2] {
    [input.v2(), input.v2()]
}

/// SMLib's status for a refusal; the port's own refusals are negative.
fn refusal_status(refusal: StepRefusal) -> i64 {
    match refusal {
        StepRefusal::Failed(fault) => status_code(fault),
        StepRefusal::ClippedStep => -1,
        StepRefusal::NotPorted(_) => -2,
        StepRefusal::SurfaceExtension => -3,
    }
}

/// A setup call whose failure is the case's answer: only a failure is printed.
fn setup<T>(out: &mut Output, result: Result<T, StepRefusal>) -> Option<T> {
    match result {
        Ok(value) => Some(value),
        Err(refusal) => {
            out.i(refusal_status(refusal));
            None
        }
    }
}

fn emit_status<T>(out: &mut Output, result: Result<T, StepRefusal>) -> Option<T> {
    match result {
        Ok(value) => {
            out.i(1000);
            Some(value)
        }
        Err(refusal) => {
            out.i(refusal_status(refusal));
            None
        }
    }
}

fn run(op: &str, input: &mut Input, out: &mut Output) {
    let (r1, r2) = (input.d(), input.d());
    let (o1, o2) = (input.i() != 0, input.i() != 0);
    let tolerance = input.d();
    let (s1, s2) = (raw_surface(input), raw_surface(input));
    let ball = RollingBall::new([&s1, &s2], [r1, r2], [o1, o2], tolerance);
    let (d1, d2) = (domain(input), domain(input));
    let mut fi = FilletIntersector::new(&ball, d1, d2);
    if op == "fi-ctor" {
        out.d(fi.approx_tolerance);
        out.d(fi.angle_tolerance);
        out.d(fi.trace_direction);
        out.flag(fi.clipped_by_through_point);
        return;
    }
    fi.approx_tolerance = input.d();
    fi.angle_tolerance = input.d();
    fi.trace_direction = input.d();
    fi.clipped_by_through_point = input.i() != 0;
    match op {
        "fi-values" => {
            let uv = uv_pair(input);
            let previous = (input.i() != 0).then(|| {
                let mut previous = blank_point();
                previous.curve_tangent = input.v3();
                previous
            });
            let mut point = blank_point();
            if emit_status(
                out,
                fi.compute_point_values(uv, &mut point, previous.as_ref()),
            )
            .is_some()
            {
                emit_point(out, &point);
            }
        }
        "fi-delta" => {
            let (uv, delta_uv) = (uv_pair(input), uv_pair(input));
            let (param, delta_param) = (input.d(), input.d());
            let (mut point, mut delta) = (blank_point(), blank_point());
            point.curve_parameter = param;
            delta.curve_parameter = delta_param;
            let result = fi
                .compute_point_values_delta(uv, delta_uv, &mut point, &mut delta)
                .map_err(StepRefusal::from);
            if emit_status(out, result).is_some() {
                emit_point(out, &point);
                emit_point(out, &delta);
            }
        }
        "fi-law" => {
            let guesses = uv_pair(input);
            let (origin, normal) = (input.v3(), input.v3());
            if let Some(found) = emit_status(out, fi.evaluate_law_point(origin, normal, guesses)) {
                // Where nothing is found SMLib returns its unset point: undefined
                // parameters and a zero curve parameter.
                let ([a, b], param) = found.unwrap_or(([Vector2::default(); 2], 0.0));
                out.v2(a);
                out.v2(b);
                out.d(param);
            }
        }
        "fi-span" | "fi-next" | "fi-next-clipped" => {
            let uv = uv_pair(input);
            let mut current = blank_point();
            if setup(out, fi.compute_point_values(uv, &mut current, None)).is_none() {
                return;
            }
            current.curve_parameter = input.d();
            if op == "fi-span" {
                let next_uv = uv_pair(input);
                let mut next = blank_point();
                if setup(
                    out,
                    fi.compute_point_values(next_uv, &mut next, Some(&current)),
                )
                .is_none()
                {
                    return;
                }
                next.curve_parameter = input.d();
                let middle = if input.i() != 0 {
                    let middle_uv = uv_pair(input);
                    let mut middle = blank_point();
                    if setup(
                        out,
                        fi.compute_point_values(middle_uv, &mut middle, Some(&current)),
                    )
                    .is_none()
                    {
                        return;
                    }
                    Some(middle)
                } else {
                    None
                };
                if let Some(span) = emit_status(
                    out,
                    fi.test_span_accuracy(&current, &mut next, middle.as_ref()),
                ) {
                    out.flag(span.good);
                    out.d(span.angle);
                    out.i(kind_number(next.kind));
                }
            } else {
                let step = input.d();
                let mut next = blank_point();
                if let Some(outcome) =
                    emit_status(out, fi.compute_next_point(&mut current, step, &mut next))
                {
                    out.flag(outcome.good);
                    out.flag(outcome.clipped);
                    out.flag(outcome.boundary_hit);
                    out.d(outcome.deviation);
                    out.d(outcome.angle);
                    out.i(kind_number(next.kind));
                    out.v2(next.uv[0]);
                    out.v2(next.uv[1]);
                    out.v2(next.uv_deriv[0]);
                    out.v2(next.uv_deriv[1]);
                    out.v3(next.curve_point);
                    out.v3(next.curve_tangent);
                    out.d(next.curve_parameter);
                    out.v2(current.uv_deriv[0]);
                    out.v2(current.uv_deriv[1]);
                    out.v3(current.curve_tangent);
                }
            }
        }
        other => panic!("unknown operation {other}"),
    }
}

const FIXTURE: &str = include_str!("fixtures/fillet-intersector/cases.txt");

fn is_case_of(line: &str, name: &str) -> bool {
    line.strip_prefix(name)
        .is_some_and(|rest| rest.starts_with(' '))
}

fn pinned_cases() -> Vec<&'static str> {
    FIXTURE
        .lines()
        .filter(|line| is_case_of(line, "fi-next-clipped"))
        .collect()
}

fn unpinned() -> String {
    FIXTURE
        .lines()
        .enumerate()
        .filter(|(index, line)| {
            !is_case_of(line, "fi-next-clipped") && !EXTENDED.contains(&(index + 1))
        })
        .map(|(_, line)| line)
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn fillet_intersector_matches_smlib() {
    let report = run_gate(&unpinned(), |_| 0, run);
    assert!(report.cases >= 500, "the gate lost cases: {}", report.cases);
    report.assert_passed(&["regular", "tangent", "touching", "singular", "boundary"]);
}

/// A step that a domain clips, and that is not a boundary hit, goes on in
/// SMLib to `LocalCurveIntersect` on an offset surface. The port refuses it
/// with `StepRefusal::ClippedStep`. SMLib's answer carries the clipped flag.
#[test]
fn clipped_steps_are_refused() {
    let cases = pinned_cases();
    assert_eq!(cases.len(), CLIPPED_CASES);
    for line in cases {
        let outputs: Vec<&str> = line
            .rsplit('|')
            .next()
            .unwrap()
            .split_whitespace()
            .collect();
        assert_eq!(
            (outputs[0], outputs[2], outputs[3]),
            ("i1000", "i1", "i0"),
            "{line}"
        );
        let report = run_gate(line, |_| u64::MAX, run);
        assert_eq!(report.candidate_bits, [1, (-1i64) as u64], "{line}");
    }
}

const CLIPPED_CASES: usize = 31;

/// Fixture lines (1-based) of five `fi-next` cases in which a Hermite
/// parameter curve leaves its surface's domain: three by rounding (a
/// parameter of `1.0000000000000002`, say), two by a little more on the
/// cubic surfaces. SMLib extends the base surface there with
/// `sm_ExtendSurface`, which is not ported, and its rail points differ from
/// the clamped ones, in the last bits for the first three (the deviation it
/// reports is `6.7e-16` where the clamped value would be `8.9e-16`). The
/// port refuses with `StepRefusal::SurfaceExtension`; SMLib answers
/// successfully, with a good point in three cases and a rejected span in two.
const EXTENDED: [usize; 5] = [550, 562, 583, 601, 615];

#[test]
fn evaluations_outside_the_domain_are_refused() {
    let lines: Vec<&str> = FIXTURE.lines().collect();
    for number in EXTENDED {
        let line = lines[number - 1];
        assert!(is_case_of(line, "fi-next"), "{line}");
        let outputs: Vec<&str> = line
            .rsplit('|')
            .next()
            .unwrap()
            .split_whitespace()
            .collect();
        assert_eq!(outputs[0], "i1000", "{line}");
        let report = run_gate(line, |_| u64::MAX, run);
        assert_eq!(report.candidate_bits, [1, (-3i64) as u64], "{line}");
    }
    // No other case refuses for this reason.
    let others: Vec<usize> = lines
        .iter()
        .enumerate()
        .filter(|(index, line)| {
            !line.starts_with('#')
                && !EXTENDED.contains(&(index + 1))
                && run_gate(line, |_| u64::MAX, run).candidate_bits == [1, (-3i64) as u64]
        })
        .map(|(index, _)| index + 1)
        .collect();
    assert!(others.is_empty(), "unpinned refusals at lines {others:?}");
}
