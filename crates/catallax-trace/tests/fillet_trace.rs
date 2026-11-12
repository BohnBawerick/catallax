#![allow(missing_docs)]
//! Differential gate for the fillet intersector's tracing loop,
//! `SmSurfaceIntersector::DoPointIntersection` as the fillet intersector runs
//! it, up to the call of `FlushCurve`, against saved outputs of the pinned C++
//! SMLib. The oracle replaces `FlushCurve` with a dump of the traced points.

mod common;

use catallax_trace::extent::Extent2;
use catallax_trace::fillet_intersector::{FilletIntersector, TraceRequest};
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

/// `SmTsectPnt()` as SMLib leaves it: every vector undefined.
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

fn run(_op: &str, input: &mut Input, out: &mut Output) {
    let (r1, r2) = (input.d(), input.d());
    let (o1, o2) = (input.i() != 0, input.i() != 0);
    let tolerance = input.d();
    let (s1, s2) = (raw_surface(input), raw_surface(input));
    let ball = RollingBall::new([&s1, &s2], [r1, r2], [o1, o2], tolerance);
    let (d1, d2) = (domain(input), domain(input));
    let mut fi = FilletIntersector::new(&ball, d1, d2);
    let extend_before_start = input.i() != 0;
    let extend_after_end = input.i() != 0;
    let start_direction = (input.i() != 0).then(|| input.v3());
    let end_direction = (input.i() != 0).then(|| input.v3());
    let approx_tolerance = (input.i() != 0).then(|| input.d());
    let angle_tolerance = (input.i() != 0).then(|| input.d());
    fi.extension_distance = input.d();
    let count = input.i() as usize;
    let mut through: Vec<TsectPoint> = (0..count)
        .map(|_| {
            let mut point = blank_point();
            point.uv = [input.v2(), input.v2()];
            point
        })
        .collect();
    let request = TraceRequest {
        extend_before_start,
        extend_after_end,
        start_direction,
        end_direction,
        approx_tolerance,
        angle_tolerance,
    };
    let result = fi.do_point_intersection(&mut through, request);
    let (status, trace) = match result {
        Ok(trace) => (1000, Some(trace)),
        Err(StepRefusal::Failed(fault)) => (status_code(fault), None),
        Err(StepRefusal::ClippedStep) => (-1, None),
        Err(StepRefusal::NotPorted(_)) => (-2, None),
    };
    out.i(status);
    if status < 0 {
        return;
    }
    // The oracle's flag for a clipped step, which the port refuses: none here.
    out.i(0);
    let traced = trace.as_ref().is_some_and(|trace| trace.traced);
    out.flag(traced);
    if let Some(trace) = trace.filter(|trace| trace.traced) {
        out.i(trace.points.len() as i64);
        for point in &trace.points {
            out.i(kind_number(point.kind));
            for s in 0..2 {
                out.v2(point.uv[s]);
                out.v2(point.uv_deriv[s]);
            }
            out.v3(point.curve_point);
            out.v3(point.curve_tangent);
            out.d(point.curve_parameter);
            out.d(point.deviation);
            out.d(point.tangent_plane_angle);
        }
        out.flag(trace.curve_is_closed);
    }
    for point in &through {
        out.i(kind_number(point.kind));
        out.d(point.curve_parameter);
    }
}

const FIXTURE: &str = include_str!("fixtures/fillet-trace/cases.txt");

fn is_case_of(line: &str, name: &str) -> bool {
    line.strip_prefix(name)
        .is_some_and(|rest| rest.starts_with(' '))
}

fn unpinned() -> String {
    FIXTURE
        .lines()
        .filter(|line| !is_case_of(line, "fi-trace-clipped"))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn fillet_trace_matches_smlib() {
    let report = run_gate(&unpinned(), |_| 0, run);
    assert!(report.cases >= 150, "the gate lost cases: {}", report.cases);
    report.assert_passed(&["regular", "tangent", "touching", "singular"]);
}

/// A trace that steps into a clipped step that is not a boundary hit goes on,
/// in SMLib, to `LocalCurveIntersect` on an offset surface, and SMLib's trace
/// succeeds. The port refuses the whole trace with `StepRefusal::ClippedStep`.
/// Every such trace runs on to a side of a domain or past its through points.
#[test]
fn traces_with_clipped_steps_are_refused() {
    let cases: Vec<&str> = FIXTURE
        .lines()
        .filter(|line| is_case_of(line, "fi-trace-clipped"))
        .collect();
    assert_eq!(cases.len(), CLIPPED_CASES);
    for line in cases {
        let outputs: Vec<&str> = line
            .rsplit('|')
            .next()
            .unwrap()
            .split_whitespace()
            .collect();
        assert_eq!((outputs[0], outputs[1]), ("i1000", "i1"), "{line}");
        let report = run_gate(line, |_| u64::MAX, run);
        assert_eq!(report.candidate_bits, [1, (-1i64) as u64], "{line}");
    }
}

const CLIPPED_CASES: usize = 131;
