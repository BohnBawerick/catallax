#![allow(missing_docs)]
use super::*;
#[path = "../common/mod.rs"]
mod common;
use common::{Input, Output};

fn evaluate(_: &str, input: &mut Input, out: &mut Output) {
    assert_eq!(input.word(), "s");
    let surface = input.surface();
    let domain = Extent2::new(input.v2(), input.v2());
    let curve = input.curve();
    let interval = input.interval();
    let uv = input.v2();
    let mut point = TracePoint {
        curve_parameter: input.d(),
        ..TracePoint::default()
    };
    let tracer = DropTracer {
        surface: &surface,
        domain,
        curve: &curve,
        interval,
        tolerance: 0.001,
        working_interval: interval,
        working_domain: domain,
        max_distance: 0.0,
        breaks: Vec::new(),
        start_points: Vec::new(),
        curves: Vec::new(),
        drops: Vec::new(),
        deviations: Vec::new(),
    };
    match tracer.compute_point_values(uv, &mut point) {
        Ok(()) => {
            out.i(common::SUCCESS);
            for p in [
                point.surface_point,
                point.du,
                point.dv,
                point.normal,
                point.curve_tangent,
            ] {
                out.d(p.x);
                out.d(p.y);
                out.d(p.z);
            }
            for p in [point.uv, point.uv_tangent] {
                out.d(p.x);
                out.d(p.y);
            }
            out.d(point.drop_distance);
            out.d(point.normal_line_distance);
        }
        Err(SolverRefusal::Failed(fault)) => out.i(common::status_code(fault)),
        Err(_) => out.i(-1),
    }
}

#[test]
fn singular_drop_starts_match_smlib() {
    let fixture = include_str!("../fixtures/drop-start/cases.txt");
    let report = common::run_gate(fixture, |_| 16, evaluate);
    assert_eq!(report.cases, 41);
    report.assert_passed(&["singular", "near-singular", "regular"]);
    common::assert_same_with_workers(fixture, |_| 16, evaluate);
}
