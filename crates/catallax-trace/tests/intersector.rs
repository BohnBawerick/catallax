#![allow(missing_docs)]
//! Differential gate for SMLib's surface-surface intersector
//! (`SmAdvSurfaceIntersector`) stage by stage, against saved outputs of the
//! pinned C++ SMLib read through a probe subclass.

mod common;

use catallax_trace::extent::Extent2;
use catallax_trace::global_solver::SolverRefusal;
use catallax_trace::spline::{BSplineCurve, BSplineSurface};
use catallax_trace::surface_intersector::{BoundaryPlaneTest, Intersector};
use catallax_trace::surface_surface::local_surface_intersect;
use catallax_trace::surface_trace::{Tracer, TsectPoint};
use catallax_trace::vector::Vector2;
use common::{Input, Output, run_gate, status_code};

fn domain(input: &mut Input) -> Extent2 {
    let (a, b, c, d) = (input.d(), input.d(), input.d(), input.d());
    Extent2::new(Vector2::new(a, b), Vector2::new(c, d))
}

fn surface(input: &mut Input) -> BSplineSurface {
    if input.word() == "s" {
        return input.surface();
    }
    let kind = input.i();
    let (origin, x, y) = (input.v3(), input.v3(), input.v3());
    let scale = input.v2();
    let bounds = domain(input);
    if kind == 1 {
        BSplineSurface::canonical_plane(origin, x, y).expect("plane")
    } else {
        BSplineSurface::plane(origin, x, y, scale, bounds).expect("plane")
    }
}

fn extent2(out: &mut Output, e: Extent2) {
    out.v2(e.min);
    out.v2(e.max);
}

fn refused(out: &mut Output, refusal: SolverRefusal) {
    match refusal {
        SolverRefusal::Failed(fault) => out.i(status_code(fault)),
        SolverRefusal::NotPorted(what) => panic!("not ported: {what:?}"),
    }
}

/// The fields of `SmTsectPnt` the probe prints.
fn point(out: &mut Output, p: &TsectPoint) {
    out.i(p.kind as i64);
    for s in 0..2 {
        out.i(p.param_types[s][0] as i64);
        out.i(p.param_types[s][1] as i64);
    }
    for s in 0..2 {
        out.v2(p.uv[s]);
        out.v2(p.uv_deriv[s]);
    }
    for s in 0..2 {
        // SrfPos, SrfDu, SrfDv, SrfDuu, SrfDuv, SrfDvv, SrfNorm.
        for cell in [0, 3, 1, 6, 4, 2, 8] {
            out.v3(p.surface[s][cell]);
        }
    }
    out.v3(p.curve_point);
    out.v3(p.curve_tangent);
    out.d(p.tangent_plane_angle);
    out.d(p.deviation);
}

/// `DoIntersection` up to `FindBoundaryStartPoints`, in the probe's order.
fn boundary(intersector: &Intersector<'_>, out: &mut Output) -> Result<(), SolverRefusal> {
    let trees = intersector.trees()?;
    let overlap = match intersector.overlap(&trees) {
        Ok(Some(overlap)) => overlap,
        Ok(None) => {
            out.i(-1);
            return Ok(());
        }
        // A grazing leaf pair: SMLib's `GrazingBoxProcess`.
        Err(SolverRefusal::NotPorted(_)) => {
            out.i(-2);
            return Ok(());
        }
        Err(refusal) => return Err(refusal),
    };
    out.i(0);
    match intersector.boundary_plane_test(overlap.bbox) {
        Ok(BoundaryPlaneTest::Unknown) => {
            out.i(common::SUCCESS);
            out.i(0);
        }
        Ok(BoundaryPlaneTest::Disjoint) => {
            out.i(common::SUCCESS);
            out.i(1);
            return Ok(());
        }
        Err(SolverRefusal::NotPorted(_)) => {
            out.i(common::SUCCESS);
            out.i(-3);
            return Ok(());
        }
        Err(refusal) => return Err(refusal),
    }
    let mut tracer = Tracer {
        trace_direction: 1.0,
        ..Tracer::default()
    };
    let mut touches = Vec::new();
    tracer.find_boundary_start_points(intersector, overlap.bbox, &mut touches)?;
    out.i(common::SUCCESS);
    out.i(tracer.start_points.len() as i64);
    for p in &tracer.start_points {
        point(out, p);
    }
    out.i(tracer.through_points.len() as i64);
    out.i(touches.len() as i64);
    for touch in &touches {
        out.i(touch.kind as i64);
        out.d(touch.deviation);
        let start = |curve: &BSplineCurve| {
            curve
                .evaluate_point(curve.natural_interval().min)
                .expect("touching curve")
        };
        out.v3(start(&touch.curve));
        for uv in &touch.uv_curves {
            let p = start(uv);
            out.d(p.x);
            out.d(p.y);
        }
    }
    Ok(())
}

/// The stored form of a curve, as the probe prints it.
fn raw(out: &mut Output, curve: &BSplineCurve) {
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

/// The whole `DoIntersection`, in the probe's order.
fn intersect(intersector: &Intersector<'_>, out: &mut Output) -> Result<(), SolverRefusal> {
    let curves = intersector.do_intersection()?;
    out.i(common::SUCCESS);
    out.i(curves.len() as i64);
    for curve in &curves {
        out.i(curve.kind as i64);
        out.d(curve.deviation);
        raw(out, &curve.curve);
        raw(out, &curve.uv_curves[0]);
        raw(out, &curve.uv_curves[1]);
    }
    Ok(())
}

fn run(op: &str, input: &mut Input, out: &mut Output) {
    let a = surface(input);
    let da = domain(input);
    let b = surface(input);
    let db = domain(input);
    let has_tolerance = input.i() != 0;
    let tolerance = input.d();
    let mut intersector = Intersector::new(&a, da, &b, db);
    if has_tolerance {
        intersector.approx_tolerance = tolerance;
    }
    if op == "intersect" {
        if let Err(refusal) = intersect(&intersector, out) {
            refused(out, refusal);
        }
        return;
    }
    if op == "irefuse" {
        match intersect(&intersector, &mut Output::default()) {
            Err(SolverRefusal::NotPorted(_)) => out.i(-1),
            other => panic!("expected a refusal, got {other:?}"),
        }
        return;
    }
    if op == "boundary" {
        if let Err(refusal) = boundary(&intersector, out) {
            refused(out, refusal);
        }
        return;
    }
    if op == "brefuse" {
        match boundary(&intersector, &mut Output::default()) {
            Err(SolverRefusal::NotPorted(_)) => out.i(-1),
            other => panic!("expected a refusal, got {other:?}"),
        }
        return;
    }
    if op == "classify" {
        let first = [input.v2(), input.v2()];
        let met = match local_surface_intersect(
            &a,
            intersector.domains[0],
            &b,
            intersector.domains[1],
            intersector.approx_tolerance,
            first[0],
            first[1],
            None,
            None,
            true,
        ) {
            Ok(met) => met,
            Err(refusal) => return refused(out, refusal),
        };
        out.i(common::SUCCESS);
        out.flag(met.is_some());
        let Some(pair) = met else {
            return;
        };
        let mut start = TsectPoint::default();
        if let Err(refusal) = intersector.compute_point_values(
            [pair.uv, pair.other_uv],
            &mut start,
            None,
            None,
            None,
            1.0,
        ) {
            return refused(out, refusal);
        }
        out.i(common::SUCCESS);
        out.i(start.kind as i64);
        out.v2(pair.uv);
        out.v2(pair.other_uv);
        match intersector.classify_boundary_start(&start) {
            Ok(degenerate) => {
                out.i(common::SUCCESS);
                out.flag(degenerate);
            }
            Err(refusal) => refused(out, refusal),
        }
        return;
    }
    if op == "trace" {
        let first = [input.v2(), input.v2()];
        let met = match local_surface_intersect(
            &a,
            intersector.domains[0],
            &b,
            intersector.domains[1],
            intersector.approx_tolerance,
            first[0],
            first[1],
            None,
            None,
            true,
        ) {
            Ok(met) => met,
            Err(refusal) => return refused(out, refusal),
        };
        out.i(common::SUCCESS);
        out.flag(met.is_some());
        let Some(pair) = met else {
            return;
        };
        let mut start = TsectPoint::default();
        if let Err(refusal) = intersector.compute_point_values(
            [pair.uv, pair.other_uv],
            &mut start,
            None,
            None,
            None,
            1.0,
        ) {
            return refused(out, refusal);
        }
        out.i(common::SUCCESS);
        let mut tracer = Tracer {
            trace_direction: 1.0,
            ..Tracer::default()
        };
        if let Err(refusal) = tracer.trace_curve(&intersector, &start, &[]) {
            return refused(out, refusal);
        }
        out.i(common::SUCCESS);
        out.i(tracer.curve_points.len() as i64);
        out.flag(tracer.curve_is_closed);
        for p in &tracer.curve_points {
            out.i(p.kind as i64);
            out.v2(p.uv[0]);
            out.v2(p.uv[1]);
            out.v3(p.curve_point);
            out.v3(p.curve_tangent);
            out.d(p.curve_parameter);
            out.d(p.deviation);
            out.d(p.tangent_plane_angle);
        }
        return;
    }
    if op == "step" {
        let first = [input.v2(), input.v2()];
        let size = input.d();
        let met = match local_surface_intersect(
            &a,
            intersector.domains[0],
            &b,
            intersector.domains[1],
            intersector.approx_tolerance,
            first[0],
            first[1],
            None,
            None,
            true,
        ) {
            Ok(met) => met,
            Err(refusal) => return refused(out, refusal),
        };
        out.i(common::SUCCESS);
        out.flag(met.is_some());
        let Some(pair) = met else {
            return;
        };
        let mut start = TsectPoint::default();
        if let Err(refusal) = intersector.compute_point_values(
            [pair.uv, pair.other_uv],
            &mut start,
            None,
            None,
            None,
            1.0,
        ) {
            return refused(out, refusal);
        }
        out.i(common::SUCCESS);
        let step = intersector.compute_step_size(size);
        out.d(step);
        let mut next = TsectPoint::default();
        match intersector.compute_next_point(&start, step, &mut next, None, 1.0) {
            Ok(outcome) => {
                out.i(common::SUCCESS);
                out.flag(outcome.good);
                out.flag(outcome.clipped);
                out.flag(outcome.boundary_hit);
                out.d(outcome.deviation);
                out.d(outcome.angle);
                if outcome.good {
                    point(out, &next);
                    out.d(next.curve_parameter);
                }
            }
            Err(refusal) => refused(out, refusal),
        }
        return;
    }
    if op == "point" {
        let two = input.i() != 0;
        let first = [input.v2(), input.v2()];
        let second = [input.v2(), input.v2()];
        let mut p = TsectPoint::default();
        match intersector.compute_point_values(first, &mut p, None, None, None, 1.0) {
            Ok(()) => {
                out.i(common::SUCCESS);
                point(out, &p);
            }
            Err(refusal) => return refused(out, refusal),
        }
        if two {
            let mut q = TsectPoint::default();
            match intersector.compute_point_values(second, &mut q, Some(&p), None, None, 1.0) {
                Ok(()) => {
                    out.i(common::SUCCESS);
                    point(out, &q);
                }
                Err(refusal) => refused(out, refusal),
            }
        }
        return;
    }
    assert!(op == "search" || op == "planes", "unknown operation {op}");
    out.i(common::SUCCESS);
    out.d(intersector.approx_tolerance);
    for i in 0..2 {
        out.flag(intersector.closed_u[i]);
        out.flag(intersector.closed_v[i]);
        extent2(out, intersector.domains[i]);
    }
    let result = intersector.trees().and_then(|trees| {
        Ok((
            trees.root_boxes_meet(intersector.approx_tolerance)?,
            intersector.search(&trees)?,
        ))
    });
    match result {
        Ok((false, _)) => out.i(-1),
        Ok((true, found)) => {
            out.i(common::SUCCESS);
            out.i(i64::from(found.is_some()));
            if let Some(found) = found {
                let p = found.bbox;
                for value in [
                    p.min.x,
                    p.min.y,
                    p.min.z,
                    p.max.x,
                    p.max.y,
                    p.max.z,
                    found.count,
                ] {
                    out.d(value);
                }
                extent2(out, found.domains[0]);
                extent2(out, found.domains[1]);
                if op == "planes" {
                    match intersector.boundary_plane_test(found.bbox) {
                        Ok(BoundaryPlaneTest::Unknown) => {
                            out.i(common::SUCCESS);
                            out.i(0);
                        }
                        Ok(BoundaryPlaneTest::Disjoint) => {
                            out.i(common::SUCCESS);
                            out.i(1);
                        }
                        // SMLib grazes here; the port refuses.
                        Err(SolverRefusal::NotPorted(_)) => {
                            out.i(common::SUCCESS);
                            out.i(-1);
                        }
                        Err(SolverRefusal::Failed(fault)) => out.i(status_code(fault)),
                    }
                }
            }
        }
        Err(SolverRefusal::Failed(fault)) => out.i(status_code(fault)),
        Err(SolverRefusal::NotPorted(what)) => panic!("not ported: {what:?}"),
    }
}

#[test]
fn intersector_search_matches_smlib() {
    let report = run_gate(include_str!("fixtures/intersector/cases.txt"), |_| 16, run);
    report.assert_passed(&[
        "tangent",
        "touching",
        "singular",
        "near-singular",
        "boundary",
    ]);
}

/// Regressions from the step 8 part 2 review.
mod review {
    use catallax_trace::extent::Extent2;
    use catallax_trace::spline::{BSplineSurface, ControlPoint, NO_WEIGHT};
    use catallax_trace::surface_intersector::Intersector;
    use catallax_trace::surface_trace::{PointKind, TsectPoint};
    use catallax_trace::vector::Vector2;

    /// A degree-(1, 1) patch over `[-1, 1]^2` with v-fast corners.
    fn bilinear(corners: [(f64, f64, f64); 4]) -> BSplineSurface {
        let knots = vec![-1.0, -1.0, 1.0, 1.0];
        let points = corners
            .iter()
            .map(|&(x, y, z)| ControlPoint {
                coordinates: [x, y, z, NO_WEIGHT],
            })
            .collect();
        BSplineSurface::new([1, 1], [knots.clone(), knots], [2, 2], points).expect("bilinear")
    }

    /// The saddle `z = u v` has no `z_uu` or `z_vv`, only a mixed term, so
    /// SMLib's planarity test (which reads `N . S_uu` and `N . S_vv` only)
    /// takes it for planar and calls its meeting with the plane `z = 0`
    /// (the crossing lines `u v = 0`) a tangent point. With the mixed term
    /// the saddle is curved in either order, and its meeting with the plane
    /// is classified from the saddle's curvatures, never as a tangent point.
    #[test]
    fn a_mixed_curvature_saddle_is_not_planar() {
        let saddle = bilinear([
            (-1.0, -1.0, 1.0),
            (-1.0, 1.0, -1.0),
            (1.0, -1.0, -1.0),
            (1.0, 1.0, 1.0),
        ]);
        let plane = bilinear([
            (-1.0, -1.0, 0.0),
            (-1.0, 1.0, 0.0),
            (1.0, -1.0, 0.0),
            (1.0, 1.0, 0.0),
        ]);
        let square = Extent2::new(Vector2::new(-1.0, -1.0), Vector2::new(1.0, 1.0));
        let origin = Vector2::new(0.0, 0.0);
        for (a, b) in [(&saddle, &plane), (&plane, &saddle)] {
            let intersector = Intersector::new(a, square, b, square);
            let mut point = TsectPoint::default();
            match intersector.compute_point_values(
                [origin, origin],
                &mut point,
                None,
                None,
                None,
                1.0,
            ) {
                // The saddle against a plane: four branches meet there.
                Ok(()) => assert_eq!(point.kind, PointKind::Singularity),
                other => panic!("expected the saddle's singularity, got {other:?}"),
            }
        }
    }
}
