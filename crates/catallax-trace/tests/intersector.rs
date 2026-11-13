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
    if intersector.grazes(&overlap) {
        out.i(-2);
        return Ok(());
    }
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
        Ok(BoundaryPlaneTest::Grazing { .. }) => {
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
            let p = start(uv.as_ref().expect("a traced touch has UV curves"));
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
        for uv in &curve.uv_curves {
            if let Some(uv) = uv {
                raw(out, uv);
            } else {
                out.i(-1);
            }
        }
    }
    Ok(())
}

fn run(op: &str, input: &mut Input, out: &mut Output) {
    let a = surface(input);
    let da = domain(input);
    // `self`: the one surface twice, which is what makes an intersector test
    // for self intersection; the second domain follows the first.
    let second = (op != "self").then(|| surface(input));
    let db = domain(input);
    let b = second.as_ref().unwrap_or(&a);
    let has_tolerance = input.i() != 0;
    let tolerance = input.d();
    let mut intersector = Intersector::new(&a, da, b, db);
    if has_tolerance {
        intersector.approx_tolerance = tolerance;
    }
    if op == "intersect" || op == "self" {
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
            b,
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
            b,
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
            b,
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
                        // The probe records a grazing stage outcome as -1.
                        Ok(BoundaryPlaneTest::Grazing { .. }) => {
                            out.i(common::SUCCESS);
                            out.i(-1);
                        }
                        Err(SolverRefusal::Failed(fault)) => out.i(status_code(fault)),
                        Err(SolverRefusal::NotPorted(operation)) => {
                            panic!("not ported: {operation:?}")
                        }
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
    common::assert_same_with_workers(include_str!("fixtures/intersector/cases.txt"), |_| 16, run);
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

    /// The quadratic patch `z = a x^2 + b y^2` over `[-1.5, 1.5]^2`
    /// (parameters `[0, 1]^2`), exact as a degree-2 net.
    fn quadratic(a: f64, b: f64) -> BSplineSurface {
        let h = 1.5;
        let (x1, x2) = ([-h, 0.0, h], [h * h, -h * h, h * h]);
        let mut points = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                points.push(ControlPoint {
                    coordinates: [x1[i], x1[j], a * x2[i] + b * x2[j], NO_WEIGHT],
                });
            }
        }
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        BSplineSurface::new([2, 2], [knots.clone(), knots], [3, 3], points).expect("quadratic")
    }

    /// At a start point (no traced curve) with a previous direction, a
    /// singularity takes the branch nearest that direction among those
    /// leaving into both domains, even when the first such branch is not the
    /// first branch found. (SMLib keeps the first branch found in that case.)
    #[test]
    fn a_start_singularity_takes_an_inside_branch() {
        let saddle = quadratic(0.3, -0.3);
        let plane = quadratic(0.0, 0.0);
        let r = std::f64::consts::FRAC_1_SQRT_2;
        let middle = Vector2::new(0.5, 0.5);
        let halves = [
            Extent2::new(Vector2::new(0.5, 0.0), Vector2::new(1.0, 1.0)),
            Extent2::new(Vector2::new(0.0, 0.0), Vector2::new(0.5, 1.0)),
            Extent2::new(Vector2::new(0.0, 0.5), Vector2::new(1.0, 1.0)),
            Extent2::new(Vector2::new(0.0, 0.0), Vector2::new(1.0, 0.5)),
        ];
        let headings = [(r, r), (r, -r), (-r, r), (-r, -r)];
        for half in halves {
            for (x, y) in headings {
                let intersector = Intersector::new(&saddle, half, &plane, half);
                let previous = TsectPoint {
                    curve_tangent: catallax_trace::vector::Vector3::new(x, y, 0.0),
                    ..TsectPoint::default()
                };
                let mut point = TsectPoint::default();
                intersector
                    .compute_point_values(
                        [middle, middle],
                        &mut point,
                        Some(&previous),
                        None,
                        None,
                        1.0,
                    )
                    .expect("singular point");
                assert_eq!(point.kind, PointKind::Singularity);
                for s in 0..2 {
                    let probe = point.uv[s] + 1e-6 * point.uv_deriv[s];
                    assert!(
                        half.contains_point(probe, 1e-11),
                        "branch {:?} leaves {half:?} (heading {x}, {y})",
                        point.curve_tangent
                    );
                }
            }
        }
    }
}

#[test]
fn grazing_natural_boundaries_keep_absent_uv_curves_explicit() {
    let a = grazing_flat_patch(0.0, 1.0, false);
    let b = grazing_flat_patch(1.0, 2.0, false);
    let curves = Intersector::new(&a, a.natural_domain(), &b, b.natural_domain())
        .do_intersection()
        .unwrap();
    assert_eq!(curves.len(), 1);
    assert!(curves[0].uv_curves.iter().all(Option::is_none));
    let domain = curves[0].curve.natural_interval();
    let start = curves[0].curve.evaluate_point(domain.min).unwrap();
    let end = curves[0].curve.evaluate_point(domain.max).unwrap();
    assert_eq!(start.x, 1.0);
    assert_eq!(end.x, 1.0);
    assert_eq!(start.z, 0.0);
    assert_eq!(end.z, 0.0);
    assert_eq!((start.y - end.y).abs(), 1.0);
}

fn grazing_flat_patch(x0: f64, x1: f64, reverse_v: bool) -> BSplineSurface {
    use catallax_trace::spline::{ControlPoint, NO_WEIGHT};
    let (lo, hi) = if reverse_v { (1.0, 0.0) } else { (0.0, 1.0) };
    BSplineSurface::new(
        [1, 1],
        [vec![0.0, 0.0, 1.0, 1.0], vec![0.0, 0.0, 1.0, 1.0]],
        [2, 2],
        [[x0, lo, 0.0], [x0, hi, 0.0], [x1, lo, 0.0], [x1, hi, 0.0]]
            .map(|p| ControlPoint {
                coordinates: [p[0], p[1], p[2], NO_WEIGHT],
            })
            .to_vec(),
    )
    .unwrap()
}

#[test]
fn grazing_requested_domains_clip_the_shared_side() {
    let run = || {
        let mut bits = Vec::new();
        let near_knot = 2.0_f64.powi(-44);
        for (reverse_a, reverse_b, va, vb, expected) in [
            (false, false, [0.25, 0.75], [0.25, 0.75], [0.25, 0.75]),
            (false, false, [0.125, 0.875], [0.25, 0.5], [0.25, 0.5]),
            (false, true, [0.25, 0.75], [0.125, 0.5], [0.5, 0.75]),
            (true, true, [0.125, 0.875], [0.25, 0.5], [0.5, 0.75]),
            (
                false,
                false,
                [near_knot, 1.0 - near_knot],
                [near_knot, 1.0 - near_knot],
                [near_knot, 1.0 - near_knot],
            ),
            (false, false, [0.25, 0.5], [0.5, 0.75], [0.5, 0.5]),
        ] {
            let a = grazing_flat_patch(0.0, 1.0, reverse_a);
            let b = grazing_flat_patch(1.0, 2.0, reverse_b);
            let da = Extent2::new(Vector2::new(0.0, va[0]), Vector2::new(1.0, va[1]));
            let db = Extent2::new(Vector2::new(0.0, vb[0]), Vector2::new(1.0, vb[1]));
            for (first, d1, second, d2) in [(&a, da, &b, db), (&b, db, &a, da)] {
                let curves = Intersector::new(first, d1, second, d2)
                    .do_intersection()
                    .unwrap();
                assert_eq!(curves.len(), 1);
                for candidate in &curves {
                    bits.push(candidate.kind as u64);
                    bits.push(candidate.deviation.to_bits());
                    bits.extend(candidate.uv_curves.iter().map(|c| u64::from(c.is_some())));
                    for spline in std::iter::once(&candidate.curve)
                        .chain(candidate.uv_curves.iter().flatten())
                    {
                        bits.extend([
                            spline.degree() as u64,
                            spline.knots().len() as u64,
                            spline.control_points().len() as u64,
                        ]);
                        bits.extend(spline.knots().iter().map(|x| x.to_bits()));
                        bits.extend(
                            spline
                                .control_points()
                                .iter()
                                .flat_map(|p| p.coordinates.map(f64::to_bits)),
                        );
                    }
                }
                let curve = &curves[0].curve;
                let interval = curve.natural_interval();
                let start = curve.evaluate_point(interval.min).unwrap();
                let end = curve.evaluate_point(interval.max).unwrap();
                assert_eq!([start.y.min(end.y), start.y.max(end.y)], expected);
                for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
                    let point = curve.evaluate_point(interval.evaluate(t)).unwrap();
                    assert_eq!(point.x, 1.0);
                    assert_eq!(point.z, 0.0);
                    assert!(point.y >= expected[0] && point.y <= expected[1]);
                }
            }
        }
        bits
    };
    let expected = run();
    for workers in [2, 4] {
        std::thread::scope(|scope| {
            let jobs: Vec<_> = (0..workers).map(|_| scope.spawn(run)).collect();
            for job in jobs {
                assert_eq!(job.join().unwrap(), expected);
            }
        });
    }
}

#[test]
fn grazing_requested_domains_exclude_missing_sides_and_disjoint_ranges() {
    let a = grazing_flat_patch(0.0, 1.0, false);
    let b = grazing_flat_patch(1.0, 2.0, false);
    let unit = a.natural_domain();
    for (da, db) in [
        (
            Extent2::new(Vector2::new(0.0, 0.0), Vector2::new(0.75, 1.0)),
            unit,
        ),
        (
            unit,
            Extent2::new(Vector2::new(0.25, 0.0), Vector2::new(1.0, 1.0)),
        ),
        (
            Extent2::new(Vector2::new(0.0, 0.125), Vector2::new(1.0, 0.25)),
            Extent2::new(Vector2::new(0.0, 0.5), Vector2::new(1.0, 0.75)),
        ),
    ] {
        for (first, d1, second, d2) in [(&a, da, &b, db), (&b, db, &a, da)] {
            assert!(
                Intersector::new(first, d1, second, d2)
                    .do_intersection()
                    .unwrap()
                    .is_empty()
            );
        }
    }
}
