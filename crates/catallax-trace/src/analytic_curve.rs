// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: output arrays become owned values; the shared
// closest-line helper returns None for singular or parallel lines, and
// objects keep the original spline identities. Tolerance tests produce
// candidates only.
// Translated from SmCurve.cpp (AnalyticCurveIntersect), SmEllipse.cpp
// (IntersectWithLine, IntersectWithEllipse), SmLine.cpp
// (SmLine::IntersectWithLine), SmGeomUtility.cpp
// (smgu_SegmentSegmentIntersect, sm_AddPoint) at
// c9c5979f04cea8c38213004d4f04df0dcfb40d05; rust-brep has no port of these.
// Endpoint tests reuse the ported segment-point distance, as the donor
// smgu_is_point_on_segment does in sm_geom_utility.rs.
//! The line and circle dispatch of `SmCurve::AnalyticCurveIntersect`.

use crate::Fault;
use crate::circle::{AnalyticArc, Circle, analytic_arc};
use crate::curve_cache::CurveCache;
use crate::curve_curve::{
    global_coincidence_check_with_drop, line_line_closest_point, line_plane_intersect,
    point_solution,
};
use crate::curve_global::{
    DroppedPoint, analytic_line, analytic_solutions, global_point_search,
    snap_to_specified_endpoint,
};
use crate::curve_solve::{local_curve_solve, local_curve_solve_with_last};
use crate::extent::{Extent1, Extent3};
use crate::geometry::{
    PlanePair, circle_coplanar_circle_intersect, intersect_two_planes, line_closest_point,
    line_coplanar_circle_intersect, point_project_to_plane, segment_point_distance,
    vector_project_to_plane,
};
use crate::global_solver::{SolutionRequested, SolverOperation, SolverRefusal};
use crate::line::Line;
use crate::scalar::{BIG_DOUBLE, EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, ZONE_TOL_3D, max};
use crate::solution::{SolutionArray, SolutionEnd, SolutionType};
use crate::spline::BSplineCurve;
use crate::vector::Vector3;

/// `sm_AddPoint`: an error adding a third distinct point is ignored by the
/// caller in SMLib, so the first two points stay in the array.
fn add_point(points: &mut Vec<Vector3>, point: Vector3, tolerance: f64) {
    if points
        .iter()
        .any(|p| point.distance_squared(*p) < tolerance * tolerance)
    {
        return;
    }
    if points.len() < 2 {
        points.push(point);
    }
}

/// `smgu_SegmentSegmentIntersect`, including its endpoint preference and
/// the near-parallel two-point case. The points lie on the first segment.
fn segment_intersections(a: [Vector3; 2], b: [Vector3; 2], tolerance: f64) -> Vec<Vector3> {
    let tolerance = max(tolerance, EFF_ZERO);
    let mut points = Vec::with_capacity(2);
    let va = a[1] - a[0];
    let vb = b[1] - b[0];
    for test in b {
        let (distance, parameter) = segment_point_distance(a[0], a[1], test);
        if distance <= tolerance {
            add_point(&mut points, a[0] + parameter * va, tolerance);
        }
    }
    for test in a {
        if segment_point_distance(b[0], b[1], test).0 <= tolerance {
            add_point(&mut points, test, tolerance);
        }
    }
    if points.len() > 1 {
        return points;
    }
    let Some((mut t1, mut t2)) = line_line_closest_point(a[0], va, b[0], vb) else {
        return points;
    };
    let mut truncated1 = false;
    if t1 < 0.0 {
        truncated1 = true;
        t1 = 0.0;
    }
    if t1 > 1.0 {
        truncated1 = true;
        t1 = 1.0;
    }
    let mut point = a[0] + t1 * va;
    if truncated1 {
        t2 = line_closest_point(b[0], vb, point);
    }
    let mut truncated2 = false;
    if t2 < 0.0 {
        truncated2 = true;
        t2 = 0.0;
    }
    if t2 > 1.0 {
        truncated2 = true;
        t2 = 1.0;
    }
    let other_point = b[0] + t2 * vb;
    if truncated2 {
        t1 = line_closest_point(a[0], va, other_point).clamp(0.0, 1.0);
        point = a[0] + t1 * va;
    }
    let distance = point.distance(other_point);
    if distance > tolerance {
        return points;
    }
    if let Some(first) = points.first_mut() {
        let end_distance = point.distance(*first);
        if end_distance > tolerance * 2.0 && end_distance < distance * 100.0 {
            points.push(point);
            return points;
        }
        *first = point;
    } else {
        points.push(point);
    }
    points
}

/// `SmLine::IntersectWithLine` on two distinct cached analytic lines.
/// The cache branch returns these answers without range refinement in SMLib.
pub(crate) fn line_line_intersect(
    first: &Line,
    interval: Extent1,
    second: &Line,
    other_interval: Extent1,
    tolerance: f64,
) -> SolutionArray {
    let a = [
        first.evaluate_point(interval.min),
        first.evaluate_point(interval.max),
    ];
    let b = [
        second.evaluate_point(other_interval.min),
        second.evaluate_point(other_interval.max),
    ];
    let mut solutions = SolutionArray::new();
    let mut bbox = Extent3::new(a[0], a[1]);
    bbox.expand_absolute(tolerance);
    if bbox.are_disjoint(Extent3::new(b[0], b[1]), 0.0) {
        return solutions;
    }
    let va = a[1] - a[0];
    let vb = b[1] - b[0];
    let points = if va.length_squared() > EFF_ZERO_SQ && vb.length_squared() > EFF_ZERO_SQ {
        segment_intersections(a, b, tolerance)
    } else {
        return solutions;
    };
    let to_end = |point: Vector3| {
        let t1 = line_closest_point(a[0], va, point);
        let t2 = line_closest_point(b[0], vb, point);
        let mut parameters = SolutionEnd::default();
        parameters.parameters[0] = interval.evaluate(t1);
        parameters.parameters[1] = other_interval.evaluate(t2);
        parameters.value = (a[0] + t1 * va).distance(b[0] + t2 * vb);
        parameters
    };
    if let Some(first) = points.first() {
        let start = to_end(*first);
        let mut solution =
            point_solution(None, start.parameters[0], start.parameters[1], start.value);
        if points.len() > 1 {
            solution.kind = SolutionType::RangeOfValues;
            solution.end = to_end(points[1]);
            if solution.start.parameters[0] > solution.end.parameters[0] {
                std::mem::swap(&mut solution.start, &mut solution.end);
            }
        }
        solutions.add(solution);
    }
    solutions
}

const INTERSECT: SolverOperation = SolverOperation::Intersect;
// SmMath.h's decimal macro, rounded once to binary64.
const SINE_10_DEG: f64 = 0.173_648_177_666_930_36;

/// The analytic dispatch may ask the general curve search to continue.
pub(crate) struct AnalyticIntersection {
    pub needs_more: bool,
    pub solutions: SolutionArray,
}

impl AnalyticIntersection {
    fn complete(solutions: SolutionArray) -> Self {
        Self {
            needs_more: false,
            solutions,
        }
    }
    fn more(solutions: SolutionArray) -> Self {
        Self {
            needs_more: true,
            solutions,
        }
    }
}

enum CachedCurve {
    Line(Line),
    Circle(Box<Circle>),
}

fn cached(curve: &BSplineCurve, interval: Extent1) -> Result<Option<CachedCurve>, SolverRefusal> {
    if let Some(line) = analytic_line(curve, interval)? {
        return Ok(Some(CachedCurve::Line(line)));
    }
    match analytic_arc(curve, interval)? {
        AnalyticArc::Circle(circle) => Ok(Some(CachedCurve::Circle(circle))),
        AnalyticArc::Spline => Ok(None),
        AnalyticArc::Rebuilt => Err(SolverRefusal::NotPorted(INTERSECT)),
    }
}

/// The direct `GlobalPointSolve` of a circle first calls its own fast drop,
/// without the refinement used for an analytic representation in a cache.
fn circle_meet(
    circle: &Circle,
    interval: Extent1,
    point: Vector3,
    tolerance: f64,
    requested: SolutionRequested,
) -> Result<SolutionArray, SolverRefusal> {
    let drop = circle.drop_point_fast(
        interval,
        true,
        point,
        tolerance,
        requested == SolutionRequested::Single,
    );
    let solutions = analytic_solutions(&drop.answers);
    if drop.complete {
        return Ok(solutions);
    }
    let curve = circle.spline();
    let cache = CurveCache::for_curve(curve, 0.5 * ZONE_TOL_3D)?;
    global_point_search(
        curve, &cache, interval, INTERSECT, point, tolerance, requested, solutions,
    )
}

fn line_meet(line: &Line, interval: Extent1, point: Vector3, tolerance: f64) -> SolutionArray {
    analytic_solutions(
        &line
            .drop_point_fast(interval, true, point, tolerance)
            .into_iter()
            .collect::<Vec<_>>(),
    )
}

fn circle_drop(
    circle: &Circle,
    interval: Extent1,
    point: Vector3,
    inward: Option<Vector3>,
    tolerance: f64,
) -> Result<DroppedPoint, SolverRefusal> {
    let solutions = circle_meet(
        circle,
        interval,
        point,
        tolerance,
        SolutionRequested::Single,
    )?;
    for solution in solutions.as_slice() {
        let mut parameter = solution.start.parameters[0];
        if snap_to_specified_endpoint(circle.spline(), interval, inward, &mut parameter, tolerance)?
            && solution.start.value < tolerance
        {
            return Ok(DroppedPoint {
                success: true,
                parameter,
                distance: solution.start.value,
            });
        }
    }
    Ok(DroppedPoint {
        success: false,
        parameter: 0.0,
        distance: 0.0,
    })
}

/// SmEllipse::IntersectWithLine, restricted to the circles the cache builds.
#[allow(clippy::too_many_lines)]
fn circle_line(
    circle: &Circle,
    interval: Extent1,
    line: &Line,
    other_interval: Extent1,
    tolerance: f64,
) -> Result<AnalyticIntersection, SolverRefusal> {
    let (center, normal, radius) = circle.geometry();
    let (start, direction) = line.geometry();
    let zero = EFF_ZERO * (1.0 + start.max_dimension());
    let dot = if direction.length_squared() < zero * zero {
        BIG_DOUBLE
    } else {
        let mut unit = direction;
        unit.unitize()?;
        normal.dot(unit)
    };
    let mut refine = false;
    let mut points = Vec::new();
    if dot.abs() > SINE_10_DEG {
        let parameter =
            line_plane_intersect(start, direction, center, normal).ok_or(Fault::Degenerate)?;
        points.push(start + parameter * direction);
        refine = true;
    } else {
        if dot.abs() > EFF_ZERO_SQRT {
            refine = true;
        }
        let projected_point = point_project_to_plane(start, center, normal);
        let projected_vector = vector_project_to_plane(direction, normal);
        if projected_vector.length_squared() < EFF_ZERO_SQ {
            return Ok(AnalyticIntersection::more(SolutionArray::new()));
        }
        let hits = line_coplanar_circle_intersect(
            projected_point,
            projected_vector,
            center,
            normal,
            radius,
            tolerance * 10.0,
        )?;
        for &t in &hits.parameters[..hits.count] {
            points.push(projected_point + t * projected_vector);
        }
    }
    let mut solutions = SolutionArray::new();
    for point in points {
        let line_solutions = line_meet(line, other_interval, point, tolerance * 10.0);
        if line_solutions.is_empty() {
            continue;
        }
        let circle_solutions = match circle_meet(
            circle,
            interval,
            point,
            tolerance * 10.0,
            SolutionRequested::All,
        ) {
            Ok(found) => found,
            Err(SolverRefusal::Failed(_)) => return Ok(AnalyticIntersection::more(solutions)),
            Err(refusal) => return Err(refusal),
        };
        for a in circle_solutions.as_slice() {
            for b in line_solutions.as_slice() {
                let mut parameters = [a.start.parameters[0], b.start.parameters[0]];
                let first = circle.spline().evaluate_point(parameters[0])?;
                let second = line.evaluate_point(parameters[1]);
                let mut distance = first.distance(second);
                let scaled_zero = EFF_ZERO * (1.0 + first.max_dimension());
                if refine && distance < tolerance * 10.0 && distance > tolerance / 100.0 {
                    if let Some(found) = local_curve_solve(
                        circle.spline(),
                        interval,
                        line,
                        other_interval,
                        SolverOperation::Minimize,
                        tolerance,
                        parameters,
                    )? {
                        if found.value <= tolerance + scaled_zero {
                            parameters = found.parameters;
                            distance = found.value;
                        }
                    }
                }
                if distance <= tolerance + scaled_zero {
                    solutions.add(point_solution(None, parameters[0], parameters[1], distance));
                }
            }
        }
    }
    Ok(AnalyticIntersection::complete(solutions))
}

/// SmEllipse::IntersectWithEllipse for circles, including the near-parallel
/// fallback, the coincidence check and the last-iterate refinement.
#[allow(clippy::too_many_lines)]
fn circle_circle(
    first: &Circle,
    interval: Extent1,
    second: &Circle,
    other_interval: Extent1,
    tolerance: f64,
) -> Result<AnalyticIntersection, SolverRefusal> {
    let (center1, normal1, radius1) = first.geometry();
    let (center2, normal2, radius2) = second.geometry();
    let sine = normal1.cross(normal2).length();
    let mut points = Vec::new();
    if sine > EFF_ZERO_SQRT {
        if sine < SINE_10_DEG {
            return Ok(AnalyticIntersection::more(SolutionArray::new()));
        }
        let PlanePair::Line { point, direction } =
            intersect_two_planes(center1, normal1, center2, normal2)?
        else {
            return Err(Fault::Degenerate.into());
        };
        let hits = line_coplanar_circle_intersect(
            point,
            direction,
            center1,
            normal1,
            radius1,
            tolerance * 10.0,
        )?;
        for &t in &hits.parameters[..hits.count] {
            points.push(point + t * direction);
        }
    } else {
        let projected = point_project_to_plane(center2, center1, normal1);
        if projected.distance(center2) > tolerance {
            return Ok(AnalyticIntersection::complete(SolutionArray::new()));
        }
        let (count, hits) = circle_coplanar_circle_intersect(
            normal1, center1, radius1, center2, radius2, tolerance,
        )?;
        if count > 2 {
            let drop = |which, domain, point, inward| {
                circle_drop(
                    if which { first } else { second },
                    domain,
                    point,
                    inward,
                    tolerance,
                )
            };
            let check = global_coincidence_check_with_drop(
                first.spline(),
                interval,
                second.spline(),
                other_interval,
                tolerance,
                false,
                &drop,
            )?;
            // SMLib reports complete even when the coincidence checker asks
            // for more work. Keep its candidates; no completeness proof follows.
            return Ok(AnalyticIntersection::complete(check.solutions));
        }
        points.extend_from_slice(&hits[..count]);
    }
    let mut solutions = SolutionArray::new();
    for point in points {
        let first_solutions = circle_meet(
            first,
            interval,
            point,
            tolerance * 10.0,
            SolutionRequested::All,
        )?;
        if first_solutions.is_empty() {
            continue;
        }
        let second_solutions = circle_meet(
            second,
            other_interval,
            point,
            tolerance * 10.0,
            SolutionRequested::All,
        )?;
        for a in first_solutions.as_slice() {
            for b in second_solutions.as_slice() {
                let mut parameters = [a.start.parameters[0], b.start.parameters[0]];
                let p1 = first.spline().evaluate_point(parameters[0])?;
                let p2 = second.spline().evaluate_point(parameters[1])?;
                let mut distance = p1.distance(p2);
                let scaled_zero = EFF_ZERO * (1.0 + p1.max_dimension());
                if distance < tolerance * 100.0 && distance > tolerance / 100.0 {
                    let (_, last) = local_curve_solve_with_last(
                        first.spline(),
                        interval,
                        second.spline(),
                        other_interval,
                        SolverOperation::Minimize,
                        tolerance,
                        parameters,
                    )?;
                    if interval.contains_value(last.parameters[0], EFF_ZERO)
                        && other_interval.contains_value(last.parameters[1], EFF_ZERO)
                    {
                        let p1 = first.spline().evaluate_point(last.parameters[0])?;
                        let p2 = second.spline().evaluate_point(last.parameters[1])?;
                        let improved = p1.distance(p2);
                        if improved < distance {
                            distance = improved;
                            parameters = last.parameters;
                        }
                    }
                }
                if distance <= tolerance + scaled_zero {
                    solutions.add(point_solution(None, parameters[0], parameters[1], distance));
                }
            }
        }
    }
    Ok(AnalyticIntersection::complete(solutions))
}

/// SmCurve::AnalyticCurveIntersect on the two cached representations.
#[allow(clippy::too_many_arguments)]
pub(crate) fn analytic_curve_intersect(
    first: &BSplineCurve,
    first_cache: Extent1,
    interval: Extent1,
    second: &BSplineCurve,
    second_cache: Extent1,
    other_interval: Extent1,
    tolerance: f64,
) -> Result<AnalyticIntersection, SolverRefusal> {
    let Some(a) = cached(first, first_cache)? else {
        return Ok(AnalyticIntersection::more(SolutionArray::new()));
    };
    let Some(b) = cached(second, second_cache)? else {
        return Ok(AnalyticIntersection::more(SolutionArray::new()));
    };
    match (&a, &b) {
        (CachedCurve::Line(a), CachedCurve::Line(b)) => Ok(AnalyticIntersection::complete(
            line_line_intersect(a, interval, b, other_interval, tolerance),
        )),
        (CachedCurve::Circle(a), CachedCurve::Line(b)) => {
            circle_line(a, interval, b, other_interval, tolerance)
        }
        (CachedCurve::Line(a), CachedCurve::Circle(b)) => {
            let mut result = circle_line(b, other_interval, a, interval, tolerance)?;
            for solution in result.solutions.as_mut_slice() {
                // Object IDs already name the original input order.
                solution.start.parameters.swap(0, 1);
                if solution.kind == SolutionType::RangeOfValues {
                    solution.end.parameters.swap(0, 1);
                }
            }
            Ok(result)
        }
        (CachedCurve::Circle(a), CachedCurve::Circle(b)) => {
            circle_circle(a, interval, b, other_interval, tolerance)
        }
    }
}
