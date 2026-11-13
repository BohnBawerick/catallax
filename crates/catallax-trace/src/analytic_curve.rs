// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: output arrays become owned values, statuses become
// Faults in the shared closest-line helper, and objects keep the original
// spline identities. Tolerance tests produce candidates only.
// Translated from SmLine.cpp (SmLine::IntersectWithLine), SmGeomUtility.cpp
// (smgu_SegmentSegmentIntersect, sm_AddPoint) at
// c9c5979f04cea8c38213004d4f04df0dcfb40d05; rust-brep has no port of these.
// Endpoint tests reuse the ported segment-point distance, as the donor
// smgu_is_point_on_segment does in sm_geom_utility.rs.
//! The line/line part of `SmCurve::AnalyticCurveIntersect`.

use crate::curve_curve::{line_line_closest_point, point_solution};
use crate::extent::{Extent1, Extent3};
use crate::geometry::{line_closest_point, segment_point_distance};
use crate::line::Line;
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ, max};
use crate::solution::{SolutionArray, SolutionEnd, SolutionType};
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
