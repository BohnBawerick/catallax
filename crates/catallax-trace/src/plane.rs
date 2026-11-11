// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: the plane is a B-spline surface carrying a
// `PlaneFrame` (see `spline.rs`), as SMLib's `SmPlane` derives from its
// B-spline surface; statuses become `Fault` and the solution container a
// returned value.
// Translated from SMLib SmPlane.cpp (DropPointFast, GlobalLineIntersect),
// SmExtent2d.cpp (IntersectWithInfiniteLine) and SmSurface.cpp
// (the DropPointFast branch of SmSurface::GlobalPointSolve); rust-brep
// crates/kernel/src/surface/sm_plane.rs, commit
// 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab, ports the same SMLib routine.
//! The analytic plane queries of step 7 that the solvers reach before their
//! B-spline searches.

use crate::Fault;
use crate::extent::Extent2;
use crate::global_solver::SolverOperation;
use crate::spline::BSplineSurface;
use crate::surface_eval::first_derivatives;
use crate::surface_solve::map_to_2d;
use crate::vector::{Vector2, Vector3};

/// `SmPlane::DropPointFast`: the projection of `test` onto the plane,
/// clamped into `domain`, which may extend past the plane's own domain: the
/// drop then reads the plane's analytic extension, as SMLib's does, while
/// `surface_eval::evaluate_point` clamps to the natural domain. A non-finite
/// projection or distance is refused.
/// The result carries its distance; `None` when the operation
/// finds nothing (a normal drop outside the domain, a meeting farther than
/// `tolerance`, or a drop beyond `target`). `Ok(None)` also stands for
/// SMLib's unsupported operations and for surfaces that are not planes,
/// whose `DropPointFast` fails.
pub fn drop_point_fast(
    surface: &BSplineSurface,
    domain: Extent2,
    operation: SolverOperation,
    test: Vector3,
    tolerance: f64,
    target: Option<f64>,
) -> Result<Option<(Vector2, f64)>, Fault> {
    let Some(frame) = surface.plane_frame() else {
        return Err(Fault::InvalidInput);
    };
    if !matches!(
        operation,
        SolverOperation::Minimize | SolverOperation::Normalize | SolverOperation::Intersect
    ) {
        return Err(Fault::InvalidInput);
    }
    let mut uv = frame.project(test);
    if !(uv.x.is_finite() && uv.y.is_finite()) {
        return Err(Fault::InvalidInput);
    }
    if operation == SolverOperation::Normalize {
        // `SmTol::MapTo2d` without a direction maps along the longer first
        // derivative.
        let [_, du, dv] =
            first_derivatives(surface, uv, [true, true]).unwrap_or([Vector3::default(); 3]);
        let direction = if du.length_squared() > dv.length_squared() {
            Vector2::new(1.0, 0.0)
        } else {
            Vector2::new(0.0, 1.0)
        };
        let tolerance_uv = map_to_2d(surface, tolerance, uv, direction);
        if !domain.contains_point(uv, tolerance_uv) {
            return Ok(None);
        }
    }
    uv = domain.clamp_point(uv);
    let distance = frame.point_fast(uv).distance(test);
    if !distance.is_finite() {
        return Err(Fault::InvalidInput);
    }
    if operation == SolverOperation::Intersect && distance > tolerance {
        return Ok(None);
    }
    if let Some(target) = target
        && matches!(
            operation,
            SolverOperation::Minimize | SolverOperation::Normalize
        )
        && distance > target
    {
        return Ok(None);
    }
    Ok(Some((uv, distance)))
}

/// `SmPlane::GlobalLineIntersect` (not as a ray): the meeting of the line
/// `point + t * vector`, over `interval` when given, with the plane over
/// `domain` within `tolerance`, as solutions `t, u, v`. A line crossing the
/// plane meets it once (or not at all outside the interval or beyond the
/// tolerance after clamping into the domain). A line parallel to the plane
/// meets it nowhere unless it lies within the tolerance over the domain,
/// where the answer is the range of the line over the domain; the degenerate
/// overlap SMLib hands to the general line search refuses with
/// [`SolverRefusal::NotPorted`]. Unlike SMLib, a nearly parallel line with
/// a trimmed end beyond the tolerance and the other end on the opposite
/// side of the plane is answered by its crossing, and one with a single end
/// within the tolerance refuses instead of meeting nothing.
pub fn global_line_intersect(
    surface: &BSplineSurface,
    domain: Extent2,
    point: Vector3,
    vector: Vector3,
    interval: Option<crate::extent::Extent1>,
    tolerance: f64,
) -> Result<crate::solution::SolutionArray, crate::global_solver::SolverRefusal> {
    use crate::extent::Extent1;
    use crate::geometry::line_closest_point;
    use crate::global_solver::SolverRefusal;
    use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ};
    use crate::solution::{Solution, SolutionArray, SolutionType};
    let frame = surface.plane_frame().ok_or(Fault::InvalidInput)?;
    if vector.length_squared() < EFF_ZERO_SQ {
        return Err(Fault::InvalidInput.into());
    }
    let mut solutions = SolutionArray::new();
    let normal = frame.z_axis();
    if normal.is_perpendicular_to(vector, 0.0001) {
        let span = match interval {
            Some(interval) => interval,
            None => {
                let mut span = Extent1::default();
                let (lo, hi) = (domain.min, domain.max);
                for uv in [lo, Vector2::new(lo.x, hi.y), Vector2::new(hi.x, lo.y), hi] {
                    span.add_value(line_closest_point(point, vector, frame.point_fast(uv)));
                }
                let grow = (span.max - span.min) * 1.2;
                Extent1::new(span.min - grow, span.max + grow)
            }
        };
        let start = point + span.min * vector;
        let end = point + span.max * vector;
        let uv_start = frame.project(start);
        let uv_vector = frame.project(end) - uv_start;
        let Some(trim) = intersect_with_infinite_line(domain, uv_start, uv_vector) else {
            return Ok(solutions);
        };
        let low = crate::scalar::max(trim.min, 0.0);
        let high = crate::scalar::min(trim.max, 1.0);
        if low > high + EFF_ZERO {
            return Ok(solutions);
        }
        let on_plane_start = frame.point_fast(uv_start + low * uv_vector);
        let on_plane_end = frame.point_fast(uv_start + high * uv_vector);
        let t_start = line_closest_point(point, vector, on_plane_start);
        let t_end = line_closest_point(point, vector, on_plane_end);
        let gap_start = point + t_start * vector - on_plane_start;
        let gap_end = point + t_end * vector - on_plane_end;
        let far_start = gap_start.length() > tolerance;
        let far_end = gap_end.length() > tolerance;
        // Deliberate deviation from SMLib, which returns no meeting as soon
        // as one end of the trimmed line is beyond the tolerance: the
        // parallel test admits lines with a slight slope, which can still
        // cross the plane between the ends. A far end with the other end on
        // the opposite side of the plane goes to the crossing below; ends on
        // one side bound the gap (it is linear along the line), so two far
        // ends meet nothing and one near end is a meeting this branch cannot
        // place, refused.
        let crossing = (far_start || far_end) && gap_start.dot(normal) * gap_end.dot(normal) < 0.0;
        if !crossing {
            if far_start && far_end {
                return Ok(solutions);
            }
            if far_start || far_end {
                // A meeting at one end only.
                return Err(SolverRefusal::NotPorted(SolverOperation::Intersect));
            }
            if crate::scalar::are_same(t_start, t_end) {
                // SMLib's general line search over the plane's patches.
                return Err(SolverRefusal::NotPorted(SolverOperation::Intersect));
            }
            // The line lies on the plane over its trimmed span: a range.
            let uv_start = frame.project(on_plane_start);
            let uv_end = frame.project(on_plane_end);
            let start_distance = (point + t_start * vector).distance(on_plane_start);
            let end_distance = (point + t_end * vector).distance(on_plane_end);
            let mut solution = Solution {
                kind: SolutionType::RangeOfValues,
                variables: 3,
                object_count: 0,
                ..Solution::default()
            };
            solution.start.value = start_distance;
            solution.start.parameters[..3].copy_from_slice(&[
                span.clamp_value(t_start),
                uv_start.x,
                uv_start.y,
            ]);
            solution.end.value = end_distance;
            solution.end.parameters[..3].copy_from_slice(&[
                span.clamp_value(t_end),
                uv_end.x,
                uv_end.y,
            ]);
            solutions.add(solution);
            return Ok(solutions);
        }
    }
    let span = interval.unwrap_or(Extent1::new(-1.0e10, 1.0e10));
    let pu = frame.x_axis * frame.uv_scale.x;
    let pv = frame.y_axis * frame.uv_scale.y;
    let lt = -vector;
    let difference = point - frame.origin;
    let determinant = pu.cross(pv).dot(lt);
    if determinant.abs() < EFF_ZERO {
        return Err(Fault::Degenerate.into());
    }
    let t = pu.cross(pv).dot(difference) / determinant;
    if !span.contains_value(t, 0.0) {
        return Ok(solutions);
    }
    let uv = Vector2::new(
        difference.cross(pv).dot(lt) / determinant,
        pu.cross(difference).dot(lt) / determinant,
    );
    let on_line = point + t * vector;
    let uv = domain.clamp_point(uv);
    let distance = frame.point_fast(uv).distance(on_line);
    if distance > tolerance {
        return Ok(solutions);
    }
    let mut solution = Solution {
        kind: SolutionType::SingleValue,
        variables: 3,
        object_count: 0,
        ..Solution::default()
    };
    solution.start.value = distance;
    solution.start.parameters[..3].copy_from_slice(&[t, uv.x, uv.y]);
    solutions.add(solution);
    Ok(solutions)
}

/// `SmExtent2d::IntersectWithInfiniteLine` with the default tolerance: the
/// line parameters where `point + s * vector` crosses the domain's sides
/// (within `100 * SM_EFF_ZERO`), if it meets the domain.
fn intersect_with_infinite_line(
    domain: Extent2,
    point: Vector2,
    vector: Vector2,
) -> Option<crate::extent::Extent1> {
    use crate::scalar::EFF_ZERO;
    let tol = EFF_ZERO * 100.0;
    let (min_x, max_x) = (domain.min.x - tol, domain.max.x + tol);
    let (min_y, max_y) = (domain.min.y - tol, domain.max.y + tol);
    let mut span = crate::extent::Extent1::default();
    if vector.x.abs() > EFF_ZERO {
        for side in [domain.min.x, domain.max.x] {
            let s = (side - point.x) / vector.x;
            let at = point.y + s * vector.y;
            if at >= min_y && at <= max_y {
                span.add_value(s);
            }
        }
    }
    if vector.y.abs() > EFF_ZERO {
        for side in [domain.min.y, domain.max.y] {
            let s = (side - point.y) / vector.y;
            let at = point.x + s * vector.x;
            if at >= min_x && at <= max_x {
                span.add_value(s);
            }
        }
    }
    (span.max > span.min - tol).then_some(span)
}
