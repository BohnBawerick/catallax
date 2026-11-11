// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline curves, surfaces and
// planes of `spline.rs` without trimming faces; statuses become `Fault`, the
// trace-point list becomes a vector and output arrays become a returned
// struct. Paths that need code from later port steps refuse with
// `SolverRefusal::NotPorted`: a curve not known to lie on the surface (the
// projected side search), the tracer's fallback when the general curve drop
// (`SmSurface::DropCurve`, in `surface_drop_general.rs`) fails or gives more
// than one curve (`SmSurfaceTracer::TraceCurve`, `FlushCurve`), a
// start point at a seam of a closed surface or on a singularity without a
// previous point, and a plane's side lines that meet the curve's box.
// Translated from SMLib SmSurfaceDropLib.cpp (SmSurface::DropAndTrimCurve,
// SmSurface::DropIsoCurve), SmSurfaceDropCurve.cpp (SmSurfaceDropCurve:
// constructor, DoTrace, FindBoundaryStartPoints, sm_AddBreak,
// sm_IsGrazingIntersection, TraceCurve2, FindInteriorCurves,
// ComputePointValues), SmSurfaceTracer.cpp (constructor, AddStartPoint),
// SmPlane.cpp (SmPlane::DropCurve), SmSurface.cpp (SnapToKnots,
// GetDomainWithContinuity), SmCurve.cpp (SnapToKnots) and SmExtent2d.cpp
// (IntersectWithInfiniteLine with tolerances); rust-brep has no port.
//! The curve drop of the range answers: the parameter curves of a 3D curve
//! lying on a surface, split where it leaves the surface's domain, as SMLib's
//! `DropAndTrimCurve` builds them for coincident curves.

use crate::Fault;
use crate::curve_cache::CurveCache;
use crate::curve_curve::global_curve_intersect;
use crate::curve_surface::{drop_vector, is_singularity};
use crate::extent::{Extent1, Extent2};
use crate::geometry::{Continuity, line_closest_point, line_point_distance};
use crate::global_solver::{SolverOperation, SolverRefusal};
use crate::line::Line;
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, PI, ZONE_TOL_3D, are_same, max};
use crate::solution::SolutionType;
use crate::spline::{BSplineCurve, BSplineSurface, ControlPoint, NO_Z, PointForm, unique_knots};
use crate::surface_cache::{CacheRefusal, SurfaceCache};
use crate::surface_eval::{evaluate, evaluate_normal, evaluate_point, first_derivatives};
use crate::vector::{Vector2, Vector3};

const INTERSECT: SolverOperation = SolverOperation::Intersect;

/// A parameter curve of a dropped curve: a B-spline in the plane `z = 0`, or
/// a plane's straight drop (`SmLine`).
#[derive(Clone, Debug)]
pub enum UvCurve {
    /// A B-spline parameter curve.
    Spline(BSplineCurve),
    /// A straight parameter curve on a plane.
    Line(Line),
}

impl UvCurve {
    /// `GetNaturalInterval`.
    pub fn natural_interval(&self) -> Extent1 {
        match self {
            Self::Spline(curve) => curve.natural_interval(),
            Self::Line(line) => line.natural_interval(),
        }
    }

    /// `EvaluatePoint`.
    pub fn evaluate_point(&self, t: f64) -> Result<Vector3, Fault> {
        match self {
            Self::Spline(curve) => curve.evaluate_point(t),
            Self::Line(line) => Ok(line.evaluate_point(t)),
        }
    }
}

/// What [`drop_and_trim_curve`] found: `rUVCurves`, `rdMaxDropToSurf` and
/// `rdMaxApproxDev`.
#[derive(Clone, Debug, Default)]
pub struct DroppedCurves {
    /// The parameter curves, one per piece of the curve over the surface.
    pub curves: Vec<UvCurve>,
    /// The largest sampled distance from the curve to the surface.
    pub max_drop: f64,
    /// The largest sampled distance from the curve to the surface normal
    /// lines.
    pub max_deviation: f64,
}

/// `SmTracePnt` as the curve drop fills it.
#[derive(Clone, Copy, Debug, Default)]
struct TracePoint {
    curve_parameter: f64,
    surface_point: Vector3,
    du: Vector3,
    dv: Vector3,
    normal: Vector3,
    curve_tangent: Vector3,
    uv: Vector2,
    uv_tangent: Vector2,
    drop_distance: f64,
    normal_line_distance: f64,
}

/// `sm_AddBreak`: inserts a break keeping the list sorted and free of
/// values within a hundred scaled zeros of each other.
fn add_break(breaks: &mut Vec<f64>, value: f64) {
    let zero = EFF_ZERO * 100.0 * (1.0 + value.abs());
    let mut insertion = breaks.len();
    for (i, &existing) in breaks.iter().enumerate() {
        if (value - existing).abs() < zero {
            return;
        }
        if value < existing {
            insertion = i;
            break;
        }
    }
    breaks.insert(insertion, value);
}

/// `SmCurve::SnapToKnots`.
fn snap_curve_to_knots(curve: &BSplineCurve, value: f64, tolerance: f64) -> f64 {
    let mut snapped = value;
    let mut jump = crate::scalar::BIG_DOUBLE;
    for (knot, _) in unique_knots(curve.knots(), None) {
        let test = (knot - value).abs();
        if test < tolerance && test < jump {
            snapped = knot;
            jump = test;
        }
    }
    snapped
}

/// The unique knots of a surface along `axis` (`GetKnots`).
fn surface_knots(surface: &BSplineSurface, axis: usize) -> Vec<f64> {
    unique_knots(surface.knots()[axis], None)
        .into_iter()
        .map(|(k, _)| k)
        .collect()
}

/// `SmSurface::SnapToKnots`.
fn snap_surface_to_knots(surface: &BSplineSurface, uv: Vector2, tolerance: f64) -> Vector2 {
    let mut snapped = uv;
    for axis in 0..2 {
        let value = if axis == 0 { uv.x } else { uv.y };
        let mut best = crate::scalar::BIG_DOUBLE;
        for knot in surface_knots(surface, axis) {
            let distance = (knot - value).abs();
            if distance < best && distance < tolerance {
                if axis == 0 {
                    snapped.x = knot;
                } else {
                    snapped.y = knot;
                }
                best = distance;
            }
        }
    }
    snapped
}

fn surface_cache(surface: &BSplineSurface) -> Result<SurfaceCache, SolverRefusal> {
    SurfaceCache::for_trimmed_surface(surface).map_err(|refusal| match refusal {
        CacheRefusal::Failed(fault) => SolverRefusal::Failed(fault),
        CacheRefusal::NotPorted(_) => SolverRefusal::NotPorted(INTERSECT),
    })
}

/// `SmSurface::GetDomainWithContinuity`: the domain around `uv` (on the
/// side `vector` points to) bounded by knot lines below `minimum`. SMLib
/// compares the end knots with the vector's components and reads the
/// unclamped `v` component at an interior `v` knot; both are kept.
fn domain_with_continuity(
    surface: &BSplineSurface,
    cache: &SurfaceCache,
    minimum: Continuity,
    uv: Vector2,
    vector: Vector2,
) -> Extent2 {
    let natural = surface.natural_domain();
    let (mut lo, mut hi) = (natural.min, natural.max);
    let mut clamped = vector;
    if clamped.x.abs() < EFF_ZERO * 1000.0 {
        clamped.x = 0.0;
    }
    if clamped.y.abs() < EFF_ZERO * 1000.0 {
        clamped.y = 0.0;
    }
    for axis in 0..2 {
        let knots = surface_knots(surface, axis);
        let conts = if axis == 0 {
            &cache.u_continuities
        } else {
            &cache.v_continuities
        };
        let (value, component, raw) = if axis == 0 {
            (uv.x, clamped.x, clamped.x)
        } else {
            (uv.y, clamped.y, vector.y)
        };
        let (low, high) = if axis == 0 {
            (&mut lo.x, &mut hi.x)
        } else {
            (&mut lo.y, &mut hi.y)
        };
        for i in 0..conts.len() {
            if conts[i] >= minimum {
                continue;
            }
            let knot = knots.get(i).copied().unwrap_or(0.0);
            if i == 0 && component < knot {
                *low = knot;
            } else if i == conts.len() - 1 && component > knot {
                *high = knot;
            } else if (knot - value).abs() < EFF_ZERO_SQRT {
                if raw <= 0.0 && i > 0 {
                    *high = knot;
                    break;
                } else if i == knots.len() - 1 {
                    *high = knot;
                } else {
                    *low = knot;
                }
            } else if knot < value {
                *low = knot;
            } else {
                *high = knot;
                break;
            }
        }
    }
    Extent2 { min: lo, max: hi }
}

/// `SmExtent2d::IntersectWithInfiniteLine` with tolerances in `u` and `v`.
fn intersect_with_infinite_line(
    domain: Extent2,
    point: Vector2,
    vector: Vector2,
    tol_u: f64,
    tol_v: f64,
) -> Option<Extent1> {
    let (min_x, max_x) = (domain.min.x - tol_u, domain.max.x + tol_u);
    let (min_y, max_y) = (domain.min.y - tol_v, domain.max.y + tol_v);
    let mut span = Extent1::default();
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
    (span.max > span.min - max(tol_u, tol_v)).then_some(span)
}

/// What a surface's general curve drop returned.
pub(crate) struct Drop {
    pub(crate) curves: Vec<UvCurve>,
    pub(crate) max_drop: f64,
    pub(crate) deviation: f64,
}

/// `SmSurface::DropCurve` (virtual, so `SmPlane::DropCurve` on a plane):
/// the parameter curves of `interval` of `curve` dropped onto `surface`
/// within `domain`, with the largest sampled drop distance and normal-line
/// deviation. Unless `keep_all`, curves that leave the tolerance are
/// discarded. The curves are unproved candidates.
pub fn drop_curve_onto(
    surface: &BSplineSurface,
    domain: Extent2,
    curve: &BSplineCurve,
    interval: Extent1,
    tolerance: f64,
    keep_all: bool,
) -> Result<DroppedCurves, SolverRefusal> {
    let dropped = drop_curve(surface, domain, curve, interval, tolerance, keep_all)?;
    Ok(DroppedCurves {
        curves: dropped.curves,
        max_drop: dropped.max_drop,
        max_deviation: dropped.deviation,
    })
}

/// `SmSurface::DropCurve`: a plane's drop (`SmPlane::DropCurve`), which
/// hands a curve off the plane to the general drop, or the general drop
/// ([`crate::surface_drop_general`]).
pub(crate) fn drop_curve(
    surface: &BSplineSurface,
    domain: Extent2,
    curve: &BSplineCurve,
    interval: Extent1,
    tolerance: f64,
    keep_all: bool,
) -> Result<Drop, SolverRefusal> {
    let Some(frame) = surface.plane_frame() else {
        return crate::surface_drop_general::general_drop_curve(
            surface, domain, curve, interval, tolerance, keep_all,
        );
    };
    let mut result = Drop {
        curves: Vec::new(),
        max_drop: 0.0,
        deviation: 0.0,
    };
    if curve.is_line(5, EFF_ZERO)? {
        let mut start3 = curve.evaluate_point(interval.min)?;
        let mut end3 = curve.evaluate_point(interval.max)?;
        let mut start2 = frame.project(start3);
        let mut end2 = frame.project(end3);
        let vector2 = end2 - start2;
        let scale = frame.uv_scale;
        let Some(trim) = intersect_with_infinite_line(
            domain,
            start2,
            vector2,
            tolerance / scale.x,
            tolerance / scale.y,
        ) else {
            return Ok(result);
        };
        let Ok(trim) = trim.intersect(Extent1::new(0.0, 1.0)) else {
            return Ok(result);
        };
        let mut final_interval = interval;
        if !trim.equals_ends(0.0, 1.0, EFF_ZERO) {
            end2 = start2 + vector2 * trim.max;
            start2 = start2 + vector2 * trim.min;
            let _ = final_interval.set_min_max(
                final_interval.evaluate(trim.min),
                final_interval.evaluate(trim.max),
            );
            start3 = curve.evaluate_point(final_interval.min)?;
            end3 = curve.evaluate_point(final_interval.max)?;
        }
        let distance1 = frame.point_fast(start2).distance(start3);
        let distance2 = frame.point_fast(end2).distance(end3);
        result.max_drop = max(distance1, distance2);
        if keep_all || result.max_drop <= tolerance {
            let mut direction = end2 - start2;
            let length = direction.length();
            if length > EFF_ZERO {
                direction = Vector2::new(direction.x / length, direction.y / length);
                let mut line = Line::new(
                    Vector3::new(start2.x, start2.y, 0.0),
                    Vector3::new(direction.x, direction.y, 0.0),
                    Extent1::new(0.0, 1.0),
                    length,
                    true,
                    false,
                );
                line.edit_parameterization(final_interval);
                result.curves.push(UvCurve::Line(line));
                result.deviation = 0.0;
            }
            return Ok(result);
        }
    }
    // `Trim` snaps the interval to nearby knots; the general drop gets the
    // snapped interval.
    let (copy, interval) = match curve.trim(interval, ZONE_TOL_3D) {
        Ok(trimmed) => (trimmed.curve, trimmed.interval),
        Err(crate::curve_edit::TrimRefusal::Failed(_)) => (curve.clone(), interval),
        Err(crate::curve_edit::TrimRefusal::ShortSpan) => {
            return Err(SolverRefusal::NotPorted(INTERSECT));
        }
    };
    let mut max_drop_sq = 0.0;
    let mut projected = Vec::with_capacity(copy.control_points().len());
    for control in copy.control_points() {
        let (point, weight) = control.get(PointForm::Euclidean);
        let uv = domain.clamp_point(frame.project(point));
        let on_surface = frame.point_fast(uv);
        let drop_sq = on_surface.distance_squared(point);
        if drop_sq > max_drop_sq {
            max_drop_sq = drop_sq;
        }
        let normal = evaluate_normal(surface, uv, [true, true]);
        let (deviation, _) = line_point_distance(on_surface, normal, point);
        if deviation > tolerance {
            return crate::surface_drop_general::general_drop_curve(
                surface, domain, curve, interval, tolerance, keep_all,
            );
        }
        projected.push((uv, weight));
    }
    if keep_all || max_drop_sq < tolerance * tolerance {
        let knots = unique_knots(copy.knots(), None);
        let values: Vec<f64> = knots.iter().map(|k| k.0).collect();
        let multiplicities: Vec<usize> = knots.iter().map(|k| k.1).collect();
        let points: Vec<Vector3> = projected
            .iter()
            .map(|(uv, _)| Vector3::new(uv.x, uv.y, 0.0))
            .collect();
        let weights: Vec<f64> = projected.iter().map(|(_, w)| *w).collect();
        let uv_curve = BSplineCurve::create_canonical(
            2,
            copy.degree(),
            &points,
            &multiplicities,
            &values,
            copy.is_rational().then_some(weights.as_slice()),
            None,
        )?;
        result.curves.push(UvCurve::Spline(uv_curve));
        result.max_drop = max_drop_sq.sqrt();
        result.deviation = 0.0;
    }
    Ok(result)
}

/// `SmSurface::DropIsoCurve`: whether `interval` of the curve runs along a
/// parameter line of the surface from `start_uv` within `tolerance` of the
/// normal lines at evenly spaced samples, and that parameter line.
#[allow(clippy::too_many_arguments)]
pub(crate) fn drop_iso_curve(
    surface: &BSplineSurface,
    domain: Extent2,
    curve: &BSplineCurve,
    interval: Extent1,
    start_uv: Vector2,
    tolerance: f64,
) -> Result<Option<(BSplineCurve, f64, f64)>, SolverRefusal> {
    let one_degree = PI / 180.0;
    let left = interval.min <= interval.mid();
    let mut pv = curve.evaluate(interval.min, 1, left, true)?;
    if pv[1].length_squared() < EFF_ZERO_SQ {
        pv = curve.evaluate(interval.evaluate(EFF_ZERO_SQRT), 1, left, false)?;
        if pv[1].length_squared() < EFF_ZERO_SQ {
            pv = curve.evaluate(interval.evaluate(1.0e-3), 1, left, false)?;
            if pv[1].length_squared() < EFF_ZERO_SQ {
                return Ok(None);
            }
        }
    }
    let tangent = pv[1];
    let mut found = None;
    let (mut du, mut dv) = (Vector3::default(), Vector3::default());
    for k in 0..4 {
        let (left_u, left_v) = match k {
            1 => (true, false),
            2 => (false, false),
            3 => (false, true),
            _ => (true, true),
        };
        let [_, u_derivative, v_derivative] =
            first_derivatives(surface, start_uv, [left_u, left_v])?;
        du = u_derivative;
        dv = v_derivative;
        let mut normal = du.cross(dv);
        let length = normal.length();
        if length > EFF_ZERO {
            normal /= length;
        } else {
            normal = evaluate_normal(surface, start_uv, [left_u, left_v]);
        }
        let projected = tangent - normal * tangent.dot(normal);
        let mut angle_u = PI / 2.0;
        let mut angle_v = PI / 2.0;
        if du.length_squared() > EFF_ZERO_SQ {
            angle_u = du.angle_between(projected)?;
        }
        if dv.length_squared() > EFF_ZERO_SQ {
            angle_v = dv.angle_between(projected)?;
        }
        let pick = if angle_u < one_degree {
            Some((true, Vector2::new(1.0, 0.0), true))
        } else if angle_u > PI - one_degree {
            Some((false, Vector2::new(-1.0, 0.0), true))
        } else if angle_v < one_degree {
            Some((true, Vector2::new(0.0, 1.0), false))
        } else if angle_v > PI - one_degree {
            Some((false, Vector2::new(0.0, -1.0), false))
        } else {
            None
        };
        if let Some((left_sides, direction, along_u)) = pick {
            let [_, u_derivative, v_derivative] =
                first_derivatives(surface, start_uv, [left_sides, left_sides])?;
            du = u_derivative;
            dv = v_derivative;
            let scale = if along_u { du.length() } else { dv.length() } / tangent.length();
            found = Some((direction, scale));
            break;
        }
    }
    let Some((direction, scale)) = found else {
        return Ok(None);
    };
    let mut knots: Vec<f64> = unique_knots(curve.knots(), None)
        .into_iter()
        .map(|(k, _)| k)
        .collect();
    while !knots.is_empty() && knots[0] < interval.min - EFF_ZERO_SQRT {
        knots.remove(0);
    }
    while knots
        .last()
        .is_some_and(|&k| k > interval.max + EFF_ZERO_SQRT)
    {
        knots.pop();
    }
    if knots.is_empty() || !are_same(knots[0], interval.min) {
        knots.insert(0, interval.min);
    }
    if knots.len() == 1 || !are_same(*knots.last().unwrap_or(&0.0), interval.max) {
        knots.push(interval.max);
    }
    let size = interval.max - interval.min;
    let (mut max_deviation, mut max_drop) = (0.0_f64, 0.0_f64);
    let count = knots.len();
    for i in 1..=count {
        let step = if i < count {
            (i as f64 / count as f64) * size
        } else {
            size
        };
        let test_uv = domain.clamp_point(start_uv + direction * step / scale);
        let parameter = interval.clamp_value(interval.min + step);
        let curve_point = curve.evaluate_point(parameter)?;
        let surface_point = evaluate_point(surface, test_uv)?;
        let normal = evaluate_normal(surface, test_uv, [true, true]);
        let (distance, mut along) = line_point_distance(surface_point, normal, curve_point);
        if distance > tolerance {
            return Ok(None);
        }
        if along < 0.0 {
            along = -along;
        }
        if max_deviation < distance {
            max_deviation = distance;
        }
        if max_drop < along {
            max_drop = along;
        }
    }
    let start = Vector3::new(start_uv.x, start_uv.y, 0.0);
    let end_uv = domain.clamp_point(start_uv + direction * size / scale);
    let end = Vector3::new(end_uv.x, end_uv.y, 0.0);
    if start.distance(end) < EFF_ZERO_SQRT
        && (du.length_squared() <= EFF_ZERO || dv.length_squared() <= EFF_ZERO)
    {
        return Ok(None);
    }
    // `SmBSplineCurve::CreateLineSegment` refuses coincident ends.
    if start.distance(end) < EFF_ZERO {
        return Err(Fault::InvalidInput.into());
    }
    // A two-dimensional curve stores no z (`NL_NOZ`).
    let flat = |p: Vector3| Vector3::new(p.x, p.y, NO_Z);
    let line = BSplineCurve::create_canonical(
        2,
        1,
        &[flat(start), flat(end)],
        &[2, 2],
        &[0.0, 1.0],
        None,
        None,
    )?;
    let line = line.edit_parameterization(interval)?;
    // Deliberate deviation from SMLib: its knot samples alone accept a
    // curved curve whose samples happen to lie on the parameter line; the
    // line is kept only when a span-wide bound proves it (see
    // `iso_line_is_proven`), else the curve takes the general path.
    if !iso_line_is_proven(surface, curve, interval, start_uv, end_uv, tolerance)? {
        return Ok(None);
    }
    Ok(Some((line, max_drop, max_deviation)))
}

/// `curve` trimmed to `interval`, or `None` where the bound below cannot be
/// formed (a failed trim or weights other than one).
fn trimmed_points(curve: &BSplineCurve, interval: Extent1) -> Option<BSplineCurve> {
    let piece = curve.trim(interval, ZONE_TOL_3D).ok()?.curve;
    // Unit weights only (a rational flag with unit weights is polynomial):
    // the hull bound compares Euclidean points.
    piece
        .control_points()
        .iter()
        .all(|p| p.get(PointForm::Euclidean).1 == 1.0)
        .then_some(piece)
}

/// `curve` with its parameter direction reversed over the same interval.
fn reversed(curve: &BSplineCurve) -> Result<BSplineCurve, Fault> {
    let knots = curve.knots();
    let (a, b) = (knots[0], knots[knots.len() - 1]);
    let reversed_knots: Vec<f64> = knots.iter().rev().map(|k| a + b - k).collect();
    let points: Vec<_> = curve.control_points().iter().rev().copied().collect();
    BSplineCurve::new(curve.degree(), reversed_knots, points)
}

/// Whether the parameter line from `start_uv` to `end_uv` over `interval`
/// provably keeps `curve` within `tolerance` at every parameter: the
/// difference between the curve and the line's image is a spline whose
/// control points bound it (convex hull). On a plane the image is affine in
/// the parameter, so its control points are its values at the Greville
/// abscissae; on a B-spline surface the image is the parameter curve,
/// compared when it shares the curve's degree and knots. Weights other than
/// one and mismatched knots are not proven.
fn iso_line_is_proven(
    surface: &BSplineSurface,
    curve: &BSplineCurve,
    interval: Extent1,
    start_uv: Vector2,
    end_uv: Vector2,
    tolerance: f64,
) -> Result<bool, SolverRefusal> {
    let Some(piece) = trimmed_points(curve, interval) else {
        return Ok(false);
    };
    let degree = piece.degree();
    let knots = piece.knots();
    let span = piece.natural_interval();
    let length = span.max - span.min;
    if length <= 0.0 {
        return Ok(false);
    }
    if let Some(frame) = surface.plane_frame() {
        for (i, control) in piece.control_points().iter().enumerate() {
            let greville = knots[i + 1..=i + degree].iter().sum::<f64>() / degree.max(1) as f64;
            let s = (greville - span.min) / length;
            let uv = start_uv + (end_uv - start_uv) * s;
            if frame.point_fast(uv).distance(control.euclidean()) > tolerance {
                return Ok(false);
            }
        }
        return Ok(true);
    }
    let along_u = (end_uv.x - start_uv.x).abs() > (end_uv.y - start_uv.y).abs();
    let (axis, value, from, to) = if along_u {
        (1, start_uv.y, start_uv.x, end_uv.x)
    } else {
        (0, start_uv.x, start_uv.y, end_uv.y)
    };
    let iso = surface.iso_curve(axis, value)?;
    let Some(iso) = trimmed_points(&iso, Extent1::new(from, to)) else {
        return Ok(false);
    };
    let iso = if to < from { reversed(&iso)? } else { iso };
    let Ok(iso) = iso.edit_parameterization(span) else {
        return Ok(false);
    };
    if iso.degree() != degree || iso.knots().len() != knots.len() {
        return Ok(false);
    }
    // Paired control points bound the difference of two splines only on
    // one basis, so knot vectors that differ (even by rounding) are first
    // refined into their common vector by knot insertion.
    let Some((iso_points, piece_points)) = common_basis(&iso, &piece)? else {
        return Ok(false);
    };
    Ok(iso_points
        .iter()
        .zip(&piece_points)
        .all(|(a, b)| a.euclidean().distance(b.euclidean()) <= tolerance))
}

/// The control points of two curves of one degree refined to their common
/// knot vector (each gets the other's knots it lacks, counted with
/// multiplicity), or `None` when their end knots differ or the refined
/// vectors are not identical.
#[allow(clippy::type_complexity)]
fn common_basis(
    a: &BSplineCurve,
    b: &BSplineCurve,
) -> Result<Option<(Vec<ControlPoint>, Vec<ControlPoint>)>, SolverRefusal> {
    use crate::curve_build::refine_parts;
    let degree = a.degree();
    let ka = a.knots();
    if degree != b.degree() || ka.len() < 2 || b.knots().len() < 2 {
        return Ok(None);
    }
    // `EditParameterization` leaves ends within the scaled zero of the
    // requested ones; such ends are mapped exactly onto `a`'s (an affine
    // relabelling, changing the curve by rounding only).
    let (a0, a1) = (ka[0], ka[ka.len() - 1]);
    let (b0, b1) = (b.knots()[0], b.knots()[b.knots().len() - 1]);
    let same_end = |x: f64, y: f64| (x - y).abs() <= EFF_ZERO * (1.0 + x.abs());
    if !same_end(a0, b0) || !same_end(a1, b1) || b1 <= b0 {
        return Ok(None);
    }
    let mapped: Vec<f64> = b
        .knots()
        .iter()
        .map(|&k| {
            if k == b0 {
                a0
            } else if k == b1 {
                a1
            } else {
                a0 + (k - b0) * ((a1 - a0) / (b1 - b0))
            }
        })
        .collect();
    let b = &BSplineCurve::new(degree, mapped, b.control_points().to_vec())?;
    let kb = b.knots();
    // The knots of `from` that `into` lacks, as a sorted multiset.
    let missing = |from: &[f64], into: &[f64]| -> Vec<f64> {
        let (mut i, mut j, mut out) = (0, 0, Vec::new());
        while i < from.len() {
            if j < into.len() && into[j] == from[i] {
                i += 1;
                j += 1;
            } else if j < into.len() && into[j] < from[i] {
                j += 1;
            } else {
                out.push(from[i]);
                i += 1;
            }
        }
        out
    };
    let refine = |curve: &BSplineCurve,
                  insert: Vec<f64>|
     -> Result<Option<(Vec<f64>, Vec<ControlPoint>)>, Fault> {
        if insert.is_empty() {
            return Ok(Some((
                curve.knots().to_vec(),
                curve.control_points().to_vec(),
            )));
        }
        match refine_parts(degree, curve.knots(), curve.control_points(), &insert) {
            Ok(parts) => Ok(Some(parts)),
            Err(_) => Ok(None),
        }
    };
    let (Some((ra, pa)), Some((rb, pb))) =
        (refine(a, missing(kb, ka))?, refine(b, missing(ka, kb))?)
    else {
        return Ok(None);
    };
    if ra != rb || pa.len() != pb.len() {
        return Ok(None);
    }
    Ok(Some((pa, pb)))
}

/// `SmSurfaceDropCurve`.
struct DropTracer<'a> {
    surface: &'a BSplineSurface,
    domain: Extent2,
    curve: &'a BSplineCurve,
    interval: Extent1,
    tolerance: f64,
    working_interval: Extent1,
    working_domain: Extent2,
    max_distance: f64,
    breaks: Vec<f64>,
    start_points: Vec<TracePoint>,
    curves: Vec<UvCurve>,
    drops: Vec<f64>,
    deviations: Vec<f64>,
}

impl DropTracer<'_> {
    /// `SmSurfaceDropCurve::ComputePointValues` for a start point (no
    /// previous point).
    fn compute_point_values(
        &self,
        uv: Vector2,
        point: &mut TracePoint,
    ) -> Result<(), SolverRefusal> {
        let left = point.curve_parameter <= self.working_interval.mid();
        let mut pv = self.curve.evaluate(point.curve_parameter, 1, left, false)?;
        if pv[1].length_squared() < EFF_ZERO_SQ {
            let t = if left {
                point.curve_parameter + EFF_ZERO_SQRT
            } else {
                point.curve_parameter - EFF_ZERO_SQRT
            };
            pv = self.curve.evaluate(t, 1, left, false)?;
        }
        let middle = self.working_domain.evaluate(0.5, 0.5);
        // `GetULeftEval` and `GetVLeftEval`.
        let sides = [uv.x <= middle.x, uv.y <= middle.y];
        let cells = evaluate(self.surface, uv, 1, 1, sides, true, false, false)?;
        point.surface_point = cells[0];
        point.du = cells[2];
        point.dv = cells[1];
        let mut normal = point.du.cross(point.dv);
        if normal.length_squared() < EFF_ZERO_SQ {
            normal = evaluate_normal(self.surface, uv, [true, true]);
            if normal.length_squared() < EFF_ZERO_SQ {
                return Err(Fault::Degenerate.into());
            }
        }
        normal.unitize()?;
        point.normal = normal;
        point.curve_tangent = pv[1] - normal * pv[1].dot(normal);
        if point.curve_tangent.length_squared() < EFF_ZERO_SQ {
            return Err(Fault::Degenerate.into());
        }
        point.uv = uv;
        let singular = is_singularity(self.surface, uv, EFF_ZERO, ZONE_TOL_3D, false)
            .map_or(true, |s| s.is_some());
        if singular {
            // The drop through a nearby point on a side is not ported.
            return Err(SolverRefusal::NotPorted(INTERSECT));
        }
        point.uv_tangent = drop_vector(point.du, point.dv, point.curve_tangent)?;
        point.drop_distance = line_closest_point(point.surface_point, point.normal, pv[0]);
        let line_point = point.surface_point + point.normal * point.drop_distance;
        point.normal_line_distance = (line_point - pv[0]).length();
        Ok(())
    }

    /// `SmSurfaceTracer::AddStartPoint` (drop points are never singular).
    fn add_start_point(&mut self, point: TracePoint) {
        if self
            .start_points
            .iter()
            .any(|p| point.uv.distance_squared(p.uv) < EFF_ZERO)
        {
            return;
        }
        self.start_points.push(point);
    }

    /// `SmSurfaceDropCurve::FindBoundaryStartPoints` for a curve known to
    /// lie on the surface.
    #[allow(clippy::too_many_lines)]
    fn find_boundary_start_points(&mut self) -> Result<(), SolverRefusal> {
        let target = Some(self.tolerance);
        self.breaks.clear();
        add_break(&mut self.breaks, self.interval.min);
        add_break(&mut self.breaks, self.interval.max);
        let cache = surface_cache(self.surface)?;
        let curve_cache = CurveCache::for_curve(self.curve, 0.5 * ZONE_TOL_3D)?;
        let tolerance = self.tolerance;
        // Side curves (u min, v min, u max, v max) and C0 parameter lines.
        let mut iso_curves: Vec<BSplineCurve> = Vec::new();
        for (axis, value) in [
            (0, self.domain.min.x),
            (1, self.domain.min.y),
            (0, self.domain.max.x),
            (1, self.domain.max.y),
        ] {
            if self.surface.plane_frame().is_some() {
                // A plane's sides are lines; they only matter when they
                // meet the curve's box.
                let natural = self.surface.natural_domain();
                let _ = natural;
                let mut box3 = self.curve.bounding_box(self.interval)?;
                box3.expand_absolute(tolerance);
                let frame = self.surface.plane_frame().ok_or(Fault::Degenerate)?;
                let (mut start, mut end) = (self.domain.min, self.domain.max);
                if axis == 0 {
                    start.x = value;
                    end.x = value;
                } else {
                    start.y = value;
                    end.y = value;
                }
                let line =
                    Line::segment(frame.point_fast(start), frame.point_fast(end), false, None)?;
                let (side_box, _) = line.bounding_boxes(line.natural_interval());
                if !side_box.are_disjoint(box3, 0.0) {
                    return Err(SolverRefusal::NotPorted(INTERSECT));
                }
                continue;
            }
            let interval = if axis == 0 {
                self.domain.v_interval()
            } else {
                self.domain.u_interval()
            };
            let iso = self.surface.iso_curve(axis, value)?;
            let curve = match iso.trim(interval, ZONE_TOL_3D) {
                Ok(trimmed) => trimmed.curve,
                Err(crate::curve_edit::TrimRefusal::Failed(_)) => iso,
                Err(crate::curve_edit::TrimRefusal::ShortSpan) => {
                    return Err(SolverRefusal::NotPorted(INTERSECT));
                }
            };
            iso_curves.push(curve);
        }
        if self.surface.plane_frame().is_none() {
            for axis in 0..2 {
                let conts = if axis == 0 {
                    &cache.u_continuities
                } else {
                    &cache.v_continuities
                };
                let knots = surface_knots(self.surface, axis);
                for i in 1..conts.len().saturating_sub(1) {
                    if conts[i] <= Continuity::C0 {
                        iso_curves.push(self.surface.iso_curve(axis, knots[i])?);
                    }
                }
            }
        }
        for iso in &iso_curves {
            if iso.is_degenerate(EFF_ZERO) {
                let pole = iso.evaluate_point(iso.natural_interval().min)?;
                let dropped = crate::curve_global::drop_point(
                    self.curve,
                    self.interval,
                    pole,
                    None,
                    tolerance,
                    None,
                    INTERSECT,
                )?;
                if dropped.success {
                    add_break(&mut self.breaks, dropped.parameter);
                }
                continue;
            }
            let solutions = global_curve_intersect(
                self.curve,
                self.interval,
                iso,
                iso.natural_interval(),
                tolerance,
                false,
            )?;
            for solution in solutions.as_slice() {
                let mut skip = false;
                if solution.kind == SolutionType::RangeOfValues {
                    if solution.start.value > tolerance / 100.0 {
                        skip = true;
                    } else if (solution.end.parameters[0] - self.interval.min).abs()
                        > EFF_ZERO_SQRT / 100.0
                        && (solution.end.parameters[0] - self.interval.max).abs()
                            > EFF_ZERO_SQRT / 100.0
                    {
                        add_break(&mut self.breaks, solution.end.parameters[0]);
                    }
                }
                let t = solution.start.parameters[0];
                if (t - self.interval.min).abs() < EFF_ZERO_SQRT / 100.0
                    || (t - self.interval.max).abs() < EFF_ZERO_SQRT / 100.0
                {
                    skip = true;
                }
                if !skip {
                    let a = self.curve.evaluate(t, 1, false, false)?;
                    let b = iso.evaluate(solution.start.parameters[1], 1, false, false)?;
                    if a[1].is_parallel_to(b[1], 1.0e-4).0 {
                        skip = true;
                    }
                }
                if !skip {
                    add_break(&mut self.breaks, t);
                }
            }
        }
        let curve_knots: Vec<f64> = unique_knots(self.curve.knots(), None)
            .into_iter()
            .map(|(k, _)| k)
            .collect();
        let conts = &curve_cache.continuities;
        for j in 1..conts.len().saturating_sub(1) {
            if conts[j] <= Continuity::C0 && self.interval.contains_value(curve_knots[j], 0.0) {
                add_break(&mut self.breaks, curve_knots[j]);
            }
        }
        let mut ii = 1;
        while ii < self.breaks.len() {
            let t0 = snap_curve_to_knots(self.curve, self.breaks[ii - 1], EFF_ZERO_SQRT / 100.0);
            let t1 = snap_curve_to_knots(self.curve, self.breaks[ii], EFF_ZERO_SQRT / 100.0);
            self.breaks[ii - 1] = t0;
            self.breaks[ii] = t1;
            if self.breaks[ii - 1] == self.breaks[ii] {
                ii += 1;
                continue;
            }
            let span = Extent1 { min: t0, max: t1 };
            let middle = self.curve.evaluate_point(span.evaluate(0.5))?;
            let found = crate::surface_global::global_point_solve_to(
                self.surface,
                self.domain,
                SolverOperation::Normalize,
                middle,
                tolerance,
                target,
            )?;
            if found.is_empty() {
                ii += 1;
                continue;
            }
            self.working_interval = span;
            let start = self.curve.evaluate_point(span.min)?;
            let mut solutions = crate::surface_global::global_point_solve_to(
                self.surface,
                self.domain,
                INTERSECT,
                start,
                tolerance,
                target,
            )?;
            if solutions.len() == 2
                && (crate::surface_solve::is_closed(self.surface, self.domain, 0, None)
                    || crate::surface_solve::is_closed(self.surface, self.domain, 1, None))
            {
                // A start on the seam of a closed surface.
                return Err(SolverRefusal::NotPorted(INTERSECT));
            }
            if solutions.len() > 1 {
                let near = self.curve.evaluate_point(span.evaluate(0.001))?;
                solutions = crate::surface_global::global_point_solve_to(
                    self.surface,
                    self.domain,
                    INTERSECT,
                    near,
                    tolerance,
                    target,
                )?;
                let (du, dv) = (self.domain.u_interval(), self.domain.v_interval());
                for solution in solutions.as_mut_slice() {
                    let ut = (solution.start.parameters[0] - du.min) / du.length();
                    let vt = (solution.start.parameters[1] - dv.min) / dv.length();
                    if ut < 0.01 {
                        solution.start.parameters[0] = du.min;
                    } else if ut > 0.99 {
                        solution.start.parameters[0] = du.max;
                    }
                    if vt < 0.01 {
                        solution.start.parameters[1] = dv.min;
                    } else if vt > 0.99 {
                        solution.start.parameters[1] = dv.max;
                    }
                }
            }
            let mut jj = 0;
            let mut last_uv = Vector2::new(0.0, 0.0);
            while jj < solutions.len() {
                let solution = solutions.as_slice()[jj];
                let mut uv =
                    Vector2::new(solution.start.parameters[0], solution.start.parameters[1]);
                if jj == 0 {
                    last_uv = uv;
                } else if last_uv.distance(uv) < EFF_ZERO_SQRT * 100.0 {
                    jj += 1;
                    continue;
                }
                if is_singularity(self.surface, uv, EFF_ZERO, ZONE_TOL_3D, false)
                    .map_or(true, |s| s.is_some())
                {
                    // The start dropped on a singularity.
                    return Err(SolverRefusal::NotPorted(INTERSECT));
                }
                uv = snap_surface_to_knots(self.surface, uv, EFF_ZERO_SQRT / 100.0);
                let iso =
                    drop_iso_curve(self.surface, self.domain, self.curve, span, uv, tolerance);
                let iso = match iso {
                    Ok(iso) => iso,
                    Err(SolverRefusal::NotPorted(what)) => {
                        return Err(SolverRefusal::NotPorted(what));
                    }
                    // SMLib reports the failure and moves on.
                    Err(SolverRefusal::Failed(_)) => {
                        jj += 1;
                        continue;
                    }
                };
                if let Some((mut line, mut drop, mut deviation)) = iso {
                    let mut kk = ii + 1;
                    while kk < self.breaks.len() {
                        let longer = Extent1 {
                            min: span.min,
                            max: self.breaks[kk],
                        };
                        kk += 1;
                        match drop_iso_curve(
                            self.surface,
                            self.domain,
                            self.curve,
                            longer,
                            uv,
                            tolerance,
                        )? {
                            Some((line2, drop2, deviation2)) => {
                                line = line2;
                                drop = drop2;
                                deviation = deviation2;
                                ii += 1;
                            }
                            None => break,
                        }
                    }
                    if solutions.len() > jj + 1 {
                        let iso_interval = line.natural_interval();
                        let a = line.evaluate_point(iso_interval.min)?;
                        let b = line.evaluate_point(iso_interval.max)?;
                        let u_range = Extent1::new(a.x.min(b.x), a.x.max(b.x));
                        let v_range = Extent1::new(a.y.min(b.y), a.y.max(b.y));
                        let mut kk = jj + 1;
                        while kk < solutions.len() {
                            let next = solutions.as_slice()[kk].start.parameters;
                            if u_range.contains_value(next[0], 0.0)
                                && v_range.contains_value(next[1], 0.0)
                            {
                                solutions.remove_at(kk, 1)?;
                            } else {
                                kk += 1;
                            }
                        }
                    }
                    self.curves.push(UvCurve::Spline(line));
                    self.drops.push(drop);
                    self.deviations.push(deviation);
                    self.max_distance = max(self.max_distance, drop);
                    jj += 1;
                    continue;
                }
                let near = self.curve.evaluate_point(span.evaluate(0.001))?;
                let near_solutions = crate::surface_global::global_point_solve_to(
                    self.surface,
                    self.domain,
                    SolverOperation::Normalize,
                    near,
                    EFF_ZERO,
                    target,
                )?;
                if let Some(first) = near_solutions.as_slice().first() {
                    let near_uv =
                        Vector2::new(first.start.parameters[0], first.start.parameters[1]);
                    self.working_domain = Extent2 {
                        min: Vector2::new(near_uv.x.min(uv.x), near_uv.y.min(uv.y)),
                        max: Vector2::new(near_uv.x.max(uv.x), near_uv.y.max(uv.y)),
                    };
                }
                if near_solutions.is_empty()
                    || (self.working_domain.x_length() < EFF_ZERO
                        && self.working_domain.y_length() < EFF_ZERO)
                {
                    jj += 1;
                    continue;
                }
                let mut point = TracePoint {
                    curve_parameter: span.min,
                    ..TracePoint::default()
                };
                self.compute_point_values(uv, &mut point)?;
                self.add_start_point(point);
                jj += 1;
            }
            ii += 1;
        }
        Ok(())
    }

    /// `SmSurfaceDropCurve::TraceCurve2`: drops the curve from a start point
    /// to the next break over the continuous part of the domain.
    fn trace_curve2(
        &mut self,
        start: TracePoint,
        cache: &SurfaceCache,
    ) -> Result<(), SolverRefusal> {
        for i in 0..self.breaks.len() {
            let (min_t, max_t) = (start.curve_parameter, self.breaks[i]);
            if max_t < min_t + EFF_ZERO {
                continue;
            }
            let working = domain_with_continuity(
                self.surface,
                cache,
                Continuity::G1,
                start.uv,
                start.uv_tangent,
            );
            if self.domain.are_disjoint(working, 0.0) {
                return Ok(());
            }
            self.working_domain = self.domain.intersect(working)?;
            self.working_interval = Extent1 {
                min: min_t,
                max: max_t,
            };
            // A failed drop falls back on the tracer, which is not ported.
            let dropped = drop_curve(
                self.surface,
                self.working_domain,
                self.curve,
                self.working_interval,
                self.tolerance,
                false,
            )
            .map_err(|_| SolverRefusal::NotPorted(INTERSECT))?;
            if dropped.curves.len() == 1 {
                let curve = dropped.curves.into_iter().next().ok_or(Fault::Degenerate)?;
                self.curves.push(curve);
                self.drops.push(dropped.max_drop);
                self.deviations.push(dropped.deviation);
            } else {
                // `SmSurfaceTracer::TraceCurve` and `FlushCurve`.
                return Err(SolverRefusal::NotPorted(INTERSECT));
            }
            return Ok(());
        }
        Ok(())
    }

    /// `SmSurfaceDropCurve::DoTrace`.
    fn do_trace(&mut self) -> Result<(), SolverRefusal> {
        self.find_boundary_start_points()?;
        let cache = surface_cache(self.surface)?;
        let starts = self.start_points.clone();
        for start in starts {
            self.trace_curve2(start, &cache)?;
        }
        // `FindInteriorCurves`.
        if self.curves.is_empty() {
            let dropped = drop_curve(
                self.surface,
                self.domain,
                self.curve,
                self.interval,
                self.tolerance,
                true,
            )?;
            for curve in dropped.curves {
                self.curves.push(curve);
                self.drops.push(dropped.max_drop);
                self.deviations.push(dropped.deviation);
            }
        }
        Ok(())
    }
}

/// `SmSurface::DropAndTrimCurve` for a curve known to lie on the surface:
/// the parameter curves of every piece of `curve` over `interval` that lies
/// over `domain`, within `tolerance` (or, for a zero `tolerance`, a
/// thousandth of the surface's size). Without `keep_all`, pieces farther
/// than the tolerance from the surface are dropped.
pub fn drop_and_trim_curve(
    surface: &BSplineSurface,
    domain: Extent2,
    curve: &BSplineCurve,
    interval: Extent1,
    tolerance: f64,
    keep_all: bool,
    on_surface: bool,
) -> Result<DroppedCurves, SolverRefusal> {
    if !on_surface {
        return Err(SolverRefusal::NotPorted(INTERSECT));
    }
    let tolerance = if tolerance == 0.0 {
        let p1 = evaluate_point(surface, domain.min)?;
        let middle = evaluate_point(surface, domain.evaluate(0.5, 0.5))?;
        let p2 = evaluate_point(surface, domain.max)?;
        (p1.distance(middle) + middle.distance(p2)) / 1000.0
    } else {
        tolerance
    };
    let mut tracer = DropTracer {
        surface,
        domain,
        curve,
        interval,
        tolerance,
        working_interval: Extent1::default(),
        working_domain: Extent2::default(),
        max_distance: 0.0,
        breaks: Vec::new(),
        start_points: Vec::new(),
        curves: Vec::new(),
        drops: Vec::new(),
        deviations: Vec::new(),
    };
    tracer.do_trace()?;
    let mut result = DroppedCurves::default();
    for ((curve, drop), deviation) in tracer
        .curves
        .into_iter()
        .zip(tracer.drops)
        .zip(tracer.deviations)
    {
        if !keep_all && drop > tolerance {
            continue;
        }
        result.curves.push(curve);
        result.max_drop = max(result.max_drop, drop);
        result.max_deviation = max(result.max_deviation, deviation);
    }
    Ok(result)
}

/// The range answers SMLib builds from a coincident curve: `DropAndTrimCurve`
/// (keeping all pieces when `keep_all`), and when that drops nothing,
/// `DropCurve` without keeping far pieces (its failure counting as no
/// curve when `drop_failure_is_empty`); one range per parameter curve, both
/// ends valued at the largest drop distance, naming no objects.
pub(crate) fn coincident_ranges(
    surface: &BSplineSurface,
    domain: Extent2,
    curve: &BSplineCurve,
    interval: Extent1,
    tolerance: f64,
    keep_all: bool,
    drop_failure_is_empty: bool,
) -> Result<Vec<crate::solution::Solution>, SolverRefusal> {
    let dropped = drop_and_trim_curve(surface, domain, curve, interval, tolerance, keep_all, true)?;
    let (curves, distance) = if dropped.curves.is_empty() {
        match drop_curve(surface, domain, curve, interval, tolerance, false) {
            Ok(drop) => (drop.curves, drop.max_drop),
            Err(SolverRefusal::Failed(_)) if drop_failure_is_empty => {
                (Vec::new(), dropped.max_drop)
            }
            Err(refusal) => return Err(refusal),
        }
    } else {
        (dropped.curves, dropped.max_drop)
    };
    let mut ranges = Vec::with_capacity(curves.len());
    for uv_curve in &curves {
        let span = uv_curve.natural_interval();
        let start = uv_curve.evaluate_point(span.min)?;
        let end = uv_curve.evaluate_point(span.max)?;
        let mut solution = crate::solution::Solution {
            kind: SolutionType::RangeOfValues,
            variables: 3,
            ..crate::solution::Solution::default()
        };
        solution.start.value = distance;
        solution.start.parameters[..3].copy_from_slice(&[span.min, start.x, start.y]);
        solution.end.value = distance;
        solution.end.parameters[..3].copy_from_slice(&[span.max, end.x, end.y]);
        ranges.push(solution);
    }
    Ok(ranges)
}

#[cfg(test)]
mod tests {
    use super::iso_line_is_proven;
    use crate::extent::Extent1;
    use crate::spline::{BSplineCurve, BSplineSurface, ControlPoint};
    use crate::vector::Vector2;

    fn point(x: f64, y: f64) -> ControlPoint {
        ControlPoint {
            coordinates: [x, y, 0.0, 1.0],
        }
    }

    /// The #458 review's witness: the surface's parameter line and the
    /// curve share their control points but differ in the middle knot by
    /// 2^-34, which moves the curve by about 1.2e-4 at t = 0.5. Paired
    /// control points prove nothing across different bases, so the line
    /// must not be proven at a tolerance of 1e-6.
    #[test]
    fn close_but_different_knots_do_not_prove_a_parameter_line() {
        let m = 2f64.powi(20);
        let surface = BSplineSurface::new(
            [1, 1],
            [vec![0.0, 0.0, 0.5, 1.0, 1.0], vec![-1.0, -1.0, 1.0, 1.0]],
            [3, 2],
            vec![
                point(0.0, -1.0),
                point(0.0, 1.0),
                point(m, -1.0),
                point(m, 1.0),
                point(0.0, -1.0),
                point(0.0, 1.0),
            ],
        )
        .expect("surface");
        let delta = 2f64.powi(-34);
        let curve = BSplineCurve::new(
            1,
            vec![0.0, 0.0, 0.5 + delta, 1.0, 1.0],
            vec![point(0.0, 0.0), point(m, 0.0), point(0.0, 0.0)],
        )
        .expect("curve");
        let proven = iso_line_is_proven(
            &surface,
            &curve,
            Extent1::new(0.0, 1.0),
            Vector2::new(0.0, 0.0),
            Vector2::new(1.0, 0.0),
            1.0e-6,
        )
        .expect("answered");
        assert!(!proven);
    }
}
