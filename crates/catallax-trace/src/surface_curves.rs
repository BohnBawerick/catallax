// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline surfaces and planes of
// `spline.rs` without trimming faces; statuses become `Fault`, the output
// arrays of `DoIntersection` one list of [`IntersectionCurve`]s. A surface
// drop that fails inside `IsPointOnCurve` is a typed failure (SMLib reads
// its unset outputs).
// Translated from SMLib SmSurfaceIntersector.cpp
// (SmSurfaceIntersector::DoIntersection, IsPointOnCurve, FlushCurve,
// ReverseCurveDirection, sm_FlushCrv, sm_IsPointUnique),
// SmAdvSurfaceIntersector.cpp (FindInteriorCurves) and SmBSplineCurve.cpp
// (SmBSplineCurve::CreateDegenerateCurve); rust-brep has no port.
//! The intersection curves of step 8: the 3D curve and the two parameter
//! curves of each meeting, built from the traced points as piecewise cubic
//! Bézier spans, or a degenerate curve where the surfaces touch.

use crate::Fault;
use crate::curve_global;
use crate::extent::PointType;
use crate::global_solver::{SolverOperation, SolverRefusal};
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ, PI};
use crate::spline::BSplineCurve;
use crate::surface_drop;
use crate::surface_eval::{evaluate_normal, evaluate_point, first_derivatives};
use crate::surface_intersector::{BoundaryPlaneTest, Intersector, IntersectorTrees, Overlap};
use crate::surface_surface::local_surface_intersect;
use crate::surface_trace::{PointKind, Tracer, TsectPoint};
use crate::vector::{Vector2, Vector3};

/// `SmTsectCurveType`, in SMLib's order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CurveKind {
    /// `SM_TC_UNKNOWN`.
    Unknown,
    /// `SM_TC_TOUCHING`: a single point where the surfaces touch.
    Touching,
    /// `SM_TC_CROSSING`: the normals are not parallel.
    Crossing,
    /// `SM_TC_TANGENT`: the surfaces touch along the curve.
    Tangent,
    /// `SM_TC_COINCIDENT`.
    Coincident,
    /// `SM_TC_NEAR_TANGENT`: a crossing whose normals stay within five
    /// degrees.
    NearTangent,
    /// `SM_TC_REGION_BOUNDARY`.
    RegionBoundary,
}

/// One answer of `DoIntersection`: the 3D curve, the parameter curve on
/// each surface (planar curves with zero `z`), the kind and the largest
/// distance from the traced points to the surfaces.
#[derive(Clone, Debug)]
pub struct IntersectionCurve {
    /// The 3D curve.
    pub curve: BSplineCurve,
    /// The parameter curves on the first and second surface.
    pub uv_curves: [BSplineCurve; 2],
    /// The kind of meeting.
    pub kind: CurveKind,
    /// The largest deviation.
    pub deviation: f64,
}

/// `SmBSplineCurve::CreateDegenerateCurve`: a degree-one curve over
/// `[0, 1]` with both control points at `point`.
pub(crate) fn degenerate_curve(dimension: usize, point: Vector3) -> Result<BSplineCurve, Fault> {
    BSplineCurve::create_canonical(
        dimension,
        1,
        &[point, point],
        &[2, 2],
        &[0.0, 1.0],
        None,
        None,
    )
}

fn flat(uv: Vector2) -> Vector3 {
    Vector3::new(uv.x, uv.y, 0.0)
}

impl IntersectionCurve {
    /// A touching curve: degenerate curves at `point` and at the
    /// parameters `uv` on each surface.
    pub(crate) fn touching(
        point: Vector3,
        uv: [Vector2; 2],
        deviation: f64,
    ) -> Result<Self, Fault> {
        Ok(Self {
            curve: degenerate_curve(3, point)?,
            uv_curves: [
                degenerate_curve(2, flat(uv[0]))?,
                degenerate_curve(2, flat(uv[1]))?,
            ],
            kind: CurveKind::Touching,
            deviation,
        })
    }
}

/// `sm_FlushCrv`: the degree-three curve with the control points of the
/// spans and each distinct knot three times inside (four at the ends).
fn flush_spans(dimension: usize, points: &[Vector3], knots: &[f64]) -> Result<BSplineCurve, Fault> {
    let mut multiplicities = vec![3; knots.len()];
    multiplicities[0] = 4;
    let last = multiplicities.len() - 1;
    multiplicities[last] = 4;
    BSplineCurve::create_canonical(dimension, 3, points, &multiplicities, knots, None, None)
}

/// `sm_IsPointUnique`: no remaining start point and no end of a curve lies
/// within `tolerance` of `point`.
fn is_point_unique(
    point: Vector3,
    starts: &[TsectPoint],
    curves: &[IntersectionCurve],
    tolerance: f64,
) -> bool {
    if starts
        .iter()
        .any(|s| point.distance(s.curve_point) < tolerance)
    {
        return false;
    }
    for entry in curves {
        let natural = entry.curve.natural_interval();
        // SMLib ignores the evaluation status and keeps the point unset.
        for t in [natural.min, natural.max] {
            if let Ok(end) = entry.curve.evaluate_point(t)
                && point.distance(end) < tolerance
            {
                return false;
            }
        }
    }
    true
}

impl Intersector<'_> {
    /// `IsPointOnCurve`: `point` lies within four tolerances of a degenerate
    /// curve, within the tolerance of another, or within ten tolerances
    /// where the surfaces still meet halfway between the point and the
    /// curve. The curve at `skip` is left out. Like SMLib, a curve that
    /// fails to evaluate or drop counts as holding the point (its status
    /// becomes a true result).
    pub fn is_point_on_curve(
        &self,
        point: Vector3,
        curves: &[IntersectionCurve],
        skip: Option<usize>,
    ) -> Result<bool, SolverRefusal> {
        let tolerance = self.approx_tolerance;
        for (index, entry) in curves.iter().enumerate() {
            if skip == Some(index) {
                continue;
            }
            let curve = &entry.curve;
            let natural = curve.natural_interval();
            if curve.is_degenerate(EFF_ZERO) {
                let Ok(start) = curve.evaluate_point(natural.min) else {
                    return Ok(true);
                };
                if start.distance(point) < tolerance * 4.0 {
                    return Ok(true);
                }
                continue;
            }
            let dropped = match curve_global::drop_point(
                curve,
                natural,
                point,
                None,
                tolerance * 10.0,
                None,
                SolverOperation::Intersect,
            ) {
                Ok(dropped) => dropped,
                Err(SolverRefusal::Failed(_)) => return Ok(true),
                Err(refusal) => return Err(refusal),
            };
            if dropped.success && dropped.distance < tolerance {
                return Ok(true);
            }
            if dropped.success && dropped.distance < tolerance * 10.0 {
                let domains = [0, 1].map(|s| self.surfaces[s].natural_domain());
                let drop = |s: usize, target: Vector3, guess: Option<Vector2>| {
                    surface_drop::drop_point(
                        self.surfaces[s],
                        domains[s],
                        target,
                        guess,
                        SolverOperation::Minimize,
                    )
                };
                let at_point = [drop(0, point, None)?, drop(1, point, None)?];
                // SMLib ignores this evaluation's status.
                let on_curve = curve.evaluate_point(dropped.parameter).unwrap_or(point);
                let at_curve = [
                    drop(0, on_curve, Some(at_point[0].uv))?,
                    drop(1, on_curve, Some(at_point[1].uv))?,
                ];
                let middle = [0, 1].map(|s| {
                    let uv = (at_point[s].uv + at_curve[s].uv) / 2.0;
                    evaluate_point(self.surfaces[s], uv).unwrap_or_default()
                });
                let gap = middle[0].distance(middle[1]);
                if gap < (at_point[0].gap + at_point[1].gap) + (at_curve[0].gap + at_curve[1].gap) {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }
}

impl Tracer {
    /// `ReverseCurveDirection`: the traced points in reverse order with
    /// their tangents negated and their parameters counted from the new
    /// start.
    pub fn reverse_curve_direction(&mut self) -> Result<(), Fault> {
        if self.curve_points.is_empty() {
            return Err(Fault::Degenerate);
        }
        let old: Vec<f64> = self
            .curve_points
            .iter()
            .map(|p| p.curve_parameter)
            .collect();
        self.curve_points.reverse();
        let count = old.len();
        let mut parameter = 0.0;
        for (k, point) in self.curve_points.iter_mut().enumerate() {
            // `i` counts down from the old end, as in SMLib.
            let i = count - k;
            point.curve_tangent = -point.curve_tangent;
            point.uv_deriv = [-point.uv_deriv[0], -point.uv_deriv[1]];
            let delta = if i > 1 { old[i - 1] - old[i - 2] } else { 0.0 };
            point.curve_parameter = parameter;
            parameter += delta;
        }
        Ok(())
    }

    /// `FlushCurve`: turns the traced points into an [`IntersectionCurve`]
    /// and empties them. One point becomes a touching curve when it is a
    /// tangent point or lies away from the start points and the ends of the
    /// curves; points all within the tolerance of their neighbours give
    /// nothing. Otherwise each span is a cubic Bézier from the end points
    /// and tangents, a last span shorter than a tenth of the one before
    /// joining it.
    pub fn flush_curve(
        &mut self,
        intersector: &Intersector<'_>,
        curves: &mut Vec<IntersectionCurve>,
    ) -> Result<(), SolverRefusal> {
        let points = std::mem::take(&mut self.curve_points);
        let tolerance = intersector.approx_tolerance;
        let Some(first) = points.first() else {
            return Ok(());
        };
        if points.len() == 1
            && (first.kind == PointKind::TangentPoint
                || is_point_unique(first.curve_point, &self.start_points, curves, tolerance))
        {
            curves.push(IntersectionCurve::touching(
                first.curve_point,
                first.uv,
                first.position(0).distance(first.position(1)),
            )?);
            return Ok(());
        }
        let too_short = points.windows(2).all(|pair| {
            pair[0]
                .curve_point
                .closer_than(tolerance, pair[1].curve_point)
        });
        if too_short {
            return Ok(());
        }
        let mut deviation = first.deviation;
        let mut angle = first.tangent_plane_angle;
        let mut kind = if first.kind == PointKind::TangentCurve {
            CurveKind::Tangent
        } else {
            CurveKind::Crossing
        };
        let mut knots = vec![first.curve_parameter];
        let mut spatial = vec![first.curve_point];
        let mut planar = [vec![flat(first.uv[0])], vec![flat(first.uv[1])]];
        let mut i = 1;
        while i < points.len() {
            let previous = &points[i - 1];
            let mut current = &points[i];
            if i == points.len() - 2 {
                let last = &points[i + 1];
                let delta = current.curve_parameter - previous.curve_parameter;
                let final_delta = last.curve_parameter - current.curve_parameter;
                if final_delta / delta < 0.1 {
                    i += 1;
                    current = last;
                }
            }
            if current.tangent_plane_angle > angle {
                angle = current.tangent_plane_angle;
            }
            if current.deviation > deviation {
                deviation = current.deviation;
            }
            knots.push(current.curve_parameter);
            let scale = current.curve_parameter - previous.curve_parameter;
            let d1 = previous.curve_tangent * scale;
            let d2 = current.curve_tangent * scale;
            spatial.push(previous.curve_point + d1 / 3.0);
            spatial.push(current.curve_point - d2 / 3.0);
            spatial.push(current.curve_point);
            for (s, curve) in planar.iter_mut().enumerate() {
                let d1 = previous.uv_deriv[s] * scale;
                let d2 = current.uv_deriv[s] * scale;
                curve.push(flat(previous.uv[s] + d1 / 3.0));
                curve.push(flat(current.uv[s] - d2 / 3.0));
                curve.push(flat(current.uv[s]));
            }
            i += 1;
        }
        let curve = flush_spans(3, &spatial, &knots)?;
        let uv_curves = [
            flush_spans(2, &planar[0], &knots)?,
            flush_spans(2, &planar[1], &knots)?,
        ];
        if angle < 5.0 * (PI / 180.0) && kind == CurveKind::Crossing && spatial.len() > 4 {
            kind = CurveKind::NearTangent;
        }
        curves.push(IntersectionCurve {
            curve,
            uv_curves,
            kind,
            deviation,
        });
        Ok(())
    }
}

/// Whether a start point's tangent leaves both domains' sides along their
/// interiors (each parameter inside, or its tangent component negligible):
/// then the curve may continue behind the start as well.
fn free_both_ways(start: &TsectPoint) -> bool {
    (0..2).all(|s| {
        (start.param_types[s][0] == PointType::Inside || start.uv_deriv[s].x.abs() <= EFF_ZERO)
            && (start.param_types[s][1] == PointType::Inside
                || start.uv_deriv[s].y.abs() <= EFF_ZERO)
    })
}

/// The angle between the surface normals at a parameter pair, folded into
/// `[0, pi/2]`, in degrees.
fn normal_angle_degrees(intersector: &Intersector<'_>, uv: [Vector2; 2]) -> Result<f64, Fault> {
    let n1 = evaluate_normal(intersector.surfaces[0], uv[0], [true, true]);
    let n2 = evaluate_normal(intersector.surfaces[1], uv[1], [true, true]);
    let mut angle = n1.angle_between(n2)?;
    if angle > PI / 2.0 {
        angle = PI - angle;
    }
    Ok(angle * (180.0 / PI))
}

impl Tracer {
    /// `SmAdvSurfaceIntersector::FindInteriorCurves`: from the closest
    /// midpoints of each leaf pair, meets the surfaces with a hundredth of
    /// the tolerance and traces from each meeting at least two degrees from
    /// tangent and off the curves found so far, both ways when the curve
    /// neither closes nor starts on a side.
    pub fn find_interior_curves(
        &mut self,
        intersector: &Intersector<'_>,
        trees: &IntersectorTrees,
        curves: &mut Vec<IntersectionCurve>,
    ) -> Result<(), SolverRefusal> {
        for guess in intersector.interior_guesses(trees)? {
            if normal_angle_degrees(intersector, guess)? < 5.0 {
                continue;
            }
            let met = local_surface_intersect(
                intersector.surfaces[0],
                intersector.domains[0],
                intersector.surfaces[1],
                intersector.domains[1],
                intersector.approx_tolerance / 100.0,
                guess[0],
                guess[1],
                None,
                None,
                false,
            )?;
            let Some(pair) = met else {
                continue;
            };
            let uvs = [pair.uv, pair.other_uv];
            let surface_point = evaluate_point(intersector.surfaces[0], uvs[0])?;
            if normal_angle_degrees(intersector, uvs)? < 2.0 {
                continue;
            }
            if intersector.is_point_on_curve(surface_point, curves, None)? {
                continue;
            }
            self.trace_direction = 1.0;
            let mut start = TsectPoint::default();
            match intersector.compute_point_values(
                uvs,
                &mut start,
                None,
                self.curve_points.last().map(|p| p.curve_tangent),
                None,
                self.trace_direction,
            ) {
                Ok(()) => {}
                Err(SolverRefusal::Failed(_)) => continue,
                Err(refusal) => return Err(refusal),
            }
            let input_direction = self.trace_direction;
            self.has_start_point = true;
            self.trace_curve(intersector, &start, curves)?;
            if !self.curve_is_closed
                && !start.is_bounded()
                && input_direction == self.trace_direction
            {
                self.reverse_curve_direction()?;
                self.trace_direction = -self.trace_direction;
                self.trace_curve(intersector, &start, curves)?;
                self.reverse_curve_direction()?;
            }
            self.flush_curve(intersector, curves)?;
            self.has_start_point = false;
        }
        Ok(())
    }
}

/// `SmVector3d::CCWAngleBetween` about `axis`: the angle in `[-pi, pi]`
/// from `start` to `end` projected across `axis`, zero when either
/// projection vanishes.
fn ccw_angle_between(axis: Vector3, start: Vector3, end: Vector3) -> Result<f64, Fault> {
    let start = axis.cross(start).cross(axis);
    if start.length_squared() < 100.0 * EFF_ZERO_SQ {
        return Ok(0.0);
    }
    let end = axis.cross(end).cross(axis);
    if end.length_squared() < 100.0 * EFF_ZERO_SQ {
        return Ok(0.0);
    }
    let angle = start.angle_between(end)?;
    Ok(if start.cross(end).dot(axis) < 0.0 {
        -angle
    } else {
        angle
    })
}

/// `SmVector3d::IsInsideSector`: whether `vector` lies in the sector swept
/// from `bound1` (with `binormal1` pointing inside) to `bound2`; vectors all
/// along one direction count as inside.
fn is_inside_sector(vector: Vector3, bound1: Vector3, binormal1: Vector3, bound2: Vector3) -> bool {
    let mut normal = bound1.cross(binormal1);
    // SMLib logs and goes on when these fail; they cannot for its callers.
    let _ = normal.unitize();
    let mut sector = ccw_angle_between(normal, bound1, bound2).unwrap_or(0.0);
    let mut this = ccw_angle_between(normal, bound1, vector).unwrap_or(0.0);
    if sector < 0.0 {
        sector += 2.0 * PI;
    }
    if this < -EFF_ZERO {
        this += 2.0 * PI;
    }
    if (this - sector).abs() >= 2.0 * PI - EFF_ZERO {
        return true;
    }
    this > -EFF_ZERO && this < sector + EFF_ZERO
}

impl Intersector<'_> {
    /// `GrazingBoxProcess`: where the meeting box is thin, looks for a pair
    /// of domain corners, one of each surface inside the meeting domains,
    /// that coincide within the tolerance. When the second surface's corner
    /// directions both point away from the first's corner sector, the
    /// surfaces only touch there: the answer is that touching curve. `None`
    /// asks for the general intersection.
    pub fn grazing_box_process(
        &self,
        overlap: &Overlap,
    ) -> Result<Option<Vec<IntersectionCurve>>, SolverRefusal> {
        let corners = [(0, 0), (0, 1), (1, 0), (1, 1)];
        for (iu1, iv1) in corners {
            let uv1 = self.domains[0].evaluate(f64::from(iu1), f64::from(iv1));
            if !overlap.domains[0].contains_point(uv1, 0.0) {
                continue;
            }
            for (iu2, iv2) in corners {
                let uv2 = self.domains[1].evaluate(f64::from(iu2), f64::from(iv2));
                if !overlap.domains[1].contains_point(uv2, 0.0) {
                    continue;
                }
                let [point1, mut du1, mut dv1] =
                    first_derivatives(self.surfaces[0], uv1, [true, true])?;
                let [point2, mut du2, mut dv2] =
                    first_derivatives(self.surfaces[1], uv2, [true, true])?;
                if iu1 == 1 {
                    du1 = -du1;
                }
                if iv1 == 1 {
                    dv1 = -dv1;
                }
                if iu2 == 1 {
                    du2 = -du2;
                }
                if iv2 == 1 {
                    dv2 = -dv2;
                }
                if point1.distance(point2) >= self.approx_tolerance {
                    continue;
                }
                if [du1, dv1, du2, dv2]
                    .iter()
                    .any(|d| d.length_squared() < EFF_ZERO)
                {
                    continue;
                }
                for d in [&mut du1, &mut dv1, &mut du2, &mut dv2] {
                    d.unitize()?;
                }
                let normal1 = du1.cross(dv1);
                if normal1.length_squared() < EFF_ZERO {
                    continue;
                }
                let binormal = normal1.cross(du1);
                let projected_du2 = normal1.cross(du2).cross(normal1);
                let projected_dv2 = normal1.cross(dv2).cross(normal1);
                if projected_du2.length_squared() < EFF_ZERO_SQ
                    && projected_dv2.length_squared() < EFF_ZERO_SQ
                {
                    return Ok(None);
                }
                let u_inside = projected_du2.length_squared() > EFF_ZERO_SQ
                    && is_inside_sector(projected_du2, du1, binormal, dv1);
                let v_inside = projected_dv2.length_squared() > EFF_ZERO_SQ
                    && is_inside_sector(projected_dv2, du1, binormal, dv1);
                if u_inside || v_inside {
                    return Ok(None);
                }
                let touching =
                    IntersectionCurve::touching(point1, [uv1, uv2], point1.distance(point2))?;
                return Ok(Some(vec![touching]));
            }
        }
        Ok(None)
    }

    /// `SmSurfaceIntersector::DoIntersection` of the advanced intersector
    /// without trimming faces: the meeting curves of the two surfaces over
    /// their domains. The overlap search and the side-plane test bound the
    /// work; start points on the sides are traced first (both ways when the
    /// curve can run behind its start), then curves inside the domains, and
    /// last a touching curve lying on another curve is dropped. Grazing
    /// meetings, self-intersection and the paths the parts refuse are
    /// [`SolverRefusal::NotPorted`].
    pub fn do_intersection(&self) -> Result<Vec<IntersectionCurve>, SolverRefusal> {
        if self.self_intersection {
            return Err(SolverRefusal::NotPorted(SolverOperation::Intersect));
        }
        let trees = self.trees()?;
        let Some(overlap) = self.overlap_box(&trees)? else {
            return Ok(Vec::new());
        };
        if self.grazes(&overlap)
            && let Some(curves) = self.grazing_box_process(&overlap)?
        {
            return Ok(curves);
        }
        if self.boundary_plane_test(overlap.bbox)? == BoundaryPlaneTest::Disjoint {
            return Ok(Vec::new());
        }
        let mut curves = Vec::new();
        let mut tracer = Tracer {
            trace_direction: 1.0,
            ..Tracer::default()
        };
        tracer.find_boundary_start_points(self, overlap.bbox, &mut curves)?;
        while let Some(start) = tracer.start_points.pop() {
            tracer.trace_direction = 1.0;
            let count = curves.len();
            let input_direction = tracer.trace_direction;
            match tracer.trace_curve(self, &start, &curves) {
                Ok(()) => {
                    if !tracer.curve_is_closed
                        && input_direction == tracer.trace_direction
                        && free_both_ways(&start)
                    {
                        tracer.reverse_curve_direction()?;
                        tracer.trace_direction = -tracer.trace_direction;
                        tracer.trace_curve(self, &start, &curves)?;
                        tracer.reverse_curve_direction()?;
                    }
                    tracer.flush_curve(self, &mut curves)?;
                }
                // SMLib keeps the points of a failed trace for the next one.
                Err(SolverRefusal::Failed(_)) => {}
                Err(refusal) => return Err(refusal),
            }
            if count == curves.len() {
                tracer.through_points.push(start);
            }
        }
        tracer.find_interior_curves(self, &trees, &mut curves)?;
        let mut index = 0;
        while index < curves.len() {
            if curves[index].curve.is_degenerate(EFF_ZERO) {
                // SMLib ignores this evaluation's status.
                let point = curves[index].curve.evaluate_point(0.0).unwrap_or_default();
                if self.is_point_on_curve(point, &curves, Some(index))? {
                    curves.remove(index);
                    continue;
                }
            }
            index += 1;
        }
        Ok(curves)
    }
}
