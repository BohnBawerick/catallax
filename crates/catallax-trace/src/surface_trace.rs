// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline surfaces and planes of
// `spline.rs`; statuses become `Fault`, the intersection point is a plain
// struct and the intersector's tracing state the `Intersector`. Paths that
// need code from later parts of step 8 refuse with
// `SolverRefusal::NotPorted`: the singularity point of nearly parallel
// surfaces (`SmAdvSurfaceIntersector::ComputeSingularityPoint`) and the step
// back from a surface singularity (`StepBackFromSingularity`).
// Translated from SMLib SmSurfaceIntersector.cpp (SmTsectPnt,
// SmTsectPnt::ClassifyPointParams, SmSurfaceIntersector::ComputePointValues);
// rust-brep has no port.
// Indexed loops follow SMLib's per-surface loops.
#![allow(clippy::needless_range_loop)]
//! The intersection points of step 8: a pair of surface parameters with the
//! surfaces' derivatives and normals, the meeting curve's point and unit
//! tangent, and the tangent projected into each parameter domain.

use crate::Fault;
use crate::curve_surface::{drop_vector, is_singularity};
use crate::extent::{Extent2, PointType};
use crate::global_solver::{SolverOperation, SolverRefusal};
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, PI, ZONE_TOL_3D, max, min};
use crate::surface_eval::{evaluate, evaluate_normal};
use crate::surface_intersector::Intersector;
use crate::surface_surface::local_surface_solve;
use crate::vector::{Vector2, Vector3};

/// `SmIntersectionPointType`, in SMLib's order (the tracer compares them).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum PointKind {
    /// `SM_IP_UNDEFINED`.
    Undefined,
    /// `SM_IP_UNKNOWN`.
    Unknown,
    /// `SM_IP_CROSSING`: the normals are not parallel.
    Crossing,
    /// `SM_IP_TANGENT_POINT`: the surfaces touch at a point.
    TangentPoint,
    /// `SM_IP_TANGENT_CURVE`: the point lies on a tangent curve.
    TangentCurve,
    /// `SM_IP_SINGULARITY`: four or six branches meet.
    Singularity,
    /// `SM_IP_COINCIDENCE`: the second fundamental forms match.
    Coincidence,
    /// `SM_IP_TOUCHING`.
    Touching,
}

/// `SmTsectPnt`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TsectPoint {
    /// `m_ePointType`.
    pub kind: PointKind,
    /// `m_eUVParamType`: each surface's u and v against its domain.
    pub param_types: [[PointType; 2]; 2],
    /// `UVPos`.
    pub uv: [Vector2; 2],
    /// `UVDeriv`: the curve tangent in each parameter domain.
    pub uv_deriv: [Vector2; 2],
    /// `m_vSurfacePV`: SMLib's evaluation matrices, cell `k * 3 + l` for
    /// the derivative `S^(k,l)`, with the unit normal in cell 8.
    pub surface: [[Vector3; 9]; 2],
    /// `CrvPos`.
    pub curve_point: Vector3,
    /// `CrvDeriv`.
    pub curve_tangent: Vector3,
    /// `m_dTangentPlaneAngleRad`.
    pub tangent_plane_angle: f64,
    /// `m_dDeviation`.
    pub deviation: f64,
    /// `m_dCurveParameter`: arc length along the traced curve.
    pub curve_parameter: f64,
}

impl Default for TsectPoint {
    /// `SmTsectPnt()`: unknown, unclassified, zeroed.
    fn default() -> Self {
        let zero = Vector3::new(0.0, 0.0, 0.0);
        Self {
            kind: PointKind::Unknown,
            param_types: [[PointType::Outside; 2]; 2],
            uv: [Vector2::new(0.0, 0.0); 2],
            uv_deriv: [Vector2::new(0.0, 0.0); 2],
            surface: [[zero; 9]; 2],
            curve_point: zero,
            curve_tangent: zero,
            tangent_plane_angle: 0.0,
            deviation: 0.0,
            curve_parameter: 0.0,
        }
    }
}

impl TsectPoint {
    /// `SrfPos`.
    pub fn position(&self, s: usize) -> Vector3 {
        self.surface[s][0]
    }
    /// `SrfDu`.
    pub fn du(&self, s: usize) -> Vector3 {
        self.surface[s][3]
    }
    /// `SrfDv`.
    pub fn dv(&self, s: usize) -> Vector3 {
        self.surface[s][1]
    }
    /// `SrfNorm`.
    pub fn normal(&self, s: usize) -> Vector3 {
        self.surface[s][8]
    }

    /// `ClassifyPointParams`: each parameter against its domain, within the
    /// 3D tolerance mapped by the first derivative (at most a tenth of the
    /// interval).
    fn classify(&mut self, uvs: [Vector2; 2], domains: [Extent2; 2], tolerance: f64) {
        for s in 0..2 {
            let (u, v) = (domains[s].u_interval(), domains[s].v_interval());
            let mut u_size = self.du(s).length();
            if u_size < EFF_ZERO {
                u_size = EFF_ZERO;
            }
            let mut v_size = self.dv(s).length();
            if v_size < EFF_ZERO {
                v_size = EFF_ZERO;
            }
            let mut tol_u = min(u.length() / 10.0, tolerance / max(0.99, u_size));
            let mut tol_v = min(v.length() / 10.0, tolerance / max(0.99, v_size));
            if tol_u < EFF_ZERO {
                tol_u = EFF_ZERO;
            }
            if tol_v < EFF_ZERO {
                tol_v = EFF_ZERO;
            }
            self.param_types[s] = [
                u.classify_point(uvs[s].x, tol_u),
                v.classify_point(uvs[s].y, tol_v),
            ];
        }
    }

    /// Whether the curve tangent leaves some side of a domain it lies on
    /// (NaN derivatives count as leaving, as in SMLib's negated tests).
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    fn points_outward(&self) -> bool {
        (0..2).any(|s| {
            let [tu, tv] = self.param_types[s];
            let d = self.uv_deriv[s];
            (tu == PointType::Start && !(d.x >= -EFF_ZERO))
                || (tu == PointType::End && !(d.x <= EFF_ZERO))
                || (tv == PointType::Start && !(d.y >= -EFF_ZERO))
                || (tv == PointType::End && !(d.y <= EFF_ZERO))
        })
    }
}

fn singular(intersector: &Intersector<'_>, s: usize, uv: Vector2) -> bool {
    // SMLib's IsSingularity reads a failed test as singular.
    is_singularity(intersector.surfaces[s], uv, EFF_ZERO, ZONE_TOL_3D, false)
        .map_or(true, |found| found.is_some())
}

fn angle_between_normals(point: &TsectPoint) -> Result<f64, Fault> {
    let mut angle = point.normal(0).angle_between(point.normal(1))?;
    if angle > PI / 2.0 {
        angle = PI - angle;
    }
    Ok(angle)
}

impl Intersector<'_> {
    /// `ComputePointValues`: evaluates both surfaces at `uvs` and fills the
    /// intersection point: positions, derivatives and unit normals, the
    /// parameter classes, the curve point and unit tangent along the
    /// normals' cross product (turned into the domains, or along the
    /// previous point's direction), its projections into the domains and the
    /// angle between the normals; near parallel normals after `previous` it
    /// refines the point by a normal solve.
    pub fn compute_point_values(
        &self,
        uvs: [Vector2; 2],
        point: &mut TsectPoint,
        previous: Option<&TsectPoint>,
        step_size: Option<&mut f64>,
        trace_direction: f64,
    ) -> Result<(), SolverRefusal> {
        for s in 0..2 {
            point.uv[s] = uvs[s];
            let cells = evaluate(
                self.surfaces[s],
                uvs[s],
                2,
                2,
                [true, true],
                true,
                true,
                true,
            )?;
            for (slot, value) in point.surface[s].iter_mut().zip(&cells) {
                *slot = *value;
            }
            let mut normal = point.du(s).cross(point.dv(s));
            if normal.length_squared() < EFF_ZERO_SQ {
                normal = evaluate_normal(self.surfaces[s], uvs[s], [true, true]);
                if normal.length_squared() < EFF_ZERO_SQ {
                    return Err(Fault::Degenerate.into());
                }
            }
            let _ = normal.unitize();
            point.surface[s][8] = normal;
        }
        point.classify(uvs, self.domains, self.approx_tolerance);
        point.curve_point = (point.position(0) + point.position(1)).divided(2.0);
        point.deviation = point.curve_point.distance(point.position(0));
        if point.kind != PointKind::Singularity {
            point.kind = PointKind::Crossing;
            point.curve_tangent = point.normal(0).cross(point.normal(1)) * trace_direction;
            point.tangent_plane_angle = angle_between_normals(point)?;
            match previous {
                None => {
                    for s in 0..2 {
                        point.uv_deriv[s] =
                            drop_vector(point.du(s), point.dv(s), point.curve_tangent)?;
                    }
                    if point.points_outward() {
                        point.uv_deriv = [-point.uv_deriv[0], -point.uv_deriv[1]];
                        point.curve_tangent = -point.curve_tangent;
                    }
                }
                Some(previous) => {
                    let same = if previous.kind == PointKind::Crossing {
                        let last = 0.0
                            < previous
                                .curve_tangent
                                .dot(previous.normal(0).cross(previous.normal(1)));
                        let this = 0.0
                            < point
                                .curve_tangent
                                .dot(point.normal(0).cross(point.normal(1)));
                        last == this
                    } else {
                        0.0 < previous.curve_tangent.dot(point.curve_tangent)
                    };
                    if !same {
                        point.curve_tangent = -point.curve_tangent;
                    }
                }
            }
            if point.tangent_plane_angle < self.start_parallel_tolerance
                || point.curve_tangent.length() < self.start_parallel_tolerance
                || previous.is_some_and(|p| p.kind == PointKind::TangentCurve)
            {
                match previous {
                    Some(previous) if singular(self, 0, uvs[0]) || singular(self, 1, uvs[1]) => {
                        point.curve_tangent = point
                            .normal(0)
                            .cross(previous.curve_tangent)
                            .cross(point.normal(0));
                        point.kind = PointKind::TangentCurve;
                    }
                    _ => {
                        // `ComputeSingularityPoint` of the advanced intersector.
                        return Err(SolverRefusal::NotPorted(SolverOperation::Intersect));
                    }
                }
            }
        }
        let tangent_length = point.curve_tangent.length();
        if tangent_length < self.start_parallel_tolerance {
            point.kind = PointKind::TangentPoint;
            return Ok(());
        }
        point.curve_tangent /= tangent_length;
        if let (Some(_), Some(_)) = (previous, step_size) {
            for s in 0..2 {
                if (point.du(s).length_squared() < EFF_ZERO_SQRT
                    || point.dv(s).length_squared() < EFF_ZERO_SQRT)
                    && singular(self, s, uvs[s])
                {
                    // `StepBackFromSingularity`.
                    return Err(SolverRefusal::NotPorted(SolverOperation::Intersect));
                }
            }
        }
        for s in 0..2 {
            point.uv_deriv[s] = drop_vector(point.du(s), point.dv(s), point.curve_tangent)?;
        }
        point.tangent_plane_angle = angle_between_normals(point)?;
        if let Some(previous) = previous
            && point.kind == PointKind::Crossing
            && point.tangent_plane_angle < 3.0 * PI / 180.0
        {
            let mut found = None;
            if previous.kind != PointKind::Unknown && previous.kind != PointKind::Singularity {
                let grow = |s: usize| {
                    let further = point.uv[s] + (point.uv[s] - previous.uv[s]);
                    Extent2::new(previous.uv[s], previous.uv[s]).union(Extent2::new(
                        self.domains[s].clamp_point(further),
                        self.domains[s].clamp_point(further),
                    ))
                };
                let (domain1, domain2) = (grow(0), grow(1));
                found = local_surface_solve(
                    self.surfaces[0],
                    domain1,
                    self.surfaces[1],
                    domain2,
                    SolverOperation::Normalize,
                    domain1.evaluate(0.5, 0.5),
                    domain2.evaluate(0.5, 0.5),
                )?;
            }
            if let Some(pair) = found {
                point.kind = PointKind::Unknown;
                self.compute_point_values(
                    [pair.uv, pair.other_uv],
                    point,
                    None,
                    None,
                    trace_direction,
                )?;
                if previous.curve_tangent.dot(point.curve_tangent) < 0.0 {
                    point.curve_tangent = -point.curve_tangent;
                    for s in 0..2 {
                        point.uv_deriv[s] =
                            drop_vector(point.du(s), point.dv(s), point.curve_tangent)?;
                    }
                }
            } else if point.kind == PointKind::Crossing
                && point.tangent_plane_angle < self.start_parallel_tolerance
            {
                point.curve_tangent = previous.curve_tangent;
                let _ = point.curve_tangent.unitize();
                for s in 0..2 {
                    point.uv_deriv[s] = drop_vector(point.du(s), point.dv(s), point.curve_tangent)?;
                }
                return Ok(());
            }
        }
        Ok(())
    }
}

/// What [`Intersector::compute_next_point`] reports besides the point:
/// SMLib's `rbFoundGoodPoint`, `rdDeviationFound`, `rdAngleFoundRad`,
/// `rbClipped` and `rbBoundaryHit`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StepOutcome {
    /// The step met the angle and distance tolerances.
    pub good: bool,
    /// The largest distance found along the span.
    pub deviation: f64,
    /// Whether the span test wrote `deviation`; SMLib's caller otherwise
    /// keeps its previous value.
    pub deviation_written: bool,
    /// The turn of the curve tangent over the span.
    pub angle: f64,
    /// A domain side shortened the step.
    pub clipped: bool,
    /// The step stopped on a domain side.
    pub boundary_hit: bool,
}

impl Intersector<'_> {
    /// `ComputeStepSize`: the old step, or a hundred tolerances to start.
    pub fn compute_step_size(&self, old_step: f64) -> f64 {
        if crate::scalar::is_zero(old_step) {
            100.0 * self.approx_tolerance
        } else {
            old_step
        }
    }

    /// `RefinePoint`: the meeting near `uvs` (looser for tangent and
    /// singular points); a tangent point is not refined and SMLib then reads
    /// an untouched solution, the parameters `(0, 0)` clamped into each
    /// domain.
    pub fn refine_point(
        &self,
        uvs: [Vector2; 2],
        kind: PointKind,
        plane_normal: Option<Vector3>,
        previous_uv: Option<Vector2>,
        boundary_testing: bool,
    ) -> Result<Option<[Vector2; 2]>, SolverRefusal> {
        let tolerance = match kind {
            PointKind::TangentCurve | PointKind::Singularity => self.approx_tolerance * 10.00,
            PointKind::Crossing => self.approx_tolerance,
            PointKind::TangentPoint => {
                let zero = Vector2::new(0.0, 0.0);
                return Ok(Some([
                    self.domains[0].clamp_point(zero),
                    self.domains[1].clamp_point(zero),
                ]));
            }
            _ => return Ok(None),
        };
        let Some(pair) = crate::surface_surface::local_surface_intersect(
            self.surfaces[0],
            self.domains[0],
            self.surfaces[1],
            self.domains[1],
            tolerance,
            uvs[0],
            uvs[1],
            plane_normal,
            previous_uv,
            boundary_testing,
        )?
        else {
            return Ok(None);
        };
        Ok(Some([
            self.domains[0].clamp_point(pair.uv),
            self.domains[1].clamp_point(pair.other_uv),
        ]))
    }

    /// `TestSpanAccuracy`: the span from `start` to `next` is good when the
    /// tangent turns less than the angle tolerance and, at five points of
    /// the Hermite curves through both ends, the two surfaces and the 3D
    /// curve stay within the tolerance of one another. A failing span
    /// demotes a singular `next` to crossing.
    pub fn test_span_accuracy(
        &self,
        start: &TsectPoint,
        next: &mut TsectPoint,
    ) -> Result<(bool, Option<f64>, f64), SolverRefusal> {
        use crate::hermite::Hermite;
        if next.curve_tangent.length() < self.stepping_parallel_tolerance {
            return Err(Fault::Degenerate.into());
        }
        let step = next.curve_parameter - start.curve_parameter;
        let angle = start.curve_tangent.angle_between(next.curve_tangent)?;
        let demote = |next: &mut TsectPoint| {
            if next.kind == PointKind::Singularity {
                next.kind = PointKind::Crossing;
            }
        };
        if angle > self.angle_tolerance && step > EFF_ZERO * 100.0 {
            demote(next);
            return Ok((false, None, angle));
        }
        let uv_curves = [0, 1].map(|s| {
            Hermite::planar(
                start.uv[s],
                start.uv_deriv[s] * step,
                next.uv[s],
                next.uv_deriv[s] * step,
            )
        });
        let curve = Hermite::new(
            start.curve_point,
            start.curve_tangent * step,
            next.curve_point,
            next.curve_tangent * step,
        );
        let tolerance_sq = self.approx_tolerance * self.approx_tolerance;
        let mut max_sq: f64 = 0.0;
        for i in 0..5 {
            let t = (i as f64 + 1.0) / (5.0 + 1.0);
            let on = |s: usize| -> Result<Vector3, SolverRefusal> {
                let q = uv_curves[s].evaluate_point(t);
                let uv = self.domains[s].clamp_point(Vector2::new(q.x, q.y));
                Ok(crate::surface_eval::evaluate_point(self.surfaces[s], uv)?)
            };
            let (p0, p1) = (on(0)?, on(1)?);
            let c = curve.evaluate_point(t);
            let distance_sq = crate::scalar::max(
                crate::scalar::max(p0.distance_squared(p1), p0.distance_squared(c)),
                p1.distance_squared(c),
            );
            if distance_sq > max_sq {
                max_sq = distance_sq;
            }
            if distance_sq > tolerance_sq && step > EFF_ZERO * 100.0 {
                demote(next);
                return Ok((false, Some(crate::scalar::sqrt(max_sq)), angle));
            }
        }
        Ok((true, Some(crate::scalar::sqrt(max_sq)), angle))
    }

    /// `ComputeNextPoint`: steps `step` along the previous point's tangent
    /// (clipped to the domains), refines the meeting there, computes its
    /// values and tests the span.
    pub fn compute_next_point(
        &self,
        previous: &TsectPoint,
        step: f64,
        next: &mut TsectPoint,
        trace_direction: f64,
    ) -> Result<StepOutcome, SolverRefusal> {
        let mut outcome = StepOutcome::default();
        if !self.domains[0].contains_point(previous.uv[0], EFF_ZERO)
            || !self.domains[1].contains_point(previous.uv[1], EFF_ZERO)
        {
            return Err(Fault::Degenerate.into());
        }
        let unclipped = step;
        let mut step = self.domains[0].clip_line(previous.uv[0], previous.uv_deriv[0], step);
        step = self.domains[1].clip_line(previous.uv[1], previous.uv_deriv[1], step);
        outcome.clipped = step < unclipped - EFF_ZERO * (1.0 + unclipped);
        if crate::scalar::is_zero(step) {
            outcome.boundary_hit = true;
            return Ok(outcome);
        }
        let guess = [0, 1]
            .map(|s| self.domains[s].clamp_point(previous.uv[s] + step * previous.uv_deriv[s]));
        let refine = |boundary: bool| {
            self.refine_point(
                guess,
                previous.kind,
                Some(previous.curve_tangent),
                Some(previous.uv[0]),
                boundary,
            )
        };
        // SMLib treats a failed refinement like a missed one.
        let Some(mut uvs) = refine(true).or_else(|r| match r {
            SolverRefusal::Failed(_) => Ok(None),
            other => Err(other),
        })?
        else {
            return Ok(outcome);
        };
        let scale = |step: &mut f64, uvs: &[Vector2; 2]| {
            let guess_step = (guess[0] - previous.uv[0]).length();
            let actual_step = (uvs[0] - previous.uv[0]).length();
            if guess_step > EFF_ZERO * actual_step {
                *step *= actual_step / guess_step;
            }
            (guess_step, actual_step)
        };
        let (guess_step, actual_step) = scale(&mut step, &uvs);
        if crate::scalar::is_zero(step) || actual_step < guess_step / 100.0 {
            let Some(again) = refine(false).or_else(|r| match r {
                SolverRefusal::Failed(_) => Ok(None),
                other => Err(other),
            })?
            else {
                return Ok(outcome);
            };
            uvs = again;
            scale(&mut step, &uvs);
            if crate::scalar::is_zero(step) {
                outcome.boundary_hit = true;
                return Ok(outcome);
            }
        }
        match self.compute_point_values(uvs, next, Some(previous), Some(&mut step), trace_direction)
        {
            Ok(()) => {}
            // SMLib reads any failure here as a boundary hit.
            Err(SolverRefusal::Failed(_)) => {
                outcome.boundary_hit = true;
                return Ok(outcome);
            }
            Err(refusal) => return Err(refusal),
        }
        if next.kind == PointKind::TangentPoint || next.kind == PointKind::Coincidence {
            return Ok(outcome);
        }
        let along = crate::geometry::line_closest_point(
            Vector3::new(previous.uv[0].x, previous.uv[0].y, 0.0),
            Vector3::new(previous.uv_deriv[0].x, previous.uv_deriv[0].y, 0.0),
            Vector3::new(next.uv[0].x, next.uv[0].y, 0.0),
        );
        if along < 0.0 {
            if next.kind == PointKind::Singularity {
                next.kind = PointKind::Crossing;
            }
            return Ok(outcome);
        }
        next.curve_parameter =
            previous.curve_parameter + next.curve_point.distance(previous.curve_point);
        let (good, deviation, angle) = self.test_span_accuracy(previous, next)?;
        outcome.good = good;
        if let Some(deviation) = deviation {
            outcome.deviation = deviation;
            outcome.deviation_written = true;
        }
        outcome.angle = angle;
        Ok(outcome)
    }
}

impl crate::curve_solve::PointSolveCurve for crate::hermite::Hermite {
    fn evaluate(
        &self,
        t: f64,
        derivatives: usize,
        _from_left: bool,
        nonzero_tangent: bool,
    ) -> Result<Vec<Vector3>, Fault> {
        crate::hermite::Hermite::evaluate(self, t, derivatives, nonzero_tangent)
    }
    fn evaluate_point(&self, t: f64) -> Result<Vector3, Fault> {
        Ok(crate::hermite::Hermite::evaluate_point(self, t))
    }
    /// `SmCurve::IsDegenerate` over the natural interval.
    fn is_degenerate(&self, tolerance: f64) -> bool {
        let interval = self.natural_interval();
        let p1 = crate::hermite::Hermite::evaluate_point(self, interval.min);
        let p2 = crate::hermite::Hermite::evaluate_point(self, interval.max);
        let scaled_zero = EFF_ZERO * (1.0 + p1.max_dimension());
        let tolerance_sq = max(tolerance * tolerance, scaled_zero * scaled_zero);
        if p1.distance_squared(p2) >= tolerance_sq {
            return false;
        }
        let lsq = interval.length() * interval.length();
        for t in [interval.min, interval.max, interval.evaluate(0.4567)] {
            let Ok(pv) = crate::hermite::Hermite::evaluate(self, t, 1, false) else {
                return false;
            };
            if pv[1].length_squared() * lsq > tolerance_sq {
                return false;
            }
            if t == interval.evaluate(0.4567) && p1.distance_squared(pv[0]) > tolerance_sq {
                return false;
            }
        }
        true
    }
}

/// `sm_IsPointOnCurve`: whether `test` lies on the Hermite span from `prev`
/// to `next`, within `tolerance` grown for shallow meetings, and either
/// against the span's direction or inside it within half the span.
pub fn is_point_on_span(
    test: &TsectPoint,
    prev: &TsectPoint,
    next: &TsectPoint,
    tolerance: f64,
) -> Result<bool, SolverRefusal> {
    use crate::hermite::Hermite;
    let mut tolerance = tolerance;
    let angle = max(prev.tangent_plane_angle, next.tangent_plane_angle);
    if angle > 0.0 {
        let mut tangent = (angle / 2.0).tan();
        if tangent < 0.001 {
            tangent = 0.001;
        }
        tolerance /= tangent;
    }
    let delta = next.curve_parameter - prev.curve_parameter;
    let (p1, p2) = (prev.curve_point, next.curve_point);
    let span = Hermite::new(
        p1,
        prev.curve_tangent * delta,
        p2,
        next.curve_tangent * delta,
    );
    let length = (p2 - p1).max_dimension();
    let (span_box, _) = span.bounding_boxes();
    let mut point_box = crate::extent::Extent3::singleton(test.curve_point);
    point_box.expand_absolute(tolerance * 2.0 + length / 4.0);
    if span_box.are_disjoint(point_box, 0.0) {
        return Ok(false);
    }
    for s in 0..2 {
        let uv_span = Hermite::planar(
            prev.uv[s],
            prev.uv_deriv[s] * delta,
            next.uv[s],
            next.uv_deriv[s] * delta,
        );
        let (mut uv_box, _) = uv_span.bounding_boxes();
        let diagonal = (uv_box.max - uv_box.min).length();
        for axis in 0..2 {
            let (lo, hi) = if axis == 0 {
                (uv_box.min.x, uv_box.max.x)
            } else {
                (uv_box.min.y, uv_box.max.y)
            };
            let mut grow = hi - lo;
            if grow < diagonal / 100.0 {
                grow = diagonal / 100.0;
            }
            if axis == 0 {
                uv_box.min.x = lo - grow;
                uv_box.max.x = hi + grow;
            } else {
                uv_box.min.y = lo - grow;
                uv_box.max.y = hi + grow;
            }
        }
        let uv = Vector3::new(test.uv[s].x, test.uv[s].y, 0.0);
        if !uv_box.contains_point(uv, 0.0) {
            return Ok(false);
        }
    }
    let interval = span.natural_interval();
    let Some(found) = crate::curve_solve::local_point_solve(
        &span,
        interval,
        SolverOperation::Minimize,
        test.curve_point,
        interval.evaluate(0.5),
    )?
    else {
        return Ok(false);
    };
    if found.distance <= tolerance {
        let pv = span.evaluate(found.parameter, 1, true)?;
        if test.curve_tangent.dot(pv[1]) < 0.0 {
            return Ok(true);
        }
    }
    if interval.contains_value(found.parameter, -EFF_ZERO) {
        return Ok(found.distance <= length / 2.0);
    }
    Ok(false)
}

/// The tracing state of `SmSurfaceIntersector`: the curve being traced, the
/// start and through points, and the trace's flags.
#[derive(Clone, Debug, Default)]
pub struct Tracer {
    /// `m_vCurvePoints`.
    pub curve_points: Vec<TsectPoint>,
    /// `m_vStartPoints`, in list order.
    pub start_points: Vec<TsectPoint>,
    /// `m_vThroughPoints`, in list order.
    pub through_points: Vec<TsectPoint>,
    /// `m_pStartPoint`, which is the first curve point when set.
    pub has_start_point: bool,
    /// `m_pEndPoint`.
    pub end_point: Option<TsectPoint>,
    /// `m_dCurveTraceDirection`.
    pub trace_direction: f64,
    /// `m_dLastThroughParmeter`.
    pub last_through_parameter: f64,
    /// `m_dExtensionDistance`.
    pub extension_distance: f64,
    /// `m_bDoingExtensionClipping`.
    pub extension_clipping: bool,
    /// `m_bCurveIsClosed`.
    pub curve_is_closed: bool,
    /// `m_bClippedByThroughPoint`.
    pub clipped_by_through_point: bool,
}

/// Which list `sm_UpdateTsectPoint` snaps to.
#[derive(Clone, Copy)]
enum Target {
    Through(usize),
    Start(usize),
    FirstCurvePoint,
    End,
}

impl Tracer {
    fn target(&self, target: Target) -> TsectPoint {
        match target {
            Target::Through(i) => self.through_points[i],
            Target::Start(i) => self.start_points[i],
            Target::FirstCurvePoint => self.curve_points[0],
            Target::End => self.end_point.unwrap_or_default(),
        }
    }

    fn target_mut(&mut self, target: Target) -> Option<&mut TsectPoint> {
        match target {
            Target::Through(i) => self.through_points.get_mut(i),
            Target::Start(i) => self.start_points.get_mut(i),
            Target::FirstCurvePoint => self.curve_points.first_mut(),
            Target::End => self.end_point.as_mut(),
        }
    }

    /// `sm_UpdateTsectPoint`: moves `next` onto the target point, keeping the
    /// curve parameter proportional to the distance from the last point.
    fn update_point(
        &mut self,
        intersector: &Intersector<'_>,
        next: &mut TsectPoint,
        target: Target,
    ) -> Result<(), SolverRefusal> {
        let prev = *self.curve_points.last().ok_or(Fault::Degenerate)?;
        let distance = prev.curve_point.distance(next.curve_point);
        let step = next.curve_parameter - prev.curve_parameter;
        let snap = self.target(target);
        intersector.compute_point_values(snap.uv, next, Some(&prev), None, self.trace_direction)?;
        let new_distance = prev.curve_point.distance(next.curve_point);
        next.curve_parameter = prev.curve_parameter + new_distance / distance * step;
        if let Some(point) = self.target_mut(target)
            && point.curve_parameter < 0.0
        {
            point.curve_parameter = next.curve_parameter;
        }
        Ok(())
    }

    /// `TestPoint`: snaps `next` to the first through point, a start point,
    /// the trace's own start (closing the curve) or its end point when it
    /// lies on the span from the last point, and reports whether the trace
    /// is done there or because `next` leaves a domain across its side.
    #[allow(clippy::too_many_lines)]
    pub fn test_point(
        &mut self,
        intersector: &Intersector<'_>,
        next: &mut TsectPoint,
    ) -> Result<bool, SolverRefusal> {
        let tolerance = intersector.approx_tolerance;
        let prev_index = self
            .curve_points
            .len()
            .checked_sub(1)
            .ok_or(Fault::Degenerate)?;
        if !self.through_points.is_empty() {
            let prev = self.curve_points[prev_index];
            if is_point_on_span(&self.through_points[0], &prev, next, tolerance)? {
                self.update_point(intersector, next, Target::Through(0))?;
                let through = self.through_points.remove(0);
                if self.extension_clipping && self.through_points.is_empty() {
                    self.last_through_parameter = through.curve_parameter;
                }
                self.clipped_by_through_point = true;
                return Ok(false);
            }
        }
        for i in (0..self.start_points.len()).rev() {
            let prev = self.curve_points[prev_index];
            if self.has_start_point && prev_index == 0 {
                continue;
            }
            let start = self.start_points[i];
            if prev.curve_point.distance(start.curve_point) < tolerance {
                let size1 = (intersector.domains[0].max - intersector.domains[0].min).length();
                let size2 = (intersector.domains[1].max - intersector.domains[1].min).length();
                if prev.uv[0].distance(start.uv[0]) < size1 / 2.0
                    && prev.uv[1].distance(start.uv[1]) < size2 / 2.0
                {
                    continue;
                }
            }
            if is_point_on_span(&start, &prev, next, tolerance)? {
                self.update_point(intersector, next, Target::Start(i))?;
                let start = self.start_points[i];
                if next.kind == PointKind::Crossing
                    && start.tangent_plane_angle > intersector.stepping_parallel_tolerance
                {
                    self.start_points.remove(i);
                }
                if start.kind == PointKind::Crossing {
                    let through = |s: usize, a: usize| {
                        let along = start.param_types[s][a];
                        let across = start.param_types[s][1 - a];
                        let d = if a == 0 {
                            start.uv_deriv[s].x
                        } else {
                            start.uv_deriv[s].y
                        };
                        (along == PointType::Start || along == PointType::End)
                            && across == PointType::Inside
                            && crate::scalar::is_zero(d)
                    };
                    let any = through(0, 0) || through(1, 0) || through(0, 1) || through(1, 1);
                    return Ok(!any);
                }
                return Ok(true);
            }
        }
        if self.has_start_point && prev_index != 0 && prev_index != 1 {
            let prev = self.curve_points[prev_index];
            if is_point_on_span(&self.curve_points[0], &prev, next, tolerance)? {
                self.update_point(intersector, next, Target::FirstCurvePoint)?;
                self.curve_is_closed = true;
                return Ok(true);
            }
        }
        if let Some(end) = self.end_point {
            let prev = self.curve_points[prev_index];
            if is_point_on_span(&end, &prev, next, tolerance)? {
                let on_seam =
                    prev_index == 0 && prev.curve_point.distance(end.curve_point) < tolerance;
                if !on_seam {
                    self.update_point(intersector, next, Target::End)?;
                    return Ok(true);
                }
            }
        } else if self.extension_clipping
            && self.last_through_parameter != 0.0
            && next.curve_parameter > self.last_through_parameter + self.extension_distance
        {
            return Ok(true);
        }
        let mut on_boundary = false;
        for s in 0..2 {
            let domain = intersector.domains[s];
            let size = max(domain.min.max_dimension(), domain.max.max_dimension());
            let tol = EFF_ZERO * (1.0 + size) * 100.0;
            let angle_tol = EFF_ZERO_SQRT * (1.0 + size) * 10.0;
            let (on, binormal) = domain.point_on_boundary(next.uv[s], tol);
            if on
                && ((binormal.x != 0.0 && next.uv_deriv[s].x.abs() > angle_tol)
                    || (binormal.y != 0.0 && next.uv_deriv[s].y.abs() > angle_tol))
            {
                on_boundary = true;
                break;
            }
        }
        if !on_boundary {
            return Ok(false);
        }
        let near = |p: &TsectPoint| {
            (0..2).all(|k| {
                let domain = intersector.domains[k];
                let size = max(domain.min.max_dimension(), domain.max.max_dimension());
                let tol = EFF_ZERO * (1.0 + size) * 100.0;
                (p.uv[k] - next.uv[k]).length_squared() <= tol
            })
        };
        self.start_points.retain(|p| !near(p));
        Ok(true)
    }

    /// `TraceCurve`: steps from `start` along the meeting curve, shrinking a
    /// step that fails the span tests and growing one that passes easily,
    /// until a side, a start or end point, or the curve's own start. Points
    /// go to `curve_points`. `existing` stands for the curves traced before
    /// (`IsPointOnCurve`), which this part does not yet consult.
    #[allow(clippy::too_many_lines)]
    pub fn trace_curve(
        &mut self,
        intersector: &Intersector<'_>,
        start: &TsectPoint,
        existing: usize,
    ) -> Result<(), SolverRefusal> {
        let not_ported = || SolverRefusal::NotPorted(SolverOperation::Intersect);
        let uv_tol = |s: usize| {
            let d = intersector.domains[s];
            EFF_ZERO * (1.0 + max(d.min.max_dimension(), d.max.max_dimension()))
        };
        if existing != 0 {
            return Err(not_ported());
        }
        if self.curve_points.is_empty() {
            let mut first = *start;
            if first.kind == PointKind::TangentPoint || first.kind == PointKind::Coincidence {
                first.curve_parameter = 0.0;
                first.deviation = 0.0;
                self.curve_points.push(first);
                self.curve_is_closed = true;
                return Ok(());
            }
            first.curve_parameter = 0.0;
            first.deviation = 0.0;
            self.curve_points.push(first);
            let probe = [0, 1].map(|s| first.uv[s] + (10.0 * EFF_ZERO_SQRT) * first.uv_deriv[s]);
            if !intersector.domains[0].contains_point(probe[0], uv_tol(0))
                || !intersector.domains[1].contains_point(probe[1], uv_tol(1))
            {
                // `ReverseCurveDirection`.
                return Err(not_ported());
            }
        }
        self.curve_is_closed = false;
        let mut count = 0usize;
        let mut step = 0.0;
        // SMLib's `dDeviation` lives across the whole trace and keeps its
        // value when a step returns before the span test writes it.
        let mut deviation = 0.0;
        let check = |p: &TsectPoint| -> Result<(), SolverRefusal> {
            if !intersector.domains[0].contains_point(p.uv[0], uv_tol(0))
                || !intersector.domains[1].contains_point(p.uv[1], uv_tol(1))
            {
                return Err(Fault::Degenerate.into());
            }
            Ok(())
        };
        loop {
            let last = *self.curve_points.last().ok_or(Fault::Degenerate)?;
            step = intersector.compute_step_size(step);
            let mut next = TsectPoint::default();
            let mut outcome = match intersector.compute_next_point(
                &last,
                step,
                &mut next,
                self.trace_direction,
            ) {
                Ok(outcome) => outcome,
                Err(SolverRefusal::Failed(_)) => break,
                Err(other) => return Err(other),
            };
            if outcome.boundary_hit && step < intersector.approx_tolerance * 100.0 {
                break;
            }
            if outcome.good {
                check(&next)?;
            }
            if outcome.deviation_written {
                deviation = outcome.deviation;
            }
            let mut decreased = false;
            while !outcome.good {
                step /= 1.5;
                if step < EFF_ZERO * 100.0 {
                    break;
                }
                decreased = true;
                outcome = match intersector.compute_next_point(
                    &last,
                    step,
                    &mut next,
                    self.trace_direction,
                ) {
                    Ok(outcome) => outcome,
                    Err(SolverRefusal::Failed(_)) => {
                        outcome.good = false;
                        break;
                    }
                    Err(other) => return Err(other),
                };
                if outcome.deviation_written {
                    deviation = outcome.deviation;
                }
                if outcome.good {
                    check(&next)?;
                }
            }
            if !outcome.good {
                break;
            }
            let last_parameter = next.curve_parameter;
            let done = self.test_point(intersector, &mut next)?;
            if done || next.curve_parameter != last_parameter {
                self.curve_points.push(next);
                if done {
                    break;
                }
                continue;
            }
            if !outcome.clipped
                && !decreased
                && deviation < intersector.approx_tolerance / 8.0
                && outcome.angle < intersector.angle_tolerance / 4.0
            {
                let mut last_parameter = next.curve_parameter;
                while outcome.good && !outcome.clipped {
                    step *= 1.5;
                    outcome = match intersector.compute_next_point(
                        &last,
                        step,
                        &mut next,
                        self.trace_direction,
                    ) {
                        Ok(outcome) => outcome,
                        Err(SolverRefusal::Failed(_)) => {
                            outcome.good = false;
                            break;
                        }
                        Err(other) => return Err(other),
                    };
                    if outcome.deviation_written {
                        deviation = outcome.deviation;
                    }
                    if !outcome.good {
                        continue;
                    }
                    check(&next)?;
                    if next.curve_parameter <= last_parameter + EFF_ZERO {
                        outcome.good = false;
                        break;
                    }
                    last_parameter = next.curve_parameter;
                    let done = self.test_point(intersector, &mut next)?;
                    if done || next.curve_parameter != last_parameter {
                        if done {
                            self.curve_points.push(next);
                            return Ok(());
                        }
                        break;
                    }
                }
                if !outcome.good {
                    step /= 1.5;
                    outcome = intersector.compute_next_point(
                        &last,
                        step,
                        &mut next,
                        self.trace_direction,
                    )?;
                    if outcome.deviation_written {
                        deviation = outcome.deviation;
                    }
                    if outcome.good {
                        check(&next)?;
                    }
                }
            }
            self.clipped_by_through_point = false;
            if step < EFF_ZERO * 100.0 {
                break;
            }
            let moved = next.curve_point.distance(last.curve_point);
            if moved < EFF_ZERO * 100.0 * (1.0 + last.curve_point.max_dimension()) {
                break;
            }
            let done = self.test_point(intersector, &mut next)?;
            next.deviation = deviation;
            self.curve_points.push(next);
            if done {
                break;
            }
            count += 1;
            if count == 5000 {
                return Err(Fault::Degenerate.into());
            }
        }
        Ok(())
    }
}

/// Whether the point lies on a side of both domains with a positive
/// deviation, as `AddStartPoint` classifies start points.
fn bounded(point: &TsectPoint, scaled_zero: f64) -> bool {
    let on_side = |t: PointType| matches!(t, PointType::Start | PointType::End | PointType::Both);
    point.deviation > scaled_zero
        && (on_side(point.param_types[0][0]) || on_side(point.param_types[0][1]))
        && (on_side(point.param_types[1][0]) || on_side(point.param_types[1][1]))
}

/// `smsurf_EvaluateNormalSection` on the point's own evaluation matrix,
/// which SMLib updates in place, along a parameter direction read as a 3D
/// vector with zero `z` (as SMLib converts it). A failure leaves zero.
fn normal_section_curvature(point: &mut TsectPoint, s: usize) -> f64 {
    let mut jet: crate::curve_surface::SurfaceJet = [[Vector3::default(); 3]; 3];
    for k in 0..3 {
        for l in 0..3 {
            jet[k][l] = point.surface[s][k * 3 + l];
        }
    }
    let d = point.uv_deriv[s];
    let result = crate::curve_surface::normal_section_of(&mut jet, Vector3::new(d.x, d.y, 0.0));
    for k in 0..3 {
        for l in 0..3 {
            point.surface[s][k * 3 + l] = jet[k][l];
        }
    }
    result.map_or(0.0, |(_, curvature)| curvature)
}

impl Tracer {
    /// `AddStartPoint`: adds `start` (crossing points last, others first)
    /// unless an existing start point matches it: the same parameters, a
    /// bounded and an unbounded point close in space and parameters with
    /// tangents within a curvature-scaled angle, or nearby parameters whose
    /// surface points meet.
    pub fn add_start_point(
        &mut self,
        intersector: &Intersector<'_>,
        mut start: TsectPoint,
    ) -> Result<(), SolverRefusal> {
        let tolerance = intersector.approx_tolerance;
        let scaled_zero = EFF_ZERO * 100.0 * (1.0 + start.position(0).max_dimension());
        let bounded_start = bounded(&start, scaled_zero);
        let mut matched = false;
        for i in 0..self.start_points.len() {
            let target = self.start_points[i];
            let bounded_target = bounded(&target, scaled_zero);
            let d1 = (start.uv[0] - target.uv[0]).length_squared();
            let d2 = (start.uv[1] - target.uv[1]).length_squared();
            if d1 < EFF_ZERO_SQ && d2 < EFF_ZERO_SQ {
                matched = true;
                break;
            }
            if bounded_start != bounded_target
                && start.position(0).distance(target.position(0)) < 2.0 * tolerance
                && start.position(1).distance(target.position(1)) < 2.0 * tolerance
            {
                let size = |s: usize| {
                    let d = intersector.surfaces[s].natural_domain();
                    max(d.min.max_dimension(), d.max.max_dimension()) / 10.0
                };
                if start.uv[0].distance(target.uv[0]) < size(0)
                    && start.uv[1].distance(target.uv[1]) < size(1)
                {
                    let c0 = normal_section_curvature(&mut start, 0);
                    let c1 = normal_section_curvature(&mut start, 1);
                    let curvature = max(c0.abs(), c1.abs());
                    let radius = if crate::scalar::is_zero(curvature) {
                        100.0 * tolerance
                    } else {
                        min(1.0 / curvature, 100.0 * tolerance)
                    };
                    let angle_tolerance = 0.01 + tolerance.atan2(radius);
                    let between = start
                        .curve_tangent
                        .angle_between(target.curve_tangent)
                        .unwrap_or(0.0);
                    if between < angle_tolerance {
                        matched = true;
                        break;
                    }
                }
            }
            if d1 < EFF_ZERO && d2 < EFF_ZERO {
                let p1 =
                    crate::surface_eval::evaluate_point(intersector.surfaces[0], target.uv[0])?;
                let p2 =
                    crate::surface_eval::evaluate_point(intersector.surfaces[1], target.uv[1])?;
                if p1.distance(p2) < tolerance {
                    matched = true;
                    break;
                }
            }
        }
        if !matched {
            if start.kind == PointKind::Crossing {
                self.start_points.push(start);
            } else {
                self.start_points.insert(0, start);
            }
        }
        Ok(())
    }
}

impl Intersector<'_> {
    /// `ClassifyBoundaryStartPointIntersection`: whether a meeting found on
    /// a side curve is only a touch, not the start of an intersection curve.
    /// Tangent points touch; points whose curve leaves both domains both ways
    /// touch; tangent, singular and coincident points start curves; and a
    /// crossing on a side tangent to the other surface is probed along the
    /// curve both ways to see whether the gap opens monotonically away from
    /// the side (a touch) or crosses over.
    #[allow(clippy::too_many_lines)]
    pub fn classify_boundary_start(&self, point: &TsectPoint) -> Result<bool, SolverRefusal> {
        if point.kind == PointKind::TangentPoint {
            return Ok(true);
        }
        // SMLib adds and subtracts the scaled derivative; the same here.
        let p_inside = (0..2).all(|s| {
            self.domains[s].contains_point(
                point.uv[s] + EFF_ZERO_SQRT * point.uv_deriv[s],
                10.0 * EFF_ZERO,
            )
        });
        let n_inside = (0..2).all(|s| {
            self.domains[s].contains_point(
                point.uv[s] - EFF_ZERO_SQRT * point.uv_deriv[s],
                10.0 * EFF_ZERO,
            )
        });
        if !p_inside && !n_inside {
            return Ok(true);
        }
        if matches!(
            point.kind,
            PointKind::TangentCurve | PointKind::Singularity | PointKind::Coincidence
        ) {
            return Ok(false);
        }
        let angle_tolerance = 2.0 * PI / 180.0;
        let gap = point.position(1) - point.position(0);
        let mut tangent_u = [false; 2];
        let mut tangent_v = [false; 2];
        let mut boundary_u = [false; 2];
        let mut boundary_v = [false; 2];
        let mut binormal_u = [Vector3::default(); 2];
        let mut binormal_v = [Vector3::default(); 2];
        for s in 0..2 {
            let uv = point.uv[s];
            let domain = self.domains[s];
            let other_normal = point.normal(1 - s);
            let (su, sv) = (point.du(s), point.dv(s));
            if crate::scalar::are_same(domain.min.x, uv.x) {
                boundary_u[s] = true;
                tangent_u[s] = other_normal.is_perpendicular_to(sv, angle_tolerance);
                binormal_u[s] = su;
            } else if crate::scalar::are_same(domain.max.x, uv.x) {
                boundary_u[s] = true;
                tangent_u[s] = other_normal.is_perpendicular_to(sv, angle_tolerance);
                binormal_u[s] = -su;
            }
            if crate::scalar::are_same(domain.min.y, uv.y) {
                boundary_v[s] = true;
                tangent_v[s] = other_normal.is_perpendicular_to(su, angle_tolerance);
                binormal_v[s] = sv;
            } else if crate::scalar::are_same(domain.max.y, uv.y) {
                boundary_v[s] = true;
                tangent_v[s] = other_normal.is_perpendicular_to(su, angle_tolerance);
                binormal_v[s] = -sv;
            }
        }
        let on_tangent =
            (0..2).any(|s| (boundary_u[s] && tangent_u[s]) || (boundary_v[s] && tangent_v[s]));
        if !on_tangent {
            return Ok(false);
        }
        let scaled_zero = EFF_ZERO * (1.0 + point.curve_point.max_dimension());
        for s in 0..2 {
            if !tangent_u[s] && !tangent_v[s] {
                continue;
            }
            let (binormal, tangent) = if tangent_u[s] {
                (binormal_u[s], point.dv(s))
            } else {
                (binormal_v[s], point.du(s))
            };
            let stride = self.approx_tolerance / tangent.length() / 2.0 * 100.0;
            let mut side_degenerate = [false; 2];
            for side in 0..2 {
                let (mut zero, mut mono, mut div) = (true, true, true);
                if (side == 0 && !p_inside) || (side == 1 && !n_inside) {
                    side_degenerate[side] = true;
                    continue;
                }
                let mut curve_point = [Vector3::default(); 3];
                let mut gaps = [Vector3::default(); 3];
                let mut gap_length = [0.0; 3];
                let mut inside = [false; 3];
                let mut sample_distance = 0.0;
                curve_point[0] = point.position(s);
                gaps[0] = Vector3::new(0.0, 0.0, 0.0);
                gap_length[0] = 0.0;
                inside[0] = true;
                let mut delta = if side == 0 { stride } else { -stride };
                let mut ii = 1usize;
                loop {
                    let tgt = if ii > 1 { 2 } else { 1 };
                    inside[tgt] = false;
                    let uv0 = point.uv[s] + delta * point.uv_deriv[s];
                    let uv1 = point.uv[1 - s] + delta * point.uv_deriv[1 - s];
                    let mut other_pv = [Vector3::default(); 3];
                    let mut drop_uv = Vector2::new(0.0, 0.0);
                    if self.domains[s].contains_point(uv0, 10.0 * EFF_ZERO)
                        && self.domains[1 - s].contains_point(uv1, 10.0 * EFF_ZERO)
                    {
                        curve_point[tgt] =
                            crate::surface_eval::evaluate_point(self.surfaces[s], uv0)
                                .unwrap_or_default();
                        sample_distance = (point.position(1) - curve_point[tgt]).length();
                        let found = crate::surface_solve::local_point_solve(
                            self.surfaces[1 - s],
                            self.domains[1 - s],
                            SolverOperation::Normalize,
                            curve_point[tgt],
                            uv1,
                            None,
                        );
                        if let Ok(found) = found
                            && found.found
                        {
                            drop_uv = found.uv;
                            if let Ok(cells) = evaluate(
                                self.surfaces[1 - s],
                                drop_uv,
                                1,
                                1,
                                [true, true],
                                true,
                                true,
                                true,
                            ) {
                                // `[0][0]`, `[1][0]` (Su) and `[0][1]` (Sv).
                                other_pv = [cells[0], cells[2], cells[1]];
                            }
                            let shift = if s == 0 { -gap } else { gap };
                            gaps[tgt] = (other_pv[0] - curve_point[tgt]) + shift;
                            gap_length[tgt] = gaps[tgt].length();
                            inside[tgt] = true;
                        }
                    }
                    if !inside[tgt] {
                        let previous = (tgt + 3 - 1) % 3;
                        side_degenerate[side] = ii < 3
                            || (mono
                                && div
                                && binormal.dot(gaps[previous]) > 0.0
                                && sample_distance < self.approx_tolerance * 10.0);
                        break;
                    } else if ii >= 2 {
                        let tgt1 = (tgt + 3 - 1) % 3;
                        let tgt2 = (tgt + 3 - 2) % 3;
                        if gap_length[tgt] > scaled_zero {
                            zero = false;
                        }
                        let bxt = binormal.cross(tangent);
                        let sc = other_pv[0] - curve_point[tgt];
                        let (su, sv) = (other_pv[1], other_pv[2]);
                        let mut a = crate::matrix::Matrix::new(3, 3);
                        a[0][0] = 2.0 * su.dot(su);
                        a[0][1] = 2.0 * su.dot(sv);
                        a[0][2] = -bxt.dot(su);
                        a[1][0] = 2.0 * sv.dot(su);
                        a[1][1] = 2.0 * sv.dot(sv);
                        a[1][2] = -bxt.dot(sv);
                        a[2][0] = bxt.dot(su);
                        a[2][1] = bxt.dot(sv);
                        a[2][2] = 0.0;
                        let b = [-2.0 * sc.dot(su), -2.0 * sc.dot(sv), -bxt.dot(sc)];
                        let mut x = Vec::new();
                        let open = match a.solve_linear_system(&b, &mut x) {
                            Err(_) => true,
                            Ok(_) => {
                                let puv = Vector2::new(drop_uv.x + x[0], drop_uv.y + x[1]);
                                let pc = other_pv[0] + x[0] * su + x[1] * sv - curve_point[tgt];
                                !self.domains[1 - s].contains_point(puv, EFF_ZERO)
                                    || pc.dot(binormal) < -EFF_ZERO
                            }
                        };
                        if !zero {
                            if gaps[tgt].dot(gaps[tgt1]) <= 0.0 {
                                mono = false;
                            }
                            if gap_length[tgt] == 0.0
                                || gap_length[tgt]
                                    < gap_length[tgt1]
                                        + 1.0001 * (gap_length[tgt1] - gap_length[tgt2])
                            {
                                div = false;
                            }
                        }
                        if !zero && mono && div && open {
                            side_degenerate[side] = true;
                            break;
                        } else if zero && mono && div {
                            curve_point[1] = curve_point[2];
                            gaps[1] = gaps[2];
                            gap_length[1] = gap_length[2];
                            inside[1] = inside[2];
                        } else {
                            side_degenerate[side] = false;
                            break;
                        }
                    }
                    ii += 1;
                    delta *= 5.0;
                }
            }
            if side_degenerate[0] && side_degenerate[1] {
                return Ok(true);
            }
        }
        Ok(false)
    }
}
