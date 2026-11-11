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
    ) -> Result<(bool, f64, f64), SolverRefusal> {
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
            return Ok((false, 0.0, angle));
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
                return Ok((false, crate::scalar::sqrt(max_sq), angle));
            }
        }
        Ok((true, crate::scalar::sqrt(max_sq), angle))
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
        outcome.deviation = deviation;
        outcome.angle = angle;
        Ok(outcome)
    }
}
