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
