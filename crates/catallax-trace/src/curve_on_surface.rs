// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: the parameter-space curve is any curve the local
// point solver can evaluate and the surface a checked B-spline surface,
// statuses become `Fault`, and only the evaluations the global solver reads
// are carried: the point, the surface parameters and the surface normal.
// Translated directly from SMLib SmCrvOnSurf.cpp (`EvaluatePoint`,
// `EvaluateUVPoint`, `EvaluateSurfaceNormal`); rust-brep has no port.
//! SMLib's `SmCrvOnSurf`: a curve in a surface's parameter space, read as a
//! curve on the surface. The solver that culls branches for
//! `SM_SO_SURFACE_PROJECTED_INTERSECT` samples it for points and normals.
//!
//! The rest of the class (derivative lifting, bounding boxes, the B-spline
//! approximation, dropping, trimming, transforms and database I/O) is not
//! carried.
//!
//! Everything returned is an unproved candidate.

use crate::{
    Fault,
    curve_solve::PointSolveCurve,
    spline::BSplineSurface,
    surface_eval::{evaluate_normal, evaluate_point},
    vector::{Vector2, Vector3},
};

/// A curve in the parameter space of a surface (`SmCrvOnSurf`).
#[derive(Clone, Copy)]
pub struct CurveOnSurface<'a> {
    uv_curve: &'a dyn PointSolveCurve,
    surface: &'a BSplineSurface,
}

impl<'a> CurveOnSurface<'a> {
    /// `SetAll`: the parameter-space curve and the surface it lies in.
    pub fn new(uv_curve: &'a dyn PointSolveCurve, surface: &'a BSplineSurface) -> Self {
        Self { uv_curve, surface }
    }

    /// `EvaluatePoint`: the surface point at the curve's parameter `t`.
    pub fn evaluate_point(&self, t: f64) -> Result<Vector3, Fault> {
        let uv = self.uv_curve.evaluate_point(t)?;
        evaluate_point(self.surface, Vector2::new(uv.x, uv.y))
    }

    /// `EvaluateUVPoint`: the surface parameters at `t`.
    pub fn evaluate_uv_point(&self, t: f64) -> Result<Vector2, Fault> {
        let uv = self.uv_curve.evaluate_point(t)?;
        Ok(Vector2::new(uv.x, uv.y))
    }

    /// `EvaluateSurfaceNormal`: the surface parameters and the parameter
    /// tangent at `t` (exact, not repaired) and the unit surface normal there.
    pub fn evaluate_surface_normal(
        &self,
        t: f64,
        from_left: bool,
    ) -> Result<([Vector2; 2], Vector3), Fault> {
        let derivatives = self.uv_curve.evaluate(t, 1, from_left, false)?;
        let uv = [
            Vector2::new(derivatives[0].x, derivatives[0].y),
            Vector2::new(derivatives[1].x, derivatives[1].y),
        ];
        let normal = evaluate_normal(self.surface, uv[0], [from_left, from_left]);
        Ok((uv, normal))
    }
}
