// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: the plane is a B-spline surface carrying a
// `PlaneFrame` (see `spline.rs`), as SMLib's `SmPlane` derives from its
// B-spline surface; statuses become `Fault` and the solution container a
// returned value.
// Translated from SMLib SmPlane.cpp (DropPointFast) and SmSurface.cpp
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
