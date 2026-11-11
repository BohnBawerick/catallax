// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline surfaces (and planes) of
// `spline.rs` without a trimming face, so the drop tolerance is SMLib's
// default intersection tolerance; statuses become `Fault`, output references
// a returned struct. The seam-side selection and the solution-array output
// are not carried: no caller in the port passes them.
// Translated from SMLib SmSurface.cpp (SmSurface::DropPoint,
// sm_GoodDropResult); rust-brep has no port.
//! The surface point drop of step 7: the nearest surface point to a point,
//! from a guess when it is good enough and otherwise by the global search.

use crate::extent::Extent2;
use crate::global_solver::{SolverOperation, SolverRefusal};
use crate::plane::drop_point_fast;
use crate::scalar::{BIG_DOUBLE, EFF_ZERO, ZONE_TOL_3D};
use crate::solution::{Solution, SolutionType};
use crate::spline::BSplineSurface;
use crate::surface_eval::{evaluate_point, first_derivatives};
use crate::surface_global::global_point_solve;
use crate::surface_solve::local_point_solve;
use crate::vector::{Vector2, Vector3};

/// `SM_XSECT_TOL_3D`: the drop tolerance of a surface without a face.
pub const XSECT_TOL_3D: f64 = 2.0 * ZONE_TOL_3D;

/// What [`drop_point`] found: SMLib's `rbSuccess`, `rResultUV`, `rdGap` and
/// `rbIsMultivalued`. Without success `uv` is the guess (or zero) and `gap`
/// is `SM_BIG_DOUBLE`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceDrop {
    /// The point dropped.
    pub success: bool,
    /// The parameters of the drop.
    pub uv: Vector2,
    /// The distance from the point to the surface there.
    pub gap: f64,
    /// The search found several answers, or a range.
    pub multivalued: bool,
}

/// `sm_GoodDropResult`: whether a local drop from `guess` converged well
/// enough to skip the global search. SMLib errs on the side of searching.
fn good_drop_result(
    surface: &BSplineSurface,
    domain: Extent2,
    point: Vector3,
    guess: Vector2,
    uv: Vector2,
    value: f64,
) -> bool {
    if value < EFF_ZERO {
        return true;
    }
    const MAX_STEP: f64 = 0.50;
    let step = uv - guess;
    if step.x.abs() > MAX_STEP * (domain.max.x - domain.min.x)
        || step.y.abs() > MAX_STEP * (domain.max.y - domain.min.y)
    {
        return false;
    }
    // SMLib ignores the evaluation status and keeps its zero-initialised
    // vectors on failure.
    let zero = Vector3::new(0.0, 0.0, 0.0);
    let [surface_point, mut su, mut sv] =
        first_derivatives(surface, uv, [true, true]).unwrap_or([zero; 3]);
    let _ = su.unitize();
    let _ = sv.unitize();
    let gap = point - surface_point;
    let fu = gap.dot(su);
    let fv = gap.dot(sv);
    if fu.abs() + fv.abs() < 0.0001 {
        return true;
    }
    for (parameter, low, high, f) in [
        (uv.x, domain.min.x, domain.max.x, fu),
        (uv.y, domain.min.y, domain.max.y, fv),
    ] {
        if parameter < low + EFF_ZERO {
            if f > EFF_ZERO {
                return false;
            }
        } else if parameter > high - EFF_ZERO {
            if f < -EFF_ZERO {
                return false;
            }
        } else if f.abs() > 0.001 {
            return false;
        }
    }
    let guess_point = evaluate_point(surface, guess).unwrap_or(zero);
    // SMLib's names are crossed: the guess's distance must not exceed the
    // result's.
    let result_distance = point.distance_squared(surface_point);
    let guess_distance = point.distance_squared(guess_point);
    guess_distance <= result_distance
}

/// `SmSurface::DropPoint` without a face: `DropPointFast` (planes only), then
/// a local solve from `guess` kept when [`good_drop_result`] accepts it,
/// then the global search. Among several global answers the one nearest the
/// guess wins, or without a guess the closest one.
pub fn drop_point(
    surface: &BSplineSurface,
    domain: Extent2,
    point: Vector3,
    guess: Option<Vector2>,
    operation: SolverOperation,
) -> Result<SurfaceDrop, SolverRefusal> {
    let mut drop = SurfaceDrop {
        success: false,
        uv: guess.unwrap_or(Vector2::new(0.0, 0.0)),
        gap: BIG_DOUBLE,
        multivalued: false,
    };
    let tolerance = XSECT_TOL_3D;
    let mut answers: Vec<Solution> = Vec::new();
    let mut local: Option<(Vector2, f64)> = None;
    if surface.plane_frame().is_some()
        && let Ok(Some((uv, distance))) =
            drop_point_fast(surface, domain, operation, point, tolerance, None)
    {
        let mut solution = Solution {
            kind: SolutionType::SingleValue,
            variables: 2,
            ..Solution::default()
        };
        solution.start.parameters[0] = uv.x;
        solution.start.parameters[1] = uv.y;
        solution.start.value = distance;
        answers.push(solution);
        drop.success = domain.contains_point(uv, 0.0);
    }
    if !drop.success
        && let Some(guess) = guess
    {
        let found = match local_point_solve(surface, domain, operation, point, guess, None) {
            Ok(found) => found.found.then_some((found.uv, found.value)),
            Err(SolverRefusal::Failed(_)) => None,
            Err(refusal) => return Err(refusal),
        };
        if let Some((uv, value)) = found
            && good_drop_result(surface, domain, point, guess, uv, value)
        {
            drop.success = true;
            local = Some((uv, value));
        }
    }
    let mut from_array = drop.success && local.is_none();
    if !drop.success {
        let found = global_point_solve(surface, domain, operation, point, tolerance, None)?;
        answers = found.as_slice().to_vec();
        if !answers.is_empty() {
            drop.success = true;
            from_array = true;
        }
    }
    if let Some((uv, value)) = local {
        drop.uv = uv;
        drop.gap = value;
    }
    if from_array {
        drop.multivalued = answers.len() > 1 || answers[0].kind == SolutionType::RangeOfValues;
        let mut best_value = BIG_DOUBLE;
        let mut best = None;
        let mut best_uv_distance = BIG_DOUBLE;
        let mut nearest = None;
        for (index, answer) in answers.iter().enumerate() {
            if answer.start.value < best_value {
                best_value = answer.start.value;
                best = Some(index);
            }
            if let Some(guess) = guess {
                let uv = Vector2::new(answer.start.parameters[0], answer.start.parameters[1]);
                let distance = (guess - uv).length_squared();
                if distance < best_uv_distance {
                    best_uv_distance = distance;
                    nearest = Some(index);
                }
            }
        }
        match nearest.or(best) {
            Some(index) => {
                let answer = answers[index];
                drop.uv = Vector2::new(answer.start.parameters[0], answer.start.parameters[1]);
                drop.gap = answer.start.value;
            }
            None => drop.success = false,
        }
    }
    Ok(drop)
}
