// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline surface candidates of
// `spline.rs`; statuses become `Fault`, output references become returned
// values, and the solver's `SmSurfaceEval` evaluations become plain
// second-order jets. The maximising and directed operations refuse with
// `SolverRefusal::NotPorted`, as the local solver they call does, and so does
// a local solve whose singular inversion drops onto a curve that the curve
// solver refuses.
// Adapted from rust-brep crates/kernel/src/surface/sm_cubic_bezier.rs
// (get_triangles_dist, well_behaved_control_net), commit
// 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab, translations of SMLib
// SmSurfacePointSolve.cpp (GetTrianglesDist, WellBehavedControlNet); the rest
// is translated directly from SMLib SmSurfacePointSolve.cpp
// (SmSurfPtGlobalSolver and SmSurface::GlobalPointSolve_1).
//! The global surface point solver of step 5: a guess from the control net
//! or a parameter grid, refined and handed to the local solver, with the
//! answer copied across closed seams.

use crate::Fault;
use crate::extent::{Extent1, Extent2};
use crate::geometry::{line_closest_point, triangle_point_distance, vec_lin_comb_two_vectors};
use crate::global_solver::{SolverOperation, SolverRefusal};
use crate::scalar::{BIG_DOUBLE, EFF_ZERO, EFF_ZERO_SQRT, ZONE_TOL_3D, is_zero};
use crate::solution::{
    ObjectId, ObjectKind, Solution, SolutionArray, SolutionEnd, SolutionObjects, SolutionType,
    SortKey,
};
use crate::spline::BSplineSurface;
use crate::surface_eval::{evaluate_out_of_bounds, evaluate_point, first_derivatives};
use crate::surface_solve::{
    Norm, Problem, Solve, Stepper, is_closed, is_on_boundary, is_on_seam, local_point_solve,
};
use crate::vector::{Vector2, Vector3};

/// The solution object identity of the solved surface.
pub const SURFACE_OBJECT: ObjectId = ObjectId(0);

/// The one surface a global point solve names in its solutions.
pub struct SurfaceObject<'a> {
    /// The surface.
    pub surface: &'a BSplineSurface,
    /// Its zone tolerance.
    pub zone: f64,
}

impl SolutionObjects for SurfaceObject<'_> {
    fn kind(&self, id: ObjectId) -> ObjectKind {
        if id == SURFACE_OBJECT {
            ObjectKind::Surface
        } else {
            ObjectKind::Other
        }
    }

    fn parameter_interval(&self, id: ObjectId, direction: usize) -> Extent1 {
        let domain = self.surface.natural_domain();
        match (id == SURFACE_OBJECT, direction) {
            (true, 0) => domain.u_interval(),
            (true, 1) => domain.v_interval(),
            _ => Extent1::default(),
        }
    }

    fn zone_tolerance(&self, id: ObjectId) -> f64 {
        if id == SURFACE_OBJECT { self.zone } else { 0.0 }
    }
}

/// Position and derivatives as `SmSurfaceEval` holds them.
#[derive(Clone, Copy)]
struct Jet {
    position: Vector3,
    su: Vector3,
    sv: Vector3,
    suu: Vector3,
    svv: Vector3,
}

/// The state of `SmSurfPtGlobalSolver`.
struct Global<'a> {
    surface: &'a BSplineSurface,
    domain: Extent2,
    test: Vector3,
    operation: SolverOperation,
    tight: f64,
    loose: f64,
    owner: Option<f64>,
    zone: f64,
    closed: [Option<bool>; 2],
}

impl Global<'_> {
    /// `SmSurfaceEval::Eval` with the out-of-bounds switch on.
    fn jet(&self, uv: Vector2) -> Result<Jet, Fault> {
        let d = evaluate_out_of_bounds(self.surface, uv, 2, 2, [true, true], true, true, true)?;
        Ok(Jet {
            position: d[0],
            su: d[3],
            sv: d[1],
            suu: d[6],
            svv: d[2],
        })
    }

    /// Control point `(i, j)`, divided by its weight as SMLib's guarded
    /// quotient does.
    fn net(&self, i: usize, j: usize) -> Vector3 {
        let nv = self.surface.counts()[1];
        let [x, y, z, w] = self.surface.control_points()[i * nv + j].coordinates;
        let mut point = Vector3::new(x, y, z);
        if self.surface.is_rational() {
            point /= w;
        }
        point
    }

    fn greville(&self, axis: usize, index: usize) -> f64 {
        self.surface.greville(axis, index).unwrap_or(0.0)
    }

    /// `ClosedSurface`: `IsClosed` over the solve domain, once per direction.
    fn closed_surface(&mut self, axis: usize) -> bool {
        if self.closed[axis].is_none() {
            self.closed[axis] = Some(is_closed(self.surface, self.domain, axis, None));
        }
        self.closed[axis] == Some(true)
    }

    /// `CalcFuncVal`: the distance to the test point, or the projected
    /// residual when normalising; `SM_BIG_DOUBLE` where evaluation fails.
    fn function_value(&self, uv: Vector2) -> f64 {
        if self.operation == SolverOperation::Normalize {
            let Ok([point, mut su, mut sv]) = first_derivatives(self.surface, uv, [false, false])
            else {
                return BIG_DOUBLE;
            };
            if point.is_undefined() {
                return BIG_DOUBLE;
            }
            let gap = point - self.test;
            let _ = su.unitize();
            let _ = sv.unitize();
            let mut value = 0.0;
            value += gap.dot(su).abs();
            value += gap.dot(sv).abs();
            value
        } else {
            let Ok(point) = evaluate_point(self.surface, uv) else {
                return BIG_DOUBLE;
            };
            if point.is_undefined() {
                return BIG_DOUBLE;
            }
            (point - self.test).length()
        }
    }

    /// `GetTrianglesDist`: the nearest approach to the two triangles of
    /// the net square at `(i, j)`, with net-local parameters.
    fn triangles_distance(&self, i: usize, j: usize) -> (f64, f64, f64) {
        let (p00, p10) = (self.net(i, j), self.net(i + 1, j));
        let (p01, p11) = (self.net(i, j + 1), self.net(i + 1, j + 1));
        let first = triangle_point_distance(p00, p10, p01, self.test);
        let second = triangle_point_distance(p11, p01, p10, self.test);
        if first.distance <= second.distance {
            (first.distance, first.u, first.v)
        } else {
            (second.distance, 1.0 - second.u, 1.0 - second.v)
        }
    }

    /// `WellBehavedControlNet`: no sharp turn at net point `(ii, jj)` and
    /// no near-parallel u and v legs.
    fn well_behaved(&self, ii: usize, jj: usize) -> bool {
        const LIMIT0: f64 = -0.61;
        const LIMIT1: f64 = -0.80;
        let [nu, nv] = self.surface.counts();
        let p1 = self.net(ii, jj);
        let turn = |p0: Vector3, p2: Vector3, at_end: bool| -> Option<bool> {
            let (v1, v2) = (p1 - p0, p2 - p1);
            let (l1, l2) = (v1.length(), v2.length());
            if l1 < EFF_ZERO || l2 < EFF_ZERO {
                return Some(at_end);
            }
            if v1.dot(v2) / (l1 * l2) < LIMIT0 {
                return Some(false);
            }
            None
        };
        if ii > 0
            && ii < nu - 1
            && let Some(answer) = turn(
                self.net(ii - 1, jj),
                self.net(ii + 1, jj),
                jj == 0 || jj >= nv - 1,
            )
        {
            return answer;
        }
        if jj > 0
            && jj < nv - 1
            && let Some(answer) = turn(
                self.net(ii, jj - 1),
                self.net(ii, jj + 1),
                ii == 0 || ii >= nu - 1,
            )
        {
            return answer;
        }
        if ii < nu - 1 && jj < nv - 1 {
            let v1 = self.net(ii + 1, jj) - p1;
            let v2 = self.net(ii, jj + 1) - p1;
            let (l1, l2) = (v1.length(), v2.length());
            if l1 < EFF_ZERO {
                return jj == 0 || jj >= nv - 1;
            } else if l2 < EFF_ZERO {
                return ii == 0 || ii >= nu - 1;
            }
            let cosine = v1.dot(v2) / (l1 * l2);
            if cosine < LIMIT1 || cosine > -LIMIT1 {
                return false;
            }
        }
        true
    }

    /// The index form of `CheckGuessForSeams`: at a closed net seam, the
    /// side of the seam the test point lies towards.
    fn check_index_for_seams(&self, best: &mut [usize; 2]) {
        let [nu, nv] = self.surface.counts();
        if best[0] > 0 && best[0] < nu - 1 && best[1] > 0 && best[1] < nv - 1 {
            return;
        }
        if nu > 2 && (best[0] == 0 || best[0] == nu - 1) {
            let seam = self.net(0, best[1]);
            if seam.closer_than(EFF_ZERO, self.net(nu - 1, best[1])) {
                let other = self.net(nu - 2, best[1]);
                let down = line_closest_point(seam, other - seam, self.test);
                let other = self.net(1, best[1]);
                let up = line_closest_point(seam, other - seam, self.test);
                best[0] = if up - down >= 0.0 { 0 } else { nu - 1 };
            }
        }
        if nv > 2 && (best[1] == 0 || best[1] == nv - 1) {
            let seam = self.net(best[0], 0);
            if seam.closer_than(EFF_ZERO, self.net(best[0], nv - 1)) {
                let other = self.net(best[0], nv - 2);
                let down = line_closest_point(seam, other - seam, self.test);
                let other = self.net(best[0], 1);
                let up = line_closest_point(seam, other - seam, self.test);
                best[1] = if up - down >= 0.0 { 0 } else { nv - 1 };
            }
        }
    }

    /// The parameter form of `CheckGuessForSeams`: on a closed seam, the
    /// side the test point lies towards.
    fn check_guess_for_seams(&mut self, uv: &mut Vector2) {
        if !is_on_seam(self.surface, *uv, self.zone).1 {
            return;
        }
        let domain = self.domain;
        for axis in 0..2 {
            let (value, low, high) = if axis == 0 {
                (uv.x, domain.min.x, domain.max.x)
            } else {
                (uv.y, domain.min.y, domain.max.y)
            };
            if (value < low + EFF_ZERO || value > high - EFF_ZERO) && self.closed_surface(axis) {
                let [point, su, sv] = first_derivatives(self.surface, *uv, [false, false])
                    .unwrap_or([Vector3::default(); 3]);
                let gap = point - self.test;
                let derivative = if axis == 0 { su } else { sv };
                let side = if gap.dot(derivative) > 0.0 { high } else { low };
                if axis == 0 {
                    uv.x = side;
                } else {
                    uv.y = side;
                }
            }
        }
    }

    /// `FindPolygonTriangleGuess`: the nearest control-net triangle, its
    /// parameters by Greville abscissae, whether the net is suspicious,
    /// the nearest distance and the mean net leg.
    fn polygon_triangle_guess(&self, best: &mut [usize; 2]) -> (Vector2, f64, f64, bool) {
        let [nu, nv] = self.surface.counts();
        let mut suspicious = self.surface.is_rational();
        let mut best_distance = BIG_DOUBLE;
        let (mut this_u, mut this_v) = (-1.0, -1.0);
        let (mut best_u, mut best_v) = (-1.0, -1.0);
        let mut accumulated = 0.0;
        let mut count = 0usize;
        for ii in 0..nu {
            for jj in 0..nv {
                if !suspicious && !self.well_behaved(ii, jj) {
                    suspicious = true;
                }
                let p0 = self.net(ii, jj);
                let mut distance = BIG_DOUBLE;
                if ii < nu - 1 && jj < nv - 1 {
                    (distance, this_u, this_v) = self.triangles_distance(ii, jj);
                }
                if distance < best_distance {
                    best_distance = distance;
                    *best = [ii, jj];
                    (best_u, best_v) = (this_u, this_v);
                }
                if ii > 0 {
                    accumulated += p0.distance(self.net(ii - 1, jj));
                    count += 1;
                }
                if jj > 0 {
                    accumulated += p0.distance(self.net(ii, jj - 1));
                    count += 1;
                }
            }
        }
        let average = accumulated / count as f64;
        self.check_index_for_seams(best);
        let mut uv = Vector2::new(self.greville(0, best[0]), self.greville(1, best[1]));
        if best_u > EFF_ZERO_SQRT && best[0] < nu - 1 {
            let next = self.greville(0, best[0] + 1);
            uv.x = uv.x + best_u * (next - uv.x);
        }
        if best_v > EFF_ZERO_SQRT && best[1] < nv - 1 {
            let next = self.greville(1, best[1] + 1);
            uv.y = uv.y + best_v * (next - uv.y);
        }
        (
            self.domain.clamp_point(uv),
            best_distance,
            average,
            suspicious,
        )
    }

    /// `RefineGuessPolygonLeg`: the Greville parameters of net point
    /// `best`, shifted towards the nearer of the neighbouring legs.
    fn polygon_leg_guess(&self, best: [usize; 2]) -> Vector2 {
        let [nu, nv] = self.surface.counts();
        let best_point = self.net(best[0], best[1]);
        let mut uv = Vector2::new(self.greville(0, best[0]), self.greville(1, best[1]));
        let shift_to = |neighbour: Vector3| {
            line_closest_point(best_point, neighbour - best_point, self.test).clamp(0.0, 0.5)
        };
        for axis in 0..2 {
            let (index, count) = (best[axis], if axis == 0 { nu } else { nv });
            let at = |k: usize| {
                if axis == 0 {
                    self.net(k, best[1])
                } else {
                    self.net(best[0], k)
                }
            };
            let down = if index > 0 {
                shift_to(at(index - 1))
            } else {
                0.0
            };
            let up = if index < count - 1 {
                shift_to(at(index + 1))
            } else {
                0.0
            };
            let mut shift = up - down;
            if shift.abs() > EFF_ZERO_SQRT {
                let other = if shift > 0.0 {
                    self.greville(axis, index + 1)
                } else {
                    shift = -shift;
                    self.greville(axis, index - 1)
                };
                if axis == 0 {
                    uv.x = (1.0 - shift) * uv.x + shift * other;
                } else {
                    uv.y = (1.0 - shift) * uv.y + shift * other;
                }
            }
        }
        uv
    }

    /// `FindPolygonGuess`: the triangle guess, or when the net is
    /// suspicious the nearest control point refined along its legs.
    fn polygon_guess(&self, best: &mut [usize; 2]) -> (Vector2, bool) {
        let (mut uv, best_distance, average, mut suspicious) = self.polygon_triangle_guess(best);
        if suspicious && best_distance < average / 100.0 {
            let point = evaluate_point(self.surface, uv).unwrap_or_default();
            if point.distance(self.test) < average / 100.0 {
                suspicious = false;
            }
        }
        if !suspicious {
            return (uv, false);
        }
        let [nu, nv] = self.surface.counts();
        let mut tied: Option<[usize; 2]> = None;
        let mut ties = 0usize;
        let mut best_distance = BIG_DOUBLE;
        for ii in 0..nu {
            for jj in 0..nv {
                let distance = self.test.distance(self.net(ii, jj));
                if distance < best_distance - EFF_ZERO {
                    best_distance = distance;
                    *best = [ii, jj];
                    tied = None;
                    ties = 0;
                } else if distance < best_distance + EFF_ZERO {
                    ties += 1;
                    tied = Some(*best);
                    *best = [ii, jj];
                }
            }
        }
        self.check_index_for_seams(best);
        let mut guess = self.polygon_leg_guess(*best);
        let mut found = true;
        let mut best_distance = BIG_DOUBLE;
        for axis in 0..2 {
            let count = if axis == 0 { nu } else { nv };
            let other = 1 - axis;
            if ties >= 2 && ties == count - 1 && tied.map(|t| t[other]) == Some(best[other]) {
                for k in 0..count {
                    let mut index = *best;
                    index[axis] = k;
                    let this_guess = self.polygon_leg_guess(index);
                    let point = evaluate_point(self.surface, this_guess).unwrap_or_default();
                    let distance = point.distance(self.test);
                    if distance < best_distance {
                        best_distance = distance;
                        guess = this_guess;
                        best[axis] = k;
                        if distance < self.loose {
                            break;
                        }
                    }
                }
                found = best_distance < BIG_DOUBLE / 2.0;
                break;
            }
        }
        if found {
            uv = guess;
        } else {
            uv = Vector2::new(self.greville(0, best[0]), self.greville(1, best[1]));
        }
        (self.domain.clamp_point(uv), suspicious)
    }

    /// `CheckBadGuessStride`: a plane step from a grid point, and a second
    /// one when the first is long; `None` where no step exists.
    fn check_bad_guess_stride(&self, uv: Vector2) -> Result<Option<(Vector2, f64)>, Fault> {
        let jet = self.jet(uv)?;
        let gap = self.test - jet.position;
        let Ok(step) = vec_lin_comb_two_vectors(jet.su, jet.sv, gap) else {
            return Ok(None);
        };
        let mut new_uv = self.domain.clamp_point(uv + step);
        let mut new_value = self.function_value(new_uv);
        let max_step_squared = 0.10 * 0.10;
        if step.length_squared() > max_step_squared * self.domain.size().length_squared() {
            let Ok(jet) = self.jet(new_uv) else {
                return Ok(Some((new_uv, new_value)));
            };
            let gap = self.test - jet.position;
            let Ok(step) = vec_lin_comb_two_vectors(jet.su, jet.sv, gap) else {
                return Ok(Some((new_uv, new_value)));
            };
            let next = self.domain.clamp_point(new_uv + step);
            let value = self.function_value(next);
            if value < new_value {
                new_uv = next;
                new_value = value;
            }
        }
        Ok(Some((new_uv, new_value)))
    }

    /// `RefinePatchByGrid`: the best grid point of `domain`, after a plane
    /// step from each, kept only when better than `best`.
    fn refine_patch_by_grid(
        &self,
        domain: Extent2,
        counts: [usize; 2],
        guess: &mut Vector2,
        best: &mut f64,
    ) -> Result<(), Fault> {
        let nu = counts[0].max(2);
        let nv = counts[1].max(2);
        for ii in 0..nu {
            let fu = ii as f64 / (nu - 1) as f64;
            for jj in 0..nv {
                let fv = jj as f64 / (nv - 1) as f64;
                let uv = domain.evaluate(fu, fv);
                let mut value = self.function_value(uv);
                let mut best_uv = uv;
                if let Some((temp_uv, temp_value)) = self.check_bad_guess_stride(uv)?
                    && temp_value < value
                {
                    value = temp_value;
                    best_uv = temp_uv;
                }
                if value < *best {
                    *guess = best_uv;
                    *best = value;
                }
            }
        }
        Ok(())
    }

    /// `FindRoughGuess`: the control-net guess, or the better of it and a
    /// shrinking parameter grid when the net is suspicious.
    fn rough_guess(&self, best: &mut [usize; 2]) -> Result<Vector2, Fault> {
        *best = [0, 0];
        let (mut uv, suspicious) = self.polygon_guess(best);
        if !suspicious {
            return Ok(uv);
        }
        let mut grid = self
            .surface
            .counts()
            .map(|count| if count < 7 && count > 0 { count } else { 7 });
        let mut size = self.domain.size();
        let mut sub = self.domain;
        let mut grid_uv = Vector2::default();
        let mut best_value = BIG_DOUBLE;
        for _ in 0..3 {
            self.refine_patch_by_grid(sub, grid, &mut grid_uv, &mut best_value)?;
            if grid[0] < 2 || grid[1] < 2 {
                break;
            }
            let fraction_u = 1.0 / (grid[0] - 1) as f64 + 0.01;
            let fraction_v = 1.0 / (grid[1] - 1) as f64 + 0.01;
            size.x *= fraction_u;
            size.y *= fraction_v;
            let low = self.domain.clamp_point(grid_uv - size);
            let high = self.domain.clamp_point(grid_uv + size);
            let _ = sub.set_min_max(low, high);
            for g in &mut grid {
                if *g > 4 {
                    *g -= 1;
                }
            }
            for g in &mut grid {
                if *g > 2 {
                    *g -= 1;
                }
            }
        }
        let polygon_value = self.function_value(uv);
        let grid_value = self.function_value(grid_uv);
        if polygon_value >= BIG_DOUBLE && grid_value >= BIG_DOUBLE {
            return Err(Fault::Degenerate);
        }
        if grid_value < polygon_value {
            uv = grid_uv;
        }
        Ok(uv)
    }

    /// `RefineGuess`: up to three plane steps from `uv` while they improve.
    fn refine_guess(&mut self, uv: &mut Vector2) {
        self.check_guess_for_seams(uv);
        let mut this_uv = *uv;
        let mut next_uv = *uv;
        let mut next_value = BIG_DOUBLE;
        let u_extent = self.domain.u_interval().length();
        let v_extent = self.domain.v_interval().length();
        let (max_u, max_v) = (u_extent * 0.10, v_extent * 0.10);
        for _ in 0..3 {
            let this_value = next_value;
            next_value = self.function_value(next_uv);
            if next_value.abs() < EFF_ZERO_SQRT {
                *uv = next_uv;
                return;
            }
            if next_value >= this_value {
                break;
            }
            let Ok([point, su, sv]) = first_derivatives(self.surface, next_uv, [false, false])
            else {
                return;
            };
            let gap = point - self.test;
            let (f, g) = (gap.dot(su), gap.dot(sv));
            let (fu, fv) = (su.dot(su), su.dot(sv));
            let gu = fv;
            let gv = sv.dot(sv);
            let det = fu * gv - fv * gu;
            let mut u_shift = -(f * gv - fv * g);
            let mut v_shift = -(fu * g - f * gu);
            if is_zero(det) || is_zero(u_shift) {
                u_shift = 0.0;
            } else if det.abs() > u_shift.abs() * EFF_ZERO_SQRT {
                u_shift /= det;
            } else {
                u_shift = 0.0;
            }
            if is_zero(det) || is_zero(v_shift) {
                v_shift = 0.0;
            } else if det.abs() > v_shift.abs() * EFF_ZERO_SQRT {
                v_shift /= det;
            } else {
                v_shift = 0.0;
            }
            if u_shift.abs() < u_extent / 1000.0 && v_shift.abs() < v_extent / 1000.0 {
                *uv = next_uv;
                return;
            }
            if u_shift.abs() > max_u {
                u_shift = if u_shift > 0.0 { max_u } else { -max_u };
            }
            if v_shift.abs() > max_v {
                v_shift = if v_shift > 0.0 { max_v } else { -max_v };
            }
            this_uv = next_uv;
            next_uv = self
                .domain
                .clamp_point(Vector2::new(this_uv.x + u_shift, this_uv.y + v_shift));
            self.check_guess_for_seams(&mut next_uv);
        }
        let this_value = next_value;
        next_value = self.function_value(next_uv);
        *uv = if next_value < this_value {
            next_uv
        } else {
            this_uv
        };
    }

    /// `RefineGuessByGrid`: three shrinking grids around the net square of
    /// `best`, keeping the best point found.
    fn refine_guess_by_grid(&self, uv: &mut Vector2, best: [usize; 2]) -> Result<(), Fault> {
        let mut best_guess = *uv;
        let mut best_value = self.function_value(best_guess);
        let [nu, nv] = self.surface.counts();
        let low_u = best[0].saturating_sub(1);
        let high_u = if best[0] < nu - 1 {
            best[0] + 1
        } else {
            nu - 1
        };
        let low_v = best[1].saturating_sub(1);
        let high_v = if best[1] < nv - 1 {
            best[1] + 1
        } else {
            nv - 1
        };
        let mut low = Vector2::new(self.greville(0, low_u), self.greville(1, low_v));
        let mut high = Vector2::new(self.greville(0, high_u), self.greville(1, high_v));
        low = self.domain.clamp_point(low);
        high = self.domain.clamp_point(high);
        let mut sub = Extent2::new(low, high);
        const GRID: usize = 6;
        let mut value = BIG_DOUBLE;
        for pass in 0..3 {
            self.refine_patch_by_grid(sub, [GRID, GRID], uv, &mut value)?;
            if pass == 2 {
                break;
            }
            if value < best_value {
                best_value = value;
                best_guess = *uv;
            }
            let delta = Vector2::new(high.x - low.x, high.y - low.y);
            low = Vector2::new(
                uv.x - delta.x / (GRID - 1) as f64,
                uv.y - delta.y / (GRID - 1) as f64,
            );
            high = Vector2::new(
                uv.x + delta.x / (GRID - 1) as f64,
                uv.y + delta.y / (GRID - 1) as f64,
            );
            low = self.domain.clamp_point(low);
            high = self.domain.clamp_point(high);
            let _ = sub.set_min_max(low, high);
        }
        if value > best_value {
            *uv = best_guess;
        }
        Ok(())
    }

    /// `IsLocalMax`: on the normal through the test point, beyond the
    /// centre of curvature along an interior direction.
    fn is_local_max(&self, jet: &Jet, uv: Vector2, tolerance: f64) -> bool {
        let gap = jet.position - self.test;
        if gap.length() <= tolerance {
            return false;
        }
        self.beyond_curvature(jet, gap, uv, tolerance)
    }

    /// `IsLocalMin`, which SMLib writes as `IsLocalMax` except at the test
    /// point itself.
    fn is_local_min(&self, jet: &Jet, uv: Vector2, tolerance: f64) -> bool {
        let gap = jet.position - self.test;
        if gap.length() <= tolerance {
            return true;
        }
        self.beyond_curvature(jet, gap, uv, tolerance)
    }

    fn beyond_curvature(&self, jet: &Jet, gap: Vector3, uv: Vector2, tolerance: f64) -> bool {
        if gap.dot(jet.su).abs() > tolerance || gap.dot(jet.sv).abs() > tolerance {
            return false;
        }
        let noise = EFF_ZERO * 1000.0;
        let flags = self.domain.point_boundaries(uv, EFF_ZERO);
        if flags & (1 | 4) == 0 {
            let f1 = jet.su.dot(jet.su);
            let f2 = gap.dot(jet.suu);
            if f1.abs() > noise && f1 + f2 < -noise {
                return true;
            }
        }
        if flags & (2 | 8) == 0 {
            let f1 = jet.sv.dot(jet.sv);
            let f2 = gap.dot(jet.svv);
            if f1.abs() > noise && f1 + f2 < -noise {
                return true;
            }
        }
        false
    }

    /// `PushingBoundary`: on a side with the plane step pointing out.
    fn pushing_boundary(&self, jet: &Jet, uv: Vector2, tolerance: f64) -> bool {
        if !self.domain.point_on_boundary(uv, tolerance.abs()).0 {
            return false;
        }
        let gap = self.test - jet.position;
        let Ok(step) = vec_lin_comb_two_vectors(jet.su, jet.sv, gap) else {
            return false;
        };
        let (u, v) = (self.domain.u_interval(), self.domain.v_interval());
        (uv.x < u.min + EFF_ZERO && step.x < -tolerance)
            || (uv.x > u.max - EFF_ZERO && step.x > tolerance)
            || (uv.y < v.min + EFF_ZERO && step.y < -tolerance)
            || (uv.y > v.max - EFF_ZERO && step.y > tolerance)
    }

    /// The local solve inside `Solve`, as an `SmSolution`.
    fn local(&self, domain: Extent2, guess: Vector2) -> Solve<(bool, Solution)> {
        let point = local_point_solve(
            self.surface,
            domain,
            self.operation,
            self.test,
            guess,
            self.owner,
        )?;
        let mut parameters = [0.0; crate::solution::MAX_PARAMETERS];
        parameters[0] = point.uv.x;
        parameters[1] = point.uv.y;
        let mut solution = Solution {
            kind: SolutionType::SingleValue,
            variables: 2,
            start: SolutionEnd {
                value: point.value,
                parameters,
            },
            object_count: 1,
            ..Solution::default()
        };
        solution.objects[0] = Some(SURFACE_OBJECT);
        Ok((point.found, solution))
    }

    /// `CheckSeams`: on a closed seam, the solution again on the other side
    /// (or both other sides and the far corner on a torus).
    fn check_seams(&mut self, solution: Solution) -> Solve<SolutionArray> {
        let mut result = SolutionArray::new();
        result.add(solution);
        let mut uv = Vector2::new(solution.start.parameters[0], solution.start.parameters[1]);
        let tolerance = 0.001 * self.domain.size().length();
        let boundary = is_on_boundary(self.surface, self.domain, uv, tolerance);
        if !boundary.on {
            return Ok(result);
        }
        let [on_u, on_v] = boundary.sides;
        let mut other_solution = solution;
        let problem = Problem {
            surface: self.surface,
            domain: self.domain,
            test: self.test,
            operation: self.operation,
            tight: self.tight,
            loose: self.loose,
        };
        let mut this_eval = Stepper::new(&problem, uv);
        let mut other_eval = this_eval.clone();
        let norm = if self.operation == SolverOperation::Normalize {
            Norm::SumAbs
        } else {
            Norm::Gap
        };
        let objects = SurfaceObject {
            surface: self.surface,
            zone: self.zone,
        };
        for axis in 0..2 {
            let on = if axis == 0 { on_u } else { on_v };
            if !(self.closed_surface(axis) && on) {
                continue;
            }
            let interval = if axis == 0 {
                self.domain.u_interval()
            } else {
                self.domain.v_interval()
            };
            let value = if axis == 0 { uv.x } else { uv.y };
            let (near, far) = if value < interval.mid() {
                (interval.min, interval.max)
            } else {
                (interval.max, interval.min)
            };
            let at = |uv: Vector2, side: f64| {
                if axis == 0 {
                    Vector2::new(side, uv.y)
                } else {
                    Vector2::new(uv.x, side)
                }
            };
            if (value - near).abs() > self.tight {
                other_eval.set_uv(at(uv, near));
                if other_eval.error_norm(norm)? <= this_eval.error_norm(norm)? {
                    if axis == 0 {
                        uv.x = near;
                    } else {
                        uv.y = near;
                    }
                    let mut first = result.get_at(0)?;
                    first.start.parameters[axis] = near;
                    first.start.value = other_eval.gap_length()?;
                    result.set_at(0, first)?;
                }
            }
            let far_uv = at(uv, far);
            other_eval.set_uv(far_uv);
            if other_eval.error_norm(norm)? < this_eval.error_norm(norm)? + EFF_ZERO_SQRT {
                other_solution.start.parameters[0] = far_uv.x;
                other_solution.start.parameters[1] = far_uv.y;
                other_solution.start.value = other_eval.gap_length()?;
                other_solution.object_count = 1;
                other_solution.objects[0] = Some(SURFACE_OBJECT);
                result.add_sorted_solution(&objects, other_solution, SortKey::ByFirstTwoParameters);
            }
        }
        if self.closed_surface(0) && self.closed_surface(1) && on_u && on_v {
            let (u, v) = (self.domain.u_interval(), self.domain.v_interval());
            let far_u = if uv.x < u.mid() { u.max } else { u.min };
            let far_v = if uv.y < v.mid() { v.max } else { v.min };
            other_eval.set_uv(Vector2::new(far_u, far_v));
            if other_eval.error_norm(norm)? < this_eval.error_norm(norm)? + EFF_ZERO_SQRT {
                other_solution.start.parameters[0] = far_u;
                other_solution.start.parameters[1] = far_v;
                other_solution.start.value = other_eval.gap_length()?;
                result.add_sorted_solution(&objects, other_solution, SortKey::ByFirstTwoParameters);
            }
        }
        Ok(result)
    }
}

/// `SmSurface::GlobalPointSolve_1`: the nearest point (minimise), the feet
/// of normals (normalise) or a coincident point (intersect) of `surface`
/// within `domain` for `test`, with `tolerance` as the solver's distance
/// tolerance and `owner_tolerance` the zone tolerance of the owning face, if
/// any. Answers on a closed seam come with their copies across it.
pub fn global_point_solve(
    surface: &BSplineSurface,
    domain: Extent2,
    operation: SolverOperation,
    test: Vector3,
    tolerance: f64,
    owner_tolerance: Option<f64>,
) -> Result<SolutionArray, SolverRefusal> {
    global_point_solve_with(
        surface,
        domain,
        operation,
        test,
        tolerance,
        owner_tolerance,
        None,
    )
}

/// [`global_point_solve`] without an owning face and with SMLib's
/// `cpdOptTargetDistance`, which only a plane's `DropPointFast` reads (the
/// spline solver ignores it).
pub fn global_point_solve_to(
    surface: &BSplineSurface,
    domain: Extent2,
    operation: SolverOperation,
    test: Vector3,
    tolerance: f64,
    target: Option<f64>,
) -> Result<SolutionArray, SolverRefusal> {
    global_point_solve_with(surface, domain, operation, test, tolerance, None, target)
}

fn global_point_solve_with(
    surface: &BSplineSurface,
    domain: Extent2,
    operation: SolverOperation,
    test: Vector3,
    tolerance: f64,
    owner_tolerance: Option<f64>,
    target: Option<f64>,
) -> Result<SolutionArray, SolverRefusal> {
    if !matches!(
        operation,
        SolverOperation::Minimize | SolverOperation::Normalize | SolverOperation::Intersect
    ) {
        return Err(SolverRefusal::NotPorted(operation));
    }
    // SMLib first tries `DropPointFast`, which only a plane answers, and
    // keeps its answers inside the domain test.
    if surface.plane_frame().is_some()
        && matches!(
            operation,
            SolverOperation::Minimize | SolverOperation::Intersect
        )
    {
        let mut solutions = SolutionArray::new();
        if let Some((uv, distance)) =
            crate::plane::drop_point_fast(surface, domain, operation, test, tolerance, target)?
            && crate::surface_solve::is_in_domain(surface, domain, uv, tolerance)
        {
            let mut solution = Solution {
                kind: SolutionType::SingleValue,
                variables: 2,
                object_count: 1,
                ..Solution::default()
            };
            solution.objects[0] = Some(SURFACE_OBJECT);
            solution.start.value = distance;
            solution.start.parameters[0] = uv.x;
            solution.start.parameters[1] = uv.y;
            solutions.add(solution);
        }
        if !solutions.is_empty() {
            return Ok(solutions);
        }
    }
    let mut solver = Global {
        surface,
        domain,
        test,
        operation,
        tight: EFF_ZERO * (1.0 + test.max_dimension()),
        loose: tolerance,
        owner: owner_tolerance,
        zone: owner_tolerance.unwrap_or(ZONE_TOL_3D),
        closed: [None, None],
    };
    let mut best = [0usize; 2];
    let rough = solver.rough_guess(&mut best)?;
    let mut solutions = Vec::new();
    let mut found = false;
    let mut solution = Solution::default();
    // GetCreaseSubdomains keeps the whole domain.
    if domain.contains_point(rough, EFF_ZERO_SQRT) {
        let mut guess = rough;
        solver.refine_guess(&mut guess);
        solver.check_guess_for_seams(&mut guess);
        (found, solution) = solver.local(domain, guess)?;
        // SMLib's SmSurfaceEval at the solution, evaluated when first read.
        let mut uv = Vector2::new(solution.start.parameters[0], solution.start.parameters[1]);
        let loose = solver.loose;
        let try_again = match operation {
            SolverOperation::Minimize => {
                if found && solver.is_local_max(&solver.jet(uv)?, uv, loose) {
                    found = false;
                }
                !found
            }
            SolverOperation::Intersect => {
                if found && solution.start.value > loose {
                    found = false;
                }
                !found
                    && !solver.is_local_min(&solver.jet(uv)?, uv, loose)
                    && !solver.pushing_boundary(&solver.jet(uv)?, uv, loose)
            }
            _ => !found,
        };
        if try_again {
            let saved = guess;
            solver.refine_guess_by_grid(&mut guess, best)?;
            if !guess.closer_than(EFF_ZERO_SQRT, saved) {
                (found, solution) = solver.local(domain, guess)?;
                if found {
                    uv = Vector2::new(solution.start.parameters[0], solution.start.parameters[1]);
                }
            }
        }
        if found {
            found = match operation {
                SolverOperation::Minimize => !solver.is_local_max(&solver.jet(uv)?, uv, loose),
                SolverOperation::Intersect => solution.start.value <= loose,
                _ => !solver.pushing_boundary(&solver.jet(uv)?, uv, loose),
            };
        }
        if found {
            if solution.variables == 0 {
                solution.variables = 2;
            }
            solutions.push(solution);
        }
    }
    let mut best_value = BIG_DOUBLE;
    for candidate in &solutions {
        if candidate.start.value < best_value {
            solution = *candidate;
            best_value = solution.start.value;
            found = true;
        }
    }
    if found {
        return solver.check_seams(solution);
    }
    Ok(SolutionArray::new())
}
