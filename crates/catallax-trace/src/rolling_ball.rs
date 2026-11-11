// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: the solver borrows two checked B-spline surfaces
// instead of SMLib's topology and fillet geometry, statuses become `Fault`,
// and SMLib's function-static first-call step is chosen per call. Translated
// directly from SMLib SmFilletSolver.cpp (`ComputeSurfaceValues`,
// `ComputePointValues`, `ComputeUVVectors`, `StepSolve`,
// `SmFilletSphereSolveENFO::Evaluate`) and SmFilletStandardSolver.cpp
// (`SmConstantRadiusFS`: the surface constructor, `SetupOffsetValues`,
// `LoadJacobian`, `GetJacobianSize`); rust-brep has no port.
//! The rolling-ball equations of a constant-radius fillet between two
//! surfaces.
//!
//! A ball of radius `r1` against surface 1 and `r2` against surface 2 touches
//! both where the two offset surfaces meet. The solver finds such points: it
//! evaluates the offset surfaces' intersection point and tangent
//! ([`RollingBall::compute_point_values`]), and steps along it by solving four
//! equations in the four surface parameters, three for the offset points to
//! coincide and one for the point to lie on the plane normal to the curve
//! ([`RollingBall::step_solve`]).
//!
//! A step that the parameter domains clip goes on to SMLib's
//! `LocalCurveIntersect` on an offset surface, which needs the surface solvers
//! generalised over offset surfaces; the port refuses it with
//! [`StepRefusal::ClippedStep`]. SMLib's fillet topology (edgeuses, fillet
//! geometry for roll-over) is outside this port: the initial values of the
//! second surface that `LoadInitialValues` takes from it are the caller's.
//!
//! Everything returned is an unproved candidate.

use crate::{
    Fault,
    curve_surface::drop_vector,
    extent::{Extent1, Extent2},
    global_solver::SolverRefusal,
    local_solve::{BoundaryHandler, EvalNFunctions, ExtentN, LocalSolveNd},
    matrix::Matrix,
    scalar::{EFF_ZERO, EFF_ZERO_SQ, PI},
    spline::BSplineSurface,
    surface_offset::OffsetSurface,
    surface_trace::{PointKind, TsectPoint},
    vector::{Vector2, Vector3},
};

/// Why a step produced no value.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StepRefusal {
    /// SMLib's `SM_ERR` or `SM_ERR_INVALID_INPUT`.
    Failed(Fault),
    /// The step was clipped by a parameter domain and SMLib re-solves it with
    /// `LocalCurveIntersect` on an offset surface, which is not ported.
    ClippedStep,
}

impl From<Fault> for StepRefusal {
    fn from(fault: Fault) -> Self {
        Self::Failed(fault)
    }
}

/// The outcome of [`RollingBall::step_solve`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepOutcome {
    /// `rbFoundGoodPoint`.
    pub found: bool,
    /// `rbClipped`: always false, a clipped step is refused.
    pub clipped: bool,
    /// `rbBoundaryHit`: the current point is on a boundary and no step is
    /// possible.
    pub boundary_hit: bool,
    /// `rdStep` after the call.
    pub step: f64,
    /// `rNextPoint`, whose surface parameters are set when `found`.
    pub next: TsectPoint,
}

/// A constant-radius rolling-ball solver (`SmConstantRadiusFS`) between two
/// surfaces.
#[derive(Clone, Copy, Debug)]
pub struct RollingBall<'a> {
    offsets: [OffsetSurface<'a>; 2],
    approximation_tolerance: f64,
    conversion_tolerance: f64,
}

impl<'a> RollingBall<'a> {
    /// `SmConstantRadiusFS` from two surfaces: balls of `radii[i]` against
    /// surface `i`, on the side its `outward[i]` orientation names. SMLib's
    /// offset-surface extension is not carried.
    pub fn new(
        surfaces: [&'a BSplineSurface; 2],
        radii: [f64; 2],
        outward: [bool; 2],
        approximation_tolerance: f64,
    ) -> Self {
        let orientation = |outward: bool| if outward { 1.0 } else { -1.0 };
        Self {
            offsets: [
                OffsetSurface::new(surfaces[0], radii[0] * orientation(outward[0])),
                OffsetSurface::new(surfaces[1], radii[1] * orientation(outward[1])),
            ],
            approximation_tolerance,
            // `SmFilletSolver`'s `m_dConversionTol`.
            conversion_tolerance: EFF_ZERO / 10.0,
        }
    }

    /// `GetSurface`: the offset surface of rail `index`.
    pub fn surface(&self, index: usize) -> &OffsetSurface<'a> {
        &self.offsets[index]
    }

    /// `SmConstantRadiusFS::GetJacobianSize`.
    pub fn jacobian_size(&self) -> usize {
        4
    }

    /// `SmSurface::Evaluate1stDerivatives` on an offset surface: the point,
    /// `Su` and `Sv`.
    fn first_derivatives(&self, index: usize, uv: Vector2) -> Result<[Vector3; 3], Fault> {
        let cells = self.offsets[index].evaluate(uv, 1, 1, [true, true], true, true)?;
        Ok([cells[0], cells[2], cells[1]])
    }

    /// `SmFilletSolver::ComputeSurfaceValues`: the offset points, derivatives
    /// and unit normals at the two parameters, the curve point and the angle
    /// between the normals folded into `[0, pi / 2]`.
    pub fn compute_surface_values(
        &self,
        uv: [Vector2; 2],
        point: &mut TsectPoint,
    ) -> Result<(), Fault> {
        for (surface, &uv) in uv.iter().enumerate() {
            let [position, du, dv] = self.first_derivatives(surface, uv)?;
            let cells = &mut point.surface[surface];
            cells[0] = position;
            cells[3] = du;
            cells[1] = dv;
            let zero = Vector3::new(0.0, 0.0, 0.0);
            cells[6] = zero;
            cells[4] = zero;
            cells[2] = zero;
            let mut normal = du.cross(dv);
            if normal.length_squared() < EFF_ZERO_SQ {
                normal = self.offsets[surface].evaluate_normal(uv, [true, true]);
                if normal.length_squared() < EFF_ZERO_SQ {
                    return Err(Fault::Degenerate);
                }
            }
            let _ = normal.unitize();
            cells[8] = normal;
            point.uv[surface] = uv;
        }
        point.curve_point = point.surface[0][0];
        let mut angle = point.normal(0).angle_between(point.normal(1))?;
        if angle > PI / 2.0 {
            angle = PI - angle;
        }
        point.tangent_plane_angle = angle;
        Ok(())
    }

    /// `SmFilletSolver::ComputeUVVectors`: each surface's parameter tangent
    /// in the direction of the curve.
    pub fn compute_uv_vectors(&self, point: &mut TsectPoint) -> Result<(), Fault> {
        for surface in 0..2 {
            point.uv_deriv[surface] =
                drop_vector(point.du(surface), point.dv(surface), point.curve_tangent).map_err(
                    |refusal| match refusal {
                        SolverRefusal::Failed(fault) => fault,
                        SolverRefusal::NotPorted(_) => Fault::Degenerate,
                    },
                )?;
        }
        Ok(())
    }

    /// `SmFilletSolver::ComputePointValues` for a constant radius: the
    /// surface values at the two parameters, the curve tangent along the cross
    /// product of the normals in the trace direction, and the parameter
    /// tangents. Where the normals are parallel the tangent is the previous
    /// point's, or found by stepping both ways.
    ///
    /// SMLib sizes the finite-difference step from a function-static
    /// initialised on its first call; the port sizes it from this point's
    /// kind on every call.
    pub fn compute_point_values(
        &self,
        uv: [Vector2; 2],
        domains: [Extent2; 2],
        trace_direction: f64,
        point: &mut TsectPoint,
        previous: Option<&TsectPoint>,
    ) -> Result<(), StepRefusal> {
        self.compute_surface_values(uv, point)?;
        point.curve_tangent = point.normal(0).cross(point.normal(1)) * trace_direction;
        if point.curve_tangent.length_squared() < 10000.0 * EFF_ZERO_SQ {
            point.kind = PointKind::TangentPoint;
            match previous {
                Some(previous) => point.curve_tangent = previous.curve_tangent,
                None => return Ok(()),
            }
        }
        point.curve_tangent.unitize()?;
        if point.kind == PointKind::TangentPoint {
            self.compute_uv_vectors(point)?;
            self.tangent_by_stepping(domains, trace_direction, point)?;
        }
        point.kind = PointKind::Crossing;
        let _ = self.compute_uv_vectors(point);
        Ok(())
    }

    /// The finite-difference branch of `ComputePointValues` at a tangent
    /// point: steps forward and back and takes the direction between the
    /// midpoints of the two offset points.
    fn tangent_by_stepping(
        &self,
        domains: [Extent2; 2],
        trace_direction: f64,
        point: &mut TsectPoint,
    ) -> Result<(), StepRefusal> {
        let mut step = 1.0e-3;
        let forward = self.step_solve(point, step, domains[0], domains[1])?;
        step = forward.step;
        let mut forward = forward;
        if !forward.found && step < 1.0e-4 {
            // SMLib retries the same step, `dBigStep` only naming the result.
            let retry = self.step_solve(point, step, domains[0], domains[1])?;
            if retry.found {
                step = 1.0e-3;
            }
            forward = retry;
        }
        let evaluate = |index: usize, uv: Vector2| self.offsets[index].evaluate_point(uv);
        let forward_points = if forward.found {
            Some((
                evaluate(0, forward.next.uv[0])?,
                evaluate(1, forward.next.uv[1])?,
            ))
        } else {
            None
        };
        step = -step;
        let backward = self.step_solve(point, step, domains[0], domains[1])?;
        step = backward.step;
        let backward_points = if backward.found {
            Some((
                evaluate(0, backward.next.uv[0])?,
                evaluate(1, backward.next.uv[1])?,
            ))
        } else {
            None
        };
        let middle = |(a, b): (Vector3, Vector3)| (a + b) / 2.0;
        let (ahead, behind) = match (forward_points, backward_points) {
            (Some(f), Some(b)) => {
                step *= 2.0;
                (middle(f), middle(b))
            }
            (Some(f), None) => (middle(f), point.curve_point),
            (None, Some(b)) => (point.curve_point, middle(b)),
            (None, None) => return Err(Fault::Degenerate.into()),
        };
        let step = step.abs();
        point.curve_tangent = (ahead - behind) / step;
        point.curve_tangent.unitize()?;
        // At a tangent point the direction is only known up to sign.
        let reference = if forward.found { &forward } else { &backward };
        let normal_0 = self.offsets[0].evaluate_normal(reference.next.uv[0], [true, true]);
        let normal_1 = self.offsets[1].evaluate_normal(reference.next.uv[1], [true, true]);
        let direction = normal_0.cross(normal_1);
        if direction.length_squared() < EFF_ZERO_SQ {
            return Err(Fault::Degenerate.into());
        }
        let direction = direction * trace_direction;
        if direction.dot(point.curve_tangent) < 0.0 {
            point.curve_tangent = -point.curve_tangent;
        }
        Ok(())
    }

    /// `SmFilletSolver::StepSolve`: from `current`, steps `step` along its
    /// parameter tangents and solves for the next point of the intersection
    /// of the offset surfaces on the plane normal to the curve there. The
    /// step is shortened to stay in the domains; a step that is clipped is
    /// refused, see the module notes.
    pub fn step_solve(
        &self,
        current: &TsectPoint,
        step: f64,
        domain_1: Extent2,
        domain_2: Extent2,
    ) -> Result<StepOutcome, StepRefusal> {
        let domains = [domain_1, domain_2];
        let mut outcome = StepOutcome {
            found: false,
            clipped: false,
            boundary_hit: false,
            step,
            next: TsectPoint::default(),
        };
        let mut old_step = step;
        let mut clipped = false;
        for (surface, domain) in domains.iter().enumerate() {
            let clipped_step =
                domain.clip_line(current.uv[surface], current.uv_deriv[surface], old_step);
            if clipped_step.abs() < old_step.abs() - EFF_ZERO * old_step.abs()
                || clipped_step * old_step < 0.0
            {
                old_step = clipped_step;
                clipped = true;
            }
        }
        let step_size = old_step;
        outcome.clipped = clipped;
        if step_size.abs() <= EFF_ZERO || step_size * step < 0.0 {
            outcome.boundary_hit = true;
            return Ok(outcome);
        }
        if clipped {
            return Err(StepRefusal::ClippedStep);
        }
        let start = [0, 1].map(|surface| {
            let target = current.uv[surface] + step_size * current.uv_deriv[surface];
            domains[surface].clamp_point(target)
        });
        let intervals = ExtentN {
            intervals: vec![
                Extent1::new(domain_1.min.x, domain_1.max.x),
                Extent1::new(domain_1.min.y, domain_1.max.y),
                Extent1::new(domain_2.min.x, domain_2.max.x),
                Extent1::new(domain_2.min.y, domain_2.max.y),
            ],
        };
        let guess = [start[0].x, start[0].y, start[1].x, start[1].y];
        let plane_abc = current.curve_tangent;
        let plane_origin = current.curve_point + step * plane_abc;
        let plane_d = -plane_origin.dot(plane_abc);
        let mut solver = LocalSolveNd::new(
            SphereSolve {
                ball: self,
                plane_abc,
                plane_d,
                plane_adjusts: false,
            },
            Some(intervals),
            Some(vec![false; 4]),
        );
        solver.set_desired_accuracy(EFF_ZERO_SQ);
        solver.set_boundary_handler(BoundaryHandler::TotalBoundaryHits, 8);
        solver.set_maximum_iterations(30);
        let mut root = solver.solve_it(&guess, self.approximation_tolerance)?;
        if !root.found {
            solver.evaluator.plane_adjusts = true;
            root = solver.solve_it(&guess, self.approximation_tolerance)?;
        }
        if root.found {
            outcome.found = true;
            outcome.next.uv[0] = Vector2::new(root.x[0], root.x[1]);
            outcome.next.uv[1] = Vector2::new(root.x[2], root.x[3]);
        }
        Ok(outcome)
    }
}

/// `SmFilletSphereSolveENFO` with `SmConstantRadiusFS::LoadJacobian`: the
/// plane equation and the three equations of the offset points coinciding.
struct SphereSolve<'a, 'b> {
    ball: &'b RollingBall<'a>,
    plane_abc: Vector3,
    plane_d: f64,
    plane_adjusts: bool,
}

impl EvalNFunctions for SphereSolve<'_, '_> {
    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        jacobian: Option<&mut Matrix>,
        found: &mut bool,
    ) -> Result<(), Fault> {
        *found = false;
        f.fill(0.0);
        let mut jacobian = jacobian;
        if let Some(jacobian) = jacobian.as_deref_mut() {
            for row in 0..jacobian.rows() {
                for column in 0..jacobian.columns() {
                    jacobian[row][column] = 0.0;
                }
            }
        }
        // `LoadJacobian`: the offset points must coincide.
        let uv_1 = Vector2::new(x[0], x[1]);
        let [point_1, du_1, dv_1] = self.ball.first_derivatives(0, uv_1)?;
        let uv_2 = Vector2::new(x[2], x[3]);
        let [point_2, du_2, dv_2] = self.ball.first_derivatives(1, uv_2)?;
        if du_1.length_squared() < EFF_ZERO_SQ
            || dv_1.length_squared() < EFF_ZERO_SQ
            || du_2.length_squared() < EFF_ZERO_SQ
            || dv_2.length_squared() < EFF_ZERO_SQ
        {
            return Err(Fault::Degenerate);
        }
        let difference = point_1 - point_2;
        if let Some(jacobian) = jacobian.as_deref_mut() {
            let rows = [
                (du_1.x, dv_1.x, -du_2.x, -dv_2.x),
                (du_1.y, dv_1.y, -du_2.y, -dv_2.y),
                (du_1.z, dv_1.z, -du_2.z, -dv_2.z),
            ];
            for (i, (a, b, c, d)) in rows.into_iter().enumerate() {
                jacobian[1 + i][0] = a;
                jacobian[1 + i][1] = b;
                jacobian[1 + i][2] = c;
                jacobian[1 + i][3] = d;
            }
        }
        f[1] = difference.x;
        f[2] = difference.y;
        f[3] = difference.z;
        let scaled_tolerance = self.ball.conversion_tolerance * (1.0 + point_1.max_dimension());
        let coincident = difference.x.abs() < scaled_tolerance
            && difference.y.abs() < scaled_tolerance
            && difference.z.abs() < scaled_tolerance;
        // The plane equation, from surface 1's first derivatives again.
        let [plane_point, du, dv] = self.ball.first_derivatives(0, uv_1)?;
        if self.plane_adjusts {
            let mut normal_1 = du.cross(dv);
            if normal_1.length_squared() < EFF_ZERO_SQ {
                normal_1 = self.ball.offsets[0].evaluate_normal(uv_1, [true, true]);
            }
            normal_1.unitize()?;
            let [point_g, du_g, dv_g] = self.ball.first_derivatives(1, uv_2)?;
            let mut normal_2 = du_g.cross(dv_g);
            if normal_2.length_squared() < EFF_ZERO_SQ {
                normal_2 = self.ball.offsets[1].evaluate_normal(uv_2, [true, true]);
            }
            normal_2.unitize()?;
            let origin = (plane_point + point_g) / 2.0;
            let mut abc = normal_1.cross(normal_2);
            if abc.dot(self.plane_abc) < 0.0 {
                abc = -abc;
            }
            self.plane_abc = abc;
            self.plane_d = -abc.dot(origin);
        }
        f[0] = self.plane_abc.dot(plane_point) + self.plane_d;
        if let Some(jacobian) = jacobian {
            jacobian[0][0] = self.plane_abc.dot(du);
            jacobian[0][1] = self.plane_abc.dot(dv);
        }
        if coincident {
            let scaled = self.ball.conversion_tolerance * (1.0 + plane_point.max_dimension());
            if f[0].abs() < scaled {
                *found = true;
            }
        }
        Ok(())
    }
}
