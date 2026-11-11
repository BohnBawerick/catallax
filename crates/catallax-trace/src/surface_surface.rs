// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline surfaces and planes of
// `spline.rs`; statuses become `Fault`, solution containers small result
// structs, and the planes SMLib builds inside the solve become their
// canonical frames. The 3D signed directed operation of `LocalSurfaceSolve`
// refuses with `SolverRefusal::NotPorted`.
// Translated from SMLib SmSurface.cpp (SmSSPIntersectENFO,
// sm_IntersectSrfBndryCrv, SmSSINormalizeENFO,
// SmSurface::LocalSurfaceIntersect, SmSSNormalizeENFO,
// SmSurface::LocalSurfaceSolve); rust-brep has no port.
//! The local surface-surface solvers of step 7: a meeting of two surfaces
//! from a pair of guesses, held on a plane across the meeting curve, and the
//! nearest, farthest or normal point pair.

use crate::Fault;
use crate::curve_solve::PointSolveCurve;
use crate::curve_surface::{local_curve_intersect, second_derivatives};
use crate::extent::{Extent1, Extent2};
use crate::global_solver::{SolverOperation, SolverRefusal};
use crate::local_solve::{
    BoundaryHandler, EvalNFunctions, ExtentN, LocalSolveNd, TerminationReason, compute_jacobian,
};
use crate::matrix::Matrix;
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, PI, ZONE_TOL_3D, are_same};
use crate::spline::{BSplineSurface, PlaneFrame};
use crate::surface_eval::{evaluate_normal, evaluate_point, first_derivatives};
use crate::surface_solve::IsoCurve;
use crate::vector::{Vector2, Vector3};

/// A point pair on two surfaces (`SmSolution` with four variables): the
/// parameters on each and the distance between the points.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfacePair {
    /// Parameters on the first surface.
    pub uv: Vector2,
    /// Parameters on the other surface.
    pub other_uv: Vector2,
    /// The distance between the two points.
    pub value: f64,
}

/// The unit normal `sABC` and point `sH` of `SmPlane(origin, normal)`
/// evaluated at `(0, 0)`, as the intersection solve reads its plane.
fn plane_through(origin: Vector3, normal: Vector3) -> Result<(Vector3, Vector3), Fault> {
    let [_, y, z] = normal.make_unit_ortho_vectors(None)?;
    let [origin, x_axis, y_axis] = PlaneFrame::canonical(origin, y, z)?;
    let frame = PlaneFrame {
        origin,
        x_axis,
        y_axis,
        uv_scale: Vector2::new(1.0, 1.0),
        domain: Extent2::default(),
    };
    let point = frame.point_fast(Vector2::new(0.0, 0.0));
    let mut abc = (frame.x_axis * frame.uv_scale.x).cross(frame.y_axis * frame.uv_scale.y);
    abc.unitize()?;
    Ok((point, abc))
}

fn intervals(domain: Extent2, other_domain: Extent2) -> ExtentN {
    ExtentN {
        intervals: vec![
            domain.u_interval(),
            domain.v_interval(),
            other_domain.u_interval(),
            other_domain.v_interval(),
        ],
    }
}

fn unit(mut v: Vector3) -> Vector3 {
    let _ = v.unitize();
    v
}

/// `SmSSPIntersectENFO`: the gap between the surface points and their
/// distance sum from a plane, the plane fixed or, in steepest descent,
/// rebuilt through the midpoint across both normals each evaluation.
struct PlaneIntersection<'a> {
    surface: &'a BSplineSurface,
    other: &'a BSplineSurface,
    steepest_descent: bool,
    plane: (Vector3, Vector3),
}

impl EvalNFunctions for PlaneIntersection<'_> {
    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        jacobian: Option<&mut Matrix>,
        found: &mut bool,
    ) -> Result<(), Fault> {
        *found = false;
        let uv = Vector2::new(x[0], x[1]);
        let uv2 = Vector2::new(x[2], x[3]);
        let [pf, duf, dvf] = first_derivatives(self.surface, uv, [true, true])?;
        let [pg, dug, dvg] = first_derivatives(self.other, uv2, [true, true])?;
        let (h, abc) = if self.steepest_descent {
            let middle = (pf + pg).divided(2.0);
            let mut n1 = duf.cross(dvf);
            let mut n2 = dug.cross(dvg);
            if n1.length_squared() < EFF_ZERO {
                n1 = evaluate_normal(self.surface, uv, [true, true]);
            }
            if n2.length_squared() < EFF_ZERO {
                n2 = evaluate_normal(self.other, uv2, [true, true]);
            }
            let mut normal = n1.cross(n2);
            if normal.length_squared() < EFF_ZERO_SQ {
                return Err(Fault::Degenerate);
            }
            let _ = normal.unitize();
            plane_through(middle, normal)?
        } else {
            self.plane
        };
        let d = -abc.dot(h);
        if duf.length_squared() < EFF_ZERO_SQ
            || dvf.length_squared() < EFF_ZERO_SQ
            || dug.length_squared() < EFF_ZERO_SQ
            || dvg.length_squared() < EFF_ZERO_SQ
        {
            return Err(Fault::Degenerate);
        }
        let diff = pf - pg;
        f[3] = abc.dot(pf) + d + abc.dot(pg) + d;
        if diff.length_squared() < EFF_ZERO_SQ && f[3] < EFF_ZERO {
            // SMLib returns before writing the gap components.
            *found = true;
            return Ok(());
        }
        if let Some(j) = jacobian {
            for (row, (a, b, c, e)) in [
                (duf.x, dvf.x, dug.x, dvg.x),
                (duf.y, dvf.y, dug.y, dvg.y),
                (duf.z, dvf.z, dug.z, dvg.z),
            ]
            .into_iter()
            .enumerate()
            {
                j[row][0] = a;
                j[row][1] = b;
                j[row][2] = -c;
                j[row][3] = -e;
            }
            j[3][0] = abc.dot(duf);
            j[3][1] = abc.dot(dvf);
            j[3][2] = abc.dot(dug);
            j[3][3] = abc.dot(dvg);
        }
        f[0] = diff.x;
        f[1] = diff.y;
        f[2] = diff.z;
        Ok(())
    }
}

/// `SmSSINormalizeENFO`: the gap held normal to a fixed plane and both
/// points on it, with SMLib's finite-difference Jacobian over the natural
/// domains (its analytic one is computed and then replaced).
struct PlaneNormalPair<'a> {
    surface: &'a BSplineSurface,
    other: &'a BSplineSurface,
    normal: Vector3,
    d: f64,
}

impl EvalNFunctions for PlaneNormalPair<'_> {
    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        jacobian: Option<&mut Matrix>,
        found: &mut bool,
    ) -> Result<(), Fault> {
        *found = false;
        let uv = Vector2::new(x[0], x[1]);
        let uv2 = Vector2::new(x[2], x[3]);
        let [mut p1, du, dv, ..] = second_derivatives(self.surface, uv, [true, true])?;
        let [p2, du2, dv2, ..] = second_derivatives(self.other, uv2, [true, true])?;
        if du.length_squared() < EFF_ZERO_SQ
            || dv.length_squared() < EFF_ZERO_SQ
            || du2.length_squared() < EFF_ZERO_SQ
            || dv2.length_squared() < EFF_ZERO_SQ
        {
            return Err(Fault::Degenerate);
        }
        let normal = self.normal;
        let mut gap = p2 - p1;
        let tolerance = EFF_ZERO * 1000.0 * (1.0 + p1.max_dimension());
        if gap.length_squared() < tolerance {
            let n1 = evaluate_normal(self.surface, uv, [true, true]);
            let n2 = evaluate_normal(self.other, uv2, [true, true]);
            let average = if n1.dot(n2) < 0.0 {
                (n1 - n2).divided(2.0)
            } else {
                (n1 + n2).divided(2.0)
            };
            let shift = normal.cross(average).cross(normal);
            p1 = if n1.dot(gap) > 0.0 { p1 - shift } else { p1 + shift };
            gap = p2 - p1;
        }
        let mut nu = du.cross(dv);
        if nu.length_squared() < EFF_ZERO_SQ {
            nu = evaluate_normal(self.surface, uv, [true, true]);
            if nu.length_squared() < EFF_ZERO_SQ {
                return Err(Fault::Degenerate);
            }
        }
        let nu = unit(nu);
        let mut n2u = du2.cross(dv2);
        if n2u.length_squared() < EFF_ZERO_SQ {
            // SMLib evaluates the first surface at the second's parameters.
            n2u = evaluate_normal(self.surface, uv2, [true, true]);
            if n2u.length_squared() < EFF_ZERO_SQ {
                return Err(Fault::Degenerate);
            }
        }
        let n2u = unit(n2u);
        let v1 = nu.cross(normal);
        let v2 = n2u.cross(normal);
        f[0] = v1.dot(gap);
        f[1] = v2.dot(gap);
        f[2] = p1.dot(normal) + self.d;
        f[3] = p2.dot(normal) + self.d;
        if let Some(j) = jacobian {
            let domains = [
                self.surface.natural_domain().u_interval(),
                self.surface.natural_domain().v_interval(),
                self.other.natural_domain().u_interval(),
                self.other.natural_domain().v_interval(),
            ];
            let mut x = x.to_vec();
            let values = f.to_vec();
            compute_jacobian(self, &mut x, &values, j, Some(&domains))?;
        }
        Ok(())
    }
}

/// `SmSSNormalizeENFO` outside the directed operation: both unit tangent
/// pairs perpendicular to the gap.
struct NormalPair<'a> {
    surface: &'a BSplineSurface,
    other: &'a BSplineSurface,
}

impl EvalNFunctions for NormalPair<'_> {
    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        jacobian: Option<&mut Matrix>,
        found: &mut bool,
    ) -> Result<(), Fault> {
        *found = false;
        let uv1 = Vector2::new(x[0], x[1]);
        let uv2 = Vector2::new(x[2], x[3]);
        let [p1, s1u, s1v, s1uv, s1uu, s1vv] = second_derivatives(self.surface, uv1, [true, true])?;
        let [p2, s2u, s2v, s2uv, s2uu, s2vv] = second_derivatives(self.other, uv2, [true, true])?;
        if s1u.length_squared() < EFF_ZERO_SQ
            || s1v.length_squared() < EFF_ZERO_SQ
            || s2u.length_squared() < EFF_ZERO_SQ
            || s2v.length_squared() < EFF_ZERO_SQ
        {
            return Err(Fault::Degenerate);
        }
        let mut gap = p2 - p1;
        let tolerance = EFF_ZERO * 1000.0 * (1.0 + p1.max_dimension());
        if gap.length_squared() < tolerance {
            let normal = evaluate_normal(self.surface, uv1, [true, true]);
            gap = if normal.dot(gap) > 0.0 {
                p2 - (p1 - normal)
            } else {
                p2 - (p1 + normal)
            };
        }
        let s1u_unit = unit(s1u);
        let s1u_unit_u = s1u.unitized_derivative(s1uu);
        let s1u_unit_v = s1u.unitized_derivative(s1uv);
        let s1v_unit = unit(s1v);
        let s1v_unit_u = s1v.unitized_derivative(s1uv);
        let s1v_unit_v = s1v.unitized_derivative(s1vv);
        let s2u_unit = unit(s2u);
        let s2u_unit_u = s2u.unitized_derivative(s2uu);
        let s2u_unit_v = s2u.unitized_derivative(s2uv);
        let s2v_unit = unit(s2v);
        let s2v_unit_u = s2v.unitized_derivative(s2uv);
        let s2v_unit_v = s2v.unitized_derivative(s2vv);
        if let Some(j) = jacobian {
            j[0][0] = s1u_unit_u.dot(gap) - s1u_unit.dot(s1u);
            j[0][1] = s1u_unit_v.dot(gap) - s1u_unit.dot(s1v);
            j[0][2] = s1u_unit.dot(s2u);
            j[0][3] = s1u_unit.dot(s2v);
            j[1][0] = s1v_unit_u.dot(gap) - s1v_unit.dot(s1u);
            j[1][1] = s1v_unit_v.dot(gap) - s1v_unit.dot(s1v);
            j[1][2] = s1v_unit.dot(s2u);
            j[1][3] = s1v_unit.dot(s2v);
            j[2][0] = -s2u_unit.dot(s1u);
            j[2][1] = -s2u_unit.dot(s1v);
            j[2][2] = s2u_unit_u.dot(gap) + s2u_unit.dot(s2u);
            j[2][3] = s2u_unit_v.dot(gap) + s2u_unit.dot(s2v);
            j[3][0] = -s2v_unit.dot(s1u);
            j[3][1] = -s2v_unit.dot(s1v);
            j[3][2] = s2v_unit_u.dot(gap) + s2v_unit.dot(s2v);
            // SMLib's entry, `S2v_Unit . S2v` where `S2u` would match row 2.
            j[3][3] = s2v_unit_v.dot(gap) + s2v_unit.dot(s2v);
        }
        f[0] = s1u_unit.dot(gap);
        f[1] = s1v_unit.dot(gap);
        f[2] = s2u_unit.dot(gap);
        f[3] = s2v_unit.dot(gap);
        Ok(())
    }
}

/// `SmSurface::LocalSurfaceSolve` outside the 3D signed directed operation:
/// a point pair whose gap is normal to both surfaces, from the guesses.
/// Minimising and maximising also keep an answer that stopped on a bound.
pub fn local_surface_solve(
    surface: &BSplineSurface,
    domain: Extent2,
    other: &BSplineSurface,
    other_domain: Extent2,
    operation: SolverOperation,
    guess: Vector2,
    other_guess: Vector2,
) -> Result<Option<SurfacePair>, SolverRefusal> {
    if !matches!(
        operation,
        SolverOperation::Minimize | SolverOperation::Maximize | SolverOperation::Normalize
    ) {
        return Err(SolverRefusal::NotPorted(operation));
    }
    let mut solver = LocalSolveNd::new(
        NormalPair { surface, other },
        Some(intervals(domain, other_domain)),
        Some(vec![false; 4]),
    );
    solver.set_boundary_handler(BoundaryHandler::TotalBoundaryHits, 8);
    let root = solver.solve_it(&[guess.x, guess.y, other_guess.x, other_guess.y], EFF_ZERO_SQRT)?;
    if !root.found
        && (solver.termination_reason() != TerminationReason::OutOfBounds
            || operation == SolverOperation::Normalize)
    {
        return Ok(None);
    }
    let uv = Vector2::new(root.x[0], root.x[1]);
    let other_uv = Vector2::new(root.x[2], root.x[3]);
    let value = evaluate_point(surface, uv)?.distance(evaluate_point(other, other_uv)?);
    Ok(Some(SurfacePair {
        uv,
        other_uv,
        value,
    }))
}

/// What `sm_IntersectSrfBndryCrv` found: the parameters on the surface and
/// along the iso-curve, the deviation, and whether the curve meets the
/// surface tangentially.
#[derive(Clone, Copy, Debug, PartialEq)]
struct BoundaryMeeting {
    uv: Vector2,
    t: f64,
    tangent: bool,
}

/// `sm_IntersectSrfBndryCrv`: the meeting of `surface` with the iso-curve
/// of `curve_surface` at `parameter` along `axis` (0 fixes u), from the
/// guesses; none when the curve is degenerate or runs along the surface at
/// the guess.
#[allow(clippy::too_many_arguments)]
fn intersect_boundary_curve(
    surface: &BSplineSurface,
    domain: Extent2,
    curve_surface: &BSplineSurface,
    axis: usize,
    parameter: f64,
    tolerance: f64,
    guess_uv: Vector2,
    guess_t: f64,
) -> Result<Option<BoundaryMeeting>, SolverRefusal> {
    let curve = IsoCurve::new(curve_surface, axis, parameter, true);
    if curve.is_degenerate(EFF_ZERO) {
        return Ok(None);
    }
    let pv = curve.evaluate(guess_t, 1, true, true)?;
    let normal = evaluate_normal(surface, guess_uv, [false, false]);
    if normal.is_perpendicular_to(pv[1], 0.05) {
        return Ok(None);
    }
    let interval: Extent1 = curve.natural_interval();
    let Some(found) = local_curve_intersect(
        surface,
        domain,
        &curve,
        interval,
        tolerance,
        (guess_t, guess_uv),
        ZONE_TOL_3D,
    )?
    else {
        return Ok(None);
    };
    let pv = curve.evaluate(found.t, 1, true, true)?;
    let normal = evaluate_normal(surface, found.uv, [true, true]);
    let angle = pv[1].angle_between(normal)?;
    Ok(Some(BoundaryMeeting {
        uv: found.uv,
        t: found.t,
        tangent: (angle - PI / 2.0).abs() < 0.01,
    }))
}

fn on_side(value: f64, interval: Extent1) -> bool {
    are_same(value, interval.min) || are_same(value, interval.max)
}

/// `SmSurface::LocalSurfaceIntersect`: a meeting point of `surface` and
/// `other` from a pair of guesses within `tolerance`. Guesses on a side
/// first meet that side's iso-curve with the other surface; the Newton
/// solve holds the pair on a plane through them across their normals (or
/// across `plane_normal`, fixed); answers near tangency are refined by the
/// normal solves; `previous` keeps an interior answer close to the last
/// point over a side answer far from it.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn local_surface_intersect(
    surface: &BSplineSurface,
    domain: Extent2,
    other: &BSplineSurface,
    other_domain: Extent2,
    tolerance: f64,
    guess: Vector2,
    other_guess: Vector2,
    plane_normal: Option<Vector3>,
    previous: Option<Vector2>,
    boundary_testing: bool,
) -> Result<Option<SurfacePair>, SolverRefusal> {
    let mut found = false;
    let mut tried_boundary = false;
    let mut this_on_boundary = false;
    let mut other_on_boundary = false;
    let mut boundary: Option<(Vector2, Vector2)> = None;
    if boundary_testing {
        let (u, v) = (domain.u_interval(), domain.v_interval());
        let (ou, ov) = (other_domain.u_interval(), other_domain.v_interval());
        for (axis, value, t) in [(0, guess.x, guess.y), (1, guess.y, guess.x)] {
            if found || !on_side(value, if axis == 0 { u } else { v }) {
                continue;
            }
            tried_boundary = true;
            this_on_boundary = true;
            if let Some(meet) = intersect_boundary_curve(
                other,
                other_domain,
                surface,
                axis,
                value,
                tolerance,
                other_guess,
                t,
            )? {
                found = true;
                let uv = if axis == 0 {
                    Vector2::new(value, meet.t)
                } else {
                    Vector2::new(meet.t, value)
                };
                boundary = Some((uv, meet.uv));
            }
        }
        for (axis, value, t) in [(0, other_guess.x, other_guess.y), (1, other_guess.y, other_guess.x)] {
            if found || !on_side(value, if axis == 0 { ou } else { ov }) {
                continue;
            }
            other_on_boundary = true;
            tried_boundary = true;
            if let Some(meet) =
                intersect_boundary_curve(surface, domain, other, axis, value, tolerance, guess, t)?
            {
                found = true;
                let other_uv = if axis == 0 {
                    Vector2::new(value, meet.t)
                } else {
                    Vector2::new(meet.t, value)
                };
                boundary = Some((meet.uv, other_uv));
            }
        }
    }
    let p1 = evaluate_point(surface, guess)?;
    let p2 = evaluate_point(other, other_guess)?;
    let origin = if this_on_boundary {
        p1
    } else if other_on_boundary {
        p2
    } else {
        (p1 + p2).divided(2.0)
    };
    let steepest_descent = plane_normal.is_none();
    let mut normal = match plane_normal {
        Some(normal) => normal,
        None => {
            let n1 = evaluate_normal(surface, guess, [true, true]);
            let n2 = evaluate_normal(other, other_guess, [true, true]);
            let normal = n1.cross(n2);
            if normal.length_squared() < EFF_ZERO_SQ {
                return Ok(None);
            }
            normal
        }
    };
    let _ = normal.unitize();
    let plane = plane_through(origin, normal)?;
    let mut solver = LocalSolveNd::new(
        PlaneIntersection {
            surface,
            other,
            steepest_descent,
            plane,
        },
        Some(intervals(domain, other_domain)),
        Some(vec![false; 4]),
    );
    solver.set_desired_accuracy(EFF_ZERO_SQ);
    solver.set_boundary_handler(BoundaryHandler::TotalBoundaryHits, 3);
    solver.set_maximum_iterations(20);
    let root = solver.solve_it(
        &[guess.x, guess.y, other_guess.x, other_guess.y],
        tolerance * 4.0,
    )?;
    found = root.found;
    let mut x = [root.x[0], root.x[1], root.x[2], root.x[3]];
    let mut uv = Vector2::new(x[0], x[1]);
    let mut other_uv = Vector2::new(x[2], x[3]);
    if boundary_testing
        && found
        && !tried_boundary
        && (domain.point_on_boundary(uv, EFF_ZERO).0
            || other_domain.point_on_boundary(other_uv, EFF_ZERO).0)
    {
        return local_surface_intersect(
            surface,
            domain,
            other,
            other_domain,
            tolerance,
            uv,
            other_uv,
            plane_normal,
            previous,
            boundary_testing,
        );
    }
    if let Some((boundary_uv, boundary_other)) = boundary {
        let mut use_boundary = true;
        if found && let Some(previous) = previous {
            let interior = previous.distance(Vector2::new(x[0], x[1]));
            let side = previous.distance(boundary_uv);
            if interior > 100.0 * side {
                use_boundary = false;
            }
        }
        if use_boundary {
            found = true;
            x = [boundary_uv.x, boundary_uv.y, boundary_other.x, boundary_other.y];
            uv = boundary_uv;
            other_uv = boundary_other;
        }
    }
    if !found && solver.termination_reason() == TerminationReason::OutOfBounds && boundary_testing {
        let (u, v) = (domain.u_interval(), domain.v_interval());
        let (ou, ov) = (other_domain.u_interval(), other_domain.v_interval());
        if on_side(x[0], u)
            && let Some(meet) = intersect_boundary_curve(
                other,
                other_domain,
                surface,
                0,
                x[0],
                tolerance,
                Vector2::new(x[2], x[3]),
                x[1],
            )?
        {
            found = true;
            x = [x[0], meet.t, meet.uv.x, meet.uv.y];
        }
        // SMLib's precedence: the max sides are tried even after a find.
        if (!found && are_same(x[1], v.min)) || are_same(x[1], v.max) {
            if let Some(meet) = intersect_boundary_curve(
                other,
                other_domain,
                surface,
                1,
                x[1],
                tolerance,
                Vector2::new(x[2], x[3]),
                x[0],
            )? {
                found = true;
                x = [meet.t, x[1], meet.uv.x, meet.uv.y];
            } else {
                found = false;
            }
        }
        if (!found && are_same(x[2], ou.min)) || are_same(x[2], ou.max) {
            if let Some(meet) = intersect_boundary_curve(
                surface,
                domain,
                other,
                0,
                x[2],
                tolerance,
                Vector2::new(x[0], x[1]),
                x[3],
            )? {
                found = true;
                x = [meet.uv.x, meet.uv.y, x[2], meet.t];
            } else {
                found = false;
            }
        }
        if (!found && are_same(x[3], ov.min)) || are_same(x[3], ov.max) {
            if let Some(meet) = intersect_boundary_curve(
                surface,
                domain,
                other,
                1,
                x[3],
                tolerance,
                Vector2::new(x[0], x[1]),
                x[2],
            )? {
                found = true;
                x = [meet.uv.x, meet.uv.y, meet.t, x[3]];
            } else {
                found = false;
            }
        }
        uv = Vector2::new(x[0], x[1]);
        other_uv = Vector2::new(x[2], x[3]);
    }
    let near_parallel = |angle: f64| angle < 3.0 * PI / 180.0 || angle > PI - 3.0 * PI / 180.0;
    if found && steepest_descent {
        let n1 = evaluate_normal(surface, uv, [true, true]);
        let n2 = evaluate_normal(other, other_uv, [true, true]);
        if near_parallel(n1.angle_between(n2)?) {
            let refined = local_surface_solve(
                surface,
                domain,
                other,
                other_domain,
                SolverOperation::Minimize,
                uv,
                other_uv,
            )?;
            let distance = evaluate_point(surface, uv)?.distance(evaluate_point(other, other_uv)?);
            if let Some(refined) = refined
                && refined.value < distance
            {
                uv = refined.uv;
                other_uv = refined.other_uv;
                x = [uv.x, uv.y, other_uv.x, other_uv.y];
            }
        }
    }
    if !found {
        return Ok(None);
    }
    let p1 = evaluate_point(surface, uv)?;
    let p2 = evaluate_point(other, other_uv)?;
    let mut distance = p1.distance(p2);
    if distance > EFF_ZERO {
        let n1 = evaluate_normal(surface, uv, [true, true]);
        let n2 = evaluate_normal(other, other_uv, [true, true]);
        if near_parallel(n1.angle_between(n2)?) {
            if plane_normal.is_none() {
                let new_normal = n1.cross(n2);
                if new_normal.length_squared() > EFF_ZERO_SQ {
                    normal = new_normal;
                }
            }
            let _ = normal.unitize();
            let d = -p1.dot(normal);
            let mut refine = LocalSolveNd::new(
                PlaneNormalPair {
                    surface,
                    other,
                    normal,
                    d,
                },
                Some(intervals(domain, other_domain)),
                Some(vec![false; 4]),
            );
            refine.set_boundary_handler(BoundaryHandler::TotalBoundaryHits, 8);
            let root = refine.solve_it(&[uv.x, uv.y, other_uv.x, other_uv.y], EFF_ZERO_SQRT)?;
            if root.found {
                let a = Vector2::new(root.x[0], root.x[1]);
                let b = Vector2::new(root.x[2], root.x[3]);
                let d2 = evaluate_point(surface, a)?.distance(evaluate_point(other, b)?);
                if d2 < distance {
                    distance = d2;
                    uv = a;
                    other_uv = b;
                    x = [a.x, a.y, b.x, b.y];
                }
            }
            if distance != 0.0 {
                // SMLib stores this solve's flag in its found flag, which the
                // distance test below no longer reads.
                if let Some(normal_pair) = local_surface_solve(
                    surface,
                    domain,
                    other,
                    other_domain,
                    SolverOperation::Normalize,
                    uv,
                    other_uv,
                )? {
                    let d2 = evaluate_point(surface, normal_pair.uv)?
                        .distance(evaluate_point(other, normal_pair.other_uv)?);
                    if d2 < distance {
                        distance = d2;
                        x = [
                            normal_pair.uv.x,
                            normal_pair.uv.y,
                            normal_pair.other_uv.x,
                            normal_pair.other_uv.y,
                        ];
                    }
                }
            }
        }
    }
    if distance < tolerance {
        return Ok(Some(SurfacePair {
            uv: Vector2::new(x[0], x[1]),
            other_uv: Vector2::new(x[2], x[3]),
            value: distance,
        }));
    }
    Ok(None)
}
