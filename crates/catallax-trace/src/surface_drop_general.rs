// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline curves, surfaces and
// planes of `spline.rs` without trimming faces or graphics; statuses become
// a small status enum and `Fault`, the drop-point stack becomes a vector of
// tuples, and output arrays become a returned struct. SMLib's tri-state
// booleans are kept where an uninitialised value is read. The offset
// surface is evaluated to first order only (`SmOffsetSurface::
// EvaluatePlainOrSimple` with its pole extrapolation); the places where
// SMLib reads a value a failed evaluation left undefined refuse with
// `SolverRefusal::NotPorted`.
// Translated from SMLib SmSurfaceDropLib.cpp (SmSurface::DropCurve and the
// SmDropPt helpers it runs: SetUninitialized, IsGoodTangentDrop,
// IsMovingAlongBoundary, TestSpanAccuracy, GetNextPoint, SetProperties,
// CheckSeamJump, FixBoundaryTangents, sm_PushPoint, sm_PopStack,
// sm_StackLook, sm_FixBoundaryCurves, sm_RemoveDuplicateDropPointSolutions,
// sm_FlushCrv, sm_AddToCurve), SmCrvOnSurf.cpp (sm_EvalDistAndDev,
// SurfaceCurveMaxDistanceBetween, EvaluateSurfaceNormal), SmSurface.cpp
// (smsurf_IsPointOnNormal, DropVectors, EvaluateDirectionalDerivs),
// SmOffsetSurface.cpp (EvaluatePlainOrSimple to first order,
// sm_FindNonsingularStepoff), SmGeomUtility.cpp
// (smgu_IsDogLeg for UV points, smgu_DirectionalDerivs) and SmTol.cpp
// (GetScaledZero for a curve); rust-brep has no port.
//! SMLib's general curve drop (`SmSurface::DropCurve`): the parameter
//! curves of a 3D curve walked onto a surface by dropping points and their
//! tangents, as piecewise cubic Hermite spans that are split until three
//! interior samples lie on the surface normals within the tolerance.
//!
//! Everything returned is an unproved candidate.

use crate::Fault;
use crate::curve::approximate_length;
use crate::curve_cache::CurveCache;
use crate::curve_global::drop_point as curve_drop_point;
use crate::curve_surface::{
    SINGULAR_U_MAX, SINGULAR_U_MIN, SINGULAR_V_MAX, SINGULAR_V_MIN, Singular,
    approx_derivative_lengths, drop_vector, get_singularities, is_singularity,
};
use crate::curve_surface_global::degenerate_on_to;
use crate::extent::{Extent1, Extent2, PointType};
use crate::geometry::{line_point_distance, vec_lin_comb_two_vectors, vector_project_to_plane};
use crate::global_solver::{SolverOperation, SolverRefusal};
use crate::hermite::Hermite;
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ, PI, ZONE_TOL_3D, degrees, is_between, max, min};
use crate::spline::{BSplineCurve, BSplineSurface};
use crate::surface_drop_curve::{Drop, UvCurve, drop_iso_curve};
use crate::surface_eval::{
    evaluate_normal_out_of_bounds, evaluate_out_of_bounds, evaluate_point, evaluate_simple,
};
use crate::surface_global::global_point_solve;
use crate::surface_solve::{
    find_degenerate_parameter, is_closed, is_on_seam, local_point_solve, map_to_2d,
};
use crate::surface_trace::principal_at_sides;
use crate::vector::{Vector2, Vector3};

const INTERSECT: SolverOperation = SolverOperation::Intersect;

/// `SM_EFF_ZERO_PARAM`: the zero for knot parameter values.
const EFF_ZERO_PARAM: f64 = 1.2e-8;

/// A place where SMLib goes on with a value a failed call left undefined.
fn unported() -> SolverRefusal {
    SolverRefusal::NotPorted(INTERSECT)
}

/// SMLib's tri-state `SmBoolean`; `UNSURE` (an uninitialised flag) reads as
/// true.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Tri {
    False,
    True,
    Unsure,
}

impl Tri {
    fn of(value: bool) -> Self {
        if value { Self::True } else { Self::False }
    }

    fn truthy(self) -> bool {
        self != Self::False
    }

    /// `b ? FALSE : TRUE`.
    fn flipped(self) -> Self {
        if self.truthy() {
            Self::False
        } else {
            Self::True
        }
    }
}

/// The status codes the walk tells apart.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Status {
    Success,
    /// A failed call: `SM_ERR` for [`Fault::Degenerate`].
    Failed(Fault),
    OutsideOfDomain,
    BadSurfacePoint,
    NotWithinTolerance,
    BadTangentDrop,
    /// `SM_BIG_ULONG`, the unset derivative status of a degenerate curve.
    Unset,
}

/// `SM_ERR`.
const ERR: Status = Status::Failed(Fault::Degenerate);

/// SMLib's `SER` inside a routine that returns a status: a failure becomes
/// the returned status, a refusal to port passes up.
macro_rules! ser {
    ($result:expr) => {
        match $result {
            Ok(value) => value,
            Err(error) => return status_of(error.into()),
        }
    };
}

fn status_of(refusal: SolverRefusal) -> Result<Status, SolverRefusal> {
    match refusal {
        SolverRefusal::Failed(fault) => Ok(Status::Failed(fault)),
        other => Err(other),
    }
}

/// `SmPointObjectContainmentType` as `IsMovingAlongBoundary` answers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Containment {
    Unknown,
    Inside,
    Outside,
    OnBoundary,
}

/// `SmSurfParamType` for closures and seams.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Closure {
    Neither,
    U,
    V,
    Both,
}

impl Closure {
    fn has_u(self) -> bool {
        matches!(self, Self::U | Self::Both)
    }

    fn has_v(self) -> bool {
        matches!(self, Self::V | Self::Both)
    }
}

/// The derivative index of `FindDegenParamForDirection`: `SM_SP_U` moves
/// `u`, anything else `v`.
fn degenerate_axis(direction: Singular) -> usize {
    if direction == Singular::U { 0 } else { 1 }
}

/// The surface evaluation SMLib runs during the drop: its out-of-bounds
/// switch is on, so parameters outside the domain extrapolate.
fn surface_evaluate(
    surface: &BSplineSurface,
    uv: Vector2,
    order: usize,
    sides: [bool; 2],
) -> Result<Vec<Vector3>, Fault> {
    evaluate_out_of_bounds(surface, uv, order, order, sides, true, true, true)
}

/// `Evaluate1stDerivatives`: position, `Su`, `Sv`.
fn first_derivatives(
    surface: &BSplineSurface,
    uv: Vector2,
    sides: [bool; 2],
) -> Result<[Vector3; 3], Fault> {
    let d = surface_evaluate(surface, uv, 1, sides)?;
    Ok([d[0], d[2], d[1]])
}

/// `sm_FindNonsingularStepoff`: a parameter step from `uv` toward the
/// domain's middle at which the base surface has a first derivative longer
/// than `tolerance` in the direction that degenerates.
fn find_nonsingular_stepoff(
    surface: &BSplineSurface,
    uv: Vector2,
    domain: Extent2,
    bu_length: f64,
    bv_length: f64,
    tolerance: f64,
) -> Result<Vector2, Fault> {
    let u_is_degenerate = bu_length < bv_length;
    let middle = Vector2::new(
        (domain.min.x + domain.max.x) / 2.0,
        (domain.min.y + domain.max.y) / 2.0,
    );
    let mut step_size = 0.001;
    while step_size < 0.2 {
        for along_v in [u_is_degenerate, !u_is_degenerate] {
            let step = if along_v {
                Vector2::new(0.0, (middle.y - uv.y) * step_size)
            } else {
                Vector2::new((middle.x - uv.x) * step_size, 0.0)
            };
            let cells = evaluate_simple(surface, uv + step, 1, 1, [true, true], true)?;
            // `sMat[2]` is the U derivative and `sMat[1]` the V derivative.
            let longer = if u_is_degenerate {
                cells[2].length() > tolerance
            } else {
                cells[1].length() > tolerance
            };
            if longer {
                return Ok(step);
            }
        }
        step_size *= 10.0;
    }
    Err(Fault::Degenerate)
}

/// `SmOffsetSurface` built by the drop (offset 1, then moved through each
/// drop point): its first-order evaluation as `EvaluatePlainOrSimple` does
/// it for a full upper-half request of one derivative each way.
struct Offset<'a> {
    base: &'a BSplineSurface,
    distance: f64,
}

impl Offset<'_> {
    /// Position, `Su` and `Sv`; `simple` skips the repair of vanishing
    /// derivatives.
    fn first_order(
        &self,
        uv_in: Vector2,
        sides: [bool; 2],
        simple: bool,
    ) -> Result<[Vector3; 3], SolverRefusal> {
        let surface = self.base;
        let offset = self.distance;
        let domain = surface.natural_domain();
        let uv = if domain.contains_point(uv_in, 0.0) {
            uv_in
        } else {
            domain.clamp_point(uv_in)
        };
        if offset == 0.0 {
            let cells = if simple {
                evaluate_simple(surface, uv, 1, 1, sides, true)?
            } else {
                crate::surface_eval::evaluate(surface, uv, 1, 1, sides, true, true, true)?
            };
            return Ok([cells[0], cells[2], cells[1]]);
        }
        let base = evaluate_simple(surface, uv, 2, 2, [true, true], true)?;
        let cell = |u: usize, v: usize| base[u * 3 + v];
        let b = cell(0, 0);
        let (mut bu, mut bv) = (cell(1, 0), cell(0, 1));
        let tolerance = 10.0 * EFF_ZERO * (1.0 + b.max_dimension());
        let (bu_length, bv_length) = (bu.length(), bv.length());
        if !simple && (bu_length < tolerance || bv_length < tolerance) {
            // At a pole: extrapolate from two steps off it.
            let step =
                find_nonsingular_stepoff(surface, uv_in, domain, bu_length, bv_length, tolerance)?;
            let first = self.first_order(uv_in + step, sides, simple)?;
            let second = self.first_order(uv_in + 2.0 * step, sides, simple)?;
            let tangent = |index: usize| -> Result<Vector3, Fault> {
                let (mut v1, mut v2) = (first[index], second[index]);
                let mut v0 = 2.0 * v1 - v2;
                if v0.length() < tolerance {
                    v1.unitize()?;
                    v2.unitize()?;
                    v0 = 2.0 * v1 - v2;
                    v0.unitize()?;
                    v0 *= tolerance * 1.1;
                }
                Ok(v0)
            };
            return Ok([2.0 * first[0] - second[0], tangent(1)?, tangent(2)?]);
        }
        let mut n = bu.cross(bv);
        let mut length = n.length();
        if length < tolerance
            && let Ok([_, du, dv]) =
                crate::surface_eval::first_derivatives(surface, uv, [true, true])
        {
            n = du.cross(dv);
            length = n.length();
            if !simple {
                bu = du;
                bv = dv;
            }
        }
        if !simple {
            let (bu_length, bv_length) = (bu.length(), bv.length());
            let smallest = length.min(bu_length).min(bv_length);
            if smallest < tolerance {
                n = evaluate_normal_out_of_bounds(surface, uv, sides);
                length = n.length();
                if bu_length < bv_length {
                    let mut du = bv.cross(n);
                    let _ = du.unitize();
                    du *= tolerance * 1.1;
                    bu = du;
                } else {
                    let mut dv = n.cross(bu);
                    let _ = dv.unitize();
                    dv *= tolerance * 1.1;
                    bv = dv;
                }
            }
        }
        if length < EFF_ZERO {
            return Err(Fault::Degenerate.into());
        }
        let c = n;
        let m = c.length();
        let n = c / m;
        let point = b + offset * n;
        let (buu, buv, bvv) = (cell(2, 0), cell(1, 1), cell(0, 2));
        let cu = buu.cross(bv) + bu.cross(buv);
        let mu = n.dot(cu);
        let nu = (cu - mu * n) / m;
        let mut su = bu + offset * nu;
        let cv = buv.cross(bv) + bu.cross(bvv);
        let mv = n.dot(cv);
        let nv = (cv - mv * n) / m;
        let mut sv = bv + offset * nv;
        if !simple {
            let tolerance = EFF_ZERO * (1.0 + point.max_dimension());
            if su.length_squared() < tolerance * tolerance
                && let Ok([_, mut du, _]) =
                    crate::surface_eval::first_derivatives(surface, uv, [true, true])
                && du.length_squared() > tolerance * tolerance
            {
                let _ = du.unitize();
                su = du * tolerance * 1.1;
            }
            if sv.length_squared() < tolerance * tolerance
                && let Ok([_, _, mut dv]) =
                    crate::surface_eval::first_derivatives(surface, uv, [true, true])
                && dv.length_squared() > tolerance * tolerance
            {
                let _ = dv.unitize();
                sv = dv * tolerance * 1.1;
            }
        }
        Ok([point, su, sv])
    }

    /// `SmSurface::DropVectors` for one vector on the offset.
    fn drop_vector(
        &self,
        uv: Vector2,
        sides: [bool; 2],
        vector: Vector3,
    ) -> Result<Vector2, SolverRefusal> {
        let [_, du, dv] = self.first_order(uv, sides, false)?;
        drop_vector(du, dv, vector)
    }

    /// `EvaluateDirectionalDerivs` to first order on the offset: the point
    /// and the derivative along the unitised `direction`.
    fn directional(
        &self,
        uv: Vector2,
        direction: Vector2,
        sides: [bool; 2],
    ) -> Result<[Vector3; 2], SolverRefusal> {
        let [point, su, sv] = self.first_order(uv, sides, true)?;
        let unit = unitized2(direction);
        Ok([point, unit.x * su + unit.y * sv])
    }
}

/// `SmVector2d::Unitize` on a copy (a vector shorter than the squared zero
/// stays as it is).
fn unitized2(mut v: Vector2) -> Vector2 {
    let length_squared = v.length_squared();
    if length_squared >= EFF_ZERO_SQ && (length_squared - 1.0).abs() > EFF_ZERO_SQ {
        let length = length_squared.sqrt();
        v = Vector2::new(v.x / length, v.y / length);
    }
    v
}

/// `smsurf_IsPointOnNormal` with a clamping domain: whether `test` lies on
/// the surface normal line at `uv` within `tolerance` (scaled by the drop
/// distance past one), the drop distance and the deviation from the line.
fn is_point_on_normal(
    surface: &BSplineSurface,
    domain: Extent2,
    uv: Vector2,
    test: Vector3,
    tangent: Vector3,
    tolerance: f64,
) -> (bool, f64, f64) {
    let clamped = domain.clamp_point(uv);
    let Ok([mut point, mut du, mut dv]) = first_derivatives(surface, clamped, [true, true]) else {
        return (false, 0.0, 0.0);
    };
    let sides = [du.dot(tangent) >= 0.0, dv.dot(tangent) >= 0.0];
    if !sides[0] || !sides[1] {
        let Ok(again) = first_derivatives(surface, clamped, sides) else {
            return (false, 0.0, 0.0);
        };
        [point, du, dv] = again;
    }
    let mut normal = du.cross(dv);
    if normal.length_squared() < EFF_ZERO_SQ {
        normal = evaluate_normal_out_of_bounds(surface, clamped, sides);
    } else {
        let _ = normal.unitize();
    }
    let (deviation, t) = line_point_distance(point, normal, test);
    let difference = (point + t * normal) - test;
    let mut to_surface = t.abs();
    if to_surface < 1.0 {
        to_surface = 1.0;
    }
    if difference.length() < tolerance * to_surface {
        (true, t.abs(), deviation)
    } else {
        (false, 0.0, deviation)
    }
}

/// `smgu_IsDogLeg` for four UV points with the default angles (35 degrees
/// of bend, 45 of offset): a `z` shape.
fn is_dog_leg(p0: Vector2, p1: Vector2, p2: Vector2, p3: Vector2) -> Result<bool, SolverRefusal> {
    let lift = |v: Vector2| Vector3::new(v.x, v.y, 0.0);
    let begin = lift(p1 - p0);
    let middle = lift(p2 - p1);
    let end = lift(p3 - p2);
    let total = lift(p3 - p0);
    // SMLib reads an unset angle when a leg has no length.
    let angle = begin.angle_between(end).map_err(|_| unported())?;
    Ok(degrees(angle) < 45.0
        && begin.is_perpendicular_to(middle, 35.0)
        && end.is_perpendicular_to(middle, 35.0)
        && total.is_parallel_to(middle, 45.0).0)
}

/// `sm_RemoveDuplicateDropPointSolutions`: drops later solutions within
/// the UV tolerances of an earlier one.
fn remove_duplicates(solutions: &mut Vec<Vector2>, tol_u: f64, tol_v: f64) {
    let mut ii = 1;
    while ii < solutions.len() {
        let uv = solutions[ii];
        let duplicate = solutions[..ii]
            .iter()
            .rev()
            .any(|other| (uv.x - other.x).abs() < tol_u && (uv.y - other.y).abs() < tol_v);
        if duplicate {
            solutions.remove(ii);
        } else {
            ii += 1;
        }
    }
}

/// `SmDropPt`: a curve point dropped to the surface and its properties.
#[derive(Clone, Copy, Debug)]
struct DropPoint {
    t: f64,
    good_drop: Tri,
    curve: [Vector3; 3],
    uv: [Vector2; 2],
    surface: [Vector3; 2],
    normal: Vector3,
    from_left: [Tri; 2],
    drop_to_surface: f64,
    approx_deviation: f64,
    normal_to_curve_deg: f64,
    surface_to_curve_deg: f64,
    speed_ratio: f64,
    on_boundary: Tri,
    on_seam: Tri,
    on_pole: Tri,
    degenerate_curve: Tri,
    leaving: Tri,
    entering: Tri,
    derivative_status: Status,
}

/// What the drop works with besides the moving points.
struct Walk<'a> {
    surface: &'a BSplineSurface,
    curve: &'a BSplineCurve,
    domain: Extent2,
    tolerance: f64,
    singularities: u32,
    closure: Closure,
    tol_u: f64,
    tol_v: f64,
}

impl DropPoint {
    /// The constructor: everything unset but the parameter.
    fn new(t: f64, degenerate_curve: bool) -> Self {
        let mut point = Self::uninitialised();
        point.t = t;
        point.degenerate_curve = Tri::of(degenerate_curve);
        point
    }

    /// `SetUninitialized`.
    fn uninitialised() -> Self {
        let undefined = crate::scalar::UNDEFINED;
        Self {
            t: undefined,
            good_drop: Tri::Unsure,
            curve: [Vector3::default(); 3],
            uv: [Vector2::default(); 2],
            surface: [Vector3::default(); 2],
            normal: Vector3::default(),
            from_left: [Tri::Unsure; 2],
            drop_to_surface: undefined,
            approx_deviation: undefined,
            normal_to_curve_deg: undefined,
            surface_to_curve_deg: undefined,
            speed_ratio: undefined,
            on_boundary: Tri::Unsure,
            on_seam: Tri::Unsure,
            on_pole: Tri::Unsure,
            degenerate_curve: Tri::Unsure,
            leaving: Tri::Unsure,
            entering: Tri::Unsure,
            derivative_status: Status::Unset,
        }
    }

    fn sides(&self) -> [bool; 2] {
        [self.from_left[0].truthy(), self.from_left[1].truthy()]
    }

    /// `IsGoodTangentDrop`.
    fn is_good_tangent_drop(&self, on_pole: bool) -> bool {
        if on_pole {
            self.surface_to_curve_deg < 15.0 && self.speed_ratio < 4.0 && self.speed_ratio > 0.25
        } else {
            self.surface_to_curve_deg < 3.0 && self.speed_ratio < 2.0 && self.speed_ratio > 0.5
        }
    }

    /// `IsMovingAlongBoundary` for one parameter direction (`axis` 0 tests
    /// the `u` sides): inside, outside or along the side the drop point
    /// lies on, from the curve tangent's angle to the side's tangent.
    fn is_moving_along_boundary(
        &self,
        surface: &BSplineSurface,
        domain: Extent2,
        tol_u: f64,
        tol_v: f64,
        limit_deg: f64,
        axis: usize,
    ) -> Result<Containment, SolverRefusal> {
        let limit = limit_deg * (PI / 180.0);
        let classify = |angle: f64, low: bool, u_side: bool| {
            let (first, second) = if u_side == low {
                (Containment::Inside, Containment::Outside)
            } else {
                (Containment::Outside, Containment::Inside)
            };
            if angle > -PI + limit && angle < 0.0 - limit {
                first
            } else if angle > 0.0 + limit && angle < PI - limit {
                second
            } else {
                Containment::OnBoundary
            }
        };
        let tangent = self.curve[1];
        let (range, value, tol) = if axis == 0 {
            (domain.u_interval(), self.uv[0].x, tol_u)
        } else {
            (domain.v_interval(), self.uv[0].y, tol_v)
        };
        if !range.is_value_on_boundary(value, tol).0 {
            return Ok(Containment::Unknown);
        }
        let low = value < range.mid();
        let sides = if axis == 0 {
            [low, false]
        } else {
            [false, low]
        };
        let [_, su, sv] = first_derivatives(surface, self.uv[0], sides).map_err(|_| unported())?;
        let normal = evaluate_normal_out_of_bounds(surface, self.uv[0], sides);
        let along = if axis == 0 { sv } else { su };
        if axis == 0
            && (along.length_squared() < EFF_ZERO_SQ || tangent.length_squared() < EFF_ZERO_SQ)
        {
            return Ok(Containment::Unknown);
        }
        let angle = crate::surface_curves::ccw_angle_between(normal, along, tangent)
            .map_err(|_| unported())?;
        let mut containment = classify(angle, low, axis == 0);
        if containment == Containment::Inside {
            let other_uv = if axis == 0 {
                Vector2::new(if low { range.max } else { range.min }, self.uv[0].y)
            } else {
                Vector2::new(self.uv[0].x, if low { range.max } else { range.min })
            };
            // SMLib evaluates the far side from the left in u for both axes.
            let other_sides = [!low, false];
            let [_, other_su, other_sv] =
                first_derivatives(surface, other_uv, other_sides).map_err(|_| unported())?;
            let other_normal = evaluate_normal_out_of_bounds(surface, other_uv, other_sides);
            let other_along = if axis == 0 { other_sv } else { other_su };
            let other_angle =
                crate::surface_curves::ccw_angle_between(other_normal, other_along, tangent)
                    .map_err(|_| unported())?;
            if classify(other_angle, !low, axis == 0) == Containment::Inside {
                let other_normal_angle = other_normal
                    .angle_between(tangent)
                    .map_err(|_| unported())?;
                containment = if (self.normal_to_curve_deg - 90.0).abs()
                    < (degrees(other_normal_angle) - 90.0).abs()
                {
                    Containment::Inside
                } else {
                    Containment::Outside
                };
            }
        }
        Ok(containment)
    }

    /// `SetProperties`: everything about the drop of `curve[0]` (with
    /// tangent `curve[1]`) to `uv[0]`; moves the offset through the curve
    /// point.
    fn set_properties(
        &mut self,
        walk: &Walk,
        offset: &mut Offset,
        adjust_pole_point: bool,
    ) -> Result<Status, SolverRefusal> {
        let surface = walk.surface;
        let tolerance = walk.tolerance;
        let mut direction = Singular::U;
        let on_pole = walk.singularities != 0
            && match is_singularity(surface, self.uv[0], tolerance, ZONE_TOL_3D, false) {
                Ok(Some(found)) => {
                    direction = found;
                    true
                }
                Ok(None) => false,
                Err(_) => return Err(unported()),
            };
        self.on_pole = Tri::of(on_pole);
        if on_pole && self.degenerate_curve != Tri::True && adjust_pole_point {
            let _ = find_degenerate_parameter(
                surface,
                &mut self.uv[0],
                self.curve[1],
                true,
                degenerate_axis(direction),
                0.0001,
            );
        }
        let [mut point, mut du, mut dv] =
            ser!(first_derivatives(surface, self.uv[0], [true, true]));
        self.surface[0] = point;
        self.from_left = [
            Tri::of(du.dot(self.curve[1]) >= 0.0),
            Tri::of(dv.dot(self.curve[1]) >= 0.0),
        ];
        if self.from_left != [Tri::True, Tri::True] {
            [point, du, dv] = ser!(first_derivatives(surface, self.uv[0], self.sides()));
            self.surface[0] = point;
        }
        self.normal = evaluate_normal_out_of_bounds(surface, self.uv[0], self.sides());
        let gap = self.curve[0] - self.surface[0];
        let (u_range, v_range) = (walk.domain.u_interval(), walk.domain.v_interval());
        let zone_u = map_to_2d(surface, ZONE_TOL_3D, self.uv[0], Vector2::new(1.0, 0.0));
        let zone_v = map_to_2d(surface, ZONE_TOL_3D, self.uv[0], Vector2::new(0.0, 1.0));
        let limit = 30.0 * (PI / 180.0);
        self.on_boundary = Tri::False;
        self.leaving = Tri::False;
        self.entering = Tri::False;
        self.on_seam = Tri::False;
        let normal = self.normal;
        let ccw = move |axis_vector: Vector3, tangent: Vector3| -> Result<f64, SolverRefusal> {
            let mut angle = crate::surface_curves::ccw_angle_between(normal, axis_vector, tangent)
                .map_err(|_| unported())?;
            if angle < 0.0 {
                angle += 2.0 * PI;
            }
            Ok(angle)
        };
        let (mut leaving, mut entering) = (self.leaving, self.entering);
        if self.uv[0].x < u_range.min + zone_u || self.uv[0].x > u_range.max - zone_u {
            self.on_boundary = Tri::True;
            if is_on_seam(surface, self.uv[0], tolerance).1 {
                self.on_seam = Tri::True;
            }
            let lengths = du.length_squared() * self.curve[1].length_squared();
            if !(on_pole && walk.singularities & (SINGULAR_V_MIN | SINGULAR_V_MAX) != 0)
                && lengths > EFF_ZERO_SQ
            {
                let angle = ccw(dv, self.curve[1])?;
                if self.uv[0].x < u_range.min + zone_u {
                    leaving = Tri::of(is_between(angle, limit, PI - limit));
                    entering = Tri::of(is_between(angle, PI + limit, 2.0 * PI - limit));
                } else {
                    leaving = Tri::of(is_between(angle, PI + limit, 2.0 * PI - limit));
                    entering = Tri::of(is_between(angle, limit, PI - limit));
                }
            }
        }
        if self.uv[0].y < v_range.min + zone_v || self.uv[0].y > v_range.max - zone_v {
            self.on_boundary = Tri::True;
            if is_on_seam(surface, self.uv[0], tolerance).1 {
                self.on_seam = Tri::True;
            }
            let lengths = dv.length_squared() * self.curve[1].length_squared();
            if !(on_pole && walk.singularities & (SINGULAR_U_MIN | SINGULAR_U_MAX) != 0)
                && lengths > EFF_ZERO_SQ
            {
                let angle = ccw(du, self.curve[1])?;
                if self.uv[0].y < v_range.min + zone_v {
                    leaving = Tri::of(is_between(angle, PI + limit, 2.0 * PI - limit));
                    entering = Tri::of(is_between(angle, limit, PI - limit));
                } else {
                    leaving = Tri::of(is_between(angle, limit, PI - limit));
                    entering = Tri::of(is_between(angle, PI + limit, 2.0 * PI - limit));
                }
            }
        }
        self.leaving = leaving;
        self.entering = entering;
        let (deviation, drop) = line_point_distance(self.surface[0], self.normal, self.curve[0]);
        self.approx_deviation = deviation;
        self.drop_to_surface = drop;
        let scaled_tolerance = (if drop < 1.0 { 1.0 } else { drop }) * tolerance;
        let mut good = deviation < scaled_tolerance;
        if self.on_boundary == Tri::False && !good && deviation < 1000.0 * scaled_tolerance {
            let (k1, k2) = match principal_at_sides(surface, self.uv[0], self.sides()) {
                Ok((k1, k2, _, _)) => (k1, k2),
                Err(_) => (0.0, 0.0),
            };
            let center1 = self.surface[0] + self.normal * k1;
            let center2 = self.surface[0] + self.normal * k2;
            let to_center1 = (center1 - self.curve[0]).length();
            let to_center2 = (center2 - self.curve[0]).length();
            let scaling = max(max(k1.abs() * to_center1, k2.abs() * to_center2), 5.0);
            good = deviation < scaling * tolerance;
        }
        self.good_drop = Tri::of(good);
        let speed = self.curve[1].length();
        if speed < EFF_ZERO && self.degenerate_curve == Tri::False {
            self.good_drop = Tri::False;
            return Ok(ERR);
        }
        if self.degenerate_curve != Tri::True {
            let angle = self
                .normal
                .angle_between(self.curve[1])
                .map_err(|_| unported())?;
            if degrees(angle) < 3.0 {
                self.good_drop = Tri::False;
                return Ok(ERR);
            }
        }
        offset.distance = gap.dot(self.normal);
        if self.degenerate_curve == Tri::False {
            self.uv[0] = walk.domain.clamp_point(self.uv[0]);
            self.uv[1] = ser!(offset.drop_vector(self.uv[0], self.sides(), self.curve[1]));
        }
        if on_pole {
            if direction == Singular::U {
                self.uv[1].x = 0.0;
            } else {
                self.uv[1].y = 0.0;
            }
        }
        if self.degenerate_curve != Tri::True {
            if self.uv[1].length_squared() >= EFF_ZERO_SQ {
                let [_, along] = offset
                    .directional(self.uv[0], self.uv[1], self.sides())
                    .map_err(|_| unported())?;
                let uv_speed = self.uv[1].length();
                let offset_tangent = uv_speed * along;
                let projected = vector_project_to_plane(self.curve[1], self.normal);
                self.normal_to_curve_deg =
                    degrees(self.curve[1].angle_between(self.normal).unwrap_or(0.0));
                self.surface_to_curve_deg =
                    degrees(projected.angle_between(offset_tangent).unwrap_or(0.0));
                let surface_speed = offset_tangent.length();
                let projected_speed = projected.length();
                self.speed_ratio = if projected_speed > EFF_ZERO {
                    surface_speed / projected_speed
                } else if surface_speed < EFF_ZERO {
                    1.0
                } else {
                    10000.0
                };
                self.surface[1] = du * self.uv[1].x + dv * self.uv[1].y;
            }
            self.derivative_status = if self.is_good_tangent_drop(on_pole) {
                Status::Success
            } else {
                Status::BadTangentDrop
            };
        }
        Ok(if self.good_drop == Tri::False {
            Status::NotWithinTolerance
        } else {
            self.derivative_status
        })
    }

    /// `CheckSeamJump`: an end that landed across a seam from this start is
    /// snapped back to the start's side when that moves it no more than the
    /// tolerance; otherwise the span is refused as leaving the domain.
    fn check_seam_jump(&self, walk: &Walk, end: &mut DropPoint) -> Result<Status, SolverRefusal> {
        if walk.closure == Closure::Neither {
            return Ok(Status::Success);
        }
        let surface = walk.surface;
        let natural = surface.natural_domain();
        for axis in 0..2 {
            let (closed, range, delta) = if axis == 0 {
                (
                    walk.closure.has_u(),
                    natural.u_interval(),
                    end.uv[0].x - self.uv[0].x,
                )
            } else {
                (
                    walk.closure.has_v(),
                    natural.v_interval(),
                    end.uv[0].y - self.uv[0].y,
                )
            };
            if !closed || delta.abs() <= range.length() / 2.0 {
                continue;
            }
            let [start3, su, sv] = ser!(first_derivatives(surface, self.uv[0], self.sides()));
            let end3 = ser!(evaluate_point(surface, end.uv[0]));
            let step = end3 - start3;
            let (across, expected) = if axis == 0 {
                (step.cross(sv), delta * su.cross(sv))
            } else {
                (step.cross(su), delta * sv.cross(su))
            };
            let numerator = across.dot(expected);
            let denominator = across.dot(across);
            if numerator < denominator * -2.0 {
                let mut snapped = end.uv[0];
                if axis == 0 {
                    snapped.x = if end.uv[0].x > range.mid() {
                        range.min
                    } else {
                        range.max
                    };
                } else {
                    snapped.y = if end.uv[0].y > range.mid() {
                        range.min
                    } else {
                        range.max
                    };
                }
                let snapped_point = ser!(evaluate_point(surface, snapped));
                if snapped_point.distance(end3) <= walk.tolerance {
                    end.uv[0] = snapped;
                } else {
                    return Ok(Status::OutsideOfDomain);
                }
            }
        }
        Ok(Status::Success)
    }

    /// `GetNextPoint`: drops the curve point at `end_t` from a Taylor guess
    /// off this start, by the local solver and then the global one.
    fn get_next_point(
        &self,
        walk: &Walk,
        offset: &mut Offset,
        end_t: f64,
        end: &mut DropPoint,
    ) -> Result<Status, SolverRefusal> {
        let surface = walk.surface;
        let domain = walk.domain;
        let mut status = ERR;
        *end = DropPoint::uninitialised();
        end.t = end_t;
        end.degenerate_curve = self.degenerate_curve;
        let zero_u = EFF_ZERO * (1.0 + domain.u_interval().length());
        let zero_v = EFF_ZERO * (1.0 + domain.v_interval().length());
        if !domain.contains_point(self.uv[0], max(zero_u, zero_v)) {
            return Ok(Status::BadSurfacePoint);
        }
        if self.on_boundary.truthy() && self.leaving.truthy() {
            return Ok(Status::OutsideOfDomain);
        }
        // SMLib does not check this evaluation.
        let d = surface_evaluate(surface, self.uv[0], 2, self.sides())
            .unwrap_or_else(|_| vec![Vector3::default(); 9]);
        let (su, sv, suv, suu, svv) = (d[3], d[1], d[4], d[6], d[2]);
        let (ut, vt) = (self.uv[1].x, self.uv[1].y);
        let second = suu * ut * ut + 2.0 * suv * ut * vt + svv * vt * vt;
        let rhs = self.curve[2] - second;
        let uvtt = vec_lin_comb_two_vectors(su, sv, rhs).unwrap_or(Vector2::new(0.0, 0.0));
        let dt = end.t - self.t;
        end.uv[0] = self.uv[0] + dt * self.uv[1];
        end.uv[0] = end.uv[0] + uvtt * (dt * dt / 2.0);
        let step = end.uv[0] - self.uv[0];
        let size = domain.size();
        if step.x.abs() > 1.5 * size.x || step.y.abs() > 1.5 * size.y {
            end.uv[0] = self.uv[0];
        }
        end.uv[0] = domain.clamp_point(end.uv[0]);
        let derivatives = ser!(walk.curve.evaluate(end.t, 2, false, true));
        end.curve = [derivatives[0], derivatives[1], derivatives[2]];
        let local = ser!(local_point_solve(
            surface,
            domain,
            SolverOperation::Minimize,
            end.curve[0],
            end.uv[0],
            None
        ));
        if !local.found {
            return Ok(ERR);
        }
        end.uv[0] = local.uv;
        let mut solution = local.uv;
        for ii in 0..2 {
            if status == Status::Success {
                break;
            }
            if ii == 1 {
                let solutions = ser!(global_point_solve(
                    surface,
                    domain,
                    SolverOperation::Minimize,
                    end.curve[0],
                    EFF_ZERO,
                    None
                ));
                let mut smallest = 0.0;
                for (jj, candidate) in solutions.as_slice().iter().enumerate() {
                    let uv =
                        Vector2::new(candidate.start.parameters[0], candidate.start.parameters[1]);
                    let difference = (uv.x - end.uv[0].x).abs() + (uv.y - end.uv[0].y).abs();
                    if jj == 0 || (difference < smallest && candidate.start.value < walk.tolerance)
                    {
                        smallest = difference;
                        solution = uv;
                    }
                }
                if end.uv[0].closer_than(EFF_ZERO, solution) {
                    break;
                }
            }
            end.uv[0] = solution;
            if let Some(direction) =
                is_singularity(surface, end.uv[0], walk.tolerance, ZONE_TOL_3D, false)
                    .map_err(|_| unported())?
            {
                let _ = find_degenerate_parameter(
                    surface,
                    &mut end.uv[0],
                    end.curve[1],
                    false,
                    degenerate_axis(direction),
                    0.0001,
                );
            }
            status = self.check_seam_jump(walk, end)?;
            if status != Status::Success {
                return Ok(status);
            }
            end.curve[1] = -end.curve[1];
            status = end.set_properties(walk, offset, false)?;
            if status != ERR {
                end.curve[1] = -end.curve[1];
                end.uv[1] = -end.uv[1];
                if end.on_boundary.truthy() {
                    if end.leaving.truthy() {
                        end.leaving = Tri::False;
                        end.entering = Tri::True;
                    } else if end.entering.truthy() {
                        end.leaving = Tri::True;
                        end.entering = Tri::False;
                    }
                }
            }
        }
        Ok(status)
    }

    /// `FixBoundaryTangents`: on a span running along a side of the domain,
    /// the tangent components across that side become zero. Returns whether
    /// anything changed.
    fn fix_boundary_tangents(
        &mut self,
        domain: Extent2,
        end: &mut DropPoint,
        tol_u: f64,
        tol_v: f64,
    ) -> bool {
        let same = |a: f64, b: f64, tol: f64| (a - b).abs() <= tol;
        let u_step = (end.uv[0].x - self.uv[0].x).abs();
        let v_step = (end.uv[0].y - self.uv[0].y).abs();
        let mut changed = false;
        for side in [domain.min.x, domain.max.x] {
            if same(self.uv[0].x, side, tol_u) && same(end.uv[0].x, side, tol_u) && v_step > u_step
            {
                self.uv[1].x = 0.0;
                end.uv[1].x = 0.0;
                changed = true;
            }
        }
        for side in [domain.min.y, domain.max.y] {
            if same(self.uv[0].y, side, tol_v) && same(end.uv[0].y, side, tol_v) && u_step > v_step
            {
                self.uv[1].y = 0.0;
                end.uv[1].y = 0.0;
                changed = true;
            }
        }
        changed
    }

    /// `TestSpanAccuracy`: whether the Hermite span from this point to
    /// `end` lies on the surface normals of the curve at three interior
    /// samples; also the largest drop distance, its parameter and the
    /// largest deviation.
    fn test_span_accuracy(
        &self,
        walk: &Walk,
        end: &DropPoint,
    ) -> Result<(bool, f64, f64, f64), SolverRefusal> {
        let surface = walk.surface;
        let mut passes = true;
        let (mut max_drop, mut max_parameter, mut max_deviation) = (0.0, self.t, 0.0);
        let interval = Extent1::new(self.t, end.t);
        let dt = end.t - self.t;
        if dt < interval.length() / 1000.0 {
            return Ok((true, max_drop, max_parameter, max_deviation));
        }
        let hermite = Hermite::planar(self.uv[0], self.uv[1] * dt, end.uv[0], end.uv[1] * dt);
        const SAMPLES: usize = 3;
        for ii in 0..SAMPLES {
            let normalized = (ii as f64 + 1.0) / (SAMPLES as f64 + 1.0);
            let h = hermite.evaluate_point(normalized);
            let mut uv = Vector2::new(h.x, h.y);
            if !walk.domain.contains_point(uv, EFF_ZERO) {
                let clamped = walk.domain.clamp_point(uv);
                let from_uv = evaluate_point(surface, uv).unwrap_or_default();
                let from_clamped = evaluate_point(surface, clamped).unwrap_or_default();
                if from_uv.distance(from_clamped) > walk.tolerance {
                    return Ok((false, max_drop, max_parameter, max_deviation));
                }
                uv = clamped;
            }
            let parameter = self.t + normalized * dt;
            let point = walk.curve.evaluate(parameter, 1, true, true)?;
            let (good, drop, deviation) =
                is_point_on_normal(surface, walk.domain, uv, point[0], point[1], walk.tolerance);
            if !good {
                passes = false;
                break;
            }
            if drop > max_drop {
                max_drop = drop;
                max_parameter = parameter;
            }
            if deviation > max_deviation {
                max_deviation = deviation;
            }
        }
        Ok((passes, max_drop, max_parameter, max_deviation))
    }
}

/// The growing parameter curve: Bezier control points and span knots.
#[derive(Default)]
struct Growing {
    points: Vec<Vector3>,
    knots: Vec<f64>,
}

impl Growing {
    fn clear(&mut self) {
        self.points.clear();
        self.knots.clear();
    }

    /// `sm_FixBoundaryCurves`: knot points clamped into the domain, and
    /// their neighbours put on a side the knot point lies on.
    fn fix_boundary_curves(&mut self, domain: Extent2, tol_u: f64, tol_v: f64) {
        let count = self.points.len();
        if count < 4 {
            return;
        }
        let spans = count / 3;
        let (u0, v0, u1, v1) = (domain.min.x, domain.min.y, domain.max.x, domain.max.y);
        for ii in 0..=spans {
            let this = 3 * ii;
            let mut knot = self.points[this];
            if knot.x < u0 {
                knot.x = u0;
            }
            if knot.x > u1 {
                knot.x = u1;
            }
            if knot.y < v0 {
                knot.y = v0;
            }
            if knot.y > v1 {
                knot.y = v1;
            }
            self.points[this] = knot;
            let on_u = knot.x < u0 + tol_u || knot.x > u1 - tol_u;
            let on_v = knot.y < v0 + tol_v || knot.y > v1 - tol_v;
            if ii > 0 && ii < spans {
                if on_u {
                    self.points[this - 1].x = knot.x;
                    self.points[this + 1].x = knot.x;
                }
                if on_v {
                    self.points[this - 1].y = knot.y;
                    self.points[this + 1].y = knot.y;
                }
            } else if ii == 0 {
                if on_u {
                    self.points[this + 1].x =
                        domain.u_interval().clamp_value(self.points[this + 1].x);
                }
                if on_v {
                    self.points[this + 1].y =
                        domain.v_interval().clamp_value(self.points[this + 1].y);
                }
            } else {
                if on_u {
                    self.points[this - 1].x =
                        domain.u_interval().clamp_value(self.points[this - 1].x);
                }
                if on_v {
                    self.points[this - 1].y =
                        domain.v_interval().clamp_value(self.points[this - 1].y);
                }
            }
        }
    }

    /// `sm_FlushCrv`: the piecewise cubic Bezier curve built so far, added
    /// to `curves`; the growing curve restarts empty.
    fn flush(
        &mut self,
        domain: Extent2,
        tol_u: f64,
        tol_v: f64,
        curves: &mut Vec<UvCurve>,
    ) -> Result<(), Fault> {
        let count = self.knots.len();
        if count == 0 {
            return Ok(());
        }
        self.fix_boundary_curves(domain, tol_u, tol_v);
        let mut multiplicities = vec![3; count];
        multiplicities[0] = 4;
        multiplicities[count - 1] = 4;
        let curve = BSplineCurve::create_canonical(
            2,
            3,
            &self.points,
            &multiplicities,
            &self.knots,
            None,
            None,
        )?;
        curves.push(UvCurve::Spline(curve));
        self.clear();
        Ok(())
    }

    /// `sm_AddToCurve`: appends the cubic span from `start` to `end`,
    /// flushing first when it does not continue the growing curve. Its
    /// failures (a failed flush, SMLib's ten-thousand-point limit) are only
    /// signalled by SMLib's caller, so they end the call quietly here.
    fn add(
        &mut self,
        walk: &Walk,
        start: &DropPoint,
        end: &DropPoint,
        curves: &mut Vec<UvCurve>,
    ) -> Result<(), SolverRefusal> {
        let scale = end.t - start.t;
        let mut p1 = start.uv[0];
        let d1 = start.uv[1] * scale;
        let mut p4 = end.uv[0];
        let d2 = end.uv[1] * scale;
        let mut p2 = p1 + d1 / 3.0;
        let mut p3 = p4 - d2 / 3.0;
        let dog_leg = if start.on_pole.truthy() || end.on_pole.truthy() {
            is_dog_leg(p1, p2, p3, p4)?
        } else {
            false
        };
        if dog_leg {
            let at = if start.on_pole.truthy() { p1 } else { p4 };
            let direction = is_singularity(walk.surface, at, walk.tolerance, ZONE_TOL_3D, false)
                .map_err(|_| unported())?;
            if direction == Some(Singular::U) {
                if start.on_pole.truthy() {
                    p1.x = p4.x;
                    p2.x = p3.x;
                } else {
                    p4.x = p1.x;
                    p3.x = p2.x;
                }
            } else if start.on_pole.truthy() {
                p1.y = p4.y;
                p2.y = p3.y;
            } else {
                p4.y = p1.y;
                p3.y = p2.y;
            }
        }
        let lift = |v: Vector2| Vector3::new(v.x, v.y, 0.0);
        if let (Some(&last_knot), Some(&last_point)) = (self.knots.last(), self.points.last())
            && (scale < EFF_ZERO
                || (last_knot - start.t).abs() > EFF_ZERO
                || lift(p1).distance(last_point) > EFF_ZERO)
            && self
                .flush(walk.domain, walk.tol_u, walk.tol_v, curves)
                .is_err()
        {
            return Ok(());
        }
        if scale < EFF_ZERO {
            return Ok(());
        }
        p1 = walk.domain.clamp_point(p1);
        p4 = walk.domain.clamp_point(p4);
        if self.knots.is_empty() {
            self.knots.push(start.t);
            self.points.push(lift(p1));
        }
        self.knots.push(end.t);
        self.points.push(lift(p2));
        self.points.push(lift(p3));
        self.points.push(lift(p4));
        Ok(())
    }
}

/// `SmTol::GetScaledZero` of a curve: the larger of the scaled zeros of its
/// middle point and of its approximate length.
fn curve_scaled_zero(curve: &BSplineCurve) -> Result<f64, Fault> {
    let interval = curve.natural_interval();
    let length = approximate_length(curve, interval, 3);
    let middle = curve.evaluate_point(interval.evaluate(0.5))?;
    Ok(max(
        EFF_ZERO * (1.0 + middle.max_dimension()),
        EFF_ZERO * (1.0 + length.abs()),
    ))
}

/// `SmCrvOnSurf::SurfaceCurveMaxDistanceBetween` for a dropped parameter
/// curve against the 3D curve, over twenty samples (SMLib switches its
/// precise refinement off): the largest distance along the surface normals
/// and the largest deviation from them. A failed sample leaves SMLib's
/// zeros.
fn surface_curve_max_distance(
    surface: &BSplineSurface,
    uv_curve: &BSplineCurve,
    curve: &BSplineCurve,
    tolerance: f64,
) -> Result<(f64, f64), SolverRefusal> {
    match sampled_max_distance(surface, uv_curve, curve, tolerance) {
        Ok(found) => Ok(found),
        Err(SolverRefusal::Failed(_)) => Ok((0.0, 0.0)),
        Err(other) => Err(other),
    }
}

fn sampled_max_distance(
    surface: &BSplineSurface,
    uv_curve: &BSplineCurve,
    curve: &BSplineCurve,
    tolerance: f64,
) -> Result<(f64, f64), SolverRefusal> {
    let interval = uv_curve.natural_interval();
    let mut other = Extent1::new(interval.min, interval.min);
    other.add_value(interval.max);
    const SAMPLES: usize = 20;
    let (mut max_distance, mut max_deviation) = (0.0, 0.0);
    for ii in 0..SAMPLES {
        let fraction = ii as f64 / (SAMPLES - 1) as f64;
        let mut t1 = interval.evaluate(fraction);
        let guess = other.evaluate(fraction);
        let mut sample =
            distance_and_deviation(surface, uv_curve, t1, curve, guess, other, tolerance)?;
        if sample.is_none() && ii == 0 {
            t1 = interval.evaluate(0.001);
            sample = distance_and_deviation(surface, uv_curve, t1, curve, guess, other, tolerance)?;
        }
        if sample.is_none() && ii == SAMPLES - 1 {
            t1 = interval.evaluate(0.999);
            sample = distance_and_deviation(surface, uv_curve, t1, curve, guess, other, tolerance)?;
        }
        if let Some((distance, deviation)) = sample {
            if distance > max_distance {
                max_distance = distance;
            }
            if deviation > max_deviation {
                max_deviation = deviation;
            }
        }
    }
    Ok((max_distance, max_deviation))
}

/// `sm_EvalDistAndDev`: the 3D curve point nearest the surface point under
/// the parameter curve at `t1` (a B-spline answers no line intersection,
/// so SMLib drops the point from the guess), its distance along the normal
/// and its deviation from it; `None` when the drop finds nothing.
fn distance_and_deviation(
    surface: &BSplineSurface,
    uv_curve: &BSplineCurve,
    t1: f64,
    curve: &BSplineCurve,
    guess: f64,
    interval: Extent1,
    tolerance: f64,
) -> Result<Option<(f64, f64)>, SolverRefusal> {
    let uv_derivatives = uv_curve.evaluate(t1, 1, true, false)?;
    let uv = Vector2::new(uv_derivatives[0].x, uv_derivatives[0].y);
    let normal = evaluate_normal_out_of_bounds(surface, uv, [true, true]);
    let surface_point = evaluate_point(surface, uv)?;
    let dropped = curve_drop_point(
        curve,
        interval,
        surface_point,
        None,
        tolerance,
        Some(guess),
        SolverOperation::Minimize,
    )?;
    if !dropped.success {
        return Ok(None);
    }
    let curve_point = curve.evaluate_point(dropped.parameter)?;
    let gap = curve_point - surface_point;
    let distance = gap.dot(normal);
    let on_normal = surface_point + distance * normal;
    let deviation = curve_point.distance(on_normal);
    Ok(Some((distance.abs(), deviation)))
}

/// `SmSurface::DropCurve`: the parameter curves of `interval` of `curve`
/// dropped onto `surface` within `domain`, walking from each drop of the
/// curve's start, with the largest sampled drop distance and normal-line
/// deviation. Unless `keep_all`, a drop with more than one curve keeps only
/// the curves within `tolerance` when the worst is beyond it.
pub(crate) fn general_drop_curve(
    surface: &BSplineSurface,
    domain: Extent2,
    curve: &BSplineCurve,
    interval: Extent1,
    tolerance: f64,
    keep_all: bool,
) -> Result<Drop, SolverRefusal> {
    let mut result = Drop {
        curves: Vec::new(),
        max_drop: 0.0,
        deviation: 0.0,
    };
    let natural = surface.natural_domain();
    let degenerate_curve = degenerate_on_to(curve, interval, curve_scaled_zero(curve)?);
    let mut start = DropPoint::new(interval.min, degenerate_curve);
    let mut end = DropPoint::new(interval.max, degenerate_curve);
    let (mut max_drop, mut max_deviation) = (0.0_f64, 0.0_f64);
    let domain = if domain.x_length() >= EFF_ZERO && domain.y_length() >= EFF_ZERO {
        domain
    } else {
        natural
    };
    // SMLib tests the u length twice.
    if domain.x_length() < EFF_ZERO {
        return Err(Fault::Degenerate.into());
    }
    if !domain.is_contained_by(natural, EFF_ZERO) {
        return Err(Fault::Degenerate.into());
    }
    let size = approx_derivative_lengths(surface, natural);
    let tol_u = tolerance / (size.x + 1.0);
    let tol_v = tolerance / (size.y + 1.0);
    let (singularities, _) = get_singularities(surface, ZONE_TOL_3D);
    let mut offset = Offset {
        base: surface,
        distance: 1.0,
    };
    let closed_u = is_closed(surface, natural, 0, Some(tolerance));
    let closed_v = is_closed(surface, natural, 1, Some(tolerance));
    let closure = match (closed_u, closed_v) {
        (true, true) => Closure::Both,
        (true, false) => Closure::U,
        (false, true) => Closure::V,
        (false, false) => Closure::Neither,
    };
    let walk = Walk {
        surface,
        curve,
        domain,
        tolerance,
        singularities,
        closure,
        tol_u,
        tol_v,
    };
    let derivatives = curve.evaluate(start.t, 2, true, true)?;
    start.curve = [derivatives[0], derivatives[1], derivatives[2]];
    let found = global_point_solve(
        surface,
        domain,
        SolverOperation::Minimize,
        start.curve[0],
        max(tolerance / 100.0, EFF_ZERO),
        None,
    )?;
    let mut solutions: Vec<Vector2> = found
        .as_slice()
        .iter()
        .map(|s| Vector2::new(s.start.parameters[0], s.start.parameters[1]))
        .collect();
    let Some(&first) = solutions.first() else {
        return Err(Fault::Degenerate.into());
    };
    start.uv[0] = first;
    // SMLib ignores the status here and reads the drop flag.
    let _ = start.set_properties(&walk, &mut offset, true)?;
    if start.good_drop == Tri::False {
        return Ok(result);
    }
    remove_duplicates(&mut solutions, tol_u, tol_v);
    start.uv[0] = solutions[0];
    let pole = if singularities != 0 {
        is_singularity(surface, start.uv[0], tolerance, ZONE_TOL_3D, false)
            .map_err(|_| unported())?
    } else {
        None
    };
    if let Some(direction) = pole {
        if degenerate_curve {
            // The point curve drops to the pole's side line.
            let middle = natural.evaluate(0.5, 0.5);
            let (lo, hi) = (natural.min, natural.max);
            let (uv_start, uv_end) = if matches!(direction, Singular::U | Singular::Both) {
                if start.uv[0].y < middle.y {
                    (Vector2::new(lo.x, lo.y), Vector2::new(hi.x, lo.y))
                } else {
                    (Vector2::new(hi.x, hi.y), Vector2::new(lo.x, hi.y))
                }
            } else if start.uv[0].x > middle.x {
                (Vector2::new(hi.x, lo.y), Vector2::new(hi.x, hi.y))
            } else {
                (Vector2::new(lo.x, hi.y), Vector2::new(lo.x, lo.y))
            };
            let flat = |p: Vector2| Vector3::new(p.x, p.y, crate::spline::NO_Z);
            let line = BSplineCurve::create_canonical(
                2,
                1,
                &[flat(uv_start), flat(uv_end)],
                &[2, 2],
                &[0.0, 1.0],
                None,
                None,
            )?
            .edit_parameterization(interval)?;
            result.curves.push(UvCurve::Spline(line));
            result.max_drop = start.drop_to_surface;
            result.deviation = start.approx_deviation;
            return Ok(result);
        }
        for solution in &mut solutions {
            start.uv[0] = *solution;
            let moved = find_degenerate_parameter(
                surface,
                &mut start.uv[0],
                start.curve[1],
                true,
                degenerate_axis(direction),
                0.0001,
            );
            if moved.is_ok() {
                if direction == Singular::U {
                    solution.x = start.uv[0].x;
                }
                if direction == Singular::V {
                    solution.y = start.uv[0].y;
                }
            }
        }
        remove_duplicates(&mut solutions, tol_u, tol_v);
        start.uv[0] = solutions[0];
    }
    let on_side = |kind: PointType| matches!(kind, PointType::Start | PointType::End);
    let u_kind = natural.u_interval().classify_point(start.uv[0].x, tol_u);
    let v_kind = natural.v_interval().classify_point(start.uv[0].y, tol_v);
    let seam = if on_side(u_kind) && closed_u && on_side(v_kind) && closed_v {
        Closure::Both
    } else if on_side(u_kind) && closed_u {
        Closure::U
    } else if on_side(v_kind) && closed_v {
        Closure::V
    } else {
        Closure::Neither
    };
    if seam != Closure::Neither {
        let original = solutions.len();
        for ii in 0..original {
            start.uv[0] = solutions[ii];
            let u_kind = domain.u_interval().classify_point(start.uv[0].x, tol_u);
            let v_kind = domain.v_interval().classify_point(start.uv[0].y, tol_v);
            let mut along_u = if seam.has_u() {
                start.is_moving_along_boundary(surface, natural, tol_u, tol_v, 5.0, 0)?
            } else {
                Containment::Unknown
            };
            let mut along_v = if seam.has_v() {
                start.is_moving_along_boundary(surface, natural, tol_u, tol_v, 5.0, 1)?
            } else {
                Containment::Unknown
            };
            if along_u == Containment::Outside {
                solutions[ii].x = if u_kind == PointType::Start {
                    natural.max.x
                } else {
                    natural.min.x
                };
            }
            if along_v == Containment::Outside {
                solutions[ii].y = if v_kind == PointType::Start {
                    natural.max.y
                } else {
                    natural.min.y
                };
            }
            if along_u == Containment::Unknown && seam.has_u() {
                along_u = Containment::OnBoundary;
            }
            if along_v == Containment::Unknown && seam.has_v() {
                along_v = Containment::OnBoundary;
            }
            if along_u == Containment::OnBoundary {
                let copy = solutions[ii];
                solutions.push(copy);
                let jj = solutions.len() - 1;
                solutions[ii].x = natural.min.x;
                solutions[jj].x = natural.max.x;
            }
            if along_v == Containment::OnBoundary {
                let copy = solutions[ii];
                solutions.push(copy);
                let jj = solutions.len() - 1;
                solutions[ii].y = natural.min.y;
                solutions[jj].y = natural.max.y;
            }
        }
        remove_duplicates(&mut solutions, tol_u, tol_v);
        start.uv[0] = solutions[0];
    }
    let mut cache: Option<CurveCache> = None;
    for &solution in &solutions {
        start.t = interval.min;
        start.uv[0] = solution;
        let uv_tolerance = EFF_ZERO_PARAM
            * (1.0 + min(domain.u_interval().length(), domain.v_interval().length()));
        if !domain.contains_point(start.uv[0], uv_tolerance) {
            continue;
        }
        if let Some((iso, drop, deviation)) =
            drop_iso_curve(surface, domain, curve, interval, start.uv[0], tolerance)?
        {
            start.drop_to_surface = drop;
            start.approx_deviation = deviation;
            result.curves.push(UvCurve::Spline(iso));
            if max_drop < drop {
                max_drop = drop;
            }
            if max_deviation < deviation {
                max_deviation = deviation;
            }
            continue;
        }
        let mut growing = Growing::default();
        let mut have_start_tangent = false;
        if cache.is_none() {
            cache = Some(CurveCache::for_curve(curve, 0.5 * ZONE_TOL_3D)?);
        }
        let spans: Vec<Extent1> = cache
            .iter()
            .flat_map(|cache| cache.spans.iter().map(|span| span.interval))
            .collect();
        'spans: for &span in &spans {
            let mut span_interval = span;
            if span_interval.max <= start.t {
                continue;
            }
            if span_interval.min > start.t {
                let _ = span_interval.set_min_max(start.t, span_interval.max);
            }
            if interval.are_disjoint(span_interval, 0.0) {
                continue;
            }
            if !span_interval.is_contained_by(interval, EFF_ZERO) {
                let Ok(part) = span_interval.intersect(interval) else {
                    continue;
                };
                if approximate_length(curve, part, 5) < tolerance {
                    continue;
                }
                span_interval = part;
            }
            start.t = span_interval.min;
            end.t = span_interval.max;
            let derivatives = curve.evaluate(start.t, 2, true, true)?;
            start.curve = [derivatives[0], derivatives[1], derivatives[2]];
            if !have_start_tangent {
                let _ = start.set_properties(&walk, &mut offset, true)?;
                if start.drop_to_surface.abs() > max_drop {
                    max_drop = start.drop_to_surface.abs();
                }
                if start.approx_deviation > max_deviation {
                    max_deviation = start.approx_deviation;
                }
            }
            let end_t = end.t;
            let mut status = start.get_next_point(&walk, &mut offset, end_t, &mut end)?;
            if status == Status::Success {
                if end.drop_to_surface.abs() > max_drop {
                    max_drop = end.drop_to_surface.abs();
                }
                if end.approx_deviation > max_deviation {
                    max_deviation = end.approx_deviation;
                }
            }
            if (end.on_pole.truthy() || end.leaving.truthy())
                && !interval.is_value_on_boundary(end.t, EFF_ZERO).0
            {
                status = Status::OutsideOfDomain;
            }
            if status != Status::Success {
                if status == Status::OutsideOfDomain
                    && start.on_boundary == Tri::True
                    && start.leaving == Tri::True
                    && growing.knots.len() > 1
                {
                    let rest = Extent1::new(start.t, interval.max);
                    if approximate_length(curve, rest, 5) < tolerance {
                        let curve_end = curve.evaluate_point(interval.max)?;
                        if start.surface[0].distance(curve_end) < tolerance
                            && let Some(last) = growing.knots.last_mut()
                        {
                            *last = interval.max;
                        }
                        break 'spans;
                    }
                }
                if status == Status::NotWithinTolerance
                    && keep_all
                    && end.approx_deviation < tolerance
                {
                    status = Status::Success;
                    end.good_drop = Tri::True;
                }
                if status != Status::Success {
                    growing.clear();
                    break 'spans;
                }
            }
            // The stack of span ends still to reach: [t, u, v, du, dv].
            let mut stack: Vec<(f64, [Vector2; 2])> = vec![(end.t, end.uv)];
            while let Some(&(t, uv)) = stack.last() {
                end.t = t;
                end.uv = uv;
                start.fix_boundary_tangents(domain, &mut end, tol_u, tol_v);
                let (mut passes, drop, _, deviation) = start.test_span_accuracy(&walk, &end)?;
                if !passes && (end.t - start.t).abs() < interval.length() / 8000.0 {
                    passes = true;
                }
                if passes {
                    if max_drop < drop {
                        max_drop = drop;
                    }
                    if max_deviation < deviation {
                        max_deviation = deviation;
                    }
                    growing.add(&walk, &start, &end, &mut result.curves)?;
                    start = end;
                    start.from_left = [end.from_left[0].flipped(), end.from_left[1].flipped()];
                    have_start_tangent = true;
                    let (t, uv) = stack.pop().ok_or(Fault::Degenerate)?;
                    end.t = t;
                    end.uv = uv;
                    continue;
                }
                let mut found_one = false;
                for _ in 0..6 {
                    let end_t = (start.t + end.t) / 2.0;
                    status = start.get_next_point(&walk, &mut offset, end_t, &mut end)?;
                    if status == Status::Success {
                        found_one = true;
                        break;
                    }
                }
                if !found_one {
                    growing.clear();
                    break 'spans;
                }
                if max_drop < end.drop_to_surface {
                    max_drop = end.drop_to_surface;
                }
                if max_deviation < end.approx_deviation {
                    max_deviation = end.approx_deviation;
                }
                stack.push((end.t, end.uv));
            }
        }
        growing.flush(domain, tol_u, tol_v, &mut result.curves)?;
    }
    if result.curves.is_empty() {
        return Ok(result);
    }
    result.max_drop = max_drop;
    result.deviation = max_deviation;
    let mut measures = Vec::with_capacity(result.curves.len());
    for uv_curve in &result.curves {
        let UvCurve::Spline(uv_curve) = uv_curve else {
            return Err(Fault::Degenerate.into());
        };
        let (distance, deviation) =
            surface_curve_max_distance(surface, uv_curve, curve, tolerance)?;
        if distance > result.max_drop {
            result.max_drop = distance;
        }
        if deviation > result.deviation {
            result.deviation = deviation;
        }
        measures.push((distance, deviation));
    }
    if !keep_all
        && (result.max_drop > tolerance || result.deviation > tolerance)
        && result.curves.len() > 1
    {
        let (mut drop, mut deviation) = (0.0, 0.0);
        let mut better = Vec::new();
        for (uv_curve, (distance, curve_deviation)) in result.curves.drain(..).zip(measures) {
            if distance < tolerance {
                better.push(uv_curve);
                if drop < distance {
                    drop = distance;
                }
                if deviation < curve_deviation {
                    deviation = curve_deviation;
                }
            }
        }
        result.curves = better;
        result.max_drop = drop;
        result.deviation = deviation;
    }
    Ok(result)
}
