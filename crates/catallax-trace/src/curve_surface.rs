// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline curve and surface
// candidates of `spline.rs`; the local solvers and singularity tests take
// any `SolverSurface` (SMLib's virtual surface methods), such as an offset
// surface; statuses become `Fault`, output references become returned
// values, and solution containers become small result structs. The directed
// operation refuses with `SolverRefusal::NotPorted`.
// Translated from SMLib SmSurface.cpp (SmCSIntersectENFO,
// SmSurface::LocalCurveIntersect, SmCSNormalizeENFO, LocalCurveSolve,
// IsSingularity, GetSingularities, ApproxDerivativeLengths); rust-brep has no
// port.
//! The local curve-surface solvers of step 6: Newton iteration to a meeting
//! of a curve and a surface, or to their nearest or normal point pair, and
//! SMLib's surface singularity tests.

use crate::Fault;
use crate::curve_solve::PointSolveCurve;
use crate::extent::{Extent1, Extent2, Extent3};
use crate::global_solver::{SolverOperation, SolverRefusal};
use crate::local_solve::{EvalNFunctions, ExtentN, LocalSolveNd};
use crate::matrix::Matrix;
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, ZONE_TOL_3D, are_same, sqrt};
use crate::spline::{BSplineCurve, BSplineSurface};
use crate::surface_eval::{
    SolverSurface, evaluate, evaluate_normal, evaluate_point, first_derivatives, first_of,
};
use crate::vector::{Vector2, Vector3};

/// The SMLib singular-side flags (`SM_SS_UMIN` and so on).
pub const SINGULAR_U_MIN: u32 = 1;
/// `SM_SS_VMIN`.
pub const SINGULAR_V_MIN: u32 = 2;
/// `SM_SS_UMAX`.
pub const SINGULAR_U_MAX: u32 = 4;
/// `SM_SS_VMAX`.
pub const SINGULAR_V_MAX: u32 = 8;

/// The parameter direction in which a surface degenerates at a point.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Singular {
    /// `SM_SP_U`: the u derivative vanishes.
    U,
    /// `SM_SP_V`: the v derivative vanishes.
    V,
    /// `SM_SP_BOTH`.
    Both,
}

/// `Evaluate2ndDerivatives`: position, `Su`, `Sv`, `Suv`, `Suu`, `Svv`,
/// with SMLib's repair of vanishing tangents.
pub(crate) fn second_derivatives<S: SolverSurface + ?Sized>(
    surface: &S,
    uv: Vector2,
    sides: [bool; 2],
) -> Result<[Vector3; 6], Fault> {
    let d = surface.evaluate_cells(uv, 2, 2, sides, true, true, true)?;
    Ok([d[0], d[3], d[1], d[4], d[6], d[2]])
}

/// `ApproxDerivativeLengths`: the mean first-derivative lengths over a
/// four-by-four sample of `domain`.
pub fn approx_derivative_lengths<S: SolverSurface + ?Sized>(
    surface: &S,
    domain: Extent2,
) -> Vector2 {
    const SAMPLES: [f64; 4] = [0.15, 0.35, 0.65, 0.85];
    let (mut u, mut v) = (0.0, 0.0);
    for su in SAMPLES {
        for sv in SAMPLES {
            let uv = domain.evaluate(su, sv);
            let [_, du, dv] =
                first_of(surface, uv, [true, true]).unwrap_or([Vector3::default(); 3]);
            u += du.length();
            v += dv.length();
        }
    }
    let count = (SAMPLES.len() * SAMPLES.len()) as f64;
    Vector2::new(u / count, v / count)
}

/// `GetSingularities` with its approximate test on: the sides along which
/// the surface collapses to a point (within the scaled zero), as SMLib's
/// flags, and the pole points of the sides that are singular or nearly so,
/// in the order u low, v low, u high, v high.
pub fn get_singularities<S: SolverSurface + ?Sized>(
    surface: &S,
    zone: f64,
) -> (u32, [Option<Vector3>; 4]) {
    let domain = surface.normal_domain();
    let (u_interval, v_interval) = (domain.u_interval(), domain.v_interval());
    let approx_tolerance = 0.5 * zone;
    let mut size = approx_derivative_lengths(surface, domain);
    size.x *= domain.max.x - domain.min.x;
    size.y *= domain.max.y - domain.min.y;
    let mut surface_size = size.length();
    if surface_size < 1.0 {
        surface_size = 1.0;
    }
    let derivative_tolerance = surface_size * EFF_ZERO_SQRT * 10.0;
    let mut singular_flags = 0;
    let mut approx_flags = 0;
    const SAMPLES: usize = 5;
    for side in 0..4 {
        let (interval, along_v, side_size, flag) = match side {
            0 => (v_interval, true, size.y, SINGULAR_U_MIN),
            1 => (u_interval, false, size.x, SINGULAR_V_MIN),
            2 => (v_interval, true, size.y, SINGULAR_U_MAX),
            _ => (u_interval, false, size.x, SINGULAR_V_MAX),
        };
        let fixed = match side {
            0 => domain.min.x,
            1 => domain.min.y,
            2 => domain.max.x,
            _ => domain.max.y,
        };
        let length = interval.length();
        let mut edge = Extent3::default();
        let mut singular = true;
        let mut approx = true;
        for j in 0..SAMPLES {
            let parameter = if j == 0 {
                interval.min
            } else if j == SAMPLES - 1 {
                interval.max
            } else {
                interval.evaluate(j as f64 / (SAMPLES - 1) as f64)
            };
            let uv = if along_v {
                Vector2::new(fixed, parameter)
            } else {
                Vector2::new(parameter, fixed)
            };
            let d = surface
                .simple_cells(uv, 1, 1, [true, true], true)
                .unwrap_or_else(|_| vec![Vector3::default(); 4]);
            let speed = if along_v {
                d[1].length()
            } else {
                d[2].length()
            };
            if length * speed > approx_tolerance {
                approx = false;
                singular = false;
            } else if length * speed > derivative_tolerance {
                singular = false;
            }
            if !singular && !approx {
                break;
            }
            edge.add_point(d[0]);
        }
        if singular || approx {
            let extent = edge.size();
            let max_length = if extent.x > extent.y {
                if extent.x > extent.z {
                    extent.x
                } else {
                    extent.z
                }
            } else if extent.y > extent.z {
                extent.y
            } else {
                extent.z
            };
            approx &= max_length <= approx_tolerance;
            if approx {
                singular &= max_length <= EFF_ZERO * (1.0 + side_size);
            } else {
                singular = false;
            }
            if singular {
                singular_flags |= flag;
            }
            if approx {
                approx_flags |= flag;
            }
        }
    }
    let pole = |flag: u32, u: f64, v: f64| {
        ((singular_flags & flag) != 0 || (approx_flags & flag) != 0)
            .then(|| surface.point_at(domain.evaluate(u, v)).unwrap_or_default())
    };
    let poles = [
        pole(SINGULAR_U_MIN, 0.0, 0.5),
        pole(SINGULAR_V_MIN, 0.5, 0.0),
        pole(SINGULAR_U_MAX, 1.0, 0.5),
        pole(SINGULAR_V_MAX, 0.5, 1.0),
    ];
    (singular_flags, poles)
}

/// `IsSingularity`: whether a derivative vanishes at `uv` (or nearly, and
/// shrinking towards the nearer side), and otherwise, unless
/// `point_test_only`, whether the point lies within `tolerance` of a pole.
pub fn is_singularity<S: SolverSurface + ?Sized>(
    surface: &S,
    uv: Vector2,
    tolerance: f64,
    zone: f64,
    point_test_only: bool,
) -> Result<Option<Singular>, Fault> {
    // `SmPlane::IsSingularity`: never.
    if surface.plane_normal().is_some() {
        return Ok(None);
    }
    let d = surface.simple_cells(uv, 1, 1, [true, true], true)?;
    let (su, sv) = (d[2], d[1]);
    let (su2, sv2) = (su.length_squared(), sv.length_squared());
    let ratio_limit = ZONE_TOL_3D;
    let absolute_limit = EFF_ZERO * 100.0;
    let natural = surface.normal_domain();
    let mut singular = [false, false];
    for (axis, flag) in singular.iter_mut().enumerate() {
        let (own2, other2) = if axis == 0 { (su2, sv2) } else { (sv2, su2) };
        if own2 >= EFF_ZERO_SQRT {
            continue;
        }
        if own2 < EFF_ZERO_SQ {
            *flag = true;
            continue;
        }
        let size = natural.size();
        let (mut own, mut other) = (sqrt(own2), sqrt(other2));
        if axis == 0 {
            own *= size.x;
            other *= size.y;
        } else {
            other *= size.x;
            own *= size.y;
        }
        if own < other * ratio_limit || own < absolute_limit {
            let m = surface.evaluate_cells(uv, 2, 2, [true, true], true, true, true)?;
            let mut unit = if axis == 0 { m[3] } else { m[1] };
            let _ = unit.unitize();
            let cross = unit.dot(m[4]);
            let (value, middle) = if axis == 0 {
                (uv.y, natural.v_interval().mid())
            } else {
                (uv.x, natural.u_interval().mid())
            };
            if (value < middle && cross > own) || (value > middle && cross < -own) {
                *flag = true;
            }
        }
    }
    let mut direction = match singular {
        [true, true] => Some(Singular::Both),
        [true, false] => Some(Singular::U),
        [false, true] => Some(Singular::V),
        [false, false] => None,
    };
    if direction.is_none() && !point_test_only {
        let (flags, poles) = get_singularities(surface, zone);
        if flags != 0 {
            let point = surface.point_at(uv).unwrap_or_default();
            let near =
                |index: usize| poles[index].is_some_and(|pole| point.closer_than(tolerance, pole));
            if flags & SINGULAR_U_MIN != 0 && near(0) {
                direction = Some(Singular::V);
            } else if flags & SINGULAR_V_MIN != 0 && near(1) {
                direction = Some(Singular::U);
            } else if flags & SINGULAR_U_MAX != 0 && near(2) {
                direction = Some(Singular::V);
            } else if flags & SINGULAR_V_MAX != 0 && near(3) {
                direction = Some(Singular::U);
            }
        }
    }
    Ok(direction)
}

/// `SmCSIntersectENFO`: the gap between surface and curve points, switching
/// to the nearest-pair functions when convergence stalls near a tangency.
struct CurveSurfaceIntersection<'a, S: ?Sized> {
    curve: &'a dyn PointSolveCurve,
    surface: &'a S,
    count: u64,
    last_sum: f64,
    minimize: bool,
}

impl<S: SolverSurface + ?Sized> EvalNFunctions for CurveSurfaceIntersection<'_, S> {
    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        jacobian: Option<&mut Matrix>,
        found: &mut bool,
    ) -> Result<(), Fault> {
        self.count += 1;
        *found = false;
        let pv = self.curve.evaluate(x[0], 1, true, true)?;
        let uv = Vector2::new(x[1], x[2]);
        let [point, mut du, mut dv] = first_of(self.surface, uv, [true, true])?;
        if du.length_squared() < EFF_ZERO_SQ
            || dv.length_squared() < EFF_ZERO_SQ
            || pv[1].length_squared() < EFF_ZERO_SQ
        {
            return Err(Fault::Degenerate);
        }
        let mut gap = point - pv[0];
        let mut jacobian = jacobian;
        if let Some(j) = jacobian.as_mut() {
            j[0][0] = -pv[1].x;
            j[0][1] = du.x;
            j[0][2] = dv.x;
            j[1][0] = -pv[1].y;
            j[1][1] = du.y;
            j[1][2] = dv.y;
            j[2][0] = -pv[1].z;
            j[2][1] = du.z;
            j[2][2] = dv.z;
        }
        f[0] = gap.x;
        f[1] = gap.y;
        f[2] = gap.z;
        let iteration_limit = 10;
        let tangent_tolerance_deg = 0.5;
        let scaled_zero = EFF_ZERO * (1.0 + point.max_dimension());
        if self.minimize || (self.count > iteration_limit && gap.length_squared() > self.last_sum) {
            let pvv = self.curve.evaluate(x[0], 2, true, true)?;
            let normal = self.surface.normal_at(uv, [true, true]);
            if normal.is_perpendicular_to(pvv[1], tangent_tolerance_deg) {
                self.minimize = true;
                let [point, su, sv, duv, duu, dvv] =
                    second_derivatives(self.surface, uv, [true, true])?;
                (du, dv) = (su, sv);
                gap = point - pvv[0];
                if gap.length_squared() < 1000.0 * scaled_zero {
                    gap = if normal.dot(gap) >= 0.0 {
                        point + normal - pvv[0]
                    } else {
                        point - normal - pvv[0]
                    };
                }
                if let Some(j) = jacobian.as_mut() {
                    j[0][0] = pvv[2].dot(gap) - pvv[1].dot(pvv[1]);
                    j[0][1] = pvv[1].dot(du);
                    j[0][2] = pvv[1].dot(dv);
                    j[1][0] = -du.dot(pvv[1]);
                    j[1][1] = duu.dot(gap) + du.dot(du);
                    j[1][2] = duv.dot(gap) + du.dot(dv);
                    j[2][0] = -dv.dot(pvv[1]);
                    j[2][1] = j[1][2];
                    j[2][2] = dvv.dot(gap) + dv.dot(dv);
                }
                f[0] = pvv[1].dot(gap);
                f[1] = du.dot(gap);
                f[2] = dv.dot(gap);
            }
        }
        if self.count > 9 {
            self.last_sum = gap.length_squared();
        }
        Ok(())
    }
}

/// `SmCSNormalizeENFO` without the directed form: the gap dotted with the
/// unit curve tangent and unit surface derivatives.
struct CurveSurfaceNormal<'a, S: ?Sized> {
    curve: &'a dyn PointSolveCurve,
    surface: &'a S,
}

impl<S: SolverSurface + ?Sized> EvalNFunctions for CurveSurfaceNormal<'_, S> {
    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        jacobian: Option<&mut Matrix>,
        found: &mut bool,
    ) -> Result<(), Fault> {
        *found = false;
        let pvv = self.curve.evaluate(x[0], 2, true, true)?;
        let uv = Vector2::new(x[1], x[2]);
        let [s, su, sv, suv, suu, svv] = second_derivatives(self.surface, uv, [true, true])?;
        let (c, ct, ctt) = (pvv[0], pvv[1], pvv[2]);
        if ct.length_squared() < EFF_ZERO_SQ
            || su.length_squared() < EFF_ZERO_SQ
            || sv.length_squared() < EFF_ZERO_SQ
        {
            return Err(Fault::Degenerate);
        }
        let mut gap = s - c;
        let scaled_zero = EFF_ZERO * 1000.0 * (1.0 + s.max_dimension());
        if gap.length_squared() < scaled_zero {
            let normal = self.surface.normal_at(uv, [true, true]);
            gap = if normal.dot(gap) >= 0.0 {
                s + normal - c
            } else {
                s - normal - c
            };
        }
        let mut ct_unit = ct;
        let _ = ct_unit.unitize();
        let ct_unit_t = ct.unitized_derivative(ctt);
        let mut su_unit = su;
        let _ = su_unit.unitize();
        let su_unit_u = su.unitized_derivative(suu);
        let su_unit_v = su.unitized_derivative(suv);
        let mut sv_unit = sv;
        let _ = sv_unit.unitize();
        let sv_unit_u = sv.unitized_derivative(suv);
        let sv_unit_v = sv.unitized_derivative(svv);
        if let Some(j) = jacobian {
            j[0][0] = ct_unit_t.dot(gap) - ct_unit.dot(ct);
            j[0][1] = ct_unit.dot(su);
            j[0][2] = ct_unit.dot(sv);
            j[1][0] = -su_unit.dot(ct);
            j[1][1] = su_unit_u.dot(gap) + su_unit.dot(su);
            j[1][2] = su_unit_v.dot(gap) + su_unit.dot(sv);
            j[2][0] = -sv_unit.dot(ct);
            j[2][1] = sv_unit_u.dot(gap) + sv_unit.dot(su);
            j[2][2] = sv_unit_v.dot(gap) + sv_unit.dot(sv);
        }
        f[0] = ct_unit.dot(gap);
        f[1] = su_unit.dot(gap);
        f[2] = sv_unit.dot(gap);
        Ok(())
    }
}

/// A curve parameter and surface parameters with their distance
/// (`SmSolution`, single value, parameters t, u, v).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurveSurfacePoint {
    /// The curve parameter.
    pub t: f64,
    /// The surface parameters.
    pub uv: Vector2,
    /// The distance between the curve and surface points.
    pub distance: f64,
}

fn intervals(interval: Extent1, domain: Extent2) -> ExtentN {
    ExtentN {
        intervals: vec![interval, domain.u_interval(), domain.v_interval()],
    }
}

/// `SmSurface::LocalCurveSolve` for minimising, maximising and
/// normalising: Newton iteration on [`CurveSurfaceNormal`] from the guess.
pub fn local_curve_solve<S: SolverSurface + ?Sized>(
    surface: &S,
    domain: Extent2,
    curve: &dyn PointSolveCurve,
    interval: Extent1,
    operation: SolverOperation,
    guess: (f64, Vector2),
) -> Result<Option<CurveSurfacePoint>, SolverRefusal> {
    if !matches!(
        operation,
        SolverOperation::Minimize | SolverOperation::Maximize | SolverOperation::Normalize
    ) {
        return Err(SolverRefusal::NotPorted(operation));
    }
    let mut solver = LocalSolveNd::new(
        CurveSurfaceNormal { curve, surface },
        Some(intervals(interval, domain)),
        Some(vec![false; 3]),
    );
    solver.set_boundary_handler(crate::local_solve::BoundaryHandler::TotalBoundaryHits, 5);
    let root = solver.solve_it(&[guess.0, guess.1.x, guess.1.y], EFF_ZERO_SQRT)?;
    if !root.found {
        return Ok(None);
    }
    let p1 = curve.evaluate_point(root.x[0])?;
    let uv = Vector2::new(root.x[1], root.x[2]);
    let p2 = surface.point_at(uv)?;
    Ok(Some(CurveSurfacePoint {
        t: root.x[0],
        uv,
        distance: p1.distance(p2),
    }))
}

/// `SmSurface::LocalCurveIntersect`: a meeting of `curve` and `surface`
/// within `tolerance`, from Newton iteration on [`CurveSurfaceIntersection`]
/// and, where the curve runs nearly tangent to the surface, the nearest
/// pair; a surface pole is moved to its side of the domain.
pub fn local_curve_intersect<S: SolverSurface + ?Sized>(
    surface: &S,
    domain: Extent2,
    curve: &dyn PointSolveCurve,
    interval: Extent1,
    tolerance: f64,
    guess: (f64, Vector2),
    zone: f64,
) -> Result<Option<CurveSurfacePoint>, SolverRefusal> {
    let mut solver = LocalSolveNd::new(
        CurveSurfaceIntersection {
            curve,
            surface,
            count: 0,
            last_sum: 0.0,
            minimize: false,
        },
        Some(intervals(interval, domain)),
        Some(vec![false; 3]),
    );
    solver.set_desired_accuracy(EFF_ZERO / 10.0);
    solver.set_step_cutting_enabled(true);
    solver.set_boundary_handler(crate::local_solve::BoundaryHandler::TotalBoundaryHits, 6);
    solver.set_maximum_iterations(20);
    let root = solver.solve_it(&[guess.0, guess.1.x, guess.1.y], tolerance * 4.0)?;
    let (mut t, mut uv) = if root.found {
        (root.x[0], Vector2::new(root.x[1], root.x[2]))
    } else {
        guess
    };
    let pv = curve.evaluate(t, 1, true, true)?;
    let point = surface.point_at(uv)?;
    let mut distance = pv[0].distance(point);
    let normal = surface.normal_at(uv, [true, true]);
    let scaled_zero = EFF_ZERO * (1.0 + point.max_dimension());
    if normal.is_perpendicular_to(pv[1], 5.0)
        && let Some(nearest) = local_curve_solve(
            surface,
            domain,
            curve,
            interval,
            SolverOperation::Minimize,
            (t, uv),
        )?
        && (nearest.distance < distance || nearest.distance < scaled_zero)
    {
        distance = nearest.distance;
        t = nearest.t;
        uv = nearest.uv;
    }
    // SMLib's IsSingularity reports a failed evaluation as singular in no
    // direction, which changes nothing here.
    if let Ok(Some(direction)) = is_singularity(surface, uv, tolerance, zone, false) {
        let natural = surface.normal_domain();
        for axis in 0..2 {
            let hits = match direction {
                Singular::Both => true,
                Singular::U => axis == 0,
                Singular::V => axis == 1,
            };
            if !hits {
                continue;
            }
            let (sub, all) = if axis == 0 {
                (domain.u_interval(), natural.u_interval())
            } else {
                (domain.v_interval(), natural.v_interval())
            };
            let low_end = are_same(sub.min, all.min);
            let high_end = are_same(sub.max, all.max);
            let value = if axis == 0 { uv.x } else { uv.y };
            let side = if low_end != high_end {
                if low_end { sub.min } else { sub.max }
            } else if value < sub.mid() {
                sub.min
            } else {
                sub.max
            };
            if axis == 0 {
                uv.x = side;
            } else {
                uv.y = side;
            }
        }
    }
    Ok((distance < tolerance).then_some(CurveSurfacePoint { t, uv, distance }))
}

/// One end of a curve-surface coincidence range.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoincidenceEnd {
    /// The curve parameter.
    pub t: f64,
    /// The surface parameters.
    pub uv: Vector2,
    /// The distance from the surface at the last accepted sample of this side.
    pub distance: f64,
}

/// `SmSurface::SimpleCoincidenceChecker`: whether `curve` over `interval`
/// lies in `surface` within `tolerance`, tested at its middle and at pairs
/// of samples spreading to its ends; the sampled range where it does. With
/// `partial_ok` a partly coincident curve also answers, after three hits.
pub fn simple_coincidence_check(
    surface: &BSplineSurface,
    domain: Extent2,
    curve: &BSplineCurve,
    interval: Extent1,
    tolerance: f64,
    partial_ok: bool,
    owner_tolerance: Option<f64>,
) -> Result<Option<[CoincidenceEnd; 2]>, SolverRefusal> {
    use crate::surface_global::global_point_solve;
    use crate::surface_solve::local_point_solve;
    let (surface_box, _) = crate::surface_cache::bounding_boxes(surface, None);
    let curve_box = curve.bounding_box(interval)?;
    if surface_box.are_disjoint(curve_box, tolerance) {
        return Ok(None);
    }
    let mut sampled = Extent1::default();
    let mut hits = 0usize;
    let middle = interval.evaluate(0.5);
    let mid = curve.evaluate(middle, 1, true, true)?;
    let solutions = global_point_solve(
        surface,
        domain,
        SolverOperation::Minimize,
        mid[0],
        tolerance,
        owner_tolerance,
    )?;
    let Some(first) = solutions.as_slice().first().copied() else {
        return Ok(None);
    };
    let natural = curve.natural_interval();
    // On the surface and along it: the unit tangent is nearly perpendicular
    // to the surface normal.
    let along = |tangent: Vector3, uv: Vector2| {
        let mut tangent = tangent;
        let length = tangent.length();
        if length > EFF_ZERO {
            tangent /= length;
        }
        let normal = evaluate_normal(surface, uv, [false, false]);
        tangent.dot(normal).abs() < 0.017
    };
    let mid_uv = Vector2::new(first.start.parameters[0], first.start.parameters[1]);
    let mut found = first.start.value <= tolerance;
    if found {
        found = along(mid[1], mid_uv);
    }
    if found {
        sampled.add_value(middle);
        hits += 1;
    } else if !partial_ok {
        return Ok(None);
    }
    let mut guess_start = mid_uv;
    let mut guess_end = mid_uv;
    let (mut uv0, mut uv1) = (mid_uv, mid_uv);
    let mut start_distance = first.start.value;
    let mut end_distance = first.start.value;
    let count = curve.degree() + 5;
    let (mut have_start, mut have_end) = (false, false);
    let mut end_points_only = true;
    for ii in 0..count {
        let mut start_fraction = 0.5 - ((ii as f64 + 1.0) / count as f64) / 2.0;
        let mut end_fraction = 0.5 + ((ii as f64 + 1.0) / count as f64) / 2.0;
        if ii == count - 1 {
            start_fraction = 0.0;
            end_fraction = 1.0;
        }
        let start_t = interval.evaluate(start_fraction);
        let start = curve.evaluate(start_t, 1, true, true)?;
        let solution = local_point_solve(
            surface,
            domain,
            SolverOperation::Minimize,
            start[0],
            guess_start,
            owner_tolerance,
        )?;
        let mut found_start = solution.found && solution.value <= tolerance;
        if found_start {
            found_start = along(start[1], solution.uv);
        }
        if found_start {
            if !natural.is_value_on_boundary(start_t, EFF_ZERO).0 {
                end_points_only = false;
            }
            sampled.add_value(start_t);
            start_distance = solution.value;
            guess_start = solution.uv;
            let previous = hits;
            hits += 1;
            if previous >= 3 && partial_ok {
                break;
            }
        } else if !partial_ok {
            return Ok(None);
        }
        let end_t = interval.evaluate(end_fraction);
        let end = curve.evaluate(end_t, 1, false, true)?;
        let solution = local_point_solve(
            surface,
            domain,
            SolverOperation::Minimize,
            end[0],
            guess_end,
            owner_tolerance,
        )?;
        let mut found_end = solution.found && solution.value <= tolerance;
        if found_end {
            found_end = along(end[1], solution.uv);
        }
        if found_end {
            if !natural.is_value_on_boundary(end_t, EFF_ZERO).0 {
                end_points_only = false;
            }
            sampled.add_value(end_t);
            end_distance = solution.value;
            guess_end = solution.uv;
            let previous = hits;
            hits += 1;
            if previous >= 3 && partial_ok {
                break;
            }
        } else if !partial_ok {
            return Ok(None);
        }
        if found_start {
            uv0 = guess_start;
            have_start = true;
            if !have_end {
                uv1 = uv0;
                have_end = true;
            }
        }
        if found_end {
            uv1 = guess_end;
            have_end = true;
            if !have_start {
                uv0 = uv1;
                have_start = true;
            }
        }
    }
    if sampled.is_init() || sampled.length() < EFF_ZERO || end_points_only {
        return Ok(None);
    }
    Ok(Some([
        CoincidenceEnd {
            t: sampled.min,
            uv: uv0,
            distance: start_distance,
        },
        CoincidenceEnd {
            t: sampled.max,
            uv: uv1,
            distance: end_distance,
        },
    ]))
}

/// SMLib's second-order surface evaluation, `aEval[3][3]`, with the
/// normal in `[2][2]`.
pub(crate) type SurfaceJet = [[Vector3; 3]; 3];

/// What `smsurf_EvaluateGeometric` derives from a surface jet.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SurfaceGeometry {
    pub(crate) gaussian: f64,
    pub(crate) mean: f64,
    pub(crate) k1: f64,
    pub(crate) k2: f64,
    pub(crate) direction1: Vector3,
    pub(crate) direction2: Vector3,
}

/// `smsurf_EvaluateGeometric`: Gaussian and mean curvature, principal
/// curvatures and directions. Like SMLib it repairs a vanishing first
/// derivative of `jet` in place.
pub(crate) fn surface_geometry(jet: &mut SurfaceJet) -> Result<SurfaceGeometry, Fault> {
    use crate::geometry::solve_quadratic;
    let mut normal = jet[1][0].cross(jet[0][1]);
    if normal.length_squared() < EFF_ZERO_SQ {
        normal = jet[2][2];
        if normal.length_squared() < EFF_ZERO_SQ {
            return Err(Fault::Degenerate);
        }
    }
    normal.unitize()?;
    if jet[1][0].length_squared() < EFF_ZERO_SQ {
        jet[1][0] = jet[0][1].cross(normal);
    }
    if jet[0][1].length_squared() < EFF_ZERO_SQ {
        jet[0][1] = normal.cross(jet[1][0]);
    }
    let (du, dv, duv, duu, dvv) = (jet[1][0], jet[0][1], jet[1][1], jet[2][0], jet[0][2]);
    let (e, f, g) = (du.dot(du), du.dot(dv), dv.dot(dv));
    let (l, m, n) = (duu.dot(normal), duv.dot(normal), dvv.dot(normal));
    let numerator1 = l * n - m * m;
    let denominator1 = e * g - f * f;
    if denominator1.abs() <= EFF_ZERO * numerator1.abs() {
        return Err(Fault::Degenerate);
    }
    let gaussian = numerator1 / denominator1;
    let numerator2 = e * n - 2.0 * f * m + g * l;
    let denominator2 = 2.0 * (e * g - f * f);
    if denominator2.abs() <= EFF_ZERO * numerator2.abs() {
        return Err(Fault::Degenerate);
    }
    let mean = numerator2 / denominator2;
    let mut hhk = mean * mean - gaussian;
    if hhk < EFF_ZERO_SQ {
        hhk = 0.0;
    }
    let root = sqrt(hhk);
    let (k1, k2) = (mean + root, mean - root);
    let scaled_zero = EFF_ZERO_SQRT * (1.0 + k1.abs() + k2.abs());
    if (k1 - k2).abs() < scaled_zero {
        let mut direction1 = du;
        direction1.unitize()?;
        let mut direction2 = direction1.cross(dv).cross(direction1);
        direction2.unitize()?;
        return Ok(SurfaceGeometry {
            gaussian,
            mean,
            k1,
            k2,
            direction1,
            direction2,
        });
    }
    let (count, roots) = solve_quadratic([e * m - f * l, e * n - g * l, f * n - g * m], EFF_ZERO);
    let (mut direction1, mut direction2);
    if count < 1 {
        if du.max_dimension() > dv.max_dimension() {
            direction1 = du;
            direction2 = normal.cross(direction1);
        } else {
            direction1 = dv;
            direction2 = direction1.cross(normal);
        }
    } else {
        direction1 = du + roots[0] * dv;
        if direction1.length_squared() < EFF_ZERO_SQ {
            direction1 = du;
        }
        direction2 = normal.cross(direction1);
    }
    {
        let w = direction1;
        let (wu, wv) = (w.dot(du), w.dot(dv));
        let a = (g * wu - f * wv) / denominator1;
        let b = (-f * wu + e * wv) / denominator1;
        let i = e * a * a + 2.0 * f * a * b + g * b * b;
        if i > EFF_ZERO_SQ {
            let kn1 = (l * a * a + 2.0 * m * a * b + n * b * b) / i;
            if (kn1 - k2).abs() < (kn1 - k1).abs() {
                std::mem::swap(&mut direction1, &mut direction2);
            }
        }
    }
    if direction1.length_squared() > EFF_ZERO_SQ {
        direction1.unitize()?;
    }
    if direction2.length_squared() > EFF_ZERO_SQ {
        direction2.unitize()?;
    }
    if direction1.cross(direction2).dot(normal) < 0.0 {
        direction2 = -direction2;
    }
    Ok(SurfaceGeometry {
        gaussian,
        mean,
        k1,
        k2,
        direction1,
        direction2,
    })
}

/// `smsurf_EvaluateNormalSection`: the unit tangent-plane projection of
/// `direction` and the normal curvature along it.
pub(crate) fn normal_section_of(
    jet: &mut SurfaceJet,
    direction: Vector3,
) -> Result<(Vector3, f64), Fault> {
    let geometry = surface_geometry(jet)?;
    let (mut v1, mut v2) = (geometry.direction1, geometry.direction2);
    if v1.length_squared() > EFF_ZERO_SQ {
        v1.unitize()?;
    }
    if v2.length_squared() > EFF_ZERO_SQ {
        v2.unitize()?;
    }
    let mut projection = v1 * v1.dot(direction) + v2 * v2.dot(direction);
    if projection.length_squared() < EFF_ZERO_SQ {
        return Ok((Vector3::new(0.0, 0.0, 0.0), 0.0));
    }
    projection.unitize()?;
    let cos_theta = projection.dot(v1);
    let sin_theta = if cos_theta >= 1.0 || cos_theta <= -1.0 {
        0.0
    } else {
        sqrt(1.0 - cos_theta * cos_theta)
    };
    Ok((
        projection,
        geometry.k1 * cos_theta * cos_theta + geometry.k2 * sin_theta * sin_theta,
    ))
}

/// The jet `Evaluate(uv, 2, 2, ..., upper, ..., bNonZeroTangents = FALSE)`
/// in SMLib's `[3][3]` layout, with the unit normal in `[2][2]`.
pub(crate) fn surface_jet(
    surface: &BSplineSurface,
    uv: Vector2,
    sides: [bool; 2],
) -> Result<SurfaceJet, Fault> {
    let d = evaluate(surface, uv, 2, 2, sides, true, false, true)?;
    let mut jet = [[Vector3::default(); 3]; 3];
    for k in 0..3 {
        for l in 0..3 - k {
            jet[k][l] = d[k * 3 + l];
        }
    }
    let mut normal = jet[1][0].cross(jet[0][1]);
    if normal.length_squared() < EFF_ZERO_SQ {
        normal = evaluate_normal(surface, uv, sides);
    } else {
        normal.unitize()?;
    }
    jet[2][2] = normal;
    Ok(jet)
}

/// `SmSurface::EvaluateNormalSection`: the tangent-plane direction, normal
/// curvature and unit normal of `surface` at `uv` along `direction`.
pub(crate) fn normal_section(
    surface: &BSplineSurface,
    uv: Vector2,
    sides: [bool; 2],
    direction: Vector3,
) -> Result<(Vector3, f64, Vector3), Fault> {
    let mut jet = surface_jet(surface, uv, sides)?;
    let normal = jet[2][2];
    let (tangent, curvature) = normal_section_of(&mut jet, direction)?;
    Ok((tangent, curvature, normal))
}

/// `smsurf_DropVectors` for one vector: its coordinates in the basis of the
/// first derivatives. SMLib's fallback through a projection plane, for a
/// basis the direct solve rejects, is not carried.
pub(crate) fn drop_vector(
    du: Vector3,
    dv: Vector3,
    vector: Vector3,
) -> Result<Vector2, SolverRefusal> {
    use crate::geometry::vec_lin_comb_two_vectors;
    let (du2, dv2) = (du.length_squared(), dv.length_squared());
    if du2 < EFF_ZERO_SQ || dv2 < EFF_ZERO_SQ {
        return Ok(Vector2::new(0.0, 0.0));
    }
    let (mut u, mut v) = (du, dv);
    if du2 < EFF_ZERO {
        let _ = u.unitize();
    }
    if dv2 < EFF_ZERO {
        let _ = v.unitize();
    }
    let normal = u.cross(v);
    if normal.length_squared() < EFF_ZERO_SQ {
        return Err(Fault::Degenerate.into());
    }
    vec_lin_comb_two_vectors(du, dv, vector)
        .map_err(|_| SolverRefusal::NotPorted(SolverOperation::Intersect))
}

/// `SmSurface::DropPointsToTangentPlane` for one point.
pub(crate) fn drop_point_to_tangent_plane(
    surface: &BSplineSurface,
    uv: Vector2,
    sides: [bool; 2],
    point: Vector3,
) -> Result<Vector2, SolverRefusal> {
    let [p, du, dv] = first_derivatives(surface, uv, sides)?;
    if du.length_squared() < EFF_ZERO_SQ || dv.length_squared() < EFF_ZERO_SQ {
        return Ok(uv);
    }
    Ok(uv + drop_vector(du, dv, point - p)?)
}

/// A curve-surface solution as `sm_FilterSolution` reads it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Candidate {
    pub(crate) t: f64,
    pub(crate) uv: Vector2,
    pub(crate) value: f64,
}

/// `sm_FilterSolution`: whether to keep a near-miss solution. Close ones
/// stay; a nearest pair across a tangency stays by the curvatures on its
/// sides; others stay only against a discontinuous end or side they point
/// out of.
#[allow(clippy::too_many_arguments)]
pub(crate) fn filter_solution(
    surface: &BSplineSurface,
    domain: Extent2,
    curve: &BSplineCurve,
    interval: Extent1,
    u_conts: [crate::geometry::Continuity; 2],
    v_conts: [crate::geometry::Continuity; 2],
    t_conts: [crate::geometry::Continuity; 2],
    candidate: Candidate,
) -> Result<bool, SolverRefusal> {
    use crate::curve::evaluate_geometric;
    use crate::geometry::Continuity;
    if candidate.value < EFF_ZERO {
        return Ok(true);
    }
    let pv = curve.evaluate(candidate.t, 1, true, true)?;
    let random = evaluate_point(surface, domain.evaluate(0.456, 0.6789))?;
    let scaled_zero =
        10.0 * EFF_ZERO * (1.0 + crate::scalar::max(pv[0].max_dimension(), random.max_dimension()));
    if candidate.value < scaled_zero {
        return Ok(true);
    }
    let uv = candidate.uv;
    let [point, du, dv] = first_derivatives(surface, uv, [true, true])?;
    let mut normal = du.cross(dv);
    if normal.length_squared() < EFF_ZERO_SQ {
        normal = evaluate_normal(surface, uv, [true, true]);
    }
    normal.unitize()?;
    let difference = pv[0] - point;
    let at_min = are_same(interval.min, candidate.t);
    let at_max = are_same(interval.max, candidate.t);
    if normal.is_parallel_to(difference, 5.0).0 && pv[1].is_perpendicular_to(difference, 5.0) {
        if (at_min && t_conts[0] >= Continuity::G1G2) || (at_max && t_conts[1] >= Continuity::G1G2)
        {
            let (_, surface_curvature, surface_normal) =
                normal_section(surface, uv, [true, true], pv[1]).unwrap_or((
                    Vector3::default(),
                    0.0,
                    Vector3::default(),
                ));
            let frenet =
                evaluate_geometric(curve, candidate.t, 2, true).unwrap_or([Vector3::default(); 4]);
            let mut curve_curvature = frenet[2].length();
            let curve_line = curve_curvature.abs() < EFF_ZERO_SQRT;
            let surface_line = surface_curvature.abs() < EFF_ZERO_SQRT;
            if curve_line && surface_line {
                return Ok(true);
            } else if curve_line {
                // SMLib keeps the point unless the product is positive.
                #[allow(clippy::neg_cmp_op_on_partial_ord)]
                let keep = !(surface_curvature * difference.dot(surface_normal) > 0.0);
                return Ok(keep);
            } else if surface_line {
                return Ok(difference.dot(frenet[2]) > 0.0);
            }
            curve_curvature = frenet[2].dot(surface_normal);
            let positive_gap = difference.dot(surface_normal) > 0.0;
            let positive_ks = surface_curvature > 0.0;
            let positive_kc = curve_curvature > 0.0;
            let kc_bigger = curve_curvature.abs() > surface_curvature.abs();
            let keep = match (positive_gap, positive_ks, positive_kc) {
                (false, true, true) | (true, false, false) => !kc_bigger,
                (false, true, false) | (true, false, true) => true,
                (false, false, true) | (true, true, false) => false,
                (false, false, false) | (true, true, true) => kc_bigger,
            };
            return Ok(keep);
        }
        return Ok(true);
    }
    let on_u_min = are_same(domain.min.x, uv.x);
    let on_u_max = are_same(domain.max.x, uv.x);
    let on_v_min = are_same(domain.min.y, uv.y);
    let on_v_max = are_same(domain.max.y, uv.y);
    let curve_dot = pv[1].dot(difference);
    let u_dot = difference.dot(du);
    let v_dot = difference.dot(dv);
    let mut curve_save = (t_conts[0] <= Continuity::C0 && at_min && curve_dot > 0.0)
        || (t_conts[1] <= Continuity::C0 && at_max && curve_dot < 0.0);
    let mut u_save = (u_conts[0] <= Continuity::C0 && on_u_min && u_dot < 0.0)
        || (u_conts[1] <= Continuity::C0 && on_u_max && u_dot > 0.0);
    let mut v_save = (v_conts[0] <= Continuity::C0 && on_v_min && v_dot < 0.0)
        || (v_conts[1] <= Continuity::C0 && on_v_max && v_dot > 0.0);
    let mut keep = curve_save || u_save || v_save;
    if keep && curve_save {
        if crate::curve_solve::is_closed(curve, curve.natural_interval(), 0.0) {
            curve_save = false;
        }
        keep = curve_save || u_save || v_save;
    }
    let natural = surface.natural_domain();
    if keep && u_save {
        if crate::surface_solve::is_closed(surface, natural, 0, None) {
            u_save = false;
        }
        keep = curve_save || u_save || v_save;
    }
    if keep && v_save {
        if crate::surface_solve::is_closed(surface, natural, 1, None) {
            v_save = false;
        }
        keep = curve_save || u_save || v_save;
    }
    Ok(keep)
}

/// `sm_ConvertToParameterDelta`.
fn parameter_delta(tangent: Vector3, delta: Vector3) -> f64 {
    let length = tangent.length();
    if length < EFF_ZERO {
        return 0.0;
    }
    let unit = tangent / length;
    unit.dot(delta) / length
}

/// `SmTsectCurveType`: how a curve meets a surface at a point.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum MeetingType {
    /// `SM_TC_CROSSING`.
    Crossing,
    /// `SM_TC_TOUCHING`.
    Touching,
    /// `SM_TC_TANGENT`.
    Tangent,
    /// `SM_TC_NEAR_TANGENT`.
    NearTangent,
    /// `SM_TC_COINCIDENT`.
    Coincident,
}

/// `sm_ClassifyCurveSurfaceXSectPoint`: near a tangency, compares the
/// curve's and the surface section's osculating circles to classify the
/// meeting and, where they suggest a second meeting nearby, a guess for it.
pub(crate) fn classify_meeting(
    surface: &BSplineSurface,
    domain: Extent2,
    curve: &BSplineCurve,
    interval: Extent1,
    uv: Vector2,
    t: f64,
) -> Result<(MeetingType, Option<(f64, Vector2)>), SolverRefusal> {
    use crate::curve::evaluate_geometric;
    use crate::geometry::{circle_coplanar_circle_intersect, line_coplanar_circle_intersect};
    use crate::scalar::{PI, asin, cos};
    let mut jet = surface_jet(surface, uv, [true, true])?;
    let surface_normal = jet[2][2];
    let Ok(frenet) = evaluate_geometric(curve, t, 2, true) else {
        return Ok((MeetingType::Touching, None));
    };
    if !surface_normal.is_perpendicular_to(frenet[1], 20.0) {
        return Ok((MeetingType::Crossing, None));
    }
    let surface_point = evaluate_point(surface, uv)?;
    let (section_tangent, section_curvature) = normal_section_of(&mut jet, frenet[1])?;
    let geometry = surface_geometry(&mut jet)?;
    #[derive(PartialEq)]
    enum Local {
        Planar,
        Elliptical,
        Hyperbolic,
        Cylindrical,
    }
    let local = if geometry.mean.abs() < EFF_ZERO_SQRT
        && geometry.k1.abs() < EFF_ZERO_SQRT
        && geometry.k2.abs() < EFF_ZERO_SQRT
    {
        Local::Planar
    } else if geometry.gaussian > EFF_ZERO_SQRT {
        Local::Elliptical
    } else if geometry.gaussian < -EFF_ZERO_SQRT {
        Local::Hyperbolic
    } else {
        Local::Cylindrical
    };
    let mut curve_normal = frenet[2];
    let curve_curvature = curve_normal.length();
    if curve_curvature.abs() > EFF_ZERO_SQRT {
        curve_normal /= curve_curvature;
    } else {
        curve_normal = surface_normal;
    }
    let mut binormal = frenet[1].cross(curve_normal);
    if binormal.length_squared() < EFF_ZERO_SQ {
        return Ok((MeetingType::Touching, None));
    }
    binormal.unitize()?;
    let angle = binormal.angle_between(surface_normal).unwrap_or(0.0);
    let cos_phi = cos(angle - PI / 2.0).abs();
    let mut curve_line;
    let mut curve_radius = 0.0;
    let (mut curve_center, mut curve_circle_normal) = (Vector3::default(), Vector3::default());
    if curve_curvature.abs() > EFF_ZERO_SQRT {
        curve_line = false;
        curve_radius = (1.0 / curve_curvature).abs();
        curve_center = frenet[0] + curve_normal * curve_radius;
        curve_circle_normal = binormal;
    } else {
        curve_line = true;
    }
    let (curve_point, curve_vector) = (frenet[0], frenet[1]);
    let (mut surface_line, mut two_lines) = (false, false);
    let mut surface_radius = 0.0;
    let (mut surface_center, mut surface_circle_normal) = (Vector3::default(), Vector3::default());
    let mut surface_vector: Vector3;
    let mut surface_vector2 = Vector3::default();
    const EFF_ZERO_RAD: f64 = 1.0e-4 * PI / 180.0;
    if cos_phi.abs() < EFF_ZERO_RAD {
        match local {
            Local::Planar => return Ok((MeetingType::Coincident, None)),
            Local::Elliptical => return Ok((MeetingType::Tangent, None)),
            Local::Hyperbolic => {
                surface_line = true;
                two_lines = true;
                let y = 10000.0;
                let mut tmp = (1.0 - geometry.k1 * y * y) / geometry.k2;
                if tmp < 0.0 {
                    tmp = 0.0;
                }
                let x = sqrt(tmp);
                surface_vector = y * geometry.direction1 + x * geometry.direction2;
                surface_vector2 = y * geometry.direction1 - x * geometry.direction2;
                let _ = surface_vector.unitize();
                let _ = surface_vector2.unitize();
            }
            Local::Cylindrical => {
                surface_line = true;
                surface_vector = if geometry.k1.abs() < EFF_ZERO_SQRT {
                    geometry.direction1
                } else {
                    geometry.direction2
                };
            }
        }
    } else if section_curvature.abs() < EFF_ZERO_SQRT {
        surface_line = true;
        surface_vector = section_tangent;
    } else {
        surface_vector = section_tangent;
        surface_radius = (1.0 / section_curvature).abs() * cos_phi;
        let sign = if section_curvature >= 0.0 { 1.0 } else { -1.0 };
        let mut center_vector = binormal.cross(sign * surface_normal).cross(binormal);
        let _ = center_vector.unitize();
        surface_center = surface_point + center_vector * surface_radius;
        surface_circle_normal = binormal;
    }
    let curve_tolerance = if curve_line {
        EFF_ZERO_SQRT
    } else {
        curve_radius / 1000.0
    };
    let surface_tolerance = if surface_line {
        EFF_ZERO_SQRT
    } else {
        surface_radius / 1000.0
    };
    let tolerance = curve_tolerance + surface_tolerance;
    if curve_line && surface_line {
        let along = (two_lines && curve_vector.is_parallel_to(surface_vector2, 2.0).0)
            || curve_vector.is_parallel_to(surface_vector, 2.0).0;
        let kind = if along {
            MeetingType::Coincident
        } else {
            MeetingType::Crossing
        };
        return Ok((kind, None));
    }
    if !curve_line && !surface_line {
        if curve_radius > 1000.0 * surface_radius {
            curve_line = true;
        } else if surface_radius > 1000.0 * curve_radius {
            surface_line = true;
        }
    }
    let mut found = Vector3::default();
    if curve_line {
        let hits = line_coplanar_circle_intersect(
            curve_point,
            curve_vector,
            surface_center,
            surface_circle_normal,
            surface_radius,
            tolerance,
        )?;
        if hits.count < 2 {
            return Ok((MeetingType::Tangent, None));
        }
        let p = hits.parameters;
        let line_t = if p[1].abs() > p[0].abs() { p[1] } else { p[0] };
        found = line_t * curve_vector;
    }
    if surface_line {
        let near = |vector: Vector3| -> Result<Option<Vector3>, SolverRefusal> {
            let hits = line_coplanar_circle_intersect(
                surface_point,
                vector,
                curve_center,
                curve_circle_normal,
                curve_radius,
                tolerance,
            )?;
            if hits.count < 2 {
                return Ok(None);
            }
            let p = hits.parameters;
            let line_t = if p[1].abs() > p[0].abs() { p[1] } else { p[0] };
            // SMLib scales the first line's vector for both lines.
            Ok(Some(line_t * surface_vector))
        };
        let first = near(surface_vector)?;
        let second = if two_lines {
            near(surface_vector2)?
        } else {
            None
        };
        match (first, second) {
            (None, None) => return Ok((MeetingType::Tangent, None)),
            (Some(a), Some(b)) => {
                found = if b.length_squared() < a.length_squared() {
                    b
                } else {
                    a
                }
            }
            (Some(a), None) => found = a,
            (None, Some(b)) => found = b,
        }
    }
    if !curve_line && !surface_line {
        let (count, points) = if curve_radius > surface_radius {
            circle_coplanar_circle_intersect(
                curve_circle_normal,
                curve_center,
                curve_radius,
                surface_center,
                surface_radius,
                tolerance,
            )?
        } else {
            circle_coplanar_circle_intersect(
                surface_circle_normal,
                surface_center,
                surface_radius,
                curve_center,
                curve_radius,
                tolerance,
            )?
        };
        match count {
            2 => {}
            3 => return Ok((MeetingType::Coincident, None)),
            _ => return Ok((MeetingType::Tangent, None)),
        }
        found = points[0] - surface_point;
        found = if found.length() < points[1].distance(surface_point) {
            points[1] - surface_point
        } else {
            points[0] - surface_point
        };
    }
    let (du, dv) = (jet[1][0], jet[0][1]);
    let sign = if frenet[1].dot(found) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    let step = if curve_line {
        found.length()
    } else {
        2.0 * curve_radius * asin(found.length() / 2.0 / curve_radius)
    };
    let new_t = t + sign * step / frenet[1].length();
    let new_uv = match crate::geometry::vec_lin_comb_two_vectors(du, dv, found) {
        Ok(delta) => uv + delta,
        Err(_) => Vector2::new(
            uv.x + parameter_delta(du, found),
            uv.y + parameter_delta(dv, found),
        ),
    };
    if !interval.contains_value(new_t, 0.0) && !domain.contains_point(new_uv, 0.0) {
        return Ok((MeetingType::NearTangent, None));
    }
    Ok((
        MeetingType::NearTangent,
        Some((interval.clamp_value(new_t), domain.clamp_point(new_uv))),
    ))
}
