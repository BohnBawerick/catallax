// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline surface candidates of
// `spline.rs` (the iso-curve also over any `SolverSurface`); statuses become
// `Fault`, output references become returned values, the lazily evaluated
// `SmSurfaceEval` becomes a cached state, and the solution container becomes
// a small result struct. Paths that need code from later port steps refuse
// with `SolverRefusal::NotPorted` instead of running: the older solver SMLib
// keeps for maximising
// (`LocalPointSolve_0`), the directed operation, and the curve drop of the
// inversion near a singular point (`InvertPointNearSingularity`) where the
// curve solver itself refuses.
// Translated from SMLib SmSurfacePointSolve.cpp (SmSurfPtSolveStepper,
// SmSurface::LocalPointSolve_1), SmSurface.cpp (SmSurfaceEval,
// SmSurface::LocalPointSolve, IsOnSeam, IsOnBoundary, IsClosed,
// EvaluateDirectionalDerivs, InvertPointNearSingularity,
// FindDegenParamForDirection, SmFindSingularParamWithGivenTangentDirEFO,
// smsurf_IsPointOnNormal), SmTol.cpp (MapTo2d), SmIsoCurve.cpp, SmCurve.cpp
// (IsDegenerate) and SmGeomUtility.cpp (smgu_DirectionalDerivs); rust-brep
// has no port.
//! The local surface point solver of step 5: Newton steps from a guess to
//! the nearest point of a B-spline surface, with SMLib's plane steps, line
//! searches and boundary handling.

use crate::Fault;
use crate::curve_global::drop_point;
use crate::curve_solve::{PointSolveCurve, local_point_solve as curve_point_solve};
use crate::extent::{Extent1, Extent2};
use crate::geometry::{
    frenet_from_derivatives, line_closest_point, line_point_distance, vec_lin_comb_two_vectors,
};
use crate::global_solver::{SolverOperation, SolverRefusal};
use crate::local_solve::{EvalFunction, LocalSolve1d, Root1d, TerminationReason};
use crate::scalar::{
    BIG_DOUBLE, EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, ZONE_TOL_3D, max, min, sqrt,
};
use crate::spline::BSplineSurface;
use crate::surface_eval::{
    SolverSurface, evaluate, evaluate_normal, evaluate_out_of_bounds, evaluate_point,
    first_derivatives,
};
use crate::vector::{Vector2, Vector3};

/// `SM_EFF_ZERO_PARAM`: zero for knot parameter values.
const EFF_ZERO_PARAM: f64 = 1.2e-8;

/// A constant-parameter curve of a surface (`SmIsoCurve`): `axis` 0 fixes
/// `u`, axis 1 fixes `v`. Evaluation follows the surface's out-of-bounds
/// switch, which SMLib turns on for the whole local point solve.
pub struct IsoCurve<'a> {
    surface: &'a BSplineSurface,
    axis: usize,
    parameter: f64,
    from_left: bool,
    out_of_bounds: bool,
}

impl<'a> IsoCurve<'a> {
    /// `SmIsoCurve(surface, param, value, bFromLeft)`, clamping `parameter`
    /// into the natural domain as SMLib does by default.
    pub fn new(surface: &'a BSplineSurface, axis: usize, parameter: f64, from_left: bool) -> Self {
        let domain = surface.natural_domain();
        let interval = if axis == 0 {
            domain.u_interval()
        } else {
            domain.v_interval()
        };
        let mut value = parameter;
        if interval.min > parameter {
            value = interval.min;
        }
        if parameter > interval.max {
            value = interval.max;
        }
        Self {
            surface,
            axis,
            parameter: value,
            from_left,
            out_of_bounds: false,
        }
    }

    /// `GetNaturalInterval`: the surface interval along the curve.
    pub fn natural_interval(&self) -> Extent1 {
        let domain = self.surface.natural_domain();
        if self.axis == 0 {
            domain.v_interval()
        } else {
            domain.u_interval()
        }
    }

    fn uv(&self, t: f64) -> Vector2 {
        if self.axis == 0 {
            Vector2::new(self.parameter, t)
        } else {
            Vector2::new(t, self.parameter)
        }
    }
}

impl PointSolveCurve for IsoCurve<'_> {
    fn evaluate(
        &self,
        t: f64,
        derivatives: usize,
        from_left: bool,
        nonzero_tangent: bool,
    ) -> Result<Vec<Vector3>, Fault> {
        let (orders, sides) = if self.axis == 0 {
            ([0, derivatives], [self.from_left, from_left])
        } else {
            ([derivatives, 0], [from_left, self.from_left])
        };
        let run = if self.out_of_bounds {
            evaluate_out_of_bounds
        } else {
            evaluate
        };
        run(
            self.surface,
            self.uv(t),
            orders[0],
            orders[1],
            sides,
            false,
            nonzero_tangent,
            true,
        )
    }

    fn evaluate_point(&self, t: f64) -> Result<Vector3, Fault> {
        evaluate_point(self.surface, self.uv(t))
    }

    fn is_degenerate(&self, tolerance: f64) -> bool {
        curve_is_degenerate(self, self.natural_interval(), tolerance)
    }
}

/// `SmCurve::IsDegenerate` over `interval`: the ends meet and the tangent
/// vanishes at both ends and at an interior sample.
fn curve_is_degenerate(curve: &dyn PointSolveCurve, interval: Extent1, tolerance: f64) -> bool {
    if interval.min == interval.max {
        return true;
    }
    let p1 = curve.evaluate_point(interval.min).unwrap_or_default();
    let p2 = curve.evaluate_point(interval.max).unwrap_or_default();
    let scaled_zero = EFF_ZERO * (1.0 + p1.max_dimension());
    let tolerance2 = max(tolerance * tolerance, scaled_zero * scaled_zero);
    if p1.distance_squared(p2) < tolerance2 {
        let length_squared = interval.length() * interval.length();
        let mut values = [Vector3::default(); 2];
        for t in [interval.min, interval.max, interval.evaluate(0.4567)] {
            if let Ok(v) = curve.evaluate(t, 1, true, false) {
                values = [v[0], v[1]];
            }
            if values[1].length_squared() * length_squared > tolerance2 {
                return false;
            }
        }
        // SMLib's negated test: a NaN distance counts as degenerate.
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        let degenerate = !(p1.distance_squared(values[0]) > tolerance2);
        return degenerate;
    }
    false
}

/// The interval of `domain` along `axis` (0 for `u`).
fn side_interval(domain: Extent2, axis: usize) -> Extent1 {
    if axis == 0 {
        domain.u_interval()
    } else {
        domain.v_interval()
    }
}

/// `SmIsoCurve` of any surface the local solvers take, such as an offset
/// surface: `axis` 0 fixes `u`, axis 1 fixes `v`. Unlike [`IsoCurve`] it has
/// no out-of-bounds switch.
pub struct SurfaceIsoCurve<'a, S: ?Sized> {
    surface: &'a S,
    axis: usize,
    parameter: f64,
    from_left: bool,
}

impl<'a, S: SolverSurface + ?Sized> SurfaceIsoCurve<'a, S> {
    /// `SmIsoCurve(surface, param, value, bFromLeft, bClampInput)`: with
    /// `clamp`, `parameter` is clamped into the natural domain.
    pub fn new(surface: &'a S, axis: usize, parameter: f64, from_left: bool, clamp: bool) -> Self {
        let interval = side_interval(surface.normal_domain(), axis);
        let mut value = parameter;
        if clamp {
            if interval.min > parameter {
                value = interval.min;
            }
            if parameter > interval.max {
                value = interval.max;
            }
        }
        Self {
            surface,
            axis,
            parameter: value,
            from_left,
        }
    }

    /// `GetNaturalInterval`: the surface interval along the curve.
    pub fn natural_interval(&self) -> Extent1 {
        side_interval(self.surface.normal_domain(), 1 - self.axis)
    }

    fn uv(&self, t: f64) -> Vector2 {
        if self.axis == 0 {
            Vector2::new(self.parameter, t)
        } else {
            Vector2::new(t, self.parameter)
        }
    }
}

impl<S: SolverSurface + ?Sized> PointSolveCurve for SurfaceIsoCurve<'_, S> {
    fn evaluate(
        &self,
        t: f64,
        derivatives: usize,
        from_left: bool,
        nonzero_tangent: bool,
    ) -> Result<Vec<Vector3>, Fault> {
        let (orders, sides) = if self.axis == 0 {
            ([0, derivatives], [self.from_left, from_left])
        } else {
            ([derivatives, 0], [from_left, self.from_left])
        };
        self.surface.evaluate_cells(
            self.uv(t),
            orders[0],
            orders[1],
            sides,
            false,
            nonzero_tangent,
            true,
        )
    }

    fn evaluate_point(&self, t: f64) -> Result<Vector3, Fault> {
        self.surface.point_at(self.uv(t))
    }

    fn is_degenerate(&self, tolerance: f64) -> bool {
        curve_is_degenerate(self, self.natural_interval(), tolerance)
    }
}

impl crate::curve_solve::SolverCurve for IsoCurve<'_> {
    fn natural_interval(&self) -> Extent1 {
        IsoCurve::natural_interval(self)
    }

    /// `SmCurve::GetDegree`, which `SmIsoCurve` does not override.
    fn degree(&self) -> usize {
        3
    }

    /// `SmCurve::GetNumberNaturalKnots`, which `SmIsoCurve` does not override.
    fn natural_knot_count(&self) -> usize {
        2
    }

    /// `SmIsoCurve::GetKnots`: the surface's distinct knots along the curve.
    fn knot_values(&self) -> Vec<f64> {
        let knots = self.surface.knots()[1 - self.axis];
        crate::spline::unique_knots(knots, None)
            .iter()
            .map(|k| k.0)
            .collect()
    }
}

/// `SmSurfParamType` values that SMLib's boundary tests return.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Bound {
    Neither,
    Both,
    UMin,
    VMin,
    UMax,
    VMax,
}

/// The problem shared by every iterate (`SmSurfacePointSolver`).
pub(crate) struct Problem<'a> {
    pub(crate) surface: &'a BSplineSurface,
    pub(crate) domain: Extent2,
    pub(crate) test: Vector3,
    pub(crate) operation: SolverOperation,
    pub(crate) tight: f64,
    pub(crate) loose: f64,
}

/// SMLib's error norms (`SmSrfPtErrorNormType`).
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Norm {
    SumAbs,
    F,
    G,
    FAbs,
    GAbs,
    Gap,
}

/// The evaluated iterate (`SmSurfaceEval` and the stepper's derived values).
#[derive(Clone, Copy)]
struct State {
    position: Vector3,
    su: Vector3,
    sv: Vector3,
    suu: Vector3,
    suv: Vector3,
    svv: Vector3,
    gap: Vector3,
    gap_length: f64,
    f_check: f64,
    g_check: f64,
    sum_check: f64,
}

/// One iterate of the local solve (`SmSurfPtSolveStepper`).
#[derive(Clone)]
pub(crate) struct Stepper<'a> {
    problem: &'a Problem<'a>,
    uv: Vector2,
    state: Option<State>,
    normal: Option<Vector3>,
    iterations: usize,
    hits: [i32; 2],
    max_step: f64,
    suspect: bool,
}

pub(crate) type Solve<T> = Result<T, SolverRefusal>;

impl<'a> Stepper<'a> {
    pub(crate) fn new(problem: &'a Problem<'a>, uv: Vector2) -> Self {
        let counts = problem.surface.counts();
        Self {
            problem,
            uv,
            state: None,
            normal: None,
            iterations: 0,
            hits: [0, 0],
            max_step: if counts == [2, 2] { 2.0 } else { 0.40 },
            suspect: false,
        }
    }

    pub(crate) fn set_uv(&mut self, uv: Vector2) {
        self.uv = uv;
        self.state = None;
        self.normal = None;
    }

    fn clamp(&mut self) {
        self.set_uv(self.problem.domain.clamp_point(self.uv));
    }

    /// `Eval`: position, derivatives, gap and the scaled residuals.
    fn eval(&mut self) -> Solve<State> {
        if let Some(state) = self.state {
            return Ok(state);
        }
        let d = evaluate_out_of_bounds(
            self.problem.surface,
            self.uv,
            2,
            2,
            [true, true],
            true,
            true,
            true,
        )?;
        let (position, sv, svv, su, suv, suu) = (d[0], d[1], d[2], d[3], d[4], d[6]);
        let gap = position - self.problem.test;
        let f = gap.dot(su);
        let g = gap.dot(sv);
        let mut su_length = su.length();
        let mut sv_length = sv.length();
        let minimum = 0.0001;
        if su_length < minimum {
            su_length = minimum;
        }
        if sv_length < minimum {
            sv_length = minimum;
        }
        let distance = 1.0;
        let f_check = f / (distance * su_length);
        let g_check = g / (distance * sv_length);
        let state = State {
            position,
            su,
            sv,
            suu,
            suv,
            svv,
            gap,
            gap_length: gap.length(),
            f_check,
            g_check,
            sum_check: f_check.abs() + g_check.abs(),
        };
        self.state = Some(state);
        Ok(state)
    }

    pub(crate) fn gap_length(&mut self) -> Solve<f64> {
        Ok(self.eval()?.gap_length)
    }

    /// `SmSurfaceEval::EvalNormal`: the unit cross product, then the
    /// cross-derivative fallback at a degenerate side, then `EvaluateNormal`.
    fn normal(&mut self) -> Solve<Vector3> {
        if let Some(normal) = self.normal {
            return Ok(normal);
        }
        let s = self.eval()?;
        let mut normal = s.su.cross(s.sv);
        if normal.length_squared() < EFF_ZERO_SQ {
            let domain = self.problem.surface.natural_domain();
            if s.sv.length_squared() < s.su.length_squared() {
                normal = s.su.cross(s.suv);
                if self.uv.y > domain.v_interval().mid() {
                    normal *= -1.0;
                }
            } else {
                normal = s.suv.cross(s.sv);
                if self.uv.x > domain.u_interval().mid() {
                    normal *= -1.0;
                }
            }
            if normal.length_squared() < EFF_ZERO_SQ {
                normal = evaluate_normal(self.problem.surface, self.uv, [true, true]);
                self.normal = Some(normal);
                return Ok(normal);
            }
        }
        normal.unitize()?;
        self.normal = Some(normal);
        Ok(normal)
    }

    pub(crate) fn error_norm(&mut self, norm: Norm) -> Solve<f64> {
        let s = self.eval()?;
        Ok(match norm {
            Norm::SumAbs => s.sum_check,
            Norm::F => s.f_check,
            Norm::G => s.g_check,
            Norm::FAbs => s.f_check.abs(),
            Norm::GAbs => s.g_check.abs(),
            Norm::Gap => s.gap_length,
        })
    }

    /// `IsSingular`: one derivative much shorter than the other and tiny.
    fn singular(&mut self) -> Solve<Option<usize>> {
        let s = self.eval()?;
        let ratio = 10000.0 * 10000.0;
        let su2 = s.su.dot(s.su);
        let sv2 = s.sv.dot(s.sv);
        let mut which = None;
        if su2 < sv2 / ratio && su2 < EFF_ZERO_PARAM {
            which = Some(0);
        }
        if sv2 < su2 / ratio && sv2 < EFF_ZERO_PARAM {
            which = Some(1);
        }
        Ok(which)
    }

    /// `Converged(double)`: on the surface, or on the surface normal through
    /// the test point.
    fn converged(&mut self, tolerance: f64) -> Solve<bool> {
        let s = self.eval()?;
        if s.gap_length < tolerance {
            return Ok(true);
        }
        if self.singular()?.is_none() {
            return Ok(s.sum_check < tolerance);
        }
        let normal = self.normal()?;
        if normal.length_squared() < EFF_ZERO_SQ {
            return Ok(false);
        }
        let t = line_closest_point(s.position, normal, self.problem.test);
        let line_point = s.position + t * normal;
        let difference = line_point - self.problem.test;
        let mut to_surface = line_point.distance(s.position);
        if to_surface < 1.0 {
            to_surface = 1.0;
        }
        Ok(difference.length() < tolerance * to_surface)
    }

    fn converged_tight(&mut self, tight: bool) -> Solve<bool> {
        let tolerance = if tight {
            self.problem.tight
        } else {
            self.problem.loose
        };
        self.converged(tolerance)
    }

    /// `PushingBoundary`: on a domain side with the gap pointing out of it.
    fn pushing_boundary(&mut self, tolerance: f64) -> Solve<Bound> {
        let s = self.eval()?;
        let gap = self.problem.test - s.position;
        let Ok(step) = vec_lin_comb_two_vectors(s.su, s.sv, gap) else {
            return Ok(Bound::Neither);
        };
        let domain = self.problem.domain;
        let mut result = Bound::Neither;
        if self.uv.x < domain.u_interval().min + EFF_ZERO && step.x < -tolerance {
            result = Bound::UMin;
        } else if self.uv.x > domain.u_interval().max - EFF_ZERO && step.x > tolerance {
            result = Bound::UMax;
        }
        if self.uv.y < domain.v_interval().min + EFF_ZERO && step.y < -tolerance {
            result = if result == Bound::Neither {
                Bound::VMin
            } else {
                Bound::Both
            };
        } else if self.uv.y > domain.v_interval().max - EFF_ZERO && step.y > tolerance {
            result = if result == Bound::Neither {
                Bound::VMax
            } else {
                Bound::Both
            };
        }
        Ok(result)
    }

    /// `ConvergedOnBoundary`: pushing a side and stationary along it.
    fn converged_on_boundary(&mut self, tolerance: f64) -> Solve<Option<Bound>> {
        let bound = self.pushing_boundary(-EFF_ZERO)?;
        if bound == Bound::Neither {
            return Ok(None);
        }
        let domain = self.problem.domain;
        if matches!(bound, Bound::UMin | Bound::UMax | Bound::Both) {
            if self.error_norm(Norm::GAbs)? < tolerance {
                return Ok(Some(bound));
            }
            let interval = domain.v_interval();
            if self.uv.y < interval.min + tolerance {
                return Ok((self.error_norm(Norm::G)? > -tolerance).then_some(bound));
            } else if self.uv.y > interval.max - tolerance {
                return Ok((self.error_norm(Norm::G)? < tolerance).then_some(bound));
            }
        }
        if matches!(bound, Bound::VMin | Bound::VMax | Bound::Both) {
            if self.error_norm(Norm::FAbs)? < tolerance {
                return Ok(Some(bound));
            }
            let interval = domain.u_interval();
            if self.uv.x < interval.min + tolerance {
                return Ok((self.error_norm(Norm::F)? > -tolerance).then_some(bound));
            } else if self.uv.x > interval.max - tolerance {
                return Ok((self.error_norm(Norm::F)? < tolerance).then_some(bound));
            }
        }
        Ok(None)
    }

    /// `SnapToBoundaries`: moves onto sides within `tolerance` when that
    /// does not worsen the error.
    fn snap_to_boundaries(&mut self, tolerance: f64) -> Solve<bool> {
        let domain = self.problem.domain;
        let mut uv = self.uv;
        if !domain.point_on_boundary(uv, tolerance).0 || domain.point_on_boundary(uv, 0.0).0 {
            return Ok(false);
        }
        let norm = if self.problem.operation == SolverOperation::Normalize {
            Norm::SumAbs
        } else {
            Norm::Gap
        };
        let current = self.error_norm(norm)?;
        let saved = uv;
        if uv.x <= domain.min.x + tolerance {
            uv.x = domain.min.x;
        } else if uv.x >= domain.max.x - tolerance {
            uv.x = domain.max.x;
        }
        if uv.y <= domain.min.y + tolerance {
            uv.y = domain.min.y;
        } else if uv.y >= domain.max.y - tolerance {
            uv.y = domain.max.y;
        }
        self.set_uv(uv);
        if self.error_norm(norm)? > current {
            self.set_uv(saved);
            return Ok(false);
        }
        Ok(true)
    }

    /// `HitBound`: counts consecutive hits on one side of `axis`.
    fn hit_bound(&mut self, axis: usize, side: i32) {
        if side == 0 {
            if self.hits[axis] != 0 {
                let interval = if axis == 0 {
                    self.problem.domain.u_interval()
                } else {
                    self.problem.domain.v_interval()
                };
                let value = if axis == 0 { self.uv.x } else { self.uv.y };
                if interval.contains_value(value, -EFF_ZERO) {
                    self.hits[axis] = 0;
                }
            }
        } else {
            self.hits[axis] = if self.hits[axis] * side > 0 {
                self.hits[axis] + side
            } else {
                side
            };
        }
    }

    /// `Step`: a Newton step, or with `plane` the step that ignores the
    /// second-derivative terms, from `previous` into `self`.
    fn step(&mut self, previous: &mut Stepper<'a>, plane: bool) -> Solve<TerminationReason> {
        self.eval()?;
        let p = previous.eval()?;
        if previous.singular()?.is_some() {
            let mut uv = previous.uv;
            let problem = self.problem;
            let natural = problem.surface.natural_domain();
            if invert_point_near_singularity(problem, natural, &mut uv)? {
                self.set_uv(uv);
                return Ok(TerminationReason::FoundAnswerConverged);
            }
        }
        let mu = p.su.length();
        let mv = p.sv.length();
        if mu < EFF_ZERO || mv < EFF_ZERO {
            return Ok(TerminationReason::UnableToConverge);
        }
        let su_n = p.su.divided(mu);
        let sv_n = p.sv.divided(mv);
        let f = p.gap.dot(su_n);
        let g = p.gap.dot(sv_n);
        let su_nu = p.su.unitized_derivative(p.suu);
        let su_nv = p.su.unitized_derivative(p.suv);
        let sv_nu = p.sv.unitized_derivative(p.suv);
        let sv_nv = p.sv.unitized_derivative(p.svv);
        let gap_u = p.su;
        let gap_v = p.sv;
        let fu1 = gap_u.dot(su_n);
        let mut fu2 = p.gap.dot(su_nu);
        let limit = 0.90;
        let operation = self.problem.operation;
        if fu2 < -fu1 * limit
            && matches!(
                operation,
                SolverOperation::Minimize | SolverOperation::Intersect
            )
        {
            self.suspect = true;
            fu2 = -fu1 * limit;
        }
        let gv1 = gap_v.dot(sv_n);
        let mut gv2 = p.gap.dot(sv_nv);
        if gv2 < -gv1 * limit
            && matches!(
                operation,
                SolverOperation::Minimize | SolverOperation::Intersect
            )
        {
            self.suspect = true;
            gv2 = -gv1 * limit;
        }
        let (fu, gv, fv, gu) = if plane {
            (fu1, gv1, gap_v.dot(su_n), gap_u.dot(sv_n))
        } else {
            (
                fu1 + fu2,
                gv1 + gv2,
                gap_v.dot(su_n) + p.gap.dot(su_nv),
                gap_u.dot(sv_n) + p.gap.dot(sv_nu),
            )
        };
        let det = fu * gv - gu * fv;
        let numerator_u = -(f * gv - g * fv);
        let numerator_v = -(fu * g - gu * f);
        if det.abs() <= EFF_ZERO * numerator_u.abs() || det.abs() <= EFF_ZERO * numerator_v.abs() {
            self.suspect = true;
            return Ok(TerminationReason::UnableToConverge);
        }
        let step = Vector2::new(numerator_u / det, numerator_v / det);
        self.set_uv(previous.uv + step);
        Ok(TerminationReason::StillIterating)
    }

    /// `CheckBounds`: whether the step from `current` to `next` leaves the
    /// domain, with the per-axis side hit, the clamped end, the end scaled
    /// back along the step, and that scale.
    fn check_bounds(
        current: Vector2,
        next: Vector2,
        domain: Extent2,
    ) -> ([i32; 2], Vector2, Vector2, f64) {
        let mut hits = [0, 0];
        let mut scale = 1.0;
        let (u_min, v_min, u_max, v_max) = (domain.min.x, domain.min.y, domain.max.x, domain.max.y);
        if next.x < u_min {
            hits[0] = -1;
            scale = if current.x >= u_min {
                (u_min - current.x) / (next.x - current.x)
            } else {
                0.0
            };
        } else if next.x > u_max {
            hits[0] = 1;
            scale = if current.x <= u_max {
                (u_max - current.x) / (next.x - current.x)
            } else {
                0.0
            };
        }
        if next.y < v_min {
            hits[1] = -1;
            let this_scale = if current.y >= v_min {
                (v_min - current.y) / (next.y - current.y)
            } else {
                0.0
            };
            if this_scale < scale {
                scale = this_scale;
            }
        } else if next.y > v_max {
            hits[1] = 1;
            let this_scale = if current.y <= v_max {
                (v_max - current.y) / (next.y - current.y)
            } else {
                0.0
            };
            if this_scale < scale {
                scale = this_scale;
            }
        }
        let mut clamped = Vector2::default();
        let mut scaled = Vector2::default();
        if hits != [0, 0] {
            scaled = current + scale * (next - current);
            clamped = next;
            if hits[0] < 0 {
                clamped.x = u_min;
            } else if hits[0] > 0 {
                clamped.x = u_max;
            }
            if hits[1] < 0 {
                clamped.y = v_min;
            } else if hits[1] > 0 {
                clamped.y = v_max;
            }
        }
        (hits, clamped, scaled, scale)
    }

    /// `AdjustStep`: limits the step from `previous` to a fraction of the
    /// domain and clips it to the domain. `Ok(None)` is SMLib's `SM_ERR`;
    /// otherwise whether the step changed and its final scale.
    fn adjust_step(&mut self, previous: &Stepper<'a>) -> Solve<Option<(bool, f64)>> {
        let domain = self.problem.domain;
        let mut step = self.uv - previous.uv;
        let original = step.length();
        if original < EFF_ZERO_SQ {
            if domain.contains_point(self.uv, 0.0) {
                return Ok(Some((false, 1.0)));
            }
            self.set_uv(domain.clamp_point(self.uv));
            return Ok(None);
        }
        let mut adjusted = false;
        let mut size_scale = 1.0;
        let max_u = self.max_step * domain.u_interval().length();
        if step.x.abs() > max_u {
            adjusted = true;
            size_scale = max_u / step.x.abs();
            step.x *= size_scale;
        }
        let max_v = self.max_step * domain.v_interval().length();
        if step.y.abs() > max_v {
            adjusted = true;
            let scale_y = min(max_v / step.y.abs(), size_scale);
            step.y *= scale_y;
            size_scale = min(size_scale, scale_y);
        }
        let current = previous.uv;
        let next = current + step;
        let (hits, clamped, scaled, bounds_scale) = Self::check_bounds(current, next, domain);
        self.hit_bound(0, hits[0]);
        self.hit_bound(1, hits[1]);
        let out = hits != [0, 0];
        if !out {
            self.set_uv(next);
        } else if bounds_scale > 1.0 - EFF_ZERO || bounds_scale < EFF_ZERO {
            self.set_uv(clamped);
        } else {
            self.set_uv(clamped);
            let mut scaled_eval = self.clone();
            scaled_eval.set_uv(scaled);
            if scaled_eval.error_norm(Norm::SumAbs)? < self.error_norm(Norm::SumAbs)? {
                *self = scaled_eval;
            }
        }
        let mut scale = 1.0;
        if size_scale < 1.0 || out {
            scale = (self.uv - previous.uv).length() / original;
            adjusted = true;
        }
        Ok(Some((adjusted, scale)))
    }

    /// `LineSearch`: a better point between `previous` and `self`.
    fn line_search(&mut self, previous: &mut Stepper<'a>) -> Solve<()> {
        if previous.error_norm(Norm::F)? * self.error_norm(Norm::F)? < 0.0
            && previous.error_norm(Norm::FAbs)? > 10.0 * previous.error_norm(Norm::GAbs)?
            && self.error_norm(Norm::FAbs)? > 10.0 * self.error_norm(Norm::GAbs)?
        {
            return self.line_search_bracketed(previous, Norm::F);
        }
        if previous.error_norm(Norm::G)? * self.error_norm(Norm::G)? < 0.0
            && previous.error_norm(Norm::GAbs)? > 10.0 * previous.error_norm(Norm::FAbs)?
            && self.error_norm(Norm::GAbs)? > 10.0 * self.error_norm(Norm::FAbs)?
        {
            return self.line_search_bracketed(previous, Norm::G);
        }
        let norm = if self.problem.operation == SolverOperation::Normalize {
            Norm::SumAbs
        } else {
            Norm::Gap
        };
        let mut third = self.clone();
        let error0 = previous.error_norm(norm)?;
        let error1 = self.error_norm(norm)?;
        let mut best_error = error0;
        let mut best_value = 0.0;
        if error1 < best_error {
            best_error = error1;
            best_value = 1.0;
        }
        let best_initial = best_error;
        for fraction in [0.1, 0.2, 0.4, 0.7] {
            third.set_uv((1.0 - fraction) * previous.uv + fraction * self.uv);
            let error = third.error_norm(norm)?;
            if error < best_error {
                if error < best_initial * 0.50 {
                    *self = third;
                    return Ok(());
                }
                best_error = error;
                best_value = fraction;
            }
        }
        if best_value < 0.01 && error0 < error1 {
            let fraction = -1.0;
            third.set_uv((1.0 - fraction) * previous.uv + fraction * self.uv);
            third.clamp();
            let error = third.error_norm(norm)?;
            if error < best_error {
                best_error = error;
                best_value = fraction;
            }
        }
        if best_value > 0.99 && error1 < error0 {
            let fraction = 2.0;
            third.set_uv((1.0 - fraction) * previous.uv + fraction * self.uv);
            third.clamp();
            let error = third.error_norm(norm)?;
            if error < best_error {
                best_value = fraction;
            }
        }
        let uv = (1.0 - best_value) * previous.uv + best_value * self.uv;
        self.set_uv(uv);
        self.clamp();
        Ok(())
    }

    /// `LineSearchBracketted`: regula falsi on a signed residual whose
    /// sign differs between `previous` and `self`.
    fn line_search_bracketed(&mut self, previous: &mut Stepper<'a>, norm: Norm) -> Solve<()> {
        let mut y0 = previous.error_norm(norm)?;
        let mut y1 = self.error_norm(norm)?;
        if y0 * y1 > 0.0 {
            return Ok(());
        }
        let best_input = min(y0.abs(), y1.abs());
        let mut max_iterations = 5;
        let (mut x0, mut x1) = (0.0, 1.0);
        let mut middle = previous.clone();
        let mut i = 0;
        while i < max_iterations {
            let mut x_mid = x0 - y0 * (x1 - x0) / (y1 - y0);
            middle.set_uv((1.0 - x_mid) * previous.uv + x_mid * self.uv);
            middle.set_uv(self.problem.domain.clamp_point(middle.uv));
            let mut y_mid = middle.error_norm(norm)?;
            if y_mid.abs() < EFF_ZERO_SQRT {
                *self = middle;
                return Ok(());
            }
            if i >= max_iterations - 1 && y_mid.abs() >= best_input {
                max_iterations += 5;
                if max_iterations > 30 {
                    break;
                }
                x_mid = (x0 + x1) / 2.0;
                middle.set_uv((1.0 - x_mid) * previous.uv + x_mid * self.uv);
                y_mid = middle.error_norm(norm)?;
            }
            if y_mid * y0 < 0.0 {
                x1 = x_mid;
                y1 = y_mid;
            } else {
                x0 = x_mid;
                y0 = y_mid;
            }
            i += 1;
        }
        *self = middle;
        Ok(())
    }

    /// `FindBoundarySolution`: when pushing a side, the nearest point on
    /// that side's iso-curve; true when it is the constrained optimum.
    fn find_boundary_solution(&mut self) -> Solve<bool> {
        let problem = self.problem;
        let tight = problem.tight;
        let mut bound = self.pushing_boundary(-tight)?;
        if bound == Bound::Neither {
            return Ok(false);
        }
        self.snap_to_boundaries(problem.loose)?;
        if bound == Bound::Both {
            if self.pushing_boundary(tight)? != Bound::Neither {
                return Ok(true);
            }
            bound = if self.uv.x < problem.domain.u_interval().mid() {
                Bound::UMin
            } else {
                Bound::UMax
            };
        }
        let natural = problem.surface.natural_domain();
        let (u, v) = (problem.domain.u_interval(), problem.domain.v_interval());
        let (axis, side, guess, iso_domain, cross_domain) = match bound {
            Bound::UMin => (0, natural.u_interval().clamp_value(u.min), self.uv.y, v, u),
            Bound::UMax => (0, natural.u_interval().clamp_value(u.max), self.uv.y, v, u),
            Bound::VMin => (1, natural.v_interval().clamp_value(v.min), self.uv.x, u, v),
            Bound::VMax => (1, natural.v_interval().clamp_value(v.max), self.uv.x, u, v),
            _ => return Ok(false),
        };
        let mut new_uv = self.uv;
        let low_side = side < cross_domain.mid();
        let mut curve = IsoCurve::new(problem.surface, axis, side, low_side);
        curve.out_of_bounds = true;
        let found = match curve_point_solve(
            &curve,
            iso_domain,
            SolverOperation::Minimize,
            problem.test,
            guess,
        ) {
            Ok(Some(found)) => found,
            _ => return Ok(false),
        };
        let mut boundary = self.clone();
        let (cross1, cross2);
        if axis == 0 {
            new_uv.y = found.parameter;
            boundary.set_uv(new_uv);
            let s = boundary.eval()?;
            (cross1, cross2) = (s.su, s.suu);
        } else {
            new_uv.x = found.parameter;
            boundary.set_uv(new_uv);
            let s = boundary.eval()?;
            (cross1, cross2) = (s.sv, s.svv);
        }
        let norm = if axis == 0 { Norm::GAbs } else { Norm::FAbs };
        if boundary.error_norm(norm)? < self.error_norm(norm)? {
            *self = boundary;
        }
        let gap = self.eval()?.gap;
        let flip = match problem.operation {
            SolverOperation::Maximize => true,
            SolverOperation::Normalize => cross1.dot(cross1) + gap.dot(cross2) < 0.0,
            _ => false,
        };
        let mut error = gap.dot(cross1);
        if flip {
            error = -error;
        }
        Ok(if low_side {
            error > -tight
        } else {
            error < tight
        })
    }

    /// `FindNextStep`: a Newton step, a plane step when Newton is suspect,
    /// boundary solutions, line searches and the choice between them.
    /// Returns SMLib's status as `false` for `SM_ERR_NOT_CONVERGING` or
    /// `SM_ERR_OUTSIDE_OF_DOMAIN`.
    fn find_next_step(
        &mut self,
        previous: &mut Stepper<'a>,
        reason: &mut TerminationReason,
    ) -> Solve<bool> {
        *reason = TerminationReason::StillIterating;
        let tight = self.problem.tight;
        let bound = self.pushing_boundary(tight)?;
        if self.iterations > 0 || bound != Bound::Neither {
            if bound == Bound::Both {
                *reason = TerminationReason::OutOfBounds;
                return Ok(true);
            }
            if self.find_boundary_solution()? {
                *reason = TerminationReason::FoundAnswerConverged;
                if self.pushing_boundary(tight)? != Bound::Neither {
                    *reason = TerminationReason::OutOfBounds;
                }
                return Ok(true);
            }
        }
        let improvement = 0.75;
        self.suspect = false;
        let mut newton = previous.clone();
        *reason = newton.step(previous, false)?;
        if *reason == TerminationReason::FoundAnswerConverged {
            *self = newton;
            return Ok(true);
        }
        let mut got_newton = *reason == TerminationReason::StillIterating;
        if got_newton {
            match newton.adjust_step(previous)? {
                Some((adjusted, scale)) => {
                    if newton.converged_tight(true)? {
                        *self = newton;
                        return Ok(true);
                    }
                    if adjusted && scale < EFF_ZERO_SQ {
                        self.suspect = true;
                    }
                }
                None => {
                    self.suspect = true;
                    got_newton = false;
                }
            }
        } else {
            self.suspect = true;
        }
        let mut suspect = self.suspect || newton.suspect;
        if !suspect {
            let step = newton.uv - self.uv;
            let factor = self.max_step / 2.0;
            if step.x.abs() > factor * self.problem.domain.u_interval().length()
                || step.y.abs() > factor * self.problem.domain.v_interval().length()
            {
                suspect = true;
            }
        }
        if !suspect
            && newton.error_norm(Norm::SumAbs)? > previous.error_norm(Norm::SumAbs)? * improvement
            && newton.error_norm(Norm::SumAbs)? > self.problem.loose
        {
            suspect = true;
        }
        self.suspect = suspect;
        if got_newton && !self.suspect {
            *self = newton;
            return Ok(true);
        }
        self.suspect = true;
        let mut plane = previous.clone();
        *reason = plane.step(previous, true)?;
        if *reason == TerminationReason::FoundAnswerConverged {
            *self = newton;
            return Ok(true);
        }
        let mut got_plane = *reason == TerminationReason::StillIterating;
        if got_plane {
            match plane.adjust_step(previous)? {
                Some(_) => {
                    if plane.converged_tight(true)? {
                        *self = plane;
                        return Ok(true);
                    }
                }
                None => got_plane = false,
            }
        }
        if !got_newton && !got_plane {
            *self = previous.clone();
            if newton.hits != [0, 0] {
                *reason = TerminationReason::OutOfBounds;
                return Ok(true);
            }
            *reason = TerminationReason::UnableToConverge;
            return Ok(false);
        }
        let norm = if self.problem.operation == SolverOperation::Normalize {
            Norm::SumAbs
        } else {
            Norm::Gap
        };
        if got_plane
            && plane.find_boundary_solution()?
            && plane.error_norm(norm)? <= previous.error_norm(norm)?
            && (!got_newton || plane.error_norm(norm)? <= newton.error_norm(norm)?)
        {
            *self = plane;
            return self.settle_on_boundary(reason);
        }
        if got_newton
            && newton.find_boundary_solution()?
            && newton.error_norm(norm)? < previous.error_norm(norm)? + EFF_ZERO
            && (!got_plane || newton.error_norm(norm)? < plane.error_norm(norm)? + EFF_ZERO)
        {
            *self = newton;
            return self.settle_on_boundary(reason);
        }
        if got_newton {
            newton.line_search(previous)?;
        }
        if got_plane {
            plane.line_search(previous)?;
        }
        if got_newton && !got_plane {
            *self = newton;
        } else if !got_newton {
            *self = plane;
        } else if newton.error_norm(norm)? < plane.error_norm(norm)? {
            *self = newton;
        } else {
            *self = plane;
        }
        if self.hits[0].abs() > 3 || self.hits[1].abs() > 3 {
            *reason = TerminationReason::OutOfBounds;
            return Ok(true);
        }
        if self.error_norm(norm)? >= previous.error_norm(norm)?
            || (self.uv - previous.uv).length() < EFF_ZERO_SQ
        {
            if self.hits[0].abs() > 1 || self.hits[1].abs() > 1 {
                *reason = TerminationReason::OutOfBounds;
                return Ok(true);
            }
            *reason = TerminationReason::UnableToConverge;
            return Ok(false);
        }
        Ok(true)
    }

    /// The end of `FindNextStep` after a boundary solution was taken.
    fn settle_on_boundary(&mut self, reason: &mut TerminationReason) -> Solve<bool> {
        if self.pushing_boundary(self.problem.tight)? != Bound::Neither {
            *reason = TerminationReason::OutOfBounds;
            return Ok(false);
        }
        *reason = TerminationReason::FoundAnswerConverged;
        Ok(true)
    }
}

/// `smsurf_IsPointOnNormal` with a domain: whether `test` lies on the
/// surface normal at `uv` (clamped into `domain`), and its distance along
/// that normal.
fn is_point_on_normal(
    surface: &BSplineSurface,
    domain: Extent2,
    uv: Vector2,
    test: Vector3,
    tolerance: f64,
) -> (bool, f64) {
    let clamped = domain.clamp_point(uv);
    let Ok([point, du, dv]) = first_derivatives(surface, clamped, [true, true]) else {
        return (false, 0.0);
    };
    let mut normal = du.cross(dv);
    if normal.length_squared() < EFF_ZERO_SQ {
        normal = evaluate_normal(surface, clamped, [true, true]);
    } else {
        let _ = normal.unitize();
    }
    let (_, t) = line_point_distance(point, normal, test);
    let difference = (point + t * normal) - test;
    let mut to_surface = t.abs();
    if to_surface < 1.0 {
        to_surface = 1.0;
    }
    if difference.length() < tolerance * to_surface {
        (true, t.abs())
    } else {
        (false, 0.0)
    }
}

/// `SmFindSingularParamWithGivenTangentDirEFO`: along the degenerate
/// parameter `s` at fixed `t`, the derivative of the unit cross tangent
/// projected on `direction`.
struct SingularTangent<'a> {
    surface: &'a BSplineSurface,
    uv: Vector2,
    direction: Vector3,
    singular: usize,
}

impl EvalFunction for SingularTangent<'_> {
    fn evaluate(
        &mut self,
        s: f64,
        value: &mut f64,
        derivative: &mut f64,
        found: &mut bool,
        _signal_errors: bool,
    ) -> Result<(), Fault> {
        *value = 0.0;
        *derivative = 0.0;
        *found = false;
        let point = if self.singular == 0 {
            Vector2::new(s, self.uv.y)
        } else {
            Vector2::new(self.uv.x, s)
        };
        let m = evaluate(self.surface, point, 3, 3, [true, true], true, true, true)?;
        let (w, ws, wss) = if self.singular == 0 {
            (m[1], m[5], m[9])
        } else {
            (m[4], m[5], m[6])
        };
        let (unit_t, unit_tt) = w.unitized_derivative2(ws, wss);
        *value = unit_t.dot(self.direction);
        *derivative = unit_tt.dot(self.direction);
        Ok(())
    }
}

/// `SmSurface::FindDegenParamForDirection`: moves the degenerate parameter
/// of `uv` to where the cross tangent points along `direction` (from the
/// pole when `from_pole`), sampling ten guesses before a local solve.
fn find_degenerate_parameter(
    surface: &BSplineSurface,
    uv: &mut Vector2,
    direction: Vector3,
    from_pole: bool,
    singular: usize,
    tolerance: f64,
) -> Result<(), Fault> {
    let mut direction = direction;
    direction.unitize()?;
    let natural = surface.natural_domain();
    let (s, t, s_interval, t_interval) = if singular == 0 {
        (uv.x, uv.y, natural.u_interval(), natural.v_interval())
    } else {
        (uv.y, uv.x, natural.v_interval(), natural.u_interval())
    };
    if (t > t_interval.mid() && from_pole) || (t <= t_interval.mid() && !from_pole) {
        direction = -direction;
    }
    let f_tolerance = 2.0 / 50.0;
    let mut sample = s_interval;
    let increment = (sample.max - sample.min) * -0.025;
    sample.min -= increment;
    sample.max += increment;
    let (mut near_best_value, mut near_best_parameter) = (-1.0, 0.5);
    let (mut best_value, mut best_parameter) = (-1.0, 0.5);
    let mut test_uv = *uv;
    const GUESSES: usize = 10;
    for i in 0..GUESSES {
        let parameter = i as f64 / (GUESSES - 1) as f64;
        if singular == 0 {
            test_uv.x = sample.evaluate(parameter);
        } else {
            test_uv.y = sample.evaluate(parameter);
        }
        let m = evaluate(surface, test_uv, 1, 1, [true, true], true, true, true)?;
        let mut w = if singular == 0 { m[1] } else { m[2] };
        let _ = w.unitize();
        let value = w.dot(direction);
        if value > best_value {
            best_value = value;
            best_parameter = parameter;
        }
        if value > near_best_value - f_tolerance
            && (value > near_best_value + f_tolerance
                || (parameter - s).abs() < (near_best_parameter - s).abs())
        {
            near_best_value = value;
            near_best_parameter = parameter;
        }
    }
    if near_best_value < 0.0 && best_value < 0.0 {
        return Err(Fault::Degenerate);
    }
    let mut solver = LocalSolve1d::with_interval(
        SingularTangent {
            surface,
            uv: *uv,
            direction,
            singular,
        },
        s_interval,
        false,
    );
    let first = solver.solve_it(s_interval.evaluate(near_best_parameter), tolerance);
    let mut found_parameter = near_best_parameter;
    let mut root = first.unwrap_or(Root1d {
        found: false,
        t: 0.0,
    });
    if !root.found {
        root = solver
            .solve_it(s_interval.evaluate(best_parameter), tolerance)
            .unwrap_or(Root1d {
                found: false,
                t: root.t,
            });
        found_parameter = best_parameter;
    }
    let value = if root.found {
        root.t
    } else {
        s_interval.evaluate(found_parameter)
    };
    if singular == 0 {
        uv.x = value;
    } else {
        uv.y = value;
    }
    if first.is_ok() && root.found {
        Ok(())
    } else {
        Err(Fault::Degenerate)
    }
}

/// `SmSurface::InvertPointNearSingularity` as the stepper calls it: `true`
/// with `uv` moved when the test point is at the pole or on its normal, or
/// once it is dropped on the iso-curve through the pole along its
/// direction; `false` when SMLib's search fails and Newton steps continue.
fn invert_point_near_singularity(
    problem: &Problem,
    domain: Extent2,
    uv: &mut Vector2,
) -> Solve<bool> {
    let surface = problem.surface;
    let (test, tolerance) = (problem.test, problem.tight);
    let [point, su, sv] = first_derivatives(surface, *uv, [false, false])?;
    if point.closer_than(tolerance, test) {
        return Ok(true);
    }
    let (singular, singular_uv) = if su.length_squared() < sv.length_squared() {
        let v = domain.v_interval();
        (
            0,
            Vector2::new(uv.x, if uv.y < v.mid() { v.min } else { v.max }),
        )
    } else {
        let u = domain.u_interval();
        (
            1,
            Vector2::new(if uv.x < u.mid() { u.min } else { u.max }, uv.y),
        )
    };
    let (good, drop) = is_point_on_normal(surface, domain, singular_uv, test, tolerance);
    if good && (drop < tolerance || drop <= EFF_ZERO) {
        *uv = singular_uv;
        return Ok(true);
    }
    let pole = evaluate_point(surface, singular_uv).unwrap_or_default();
    let difference = test - pole;
    let mut moved = singular_uv;
    if find_degenerate_parameter(surface, &mut moved, difference, true, singular, tolerance)
        .is_err()
    {
        return Ok(false);
    }
    let (degenerate, good_domain) = if singular == 0 {
        (moved.x, domain.v_interval())
    } else {
        (moved.y, domain.u_interval())
    };
    let Ok(iso_curve) = surface.iso_curve(singular, degenerate) else {
        return Ok(false);
    };
    // SMLib lets the iso-curve evaluate outside its natural interval; the
    // port's curve clamps, which agrees only inside it.
    if !good_domain.is_contained_by(iso_curve.natural_interval(), 0.0) {
        return Err(SolverRefusal::NotPorted(problem.operation));
    }
    let [_, su, sv] = first_derivatives(surface, moved, [false, false])?;
    let good = if singular == 0 { sv } else { su };
    let guess = good_domain.clamp_value(difference.dot(good) / good.dot(good));
    let dropped = match drop_point(
        &iso_curve,
        good_domain,
        test,
        None,
        tolerance,
        Some(guess),
        SolverOperation::Minimize,
    ) {
        Ok(dropped) => dropped,
        Err(SolverRefusal::Failed(_)) => return Ok(false),
        Err(refusal) => return Err(refusal),
    };
    *uv = if singular == 0 {
        Vector2::new(degenerate, dropped.parameter)
    } else {
        Vector2::new(dropped.parameter, degenerate)
    };
    Ok(true)
}

/// `SmTol::MapTo2d`: a 3D tolerance as a parameter step along `direction`.
/// `SmSurface::IsInDomain`: `uv` lies in `domain` within the scaled zero,
/// or in the natural domain within `tolerance` mapped to each parameter
/// direction (SMLib tests the natural domain there, not `domain`).
pub(crate) fn is_in_domain(
    surface: &BSplineSurface,
    domain: Extent2,
    uv: Vector2,
    tolerance: f64,
) -> bool {
    if domain.contains_point(uv, EFF_ZERO) {
        return true;
    }
    let tolerance_u = map_to_2d(surface, tolerance, uv, Vector2::new(1.0, 0.0));
    let tolerance_v = map_to_2d(surface, tolerance, uv, Vector2::new(0.0, 1.0));
    let natural = surface.natural_domain();
    natural.u_interval().contains_value(uv.x, tolerance_u)
        && natural.v_interval().contains_value(uv.y, tolerance_v)
}

pub(crate) fn map_to_2d(
    surface: &BSplineSurface,
    tolerance: f64,
    uv: Vector2,
    direction: Vector2,
) -> f64 {
    let [_, du, dv] =
        first_derivatives(surface, uv, [true, true]).unwrap_or([Vector3::default(); 3]);
    let domain = surface.natural_domain();
    let (du_size, dv_size) = (du.length(), dv.length());
    let scaled_u = if du_size < EFF_ZERO {
        domain.u_interval().length() / 100.0
    } else {
        tolerance / du_size
    };
    let scaled_v = if dv_size < EFF_ZERO {
        domain.v_interval().length() / 100.0
    } else {
        tolerance / dv_size
    };
    let size = sqrt(direction.x * direction.x + direction.y * direction.y);
    if size < EFF_ZERO {
        min(scaled_u, scaled_v)
    } else {
        direction.x.abs() / size * scaled_u + direction.y.abs() / size * scaled_v
    }
}

/// `SmSurface::IsClosed` without the continuity output: the two sides of
/// `axis` meet at six samples, and a side that collapses to a point counts
/// only when it is a genuine pole. Without a `tolerance` each sample uses
/// a thousand scaled zeros of its far side point.
pub fn is_closed(
    surface: &BSplineSurface,
    domain: Extent2,
    axis: usize,
    tolerance: Option<f64>,
) -> bool {
    // `SmPlane::IsClosed`: planes are never closed.
    if surface.plane_frame().is_some() {
        return false;
    }
    const SAMPLES: usize = 5;
    let direction = if axis == 0 {
        Vector2::new(1.0, 0.0)
    } else {
        Vector2::new(0.0, 1.0)
    };
    let mut middle = [Vector3::default(); 4];
    let mut frame = [Vector3::default(); 4];
    for j in 0..=SAMPLES {
        let t = (j as f64 + 1.0) / (2.0 + SAMPLES as f64);
        let uv = domain.evaluate(t, t);
        let (low, high, mid) = if axis == 0 {
            let m = (domain.min.x + domain.max.x) / 2.0;
            (
                Vector2::new(domain.min.x, uv.y),
                Vector2::new(domain.max.x, uv.y),
                Vector2::new(m, uv.y),
            )
        } else {
            let m = (domain.min.y + domain.max.y) / 2.0;
            (
                Vector2::new(uv.x, domain.min.y),
                Vector2::new(uv.x, domain.max.y),
                Vector2::new(uv.x, m),
            )
        };
        let p_low = evaluate_point(surface, low).unwrap_or_default();
        let p_high = evaluate_point(surface, high).unwrap_or_default();
        let tolerance = tolerance.unwrap_or(EFF_ZERO * 1000.0 * (1.0 + p_high.max_dimension()));
        if (p_low - p_high).length_squared() > tolerance * tolerance {
            return false;
        }
        // EvaluateDirectionalDerivs to second order along the axis.
        if let Ok(d) = evaluate(surface, mid, 2, 2, [true, true], true, false, true) {
            let (x, y) = (direction.x, direction.y);
            middle[0] = d[0];
            middle[1] = x * d[3] + y * d[1];
            middle[2] = x * x * d[6] + 2.0 * x * y * d[4] + y * y * d[2];
        }
        let gap1 = (middle[0] - p_low).length_squared();
        let gap2 = (middle[0] - p_high).length_squared();
        if max(gap1, gap2) < 4.0 * tolerance * tolerance {
            let _ = frenet_from_derivatives(&middle, 2, &mut frame);
            let curvature2 = frame[2].length_squared();
            let radius = if curvature2 > EFF_ZERO_SQ {
                1.0 / sqrt(curvature2)
            } else {
                BIG_DOUBLE
            };
            if radius > tolerance {
                return false;
            }
        }
    }
    true
}

/// What `SmSurface::IsOnBoundary` finds.
pub(crate) struct OnBoundary {
    /// Near a side of the domain and inside the natural domain.
    pub(crate) on: bool,
    /// Near a u side and near a v side.
    pub(crate) sides: [bool; 2],
    /// The opposite u and v side values, where near.
    pub(crate) other: Vector2,
    /// The zone tolerance mapped to u and v steps.
    pub(crate) tolerances: [f64; 2],
}

/// `SmSurface::IsOnBoundary`: whether `uv` is within the 3D `zone`
/// tolerance, mapped to parameters at `uv`, of a side of `domain`, and
/// within it of the natural domain.
pub(crate) fn is_on_boundary(
    surface: &BSplineSurface,
    domain: Extent2,
    uv: Vector2,
    zone: f64,
) -> OnBoundary {
    let natural = surface.natural_domain();
    let tolerance_u = map_to_2d(surface, zone, uv, Vector2::new(1.0, 0.0));
    let tolerance_v = map_to_2d(surface, zone, uv, Vector2::new(0.0, 1.0));
    let (on_u, other_u) = domain.u_interval().is_value_on_boundary(uv.x, tolerance_u);
    let (on_v, other_v) = domain.v_interval().is_value_on_boundary(uv.y, tolerance_v);
    let on = (on_u || on_v)
        && natural.u_interval().contains_value(uv.x, tolerance_u)
        && natural.v_interval().contains_value(uv.y, tolerance_v);
    OnBoundary {
        on,
        sides: [on_u, on_v],
        other: Vector2::new(other_u, other_v),
        tolerances: [tolerance_u, tolerance_v],
    }
}

/// `SmSurface::IsClosed` with its continuity output: like [`is_closed`],
/// evaluating second derivatives at the paired sides, and also the lowest
/// continuity class across them.
pub(crate) fn is_closed_with_continuity(
    surface: &BSplineSurface,
    domain: Extent2,
    axis: usize,
) -> (bool, crate::geometry::Continuity) {
    use crate::geometry::{Continuity, CurveJet, evaluate_curve_continuity};
    if surface.plane_frame().is_some() {
        return (false, Continuity::Discontinuous);
    }
    const SAMPLES: usize = 5;
    let direction = if axis == 0 {
        Vector2::new(1.0, 0.0)
    } else {
        Vector2::new(0.0, 1.0)
    };
    let zero = Vector3::new(0.0, 0.0, 0.0);
    let mut minimum = Continuity::CInfinity;
    let mut middle = [Vector3::default(); 4];
    let mut frame = [Vector3::default(); 4];
    let (mut p_low, mut p_high) = (Vector3::default(), Vector3::default());
    let (mut low_d, mut high_d) = ([zero; 2], [zero; 2]);
    for j in 0..=SAMPLES {
        let t = (j as f64 + 1.0) / (2.0 + SAMPLES as f64);
        let uv = domain.evaluate(t, t);
        let (low, high, mid) = if axis == 0 {
            let m = (domain.min.x + domain.max.x) / 2.0;
            (
                Vector2::new(domain.min.x, uv.y),
                Vector2::new(domain.max.x, uv.y),
                Vector2::new(m, uv.y),
            )
        } else {
            let m = (domain.min.y + domain.max.y) / 2.0;
            (
                Vector2::new(uv.x, domain.min.y),
                Vector2::new(uv.x, domain.max.y),
                Vector2::new(uv.x, m),
            )
        };
        let second = |at: Vector2| evaluate(surface, at, 2, 2, [true, true], true, true, true);
        if let Ok(d) = second(low) {
            p_low = d[0];
            low_d = if axis == 0 {
                [d[3], d[6]]
            } else {
                [d[1], d[2]]
            };
        }
        if let Ok(d) = second(high) {
            p_high = d[0];
            high_d = if axis == 0 {
                [d[3], d[6]]
            } else {
                [d[1], d[2]]
            };
        }
        let tolerance = EFF_ZERO * 1000.0 * (1.0 + p_high.max_dimension());
        if (p_low - p_high).length_squared() > tolerance * tolerance {
            return (false, Continuity::Discontinuous);
        }
        let jet = |position: Vector3, d: [Vector3; 2]| CurveJet {
            position,
            first: d[0],
            second: d[1],
            third: None,
        };
        let class = evaluate_curve_continuity(
            &jet(p_low, low_d),
            &jet(p_high, high_d),
            crate::scalar::ANGLE_TOL_DEG,
        );
        if minimum > class {
            minimum = class;
        }
        if minimum > Continuity::Discontinuous {
            if let Ok(d) = evaluate(surface, mid, 2, 2, [true, true], true, false, true) {
                let (x, y) = (direction.x, direction.y);
                middle[0] = d[0];
                middle[1] = x * d[3] + y * d[1];
                middle[2] = x * x * d[6] + 2.0 * x * y * d[4] + y * y * d[2];
            }
            let gap1 = (middle[0] - p_low).length_squared();
            let gap2 = (middle[0] - p_high).length_squared();
            if max(gap1, gap2) < 4.0 * tolerance * tolerance {
                let _ = frenet_from_derivatives(&middle, 2, &mut frame);
                let curvature2 = frame[2].length_squared();
                let radius = if curvature2 > EFF_ZERO_SQ {
                    1.0 / sqrt(curvature2)
                } else {
                    BIG_DOUBLE
                };
                if radius > tolerance {
                    return (false, Continuity::Discontinuous);
                }
            }
        }
    }
    (true, minimum)
}

/// `SmSurface::IsOnSeam` over the natural domain: the closed directions of
/// the sides `uv` lies on, and the matching points across each seam.
pub fn is_on_seam(
    surface: &BSplineSurface,
    uv: Vector2,
    zone: f64,
) -> (Option<usize>, bool, Vec<Vector2>) {
    // `SmPlane::IsOnSeam`: a plane has no seam.
    if surface.plane_frame().is_some() {
        return (None, false, Vec::new());
    }
    let domain = surface.natural_domain();
    let boundary = is_on_boundary(surface, domain, uv, zone);
    if !boundary.on {
        return (None, false, Vec::new());
    }
    let [on_u, on_v] = boundary.sides;
    let [tolerance_u, tolerance_v] = boundary.tolerances;
    let (other_u, other_v) = (boundary.other.x, boundary.other.y);
    let closed_u = on_u && is_closed(surface, domain, 0, Some(tolerance_u));
    let closed_v = on_v && is_closed(surface, domain, 1, Some(tolerance_v));
    let direction = match (closed_u, closed_v) {
        (true, true) => Some(2),
        (true, false) => Some(0),
        (false, true) => Some(1),
        (false, false) => None,
    };
    let mut across = Vec::new();
    if closed_u {
        across.push(Vector2::new(other_u, uv.y));
    }
    if closed_v {
        across.push(Vector2::new(uv.x, other_v));
    }
    if on_u && on_v && closed_u && closed_v {
        across.push(Vector2::new(other_u, other_v));
    }
    (direction, closed_u || closed_v, across)
}

/// The answer of [`local_point_solve`] (`SmSolution`, single value, with
/// `LocalPointSolve_1`'s termination reason and iteration count).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfacePoint {
    /// The final parameters.
    pub uv: Vector2,
    /// The distance to the test point at `uv`.
    pub value: f64,
    /// Whether SMLib reports a solution.
    pub found: bool,
    /// Why iteration stopped.
    pub reason: TerminationReason,
    /// Completed iterations.
    pub iterations: usize,
}

/// `SmSurface::LocalPointSolve` (through `LocalPointSolve_1`): refines
/// `guess` to the nearest point (minimise), the foot of a normal
/// (normalise) or a coincident point (intersect) of `surface` within
/// `domain`. `owner_tolerance` is the zone tolerance of the face that owns
/// the surface, if any. Minimising still answers when it stops against the
/// domain boundary.
pub fn local_point_solve(
    surface: &BSplineSurface,
    domain: Extent2,
    operation: SolverOperation,
    test: Vector3,
    guess: Vector2,
    owner_tolerance: Option<f64>,
) -> Result<SurfacePoint, SolverRefusal> {
    if !matches!(
        operation,
        SolverOperation::Minimize | SolverOperation::Normalize | SolverOperation::Intersect
    ) {
        return Err(SolverRefusal::NotPorted(operation));
    }
    let guess = domain.clamp_point(guess);
    let loose = owner_tolerance.unwrap_or(if operation == SolverOperation::Intersect {
        0.01
    } else {
        EFF_ZERO_SQRT * 100.0
    });
    let zone = owner_tolerance.unwrap_or(ZONE_TOL_3D);
    let problem = Problem {
        surface,
        domain,
        test,
        operation,
        tight: EFF_ZERO * (1.0 + test.max_dimension()),
        loose,
    };
    let mut current = Stepper::new(&problem, guess);
    let start_value = current.gap_length()?;
    if current.converged(EFF_ZERO / 1000.0)?
        && (current.singular()?.is_none() || current.gap_length()? < problem.tight)
    {
        return Ok(SurfacePoint {
            uv: guess,
            value: start_value,
            found: true,
            reason: TerminationReason::FoundAnswerConverged,
            iterations: 0,
        });
    }
    const MAX_ITERATIONS: usize = 35;
    let mut reason = TerminationReason::StillIterating;
    let mut next = current.clone();
    let mut iterations = 0;
    while iterations < MAX_ITERATIONS {
        let ok = next.find_next_step(&mut current, &mut reason)?;
        if next.converged_tight(true)? {
            reason = TerminationReason::FoundAnswerConverged;
            break;
        } else if let Some(bound) = next.converged_on_boundary(problem.tight)? {
            next.snap_to_boundaries(problem.loose)?;
            if next.gap_length()? > problem.loose {
                let mut new_point = next.uv;
                let (seam, on_seam, across) = is_on_seam(surface, next.uv, zone);
                let which = match bound {
                    Bound::UMin | Bound::UMax => 0,
                    Bound::VMin | Bound::VMax => 1,
                    _ => 2,
                };
                if on_seam {
                    let seam = seam.unwrap_or(2);
                    new_point = if seam != 2 && (seam == which || which == 2) {
                        across[0]
                    } else if seam != 2 {
                        next.uv
                    } else if which != 2 {
                        across[which]
                    } else {
                        across[2]
                    };
                }
                if new_point.differs(next.uv) {
                    next.set_uv(new_point);
                } else {
                    reason = TerminationReason::OutOfBounds;
                    break;
                }
            } else {
                reason = TerminationReason::OutOfBounds;
                break;
            }
        }
        if !ok || reason != TerminationReason::StillIterating {
            break;
        }
        current = next.clone();
        current.iterations = iterations + 1;
        iterations += 1;
    }
    let uv = next.uv;
    let value = next.error_norm(Norm::Gap)?;
    let mut found = next.converged_tight(false)?;
    if !found {
        if operation == SolverOperation::Minimize {
            found = true;
        }
        if iterations >= MAX_ITERATIONS && reason != TerminationReason::OutOfBounds {
            reason = TerminationReason::UnableToConverge;
        }
    }
    Ok(SurfacePoint {
        uv,
        value,
        found,
        reason,
        iterations,
    })
}
