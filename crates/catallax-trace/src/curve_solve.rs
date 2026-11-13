// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline candidates of `spline.rs`;
// statuses become `Fault`, output references become returned values and
// solution containers become small result structs. Operations and fallbacks
// that need code from later port steps refuse with
// `SolverRefusal::NotPorted` instead of running.
// Translated from SMLib SmCurve.cpp (LocalPointSolve, SmFindPCExtremaEFO,
// LocalCurveIntersect, SmCCIntersectENFO, SmCCTangentIntENFO,
// LocalCurveSolve, SmFindCCNormalENFO, SimpleCoincidenceChecker,
// CurveMaxDistanceBetween, IsClosed, IsPeriodic); rust-brep has no port.
//! Local curve solvers of step 4: the nearest curve point to a point, and
//! the local meeting of two curves with its coincidence check.

use crate::Fault;
use crate::curve::{evaluate_geometric, unitize_tangent};
use crate::extent::Extent1;
use crate::global_solver::{SolverOperation, SolverRefusal};
use crate::local_solve::{
    BoundaryHandler, EvalFunction, EvalNFunctions, ExtentN, LocalSolve1d, LocalSolveNd,
    TerminationReason,
};
use crate::matrix::Matrix;
use crate::pseudo_box::PseudoBox;
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, are_same, max};
use crate::spline::{BSplineCurve, unique_knots};
use crate::vector::Vector3;

/// A point found on one curve (`SmSolution`, single value): the parameter
/// and the distance to the target.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvePoint {
    /// Curve parameter.
    pub parameter: f64,
    /// Distance from the target point.
    pub distance: f64,
}

/// A parameter pair on two curves with its solution value.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvePair {
    /// Parameters on the first and second curve.
    pub parameters: [f64; 2],
    /// The solution value: the distance between the points, or the
    /// operation's own measure.
    pub value: f64,
}

/// `IsClosed`: the ends meet and the middle is away from them.
pub fn is_closed(curve: &dyn SolverCurve, interval: Extent1, tolerance: f64) -> bool {
    let p1 = curve.evaluate_point(interval.min).unwrap_or_default();
    let p2 = curve.evaluate_point(interval.max).unwrap_or_default();
    let mid = curve
        .evaluate_point(interval.evaluate(0.5))
        .unwrap_or_default();
    let tol = if tolerance != 0.0 {
        tolerance
    } else {
        EFF_ZERO * (1.0 + p1.max_dimension())
    };
    p1.distance(p2) < tol && mid.distance(p1) > tol
}

/// `IsPeriodic`: closed with end tangents within one degree.
pub fn is_periodic(curve: &dyn SolverCurve, interval: Extent1) -> bool {
    if !is_closed(curve, interval, 0.0) {
        return false;
    }
    let start = curve.evaluate(interval.min, 1, true, true);
    let end = curve.evaluate(interval.max, 1, true, true);
    let (Ok(start), Ok(end)) = (start, end) else {
        return false;
    };
    match start[1].angle_between(end[1]) {
        Ok(angle) => angle < 1.0 * (crate::scalar::PI / 180.0),
        Err(_) => false,
    }
}

/// The `SmCurve` virtuals that SMLib's local point solver calls, so that it
/// runs on B-spline curves and on surface iso-curves alike.
pub trait PointSolveCurve {
    /// `Evaluate`: position and the first `derivatives` derivatives.
    fn evaluate(
        &self,
        t: f64,
        derivatives: usize,
        from_left: bool,
        nonzero_tangent: bool,
    ) -> Result<Vec<Vector3>, Fault>;
    /// `EvaluatePoint`.
    fn evaluate_point(&self, t: f64) -> Result<Vector3, Fault>;
    /// `IsDegenerate` over the natural interval.
    fn is_degenerate(&self, tolerance: f64) -> bool;
}

impl PointSolveCurve for BSplineCurve {
    fn evaluate(
        &self,
        t: f64,
        derivatives: usize,
        from_left: bool,
        nonzero_tangent: bool,
    ) -> Result<Vec<Vector3>, Fault> {
        BSplineCurve::evaluate(self, t, derivatives, from_left, nonzero_tangent)
    }
    fn evaluate_point(&self, t: f64) -> Result<Vector3, Fault> {
        BSplineCurve::evaluate_point(self, t)
    }
    fn is_degenerate(&self, tolerance: f64) -> bool {
        BSplineCurve::is_degenerate(self, tolerance)
    }
}

/// The further `SmCurve` virtuals the curve-curve solvers call, so that
/// they run on B-spline curves and surface iso-curves alike.
pub trait SolverCurve: PointSolveCurve {
    /// `GetNaturalInterval`.
    fn natural_interval(&self) -> Extent1;
    /// `GetDegree` (SMLib's default for curves without one is 3).
    fn degree(&self) -> usize;
    /// `GetNumberNaturalKnots` (SMLib's default is 2).
    fn natural_knot_count(&self) -> usize;
    /// `GetKnots`: the distinct knot values.
    fn knot_values(&self) -> Vec<f64>;
}

impl SolverCurve for BSplineCurve {
    fn natural_interval(&self) -> Extent1 {
        BSplineCurve::natural_interval(self)
    }
    fn degree(&self) -> usize {
        BSplineCurve::degree(self)
    }
    fn natural_knot_count(&self) -> usize {
        self.knots().len()
    }
    fn knot_values(&self) -> Vec<f64> {
        unique_knots(self.knots(), None)
            .iter()
            .map(|k| k.0)
            .collect()
    }
}

/// `SmFindPCExtremaEFO`: `F(t) = (C(t) - P) · C'(t)`, zero where the gap
/// is perpendicular to the curve. Large tangents are scaled down by powers
/// of ten, fixed at the first evaluation.
struct PointCurveExtrema<'a, C: PointSolveCurve + ?Sized> {
    curve: &'a C,
    test: Vector3,
    domain: Extent1,
    squared_tolerance: f64,
    maximize: bool,
    have_scale: bool,
    scale_steps: u32,
}

impl<C: PointSolveCurve + ?Sized> EvalFunction for PointCurveExtrema<'_, C> {
    fn evaluate(
        &mut self,
        t: f64,
        value: &mut f64,
        derivative: &mut f64,
        found: &mut bool,
        _signal_errors: bool,
    ) -> Result<(), Fault> {
        *value = 0.0;
        *derivative = 0.0;
        *found = false;
        let side = t < self.domain.mid();
        let pv = self.curve.evaluate(t, 2, side, true)?;
        let mut v = pv[0] - self.test;
        let mut pv = [pv[0], pv[1], pv[2]];
        if self.have_scale {
            let mut scale = 1.0;
            for _ in 0..self.scale_steps {
                scale *= 10.0;
            }
            if scale > 1.0 {
                v = v.divided(scale);
                for p in &mut pv {
                    *p = p.divided(scale);
                }
            }
        } else {
            while pv[1].length_squared() > 1.0e4 {
                self.scale_steps += 1;
                v = v.divided(10.0);
                for p in &mut pv {
                    *p = p.divided(10.0);
                }
            }
            self.have_scale = true;
        }
        let tan_dot_tan = pv[1].dot(pv[1]);
        if tan_dot_tan < EFF_ZERO_SQ {
            if self.curve.is_degenerate(EFF_ZERO) {
                *value = 0.0;
                *found = true;
                return Ok(());
            }
            return Err(Fault::Degenerate);
        }
        let v_dot_v = v.dot(v);
        let v_dot_tan = v.dot(pv[1]);
        let mut v_dot_second = v.dot(pv[2]);
        *value = v_dot_tan;
        let limit = 0.95;
        if v_dot_second < -tan_dot_tan * limit {
            v_dot_second = -tan_dot_tan * limit;
        }
        *derivative = tan_dot_tan + v_dot_second;
        if self.maximize {
            *derivative = -*derivative;
        }
        if v_dot_v < self.squared_tolerance {
            if self.maximize {
                return Err(Fault::Degenerate);
            }
            *found = true;
            return Ok(());
        }
        if (v_dot_tan * v_dot_tan) / (v_dot_v * tan_dot_tan) < self.squared_tolerance {
            *found = true;
        }
        Ok(())
    }
}

/// `LocalPointSolve` for minimising, maximising, normalising and
/// intersecting: Newton iteration on [`PointCurveExtrema`] from `guess`.
/// Hitting the interval bounds still answers, except when normalising.
pub fn local_point_solve<C: PointSolveCurve + ?Sized>(
    curve: &C,
    interval: Extent1,
    operation: SolverOperation,
    test: Vector3,
    guess: f64,
) -> Result<Option<CurvePoint>, SolverRefusal> {
    if !matches!(
        operation,
        SolverOperation::Minimize
            | SolverOperation::Maximize
            | SolverOperation::Normalize
            | SolverOperation::Intersect
    ) {
        return Err(SolverRefusal::NotPorted(operation));
    }
    let guess = interval.clamp_value(guess);
    let mut solver = LocalSolve1d::with_interval(
        PointCurveExtrema {
            curve,
            test,
            domain: interval,
            squared_tolerance: EFF_ZERO_SQ,
            maximize: operation == SolverOperation::Maximize,
            have_scale: false,
            scale_steps: 0,
        },
        interval,
        false,
    );
    let Ok(root) = solver.solve_it(guess, 1.0e-8) else {
        return Ok(None);
    };
    if !root.found
        && (solver.termination_reason() != TerminationReason::OutOfBounds
            || operation == SolverOperation::Normalize)
    {
        return Ok(None);
    }
    let point = curve.evaluate_point(root.t).unwrap_or_default();
    Ok(Some(CurvePoint {
        parameter: root.t,
        distance: point.distance(test),
    }))
}

fn evaluate2(curve: &dyn SolverCurve, t: f64) -> Result<[Vector3; 3], Fault> {
    let v = curve.evaluate(t, 2, true, true)?;
    Ok([v[0], v[1], v[2]])
}

/// `SmCCTangentIntENFO::StaticEvaluate`: parallel unit tangents and a zero
/// squared gap, for a tangent meeting of two curves.
fn tangent_evaluate(
    a: &dyn SolverCurve,
    b: &dyn SolverCurve,
    x: &[f64],
    f: &mut [f64],
    jacobian: Option<&mut Matrix>,
) -> Result<(), Fault> {
    let p1 = evaluate2(a, x[0])?;
    let p2 = evaluate2(b, x[1])?;
    if p2[1].length_squared() < EFF_ZERO || p1[1].length_squared() < EFF_ZERO {
        return Err(Fault::Degenerate);
    }
    let mut ufp = p1[1];
    ufp.unitize()?;
    let mut ugp = p2[1];
    ugp.unitize()?;
    let ufpp = p1[1].unitized_derivative(p1[2]);
    let ugpp = p2[1].unitized_derivative(p2[2]);
    let gap = p1[0] - p2[0];
    if let Some(j) = jacobian {
        j[0][0] = ufpp.dot(ugp);
        j[0][1] = ufp.dot(ugpp);
        j[1][0] = 2.0 * p1[1].dot(gap);
        j[1][1] = -2.0 * p2[1].dot(gap);
    }
    let sign = if ufp.dot(ugp) > 0.0 { 1.0 } else { -1.0 };
    f[0] = ufp.dot(ugp) - sign;
    f[1] = gap.dot(gap);
    Ok(())
}

/// `SmCCTangentIntENFO`.
struct TangentIntersection<'a> {
    a: &'a dyn SolverCurve,
    b: &'a dyn SolverCurve,
}

impl EvalNFunctions for TangentIntersection<'_> {
    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        jacobian: Option<&mut Matrix>,
        found: &mut bool,
    ) -> Result<(), Fault> {
        *found = false;
        tangent_evaluate(self.a, self.b, x, f, jacobian)
    }
}

/// `SmCCIntersectENFO`: the gap dotted with each unit tangent, switching to
/// the tangent functions near a tangent meeting.
struct CurveIntersection<'a> {
    a: &'a dyn SolverCurve,
    b: &'a dyn SolverCurve,
    solution_size: f64,
    scaled_zero: f64,
}

impl EvalNFunctions for CurveIntersection<'_> {
    fn solution_size(&self) -> f64 {
        self.solution_size
    }

    fn scaled_zero(&self) -> f64 {
        self.scaled_zero
    }

    fn has_intervals(&mut self) -> bool {
        true
    }

    fn set_intervals(
        &mut self,
        _intervals: Option<&ExtentN>,
        _periodic: Option<&[bool]>,
    ) -> Result<(), Fault> {
        Ok(())
    }

    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        jacobian: Option<&mut Matrix>,
        found: &mut bool,
    ) -> Result<(), Fault> {
        *found = false;
        let p1 = evaluate2(self.a, x[0])?;
        let p2 = evaluate2(self.b, x[1])?;
        let gap = p1[0] - p2[0];
        self.solution_size = gap.length();
        self.scaled_zero = EFF_ZERO * (1.0 + max(p1[0].max_dimension(), p2[0].max_dimension()));
        let tangent = gap.length_squared() < 0.01 && p1[1].is_parallel_to(p2[1], 0.001).0;
        if tangent {
            return tangent_evaluate(self.a, self.b, x, f, jacobian);
        }
        let (unit1, _) = unitize_tangent(&p1, -gap);
        let (unit2, _) = unitize_tangent(&p2, gap);
        f[0] = gap.dot(unit1);
        f[1] = gap.dot(unit2);
        if let Some(j) = jacobian {
            let unit1_t = p1[1].unitized_derivative(p1[2]);
            let unit2_t = p2[1].unitized_derivative(p2[2]);
            let term11 = p1[1].dot(unit1);
            let mut term12 = gap.dot(unit1_t);
            let term21 = -p2[1].dot(unit2);
            let mut term22 = gap.dot(unit2_t);
            if term12 < -term11 / 2.0 {
                term12 = -term11 / 2.0;
            }
            if term22 > -term21 / 2.0 {
                term22 = -term21 / 2.0;
            }
            j[0][0] = term11 + term12;
            j[0][1] = -p2[1].dot(unit1);
            j[1][0] = p1[1].dot(unit2);
            j[1][1] = term21 + term22;
        }
        Ok(())
    }
}

/// `SmFindCCNormalENFO` for minimising, maximising and normalising: the
/// gap dotted with each tangent. At the first evaluation, a gap below the
/// tolerance is stepped back along a curvature normal so the Newton step
/// can see which way to go.
struct CurveNormal<'a> {
    a: &'a dyn SolverCurve,
    b: &'a dyn SolverCurve,
    squared_tolerance: f64,
    count: u32,
    step_back: bool,
}

impl EvalNFunctions for CurveNormal<'_> {
    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        jacobian: Option<&mut Matrix>,
        found: &mut bool,
    ) -> Result<(), Fault> {
        if jacobian.is_some() {
            self.count += 1;
        }
        *found = false;
        let p1 = evaluate2(self.a, x[0])?;
        let p2 = evaluate2(self.b, x[1])?;
        let mut gap = p1[0] - p2[0];
        if p1[1].length_squared() < EFF_ZERO_SQ || p2[1].length_squared() < EFF_ZERO_SQ {
            return Err(Fault::Degenerate);
        }
        if (self.count == 1 || self.step_back) && gap.length_squared() < self.squared_tolerance {
            self.step_back = true;
            let mut shift = Vector3::new(0.0, 0.0, 0.0);
            let frame = evaluate_geometric(self.a, x[0], 2, true)?;
            if frame[2].length_squared() > EFF_ZERO_SQ {
                let mut normal = frame[2];
                let _ = normal.unitize();
                shift = normal;
            } else {
                let frame = evaluate_geometric(self.b, x[1], 2, true)?;
                if frame[2].length_squared() > EFF_ZERO_SQ {
                    let mut normal = frame[2];
                    let _ = normal.unitize();
                    shift = normal;
                }
            }
            gap = if shift.dot(gap) >= 0.0 {
                p1[0] + shift - p2[0]
            } else {
                p1[0] - shift - p2[0]
            };
        }
        if let Some(j) = jacobian {
            j[0][0] = p1[2].dot(gap) + p1[1].dot(p1[1]);
            j[0][1] = -p1[1].dot(p2[1]);
            j[1][0] = -j[0][1];
            j[1][1] = p2[2].dot(gap) - p2[1].dot(p2[1]);
        }
        f[0] = p1[1].dot(gap);
        f[1] = p2[1].dot(gap);
        if gap.length_squared() < self.squared_tolerance {
            *found = true;
            return Ok(());
        }
        let projected1 = (f[0] * f[0]) / gap.length_squared();
        let projected2 = (f[1] * f[1]) / gap.length_squared();
        if projected1 < self.squared_tolerance && projected2 < self.squared_tolerance {
            *found = true;
        }
        Ok(())
    }
}

fn intervals(a: Extent1, b: Extent1) -> ExtentN {
    ExtentN {
        intervals: vec![a, b],
    }
}

/// `LocalCurveSolve` for minimising, maximising and normalising: Newton
/// iteration on [`CurveNormal`] from the guesses; the value is the distance
/// between the points reached. Hitting the bounds still answers, except when
/// normalising.
#[allow(clippy::too_many_arguments)]
pub fn local_curve_solve(
    a: &dyn SolverCurve,
    interval_a: Extent1,
    b: &dyn SolverCurve,
    interval_b: Extent1,
    operation: SolverOperation,
    tolerance: f64,
    guesses: [f64; 2],
) -> Result<Option<CurvePair>, SolverRefusal> {
    let (found, pair) =
        local_curve_solve_with_last(a, interval_a, b, interval_b, operation, tolerance, guesses)?;
    Ok(found.then_some(pair))
}

/// `LocalCurveSolve` also fills its last iterate when its found flag is
/// false. Analytic circle intersections check that pair for improvement.
#[allow(clippy::too_many_arguments)]
pub(crate) fn local_curve_solve_with_last(
    a: &dyn SolverCurve,
    interval_a: Extent1,
    b: &dyn SolverCurve,
    interval_b: Extent1,
    operation: SolverOperation,
    tolerance: f64,
    guesses: [f64; 2],
) -> Result<(bool, CurvePair), SolverRefusal> {
    if !matches!(
        operation,
        SolverOperation::Minimize | SolverOperation::Maximize | SolverOperation::Normalize
    ) {
        return Err(SolverRefusal::NotPorted(operation));
    }
    let mut solver = LocalSolveNd::new(
        CurveNormal {
            a,
            b,
            squared_tolerance: EFF_ZERO * EFF_ZERO,
            count: 0,
            step_back: false,
        },
        Some(intervals(interval_a, interval_b)),
        Some(vec![false, false]),
    );
    solver.set_boundary_handler(BoundaryHandler::TotalBoundaryHits, 8);
    let root = solver.solve_it(&guesses, tolerance)?;
    let p1 = a.evaluate_point(root.x[0])?;
    let p2 = b.evaluate_point(root.x[1])?;
    let found = root.found
        || (solver.termination_reason() == TerminationReason::OutOfBounds
            && operation != SolverOperation::Normalize);
    Ok((
        found,
        CurvePair {
            parameters: [root.x[0], root.x[1]],
            value: p1.distance(p2),
        },
    ))
}

fn intersect_solver<'a>(
    a: &'a dyn SolverCurve,
    b: &'a dyn SolverCurve,
    interval_a: Extent1,
    interval_b: Extent1,
) -> LocalSolveNd<CurveIntersection<'a>> {
    let mut solver = LocalSolveNd::new(
        CurveIntersection {
            a,
            b,
            solution_size: 0.0,
            scaled_zero: 1.0,
        },
        Some(intervals(interval_a, interval_b)),
        Some(vec![false, false]),
    );
    solver.set_desired_accuracy(EFF_ZERO / 10.0);
    solver.set_boundary_handler(BoundaryHandler::HitsAllBoundaries, 3);
    solver
}

fn first_derivatives(
    a: &dyn SolverCurve,
    ta: f64,
    b: &dyn SolverCurve,
    tb: f64,
) -> Result<([Vector3; 2], [Vector3; 2]), Fault> {
    let p1 = a.evaluate(ta, 1, true, true)?;
    let p2 = b.evaluate(tb, 1, true, true)?;
    Ok(([p1[0], p1[1]], [p2[0], p2[1]]))
}

/// `LocalCurveIntersect`: Newton iteration for a meeting of two curves near
/// the guesses, a second start from the first answer, a closest-approach
/// solve and a tangent refinement. Returns the parameters and their gap when
/// the gap is below `tolerance`.
pub fn local_curve_intersect(
    a: &dyn SolverCurve,
    interval_a: Extent1,
    b: &dyn SolverCurve,
    interval_b: Extent1,
    tolerance: f64,
    guesses: [f64; 2],
) -> Result<Option<CurvePair>, SolverRefusal> {
    let mut solver = intersect_solver(a, b, interval_a, interval_b);
    let root = solver.solve_it(&guesses, tolerance * 2.0)?;
    if !root.found {
        return Ok(None);
    }
    let mut solution = [root.x[0], root.x[1]];
    let (mut pv1, mut pv2) = first_derivatives(a, solution[0], b, solution[1])?;
    let mut deviation = pv1[0].distance(pv2[0]);
    if deviation < EFF_ZERO_SQ {
        return Ok(Some(CurvePair {
            parameters: solution,
            value: deviation,
        }));
    }
    // Negated as in SMLib, so a NaN gap is not an intersection.
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    if !(deviation < tolerance) {
        return Ok(None);
    }
    if deviation > EFF_ZERO {
        let old = solution;
        let same = |p: f64, q: f64| (p - q).abs() <= EFF_ZERO_SQRT;
        if !same(old[0], guesses[0]) || !same(old[1], guesses[1]) {
            let mut again = intersect_solver(a, b, interval_a, interval_b);
            let root = again.solve_it(&solution, tolerance * 2.0)?;
            if root.found {
                (pv1, pv2) = first_derivatives(a, root.x[0], b, root.x[1])?;
                solution = [root.x[0], root.x[1]];
                let new_deviation = pv1[0].distance(pv2[0]);
                if new_deviation < deviation {
                    deviation = new_deviation;
                } else {
                    solution = old;
                }
            } else {
                solution = old;
            }
        }
    }
    if deviation > EFF_ZERO
        && !are_same(interval_a.min, interval_a.max)
        && !are_same(interval_b.min, interval_b.max)
        && let Some(closest) = local_curve_solve(
            a,
            interval_a,
            b,
            interval_b,
            SolverOperation::Minimize,
            tolerance,
            solution,
        )?
        && closest.value < deviation
    {
        deviation = closest.value;
        solution = closest.parameters;
    }
    if deviation > EFF_ZERO && pv1[1].is_parallel_to(pv2[1], 2.0).0 {
        let mut tangent = LocalSolveNd::new(
            TangentIntersection { a, b },
            Some(intervals(interval_a, interval_b)),
            Some(vec![is_periodic(a, interval_a), is_periodic(b, interval_b)]),
        );
        tangent.set_desired_accuracy(EFF_ZERO_SQ);
        let root = tangent.solve_it(&solution, tolerance * 2.0)?;
        let (q1, q2) = first_derivatives(a, root.x[0], b, root.x[1])?;
        let new_deviation = q1[0].distance(q2[0]);
        if new_deviation < deviation {
            return Ok(Some(CurvePair {
                parameters: [root.x[0], root.x[1]],
                value: new_deviation,
            }));
        }
    }
    Ok(Some(CurvePair {
        parameters: solution,
        value: deviation,
    }))
}

/// What [`curve_max_distance_between`] found.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaxDistance {
    /// The largest gap seen.
    pub distance: f64,
    /// Parameter on this curve where it was seen.
    pub this_parameter: f64,
    /// Parameter on the other curve where it was seen.
    pub other_parameter: f64,
    /// Whether a nearest point was found out of order along this curve;
    /// SMLib does not report it when the search stops early at `needed`.
    pub incorrect_direction: Option<bool>,
}

/// `CurveMaxDistanceBetween`: samples the other curve between two of its
/// parameters, drops each sample onto this curve and keeps the largest gap,
/// stopping once it exceeds `needed`. With no sample count it samples
/// densely and refines the largest gap. A sample that does not drop locally
/// needs SMLib's global point search, which a later port step brings: that
/// refuses.
#[allow(clippy::too_many_arguments)]
pub fn curve_max_distance_between(
    this: &dyn SolverCurve,
    interval: Extent1,
    other: &dyn SolverCurve,
    other_start: f64,
    other_end: f64,
    samples: usize,
    needed: Option<f64>,
) -> Result<MaxDistance, SolverRefusal> {
    if std::ptr::addr_eq(this, other)
        && ((are_same(interval.min, other_start) && are_same(interval.max, other_end))
            || (are_same(interval.min, other_end) && are_same(interval.max, other_start)))
    {
        return Ok(MaxDistance {
            distance: 0.0,
            this_parameter: interval.mid(),
            other_parameter: interval.mid(),
            incorrect_direction: None,
        });
    }
    let mut other_interval = Extent1::new(other_start, other_start);
    other_interval.add_value(other_end);
    let same_direction = other_end >= other_start;
    let mut samples = samples;
    let mut precise = false;
    if samples == 0 {
        precise = true;
        samples = 22 + this.knot_values().len() + other.knot_values().len();
    }
    if this.degree() == 1
        && this.natural_knot_count() == 4
        && other.degree() == 1
        && other.natural_knot_count() == 4
    {
        precise = false;
        samples = 2;
    }
    let mut this_max = interval.evaluate(0.5);
    let mut other_max = other_interval.evaluate(0.5);
    let mut max_distance = 0.0;
    let mut incorrect = false;
    let mut previous: Option<f64> = None;
    for ii in 0..samples {
        let fraction = ii as f64 / (samples - 1) as f64;
        let guess1 = interval.evaluate(fraction);
        let mut guess2 = if same_direction {
            other_interval.evaluate(fraction)
        } else {
            other_interval.evaluate(1.0 - fraction)
        };
        let mut point = other.evaluate_point(guess2)?;
        let mut found =
            local_point_solve(this, interval, SolverOperation::Minimize, point, guess1)?;
        if found.is_none() && ii == 0 {
            guess2 = other_interval.evaluate(0.001);
            point = other.evaluate_point(guess2)?;
            found = local_point_solve(this, interval, SolverOperation::Minimize, point, guess1)?;
        }
        if found.is_none() && ii == samples - 1 {
            guess2 = other_interval.evaluate(0.999);
            point = other.evaluate_point(guess2)?;
            found = local_point_solve(this, interval, SolverOperation::Minimize, point, guess1)?;
        }
        let Some(hit) = found else {
            return Err(SolverRefusal::NotPorted(SolverOperation::Minimize));
        };
        let distance = this.evaluate_point(hit.parameter)?.distance(point);
        if distance > max_distance {
            this_max = hit.parameter;
            other_max = guess2;
            max_distance = distance;
        }
        if let Some(previous) = previous
            && (if same_direction {
                previous > hit.parameter
            } else {
                previous < hit.parameter
            })
        {
            incorrect = true;
        }
        previous = Some(hit.parameter);
        if let Some(needed) = needed
            && needed < max_distance
        {
            return Ok(MaxDistance {
                distance: max_distance,
                this_parameter: this_max,
                other_parameter: other_max,
                incorrect_direction: None,
            });
        }
    }
    let point = this.evaluate_point(this_max)?;
    if precise
        && max_distance > EFF_ZERO * (1.0 + point.max_dimension())
        && let Some(far) = local_curve_solve(
            this,
            interval,
            other,
            other_interval,
            SolverOperation::Maximize,
            EFF_ZERO_SQRT,
            [this_max, other_max],
        )?
    {
        let [t1, t2] = far.parameters;
        let distance = this.evaluate_point(t1)?.distance(other.evaluate_point(t2)?);
        if distance > max_distance {
            this_max = t1;
            other_max = t2;
            max_distance = distance;
        }
    }
    Ok(MaxDistance {
        distance: max_distance,
        this_parameter: this_max,
        other_parameter: other_max,
        incorrect_direction: Some(incorrect),
    })
}

/// What [`simple_coincidence_check`] found between two curve pieces.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Coincidence {
    /// No end of either piece lies on the other.
    None,
    /// One point of contact.
    Point(CurvePair),
    /// Two separate points of contact.
    Points([CurvePair; 2]),
    /// The pieces coincide between these two parameter pairs.
    Range([CurvePair; 2]),
}

/// `SimpleCoincidenceChecker`: drops each piece's ends onto the other piece,
/// then decides between one contact, two contacts and a coincident range by
/// the largest gap between the contacts. Assumes neither piece is closed and
/// both are flat enough to drop points locally.
#[allow(clippy::too_many_arguments)]
pub fn simple_coincidence_check(
    this: &dyn SolverCurve,
    interval: Extent1,
    other: &dyn SolverCurve,
    other_interval: Extent1,
    tolerance: f64,
    this_box: Option<&PseudoBox>,
    other_box: Option<&PseudoBox>,
) -> Result<Coincidence, SolverRefusal> {
    let start1 = this.evaluate_point(interval.min)?;
    let end1 = this.evaluate_point(interval.max)?;
    let start2 = other.evaluate_point(other_interval.min)?;
    let end2 = other.evaluate_point(other_interval.max)?;
    let near = |pseudo: Option<&PseudoBox>, p: Vector3| {
        pseudo.is_none_or(|b| b.distance_to_point(p) < tolerance)
    };
    let mut found: Vec<CurvePair> = Vec::with_capacity(4);
    let keep = |hit: Option<CurvePoint>| hit.filter(|h| h.distance < tolerance);
    if near(other_box, start1)
        && let Some(hit) = keep(local_point_solve(
            other,
            other_interval,
            SolverOperation::Minimize,
            start1,
            other_interval.evaluate(0.5),
        )?)
    {
        found.push(CurvePair {
            parameters: [interval.min, hit.parameter],
            value: hit.distance,
        });
    }
    if near(this_box, start2)
        && let Some(hit) = keep(local_point_solve(
            this,
            interval,
            SolverOperation::Minimize,
            start2,
            interval.evaluate(0.5),
        )?)
    {
        found.push(CurvePair {
            parameters: [interval.clamp_value(hit.parameter), other_interval.min],
            value: hit.distance,
        });
    }
    if near(this_box, end2)
        && let Some(hit) = keep(local_point_solve(
            this,
            interval,
            SolverOperation::Minimize,
            end2,
            interval.evaluate(0.5),
        )?)
    {
        let pair = CurvePair {
            parameters: [interval.clamp_value(hit.parameter), other_interval.max],
            value: hit.distance,
        };
        if let Some(last) = found.last()
            && pair.parameters[0] < last.parameters[0]
        {
            let index = found.len() - 1;
            found.insert(index, pair);
        } else {
            found.push(pair);
        }
    }
    if near(other_box, end1)
        && let Some(hit) = keep(local_point_solve(
            other,
            other_interval,
            SolverOperation::Minimize,
            end1,
            other_interval.evaluate(0.5),
        )?)
    {
        found.push(CurvePair {
            parameters: [interval.max, hit.parameter],
            value: hit.distance,
        });
    }
    let Some(&first) = found.first() else {
        return Ok(Coincidence::None);
    };
    let last = found[found.len() - 1];
    if found.len() == 1
        || (last.parameters[0] - first.parameters[0]).abs() < EFF_ZERO
        || (last.parameters[1] - first.parameters[1]).abs() < EFF_ZERO
    {
        let pick = if found.len() == 1 || first.value < last.value {
            first
        } else {
            last
        };
        return Ok(Coincidence::Point(pick));
    }
    let span1 = Extent1::new(first.parameters[0], last.parameters[0]);
    let mut span2 = Extent1::new(first.parameters[1], first.parameters[1]);
    span2.add_value(last.parameters[1]);
    let coincident = if are_same(span1.max, span1.min) || are_same(span2.max, span2.min) {
        false
    } else {
        let degree = this.degree().max(other.degree());
        let gap = curve_max_distance_between(
            this,
            span1,
            other,
            first.parameters[1],
            last.parameters[1],
            degree + 3,
            Some(tolerance),
        )?;
        gap.distance <= tolerance
    };
    Ok(if coincident {
        Coincidence::Range([first, last])
    } else {
        Coincidence::Points([first, last])
    })
}
