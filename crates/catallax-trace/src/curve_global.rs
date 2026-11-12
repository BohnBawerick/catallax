// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline candidates of `spline.rs`;
// statuses become `Fault`, the point tree becomes one node appended to the
// curve cache's arena, and solution containers keep SMLib's `SmSolution`
// layout. Operations other than minimising, maximising, normalising and
// intersecting refuse with `SolverRefusal::NotPorted`, and so do curves
// SMLib hands to an analytic circle (`SmCircle::DropPointFast`) or to the
// sampled line of a rational two-point curve. A polynomial two-point line
// goes to its analytic `SmLine` (`line.rs`), as in SMLib.
// Translated from SMLib SmCurve.cpp (SmCurve::GlobalPointSolve,
// SmCPGlobalSolver, SmCurve::DropPoint, sm_SnapSolToSpecifiedEndpoint) and
// SmBSplineCurve.cpp (CreateAnalyticCurve's line, and the analytic tests of
// its arc); rust-brep
// has no port.
//! The global point-curve solver of step 4: every nearest, farthest or
//! normal point of a curve to a test point, searched over the curve cache,
//! and the curve point drop built on it.

use crate::Fault;
use crate::curve::{CurveSide, evaluate_continuity};
use crate::curve_cache::{CurveCache, NodeData};
use crate::curve_solve::{CurvePoint, is_closed, local_point_solve};
use crate::extent::{Extent1, Extent3};
use crate::geometry::{Continuity, segment_point_angle_deg_to_perp};
use crate::global_solver::{
    AuxData, Branch, GlobalSolver, SearchTrees, SolutionRequested, SolverHooks, SolverOperation,
    SolverRefusal, TreeNode,
};
use crate::line::Line;
use crate::scalar::{
    ANGLE_TOL_DEG, BIG_DOUBLE, EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, PI, ZONE_TOL_3D, are_same,
    is_zero,
};
use crate::solution::{
    MAX_PARAMETERS, ObjectId, ObjectKind, Solution, SolutionArray, SolutionEnd, SolutionObjects,
    SolutionType, SortKey,
};
use crate::spline::BSplineCurve;
use crate::vector::Vector3;

/// The solution identity of the curve.
pub const CURVE_OBJECT: ObjectId = ObjectId(0);

/// `SM_EFF_ZERO_RAD`.
const EFF_ZERO_RAD: f64 = 1.0e-4 * PI / 180.0;

/// The curve a point-curve solution names.
pub struct CurveObject<'a> {
    /// The curve.
    pub curve: &'a BSplineCurve,
    /// Its zone tolerance (SMLib's global one without an owner).
    pub zone: f64,
}

impl SolutionObjects for CurveObject<'_> {
    fn kind(&self, id: ObjectId) -> ObjectKind {
        if id == CURVE_OBJECT {
            ObjectKind::Curve
        } else {
            ObjectKind::Other
        }
    }

    fn parameter_interval(&self, id: ObjectId, _direction: usize) -> Extent1 {
        if id == CURVE_OBJECT {
            self.curve.natural_interval()
        } else {
            Extent1::default()
        }
    }

    fn zone_tolerance(&self, id: ObjectId) -> f64 {
        if id == CURVE_OBJECT { self.zone } else { 0.0 }
    }
}

fn ported(operation: SolverOperation) -> Result<(), SolverRefusal> {
    if matches!(
        operation,
        SolverOperation::Minimize
            | SolverOperation::Maximize
            | SolverOperation::Normalize
            | SolverOperation::Intersect
    ) {
        Ok(())
    } else {
        Err(SolverRefusal::NotPorted(operation))
    }
}

/// `SM_IS_ZERO_TO_TOL`.
fn is_zero_to_tol(a: f64, tol: f64) -> bool {
    a.abs() <= tol
}

/// `SmCPGlobalSolver`: the point tree is tree 0, the curve tree tree 1.
struct CurvePointHooks<'a> {
    curve: &'a BSplineCurve,
    cache: &'a CurveCache,
    interval: Extent1,
    test: Vector3,
}

impl CurvePointHooks<'_> {
    /// Whether answers sort by value rather than by parameter.
    fn by_value(solver: &GlobalSolver<'_>) -> SortKey {
        let operation = solver.operation;
        if matches!(
            operation,
            SolverOperation::Minimize | SolverOperation::Maximize
        ) || (operation == SolverOperation::Normalize
            && solver.requested == SolutionRequested::Single)
        {
            SortKey::ByValueWithParameter
        } else {
            SortKey::ByFirstParameter
        }
    }

    fn solution(branch: &Branch, kind: SolutionType, start: SolutionEnd) -> Solution {
        let mut solution = Solution {
            kind,
            variables: 1,
            start,
            object_count: 1,
            ..Solution::default()
        };
        solution.objects[0] = Some(CURVE_OBJECT);
        solution.nodes[0] = branch[0];
        solution.nodes[1] = branch[1];
        solution
    }
}

fn end(t: f64, value: f64) -> SolutionEnd {
    let mut parameters = [0.0; MAX_PARAMETERS];
    parameters[0] = t;
    SolutionEnd { value, parameters }
}

impl SolverHooks for CurvePointHooks<'_> {
    fn branch_may_contain_answers(
        &mut self,
        solver: &GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<bool, SolverRefusal> {
        if solver.operation != SolverOperation::Normalize {
            return solver.base_branch_may_contain_answers(branch);
        }
        let Some(curve_id) = branch[1] else {
            return Err(Fault::Degenerate.into());
        };
        let node = self
            .cache
            .trees
            .get(curve_id)
            .ok_or(SolverRefusal::Failed(Fault::Degenerate))?;
        if node.aux != AuxData::BezierCurve && node.aux != AuxData::CurveData {
            return Ok(true);
        }
        let span = match self.cache.node_data[curve_id.0 as usize] {
            NodeData::Leaf(index) => self.cache.spans[index],
            NodeData::Parent(index) => self.cache.parents[index],
            NodeData::None => return Ok(true),
        };
        if span.max_turning_angle_deg > 90.0 {
            return Ok(true);
        }
        let start = self
            .curve
            .evaluate_point(span.interval.min)
            .unwrap_or_default();
        let finish = self
            .curve
            .evaluate_point(span.interval.max)
            .unwrap_or_default();
        if self.test.closer_than(solver.tolerance, start)
            || self.test.closer_than(solver.tolerance, finish)
        {
            return Ok(true);
        }
        let (angle, _) = segment_point_angle_deg_to_perp(start, finish, self.test);
        Ok(angle * 0.99 - EFF_ZERO <= span.max_turning_angle_deg)
    }

    #[allow(clippy::too_many_lines)]
    fn local_solve(
        &mut self,
        solver: &mut GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<bool, SolverRefusal> {
        let Some(curve_id) = branch[1] else {
            return Err(Fault::Degenerate.into());
        };
        let NodeData::Leaf(span_index) = self.cache.node_data[curve_id.0 as usize] else {
            return Err(Fault::Degenerate.into());
        };
        let span = self.cache.spans[span_index];
        let span_box = self.cache.pseudo_boxes[span_index];
        let (curve, test, interval) = (self.curve, self.test, self.interval);
        let operation = solver.operation;
        let tolerance = solver.tolerance;
        let single = solver.requested == SolutionRequested::Single;
        let mut bez = span.interval;
        let mut conts = span.continuities;

        if bez.are_disjoint(interval, 0.0) {
            return Ok(false);
        }
        if bez.min <= interval.min {
            conts[0] = Continuity::Discontinuous;
        }
        if bez.max >= interval.max {
            conts[1] = Continuity::Discontinuous;
        }
        if !bez.is_contained_by(interval, EFF_ZERO) {
            if let Ok(meet) = bez.intersect(interval) {
                bez = meet;
            }
            if are_same(bez.max, bez.min) {
                return Ok(false);
            }
        }

        let natural = curve.natural_interval();
        if (is_zero(bez.min - natural.min) || is_zero(bez.max - natural.max))
            && (is_zero(interval.min - natural.min) && is_zero(interval.max - natural.max))
            && is_closed(curve, natural, 0.0)
        {
            let class = evaluate_continuity(
                CurveSide {
                    curve,
                    parameter: natural.min,
                    from_left: false,
                },
                CurveSide {
                    curve,
                    parameter: natural.max,
                    from_left: true,
                },
                ANGLE_TOL_DEG,
            )
            .map_or(Continuity::Discontinuous, |(class, _, _)| class);
            if is_zero(bez.min - natural.min) {
                conts[0] = class;
            }
            if is_zero(bez.max - natural.max) {
                conts[1] = class;
            }
        }

        if operation == SolverOperation::Minimize
            || (operation == SolverOperation::Normalize && single)
        {
            let mut distance = span_box.distance_to_point(test) - tolerance;
            if distance < 0.0 {
                distance = 0.0;
            }
            if distance * distance > solver.best_answer_so_far_sq {
                return Ok(false);
            }
        }
        if operation == SolverOperation::Intersect && span_box.distance_to_point(test) > tolerance {
            return Ok(false);
        }

        let do_line_test = span.max_turning_angle_deg < 90.0
            && (operation == SolverOperation::Normalize
                || (matches!(
                    operation,
                    SolverOperation::Minimize | SolverOperation::Maximize
                ) && conts[0] > Continuity::C0
                    && conts[1] > Continuity::C0));

        let start = curve.evaluate(bez.min, 1, true, true)?;
        let finish = curve.evaluate(bez.max, 1, true, true)?;
        let (start_point, start_inward) = (start[0], -start[1]);
        let (end_point, end_tangent) = (finish[0], finish[1]);
        let chord = end_point - start_point;
        let key = Self::by_value(solver);

        if chord.length_squared() < EFF_ZERO_SQ {
            if operation == SolverOperation::Normalize && !single {
                return Ok(false);
            }
            let distance = end_point.distance(test);
            if operation == SolverOperation::Intersect && distance > tolerance {
                return Ok(false);
            }
            let solution =
                Self::solution(branch, SolutionType::SingleValue, end(bez.min, distance));
            solver.add_sorted_solution(&*self, &solution, key, true);
            return Ok(false);
        }

        let (angle_from_perp, dt) = segment_point_angle_deg_to_perp(start_point, end_point, test);
        if do_line_test
            && angle_from_perp * 0.99 - EFF_ZERO > span.max_turning_angle_deg
            && !test.closer_than(tolerance, start_point)
            && !test.closer_than(tolerance, end_point)
        {
            return Ok(false);
        }

        let mut guesses = vec![bez.evaluate(dt)];
        for at in [0.5, 0.0, 1.0] {
            if !is_zero_to_tol(dt - at, EFF_ZERO_SQRT) {
                guesses.push(bez.evaluate(at));
            }
        }
        let local_operation = if operation == SolverOperation::Intersect {
            guesses.truncate(1);
            SolverOperation::Minimize
        } else {
            operation
        };

        // SMLib keeps one `SmSolution` across the guesses; a guess that finds
        // nothing leaves the previous answer in it.
        let mut current = CurvePoint {
            parameter: 0.0,
            distance: 0.0,
        };
        let mut last_found = false;
        let mut coincident: Vec<CurvePoint> = Vec::new();
        let mut first_value = BIG_DOUBLE;
        let (mut last_deviation, mut last_parameter) = (0.0, 0.0);

        for (i, &guess) in guesses.iter().enumerate() {
            let found = match local_point_solve(curve, bez, local_operation, test, guess)? {
                Some(point) => {
                    current = point;
                    true
                }
                None => false,
            };
            if i == 0 && !found {
                if operation == SolverOperation::Normalize {
                    continue;
                }
                return Ok(false);
            }
            if last_found
                && found
                && is_zero(last_deviation - current.distance)
                && is_zero(last_parameter - current.parameter)
            {
                continue;
            }
            if found {
                last_found = true;
                last_deviation = current.distance;
                last_parameter = current.parameter;
            }
            if i == 0 {
                first_value = current.distance;
            }
            if operation == SolverOperation::Intersect && current.distance > tolerance {
                return Ok(false);
            }

            let mut boundary = None;
            if is_zero(current.parameter - bez.min) && conts[0] >= Continuity::G1 {
                boundary = Some((start_point, start_inward));
            }
            if is_zero(current.parameter - bez.max) && conts[1] >= Continuity::G1 {
                boundary = Some((end_point, end_tangent));
            }
            if let Some((point, direction)) = boundary
                && local_operation != SolverOperation::Normalize
                && !is_zero_to_tol(current.distance, EFF_ZERO * 100.0)
            {
                let gap = test - point;
                let mut angle = 0.0;
                let mut outside = false;
                if matches!(
                    local_operation,
                    SolverOperation::Minimize | SolverOperation::Maximize
                ) {
                    if gap.length_squared() > EFF_ZERO_SQ {
                        angle = direction.angle_between(gap).unwrap_or(0.0);
                    }
                    angle = (angle - PI / 2.0).abs();
                    outside = direction.dot(gap) > 0.0;
                }
                if !is_zero_to_tol(angle, EFF_ZERO_RAD) && outside {
                    continue;
                }
            }

            if is_zero(current.parameter - guess) && is_zero(current.distance - first_value) {
                coincident.push(current);
                continue;
            }
            if !found {
                continue;
            }
            let solution = Self::solution(
                branch,
                SolutionType::SingleValue,
                end(current.parameter, current.distance),
            );
            solver.add_sorted_solution(&*self, &solution, key, true);
        }

        if coincident.len() < 3 {
            for point in &coincident {
                let solution = Self::solution(
                    branch,
                    SolutionType::SingleValue,
                    end(point.parameter, point.distance),
                );
                solver.add_sorted_solution(&*self, &solution, key, true);
            }
        } else {
            let mut range =
                Self::solution(branch, SolutionType::RangeOfValues, SolutionEnd::default());
            let (mut done, mut first, mut last) = (0, false, false);
            for point in &coincident {
                if is_zero(point.parameter - bez.min) {
                    done += 1;
                    first = true;
                    range.start = end(point.parameter, point.distance);
                }
                if is_zero(point.parameter - bez.max) {
                    done += 1;
                    last = true;
                    range.end = end(point.parameter, point.distance);
                }
            }
            if done == 2 && first && last {
                solver.add_sorted_solution(&*self, &range, key, true);
            }
        }
        Ok(false)
    }
}

/// Whether SMLib's curve cache would hand `curve` to an analytic curve
/// (`SmCurveCache::GetAnalyticCurve`, `SmBSplineCurve::CreateAnalyticCurve`):
/// a two-point line, or a rational quadratic that may be a circular arc.
/// The arc test is conservative: every rational quadratic that is not
/// degenerate counts, so a curve SMLib keeps as a spline may refuse but none
/// SMLib solves analytically runs the spline search.
pub(crate) fn may_be_analytic(curve: &BSplineCurve, interval: Extent1) -> Result<bool, Fault> {
    let middle = curve.evaluate_point(interval.evaluate(0.456_789))?;
    let scaled_zero = EFF_ZERO * (1.0 + middle.max_dimension());
    if curve.degree() == 1 && curve.control_points().len() == 2 {
        return curve.is_line(10, scaled_zero);
    }
    Ok(curve.degree() == 2 && curve.is_rational() && !curve.is_degenerate(scaled_zero))
}

/// What `SmBSplineCurve::CreateAnalyticCurve` makes of a polynomial
/// two-point line over `interval`: the `SmLine` through its control points
/// with `interval` as its analytic domain, or `None` when its ends do not
/// reproduce the curve's within the scaled zero. Other curves are `None`.
pub(crate) fn analytic_line(
    curve: &BSplineCurve,
    interval: Extent1,
) -> Result<Option<Line>, Fault> {
    if curve.degree() != 1 || curve.control_points().len() != 2 || curve.is_rational() {
        return Ok(None);
    }
    let middle = curve.evaluate_point(interval.evaluate(0.456_789))?;
    let scaled_zero = EFF_ZERO * (1.0 + middle.max_dimension());
    if !curve.is_line(10, scaled_zero)? {
        return Ok(None);
    }
    let [(point, _), (end, _)] = curve.ends();
    let vector = end - point;
    let origin = point - interval.min * (vector / interval.length());
    let direction = vector / vector.length();
    let scale = vector.length() / interval.length();
    let line = Line::new(
        origin,
        direction,
        interval,
        scale,
        curve.dimension() == 2,
        false,
    );
    for t in [interval.min, interval.max] {
        if curve.evaluate_point(t)?.distance(line.evaluate_point(t)) > scaled_zero {
            return Ok(None);
        }
    }
    Ok(Some(line))
}

/// The analytic branch of `SmCurve::GlobalPointSolve` for a line: the
/// line's drop, refined by a local solve on the curve from its parameter
/// when it is not on the curve and the refinement comes closer.
fn analytic_line_solve(
    curve: &BSplineCurve,
    line: &Line,
    interval: Extent1,
    operation: SolverOperation,
    test: Vector3,
    tolerance: f64,
) -> Result<SolutionArray, SolverRefusal> {
    let mut solutions = SolutionArray::new();
    let meeting = operation == SolverOperation::Intersect;
    let Some((parameter, distance)) = line.drop_point_fast(interval, meeting, test, tolerance)
    else {
        return Ok(solutions);
    };
    let mut answer = (parameter, distance);
    if distance > EFF_ZERO
        && let Some(found) = local_point_solve(curve, interval, operation, test, parameter)?
        && found.distance < distance
    {
        answer = (found.parameter, found.distance);
    }
    let mut solution = Solution {
        kind: SolutionType::SingleValue,
        variables: 1,
        start: end(answer.0, answer.1),
        object_count: 1,
        ..Solution::default()
    };
    solution.objects[0] = Some(CURVE_OBJECT);
    solutions.add(solution);
    Ok(solutions)
}

/// `SmCurve::GlobalPointSolve`: the points of `curve` over `interval` that
/// minimise, maximise or are normal to the distance from `test`, or meet it
/// within `tolerance`, as single values sorted by value or parameter.
///
/// SMLib's search status is ignored, so a refusal inside the search returns
/// the answers found before it, as SMLib does; an unported operation or an
/// analytic curve refuses up front.
pub fn global_point_solve(
    curve: &BSplineCurve,
    interval: Extent1,
    operation: SolverOperation,
    test: Vector3,
    tolerance: f64,
    requested: SolutionRequested,
) -> Result<SolutionArray, SolverRefusal> {
    ported(operation)?;
    let cache = CurveCache::for_curve(curve, 0.5 * ZONE_TOL_3D)?;
    if matches!(
        operation,
        SolverOperation::Minimize | SolverOperation::Intersect
    ) {
        if let Some(line) = analytic_line(curve, cache.interval)? {
            return analytic_line_solve(curve, &line, interval, operation, test, tolerance);
        }
        // Circular arcs (`SmCircle::DropPointFast`) and rational two-point
        // lines (`SmCurve::IsLine`'s sampled line) are not ported.
        if may_be_analytic(curve, cache.interval)? && !(curve.degree() == 1 && !curve.is_rational())
        {
            return Err(SolverRefusal::NotPorted(operation));
        }
    }
    let mut trees: SearchTrees = cache.trees.clone();
    let point = trees.add(TreeNode {
        bbox: Extent3::singleton(test),
        children: None,
        aux: AuxData::None,
    });
    let objects = CurveObject {
        curve,
        zone: ZONE_TOL_3D,
    };
    let mut solver =
        GlobalSolver::with_two_trees(&trees, &objects, Some(point), Some(cache.top), 1);
    let mut hooks = CurvePointHooks {
        curve,
        cache: &cache,
        interval,
        test,
    };
    let best = match operation {
        SolverOperation::Maximize => 0.0,
        _ => BIG_DOUBLE,
    };
    let mut solutions = SolutionArray::new();
    match solver.solve_it(
        &mut hooks,
        operation,
        requested,
        tolerance,
        best,
        &[],
        &mut solutions,
    ) {
        Err(SolverRefusal::NotPorted(what)) => Err(SolverRefusal::NotPorted(what)),
        _ => Ok(solutions),
    }
}

/// `sm_SnapSolToSpecifiedEndpoint`: with `inward`, an answer within
/// `snap_tolerance` of an interval end is kept and snapped there when the
/// curve leaves that end along `inward`, and dropped otherwise.
fn snap_to_specified_endpoint(
    curve: &BSplineCurve,
    interval: Extent1,
    inward: Option<Vector3>,
    parameter: &mut f64,
    snap_tolerance: f64,
) -> Result<bool, Fault> {
    let Some(inward) = inward else {
        return Ok(true);
    };
    let mut good = true;
    let point = curve.evaluate_point(*parameter)?;
    let low = curve.evaluate(interval.min, 1, true, true)?;
    let mut at_min = false;
    if low[0].distance(point) < snap_tolerance {
        at_min = true;
        if inward.dot(low[1]) < 0.0 {
            good = false;
        } else {
            good = true;
            *parameter = interval.min;
        }
    }
    let high = curve.evaluate(interval.max, 1, true, true)?;
    if high[0].distance(point) < snap_tolerance && (!at_min || !good) {
        if inward.dot(-high[1]) < 0.0 {
            good = false;
        } else {
            good = true;
            *parameter = interval.max;
        }
    }
    Ok(good)
}

/// What [`drop_point`] found: SMLib's `rbSuccess`, `rdDroppedParameter` and
/// `rdDistanceToCurve`, which keep the local answer even when it is not a
/// success.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DroppedPoint {
    /// The drop is an answer.
    pub success: bool,
    /// Curve parameter.
    pub parameter: f64,
    /// Distance from the dropped point.
    pub distance: f64,
}

/// `SmCurve::DropPoint`: the nearest (or, normalising, a normal) point of
/// `curve` over `interval` to `point`, from `guess` with the local solver and
/// otherwise with the global one; `inward` names an end to snap to.
pub fn drop_point(
    curve: &BSplineCurve,
    interval: Extent1,
    point: Vector3,
    inward: Option<Vector3>,
    tolerance: f64,
    guess: Option<f64>,
    operation: SolverOperation,
) -> Result<DroppedPoint, SolverRefusal> {
    let operation = match operation {
        SolverOperation::Minimize | SolverOperation::Intersect | SolverOperation::Normalize => {
            operation
        }
        _ => SolverOperation::Minimize,
    };
    let mut dropped = DroppedPoint {
        success: false,
        parameter: 0.0,
        distance: 0.0,
    };
    let intersect_tolerance = if operation == SolverOperation::Intersect {
        tolerance
    } else {
        BIG_DOUBLE
    };
    let mut converged = false;
    let mut local_value = 0.0;
    let local = match guess {
        Some(guess) => local_point_solve(curve, interval, operation, point, guess)?,
        None => None,
    };
    if let Some(found) = local {
        converged = true;
        local_value = found.distance;
        dropped.parameter = found.parameter;
        dropped.distance = found.distance;
        if found.distance < intersect_tolerance {
            let mut parameter = found.parameter;
            let keep =
                snap_to_specified_endpoint(curve, interval, inward, &mut parameter, tolerance)?;
            if keep && found.distance < intersect_tolerance {
                dropped = DroppedPoint {
                    success: true,
                    parameter,
                    distance: found.distance,
                };
            }
        }
    }
    if !converged || (operation == SolverOperation::Intersect && local_value > intersect_tolerance)
    {
        let solutions = global_point_solve(
            curve,
            interval,
            operation,
            point,
            tolerance,
            SolutionRequested::Single,
        )?;
        for solution in solutions.as_slice() {
            let mut parameter = solution.start.parameters[0];
            let keep =
                snap_to_specified_endpoint(curve, interval, inward, &mut parameter, tolerance)?;
            if keep && solution.start.value < intersect_tolerance {
                dropped = DroppedPoint {
                    success: true,
                    parameter,
                    distance: solution.start.value,
                };
                break;
            }
        }
    }
    Ok(dropped)
}
