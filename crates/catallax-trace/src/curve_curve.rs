// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline curves of `spline.rs`;
// statuses become `Fault`, tree-node pointers become identities in one arena
// holding both curve trees, and solution containers keep SMLib's
// `SmSolution` layout. Paths no caller in the port reaches refuse with
// `SolverRefusal::NotPorted`: a curve intersected with itself, two curves
// that SMLib's caches hand to analytic lines or arcs
// (`AnalyticCurveIntersect`), and a zero tangent where SMLib reads past its
// two-vector array (`sm_UnitizeTangent`).
// Translated from SMLib SmCurve.cpp (SmCurve::GlobalCurveIntersect,
// SmCCIGlobalSolver::LocalSolve, sm_FilterSolution, sm_UnitizeTangent,
// sm_FindAdditionalCcTangentGuess, sm_RefineCCSolutions,
// sm_IsPointSolutionAcrossClosedSeamBoundary,
// SmCurve::GlobalCoincidenceChecker, SmCurve::ClassifyIntersectingRangeSolution,
// SmCurve::ApproximatePolarBox) and SmGeomUtility.cpp
// (smgu_LineLineClosestPoint, smgu_LinePlaneIntersect); rust-brep has no port.
//! The global curve-curve intersector: both curve trees searched together,
//! local intersections in each pair of spans, the coincidence check of the
//! whole curves, and the combination of the answers into points and
//! coincident ranges.
// Indexed loops follow SMLib's loops over span ends.
#![allow(clippy::needless_range_loop)]

use crate::Fault;
use crate::curve::evaluate_geometric;
use crate::curve_cache::{CurveCache, NodeData};
use crate::curve_global::{drop_point, may_be_analytic};
use crate::curve_solve::{
    Coincidence, CurvePair, curve_max_distance_between, is_closed, local_curve_intersect,
    local_curve_solve, local_point_solve, simple_coincidence_check,
};
use crate::curve_surface_global::degenerate_on_to;
use crate::extent::{Extent1, Extent3};
use crate::geometry::{
    Continuity, circle_coplanar_circle_intersect, line_closest_point,
    line_coplanar_circle_intersect,
};
use crate::global_solver::{
    Branch, GlobalSolver, SearchTrees, SolutionRequested, SolverHooks, SolverOperation,
    SolverRefusal,
};
use crate::polar_box::{PolarBox, PolarState};
use crate::scalar::{
    BIG_DOUBLE, EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, PI, ZONE_TOL_3D, are_same, is_zero, max, min,
    scaled_zero,
};
use crate::solution::{
    MAX_PARAMETERS, NodeId, ObjectId, ObjectKind, Solution, SolutionArray, SolutionEnd,
    SolutionObjects, SolutionType, SortKey,
};
use crate::spline::BSplineCurve;
use crate::vector::Vector3;

/// The first curve's object identity in solutions.
pub const FIRST_CURVE: ObjectId = ObjectId(0);
/// The second curve's object identity in solutions.
pub const SECOND_CURVE: ObjectId = ObjectId(1);

const INTERSECT: SolverOperation = SolverOperation::Intersect;

/// The two curves a curve-curve solution names.
pub struct CurvePairObjects<'a> {
    /// The first curve.
    pub first: &'a BSplineCurve,
    /// The second curve.
    pub second: &'a BSplineCurve,
}

impl SolutionObjects for CurvePairObjects<'_> {
    fn kind(&self, id: ObjectId) -> ObjectKind {
        match id {
            FIRST_CURVE | SECOND_CURVE => ObjectKind::Curve,
            _ => ObjectKind::Other,
        }
    }

    fn parameter_interval(&self, id: ObjectId, _direction: usize) -> Extent1 {
        match id {
            FIRST_CURVE => self.first.natural_interval(),
            SECOND_CURVE => self.second.natural_interval(),
            _ => Extent1::default(),
        }
    }

    fn zone_tolerance(&self, id: ObjectId) -> f64 {
        match id {
            FIRST_CURVE | SECOND_CURVE => ZONE_TOL_3D,
            _ => 0.0,
        }
    }
}

fn end(t1: f64, t2: f64, value: f64) -> SolutionEnd {
    let mut parameters = [0.0; MAX_PARAMETERS];
    parameters[0] = t1;
    parameters[1] = t2;
    SolutionEnd { value, parameters }
}

/// A single-value curve-curve solution naming both curves.
fn point_solution(branch: Option<&Branch>, t1: f64, t2: f64, value: f64) -> Solution {
    let mut solution = Solution {
        kind: SolutionType::SingleValue,
        variables: 2,
        start: end(t1, t2, value),
        object_count: 2,
        ..Solution::default()
    };
    solution.objects[0] = Some(FIRST_CURVE);
    solution.objects[1] = Some(SECOND_CURVE);
    if let Some(branch) = branch {
        solution.nodes[0] = branch[0];
        solution.nodes[1] = branch[1];
    }
    solution
}

/// `smgu_LinePlaneIntersect` (point and normal form).
fn line_plane_intersect(
    line_point: Vector3,
    line_vector: Vector3,
    plane_point: Vector3,
    normal: Vector3,
) -> Option<f64> {
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    if !(line_vector.length_squared() >= EFF_ZERO_SQ) || !(normal.length_squared() >= EFF_ZERO_SQ) {
        return None;
    }
    let d = -normal.dot(plane_point);
    let tolerance = EFF_ZERO * line_vector.max_dimension() + EFF_ZERO_SQ;
    let denominator = line_vector.dot(normal);
    if denominator.abs() < tolerance {
        return None;
    }
    Some(-(d + line_point.dot(normal)) / denominator)
}

/// `smgu_LineLineClosestPoint`: the parameters of the closest points of two
/// lines, `None` for a zero or parallel line.
pub(crate) fn line_line_closest_point(
    p1: Vector3,
    v1: Vector3,
    p2: Vector3,
    v2: Vector3,
) -> Option<(f64, f64)> {
    if v1.length_squared() < EFF_ZERO_SQ || v2.length_squared() < EFF_ZERO_SQ {
        return None;
    }
    let mut cross = v1.cross(v2);
    if cross.length_squared() < EFF_ZERO_SQ {
        return None;
    }
    cross.unitize().ok()?;
    let mut plane1 = cross.cross(v1);
    let mut plane2 = cross.cross(v2);
    plane1.unitize().ok()?;
    plane2.unitize().ok()?;
    let t1 = line_plane_intersect(p1, v1, p2, plane2)?;
    let t2 = line_plane_intersect(p2, v2, p1, plane1)?;
    Some((t1, t2))
}

/// `sm_UnitizeTangent` for a nonzero tangent. SMLib falls back to a
/// second derivative its callers here never pass (it reads past their
/// two-vector arrays), so a tangent no longer than `SM_EFF_ZERO` refuses.
fn unitize_tangent(tangent: Vector3) -> Result<Vector3, SolverRefusal> {
    let length = tangent.length();
    if length > EFF_ZERO {
        Ok(tangent / length)
    } else {
        Err(SolverRefusal::NotPorted(INTERSECT))
    }
}

/// `sm_FilterSolution` (curve-curve): whether to keep an answer found on a
/// span. Answers at a span end that is not a curve end are kept only when
/// the curves do not continue past each other there, at a point of closest
/// approach, or at a C0 end.
#[allow(clippy::too_many_arguments)]
fn filter_solution(
    curve1: &BSplineCurve,
    interval1: Extent1,
    curve2: &BSplineCurve,
    interval2: Extent1,
    conts1: [Continuity; 2],
    conts2: [Continuity; 2],
    solution: &Solution,
) -> Result<bool, SolverRefusal> {
    if solution.start.value < EFF_ZERO {
        return Ok(true);
    }
    let (t1, t2) = (solution.start.parameters[0], solution.start.parameters[1]);
    let pv1 = curve1.evaluate(t1, 1, true, true)?;
    let pv2 = curve2.evaluate(t2, 1, true, true)?;
    let scaled = EFF_ZERO
        * (1.0
            + max(
                max(
                    max(pv1[0].max_dimension(), pv1[1].max_dimension()),
                    max(pv2[0].max_dimension(), pv2[1].max_dimension()),
                ),
                max(t1, t2),
            ));
    if solution.start.value < scaled {
        return Ok(true);
    }
    let start1 = are_same(interval1.min, t1);
    let end1 = are_same(interval1.max, t1);
    let mid1 = !start1 && !end1;
    let next1 =
        mid1 || (start1 && conts1[0] >= Continuity::G1) || (end1 && conts1[1] >= Continuity::G1);
    let start2 = are_same(interval2.min, t2);
    let end2 = are_same(interval2.max, t2);
    let mid2 = !start2 && !end2;
    let next2 =
        mid2 || (start2 && conts2[0] >= Continuity::G1) || (end2 && conts2[1] >= Continuity::G1);
    let gap = pv2[0] - pv1[0];
    if (!mid1 && next1) || (!mid2 && next2) {
        let mut tangent1 = unitize_tangent(pv1[1])?;
        if start1 {
            tangent1 = -tangent1;
        }
        let mut tangent2 = unitize_tangent(pv2[1])?;
        if start2 {
            tangent2 = -tangent2;
        }
        let gap1 = tangent1.dot(gap);
        let gap2 = tangent2.dot(gap);
        let limit = EFF_ZERO_SQRT;
        if (!mid1 && next1 && gap1 > limit) || (!mid2 && next2 && gap2 < -limit) {
            return Ok(false);
        }
        if !mid1 && next1 && (-tangent2 - tangent1).dot(gap) < -limit {
            return Ok(false);
        }
        if !mid1 && next1 && next2 && (tangent2 - tangent1).dot(gap) < -limit {
            return Ok(false);
        }
        if !mid2 && next2 && (tangent2 + tangent1).dot(gap) < -limit {
            return Ok(false);
        }
        if !mid2 && next2 && next1 && (tangent2 - tangent1).dot(gap) < -limit {
            return Ok(false);
        }
    }
    if gap.is_perpendicular_to(pv1[1], 1.0) && gap.is_perpendicular_to(pv2[1], 1.0) {
        return Ok(true);
    }
    Ok((conts1[0] <= Continuity::C0 && start1)
        || (conts1[1] <= Continuity::C0 && end1)
        || (conts2[0] <= Continuity::C0 && start2)
        || (conts2[1] <= Continuity::C0 && end2))
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

/// `sm_FindAdditionalCcTangentGuess`: near a tangency, where the curves'
/// osculating circles (or a line and a circle) cross a second time, a
/// guess for a second meeting.
fn additional_tangent_guess(
    curve1: &BSplineCurve,
    interval1: Extent1,
    curve2: &BSplineCurve,
    interval2: Extent1,
    t1: f64,
    t2: f64,
) -> Result<Option<(f64, f64)>, SolverRefusal> {
    let pv1 = curve1.evaluate(t1, 1, true, true)?;
    let pv2 = curve2.evaluate(t2, 1, true, true)?;
    if !pv1[1].is_parallel_to(pv2[1], 20.0).0 {
        return Ok(None);
    }
    let geometric1 = evaluate_geometric(curve1, t1, 2, true)?;
    let mut curvature_vector1 = geometric1[2];
    let curvature1 = curvature_vector1.length();
    if curvature1 > EFF_ZERO {
        curvature_vector1.unitize()?;
    }
    let geometric2 = evaluate_geometric(curve2, t2, 2, true)?;
    let mut curvature_vector2 = geometric2[2];
    let curvature2 = curvature_vector2.length();
    if curvature2 > EFF_ZERO {
        curvature_vector2.unitize()?;
    }
    let mut line_circle = false;
    let radius1;
    let mut radius2 = 0.0;
    let (mut line_point, mut line_vector) = (Vector3::default(), Vector3::default());
    let center1;
    let mut center2 = Vector3::default();
    let mut normal1;
    let scaled_tolerance = 1.0e-4 * (1.0 + curvature1.abs() + curvature2.abs());
    if curvature1.abs() < scaled_tolerance {
        if curvature2.abs() < scaled_tolerance {
            return Ok(None);
        }
        radius1 = (1.0 / curvature2).abs();
        center1 = pv2[0] + curvature_vector2 * radius1;
        normal1 = geometric2[1].cross(geometric2[2]);
        if normal1.length_squared() < EFF_ZERO_SQ {
            return Ok(None);
        }
        normal1.unitize()?;
        line_point = pv1[0];
        line_vector = pv1[1];
        line_vector.unitize()?;
        line_circle = true;
    } else {
        if curvature2.abs() < scaled_tolerance {
            line_point = pv2[0];
            line_vector = pv2[1];
            line_vector.unitize()?;
            line_circle = true;
        } else {
            radius2 = (1.0 / curvature2).abs();
            center2 = pv2[0] + curvature_vector2 * radius2;
            let mut normal2 = geometric2[1].cross(geometric2[2]);
            if normal2.length_squared() < EFF_ZERO_SQ {
                return Ok(None);
            }
            normal2.unitize()?;
        }
        radius1 = (1.0 / curvature1).abs();
        center1 = pv1[0] + curvature_vector1 * radius1;
        normal1 = geometric1[1].cross(geometric1[2]);
        if normal1.length_squared() < EFF_ZERO_SQ {
            return Ok(None);
        }
        normal1.unitize()?;
    }
    let found = if line_circle {
        let hits = line_coplanar_circle_intersect(
            line_point,
            line_vector,
            center1,
            normal1,
            radius1,
            EFF_ZERO,
        )?;
        if hits.count != 2 {
            return Ok(None);
        }
        let mut t = hits.parameters[0];
        if hits.parameters[1].abs() > hits.parameters[0].abs() {
            t = hits.parameters[1];
        }
        t * line_vector
    } else {
        let (count, points) = circle_coplanar_circle_intersect(
            normal1, center1, radius1, center2, radius2, EFF_ZERO,
        )?;
        if count != 2 {
            return Ok(None);
        }
        let mut found = points[0] - pv1[0];
        if found.length() < points[1].distance(pv1[0]) {
            found = points[1] - pv1[0];
        }
        found
    };
    let guess1 = t1 + parameter_delta(pv1[1], found);
    let guess2 = t2 + parameter_delta(pv2[1], found);
    if !interval1.contains_value(guess1, 0.0) && !interval2.contains_value(guess2, 0.0) {
        return Ok(None);
    }
    Ok(Some((
        interval1.clamp_value(guess1),
        interval2.clamp_value(guess2),
    )))
}

/// `SmCCIGlobalSolver`: the global solver over both curves' trees.
struct CurveCurveHooks<'a> {
    first: &'a BSplineCurve,
    second: &'a BSplineCurve,
    first_cache: &'a CurveCache,
    second_cache: &'a CurveCache,
    second_offset: u32,
    interval1: Extent1,
    interval2: Extent1,
}

/// One leaf span of a curve tree.
#[derive(Clone, Copy)]
struct Leaf {
    index: usize,
    bbox: Extent3,
}

impl CurveCurveHooks<'_> {
    fn leaf(&self, id: NodeId, second: bool) -> Result<Leaf, SolverRefusal> {
        let (cache, local) = if second {
            (self.second_cache, NodeId(id.0 - self.second_offset))
        } else {
            (self.first_cache, id)
        };
        let node = cache.trees.get(local).ok_or(Fault::Degenerate)?;
        let NodeData::Leaf(index) = cache.node_data[local.0 as usize] else {
            return Err(Fault::Degenerate.into());
        };
        Ok(Leaf {
            index,
            bbox: node.bbox,
        })
    }
}

/// The end points of a span within the solved intervals (`sLinePt`,
/// `sEnd`) and the chord between them.
struct Chord {
    start: Vector3,
    end: Vector3,
    vector: Vector3,
}

fn chord(curve: &BSplineCurve, interval: Extent1) -> Result<Chord, SolverRefusal> {
    let start = curve.evaluate_point(interval.min)?;
    let end = curve.evaluate_point(interval.max)?;
    Ok(Chord {
        start,
        end,
        vector: end - start,
    })
}

impl SolverHooks for CurveCurveHooks<'_> {
    #[allow(clippy::too_many_lines)]
    fn local_solve(
        &mut self,
        solver: &mut GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<bool, SolverRefusal> {
        let tolerance = solver.tolerance;
        let (Some(id1), Some(id2)) = (branch[0], branch[1]) else {
            return Err(Fault::Degenerate.into());
        };
        let leaf1 = self.leaf(id1, false)?;
        let leaf2 = self.leaf(id2, true)?;
        let span1 = self.first_cache.spans[leaf1.index];
        let span2 = self.second_cache.spans[leaf2.index];
        let box1 = self.first_cache.pseudo_boxes[leaf1.index];
        let box2 = self.second_cache.pseudo_boxes[leaf2.index];
        let (curve1, curve2) = (self.first, self.second);
        let mut interval1 = span1.interval;
        let mut interval2 = span2.interval;
        let mut expanded = box1;
        expanded.expand_absolute(tolerance);
        if expanded.are_disjoint(&box2) {
            return Ok(false);
        }
        if solver.operation == SolverOperation::IntersectionTest {
            let mut bbox = leaf1.bbox;
            bbox.expand_absolute(tolerance);
            let meet = bbox.intersect(leaf2.bbox).unwrap_or_default();
            if solver.solutions.is_empty() {
                let mut solution = Solution {
                    kind: SolutionType::SingleValue,
                    variables: 7,
                    object_count: 2,
                    ..Solution::default()
                };
                solution.nodes[0] = branch[0];
                solution.nodes[1] = branch[1];
                solution.start.parameters[6] = 1.0;
                solution.start.parameters[0] = meet.min.x;
                solution.start.parameters[1] = meet.min.y;
                solution.start.parameters[2] = meet.min.z;
                solution.start.parameters[3] = meet.max.x;
                solution.start.parameters[4] = meet.max.y;
                solution.start.parameters[5] = meet.max.z;
                solution.end.parameters[0] = interval1.min;
                solution.end.parameters[1] = interval1.max;
                solution.end.parameters[2] = interval2.min;
                solution.end.parameters[3] = interval2.max;
                solution.objects[0] = Some(FIRST_CURVE);
                solution.objects[1] = Some(SECOND_CURVE);
                if are_same(interval1.min, self.interval1.min)
                    && are_same(interval1.max, self.interval1.max)
                    && are_same(interval2.min, self.interval2.min)
                    && are_same(interval2.max, self.interval2.max)
                {
                    solution.start.parameters[6] = 2.0;
                }
                solver.solutions.add(solution);
            } else {
                let solution = &mut solver.solutions.as_mut_slice()[0];
                solution.start.parameters[6] += 1.0;
                let p = solution.start.parameters;
                let old = Extent3::new(
                    Vector3::new(p[0], p[1], p[2]),
                    Vector3::new(p[3], p[4], p[5]),
                );
                let union = meet.union(old);
                solution.start.parameters[0] = union.min.x;
                solution.start.parameters[1] = union.min.y;
                solution.start.parameters[2] = union.min.z;
                solution.start.parameters[3] = union.max.x;
                solution.start.parameters[4] = union.max.y;
                solution.start.parameters[5] = union.max.z;
                let e = solution.end.parameters;
                let old1 = Extent1 {
                    min: e[0],
                    max: e[1],
                };
                let old2 = Extent1 {
                    min: e[2],
                    max: e[3],
                };
                let union1 = Extent1 {
                    min: min(old1.min, interval1.min),
                    max: max(old1.max, interval1.max),
                };
                let union2 = Extent1 {
                    min: min(old2.min, interval2.min),
                    max: max(old2.max, interval2.max),
                };
                solution.end.parameters[0] = union1.min;
                solution.end.parameters[1] = union1.max;
                solution.end.parameters[2] = union2.min;
                solution.end.parameters[3] = union2.max;
            }
            return Ok(false);
        }
        if self.interval1.are_disjoint(interval1, 0.0)
            || self.interval2.are_disjoint(interval2, 0.0)
        {
            return Ok(false);
        }
        let mut conts1 = span1.continuities;
        let mut conts2 = span2.continuities;
        if interval1.min <= self.interval1.min {
            conts1[0] = Continuity::Discontinuous;
        }
        if interval1.max >= self.interval1.max {
            conts1[1] = Continuity::Discontinuous;
        }
        if !interval1.is_contained_by(self.interval1, EFF_ZERO) {
            interval1 = interval1.intersect(self.interval1).unwrap_or_default();
            if are_same(interval1.max, interval1.min) {
                return Ok(false);
            }
        }
        let line1 = chord(curve1, interval1)?;
        if interval2.min <= self.interval2.min {
            conts2[0] = Continuity::Discontinuous;
        }
        if interval2.max >= self.interval2.max {
            conts2[1] = Continuity::Discontinuous;
        }
        if !interval2.is_contained_by(self.interval2, EFF_ZERO) {
            interval2 = interval2.intersect(self.interval2).unwrap_or_default();
            if are_same(interval2.max, interval2.min) {
                return Ok(false);
            }
        }
        let line2 = chord(curve2, interval2)?;
        let natural1 = curve1.natural_interval();
        let natural2 = curve2.natural_interval();
        for (interval, whole, natural, conts, curve) in [
            (interval1, self.interval1, natural1, &mut conts1, curve1),
            (interval2, self.interval2, natural2, &mut conts2, curve2),
        ] {
            if (is_zero(interval.min - natural.min) || is_zero(interval.max - natural.max))
                && is_zero(whole.min - natural.min)
                && is_zero(whole.max - natural.max)
                && is_closed(curve, natural, 0.0)
            {
                let continuity = crate::curve::evaluate_continuity(
                    crate::curve::CurveSide {
                        curve,
                        parameter: natural.min,
                        from_left: false,
                    },
                    crate::curve::CurveSide {
                        curve,
                        parameter: natural.max,
                        from_left: true,
                    },
                    crate::scalar::ANGLE_TOL_DEG,
                )
                .map_or(Continuity::Discontinuous, |(class, _, _)| class);
                if is_zero(interval.min - natural.min) {
                    conts[0] = continuity;
                }
                if is_zero(interval.max - natural.max) {
                    conts[1] = continuity;
                }
            }
        }
        let mut angle = if line1.vector.length_squared() < EFF_ZERO_SQ
            || line2.vector.length_squared() < EFF_ZERO_SQ
        {
            PI / 2.0
        } else {
            line1.vector.angle_between(line2.vector)?
        };
        if angle > PI / 2.0 {
            angle = PI - angle;
        }
        let angle_deg = angle * (180.0 / PI);
        let combined_turn = span1.max_turning_angle_deg + span2.max_turning_angle_deg;
        let mut only_one = if combined_turn < EFF_ZERO * 100.0 {
            let length = min(line1.vector.length(), line2.vector.length());
            length * angle > tolerance
        } else {
            let scaled = 100.0
                * EFF_ZERO
                * (1.0 + max(line1.vector.max_dimension(), line2.vector.max_dimension()));
            angle_deg * 0.99 > combined_turn + scaled
        };
        if line1.vector.length_squared() >= EFF_ZERO_SQ
            && line2.vector.length_squared() >= EFF_ZERO_SQ
            && (conts1[0] == Continuity::Discontinuous || conts1[1] == Continuity::Discontinuous)
            && (conts2[0] == Continuity::Discontinuous || conts2[1] == Continuity::Discontinuous)
        {
            for j1 in 0..2 {
                for j2 in 0..2 {
                    if conts1[j1] != Continuity::Discontinuous
                        || conts2[j2] != Continuity::Discontinuous
                    {
                        continue;
                    }
                    let (point1, t1, mut direction1) = if j1 == 0 {
                        (line1.start, interval1.min, line1.vector)
                    } else {
                        (line1.end, interval1.max, -line1.vector)
                    };
                    let (point2, t2, mut direction2) = if j2 == 0 {
                        (line2.start, interval2.min, line2.vector)
                    } else {
                        (line2.end, interval2.max, -line2.vector)
                    };
                    let value = point1.distance(point2);
                    if value >= tolerance {
                        continue;
                    }
                    let mut test_one = false;
                    if only_one || (direction1.dot(direction2) < -EFF_ZERO && combined_turn < 90.0)
                    {
                        test_one = true;
                    } else if curve1.degree() <= 2 && curve2.degree() <= 2 && !only_one {
                        let pv1 = curve1.evaluate(t1, 1, true, true)?;
                        let pv2 = curve2.evaluate(t2, 1, true, true)?;
                        let (mut d1, mut d2) = (pv1[1], pv2[1]);
                        d1.unitize()?;
                        d2.unitize()?;
                        let cross = d1.cross(d2);
                        direction1.unitize()?;
                        direction2.unitize()?;
                        let cross2 = direction1.cross(direction2);
                        if cross.length_squared() < EFF_ZERO_SQ
                            && cross2.length_squared() > EFF_ZERO_SQ
                        {
                            test_one = true;
                        }
                    }
                    if !test_one {
                        continue;
                    }
                    let solution = point_solution(Some(branch), t1, t2, value);
                    if let Some((p1, p2)) = line_line_closest_point(
                        line1.start,
                        line1.vector,
                        line2.start,
                        line2.vector,
                    ) {
                        if (p1 < EFF_ZERO_SQRT || p1 > 1.0 - EFF_ZERO_SQRT)
                            && (p2 < EFF_ZERO_SQRT || p2 > 1.0 - EFF_ZERO_SQRT)
                        {
                            solver.add_sorted_solution(
                                &*self,
                                &solution,
                                SortKey::ByAllParameters,
                                true,
                            );
                            return Ok(false);
                        }
                    } else if direction1.dot(direction2) < EFF_ZERO {
                        let p1 = line_closest_point(line1.start, line1.vector, point2);
                        if p1 < EFF_ZERO || p1 > 1.0 - EFF_ZERO {
                            solver.add_sorted_solution(
                                &*self,
                                &solution,
                                SortKey::ByAllParameters,
                                true,
                            );
                            return Ok(false);
                        }
                    }
                }
            }
        }
        if !only_one {
            let found = simple_coincidence_check(
                curve1,
                interval1,
                curve2,
                interval2,
                tolerance,
                Some(&box1),
                Some(&box2),
            )?;
            let (pairs, range): (Vec<CurvePair>, bool) = match found {
                Coincidence::None => (Vec::new(), false),
                Coincidence::Point(pair) => (vec![pair], false),
                Coincidence::Points(pairs) => (pairs.to_vec(), false),
                Coincidence::Range(pairs) => (pairs.to_vec(), true),
            };
            if range {
                let mut solution = point_solution(
                    Some(branch),
                    pairs[0].parameters[0],
                    pairs[0].parameters[1],
                    pairs[0].value,
                );
                solution.kind = SolutionType::RangeOfValues;
                solution.end = end(
                    pairs[1].parameters[0],
                    pairs[1].parameters[1],
                    pairs[1].value,
                );
                solver.add_sorted_solution(&*self, &solution, SortKey::ByAllParameters, true);
                return Ok(false);
            }
            for pair in pairs {
                let mut solution = point_solution(
                    Some(branch),
                    pair.parameters[0],
                    pair.parameters[1],
                    pair.value,
                );
                if solution.start.value > tolerance / 100.0
                    && let Some(better) = local_curve_intersect(
                        curve1,
                        interval1,
                        curve2,
                        interval2,
                        tolerance,
                        [pair.parameters[0], pair.parameters[1]],
                    )?
                    && better.value < solution.start.value
                {
                    solution.start = end(better.parameters[0], better.parameters[1], better.value);
                }
                if filter_solution(
                    curve1, interval1, curve2, interval2, conts1, conts2, &solution,
                )? {
                    solver.add_sorted_solution(&*self, &solution, SortKey::ByAllParameters, true);
                }
            }
        }
        if !only_one {
            if leaf1.bbox.are_disjoint(leaf2.bbox, 0.0) {
                only_one = false;
            } else if let Ok(meet) = leaf1.bbox.intersect(leaf2.bbox) {
                let zero = EFF_ZERO * (1.0 + meet.max.max_dimension());
                let mut zero_dimensions = 0;
                if meet.max.x - meet.min.x < zero {
                    zero_dimensions += 1;
                }
                if meet.max.y - meet.min.y < zero {
                    zero_dimensions += 1;
                }
                if meet.max.z - meet.min.z < zero {
                    zero_dimensions += 1;
                }
                if zero_dimensions >= 3 {
                    only_one = true;
                }
            }
        }
        let mut found_end = false;
        let mut end_solution = point_solution(Some(branch), 0.0, 0.0, 0.0);
        let mut end_solutions: Vec<Solution> = Vec::new();
        // The four span ends against the other span. SMLib tests the first
        // span's start and the second's start for C0 or less, and the ends
        // for a discontinuity only.
        let ends = [
            (
                true,
                conts1[0] <= Continuity::C0,
                line1.start,
                interval1.min,
            ),
            (
                true,
                conts1[1] <= Continuity::Discontinuous,
                line1.end,
                interval1.max,
            ),
            (
                false,
                conts2[0] <= Continuity::C0,
                line2.start,
                interval2.min,
            ),
            (
                false,
                conts2[1] <= Continuity::Discontinuous,
                line2.end,
                interval2.max,
            ),
        ];
        for (on_first, open, point, t) in ends {
            let (target, target_interval, target_box, own_line) = if on_first {
                (curve2, interval2, &box2, &line1)
            } else {
                (curve1, interval1, &box1, &line2)
            };
            if !open || target_box.distance_to_point(point) >= tolerance {
                continue;
            }
            let Some(found) = local_point_solve(
                target,
                target_interval,
                SolverOperation::Minimize,
                point,
                target_interval.evaluate(0.5),
            )?
            else {
                continue;
            };
            if found.distance >= tolerance {
                continue;
            }
            if !(found.distance < tolerance / 100.0
                || !target_interval
                    .is_value_on_boundary(found.parameter, EFF_ZERO)
                    .0)
            {
                continue;
            }
            end_solution.start = if on_first {
                end(t, found.parameter, found.distance)
            } else {
                end(found.parameter, t, found.distance)
            };
            found_end = true;
            if only_one {
                let target_point = target.evaluate_point(found.parameter)?;
                if own_line.vector.length_squared() < EFF_ZERO_SQ {
                    solver.add_sorted_solution(
                        &*self,
                        &end_solution,
                        SortKey::ByAllParameters,
                        true,
                    );
                    return Ok(false);
                }
                let parameter = line_closest_point(own_line.start, own_line.vector, target_point);
                if parameter < EFF_ZERO || parameter > 1.0 - EFF_ZERO {
                    solver.add_sorted_solution(
                        &*self,
                        &end_solution,
                        SortKey::ByAllParameters,
                        true,
                    );
                    return Ok(false);
                }
            }
            end_solutions.push(end_solution);
        }
        let mut guesses: Vec<(f64, f64)> = if only_one {
            let (t1, t2) =
                line_line_closest_point(line1.start, line1.vector, line2.start, line2.vector)
                    .unwrap_or((0.5, 0.5));
            let unit = Extent1::new(0.0, 1.0);
            vec![(
                interval1.evaluate(unit.clamp_value(t1)),
                interval2.evaluate(unit.clamp_value(t2)),
            )]
        } else {
            vec![
                (interval1.evaluate(0.5), interval2.evaluate(0.5)),
                (interval1.evaluate(0.0), interval2.evaluate(0.0)),
                (interval1.evaluate(0.0), interval2.evaluate(1.0)),
                (interval1.evaluate(1.0), interval2.evaluate(0.0)),
                (interval1.evaluate(1.0), interval2.evaluate(1.0)),
            ]
        };
        if curve1.is_degenerate(EFF_ZERO) || curve2.is_degenerate(EFF_ZERO) {
            // SMLib adds the last end answer once for each end answer.
            for _ in 0..end_solutions.len() {
                solver.add_sorted_solution(&*self, &end_solution, SortKey::ByAllParameters, true);
            }
            return Ok(false);
        }
        let mut added = false;
        let mut i = 0;
        while i < guesses.len() {
            let (guess1, guess2) = guesses[i];
            let found = local_curve_intersect(
                curve1,
                interval1,
                curve2,
                interval2,
                tolerance,
                [guess1, guess2],
            )?;
            if let Some(found) = found {
                let (p1, p2) = (found.parameters[0], found.parameters[1]);
                let solution = point_solution(Some(branch), p1, p2, found.value);
                let keep = filter_solution(
                    curve1, interval1, curve2, interval2, conts1, conts2, &solution,
                )?;
                let skip = interval1.is_value_on_boundary(p1, EFF_ZERO).0
                    && interval2.is_value_on_boundary(p2, EFF_ZERO).0
                    && found_end
                    && found.value > tolerance / 100.0;
                if keep && !skip {
                    solver.add_sorted_solution(&*self, &solution, SortKey::ByAllParameters, true);
                    added = true;
                    if i < 5
                        && !only_one
                        && let Some(guess) =
                            additional_tangent_guess(curve1, interval1, curve2, interval2, p1, p2)?
                    {
                        guesses.push(guess);
                    }
                }
            }
            i += 1;
        }
        if !added || !only_one {
            for solution in &end_solutions {
                solver.add_sorted_solution(&*self, solution, SortKey::ByAllParameters, true);
            }
        }
        Ok(false)
    }
}

/// `sm_RefineCCSolutions`: moves the ends of each coincident range onto a
/// tighter meeting when it lies within a five-hundredth of the range.
fn refine_ranges(
    first: &BSplineCurve,
    second: &BSplineCurve,
    solutions: &mut SolutionArray,
) -> Result<(), SolverRefusal> {
    let natural1 = first.natural_interval();
    let natural2 = second.natural_interval();
    for solution in solutions.as_mut_slice() {
        if solution.kind != SolutionType::RangeOfValues {
            continue;
        }
        if solution.variables != 2 {
            return Err(Fault::InvalidInput.into());
        }
        let mut range1 = Extent1::default();
        range1.add_value(solution.start.parameters[0]);
        range1.add_value(solution.end.parameters[0]);
        let reverse1 = solution.start.parameters[0] > solution.end.parameters[0];
        let mut range2 = Extent1::default();
        range2.add_value(solution.start.parameters[1]);
        range2.add_value(solution.end.parameters[1]);
        let reverse2 = solution.start.parameters[1] > solution.end.parameters[1];
        for which in [0, 1] {
            let point = if which == 0 {
                solution.start
            } else {
                solution.end
            };
            if point.value <= 0.000_001 {
                continue;
            }
            let Some(found) = local_curve_intersect(
                first,
                natural1,
                second,
                natural2,
                1.0e-10,
                [point.parameters[0], point.parameters[1]],
            )?
            else {
                continue;
            };
            if found.value >= point.value {
                continue;
            }
            let (Some(mut f1), Some(mut f2)) = (
                inversion(range1, found.parameters[0]),
                inversion(range2, found.parameters[1]),
            ) else {
                // SMLib's failed inversion leaves the fraction undefined,
                // which fails the shift test.
                continue;
            };
            if reverse1 {
                f1 = 1.0 - f1;
            }
            if reverse2 {
                f2 = 1.0 - f2;
            }
            if f1.abs() < 0.002 && f2.abs() < 0.002 {
                let target = if which == 0 {
                    &mut solution.start
                } else {
                    &mut solution.end
                };
                target.value = found.value;
                target.parameters[0] = found.parameters[0];
                target.parameters[1] = found.parameters[1];
            }
        }
    }
    Ok(())
}

/// `SmExtent1d::Inversion`'s fraction; SMLib computes it for values outside
/// the interval too, failing only for a zero-length interval.
fn inversion(interval: Extent1, value: f64) -> Option<f64> {
    let numerator = value - interval.min;
    let denominator = interval.max - interval.min;
    if denominator < numerator.abs() * EFF_ZERO {
        None
    } else {
        Some(numerator / denominator)
    }
}

/// `SmCurve::ApproximatePolarBox`: a polar box of the tangents at the
/// middle and at `count` (made odd) evenly spaced parameters.
fn approximate_polar_box(
    curve: &BSplineCurve,
    interval: Extent1,
    count: usize,
) -> Result<PolarBox, Fault> {
    let count = if count.is_multiple_of(2) {
        count + 1
    } else {
        count
    };
    let natural = curve.natural_interval();
    let interval = Extent1 {
        min: natural.clamp_value(interval.min),
        max: natural.clamp_value(interval.max),
    };
    let mut polar = PolarBox::default();
    let middle = curve.evaluate(interval.evaluate(0.5), 1, true, true)?;
    polar.add_vector(middle[1], None);
    let mut parameter = interval.min;
    let step = interval.length() / (count - 1) as f64;
    for _ in 0..count {
        let pv = curve.evaluate(parameter, 1, true, true)?;
        polar.add_vector(pv[1], None);
        parameter += step;
    }
    Ok(polar)
}

/// Whether a curve needs splitting before the end-point coincidence test:
/// closed over `interval`, or turning by half a circle or more.
fn needs_split(curve: &BSplineCurve, interval: Extent1, tolerance: f64) -> Result<bool, Fault> {
    if is_closed(curve, interval, tolerance) {
        return Ok(true);
    }
    if interval.length() > ZONE_TOL_3D / 10.0 {
        let polar = approximate_polar_box(curve, interval, 5)?;
        let size = max(polar.max.x - polar.min.x, polar.max.y - polar.min.y);
        if polar.state == PolarState::Unbounded
            || (matches!(polar.state, PolarState::Arc | PolarState::Region)
                && size >= PI - EFF_ZERO)
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// `sm_IsPointSolutionAcrossClosedSeamBoundary`: a point answer at an end of
/// a closed curve whose tangents (turned inward) point apart.
fn across_closed_seam(
    solution: &Solution,
    this: &BSplineCurve,
    interval: Extent1,
    other: &BSplineCurve,
    other_interval: Extent1,
    this_is_first: bool,
) -> bool {
    let (start1, start2) = if this_is_first {
        (solution.start.parameters[0], solution.start.parameters[1])
    } else {
        (solution.start.parameters[1], solution.start.parameters[0])
    };
    if solution.kind == SolutionType::SingleValue
        && is_closed(this, interval, 0.0)
        && interval.is_value_on_boundary(start1, EFF_ZERO).0
    {
        let (Ok(pv1), Ok(pv2)) = (
            this.evaluate(start1, 1, true, true),
            other.evaluate(start2, 1, true, true),
        ) else {
            return false;
        };
        let mut d1 = pv1[1];
        let mut d2 = pv2[1];
        if are_same(start1, interval.max) {
            d1 = -d1;
        }
        if are_same(start2, other_interval.max) {
            d2 = -d2;
        }
        if d1.dot(d2) < 0.0 {
            return true;
        }
    }
    false
}

/// What [`global_coincidence_check`] found.
#[derive(Clone, Debug, Default)]
pub struct CoincidenceCheck {
    /// The meetings at curve ends and the coincident range, if any.
    pub solutions: SolutionArray,
    /// `rbAdditionalWorkNeeded`: whether other meetings may remain.
    pub more_work: bool,
}

/// `SmCurve::GlobalCoincidenceChecker` (the version without
/// `SM_TOLERANT_XSECT`): whether the curves coincide between their end
/// points dropped onto each other, splitting closed or strongly turning
/// curves in three first.
#[allow(clippy::too_many_lines)]
pub fn global_coincidence_check(
    this: &BSplineCurve,
    interval: Extent1,
    other: &BSplineCurve,
    other_interval: Extent1,
    tolerance: f64,
    snap_to_ends: bool,
) -> Result<CoincidenceCheck, SolverRefusal> {
    let mut result = CoincidenceCheck::default();
    if degenerate_on_to(this, interval, tolerance / 2.0)
        || degenerate_on_to(other, other_interval, tolerance / 2.0)
    {
        return Ok(result);
    }
    if let (Ok(this_box), Ok(other_box)) = (
        this.bounding_box(interval),
        other.bounding_box(other_interval),
    ) && this_box.are_disjoint(other_box, tolerance)
    {
        return Ok(result);
    }
    let objects = CurvePairObjects {
        first: this,
        second: other,
    };
    for split_other in [false, true] {
        let split = if split_other {
            needs_split(other, other_interval, tolerance)?
        } else {
            needs_split(this, interval, tolerance)?
        };
        if !split {
            continue;
        }
        let whole = if split_other {
            other_interval
        } else {
            interval
        };
        let pieces = [
            Extent1 {
                min: whole.min,
                max: whole.evaluate(0.375),
            },
            Extent1 {
                min: whole.evaluate(0.375),
                max: whole.evaluate(0.625),
            },
            Extent1 {
                min: whole.evaluate(0.625),
                max: whole.max,
            },
        ];
        let empty = SearchTrees::new();
        let mut gather = GlobalSolver::with_two_trees(&empty, &objects, None, None, 2);
        for piece in pieces {
            // SMLib's recursion snaps to the ends (the default argument).
            let part = if split_other {
                global_coincidence_check(this, interval, other, piece, tolerance, true)?
            } else {
                global_coincidence_check(this, piece, other, other_interval, tolerance, true)?
            };
            for solution in part.solutions.as_slice() {
                let skip = if split_other {
                    across_closed_seam(solution, this, interval, other, other_interval, true)
                } else {
                    across_closed_seam(solution, other, other_interval, this, interval, false)
                };
                if skip {
                    continue;
                }
                gather.add_sorted_solution(
                    &crate::global_solver::BaseHooks,
                    solution,
                    SortKey::ByAllParameters,
                    true,
                );
            }
        }
        result.solutions = std::mem::take(&mut gather.solutions);
        let mut check1 = Extent1::default();
        let mut check2 = Extent1::default();
        for solution in result.solutions.as_slice() {
            if solution.kind != SolutionType::RangeOfValues {
                continue;
            }
            let negative = check1.min > check1.max;
            if (negative
                && !are_same(solution.start.parameters[0], interval.min)
                && !(other_interval
                    .is_value_on_boundary(solution.start.parameters[1], EFF_ZERO)
                    .0
                    || other_interval
                        .is_value_on_boundary(solution.end.parameters[1], EFF_ZERO)
                        .0))
                || (!negative && !are_same(solution.start.parameters[0], check1.max))
            {
                break;
            }
            check1.add_value(solution.start.parameters[0]);
            check1.add_value(solution.end.parameters[0]);
            check2.add_value(solution.start.parameters[1]);
            check2.add_value(solution.end.parameters[1]);
        }
        result.more_work = !((check1.length() - interval.length()).abs() < EFF_ZERO
            || (check2.length() - other_interval.length()).abs() < EFF_ZERO);
        return Ok(result);
    }
    let start1 = this.evaluate(interval.min, 1, true, true)?;
    let mut end1 = this.evaluate(interval.max, 1, true, true)?;
    end1[1] = -end1[1];
    let start2 = other.evaluate(other_interval.min, 1, true, true)?;
    let mut end2 = other.evaluate(other_interval.max, 1, true, true)?;
    end2[1] = -end2[1];
    for (a, ta, b, tb) in [
        (&start1, interval.min, &start2, other_interval.min),
        (&start1, interval.min, &end2, other_interval.max),
        (&end1, interval.max, &start2, other_interval.min),
        (&end1, interval.max, &end2, other_interval.max),
    ] {
        let distance = a[0].distance(b[0]);
        if distance < tolerance && a[1].dot(b[1]) < 0.0 {
            result.solutions.add(point_solution(None, ta, tb, distance));
        }
    }
    let mut t1: Vec<f64> = Vec::new();
    let mut t2: Vec<f64> = Vec::new();
    let mut deviations: Vec<f64> = Vec::new();
    let snap = |v: Vector3| if snap_to_ends { Some(v) } else { None };
    let drop =
        |curve: &BSplineCurve, interval: Extent1, point: Vector3, inward: Option<Vector3>| {
            drop_point(curve, interval, point, inward, tolerance, None, INTERSECT)
        };
    let dropped = drop(other, other_interval, start1[0], snap(start1[1]))?;
    if dropped.success && dropped.distance < tolerance {
        t1.push(interval.min);
        t2.push(dropped.parameter);
        deviations.push(dropped.distance);
    }
    let dropped = drop(this, interval, start2[0], snap(start2[1]))?;
    if dropped.success && dropped.distance < tolerance {
        t1.push(dropped.parameter);
        t2.push(other_interval.min);
        deviations.push(dropped.distance);
    }
    let dropped = drop(this, interval, end2[0], snap(end2[1]))?;
    if dropped.success && dropped.distance < tolerance {
        let n = t1.len();
        if n > 0 && dropped.parameter < t1[n - 1] {
            t1.insert(n - 1, dropped.parameter);
            t2.insert(n - 1, other_interval.max);
            deviations.insert(n - 1, dropped.distance);
        } else {
            t1.push(dropped.parameter);
            t2.push(other_interval.max);
            deviations.push(dropped.distance);
        }
    }
    if t1.is_empty() {
        result.more_work = true;
        return Ok(result);
    }
    let dropped = drop(other, other_interval, end1[0], snap(end1[1]))?;
    if dropped.success && dropped.distance < tolerance {
        t1.push(interval.max);
        t2.push(dropped.parameter);
        deviations.push(dropped.distance);
    }
    let n = t1.len();
    let zero1 = EFF_ZERO * (1.0 + max(t1[n - 1].abs(), t1[0].abs()));
    let zero2 = EFF_ZERO * (1.0 + max(t2[n - 1].abs(), t2[0].abs()));
    if n == 1 || ((t1[n - 1] - t1[0]).abs() < zero1 || (t2[n - 1] - t2[0]).abs() < zero2) {
        result.more_work = true;
        return Ok(result);
    }
    let range1 = Extent1 {
        min: t1[0],
        max: t1[n - 1],
    };
    let degree = this.degree().max(other.degree());
    let between = curve_max_distance_between(
        this,
        range1,
        other,
        t2[0],
        t2[n - 1],
        degree + 3,
        Some(tolerance),
    )?;
    if between.distance > tolerance {
        result.more_work = true;
        return Ok(result);
    }
    let mut range = point_solution(None, t1[0], t2[0], deviations[0]);
    range.kind = SolutionType::RangeOfValues;
    range.end = end(t1[n - 1], t2[n - 1], deviations[n - 1]);
    result.solutions.add(range);
    result.more_work = true;
    let low = min(t1[0], t1[n - 1]);
    let high = max(t1[0], t1[n - 1]);
    let delta0 = (low - interval.min).abs();
    let delta1 = (high - interval.max).abs();
    let d0 = if delta0 > 0.0 {
        delta0 * start1[1].length()
    } else {
        0.0
    };
    let d1 = if delta1 > 0.0 {
        delta1 * end1[1].length()
    } else {
        0.0
    };
    if d0 < tolerance && d1 < tolerance {
        result.more_work = false;
    }
    if result.more_work {
        let low = min(t2[0], t2[n - 1]);
        let high = max(t2[0], t2[n - 1]);
        let delta0 = (low - other_interval.min).abs();
        let delta1 = (high - other_interval.max).abs();
        let d0 = if delta0 > 0.0 {
            delta0 * start2[1].length()
        } else {
            0.0
        };
        let d1 = if delta1 > 0.0 {
            delta1 * end2[1].length()
        } else {
            0.0
        };
        if d0 < tolerance && d1 < tolerance {
            result.more_work = false;
        }
    }
    Ok(result)
}

/// `SmCurve::ClassifyIntersectingRangeSolution`: whether a range answer is
/// really a tangent meeting at one point (the gap between the curves grows
/// steadily away from the range's tighter end), and that point.
fn classify_range(
    this: &BSplineCurve,
    other: &BSplineCurve,
    solution: &Solution,
) -> Result<Option<SolutionEnd>, SolverRefusal> {
    if solution.kind == SolutionType::SingleValue {
        return Ok(Some(solution.start));
    }
    let natural = this.natural_interval();
    let other_natural = other.natural_interval();
    let (start_deviation, end_deviation, start, delta) =
        if solution.start.value < solution.end.value {
            (
                solution.start.value,
                solution.end.value,
                [solution.start.parameters[0], solution.start.parameters[1]],
                [
                    solution.end.parameters[0] - solution.start.parameters[0],
                    solution.end.parameters[1] - solution.start.parameters[1],
                ],
            )
        } else {
            (
                solution.end.value,
                solution.start.value,
                [solution.end.parameters[0], solution.end.parameters[1]],
                [
                    solution.start.parameters[0] - solution.end.parameters[0],
                    solution.start.parameters[1] - solution.end.parameters[1],
                ],
            )
        };
    let _ = end_deviation;
    let mut this_range = Extent1::default();
    this_range.add_value(solution.start.parameters[0]);
    this_range.add_value(solution.end.parameters[0]);
    let mut other_range = Extent1::default();
    other_range.add_value(solution.start.parameters[1]);
    other_range.add_value(solution.end.parameters[1]);
    let Some(minimum) = local_curve_solve(
        this,
        this_range,
        other,
        other_range,
        SolverOperation::Minimize,
        1.0,
        [this_range.evaluate(0.5), other_range.evaluate(0.5)],
    )?
    else {
        return Ok(None);
    };
    if minimum.value > start_deviation / 10.0 {
        return Ok(None);
    }
    let point0 = this.evaluate_point(start[0])?;
    let point1 = other.evaluate_point(start[1])?;
    let meeting_gap = point1 - point0;
    let zero = EFF_ZERO * (1.0 + point0.max_dimension());
    let mut gaps = [Vector3::default(); 3];
    let mut lengths = [0.0; 3];
    let (mut is_zero_gap, mut monotone, mut diverging) = (true, true, true);
    let mut degenerate = false;
    let mut done = false;
    let mut step = 1.0 / 16.0001;
    let mut ii = 1;
    while !done && step < 1.0 {
        let tgt = if ii > 1 { 2 } else { 1 };
        let t1 = start[0] + step * delta[0];
        let t2 = start[1] + step * delta[1];
        let curve_point = this.evaluate_point(t1)?;
        match local_point_solve(
            other,
            other_natural,
            SolverOperation::Normalize,
            curve_point,
            t2,
        )? {
            Some(found) => {
                let other_point = other.evaluate_point(found.parameter)?;
                gaps[tgt] = (other_point - curve_point) - meeting_gap;
                lengths[tgt] = gaps[tgt].length();
            }
            None => {
                degenerate = true;
                done = true;
                ii += 1;
                step *= 2.0;
                continue;
            }
        }
        if ii >= 2 {
            let tgt1 = (tgt + 3 - 1) % 3;
            let tgt2 = (tgt + 3 - 2) % 3;
            if lengths[tgt] > zero {
                is_zero_gap = false;
            }
            if gaps[tgt].dot(gaps[tgt1]) <= 0.0 {
                monotone = false;
            }
            if lengths[tgt] == 0.0
                || lengths[tgt] < lengths[tgt1] + 1.0001 * (lengths[tgt1] - lengths[tgt2])
            {
                diverging = false;
            }
            if !is_zero_gap && monotone && diverging {
                degenerate = true;
                done = true;
            } else if is_zero_gap && monotone && diverging {
                done = false;
                gaps[1] = gaps[2];
                lengths[1] = lengths[2];
            } else {
                degenerate = false;
                done = true;
            }
        }
        ii += 1;
        step *= 2.0;
    }
    if !done && diverging {
        degenerate = delta[0].abs() < natural.length() / 10.0
            && delta[1].abs() < other_natural.length() / 10.0;
    }
    if degenerate {
        Ok(Some(end(start[0], start[1], start_deviation)))
    } else {
        Ok(None)
    }
}

/// `SmCurve::GlobalCurveIntersect` for two B-spline curves: their meetings
/// within `tolerance` over the intervals, as points (`t1`, `t2`, gap) and
/// coincident ranges, sorted along the first curve. With
/// `skip_coincidence`, the whole-curve coincidence check is skipped.
#[allow(clippy::too_many_lines)]
pub fn global_curve_intersect(
    this: &BSplineCurve,
    interval: Extent1,
    other: &BSplineCurve,
    other_interval: Extent1,
    tolerance: f64,
    skip_coincidence: bool,
) -> Result<SolutionArray, SolverRefusal> {
    if std::ptr::eq(this, other) {
        return Err(SolverRefusal::NotPorted(INTERSECT));
    }
    let cache1 = CurveCache::for_curve(this, 0.5 * ZONE_TOL_3D)?;
    let cache2 = CurveCache::for_curve(other, 0.5 * ZONE_TOL_3D)?;
    let mut trees = SearchTrees::new();
    trees.append(&cache1.trees);
    let offset = trees.append(&cache2.trees);
    let top1 = cache1.trees.get(cache1.top).ok_or(Fault::Degenerate)?.bbox;
    let top2 = cache2.trees.get(cache2.top).ok_or(Fault::Degenerate)?.bbox;
    let mut root = top1;
    root.expand_absolute(tolerance);
    if root.are_disjoint(top2, 0.0) {
        return Ok(SolutionArray::new());
    }
    if may_be_analytic(this, cache1.interval)? && may_be_analytic(other, cache2.interval)? {
        // `AnalyticCurveIntersect` is not ported.
        return Err(SolverRefusal::NotPorted(INTERSECT));
    }
    let objects = CurvePairObjects {
        first: this,
        second: other,
    };
    let tops = (Some(cache1.top), Some(NodeId(cache2.top.0 + offset)));
    let mut hooks = CurveCurveHooks {
        first: this,
        second: other,
        first_cache: &cache1,
        second_cache: &cache2,
        second_offset: offset,
        interval1: interval,
        interval2: other_interval,
    };
    let mut solutions = SolutionArray::new();
    {
        let mut test = GlobalSolver::with_two_trees(&trees, &objects, tops.0, tops.1, 2);
        test.solve_it(
            &mut hooks,
            SolverOperation::IntersectionTest,
            SolutionRequested::All,
            tolerance,
            BIG_DOUBLE,
            &[],
            &mut solutions,
        )?;
        if solutions.is_empty() {
            return Ok(solutions);
        }
        let first = solutions.as_slice()[0];
        if first.start.parameters[6] < 1.00001 {
            test.operation = INTERSECT;
            test.solutions = SolutionArray::new();
            let branch: Branch = {
                let mut b: Branch = [None; crate::global_solver::MAX_TREES];
                b[0] = first.nodes[0];
                b[1] = first.nodes[1];
                b
            };
            hooks.local_solve(&mut test, &branch)?;
            solutions = std::mem::take(&mut test.solutions);
            let mut finished = true;
            if solutions.len() == 1 {
                let only = &mut solutions.as_mut_slice()[0];
                only.object_count = 2;
                only.objects[0] = Some(FIRST_CURVE);
                only.objects[1] = Some(SECOND_CURVE);
                if only.kind == SolutionType::SingleValue && only.start.value > tolerance / 1000.0 {
                    let pv1 = this.evaluate(only.start.parameters[0], 1, true, true)?;
                    let pv2 = other.evaluate(only.start.parameters[1], 1, true, true)?;
                    let degrees = pv1[1]
                        .angle_between(pv2[1])
                        .map_or(90.0, |a| a * 180.0 / PI);
                    if !(20.0..=160.0).contains(&degrees) {
                        finished = false;
                    }
                }
            } else if solutions.len() > 1 {
                finished = false;
            }
            if finished {
                return Ok(solutions);
            }
        }
        solutions = SolutionArray::new();
    }
    if !skip_coincidence {
        let check =
            global_coincidence_check(this, interval, other, other_interval, tolerance, false)?;
        solutions = check.solutions;
        if !check.more_work {
            return Ok(solutions);
        }
    }
    let mut interval1 = interval;
    let mut interval2 = other_interval;
    for solution in solutions.as_slice() {
        if solution.kind != SolutionType::RangeOfValues {
            continue;
        }
        if are_same(interval1.min, solution.start.parameters[0])
            && solution.end.parameters[0] < interval1.max
        {
            let _ = interval1.set_min_max(solution.end.parameters[0], interval1.max);
            if interval1.length().abs() < EFF_ZERO {
                return Ok(solutions);
            }
        }
        if are_same(interval1.max, solution.end.parameters[0])
            && interval1.min < solution.start.parameters[0]
        {
            let _ = interval1.set_min_max(interval1.min, solution.start.parameters[0]);
            if interval1.length().abs() < EFF_ZERO {
                return Ok(solutions);
            }
        }
        let mut other_range = Extent1 {
            min: solution.start.parameters[1],
            max: solution.start.parameters[1],
        };
        other_range.add_value(solution.end.parameters[1]);
        if are_same(interval2.min, other_range.min) && other_range.max < interval2.max {
            let _ = interval2.set_min_max(other_range.max, interval2.max);
            if interval2.length().abs() < EFF_ZERO {
                return Ok(solutions);
            }
        }
        if are_same(interval2.max, other_range.max) && interval2.min < other_range.min {
            // SMLib sets the first curve's interval here.
            let _ = interval1.set_min_max(interval2.min, other_range.min);
            if interval2.length().abs() < EFF_ZERO {
                return Ok(solutions);
            }
        }
    }
    let kept: Vec<Solution> = solutions
        .as_slice()
        .iter()
        .copied()
        .filter(|s| !(s.kind == SolutionType::SingleValue && s.start.value > EFF_ZERO))
        .collect();
    solutions = SolutionArray::new();
    for solution in kept {
        solutions.add(solution);
    }
    hooks.interval1 = interval1;
    hooks.interval2 = interval2;
    let mut solver = GlobalSolver::with_two_trees(&trees, &objects, tops.0, tops.1, 2);
    solver.solve_it(
        &mut hooks,
        INTERSECT,
        SolutionRequested::All,
        tolerance,
        BIG_DOUBLE,
        &[],
        &mut solutions,
    )?;
    refine_ranges(this, other, &mut solutions)?;
    combine(
        this,
        interval,
        other,
        other_interval,
        interval2,
        tolerance,
        &hooks,
        &mut solutions,
    )?;
    Ok(solutions)
}

/// The interval of the span a solution's tree node holds, or the point.
fn node_interval(
    hooks: &CurveCurveHooks<'_>,
    node: Option<NodeId>,
    second: bool,
    value: f64,
) -> Extent1 {
    node.and_then(|id| hooks.leaf(id, second).ok())
        .map(|leaf| {
            if second {
                hooks.second_cache.spans[leaf.index].interval
            } else {
                hooks.first_cache.spans[leaf.index].interval
            }
        })
        .unwrap_or(Extent1 {
            min: value,
            max: value,
        })
}

/// The tail of `GlobalCurveIntersect`: merges neighbouring answers of the
/// same type, culls points inside or next to ranges, drops the worse of
/// close point pairs and turns short or tangent ranges into points.
#[allow(clippy::too_many_lines, clippy::too_many_arguments)]
fn combine(
    this: &BSplineCurve,
    interval: Extent1,
    other: &BSplineCurve,
    other_interval: Extent1,
    trimmed_other: Extent1,
    tolerance: f64,
    hooks: &CurveCurveHooks<'_>,
    solutions: &mut SolutionArray,
) -> Result<(), SolverRefusal> {
    const NONE: usize = usize::MAX;
    let mut items: Vec<Solution> = solutions.as_slice().to_vec();
    let (mut last_range, mut last_point) =
        if !items.is_empty() && items[0].kind == SolutionType::RangeOfValues {
            (0, NONE)
        } else {
            (NONE, 0)
        };
    let this_whole = this.natural_interval();
    let other_whole = other.natural_interval();
    for ii in 1..items.len() {
        let current = items[ii];
        let is_range = current.kind == SolutionType::RangeOfValues;
        let last_index = if is_range { last_range } else { last_point };
        let combine = if last_index == NONE {
            false
        } else {
            let last = items[last_index];
            let last_end = if last.kind == SolutionType::RangeOfValues {
                last.end
            } else {
                last.start
            };
            let le = last_end.parameters;
            let cs = current.start.parameters;
            // Answers at one parameter of one curve far apart on the other
            // stay apart (SMLib's two separate tests).
            if ((le[0] - cs[0]).abs() < EFF_ZERO_SQRT
                && (le[1] - cs[1]).abs() * 10.0 > other_interval.max - other_interval.min)
                || ((le[1] - cs[1]).abs() < EFF_ZERO_SQRT
                    && (le[0] - cs[0]).abs() * 10.0 > interval.max - interval.min)
            {
                false
            } else if are_same(le[0], cs[0]) || are_same(le[1], cs[1]) {
                let mut combine = true;
                if is_range {
                    let (start1, end1) = (last.start.value, last.end.value);
                    if start1 > EFF_ZERO && start1 > end1 * 100_000.0 {
                        combine = false;
                    }
                    let (start2, end2) = (current.start.value, current.end.value);
                    if end2 > EFF_ZERO && end2 > start2 * 100_000.0 {
                        combine = false;
                    }
                    let last_length =
                        max(EFF_ZERO, last.end.parameters[0] - last.start.parameters[0]);
                    let current_length = max(
                        EFF_ZERO,
                        current.end.parameters[0] - current.start.parameters[0],
                    );
                    if last_length > current_length * 1000.0
                        || current_length > last_length * 1000.0
                    {
                        combine = true;
                    }
                }
                combine
            } else if !is_range {
                let p1 = this.evaluate_point(cs[0])?;
                let p2 = this.evaluate_point(le[0])?;
                let mut combine = false;
                if p1.distance_squared(p2) < tolerance * tolerance {
                    let across = (cs[0] - le[0]).abs() > 0.9 * this_whole.length()
                        || (cs[1] - le[1]).abs() > 0.9 * other_whole.length();
                    if !across {
                        combine = true;
                    }
                }
                combine
            } else if le[0] >= cs[0] {
                is_between(cs[1], last.start.parameters[1], last.end.parameters[1])
            } else {
                let mut range1 = Extent1::default();
                range1.add_value(le[0]);
                range1.add_value(cs[0]);
                let mut range2 = Extent1::default();
                range2.add_value(le[1]);
                range2.add_value(cs[1]);
                let on1 = this.evaluate_point(range1.evaluate(0.45))?;
                let dropped = drop_point(
                    other,
                    range2,
                    on1,
                    None,
                    tolerance,
                    Some(range2.evaluate(0.5)),
                    SolverOperation::Minimize,
                )?;
                if !dropped.success {
                    false
                } else {
                    let on2 = other.evaluate_point(trimmed_other.evaluate(0.38))?;
                    let dropped = drop_point(
                        this,
                        range1,
                        on2,
                        None,
                        tolerance,
                        Some(range1.evaluate(0.5)),
                        SolverOperation::Minimize,
                    )?;
                    if !dropped.success {
                        false
                    } else if range1.length() < interval.length() / 1000.0
                        && range2.length() < other_interval.length() / 1000.0
                    {
                        true
                    } else {
                        let degree = this.degree().max(other.degree());
                        let between = curve_max_distance_between(
                            this,
                            range1,
                            other,
                            le[1],
                            cs[1],
                            degree + 3,
                            Some(tolerance),
                        )?;
                        between.distance < tolerance
                    }
                }
            }
        };
        if combine {
            let last = &mut items[last_index];
            if is_range {
                if last.end.parameters[0] < current.end.parameters[0] {
                    last.end = current.end;
                }
            } else {
                let last_end = if last.kind == SolutionType::RangeOfValues {
                    last.end
                } else {
                    last.start
                };
                if current.start.value < last_end.value {
                    *last = current;
                }
            }
        } else {
            let next = if last_range == NONE {
                last_point.wrapping_add(1)
            } else if last_point == NONE {
                last_range + 1
            } else {
                last_range.max(last_point) + 1
            };
            if is_range {
                last_range = next;
            } else {
                last_point = next;
            }
            if next != ii {
                items[next] = current;
            }
        }
    }
    if !items.is_empty() {
        let last = if last_range == NONE {
            last_point
        } else if last_point == NONE {
            last_range
        } else {
            last_range.max(last_point)
        };
        items.truncate(last + 1);
    }
    // Cull points inside or beside ranges, and ranges too short to keep.
    let count = items.len();
    let mut cull = vec![false; count];
    let mut compress = false;
    for ii in 0..count {
        let range = items[ii];
        if range.kind != SolutionType::RangeOfValues {
            continue;
        }
        let mut span1 = Extent1::default();
        span1.add_value(range.start.parameters[0]);
        span1.add_value(range.end.parameters[0]);
        let mut span2 = Extent1::default();
        span2.add_value(range.start.parameters[1]);
        span2.add_value(range.end.parameters[1]);
        for jj in 0..count {
            if ii == jj {
                continue;
            }
            let point = items[jj];
            if point.kind != SolutionType::SingleValue {
                continue;
            }
            let p = point.start.parameters;
            let node1 = node_interval(hooks, point.nodes[0], false, p[0]);
            let node2 = node_interval(hooks, point.nodes[1], true, p[1]);
            if node1.intersect(span1).is_err() || node2.intersect(span2).is_err() {
                continue;
            }
            let close = interval.length() / 1000.0;
            if span1.length() < close {
                cull[ii] = true;
                compress = true;
            } else if p[0] < span1.min {
                if span1.min - p[0] < close && point.start.value > 2.0 * range.start.value {
                    cull[jj] = true;
                    compress = true;
                }
            } else if p[0] > span1.max {
                if p[0] - span1.max < close && point.start.value > 2.0 * range.end.value {
                    cull[jj] = true;
                    compress = true;
                }
            } else {
                cull[jj] = true;
                compress = true;
            }
            if cull[ii] {
                break;
            }
            if p[0] > span1.max || p[1] > span2.max {
                break;
            }
        }
    }
    if compress {
        let mut kept: Vec<Solution> = Vec::new();
        let mut ii = 0;
        while ii < items.len() {
            if cull[ii] {
                if items[ii].kind == SolutionType::RangeOfValues {
                    let (start, finish) = (items[ii].start.value, items[ii].end.value);
                    if start < EFF_ZERO && finish < EFF_ZERO {
                        let mut copy = items[ii];
                        copy.kind = SolutionType::SingleValue;
                        items[ii].kind = SolutionType::SingleValue;
                        items[ii].start = items[ii].end;
                        items.insert(ii + 1, copy);
                        cull[ii] = false;
                        cull.insert(ii + 1, false);
                    } else if start < EFF_ZERO {
                        items[ii].kind = SolutionType::SingleValue;
                        cull[ii] = false;
                    } else if finish < EFF_ZERO {
                        items[ii].kind = SolutionType::SingleValue;
                        items[ii].start = items[ii].end;
                        cull[ii] = false;
                    } else {
                        ii += 1;
                        continue;
                    }
                } else {
                    ii += 1;
                    continue;
                }
            }
            kept.push(items[ii]);
            ii += 1;
        }
        items = kept;
    }
    // Of two close point answers along both curves whose gaps differ a
    // hundredfold, drop the worse when the curves part between them.
    let delta_a = this_whole.length() / 100.0;
    let delta_b = other_whole.length() / 100.0;
    let count = items.len();
    for ii in 0..count.saturating_sub(1) {
        if items[ii].kind != SolutionType::SingleValue || items[ii].start.value > BIG_DOUBLE / 2.0 {
            continue;
        }
        let (a1, b1) = (items[ii].start.parameters[0], items[ii].start.parameters[1]);
        for jj in ii + 1..count {
            if items[jj].kind != SolutionType::SingleValue
                || items[jj].start.value > BIG_DOUBLE / 2.0
            {
                continue;
            }
            let (a2, b2) = (items[jj].start.parameters[0], items[jj].start.parameters[1]);
            if (a2 - a1).abs() > delta_a || (b2 - b1).abs() > delta_b {
                break;
            }
            let (value_a, value_b) = (items[ii].start.value, items[jj].start.value);
            if value_a < value_b * 100.0 && value_a > value_b / 100.0 {
                continue;
            }
            let distance1 = this.evaluate_point(a1)?.distance(other.evaluate_point(b1)?);
            let distance2 = this.evaluate_point(a2)?.distance(other.evaluate_point(b2)?);
            let middle = this
                .evaluate_point((a1 + a2) / 2.0)?
                .distance(other.evaluate_point((b1 + b2) / 2.0)?);
            if distance1 < distance2 && middle < distance2 {
                items[jj].start.value = BIG_DOUBLE;
            } else if distance2 < distance1 && middle < distance1 {
                items[ii].start.value = BIG_DOUBLE;
            }
        }
    }
    items.retain(|s| s.start.value <= BIG_DOUBLE / 2.0);
    // Ranges that are short, or a tangent meeting at one point, become
    // points.
    let mut this_length = BIG_DOUBLE;
    let mut other_length = BIG_DOUBLE;
    // SMLib declares the flag once for the loop: a short range keeps the
    // previous range's classification.
    let mut degenerate = false;
    for current in &mut items {
        if current.kind != SolutionType::RangeOfValues {
            continue;
        }
        let (start_deviation, end_deviation) = (current.start.value, current.end.value);
        let mut range = Extent1 {
            min: current.start.parameters[0],
            max: current.start.parameters[0],
        };
        range.add_value(current.end.parameters[0]);
        let mut other_range = Extent1 {
            min: current.start.parameters[1],
            max: current.start.parameters[1],
        };
        other_range.add_value(current.end.parameters[1]);
        let (length, other_range_length) = (range.length(), other_range.length());
        let point = this.evaluate_point(current.start.parameters[0])?;
        if this_length == BIG_DOUBLE {
            this_length = crate::curve::approximate_length(this, interval, 5);
        }
        if other_length == BIG_DOUBLE {
            other_length = crate::curve::approximate_length(other, other_interval, 5);
        }
        let segment = crate::curve::approximate_length(this, range, 5);
        let zero = EFF_ZERO * (1.0 + point.max_dimension());
        let tolerance1 = 10000.0 * EFF_ZERO * (1.0 + interval.max_dimension());
        let tolerance2 = 10000.0 * EFF_ZERO * (1.0 + other_interval.max_dimension());
        if start_deviation < zero
            && end_deviation < zero
            && length > tolerance1
            && other_range_length > tolerance2
        {
            continue;
        }
        let large = segment > min(this_length, other_length) / 10.0;
        let small = length < tolerance1 || other_range_length < tolerance2;
        if !small {
            degenerate = false;
            if let Some(point) = classify_range(this, other, current)? {
                degenerate = true;
                current.start = point;
            }
        }
        if large && !degenerate {
            continue;
        }
        if small || degenerate {
            if small && current.start.value < EFF_ZERO && current.end.value < EFF_ZERO {
                let s = current.start.parameters;
                let e = current.end.parameters;
                let delta_start0 = s[0] - interval.min;
                let delta_end0 = interval.max - e[0];
                let start_begin1 = s[1] - other_interval.min;
                let start_end1 = other_interval.max - s[1];
                let end_begin1 = e[1] - other_interval.min;
                let end_end1 = other_interval.max - e[1];
                let least = min(min(start_begin1, start_end1), min(end_begin1, end_end1));
                let at_start = delta_start0 < EFF_ZERO;
                let at_end = delta_end0 < EFF_ZERO;
                let use_end = !at_start && (at_end || end_begin1 == least || end_end1 == least);
                if use_end {
                    current.start = current.end;
                }
            } else if current.start.value > current.end.value {
                current.start = current.end;
            }
            current.kind = SolutionType::SingleValue;
        }
    }
    let mut combined = SolutionArray::new();
    for item in items {
        combined.add(item);
    }
    *solutions = combined;
    Ok(())
}

/// `SM_IS_BETWEEN(c, a, b)`: strictly inside either ordering of `a`, `b`
/// by their scaled zero.
fn is_between(c: f64, a: f64, b: f64) -> bool {
    let zero = max(scaled_zero(a), scaled_zero(b));
    (a + zero < c && c < b - zero) || (b + zero < c && c < a - zero)
}
