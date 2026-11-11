// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline curve and surface
// candidates of `spline.rs` without trimming faces; statuses become `Fault`,
// tree-node pointers become identities in one arena holding the curve and
// surface trees, and solution containers keep SMLib's `SmSolution` layout.
// Paths that need code from later port steps refuse with
// `SolverRefusal::NotPorted`: coincident curves (`DropAndTrimCurve`, `DropPoint` and the range
// classification of the solution combiner) and the polar-box cone test.
// Three answers deliberately differ from SMLib: a failed point search for a
// degenerate curve is a typed failure (SMLib ignores its status), its
// answers declare the three variables they carry (SMLib leaves two), and the
// end-point fallback drops points outside the requested domain (SMLib keeps
// them).
// Translated from SMLib SmSurface.cpp (SmSurface::GlobalCurveIntersect,
// SmCSIGlobalSolver, sm_CsArePerpenDisjoint,
// sm_GlobalCurveIntersectCombineSolutions); rust-brep has no port.
//! The global curve-surface intersector of step 6: the curve and surface
//! cache trees searched together, local intersections in each pair of
//! leaves, checks along discontinuous ends and sides, and the combination
//! of the answers.

use crate::Fault;
use crate::curve::{CurveSide, evaluate_continuity};
use crate::curve_cache::{CurveCache, NodeData};
use crate::curve_solve::{Coincidence, CurvePair, is_closed, simple_coincidence_check};
use crate::curve_surface::{
    Candidate, MeetingType, classify_meeting, drop_point_to_tangent_plane, filter_solution,
    is_singularity, local_curve_intersect,
};
use crate::extent::{Extent1, Extent2};
use crate::geometry::Continuity;
use crate::global_solver::{
    Branch, GlobalSolver, SearchTrees, SolutionRequested, SolverHooks, SolverOperation,
    SolverRefusal,
};
use crate::scalar::{ANGLE_TOL_DEG, BIG_DOUBLE, EFF_ZERO, PI, ZONE_TOL_3D, are_same, is_zero};
use crate::solution::{
    MAX_PARAMETERS, NodeId, ObjectId, ObjectKind, Preferred, Solution, SolutionArray, SolutionEnd,
    SolutionObjects, SolutionPair, SolutionType, SortKey,
};
use crate::spline::{BSplineCurve, BSplineSurface};
use crate::surface_cache::{CacheRefusal, SurfaceCache, SurfaceNodeData};
use crate::surface_eval::evaluate_point;
use crate::surface_solve::{IsoCurve, is_closed_with_continuity, is_on_seam};
use crate::vector::{Vector2, Vector3};

/// The solution identity of the curve.
pub const CURVE_OBJECT: ObjectId = ObjectId(0);
/// The solution identity of the surface.
pub const SURFACE_OBJECT: ObjectId = ObjectId(1);

const INTERSECT: SolverOperation = SolverOperation::Intersect;

/// The curve and the surface a curve-surface solution names.
pub struct CurveSurfaceObjects<'a> {
    /// The curve.
    pub curve: &'a BSplineCurve,
    /// The surface.
    pub surface: &'a BSplineSurface,
    /// The zone tolerance of both (SMLib's global one without an owner).
    pub zone: f64,
}

impl SolutionObjects for CurveSurfaceObjects<'_> {
    fn kind(&self, id: ObjectId) -> ObjectKind {
        match id {
            CURVE_OBJECT => ObjectKind::Curve,
            SURFACE_OBJECT => ObjectKind::Surface,
            _ => ObjectKind::Other,
        }
    }

    fn parameter_interval(&self, id: ObjectId, direction: usize) -> Extent1 {
        match (id, direction) {
            (CURVE_OBJECT, _) => self.curve.natural_interval(),
            (SURFACE_OBJECT, 0) => self.surface.natural_domain().u_interval(),
            (SURFACE_OBJECT, 1) => self.surface.natural_domain().v_interval(),
            _ => Extent1::default(),
        }
    }

    fn zone_tolerance(&self, id: ObjectId) -> f64 {
        if id == CURVE_OBJECT || id == SURFACE_OBJECT {
            self.zone
        } else {
            0.0
        }
    }
}

fn start(t: f64, u: f64, v: f64, value: f64) -> SolutionEnd {
    let mut parameters = [0.0; MAX_PARAMETERS];
    parameters[..3].copy_from_slice(&[t, u, v]);
    SolutionEnd { value, parameters }
}

/// A single-value curve-surface solution from the global solver.
fn point_solution(branch: &Branch, t: f64, uv: Vector2, value: f64) -> Solution {
    let mut solution = Solution {
        kind: SolutionType::SingleValue,
        variables: 3,
        start: start(t, uv.x, uv.y, value),
        object_count: 2,
        ..Solution::default()
    };
    solution.objects[0] = Some(CURVE_OBJECT);
    solution.objects[1] = Some(SURFACE_OBJECT);
    solution.nodes[0] = branch[0];
    solution.nodes[1] = branch[1];
    solution
}

/// `SmCSIGlobalSolver`.
struct CurveSurfaceHooks<'a> {
    curve: &'a BSplineCurve,
    surface: &'a BSplineSurface,
    domain: Extent2,
    interval: Extent1,
    coincidence: bool,
    curve_cache: &'a CurveCache,
    surface_cache: &'a SurfaceCache,
    surface_offset: u32,
    zone: f64,
    owner: Option<f64>,
    tolerance: f64,
}

impl CurveSurfaceHooks<'_> {
    fn candidate(solution: &Solution) -> Candidate {
        Candidate {
            t: solution.start.parameters[0],
            uv: Vector2::new(solution.start.parameters[1], solution.start.parameters[2]),
            value: solution.start.value,
        }
    }

    /// `sm_CsArePerpenDisjoint`: whether the curve node lies in a cone of
    /// directions the surface node's normals avoid.
    fn are_perpendicular_disjoint(
        &self,
        solver: &GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<bool, SolverRefusal> {
        let _ = solver;
        let (Some(curve_id), Some(surface_id)) = (branch[0], branch[1]) else {
            return Err(Fault::Degenerate.into());
        };
        let trees = self.trees_of(curve_id, surface_id)?;
        let Ok((_, cone_angle)) = trees.1.bbox.cone_bound_to(trees.0.bbox) else {
            // SMLib's SER returns the error status, read as true.
            return Ok(true);
        };
        if cone_angle > PI / 2.0 {
            return Ok(false);
        }
        Err(SolverRefusal::NotPorted(INTERSECT))
    }

    fn trees_of(
        &self,
        curve_id: NodeId,
        surface_id: NodeId,
    ) -> Result<
        (
            crate::global_solver::TreeNode,
            crate::global_solver::TreeNode,
        ),
        SolverRefusal,
    > {
        let curve = self
            .curve_cache
            .trees
            .get(curve_id)
            .copied()
            .ok_or(SolverRefusal::Failed(Fault::Degenerate))?;
        let surface = self
            .surface_cache
            .trees
            .get(NodeId(surface_id.0 - self.surface_offset))
            .copied()
            .ok_or(SolverRefusal::Failed(Fault::Degenerate))?;
        Ok((curve, surface))
    }
}

impl SolverHooks for CurveSurfaceHooks<'_> {
    #[allow(clippy::too_many_lines)]
    fn local_solve(
        &mut self,
        solver: &mut GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<bool, SolverRefusal> {
        let tolerance = solver.tolerance;
        let (Some(curve_id), Some(surface_id)) = (branch[0], branch[1]) else {
            return Err(Fault::Degenerate.into());
        };
        let (curve_node, surface_node) = self.trees_of(curve_id, surface_id)?;
        let NodeData::Leaf(span_index) = self.curve_cache.node_data[curve_id.0 as usize] else {
            return Err(Fault::Degenerate.into());
        };
        let SurfaceNodeData::Patch(patch_index) =
            self.surface_cache.node_data[(surface_id.0 - self.surface_offset) as usize]
        else {
            return Err(Fault::Degenerate.into());
        };
        let span = self.curve_cache.spans[span_index];
        let span_box = self.curve_cache.pseudo_boxes[span_index];
        let patch = &self.surface_cache.patches[patch_index];
        let (curve, surface) = (self.curve, self.surface);
        let mut surface_domain = patch.uv_domain;
        let (patch_u, patch_v) = (surface_domain.u_interval(), surface_domain.v_interval());
        let mut curve_interval = span.interval;
        if surface_domain.are_disjoint(self.domain, 0.0) {
            return Ok(false);
        }
        if self.interval.are_disjoint(curve_interval, 0.0) {
            return Ok(false);
        }
        let mut expanded = span_box;
        expanded.expand_absolute(tolerance);
        if expanded.are_disjoint(&patch.pseudo_box) {
            return Ok(false);
        }
        let curve_size = curve_node.bbox.max - curve_node.bbox.min;
        let surface_size = surface_node.bbox.max - surface_node.bbox.min;
        let scaled_zero = EFF_ZERO * (1.0 + curve_node.bbox.max.max_dimension());
        let curve_bigger = curve_size.length() > surface_size.length() * 2.0;
        let mut t_conts = span.continuities;
        let mut u_conts = patch.u_continuities;
        let mut v_conts = patch.v_continuities;
        let domain = self.domain;
        if domain.min.x >= surface_domain.min.x {
            u_conts[0] = Continuity::Discontinuous;
        }
        if domain.max.x <= surface_domain.max.x {
            u_conts[1] = Continuity::Discontinuous;
        }
        if domain.min.y >= surface_domain.min.y {
            v_conts[0] = Continuity::Discontinuous;
        }
        if domain.max.y <= surface_domain.max.y {
            v_conts[1] = Continuity::Discontinuous;
        }
        if !surface_domain.is_contained_by(domain, EFF_ZERO) {
            if let Ok(meet) = surface_domain.intersect(domain) {
                surface_domain = meet;
            }
            if are_same(surface_domain.max.x, surface_domain.min.x)
                || are_same(surface_domain.max.y, surface_domain.min.y)
            {
                return Ok(false);
            }
        }
        if self.interval.min >= curve_interval.min {
            t_conts[0] = Continuity::Discontinuous;
        }
        if self.interval.max <= curve_interval.max {
            t_conts[1] = Continuity::Discontinuous;
        }
        if !curve_interval.is_contained_by(self.interval, EFF_ZERO) {
            curve_interval = curve_interval.intersect(self.interval)?;
            if are_same(curve_interval.max, curve_interval.min) {
                return Ok(false);
            }
        }
        let natural_t = curve.natural_interval();
        let natural = surface.natural_domain();
        let (natural_u, natural_v) = (natural.u_interval(), natural.v_interval());
        if (is_zero(curve_interval.min - natural_t.min)
            || is_zero(curve_interval.max - natural_t.max))
            && (is_zero(self.interval.min - natural_t.min)
                && is_zero(self.interval.max - natural_t.max))
            && is_closed(curve, natural_t, 0.0)
        {
            let class = evaluate_continuity(
                CurveSide {
                    curve,
                    parameter: natural_t.min,
                    from_left: false,
                },
                CurveSide {
                    curve,
                    parameter: natural_t.max,
                    from_left: true,
                },
                ANGLE_TOL_DEG,
            )
            .map_or(Continuity::Discontinuous, |(class, _, _)| class);
            if is_zero(curve_interval.min - natural_t.min) {
                t_conts[0] = class;
            }
            if is_zero(curve_interval.max - natural_t.max) {
                t_conts[1] = class;
            }
        }
        for (axis, conts, side, natural_side, low, high) in [
            (
                0,
                &mut u_conts,
                patch_u,
                natural_u,
                domain.min.x,
                domain.max.x,
            ),
            (
                1,
                &mut v_conts,
                patch_v,
                natural_v,
                domain.min.y,
                domain.max.y,
            ),
        ] {
            if (is_zero(side.min - natural_side.min) || is_zero(side.max - natural_side.max))
                && (is_zero(low - natural_side.min) && is_zero(high - natural_side.max))
            {
                let (closed, class) = is_closed_with_continuity(surface, natural, axis);
                if closed {
                    if is_zero(side.min - natural_side.min) {
                        conts[0] = class;
                    }
                    if is_zero(side.max - natural_side.max) {
                        conts[1] = class;
                    }
                }
            }
        }
        let center = surface_domain.evaluate(0.5, 0.5);
        let mut uv_guesses = vec![center];
        let mut t_guesses = vec![curve_interval.evaluate(0.5)];
        if curve_bigger {
            let drop = evaluate_point(surface, uv_guesses[0])?;
            if let Some(found) = crate::curve_solve::local_point_solve(
                curve,
                curve_interval,
                SolverOperation::Minimize,
                drop,
                t_guesses[0],
            )? {
                t_guesses[0] = found.parameter;
            }
        }
        let drop = curve.evaluate_point(t_guesses[0])?;
        uv_guesses[0] = drop_point_to_tangent_plane(surface, uv_guesses[0], [true, true], drop)?;
        if !surface_domain.contains_point(uv_guesses[0], 0.0) {
            if curve_bigger {
                uv_guesses[0] = surface_domain.clamp_point(uv_guesses[0]);
                let drop = evaluate_point(surface, uv_guesses[0])?;
                let mut second_t = curve_interval.evaluate(0.5);
                if let Some(found) = crate::curve_solve::local_point_solve(
                    curve,
                    curve_interval,
                    SolverOperation::Minimize,
                    drop,
                    second_t,
                )? {
                    second_t = found.parameter;
                    let drop = curve.evaluate_point(second_t)?;
                    let second_uv =
                        drop_point_to_tangent_plane(surface, uv_guesses[0], [true, true], drop)?;
                    uv_guesses.push(surface_domain.clamp_point(second_uv));
                    t_guesses.push(second_t);
                }
            } else {
                let mut placed = false;
                for t in [curve_interval.min, curve_interval.max] {
                    t_guesses[0] = t;
                    let drop = curve.evaluate_point(t)?;
                    uv_guesses[0] =
                        drop_point_to_tangent_plane(surface, center, [true, true], drop)?;
                    if surface_domain.contains_point(uv_guesses[0], 0.0) {
                        placed = true;
                        break;
                    }
                }
                if !placed {
                    t_guesses[0] = curve_interval.evaluate(0.5);
                    let drop = curve.evaluate_point(t_guesses[0])?;
                    uv_guesses[0] =
                        drop_point_to_tangent_plane(surface, center, [true, true], drop)?;
                    uv_guesses[0] = surface_domain.clamp_point(uv_guesses[0]);
                    for fraction in [0.5, 0.0, 1.0] {
                        uv_guesses.push(center);
                        t_guesses.push(curve_interval.evaluate(fraction));
                    }
                }
            }
        }
        if self.coincidence {
            uv_guesses.clear();
            t_guesses.clear();
        }
        let mut found_intersection = false;
        let mut deviation = 0.0;
        let mut i = 0;
        while i < uv_guesses.len() {
            let found = local_curve_intersect(
                surface,
                surface_domain,
                curve,
                curve_interval,
                tolerance,
                (t_guesses[i], uv_guesses[i]),
                self.zone,
            )?;
            found_intersection = found.is_some();
            if let Some(point) = found {
                deviation = point.distance;
                let solution = point_solution(branch, point.t, point.uv, point.distance);
                let keep = filter_solution(
                    surface,
                    surface_domain,
                    curve,
                    curve_interval,
                    u_conts,
                    v_conts,
                    t_conts,
                    Self::candidate(&solution),
                )?;
                if keep {
                    solver.add_sorted_solution(&*self, &solution, SortKey::ByAllParameters, true);
                    if self.are_perpendicular_disjoint(solver, branch)? {
                        return Ok(false);
                    }
                    if i < 1
                        && let Ok((_, Some((t, uv)))) = classify_meeting(
                            surface,
                            surface_domain,
                            curve,
                            curve_interval,
                            point.uv,
                            point.t,
                        )
                    {
                        uv_guesses.push(uv);
                        t_guesses.push(t);
                    }
                }
            }
            i += 1;
        }
        let pseudo = expanded;
        for j in 0..2 {
            if t_conts[j] <= Continuity::C0 && !self.coincidence {
                let t = if j == 0 {
                    curve_interval.min
                } else {
                    curve_interval.max
                };
                let point = curve.evaluate_point(t)?;
                if pseudo.distance_to_point(point) < EFF_ZERO {
                    let found = crate::surface_solve::local_point_solve(
                        surface,
                        natural,
                        SolverOperation::Minimize,
                        point,
                        surface_domain.evaluate(0.5, 0.5),
                        self.owner,
                    )?;
                    // SMLib solves over the natural domain and can keep a
                    // point outside the requested one; the port drops it.
                    if found.found
                        && found.value < tolerance
                        && self.domain.contains_point(found.uv, 0.0)
                    {
                        let solution = point_solution(branch, t, found.uv, found.value);
                        let keep = filter_solution(
                            surface,
                            natural,
                            curve,
                            curve_interval,
                            u_conts,
                            v_conts,
                            t_conts,
                            Self::candidate(&solution),
                        )?;
                        if keep {
                            solver.add_sorted_solution(
                                &*self,
                                &solution,
                                SortKey::ByAllParameters,
                                true,
                            );
                        }
                    }
                }
            }
            for axis in 0..2 {
                let conts = if axis == 0 { u_conts } else { v_conts };
                if conts[j] > Continuity::C0 {
                    continue;
                }
                let (fixed, from_left) = if axis == 0 {
                    if j == 0 {
                        (surface_domain.min.x, true)
                    } else {
                        (surface_domain.max.x, false)
                    }
                } else if j == 0 {
                    (surface_domain.min.y, true)
                } else {
                    (surface_domain.max.y, false)
                };
                let iso = IsoCurve::new(surface, axis, fixed, from_left);
                if crate::curve_solve::PointSolveCurve::is_degenerate(&iso, EFF_ZERO) {
                    continue;
                }
                let iso_interval = if axis == 0 {
                    surface_domain.v_interval()
                } else {
                    surface_domain.u_interval()
                };
                let coincidence = simple_coincidence_check(
                    curve,
                    curve_interval,
                    &iso,
                    iso_interval,
                    tolerance,
                    Some(&span_box),
                    Some(&patch.pseudo_box),
                )?;
                let as_solution = |pair: &CurvePair| {
                    let (u, v) = if axis == 0 {
                        (fixed, pair.parameters[1])
                    } else {
                        (pair.parameters[1], fixed)
                    };
                    point_solution(branch, pair.parameters[0], Vector2::new(u, v), pair.value)
                };
                let found: Vec<Solution> = match coincidence {
                    Coincidence::None => Vec::new(),
                    Coincidence::Point(p) => vec![as_solution(&p)],
                    Coincidence::Points([p, q]) => vec![as_solution(&p), as_solution(&q)],
                    Coincidence::Range([p, q]) => {
                        let mut range = as_solution(&p);
                        range.kind = SolutionType::RangeOfValues;
                        range.end = as_solution(&q).start;
                        vec![range]
                    }
                };
                let mut max_error: f64 = 0.0;
                let mut found_range = false;
                for solution in &found {
                    let skip = if axis == 0 {
                        found_intersection && solution.start.value + EFF_ZERO > deviation
                    } else {
                        found_intersection && solution.start.value > deviation + EFF_ZERO
                    };
                    if skip {
                        continue;
                    }
                    max_error = crate::scalar::max(max_error, solution.start.value);
                    if solution.kind == SolutionType::RangeOfValues {
                        found_range = true;
                        solver.add_sorted_solution(
                            &*self,
                            solution,
                            SortKey::ByAllParameters,
                            true,
                        );
                    } else if filter_solution(
                        surface,
                        surface_domain,
                        curve,
                        curve_interval,
                        u_conts,
                        v_conts,
                        t_conts,
                        Self::candidate(solution),
                    )? {
                        solver.add_sorted_solution(
                            &*self,
                            solution,
                            SortKey::ByAllParameters,
                            true,
                        );
                    }
                }
                if found_range {
                    return Ok(false);
                }
                let limit = if axis == 0 { EFF_ZERO } else { scaled_zero };
                if max_error > limit || found.is_empty() {
                    let meet = crate::curve_solve::local_curve_intersect(
                        &iso,
                        iso_interval,
                        curve,
                        curve_interval,
                        tolerance,
                        [iso_interval.evaluate(0.5), curve_interval.evaluate(0.5)],
                    )?;
                    if let Some(meet) = meet {
                        let (u, v) = if axis == 0 {
                            (fixed, meet.parameters[0])
                        } else {
                            (meet.parameters[0], fixed)
                        };
                        let solution = point_solution(
                            branch,
                            meet.parameters[1],
                            Vector2::new(u, v),
                            meet.value,
                        );
                        if filter_solution(
                            surface,
                            surface_domain,
                            curve,
                            curve_interval,
                            u_conts,
                            v_conts,
                            t_conts,
                            Self::candidate(&solution),
                        )? {
                            solver.add_sorted_solution(
                                &*self,
                                &solution,
                                SortKey::ByAllParameters,
                                true,
                            );
                        }
                    }
                }
            }
        }
        Ok(false)
    }

    fn identical_solutions(
        &self,
        objects: &dyn SolutionObjects,
        first: &Solution,
        second: &Solution,
    ) -> (SolutionPair, Preferred) {
        let (pair, preferred) = first.compare_for_identity(objects, second);
        if pair == SolutionPair::Identical {
            return (pair, preferred);
        }
        if pair == SolutionPair::DifferentParameters
            && first.kind == SolutionType::SingleValue
            && second.kind == SolutionType::SingleValue
        {
            let uv1 = Vector2::new(first.start.parameters[1], first.start.parameters[2]);
            let uv2 = Vector2::new(second.start.parameters[1], second.start.parameters[2]);
            let p1 = evaluate_point(self.surface, uv1).unwrap_or_default();
            let p2 = evaluate_point(self.surface, uv2).unwrap_or_default();
            let tolerance = self.tolerance;
            if !p1.closer_than(tolerance, p2) {
                return (pair, preferred);
            }
            // A failed evaluation reads as singular, as SMLib's status does.
            let singular = |uv: Vector2| {
                is_singularity(self.surface, uv, EFF_ZERO, self.zone, false)
                    .map_or(true, |d| d.is_some())
            };
            if !singular(uv1) || !singular(uv2) {
                return (pair, preferred);
            }
            let on1 = is_on_seam(self.surface, uv1, self.zone).1;
            let on2 = is_on_seam(self.surface, uv2, self.zone).1;
            match (on1, on2) {
                (true, false) => return (SolutionPair::Identical, Preferred::First),
                (false, true) => return (SolutionPair::Identical, Preferred::Second),
                (false, false) => {
                    let pick = if first.start.value < second.start.value {
                        Preferred::First
                    } else {
                        Preferred::Second
                    };
                    return (SolutionPair::Identical, pick);
                }
                (true, true) => {}
            }
        }
        (pair, preferred)
    }
}

/// `SmCurve::IsDegenerate` over `interval`: SMLib's B-spline test over the
/// whole curve, or its generic sampled test over a part.
pub(crate) fn degenerate_on(curve: &BSplineCurve, interval: Extent1) -> bool {
    degenerate_on_to(curve, interval, EFF_ZERO)
}

/// [`degenerate_on`] with SMLib's `d3DTol` argument.
pub(crate) fn degenerate_on_to(curve: &BSplineCurve, interval: Extent1, tolerance: f64) -> bool {
    if interval.are_equal(curve.natural_interval(), 0.0) {
        return curve.is_degenerate(tolerance);
    }
    if interval.min == interval.max {
        return true;
    }
    let p1 = curve.evaluate_point(interval.min).unwrap_or_default();
    let p2 = curve.evaluate_point(interval.max).unwrap_or_default();
    let scaled_zero = EFF_ZERO * (1.0 + p1.max_dimension());
    let tolerance2 = crate::scalar::max(tolerance * tolerance, scaled_zero * scaled_zero);
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
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        let degenerate = !(p1.distance_squared(values[0]) > tolerance2);
        return degenerate;
    }
    false
}

/// `sm_GlobalCurveIntersectCombineSolutions` for single-value answers:
/// merges answers that coincide or lie within the tolerance along the
/// curve, keeping the closer. Range answers and coincident meetings need
/// the surface point drop (`SmSurface::DropPoint`) and refuse.
#[allow(clippy::too_many_lines)]
fn combine_solutions(
    surface: &BSplineSurface,
    domain: Extent2,
    curve: &BSplineCurve,
    interval: Extent1,
    tolerance: f64,
    solutions: &mut SolutionArray,
) -> Result<(), SolverRefusal> {
    let mut items: Vec<Solution> = solutions.as_slice().to_vec();
    if items.iter().any(|s| s.kind == SolutionType::RangeOfValues) {
        return Err(SolverRefusal::NotPorted(INTERSECT));
    }
    let uv_size = surface.natural_domain().size();
    let t_size = curve.natural_interval().length();
    let mut coincidence = false;
    let classify_first = |items: &[Solution]| {
        let first = items[0].start.parameters;
        classify_meeting(
            surface,
            domain,
            curve,
            interval,
            Vector2::new(first[1], first[2]),
            first[0],
        )
        .is_ok_and(|(kind, _)| kind == MeetingType::Coincident)
    };
    if !items.is_empty() && classify_first(&items) {
        coincidence = true;
    }
    let mut last = 0usize;
    for ii in 1..items.len() {
        let last_end = items[last].start;
        let mut closest = ii;
        let scaled_zero = EFF_ZERO * (1.0 + items[ii].start.parameters[0].abs());
        for jj in ii + 1..items.len() {
            if (items[closest].start.parameters[0] - items[jj].start.parameters[0]).abs()
                < scaled_zero
            {
                let d1 = (items[closest].start.parameters[1] - last_end.parameters[1]).abs();
                let e1 = (items[closest].start.parameters[2] - last_end.parameters[2]).abs();
                let d2 = (items[jj].start.parameters[1] - last_end.parameters[1]).abs();
                let e2 = (items[jj].start.parameters[2] - last_end.parameters[2]).abs();
                if d2 + e2 < d1 + e1 {
                    closest = jj;
                }
            } else {
                break;
            }
        }
        if closest != ii {
            items.swap(closest, ii);
        }
        let current = items[ii];
        let c = current.start.parameters;
        let l = last_end.parameters;
        let uv_diff = Vector2::new((c[1] - l[1]).abs(), (c[2] - l[2]).abs());
        let t_diff = (c[0] - l[0]).abs();
        let far = uv_diff.x > uv_size.x * 0.90 || uv_diff.y > uv_size.y * 0.90;
        let combine = if far || t_diff > t_size * 0.90 {
            false
        } else if are_same(l[0], c[0]) || (are_same(l[1], c[1]) && are_same(l[2], c[2])) {
            true
        } else {
            let p1 = curve.evaluate_point(c[0]).unwrap_or_default();
            let p2 = curve.evaluate_point(l[0]).unwrap_or_default();
            if p1.distance_squared(p2) < tolerance * tolerance {
                true
            } else {
                if classify_first(&items) {
                    coincidence = true;
                }
                false
            }
        };
        if combine {
            if current.start.value < last_end.value {
                items[last] = current;
            }
        } else {
            last += 1;
            if last != ii {
                items[last] = current;
            }
        }
    }
    if !items.is_empty() {
        items.truncate(last + 1);
    }
    if coincidence {
        return Err(SolverRefusal::NotPorted(INTERSECT));
    }
    let mut combined = SolutionArray::new();
    for item in items {
        combined.add(item);
    }
    *solutions = combined;
    Ok(())
}

/// The straight-curve branch of `SmSurface::GlobalCurveIntersect`: the
/// line through the curve's point and exact tangent `pv` at the interval's
/// start, met by `GlobalLineIntersect` (analytic for a plane, else the line
/// search), its answers naming the curve and the surface.
pub(crate) fn line_intersect(
    surface: &BSplineSurface,
    domain: Extent2,
    pv: Vec<Vector3>,
    interval: Extent1,
    tolerance: f64,
) -> Result<SolutionArray, SolverRefusal> {
    let point = pv[0] - interval.min * pv[1];
    let mut solutions = if surface.plane_frame().is_some() {
        crate::plane::global_line_intersect(
            surface,
            domain,
            point,
            pv[1],
            Some(interval),
            tolerance,
        )?
    } else {
        crate::line_surface::global_line_intersect(
            surface,
            domain,
            point,
            pv[1],
            Some(interval),
            tolerance,
        )?
    };
    for index in 0..solutions.len() {
        let solution = &mut solutions[index];
        solution.object_count = 2;
        solution.objects[0] = Some(CURVE_OBJECT);
        solution.objects[1] = Some(SURFACE_OBJECT);
    }
    Ok(solutions)
}

/// `SmSurface::GlobalCurveIntersect`: every meeting of `curve` over
/// `interval` with `surface` over `domain` within `tolerance`, sorted along
/// the curve, as `SmSolution`s with parameters `t, u, v`. A degenerate curve
/// is solved as a point.
pub fn global_curve_intersect(
    surface: &BSplineSurface,
    domain: Extent2,
    curve: &BSplineCurve,
    interval: Extent1,
    tolerance: f64,
    owner_tolerance: Option<f64>,
) -> Result<SolutionArray, SolverRefusal> {
    let zone = owner_tolerance.unwrap_or(ZONE_TOL_3D);
    if degenerate_on(curve, interval) {
        let point = curve.evaluate_point(interval.min)?;
        // SMLib ignores the search status here and keeps whatever it found;
        // the port reports a failed search as a typed failure instead.
        let found = crate::surface_global::global_point_solve(
            surface,
            domain,
            INTERSECT,
            point,
            tolerance,
            owner_tolerance,
        )?;
        let mut solutions = SolutionArray::new();
        for mut solution in found.as_slice().iter().copied() {
            // SMLib leaves the surface solver's two variables; the answer
            // carries `t, u, v`, so the port declares three.
            solution.variables = 3;
            solution.object_count = 2;
            solution.objects[0] = Some(CURVE_OBJECT);
            solution.objects[1] = Some(SURFACE_OBJECT);
            solution.start.parameters[2] = solution.start.parameters[1];
            solution.start.parameters[1] = solution.start.parameters[0];
            solution.start.parameters[0] = interval.min;
            solutions.add(solution);
        }
        return Ok(solutions);
    }
    if curve.is_line(10, EFF_ZERO)? {
        let pv = curve.evaluate(interval.min, 1, true, false)?;
        return line_intersect(surface, domain, pv, interval, tolerance);
    }
    let coincident = crate::curve_surface::simple_coincidence_check(
        surface,
        domain,
        curve,
        interval,
        tolerance,
        true,
        owner_tolerance,
    )?;
    if coincident.is_some() {
        return Err(SolverRefusal::NotPorted(INTERSECT));
    }
    let surface_cache =
        SurfaceCache::for_trimmed_surface(surface).map_err(|refusal| match refusal {
            CacheRefusal::Failed(fault) => SolverRefusal::Failed(fault),
            CacheRefusal::NotPorted(_) => SolverRefusal::NotPorted(INTERSECT),
        })?;
    let curve_cache = CurveCache::for_curve(curve, 0.5 * zone)?;
    let mut trees = SearchTrees::new();
    trees.append(&curve_cache.trees);
    let offset = trees.append(&surface_cache.trees);
    let objects = CurveSurfaceObjects {
        curve,
        surface,
        zone,
    };
    let mut solver = GlobalSolver::with_two_trees(
        &trees,
        &objects,
        Some(curve_cache.top),
        Some(NodeId(surface_cache.top.0 + offset)),
        3,
    );
    let mut hooks = CurveSurfaceHooks {
        curve,
        surface,
        domain,
        interval,
        coincidence: false,
        curve_cache: &curve_cache,
        surface_cache: &surface_cache,
        surface_offset: offset,
        zone,
        owner: owner_tolerance,
        tolerance,
    };
    let mut solutions = SolutionArray::new();
    solver.solve_it(
        &mut hooks,
        INTERSECT,
        SolutionRequested::All,
        tolerance,
        BIG_DOUBLE,
        &[],
        &mut solutions,
    )?;
    combine_solutions(surface, domain, curve, interval, tolerance, &mut solutions)?;
    Ok(solutions)
}
