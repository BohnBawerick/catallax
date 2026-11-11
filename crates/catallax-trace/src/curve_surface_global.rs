// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline curve and surface
// candidates of `spline.rs` without trimming faces; statuses become `Fault`,
// tree-node pointers become identities in one arena holding the curve and
// surface trees, and solution containers keep SMLib's `SmSolution` layout.
// Paths that need code from later port steps refuse with
// `SolverRefusal::NotPorted`: a coincident curve that needs SMLib's general
// curve drop, and the polar-box cone test. A coincident curve's ranges come
// from `DropAndTrimCurve`. The solution combiner keeps three SMLib slips
// that can change answers: the join test of overlapping ranges reads the
// second parameter twice, the point-beside-range test leaves its third span
// empty (so it skips every point whose tree nodes carry data), and a later
// answer's cull mark is cleared when the loop reaches it.
// Three answers deliberately differ from SMLib: a failed point search for a
// degenerate curve is a typed failure (SMLib ignores its status), its
// answers declare the three variables they carry (SMLib leaves two), and the
// end-point fallback drops points outside the requested domain (SMLib keeps
// them).
// Translated from SMLib SmSurface.cpp (SmSurface::GlobalCurveIntersect,
// SmCSIGlobalSolver, sm_CsArePerpenDisjoint,
// sm_GlobalCurveIntersectCombineSolutions,
// SmSurface::ClassifyIntersectingRangeSolution); rust-brep has no port.
//! The global curve-surface intersector of step 6: the curve and surface
//! cache trees searched together, local intersections in each pair of
//! leaves, checks along discontinuous ends and sides, and the combination
//! of the answers.

use crate::Fault;
use crate::curve::{CurveSide, approximate_length, evaluate_continuity};
use crate::curve_cache::{CurveCache, NodeData};
use crate::curve_solve::{Coincidence, CurvePair, is_closed, simple_coincidence_check};
use crate::curve_surface::{
    Candidate, MeetingType, classify_meeting, drop_point_to_tangent_plane, filter_solution,
    is_singularity, local_curve_intersect,
};
use crate::extent::{Extent1, Extent2};
use crate::geometry::Continuity;
use crate::global_solver::{
    AuxData, Branch, GlobalSolver, SearchTrees, SolutionRequested, SolverHooks, SolverOperation,
    SolverRefusal,
};
use crate::scalar::{
    ANGLE_TOL_DEG, BIG_DOUBLE, EFF_ZERO, PI, ZONE_TOL_3D, are_same, is_between, is_zero,
};
use crate::solution::{
    MAX_PARAMETERS, NodeId, ObjectId, ObjectKind, Preferred, Solution, SolutionArray, SolutionEnd,
    SolutionObjects, SolutionPair, SolutionType, SortKey,
};
use crate::spline::{BSplineCurve, BSplineSurface};
use crate::surface_cache::{CacheRefusal, SurfaceCache, SurfaceNodeData};
use crate::surface_drop::drop_point;
use crate::surface_drop_curve::{curve_ranges, drop_and_trim_curve};
use crate::surface_eval::{evaluate_normal, evaluate_point};
use crate::surface_solve::{IsoCurve, is_closed_with_continuity, is_on_seam, local_point_solve};
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

/// `SmSurface::ClassifyIntersectingRangeSolution` for a range answer:
/// whether the curve only touches the surface there. From the end with the
/// smaller deviation, points at geometrically growing fractions towards the
/// other end are dropped normally to the surface; a gap (beyond the one at
/// that end) that grows monotonically and diverges marks a touching point.
/// Returns that end, as SMLib's replacement start, when it does.
fn classify_intersecting_range(
    surface: &BSplineSurface,
    curve: &BSplineCurve,
    range: &Solution,
) -> Result<Option<SolutionEnd>, SolverRefusal> {
    let interval = curve.natural_interval();
    let region = surface.natural_domain();
    let (near, far) = if range.start.value < range.end.value {
        (range.start, range.end)
    } else {
        (range.end, range.start)
    };
    let start = near.parameters;
    let delta = [
        far.parameters[0] - start[0],
        far.parameters[1] - start[1],
        far.parameters[2] - start[2],
    ];
    let curve_point = curve.evaluate_point(start[0]).unwrap_or_default();
    let surface_point =
        evaluate_point(surface, Vector2::new(start[1], start[2])).unwrap_or_default();
    let crossing_gap = surface_point - curve_point;
    let scaled_zero = EFF_ZERO * (1.0 + curve_point.max_dimension());
    // The gap at the near end is exactly zero once the crossing gap is
    // taken off.
    let mut gaps = [Vector3::default(); 3];
    let mut lengths = [0.0f64; 3];
    let (mut zero, mut monotonic, mut diverging) = (true, true, true);
    let mut degenerate = false;
    let mut done = false;
    let mut fraction = 1.0 / 16.0001;
    let mut ii = 1usize;
    while !done && fraction < 1.0 {
        let target = if ii > 1 { 2 } else { 1 };
        let t = start[0] + fraction * delta[0];
        let guess = Vector2::new(
            start[1] + fraction * delta[1],
            start[2] + fraction * delta[2],
        );
        let point = curve.evaluate_point(t).unwrap_or_default();
        let found = match local_point_solve(
            surface,
            region,
            SolverOperation::Normalize,
            point,
            guess,
            None,
        ) {
            Ok(found) => found.found.then_some(found.uv),
            Err(SolverRefusal::Failed(_)) => None,
            Err(refusal) => return Err(refusal),
        };
        let Some(uv) = found else {
            degenerate = true;
            done = true;
            ii += 1;
            fraction *= 2.0;
            continue;
        };
        let on_surface = evaluate_point(surface, uv).unwrap_or_default();
        gaps[target] = (on_surface - point) - crossing_gap;
        lengths[target] = gaps[target].length();
        if ii >= 2 {
            let previous = (target + 2) % 3;
            let before = (target + 1) % 3;
            if lengths[target] > scaled_zero {
                zero = false;
            }
            if gaps[target].dot(gaps[previous]) <= 0.0 {
                monotonic = false;
            }
            if lengths[target] == 0.0
                || lengths[target]
                    < lengths[previous] + 1.0001 * (lengths[previous] - lengths[before])
            {
                diverging = false;
            }
            if !zero && monotonic && diverging {
                degenerate = true;
                done = true;
            } else if zero && monotonic && diverging {
                gaps[1] = gaps[2];
                lengths[1] = lengths[2];
            } else {
                degenerate = false;
                done = true;
            }
        }
        ii += 1;
        fraction *= 2.0;
    }
    if !done && diverging {
        let region_size = region.size().length();
        degenerate = delta[0].abs() < interval.length() / 10.0
            && delta[1].abs() < region_size / 10.0
            && delta[2].abs() < region_size / 10.0;
    }
    if !degenerate {
        return Ok(None);
    }
    let mut point = range.start;
    point.value = near.value;
    point.parameters[..3].copy_from_slice(&start[..3]);
    Ok(Some(point))
}

/// Whether the curve point at `fraction` of `gap` drops to the surface
/// within `tolerance`, from the matching fraction of the surface gaps.
fn gap_point_drops(
    surface: &BSplineSurface,
    domain: Extent2,
    curve: &BSplineCurve,
    tolerance: f64,
    gap: [Extent1; 3],
    fraction: f64,
) -> Result<bool, SolverRefusal> {
    let point = curve.evaluate_point(gap[0].evaluate(fraction))?;
    let guess = Vector2::new(gap[1].evaluate(fraction), gap[2].evaluate(fraction));
    let drop = drop_point(
        surface,
        domain,
        point,
        Some(guess),
        SolverOperation::Minimize,
    )?;
    Ok(drop.success && drop.gap <= tolerance)
}

/// Whether two consecutive range answers join: SMLib's test of the gap
/// between the end of `last` and the start of `current`, along the curve
/// and over the surface.
fn ranges_join(
    surface: &BSplineSurface,
    domain: Extent2,
    curve: &BSplineCurve,
    tolerance: f64,
    last: &Solution,
    current: &Solution,
) -> Result<bool, SolverRefusal> {
    let l = last.end.parameters;
    let c = current.start.parameters;
    if l[0] >= c[0] {
        // SMLib tests the second parameter twice and never the third.
        return Ok(is_between(
            c[1],
            last.start.parameters[1],
            last.end.parameters[1],
        ));
    }
    let gap = [
        Extent1::new(l[0], c[0]),
        Extent1::new(l[1], c[1]),
        Extent1::new(l[2], c[2]),
    ];
    for jj in 0..2 {
        if !gap_point_drops(
            surface,
            domain,
            curve,
            tolerance,
            gap,
            0.33 * f64::from(jj + 1),
        )? {
            return Ok(false);
        }
    }
    let uv_size = surface.natural_domain().size();
    let t_size = curve.natural_interval().length();
    if gap[0].length() < t_size / 1000.0
        && (gap[1].length() < uv_size.x / 1000.0 || gap[2].length() < uv_size.y / 1000.0)
    {
        return Ok(true);
    }
    for jj in 0..3 {
        if !gap_point_drops(
            surface,
            domain,
            curve,
            tolerance,
            gap,
            0.25 * f64::from(jj + 1),
        )? {
            return Ok(false);
        }
    }
    Ok(true)
}

/// `sm_GlobalCurveIntersectCombineSolutions`. Point answers and range
/// answers are merged separately along the curve: points that coincide or
/// lie within the tolerance keep the closer, ranges that touch or whose gap
/// drops to the surface join. Then ranges too short to classify, or whose
/// classification finds a touching point, become points. When the first
/// point is a coincident meeting, or a range covers more than 2% of the
/// interval, the curve's `DropAndTrimCurve` pieces replace every range and
/// every point at which the curve is tangent to the surface, unless the
/// curve lies farther than the tolerance from the surface.
#[allow(clippy::too_many_lines, clippy::too_many_arguments)]
fn combine_solutions(
    surface: &BSplineSurface,
    domain: Extent2,
    curve: &BSplineCurve,
    interval: Extent1,
    tolerance: f64,
    trees: &SearchTrees,
    solutions: &mut SolutionArray,
) -> Result<(), SolverRefusal> {
    const NONE: usize = usize::MAX;
    let mut items: Vec<Solution> = solutions.as_slice().to_vec();
    let uv_size = surface.natural_domain().size();
    let t_size = curve.natural_interval().length();
    let first_is_coincident = |items: &[Solution]| {
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
    let mut coincidence = false;
    // The last kept range and the last kept point; either may be a range
    // after a point and a range swap places below.
    let (mut last_range, mut last_point) = (NONE, NONE);
    if let Some(first) = items.first() {
        if first.kind == SolutionType::RangeOfValues {
            last_range = 0;
        } else {
            last_point = 0;
            coincidence = first_is_coincident(&items);
        }
    }
    for ii in 1..items.len() {
        let range_slot = items[ii].kind == SolutionType::RangeOfValues;
        let last_index = if range_slot { last_range } else { last_point };
        let mut combine = false;
        if last_index != NONE {
            let last = items[last_index];
            let last_end = if last.kind == SolutionType::RangeOfValues {
                last.end
            } else {
                last.start
            };
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
            let is_range = current.kind == SolutionType::RangeOfValues;
            let far = (c[1] - l[1]).abs() > uv_size.x * 0.90
                || (c[2] - l[2]).abs() > uv_size.y * 0.90
                || (c[0] - l[0]).abs() > t_size * 0.90;
            if far {
                combine = false;
            } else if are_same(l[0], c[0]) || (are_same(l[1], c[1]) && are_same(l[2], c[2])) {
                combine = true;
                if is_range {
                    if last.start.value > EFF_ZERO && last.start.value > last.end.value * 100_000.0
                    {
                        combine = false;
                    }
                    if current.end.value > EFF_ZERO
                        && current.end.value > current.start.value * 100_000.0
                    {
                        combine = false;
                    }
                    let last_length = crate::scalar::max(
                        EFF_ZERO,
                        last.end.parameters[0] - last.start.parameters[0],
                    );
                    let current_length = crate::scalar::max(
                        EFF_ZERO,
                        current.end.parameters[0] - current.start.parameters[0],
                    );
                    if last_length > current_length * 1000.0
                        || current_length > last_length * 1000.0
                    {
                        combine = true;
                    }
                }
            } else if is_range {
                combine = ranges_join(surface, domain, curve, tolerance, &last, &current)?;
            } else {
                let p1 = curve.evaluate_point(c[0])?;
                let p2 = curve.evaluate_point(l[0])?;
                combine = p1.distance_squared(p2) < tolerance * tolerance;
                if !combine && first_is_coincident(&items) {
                    coincidence = true;
                }
            }
            if combine {
                if is_range {
                    if items[last_index].end.parameters[0] < current.end.parameters[0] {
                        items[last_index].end = current.end;
                    }
                } else if current.start.value < last_end.value {
                    items[last_index] = current;
                }
            }
        }
        if !combine {
            let next = if last_range == NONE {
                last_point + 1
            } else if last_point == NONE {
                last_range + 1
            } else {
                last_range.max(last_point) + 1
            };
            if range_slot {
                last_range = next;
            } else {
                last_point = next;
            }
            if next != ii {
                items[next] = items[ii];
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

    // Points next to a range: SMLib meant to drop a point inside a range, or
    // close to an end and much farther from the surface than that end, and a
    // range too short to keep beside a point. It puts the third parameters
    // into the second span and leaves the third empty, so every point whose
    // tree nodes carry data is skipped and the first other point stops the
    // search. Each answer's mark is cleared when the outer loop reaches it,
    // dropping marks an earlier range set on a later point.
    let mut cull = vec![false; items.len()];
    let mut compress = false;
    for ii in 0..items.len() {
        cull[ii] = false;
        if items[ii].kind != SolutionType::RangeOfValues {
            continue;
        }
        let range = items[ii];
        let span1 = Extent1::new(range.start.parameters[0], range.end.parameters[0]);
        let mut span2 = Extent1::new(range.start.parameters[1], range.end.parameters[1]);
        span2.add_value(range.start.parameters[2]);
        span2.add_value(range.end.parameters[2]);
        let span3 = Extent1::default();
        for jj in 0..items.len() {
            if ii == jj || items[jj].kind != SolutionType::SingleValue {
                continue;
            }
            let point = items[jj];
            let aux = |node: Option<NodeId>| node.and_then(|id| trees.get(id)).map(|n| n.aux);
            let curve_node = matches!(
                aux(point.nodes[0]),
                Some(AuxData::BezierCurve | AuxData::CurveData)
            );
            let surface_node = matches!(
                aux(point.nodes[1]),
                Some(AuxData::BezierSurface | AuxData::BezierSurfaceHead | AuxData::AuxData)
            );
            if curve_node && surface_node {
                continue;
            }
            let close = interval.length() / 1000.0;
            if span1.length() < close {
                cull[ii] = true;
                compress = true;
            } else {
                let t = point.start.parameters[0];
                if t < span1.min {
                    if span1.min - t < close && point.start.value > 2.0 * range.start.value {
                        cull[jj] = true;
                        compress = true;
                    }
                } else if t > span1.max {
                    if t - span1.max < close && point.start.value > 2.0 * range.end.value {
                        cull[jj] = true;
                        compress = true;
                    }
                } else {
                    cull[jj] = true;
                    compress = true;
                }
            }
            if cull[ii] {
                break;
            }
            let p = point.start.parameters;
            if p[0] > span1.max || p[1] > span2.max || p[2] > span3.max {
                break;
            }
        }
    }
    if compress {
        let mut keep = cull.iter().map(|culled| !culled);
        items.retain(|_| keep.next().unwrap_or(true));
    }

    // Ranges too short to classify, or that only touch the surface, become
    // points at their end nearer the surface. Like SMLib the touching flag
    // carries over a range too short to classify, which becomes a point
    // either way.
    let mut length = BIG_DOUBLE;
    let mut degenerate = false;
    for item in &mut items {
        if item.kind != SolutionType::RangeOfValues {
            continue;
        }
        let start_value = item.start.value;
        let end_value = item.end.value;
        let span = Extent1::new(item.start.parameters[0], item.end.parameters[0]);
        let region = Extent2::new(
            Vector2::new(item.start.parameters[1], item.start.parameters[2]),
            Vector2::new(item.end.parameters[1], item.end.parameters[2]),
        );
        let span_length = span.length();
        let region_length = region.size().length();
        let point = curve.evaluate_point(item.start.parameters[0])?;
        if length == BIG_DOUBLE {
            length = approximate_length(curve, interval, 5);
        }
        let piece = approximate_length(curve, span, 5);
        let scaled_zero = EFF_ZERO * (1.0 + point.max_dimension());
        let span_tolerance = 100.0 * EFF_ZERO * (1.0 + interval.max_dimension());
        let region_tolerance = 100.0 * EFF_ZERO * (1.0 + domain.max_dimension());
        if start_value < scaled_zero
            && end_value < scaled_zero
            && span_length > span_tolerance
            && region_length > region_tolerance
        {
            continue;
        }
        let large = piece > length / 10.0;
        let small = span_length < span_tolerance || region_length < region_tolerance;
        if !small {
            let touching = classify_intersecting_range(surface, curve, item)?;
            degenerate = touching.is_some();
            if let Some(start) = touching {
                item.start = start;
            }
        }
        if large && !degenerate {
            continue;
        }
        if small || degenerate {
            if item.start.value > item.end.value {
                item.start = item.end;
            }
            item.kind = SolutionType::SingleValue;
        }
    }

    if !coincidence {
        coincidence = items.iter().any(|item| {
            item.kind == SolutionType::RangeOfValues
                && item.end.parameters[0] - item.start.parameters[0] > 0.02 * interval.length()
        });
    }
    if coincidence {
        let mut kept = Vec::new();
        for item in &items {
            if item.kind != SolutionType::SingleValue {
                continue;
            }
            let p = item.start.parameters;
            let pv = curve.evaluate(p[0], 1, true, true)?;
            let normal = evaluate_normal(surface, Vector2::new(p[1], p[2]), [true, true]);
            if !pv[1].is_perpendicular_to(normal, 2.0) {
                kept.push(*item);
            }
        }
        let dropped = drop_and_trim_curve(surface, domain, curve, interval, tolerance, true, true)?;
        if dropped.max_drop <= tolerance && !dropped.curves.is_empty() {
            for mut range in curve_ranges(&dropped.curves, dropped.max_drop)? {
                range.object_count = 2;
                range.objects[0] = Some(CURVE_OBJECT);
                range.objects[1] = Some(SURFACE_OBJECT);
                kept.push(range);
            }
            items = kept;
        }
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
        // The curve's pieces over the surface, as ranges naming both.
        let ranges = crate::surface_drop_curve::coincident_ranges(
            surface, domain, curve, interval, tolerance, true, true,
        )?;
        let mut solutions = SolutionArray::new();
        for mut range in ranges {
            range.object_count = 2;
            range.objects[0] = Some(CURVE_OBJECT);
            range.objects[1] = Some(SURFACE_OBJECT);
            solutions.add(range);
        }
        return Ok(solutions);
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
    combine_solutions(
        surface,
        domain,
        curve,
        interval,
        tolerance,
        &trees,
        &mut solutions,
    )?;
    Ok(solutions)
}
