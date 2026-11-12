// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline surfaces of `spline.rs`
// without trimming faces, so every node and point is valid and the surface
// to evaluate is the cache's own; statuses become `Fault`, the solver's
// virtuals are `SolverHooks`. A ray refuses with
// `SolverRefusal::NotPorted`, and so does a line lying in a B-spline surface
// where SMLib's general curve drop is needed. Like the
// curve-surface port, a line end solved over the natural domain is dropped
// when it lies outside the requested domain. A box clip that cannot be
// computed refuses instead of pruning, where SMLib ignores its status.
// Translated from SMLib SmSurface.cpp (SmLSIGlobalSolver,
// BranchMayContainAnswers, LocalSolve, IdenticalSolutions,
// SmSurface::GlobalLineIntersect, SmSurfaceEval::Eval, EvalNormal,
// EvalNormalSection); rust-brep has no port.
//! The global line-surface intersection (`SmSurface::GlobalLineIntersect`)
//! that `GlobalCurveIntersect` uses for straight curves: a subdivision
//! search of the surface's patches along the line with Newton steps from
//! guesses, the line's meetings with the patches' C0 sides, and the merge of
//! duplicate and adjacent answers.
#![allow(clippy::needless_range_loop)]

use crate::Fault;
use crate::curve_solve::{Coincidence, CurvePair, PointSolveCurve};
use crate::curve_surface::{
    Candidate, SurfaceJet, classify_meeting, drop_point_to_tangent_plane, filter_solution,
    is_singularity, local_curve_intersect, normal_section_of,
};
use crate::extent::{Extent1, Extent2, Extent3};
use crate::geometry::{Continuity, line_closest_point};
use crate::global_solver::{
    Branch, GlobalSolver, SolutionRequested, SolverHooks, SolverOperation, SolverRefusal,
};
use crate::line::Line;
use crate::polar_box::PolarBoxes;
use crate::scalar::{BIG_DOUBLE, EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, are_same, max};
use crate::solution::{
    MAX_PARAMETERS, ObjectId, ObjectKind, Preferred, Solution, SolutionArray, SolutionEnd,
    SolutionObjects, SolutionPair, SolutionType, SortKey,
};
use crate::spline::BSplineSurface;
use crate::surface_cache::{CacheRefusal, SurfaceCache, SurfaceNodeData};
use crate::surface_eval::{evaluate, evaluate_normal, evaluate_point};
use crate::surface_solve::{IsoCurve, is_on_seam};
use crate::vector::{Vector2, Vector3};

const INTERSECT: SolverOperation = SolverOperation::Intersect;

/// The surface's identity in line-surface solutions (the line has none).
pub const SURFACE_OBJECT: ObjectId = ObjectId(1);

struct LineSurfaceObjects<'a> {
    surface: &'a BSplineSurface,
}

impl SolutionObjects for LineSurfaceObjects<'_> {
    fn kind(&self, id: ObjectId) -> ObjectKind {
        if id == SURFACE_OBJECT {
            ObjectKind::Surface
        } else {
            ObjectKind::Other
        }
    }
    fn parameter_interval(&self, id: ObjectId, direction: usize) -> Extent1 {
        match (id, direction) {
            (SURFACE_OBJECT, 0) => self.surface.natural_domain().u_interval(),
            (SURFACE_OBJECT, 1) => self.surface.natural_domain().v_interval(),
            _ => Extent1::default(),
        }
    }
    fn zone_tolerance(&self, id: ObjectId) -> f64 {
        if id == SURFACE_OBJECT {
            crate::scalar::ZONE_TOL_3D
        } else {
            0.0
        }
    }
}

fn single(branch: &Branch, t: f64, uv: Vector2, value: f64) -> Solution {
    let mut parameters = [0.0; MAX_PARAMETERS];
    parameters[..3].copy_from_slice(&[t, uv.x, uv.y]);
    let mut solution = Solution {
        kind: SolutionType::SingleValue,
        variables: 3,
        start: SolutionEnd { value, parameters },
        object_count: 2,
        ..Solution::default()
    };
    solution.objects[1] = Some(SURFACE_OBJECT);
    solution.nodes[0] = branch[0];
    solution.nodes[1] = branch[1];
    solution
}

fn candidate(solution: &Solution) -> Candidate {
    let p = solution.start.parameters;
    Candidate {
        t: p[0],
        uv: Vector2::new(p[1], p[2]),
        value: solution.start.value,
    }
}

/// `SmLSIGlobalSolver`.
struct LineSurfaceHooks<'a> {
    surface: &'a BSplineSurface,
    cache: &'a SurfaceCache,
    polar: PolarBoxes<'a>,
    domain: Extent2,
    point: Vector3,
    vector: Vector3,
    interval: Option<Extent1>,
    tolerance: f64,
    line: Line,
}

impl LineSurfaceHooks<'_> {
    fn keep(
        &self,
        solver: &mut GlobalSolver<'_>,
        solution: &Solution,
        domain: Extent2,
        line_interval: Extent1,
        conts: [[Continuity; 2]; 3],
    ) -> Result<bool, SolverRefusal> {
        let [u_conts, v_conts, t_conts] = conts;
        let keep = filter_solution(
            self.surface,
            domain,
            &self.line,
            line_interval,
            u_conts,
            v_conts,
            t_conts,
            candidate(solution),
        )?;
        if keep {
            solver.add_sorted_solution(self, solution, SortKey::ByAllParameters, true);
        }
        Ok(keep)
    }
}

impl SolverHooks for LineSurfaceHooks<'_> {
    /// `SmLSIGlobalSolver::BranchMayContainAnswers`: the line passes within
    /// the tolerance of the node box's bounding sphere (and of its interval
    /// ends), and through the box or near its inner sphere.
    fn branch_may_contain_answers(
        &mut self,
        solver: &GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<bool, SolverRefusal> {
        let _ = solver;
        let id = branch[0].ok_or(Fault::Degenerate)?;
        let bbox = self.cache.trees.get(id).ok_or(Fault::Degenerate)?.bbox;
        let (center, radius, min_radius) = bbox.sphere_bound();
        let reach = (radius + self.tolerance) * (radius + self.tolerance);
        let t = line_closest_point(self.point, self.vector, center);
        let on_line = self.point + t * self.vector;
        let distance_sq = on_line.distance_squared(center);
        if distance_sq > reach {
            return Ok(false);
        }
        if let Some(interval) = self.interval {
            for (beyond, end) in [
                (t < interval.min, interval.min),
                (t > interval.max, interval.max),
            ] {
                if beyond && (self.point + end * self.vector).distance_squared(center) > reach {
                    return Ok(false);
                }
            }
        }
        let inner = (min_radius + self.tolerance) * (min_radius + self.tolerance);
        if distance_sq <= inner || bbox.contains_point(on_line, self.tolerance) {
            return Ok(true);
        }
        // SMLib ignores the clip's status, so a failure finds nothing there;
        // here it refuses rather than prune a branch it could not settle.
        Ok(bbox.intersect_line(on_line, self.vector)?.is_some())
    }

    /// `SmLSIGlobalSolver::LocalSolve` on a leaf patch.
    #[allow(clippy::too_many_lines)]
    fn local_solve(
        &mut self,
        solver: &mut GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<bool, SolverRefusal> {
        let id = branch[0].ok_or(Fault::Degenerate)?;
        let SurfaceNodeData::Patch(index) = self.cache.node_data[id.0 as usize] else {
            return Err(Fault::Degenerate.into());
        };
        let patch = &self.cache.patches[index];
        let surface = self.surface;
        let mut uv_domain = patch.uv_domain;
        if uv_domain.are_disjoint(self.domain, 0.0) {
            return Ok(false);
        }
        let t_conts = [Continuity::Discontinuous; 2];
        let mut u_conts = patch.u_continuities;
        let mut v_conts = patch.v_continuities;
        if self.domain.min.x >= uv_domain.min.x {
            u_conts[0] = Continuity::Discontinuous;
        }
        if self.domain.max.x <= uv_domain.max.x {
            u_conts[1] = Continuity::Discontinuous;
        }
        if self.domain.min.y >= uv_domain.min.y {
            v_conts[0] = Continuity::Discontinuous;
        }
        if self.domain.max.y <= uv_domain.max.y {
            v_conts[1] = Continuity::Discontinuous;
        }
        if !uv_domain.is_contained_by(self.domain, EFF_ZERO) {
            if let Ok(meet) = uv_domain.intersect(self.domain) {
                uv_domain = meet;
            }
            if are_same(uv_domain.max.x, uv_domain.min.x)
                || are_same(uv_domain.max.y, uv_domain.min.y)
            {
                return Ok(false);
            }
        }
        let tolerance = solver.tolerance;
        let mut pseudo = patch.pseudo_box;
        pseudo.expand_absolute(self.tolerance);
        let mut grown = self.cache.trees.get(id).ok_or(Fault::Degenerate)?.bbox;
        grown.expand_absolute(tolerance);
        // A clip that cannot be computed refuses instead of pruning the leaf.
        let Some(box_span) = grown.intersect_line(self.point, self.vector)? else {
            return Ok(false);
        };
        let Some(pseudo_span) = pseudo.intersect_line(self.point, self.vector)? else {
            return Ok(false);
        };
        if max(box_span.0, pseudo_span.0) > crate::scalar::min(box_span.1, pseudo_span.1) {
            return Ok(false);
        }
        let mut bound = Extent1::new(
            max(box_span.0, pseudo_span.0),
            crate::scalar::min(box_span.1, pseudo_span.1),
        );
        bound.expand_absolute(2.0 * tolerance);
        if let Some(interval) = self.interval {
            if let Ok(meet) = bound.intersect(interval) {
                bound = meet;
            }
            if interval.are_disjoint(bound, 0.0) {
                return Ok(false);
            }
        }
        let span = match self.interval {
            Some(interval) => interval.intersect(bound).unwrap_or(bound),
            None => bound,
        };
        let center = uv_domain.evaluate(0.5, 0.5);
        let mut t_guesses = vec![span.evaluate(0.5)];
        let mut uv_guesses = vec![center];
        let at = |t: f64| self.point + t * self.vector;
        uv_guesses[0] =
            drop_point_to_tangent_plane(surface, center, [true, true], at(t_guesses[0]))?;
        let eps = EFF_ZERO * (1.0 + uv_domain.size().max_dimension());
        if !uv_domain.contains_point(uv_guesses[0], eps) {
            t_guesses[0] = span.evaluate(0.1);
            uv_guesses[0] =
                drop_point_to_tangent_plane(surface, center, [true, true], at(t_guesses[0]))?;
            if !uv_domain.contains_point(uv_guesses[0], 0.0) {
                t_guesses[0] = span.evaluate(0.9);
                uv_guesses[0] =
                    drop_point_to_tangent_plane(surface, center, [true, true], at(t_guesses[0]))?;
                if !uv_domain.contains_point(uv_guesses[0], 0.0) {
                    uv_guesses[0] = uv_domain.clamp_point(uv_guesses[0]);
                    for fraction in [0.5, 0.1, 0.9] {
                        uv_guesses.push(center);
                        t_guesses.push(span.evaluate(fraction));
                    }
                }
            }
        }
        let conts = [u_conts, v_conts, t_conts];
        let mut unit = self.vector;
        let mut i = 0;
        while i < uv_guesses.len() {
            let found = local_curve_intersect(
                surface,
                uv_domain,
                &self.line,
                span,
                tolerance,
                (t_guesses[i], uv_guesses[i]),
                crate::scalar::ZONE_TOL_3D,
            )?;
            if let Some(found) = found {
                let solution = single(branch, found.t, found.uv, found.distance);
                if !self.keep(solver, &solution, uv_domain, span, conts)? {
                    i += 1;
                    continue;
                }
                unit.unitize()?;
                if !self
                    .polar
                    .get(id.0 as usize)
                    .has_perpendicular_to_vector(unit)
                {
                    return Ok(false);
                }
                if i < 1
                    && let Ok((_, Some((t, uv)))) =
                        classify_meeting(surface, uv_domain, &self.line, span, found.uv, found.t)
                {
                    uv_guesses.push(uv);
                    t_guesses.push(t);
                }
            }
            i += 1;
        }
        let natural = surface.natural_domain();
        for j in 0..2 {
            // The line's own ends (its continuity is always broken there).
            let t = if j == 0 { span.min } else { span.max };
            let end = self.line.evaluate_point(t);
            if pseudo.distance_to_point(end) < EFF_ZERO {
                let found = crate::surface_solve::local_point_solve(
                    surface,
                    natural,
                    SolverOperation::Minimize,
                    end,
                    uv_domain.evaluate(0.5, 0.5),
                    None,
                )?;
                if found.found
                    && found.value < tolerance
                    && self.domain.contains_point(found.uv, 0.0)
                {
                    let solution = single(branch, t, found.uv, found.value);
                    let keep = filter_solution(
                        surface,
                        natural,
                        &self.line,
                        span,
                        u_conts,
                        v_conts,
                        t_conts,
                        candidate(&solution),
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
            for axis in 0..2 {
                let side_conts = if axis == 0 { u_conts } else { v_conts };
                if side_conts[j] > Continuity::C0 {
                    continue;
                }
                let (fixed, from_left) = if axis == 0 {
                    (
                        if j == 0 {
                            uv_domain.min.x
                        } else {
                            uv_domain.max.x
                        },
                        j == 0,
                    )
                } else {
                    (
                        if j == 0 {
                            uv_domain.min.y
                        } else {
                            uv_domain.max.y
                        },
                        j == 0,
                    )
                };
                let iso = IsoCurve::new(surface, axis, fixed, from_left);
                if PointSolveCurve::is_degenerate(&iso, EFF_ZERO) {
                    continue;
                }
                let iso_interval = if axis == 0 {
                    uv_domain.v_interval()
                } else {
                    uv_domain.u_interval()
                };
                let coincidence = crate::curve_solve::simple_coincidence_check(
                    &self.line,
                    span,
                    &iso,
                    iso_interval,
                    tolerance,
                    None,
                    Some(&patch.pseudo_box),
                )?;
                let as_solution = |pair: &CurvePair| {
                    let uv = if axis == 0 {
                        Vector2::new(fixed, pair.parameters[1])
                    } else {
                        Vector2::new(pair.parameters[1], fixed)
                    };
                    single(branch, pair.parameters[0], uv, pair.value)
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
                let mut coincident = false;
                for solution in &found {
                    if solution.kind == SolutionType::RangeOfValues {
                        coincident = true;
                        solver.add_sorted_solution(
                            &*self,
                            solution,
                            SortKey::ByAllParameters,
                            true,
                        );
                    } else if self.keep(solver, solution, uv_domain, span, conts)?
                        && axis == 1
                        && (span.length() < EFF_ZERO || iso_interval.length() < EFF_ZERO)
                    {
                        coincident = true;
                    }
                }
                if coincident {
                    return Ok(false);
                }
                if found.is_empty() {
                    let meet = crate::curve_solve::local_curve_intersect(
                        &iso,
                        iso_interval,
                        &self.line,
                        span,
                        tolerance,
                        [iso_interval.evaluate(0.5), span.evaluate(0.5)],
                    )?;
                    if let Some(meet) = meet {
                        let uv = if axis == 0 {
                            Vector2::new(fixed, meet.parameters[0])
                        } else {
                            Vector2::new(meet.parameters[0], fixed)
                        };
                        let solution = single(branch, meet.parameters[1], uv, meet.value);
                        self.keep(solver, &solution, uv_domain, span, conts)?;
                    }
                }
            }
        }
        Ok(false)
    }

    /// `SmLSIGlobalSolver::IdenticalSolutions`: SMLib's own comparison, and
    /// two points that meet at a singular point of the surface are one (the
    /// one off a seam, or the closer).
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
        if pair != SolutionPair::DifferentParameters
            || first.kind != SolutionType::SingleValue
            || second.kind != SolutionType::SingleValue
        {
            return (pair, Preferred::Neither);
        }
        let uv1 = Vector2::new(first.start.parameters[1], first.start.parameters[2]);
        let uv2 = Vector2::new(second.start.parameters[1], second.start.parameters[2]);
        let p1 = evaluate_point(self.surface, uv1).unwrap_or_default();
        let p2 = evaluate_point(self.surface, uv2).unwrap_or_default();
        if !p1.closer_than(self.tolerance_3d(), p2) {
            return (SolutionPair::DifferentParameters, Preferred::Neither);
        }
        let singular = |uv: Vector2| {
            is_singularity(
                self.surface,
                uv,
                EFF_ZERO,
                crate::scalar::ZONE_TOL_3D,
                false,
            )
            .map_or(true, |found| found.is_some())
        };
        if !singular(uv1) || !singular(uv2) {
            return (SolutionPair::DifferentParameters, Preferred::Neither);
        }
        let seam1 = is_on_seam(self.surface, uv1, crate::scalar::ZONE_TOL_3D).1;
        let seam2 = is_on_seam(self.surface, uv2, crate::scalar::ZONE_TOL_3D).1;
        let preferred = match (seam1, seam2) {
            (true, false) => Preferred::First,
            (false, true) => Preferred::Second,
            (false, false) => {
                if first.start.value < second.start.value {
                    Preferred::First
                } else {
                    Preferred::Second
                }
            }
            (true, true) => return (SolutionPair::DifferentParameters, Preferred::Neither),
        };
        (SolutionPair::Identical, preferred)
    }
}

impl LineSurfaceHooks<'_> {
    /// `m_d3dTolerance`, which `IdenticalSolutions` reads; the search runs
    /// with the distance tolerance.
    fn tolerance_3d(&self) -> f64 {
        self.tolerance
    }
}

/// `SmSurfaceEval` at `uv` along `direction`: the unit normal (with SMLib's
/// fallbacks at a vanishing first derivative) and the normal curvature
/// along the direction.
fn normal_and_section(
    surface: &BSplineSurface,
    uv: Vector2,
    direction: Vector3,
) -> Result<(Vector3, f64), Fault> {
    let d = evaluate(surface, uv, 2, 2, [true, true], true, true, true)?;
    let mut jet: SurfaceJet = [[Vector3::default(); 3]; 3];
    for k in 0..3 {
        for l in 0..3 - k {
            jet[k][l] = d[k * 3 + l];
        }
    }
    let mut normal = jet[1][0].cross(jet[0][1]);
    if normal.length_squared() < EFF_ZERO_SQ {
        let natural = surface.natural_domain();
        let along_u = jet[0][1].length_squared() >= jet[1][0].length_squared();
        normal = if along_u {
            let n = jet[1][1].cross(jet[0][1]);
            if uv.x > natural.u_interval().mid() {
                -n
            } else {
                n
            }
        } else {
            let n = jet[1][0].cross(jet[1][1]);
            if uv.y > natural.v_interval().mid() {
                -n
            } else {
                n
            }
        };
        if normal.length_squared() < EFF_ZERO_SQ {
            normal = evaluate_normal(surface, uv, [true, true]);
        } else {
            normal.unitize()?;
        }
    } else {
        normal.unitize()?;
    }
    jet[2][2] = normal;
    let (_, curvature) = normal_section_of(&mut jet, direction)?;
    Ok((normal, curvature))
}

/// `SmSurface::GlobalLineIntersect` (not as a ray) for a B-spline surface:
/// the meetings of `point + t * vector`, over `interval` when given, with
/// `surface` over `domain` within `tolerance`, as solutions `t, u, v` with
/// one object, the surface. Duplicates within a parameter tolerance merge.
/// Where every answer could be a point of a line lying in the surface, the
/// answers are the line's pieces over the surface (`DropAndTrimCurve`) as
/// ranges (a line without an interval refuses there); neighbouring answers the line joins within the surface merge
/// into ranges.
#[allow(clippy::too_many_lines)]
pub fn global_line_intersect(
    surface: &BSplineSurface,
    domain: Extent2,
    point: Vector3,
    vector: Vector3,
    interval: Option<Extent1>,
    tolerance: f64,
) -> Result<SolutionArray, SolverRefusal> {
    if vector.length_squared() < EFF_ZERO_SQ {
        return Err(Fault::InvalidInput.into());
    }
    let cache = SurfaceCache::for_trimmed_surface(surface).map_err(|refusal| match refusal {
        CacheRefusal::Failed(fault) => SolverRefusal::Failed(fault),
    })?;
    let tree_box: Extent3 = cache.trees.get(cache.top).ok_or(Fault::Degenerate)?.bbox;
    let scale = max(vector.length(), 0.1);
    let unit = vector / scale;
    let scaled = interval.map(|i| Extent1::new(i.min * scale, i.max * scale));
    let objects = LineSurfaceObjects { surface };
    let mut solver = GlobalSolver::with_one_tree(&cache.trees, &objects, Some(cache.top), 2);
    let length = unit.length();
    let mut hooks = LineSurfaceHooks {
        surface,
        cache: &cache,
        polar: PolarBoxes::new(surface, &cache),
        domain,
        point,
        vector: unit,
        interval: scaled,
        tolerance,
        line: Line::new(
            point,
            unit,
            Extent1::new(-1234567.0, 1234567.0),
            length,
            false,
            false,
        ),
    };
    let mut found = SolutionArray::new();
    solver.solve_it(
        &mut hooks,
        INTERSECT,
        SolutionRequested::All,
        tolerance,
        BIG_DOUBLE,
        &[],
        &mut found,
    )?;
    let mut solutions: Vec<Solution> = found.as_slice().to_vec();
    for solution in &mut solutions {
        solution.start.parameters[0] /= scale;
        if solution.kind == SolutionType::RangeOfValues {
            solution.end.parameters[0] /= scale;
        }
    }
    let line_max = interval.map_or(1.0, |i| i.max);
    let line_length = interval.map_or(1.0, Extent1::length);
    let tolerances = [
        max(EFF_ZERO_SQRT * (1.0 + line_max).abs(), line_length * 0.0001),
        max(
            EFF_ZERO_SQRT * (1.0 + domain.max.x).abs(),
            (domain.max.x - domain.min.x) * 0.0001,
        ),
        max(
            EFF_ZERO_SQRT * (1.0 + domain.max.y).abs(),
            (domain.max.y - domain.min.y) * 0.0001,
        ),
    ];
    let mut ii = 0;
    while ii < solutions.len() {
        let mut jj = ii + 1;
        while jj < solutions.len() {
            let (a, b) = (&solutions[ii], &solutions[jj]);
            let mut same = a.kind == b.kind;
            for k in 0..3 {
                if !same {
                    break;
                }
                if (a.start.parameters[k] - b.start.parameters[k]).abs() > tolerances[k]
                    || (a.kind == SolutionType::RangeOfValues
                        && (a.end.parameters[k] - b.end.parameters[k]).abs() > tolerances[k])
                {
                    same = false;
                }
            }
            if same {
                if solutions[ii].start.value > solutions[jj].start.value {
                    solutions[ii] = solutions[jj];
                }
                solutions.remove(jj);
                continue;
            }
            jj += 1;
        }
        ii += 1;
    }
    let mut might_be_coincident = !solutions.is_empty();
    let size = 2.0 * tree_box.max_dimension();
    for solution in &solutions {
        let ends: &[SolutionEnd] = if solution.kind == SolutionType::RangeOfValues {
            &[solution.start, solution.end]
        } else {
            &[solution.start]
        };
        for end in ends {
            let uv = Vector2::new(end.parameters[1], end.parameters[2]);
            // SMLib reads whatever a failed evaluation leaves; the port
            // takes it as not coincident.
            match normal_and_section(surface, uv, vector) {
                Ok((normal, curvature))
                    if vector.is_perpendicular_to(normal, 0.0001)
                        && curvature.abs() <= 1.0 / size => {}
                _ => might_be_coincident = false,
            }
        }
    }
    if might_be_coincident {
        let (start, end, line_interval) = match interval {
            Some(i) => (point + i.min * vector, point + i.max * vector, i),
            None => (point, point + vector, Extent1::new(0.0, 1.0)),
        };
        // `SmBSplineCurve::CreateLineSegment`, then `EditParameterization`.
        if start.distance(end) < EFF_ZERO {
            return Ok(array(solutions));
        }
        let segment = crate::spline::BSplineCurve::create_canonical(
            3,
            1,
            &[start, end],
            &[2, 2],
            &[line_interval.min, line_interval.max],
            None,
            None,
        )?;
        let coincident = crate::curve_surface::simple_coincidence_check(
            surface,
            domain,
            &segment,
            line_interval,
            tolerance,
            true,
            None,
        )?;
        if coincident.is_some() {
            // Deliberate deviation from SMLib: without an interval it probes
            // and drops the unit segment only, so its ranges could omit
            // overlap outside `[0, 1]`. An unlimited coincident line refuses.
            if interval.is_none() {
                return Err(SolverRefusal::NotPorted(INTERSECT));
            }
            // The line's pieces over the surface, as ranges.
            let ranges = crate::surface_drop_curve::coincident_ranges(
                surface,
                domain,
                &segment,
                line_interval,
                tolerance,
                false,
                false,
            )?;
            let mut found = SolutionArray::new();
            for range in ranges {
                found.add(range);
            }
            return Ok(found);
        }
    }
    // Merge answers that the line joins within the surface into ranges:
    // SMLib's last pass.
    let mut processed = 0;
    let mut ii = 1;
    while ii < solutions.len() {
        let last_end = if solutions[processed].kind == SolutionType::RangeOfValues {
            solutions[processed].end
        } else {
            solutions[processed].start
        };
        for jj in ii + 1..solutions.len() {
            if (solutions[ii].start.parameters[0] - solutions[jj].start.parameters[0]).abs()
                < EFF_ZERO
            {
                let du2 = (solutions[jj].start.parameters[1] - last_end.parameters[1]).abs();
                let dv2 = (solutions[jj].start.parameters[2] - last_end.parameters[2]).abs();
                let du1 = (solutions[ii].start.parameters[1] - last_end.parameters[1]).abs();
                let dv1 = (solutions[ii].start.parameters[2] - last_end.parameters[2]).abs();
                if du2 + dv2 < du1 + dv1 {
                    solutions.swap(ii, jj);
                }
            } else {
                break;
            }
        }
        let current = solutions[ii];
        let combine = combines(
            surface,
            domain,
            point,
            vector,
            tolerance,
            &last_end,
            &current.start,
        )?;
        if combine {
            let last = &mut solutions[processed];
            let original = last.kind;
            last.kind = SolutionType::RangeOfValues;
            if current.kind == SolutionType::RangeOfValues {
                if last_end.parameters[0] < current.start.parameters[0] {
                    last.end = current.end;
                }
            } else {
                let short = interval.map_or(EFF_ZERO * 100.0, |i| EFF_ZERO * 100.0 * i.length());
                if (last.start.parameters[0] - current.start.parameters[0]).abs() < short {
                    last.kind = original;
                } else {
                    last.end = current.start;
                    if last.start.parameters[0] > last.end.parameters[0] {
                        std::mem::swap(&mut last.start, &mut last.end);
                    }
                }
            }
        } else {
            processed += 1;
            if ii != processed {
                solutions[processed] = current;
            }
        }
        ii += 1;
    }
    if solutions.len() > processed + 1 {
        solutions.truncate(processed + 1);
    }
    Ok(array(solutions))
}

/// Whether the line between two answers stays on the surface (SMLib's
/// `bDoCombine`): the answers follow each other along the line, do not
/// straddle a seam, and the line's middle between them drops onto the
/// surface within the tolerance.
fn combines(
    surface: &BSplineSurface,
    domain: Extent2,
    point: Vector3,
    vector: Vector3,
    tolerance: f64,
    last: &SolutionEnd,
    current: &SolutionEnd,
) -> Result<bool, SolverRefusal> {
    let (lp, cp) = (last.parameters, current.parameters);
    if lp[0] > cp[0] + EFF_ZERO_SQRT {
        return Ok(false);
    }
    let mut between = Extent1::new(lp[0], lp[0]);
    between.add_value(cp[0]);
    let last_uv = Vector2::new(lp[1], lp[2]);
    let current_uv = Vector2::new(cp[1], cp[2]);
    if between.length() < EFF_ZERO_SQRT
        && ((last_uv.x - current_uv.x).abs() > 0.99 * domain.x_length()
            || (last_uv.y - current_uv.y).abs() > 0.99 * domain.y_length())
    {
        return Ok(false);
    }
    let Ok(uv_between) = Extent2::new(last_uv, current_uv).intersect(domain) else {
        return Ok(false);
    };
    let largest = max(max(lp[1].abs(), cp[1].abs()), max(lp[2].abs(), cp[2].abs()));
    let grow = 1.1 * EFF_ZERO * (1.0 + largest);
    let grown = Extent2::new(
        Vector2::new(uv_between.min.x - grow, uv_between.min.y - grow),
        Vector2::new(uv_between.max.x + grow, uv_between.max.y + grow),
    );
    let Ok(uv_between) = grown.intersect(domain) else {
        return Ok(false);
    };
    let middle = point + between.evaluate(0.5) * vector;
    let found = crate::surface_solve::local_point_solve(
        surface,
        uv_between,
        SolverOperation::Minimize,
        middle,
        uv_between.evaluate(0.5, 0.5),
        None,
    )?;
    Ok(found.found && found.value <= tolerance)
}

fn array(solutions: Vec<Solution>) -> SolutionArray {
    let mut array = SolutionArray::new();
    for solution in solutions {
        array.add(solution);
    }
    array
}
