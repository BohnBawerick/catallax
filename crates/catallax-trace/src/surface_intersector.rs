// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline surfaces and planes of
// `spline.rs` without trimming faces; statuses become `Fault`, the two
// surface trees share one arena, and the intersector's state is a struct
// whose stages are methods. Grazing side planes refuse with
// `SolverRefusal::NotPorted`.
// Translated from SMLib SmSurfaceIntersector.cpp (the
// SmSurfaceIntersector constructor and the search stage of DoIntersection)
// and SmAdvSurfaceIntersector.cpp (the constructor, SolveBranch, LocalSolve
// for SM_SO_INTERSECTION_TEST and SM_SO_INTERSECT, sm_FindLargestNonLeaf and
// FindBestSubdivisionNode); rust-brep has no port.
//! The surface-surface intersector of step 8 (Gate B): its set-up and the
//! search for the region where the two surfaces may meet.

use crate::Fault;
use crate::extent::{Extent2, Extent3};
use crate::global_solver::{
    Branch, GlobalSolver, SearchTrees, SolutionRequested, SolverHooks, SolverOperation,
    SolverRefusal, Subdivision,
};
use crate::scalar::{BIG_DOUBLE, EFF_ZERO, PI};
use crate::solution::{NoObjects, NodeId, Solution, SolutionArray, SolutionType};
use crate::spline::BSplineSurface;
use crate::surface_cache::{CacheRefusal, SurfaceCache, SurfaceNodeData};
use crate::surface_eval::evaluate_point;
use crate::surface_solve::is_closed;
use crate::vector::Vector3;

const INTERSECT: SolverOperation = SolverOperation::Intersect;

/// The state of `SmAdvSurfaceIntersector` after construction.
pub struct Intersector<'a> {
    /// `m_cpSurface`.
    pub surfaces: [&'a BSplineSurface; 2],
    /// `m_vUVDomain`: each given domain, cut to the natural one.
    pub domains: [Extent2; 2],
    /// `m_bClosedU`.
    pub closed_u: [bool; 2],
    /// `m_bClosedV`.
    pub closed_v: [bool; 2],
    /// `m_bUseSurfaceEdges`.
    pub use_surface_edges: [bool; 2],
    /// `m_dThisAngTolRad`.
    pub angle_tolerance: f64,
    /// `m_dThisApproxTol3d`.
    pub approx_tolerance: f64,
    /// `m_dSteppingParallelTolRadians`.
    pub stepping_parallel_tolerance: f64,
    /// `m_dStartParallelTolRadians`.
    pub start_parallel_tolerance: f64,
    /// `m_bDoingSelfIntersection`.
    pub self_intersection: bool,
}

/// What the search stage of `DoIntersection` found: the union of the
/// meeting boxes of the leaf pairs (`m_vStart[0..6]`), how many pairs met
/// (`m_vStart[6]`) and the unions of their parameter domains (`m_vEnd`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Overlap {
    /// The union of the leaf-pair meeting boxes.
    pub bbox: Extent3,
    /// The number of meeting leaf pairs.
    pub count: f64,
    /// The union of the meeting patches' domains on each surface.
    pub domains: [Extent2; 2],
}

/// The surface caches of both surfaces in one arena: the first surface's
/// tree, then the second's from `offset`.
pub struct IntersectorTrees {
    /// The first surface's cache.
    pub first: SurfaceCache,
    /// The second surface's cache.
    pub second: SurfaceCache,
    /// Both trees.
    pub trees: SearchTrees,
    /// Where the second tree starts.
    pub offset: u32,
}

impl IntersectorTrees {
    fn patch(&self, id: NodeId) -> Option<&crate::surface_cache::PatchData> {
        let (cache, index) = if id.0 < self.offset {
            (&self.first, id.0)
        } else {
            (&self.second, id.0 - self.offset)
        };
        match cache.node_data.get(index as usize)? {
            SurfaceNodeData::Patch(patch) => cache.patches.get(*patch),
            SurfaceNodeData::Head(_) => None,
        }
    }
}

fn cache_of(surface: &BSplineSurface) -> Result<SurfaceCache, SolverRefusal> {
    SurfaceCache::for_trimmed_surface(surface).map_err(|refusal| match refusal {
        CacheRefusal::Failed(fault) => SolverRefusal::Failed(fault),
        CacheRefusal::NotPorted(_) => SolverRefusal::NotPorted(INTERSECT),
    })
}

impl<'a> Intersector<'a> {
    /// `SmSurfaceIntersector` and `SmAdvSurfaceIntersector` constructors
    /// (not from the fillet intersector).
    pub fn new(
        surface1: &'a BSplineSurface,
        domain1: Extent2,
        surface2: &'a BSplineSurface,
        domain2: Extent2,
    ) -> Self {
        let clip = |surface: &BSplineSurface, domain: Extent2| {
            let natural = surface.natural_domain();
            if domain.is_contained_by(natural, EFF_ZERO) {
                domain
            } else {
                // SMLib warns and keeps the natural domain when the two do
                // not overlap.
                natural.intersect(domain).unwrap_or(natural)
            }
        };
        let closed = |surface: &BSplineSurface, domain: Extent2, axis: usize| {
            is_closed(surface, domain, axis, None)
        };
        let p1 = evaluate_point(surface1, domain1.min).unwrap_or_default();
        let middle = evaluate_point(surface1, domain1.evaluate(0.5, 0.5)).unwrap_or_default();
        let p2 = evaluate_point(surface1, domain1.max).unwrap_or_default();
        let size = p1.distance(middle) + middle.distance(p2);
        Self {
            surfaces: [surface1, surface2],
            domains: [clip(surface1, domain1), clip(surface2, domain2)],
            closed_u: [closed(surface1, domain1, 0), closed(surface2, domain2, 0)],
            closed_v: [closed(surface1, domain1, 1), closed(surface2, domain2, 1)],
            use_surface_edges: [true, true],
            angle_tolerance: 20.0 * PI / 180.0,
            approx_tolerance: size / 1000.0,
            stepping_parallel_tolerance: EFF_ZERO,
            start_parallel_tolerance: 1.0e-4 * PI / 180.0,
            self_intersection: std::ptr::eq(surface1, surface2),
        }
    }

    /// Both surfaces' caches in one arena.
    pub fn trees(&self) -> Result<IntersectorTrees, SolverRefusal> {
        let first = cache_of(self.surfaces[0])?;
        let second = cache_of(self.surfaces[1])?;
        let mut trees = SearchTrees::new();
        trees.append(&first.trees);
        let offset = trees.append(&second.trees);
        Ok(IntersectorTrees {
            first,
            second,
            trees,
            offset,
        })
    }

    /// The search stage of `DoIntersection`: `None` when the root boxes are
    /// apart or no leaf pair meets; otherwise the [`Overlap`]. A grazing
    /// meeting (one leaf pair in a box thinner than twice the tolerance)
    /// needs SMLib's `GrazingBoxProcess` and refuses.
    pub fn overlap(&self, trees: &IntersectorTrees) -> Result<Option<Overlap>, SolverRefusal> {
        let Some(overlap) = self.overlap_box(trees)? else {
            return Ok(None);
        };
        if self.grazes(&overlap) {
            return Err(SolverRefusal::NotPorted(INTERSECT));
        }
        Ok(Some(overlap))
    }

    /// The search stage of `DoIntersection` without the grazing test:
    /// `None` when the root boxes are apart or no leaf pair meets.
    pub fn overlap_box(&self, trees: &IntersectorTrees) -> Result<Option<Overlap>, SolverRefusal> {
        if !trees.root_boxes_meet(self.approx_tolerance)? {
            return Ok(None);
        }
        self.search(trees)
    }

    /// Whether one leaf pair meets in a box thinner than twice the
    /// tolerance, where `DoIntersection` first tries `GrazingBoxProcess`.
    pub fn grazes(&self, overlap: &Overlap) -> bool {
        if overlap.count >= 2.0 {
            return false;
        }
        let size = overlap.bbox.max - overlap.bbox.min;
        let limit = 2.0 * self.approx_tolerance;
        size.x < limit || size.y < limit || size.z < limit
    }

    /// `SolveIt(SM_SO_INTERSECTION_TEST)` over both trees with the
    /// advanced intersector's overrides: the first solution, if any.
    pub fn search(&self, trees: &IntersectorTrees) -> Result<Option<Overlap>, SolverRefusal> {
        let objects = NoObjects;
        let mut solver = GlobalSolver::with_two_trees(
            &trees.trees,
            &objects,
            Some(trees.first.top),
            Some(trees.second_top()),
            4,
        );
        let mut hooks = AdvancedHooks {
            trees,
            surfaces: self.surfaces,
            domains: self.domains,
            self_intersection: self.self_intersection,
        };
        let mut solutions = SolutionArray::new();
        solver.solve_it(
            &mut hooks,
            SolverOperation::IntersectionTest,
            SolutionRequested::All,
            self.approx_tolerance,
            BIG_DOUBLE,
            &[],
            &mut solutions,
        )?;
        let Ok(first) = solutions.get_at(0) else {
            return Ok(None);
        };
        let p = first.start.parameters;
        let e = first.end.parameters;
        let v2 = crate::vector::Vector2::new;
        Ok(Some(Overlap {
            bbox: Extent3::new(
                Vector3::new(p[0], p[1], p[2]),
                Vector3::new(p[3], p[4], p[5]),
            ),
            count: p[6],
            domains: [
                Extent2::new(v2(e[0], e[1]), v2(e[2], e[3])),
                Extent2::new(v2(e[4], e[5]), v2(e[6], e[7])),
            ],
        }))
    }
}

/// `SmBoundaryPlaneTestResult` without the grazing outcome, which refuses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundaryPlaneTest {
    /// `SM_BR_UNKNOWN`: a meeting is possible; more work is needed.
    Unknown,
    /// `SM_BR_DISJOINT`: a pair of opposed side planes separates the
    /// surfaces.
    Disjoint,
}

impl Intersector<'_> {
    /// `DoBoundaryPlaneTest`: compares each pair of opposed side planes of
    /// the two surfaces. Surfaces wholly beyond one pair are disjoint; a pair
    /// on which they only graze needs the curve-curve intersection of the
    /// side curves (`SM_BR_GRAZE`), which refuses.
    pub fn boundary_plane_test(&self, bbox: Extent3) -> Result<BoundaryPlaneTest, SolverRefusal> {
        let size = (bbox.max - bbox.min).length();
        let angle_tolerance = 1e-4;
        let distance_tolerance = crate::scalar::min(1e-2, self.approx_tolerance / size) / 100.0;
        let planes = |surface: &BSplineSurface| {
            crate::surface_cache::boundary_planes(surface, surface.natural_domain())
                .map(|planes| planes.planes)
                .unwrap_or([None; 4])
        };
        let (planes1, planes2) = (planes(self.surfaces[0]), planes(self.surfaces[1]));
        for plane1 in planes1.iter().flatten() {
            for plane2 in planes2.iter().flatten() {
                if plane1.normal.dot(plane2.normal) > -1.0 + angle_tolerance {
                    continue;
                }
                let mut disjoint = true;
                let mut overlap = false;
                let largest =
                    crate::scalar::max(plane1.base.max_dimension(), plane2.base.max_dimension());
                let tolerance_sum = crate::scalar::max(
                    distance_tolerance / 1000.0,
                    plane1.tolerance + plane2.tolerance,
                ) * (1.0 + largest);
                for (sample1, sample2) in plane1.samples.iter().zip(&plane2.samples) {
                    let distance2 = plane2.normal.dot(*sample1 - plane2.base);
                    let distance1 = plane1.normal.dot(*sample2 - plane1.base);
                    if distance2 <= tolerance_sum || distance1 <= tolerance_sum {
                        disjoint = false;
                    }
                    if distance1 < -distance_tolerance || distance2 < -distance_tolerance {
                        overlap = true;
                        break;
                    }
                }
                if disjoint {
                    return Ok(BoundaryPlaneTest::Disjoint);
                }
                if overlap {
                    continue;
                }
                return Err(SolverRefusal::NotPorted(INTERSECT));
            }
        }
        Ok(BoundaryPlaneTest::Unknown)
    }
}

impl IntersectorTrees {
    /// The second tree's root in the shared arena.
    pub fn second_top(&self) -> NodeId {
        NodeId(self.second.top.0 + self.offset)
    }

    /// Whether the first root box, grown by `tolerance`, meets the second.
    pub fn root_boxes_meet(&self, tolerance: f64) -> Result<bool, SolverRefusal> {
        let root = |id: NodeId| {
            self.trees
                .get(id)
                .map(|node| node.bbox)
                .ok_or(SolverRefusal::Failed(Fault::Degenerate))
        };
        let mut bbox = root(self.first.top)?;
        bbox.expand_absolute(tolerance);
        Ok(!bbox.are_disjoint(root(self.second_top())?, 0.0))
    }
}

/// The overrides of `SmAdvSurfaceIntersector` that the box search and the
/// interior search use.
struct AdvancedHooks<'a> {
    trees: &'a IntersectorTrees,
    surfaces: [&'a BSplineSurface; 2],
    domains: [Extent2; 2],
    self_intersection: bool,
}

impl AdvancedHooks<'_> {
    /// `SmAdvSurfaceIntersector::LocalSolve` and the answer that its
    /// `SolveBranch` saves for `SM_SO_INTERSECT`: for a leaf pair whose
    /// pseudo boxes meet, the patch midpoints, the larger patch's moved
    /// onto its tangent plane at the smaller's when that is closer.
    fn closest_midpoints(
        &self,
        solver: &mut GlobalSolver<'_>,
        ids: [NodeId; 2],
        patches: [&crate::surface_cache::PatchData; 2],
    ) -> Result<bool, SolverRefusal> {
        let mut box1 = patches[0].pseudo_box;
        box1.expand_absolute(solver.tolerance);
        if box1.are_disjoint(&patches[1].pseudo_box) {
            return Ok(false);
        }
        let node = |id: NodeId| {
            self.trees
                .trees
                .get(id)
                .ok_or(SolverRefusal::Failed(Fault::Degenerate))
        };
        let size = |id: NodeId| -> Result<f64, SolverRefusal> {
            let bbox = node(id)?.bbox;
            Ok((bbox.max - bbox.min).length())
        };
        let (size1, size2) = (size(ids[0])?, size(ids[1])?);
        let mut uv = [0, 1].map(|s| patches[s].uv_domain.evaluate(0.5, 0.5));
        let points = [
            evaluate_point(self.surfaces[0], uv[0])?,
            evaluate_point(self.surfaces[1], uv[1])?,
        ];
        let mut distance = points[0].distance(points[1]);
        let mover = if size1 > 4.0 * size2 {
            Some(0)
        } else if size2 > 4.0 * size1 {
            Some(1)
        } else {
            None
        };
        if let Some(s) = mover {
            let other = 1 - s;
            let moved = crate::curve_surface::drop_point_to_tangent_plane(
                self.surfaces[s],
                uv[s],
                [true, true],
                points[other],
            )?;
            let moved = patches[s].uv_domain.clamp_point(moved);
            let new_distance = evaluate_point(self.surfaces[s], moved)?.distance(points[other]);
            if new_distance < distance {
                uv[s] = moved;
                distance = new_distance;
            }
        }
        if distance < BIG_DOUBLE {
            let mut solution = Solution {
                kind: SolutionType::SingleValue,
                variables: 4,
                ..Solution::default()
            };
            solution.start.parameters[..4].copy_from_slice(&[uv[0].x, uv[0].y, uv[1].x, uv[1].y]);
            solution.nodes[0] = Some(ids[0]);
            solution.nodes[1] = Some(ids[1]);
            solver.solutions.add(solution);
        }
        Ok(false)
    }
}

impl Intersector<'_> {
    /// `SolveIt(SM_SO_INTERSECT)` over both trees with the advanced
    /// intersector's overrides: one parameter pair per leaf pair whose
    /// pseudo boxes meet, in search order (`FindInteriorCurves`' guesses).
    pub fn interior_guesses(
        &self,
        trees: &IntersectorTrees,
    ) -> Result<Vec<[crate::vector::Vector2; 2]>, SolverRefusal> {
        let objects = NoObjects;
        let mut solver = GlobalSolver::with_two_trees(
            &trees.trees,
            &objects,
            Some(trees.first.top),
            Some(trees.second_top()),
            4,
        );
        let mut hooks = AdvancedHooks {
            trees,
            surfaces: self.surfaces,
            domains: self.domains,
            self_intersection: self.self_intersection,
        };
        let mut solutions = SolutionArray::new();
        solver.solve_it(
            &mut hooks,
            SolverOperation::Intersect,
            SolutionRequested::All,
            self.approx_tolerance,
            BIG_DOUBLE,
            &[],
            &mut solutions,
        )?;
        let v2 = crate::vector::Vector2::new;
        Ok(solutions
            .as_slice()
            .iter()
            .map(|s| {
                let p = s.start.parameters;
                [v2(p[0], p[1]), v2(p[2], p[3])]
            })
            .collect())
    }
}

impl SolverHooks for AdvancedHooks<'_> {
    /// `sm_FindLargestNonLeaf` as `SmAdvSurfaceIntersector` uses it: split
    /// the second surface's node once the first's is a leaf, else the
    /// first's.
    fn find_best_subdivision_node(
        &self,
        _solver: &GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<Subdivision, SolverRefusal> {
        let is_leaf = |id: Option<NodeId>| -> Result<bool, SolverRefusal> {
            let node = id
                .and_then(|id| self.trees.trees.get(id))
                .ok_or(SolverRefusal::Failed(Fault::Degenerate))?;
            let Some([child, _]) = node.children else {
                return Ok(true);
            };
            let child = self
                .trees
                .trees
                .get(child)
                .ok_or(SolverRefusal::Failed(Fault::Degenerate))?;
            Ok(node.aux == crate::global_solver::AuxData::BezierSurface
                && child.aux != crate::global_solver::AuxData::BezierSurface)
        };
        let tree = if is_leaf(branch[0])? { 1 } else { 0 };
        let _ = is_leaf(branch[1])?;
        Ok(Subdivision {
            tree,
            swap_traversal_order: false,
        })
    }

    /// `SmAdvSurfaceIntersector::LocalSolve` for
    /// `SM_SO_INTERSECTION_TEST`: the union of the meeting boxes and patch
    /// domains of the leaf pairs whose pseudo boxes meet.
    fn local_solve(
        &mut self,
        solver: &mut GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<bool, SolverRefusal> {
        let (Some(id1), Some(id2)) = (branch[0], branch[1]) else {
            return Err(Fault::Degenerate.into());
        };
        let (Some(patch1), Some(patch2)) = (self.trees.patch(id1), self.trees.patch(id2)) else {
            return Err(Fault::Degenerate.into());
        };
        let (domain1, domain2) = (patch1.uv_domain, patch2.uv_domain);
        if domain1.are_disjoint(self.domains[0], 0.0) || domain2.are_disjoint(self.domains[1], 0.0)
        {
            return Ok(false);
        }
        if self.self_intersection && !domain1.are_disjoint(domain2, 0.0) {
            return Ok(false);
        }
        if solver.operation == SolverOperation::Intersect {
            return self.closest_midpoints(solver, [id1, id2], [patch1, patch2]);
        }
        if solver.operation != SolverOperation::IntersectionTest {
            return Err(SolverRefusal::NotPorted(solver.operation));
        }
        let tolerance = solver.tolerance;
        let mut box1 = patch1.pseudo_box;
        box1.expand_absolute(tolerance);
        if box1.are_disjoint(&patch2.pseudo_box) {
            return Ok(false);
        }
        let node1 = self
            .trees
            .trees
            .get(id1)
            .ok_or(SolverRefusal::Failed(Fault::Degenerate))?;
        let node2 = self
            .trees
            .trees
            .get(id2)
            .ok_or(SolverRefusal::Failed(Fault::Degenerate))?;
        let mut bbox = node1.bbox;
        bbox.expand_absolute(tolerance);
        // SMLib ignores the status; disjoint boxes leave an empty result.
        let meeting = bbox
            .intersect(node2.bbox)
            .map_err(|_| SolverRefusal::NotPorted(INTERSECT))?;
        let mut solutions = std::mem::take(&mut solver.solutions);
        if solutions.is_empty() {
            let mut solution = Solution {
                kind: SolutionType::SingleValue,
                variables: 7,
                ..Solution::default()
            };
            solution.nodes[0] = Some(id1);
            solution.nodes[1] = Some(id2);
            let p = &mut solution.start.parameters;
            p[..3].copy_from_slice(&[meeting.min.x, meeting.min.y, meeting.min.z]);
            p[3..6].copy_from_slice(&[meeting.max.x, meeting.max.y, meeting.max.z]);
            p[6] = 1.0;
            let e = &mut solution.end.parameters;
            e[..4].copy_from_slice(&[domain1.min.x, domain1.min.y, domain1.max.x, domain1.max.y]);
            e[4..8].copy_from_slice(&[domain2.min.x, domain2.min.y, domain2.max.x, domain2.max.y]);
            solutions.add(solution);
        } else {
            let mut solution = solutions.get_at(0)?;
            let p = solution.start.parameters;
            let old = Extent3::new(
                Vector3::new(p[0], p[1], p[2]),
                Vector3::new(p[3], p[4], p[5]),
            );
            let merged = meeting.union(old);
            let e = solution.end.parameters;
            let v2 = crate::vector::Vector2::new;
            let old1 = Extent2::new(v2(e[0], e[1]), v2(e[2], e[3]));
            let old2 = Extent2::new(v2(e[4], e[5]), v2(e[6], e[7]));
            let d1 = old1.union(domain1);
            let d2 = old2.union(domain2);
            let p = &mut solution.start.parameters;
            p[6] += 1.0;
            p[..3].copy_from_slice(&[merged.min.x, merged.min.y, merged.min.z]);
            p[3..6].copy_from_slice(&[merged.max.x, merged.max.y, merged.max.z]);
            let e = &mut solution.end.parameters;
            e[..4].copy_from_slice(&[d1.min.x, d1.min.y, d1.max.x, d1.max.y]);
            e[4..8].copy_from_slice(&[d2.min.x, d2.min.y, d2.max.x, d2.max.y]);
            solutions.set_at(0, solution)?;
        }
        solver.solutions = solutions;
        Ok(false)
    }
}
