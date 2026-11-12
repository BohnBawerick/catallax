// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: translated from SMLib SmGlobalSolver.h and
// SmGlobalSolver.cpp, SmTree.h (the tree node fields the solver reads) and
// SmCoreTypes.h (the operation, category, request and sort enumerations).
// SMLib's virtual methods become the `SolverHooks` trait with SMLib's base
// behaviour as the default; trees become an index arena; raw pointers to
// trees, nodes, objects and the output array become owned or borrowed Rust
// values; debug drawing and `Dump` are not carried over; failures that SMLib
// meets through a null pointer return a typed refusal. The branch of
// `BranchMayContainAnswers` that samples a curve on a surface
// (`SM_SO_SURFACE_PROJECTED_INTERSECT`) reads the second tree's owner and
// spans from a `SurfaceCurveSpans` the caller gives the solver, where SMLib
// reaches them through the tree node's pointers.
//! SMLib's global search solver: traversal of two or three bounding-box
//! trees, culling by operation, and the accumulation of answers.
//!
//! The solver proposes candidates. Its tolerances decide what to search, never
//! what Catallax accepts.

use crate::Fault;
use crate::curve_on_surface::CurveOnSurface;
use crate::extent::{Extent1, Extent3};
use crate::geometry::{line_closest_point, point_project_to_plane};
use crate::scalar::{BIG_DOUBLE, EFF_ZERO, EFF_ZERO_SQRT, UNDEFINED, max, sqrt};
use crate::search_math::{
    maximum_distance_squared, minimum_distance, minimum_distance_squared,
    perspective_project_point_to_plane, perspective_project_vector_to_plane,
    rotate_project_point_to_plane,
};
use crate::solution::{
    NodeId, Preferred, Solution, SolutionArray, SolutionObjects, SolutionPair, SolutionType,
    SortKey,
};
use crate::vector::Vector3;

/// `SM_GS_MAX_TREES`: the most trees one solver walks together.
pub const MAX_TREES: usize = 3;

/// One tree node from each tree, in tree order; unused slots are `None`.
pub type Branch = [Option<NodeId>; MAX_TREES];

/// `SmAuxDataType`: what a tree node's auxiliary data holds. The
/// discriminants are SMLib's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AuxData {
    /// `SM_AD_NONE`.
    None = 0,
    /// `SM_AD_BEZIER_CURVE`: a Bezier span, on curve-cache leaves.
    BezierCurve = 1,
    /// `SM_AD_CURVE_DATA`: curve-cache internal nodes.
    CurveData = 2,
    /// `SM_AD_BEZIER_SURFACE`: a Bezier patch, marking decomposition nodes.
    BezierSurface = 3,
    /// `SM_AD_BEZIER_SURFACE_HEAD`: organising ancestors of decomposition nodes.
    BezierSurfaceHead = 4,
    /// `SM_AD_SURFACE_DATA`.
    SurfaceData = 5,
    /// `SM_AD_AUX_DATA`: non-decomposition descendants of decomposition nodes.
    AuxData = 6,
    /// `SM_AD_CUBEAUX_DATA`.
    CubeAuxData = 7,
    /// `SM_AD_OBJECT_LIST`: an object list, on point and spatial trees.
    ObjectList = 8,
}

impl AuxData {
    /// The auxiliary data type with this SMLib discriminant.
    pub fn from_code(code: usize) -> Option<Self> {
        [
            Self::None,
            Self::BezierCurve,
            Self::CurveData,
            Self::BezierSurface,
            Self::BezierSurfaceHead,
            Self::SurfaceData,
            Self::AuxData,
            Self::CubeAuxData,
            Self::ObjectList,
        ]
        .get(code)
        .copied()
    }
}

/// The fields of `SmTreeNode` the solver reads.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TreeNode {
    /// `m_sBBox`: the bounding box of this node.
    pub bbox: Extent3,
    /// `m_pChild1` and `m_pChild2`: both children, or none for a leaf.
    pub children: Option<[NodeId; 2]>,
    /// `m_eAuxDataType`.
    pub aux: AuxData,
}

/// The nodes of every tree one solver walks, as an index arena.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SearchTrees {
    nodes: Vec<TreeNode>,
}

impl SearchTrees {
    /// An empty arena.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a node and returns its identity.
    ///
    /// # Panics
    ///
    /// Panics when the arena would hold more than `u32::MAX` nodes.
    pub fn add(&mut self, node: TreeNode) -> NodeId {
        self.nodes.push(node);
        NodeId(u32::try_from(self.nodes.len() - 1).expect("fewer than 2^32 tree nodes"))
    }

    /// Appends every node of `other`, renumbering its child links, and
    /// returns the offset added to `other`'s identities.
    ///
    /// # Panics
    ///
    /// Panics when the arena would hold more than `u32::MAX` nodes.
    pub fn append(&mut self, other: &Self) -> u32 {
        let offset = u32::try_from(self.nodes.len()).expect("fewer than 2^32 tree nodes");
        self.nodes.extend(other.nodes.iter().map(|node| {
            TreeNode {
                children: node
                    .children
                    .map(|[a, b]| [NodeId(a.0 + offset), NodeId(b.0 + offset)]),
                ..*node
            }
        }));
        offset
    }

    /// The node with this identity.
    pub fn get(&self, id: NodeId) -> Option<&TreeNode> {
        self.nodes.get(id.0 as usize)
    }

    /// How many nodes the arena holds.
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Whether the arena holds no node.
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }
}

/// `SmSolverOperationType`: what the local solver is asked to do. The
/// discriminants are SMLib's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolverOperation {
    /// `SM_SO_INTERSECT`.
    Intersect = 0,
    /// `SM_SO_INTERSECT_WIREFRAME`.
    IntersectWireframe = 1,
    /// `SM_SO_INTERSECTION_TEST`.
    IntersectionTest = 2,
    /// `SM_SO_MINIMIZE`.
    Minimize = 3,
    /// `SM_SO_MAXIMIZE`.
    Maximize = 4,
    /// `SM_SO_NORMALIZE`.
    Normalize = 5,
    /// `SM_SO_AT_DISTANCE`.
    AtDistance = 6,
    /// `SM_SO_FIND`.
    Find = 7,
    /// `SM_SO_RAYFIRE`.
    Rayfire = 8,
    /// `SM_SO_3D_SIGNED_DIRECTED_MINIMIZE`.
    SignedDirectedMinimize3d = 9,
    /// `SM_SO_3D_SIGNED_DIRECTED_MAXIMIZE`.
    SignedDirectedMaximize3d = 10,
    /// `SM_SO_SIGNED_DIRECTED_MINIMIZE`.
    SignedDirectedMinimize = 11,
    /// `SM_SO_DIRECTED_MINIMIZE`.
    DirectedMinimize = 12,
    /// `SM_SO_DIRECTED_MAXIMIZE`.
    DirectedMaximize = 13,
    /// `SM_SO_PROJECTED_MINIMIZE`.
    ProjectedMinimize = 14,
    /// `SM_SO_PROJECTED_MAXIMIZE`.
    ProjectedMaximize = 15,
    /// `SM_SO_PROJECTED_INTERSECT`.
    ProjectedIntersect = 16,
    /// `SM_SO_ROTATED_PROJECTED_INTERSECT`.
    RotatedProjectedIntersect = 17,
    /// `SM_SO_PERSPECTIVE_INTERSECT`.
    PerspectiveIntersect = 18,
    /// `SM_SO_SURFACE_PROJECTED_INTERSECT`.
    SurfaceProjectedIntersect = 19,
    /// `SM_SO_ANGLE_MINIMIZE`.
    AngleMinimize = 20,
    /// `SM_SO_SIGNED_ANGLE_MINIMIZE`.
    SignedAngleMinimize = 21,
    /// `SM_SO_PROJECTED_TANGENCY`.
    ProjectedTangency = 22,
    /// `SM_SO_SIGNED_PIVOT_MINIMIZE`.
    SignedPivotMinimize = 23,
    /// `SM_SO_PROJECTED_TANGENT_THROUGH_POINT`.
    ProjectedTangentThroughPoint = 24,
}

impl SolverOperation {
    /// Every operation, indexed by its SMLib discriminant.
    pub const ALL: [Self; 25] = [
        Self::Intersect,
        Self::IntersectWireframe,
        Self::IntersectionTest,
        Self::Minimize,
        Self::Maximize,
        Self::Normalize,
        Self::AtDistance,
        Self::Find,
        Self::Rayfire,
        Self::SignedDirectedMinimize3d,
        Self::SignedDirectedMaximize3d,
        Self::SignedDirectedMinimize,
        Self::DirectedMinimize,
        Self::DirectedMaximize,
        Self::ProjectedMinimize,
        Self::ProjectedMaximize,
        Self::ProjectedIntersect,
        Self::RotatedProjectedIntersect,
        Self::PerspectiveIntersect,
        Self::SurfaceProjectedIntersect,
        Self::AngleMinimize,
        Self::SignedAngleMinimize,
        Self::ProjectedTangency,
        Self::SignedPivotMinimize,
        Self::ProjectedTangentThroughPoint,
    ];

    /// The operation with this SMLib discriminant.
    pub fn from_code(code: usize) -> Option<Self> {
        Self::ALL.get(code).copied()
    }

    /// Whether `SolveIt` marks the operation as projected
    /// (`m_bProjectedOperation`).
    pub fn is_projected(self) -> bool {
        matches!(
            self,
            Self::SignedDirectedMinimize
                | Self::DirectedMinimize
                | Self::ProjectedMinimize
                | Self::SignedAngleMinimize
                | Self::DirectedMaximize
                | Self::ProjectedMaximize
                | Self::SignedPivotMinimize
                | Self::SignedDirectedMaximize3d
        )
    }

    /// The category `SolveIt` files the operation under.
    pub fn category(self) -> OperationCategory {
        match self {
            Self::DirectedMinimize
            | Self::ProjectedMinimize
            | Self::Minimize
            | Self::AngleMinimize
            | Self::Rayfire => OperationCategory::Minimize,
            Self::DirectedMaximize | Self::ProjectedMaximize | Self::Maximize => {
                OperationCategory::Maximize
            }
            Self::SignedDirectedMinimize
            | Self::SignedAngleMinimize
            | Self::SignedPivotMinimize
            | Self::SignedDirectedMinimize3d
            | Self::SignedDirectedMaximize3d => OperationCategory::Signed,
            _ => OperationCategory::Other,
        }
    }

    fn is_maximizing(self) -> bool {
        matches!(
            self,
            Self::Maximize | Self::DirectedMaximize | Self::ProjectedMaximize
        )
    }
}

/// `SmOperationCategory`. The discriminants are SMLib's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OperationCategory {
    /// `SM_OPERATION_MINIMIZE`.
    Minimize = 0,
    /// `SM_OPERATION_MAXIMIZE`.
    Maximize = 1,
    /// `SM_OPERATION_SIGNED`.
    Signed = 2,
    /// `SM_OPERATION_OTHER`.
    Other = 3,
}

impl OperationCategory {
    /// The category with this SMLib discriminant.
    pub fn from_code(code: usize) -> Option<Self> {
        [Self::Minimize, Self::Maximize, Self::Signed, Self::Other]
            .get(code)
            .copied()
    }
}

/// `SmSolutionRequestedType`. The discriminants are SMLib's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolutionRequested {
    /// `SM_SR_SINGLE`: only the best solution.
    Single = 0,
    /// `SM_SR_ALL`: every solution within tolerance.
    All = 1,
    /// `SM_SR_FIND_AMBIGUITIES`: not used by SMLib.
    FindAmbiguities = 2,
    /// `SM_SR_NODES`: only the tree nodes where solutions may exist.
    Nodes = 3,
}

impl SolutionRequested {
    /// The request with this SMLib discriminant.
    pub fn from_code(code: usize) -> Option<Self> {
        [Self::Single, Self::All, Self::FindAmbiguities, Self::Nodes]
            .get(code)
            .copied()
    }
}

/// Why a ported solver routine produced no value.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolverRefusal {
    /// SMLib's `SM_ERR` or `SM_ERR_INVALID_INPUT`, including the places where
    /// it would read through a null pointer: a missing node, a missing child,
    /// or a missing optional vector.
    Failed(Fault),
    /// The operation needs SMLib code that a later port step brings.
    NotPorted(SolverOperation),
}

impl From<Fault> for SolverRefusal {
    fn from(fault: Fault) -> Self {
        Self::Failed(fault)
    }
}

/// The node to subdivide next and whether to visit its second child first.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Subdivision {
    /// `rlBestSubdivisionIndex`: the tree whose node splits.
    pub tree: usize,
    /// `rbSwapTraversalOrder`: visit child 2 before child 1.
    pub swap_traversal_order: bool,
}

/// The behaviour that SMLib's derived solvers override.
///
/// Every method defaults to `SmGlobalSolver`'s own implementation, so a
/// solver that overrides nothing is [`BaseHooks`]. A derived solver
/// overrides `local_solve` at least.
pub trait SolverHooks {
    /// `BranchMayContainAnswers`: culls a branch that cannot hold an answer.
    fn branch_may_contain_answers(
        &mut self,
        solver: &GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<bool, SolverRefusal> {
        solver.base_branch_may_contain_answers(branch)
    }

    /// `AreReadyForLocalSolve`: whether the branch is ready for `local_solve`.
    fn are_ready_for_local_solve(
        &self,
        solver: &GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<bool, SolverRefusal> {
        solver.base_are_ready_for_local_solve(branch)
    }

    /// `FindBestSubdivisionNode`: which tree's node to split next.
    fn find_best_subdivision_node(
        &self,
        solver: &GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<Subdivision, SolverRefusal> {
        solver.base_find_best_subdivision_node(branch)
    }

    /// `Subdivide`: splits a branch into two sub-branches.
    fn subdivide(
        &mut self,
        solver: &GlobalSolver<'_>,
        original: &Branch,
    ) -> Result<(Branch, Branch), SolverRefusal> {
        solver.base_subdivide(&*self, original)
    }

    /// `LocalSolve`: solves a ready branch and reports whether it still needs
    /// more subdivision.
    fn local_solve(
        &mut self,
        solver: &mut GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<bool, SolverRefusal> {
        solver.base_local_solve(branch)
    }

    /// `IdenticalSolutions`: whether two solutions are the same answer.
    fn identical_solutions(
        &self,
        objects: &dyn SolutionObjects,
        first: &Solution,
        second: &Solution,
    ) -> (SolutionPair, Preferred) {
        first.compare_for_identity(objects, second)
    }
}

/// The solver that overrides nothing: `SmGlobalSolver` itself, whose local
/// solve records the leaf node pairs.
#[derive(Clone, Copy, Debug, Default)]
pub struct BaseHooks;

impl SolverHooks for BaseHooks {}

/// `SmGlobalSolver`.
pub struct GlobalSolver<'a> {
    nodes: &'a SearchTrees,
    objects: &'a dyn SolutionObjects,
    tops: [Option<NodeId>; MAX_TREES],
    /// `m_eSolverOperation`.
    pub operation: SolverOperation,
    /// `m_eOperationCategory`.
    pub category: OperationCategory,
    /// `m_bProjectedOperation`.
    pub projected: bool,
    /// `m_eSolutionRequested`.
    pub requested: SolutionRequested,
    /// `m_lNumTrees`.
    pub num_trees: usize,
    /// `m_lNumVariables`.
    pub num_variables: usize,
    /// `m_d3dTolerance`.
    pub tolerance: f64,
    /// `m_dBestAnswerSoFarSq`.
    pub best_answer_so_far_sq: f64,
    /// `m_dAtDistance`.
    pub at_distance: f64,
    /// `m_dAtDistanceSq`.
    pub at_distance_sq: f64,
    /// `m_cpOptVectors`: the optional vectors some operations read; empty
    /// stands for SMLib's null pointer.
    pub opt_vectors: Vec<Vector3>,
    /// `m_pSolutions`: the answers found so far.
    pub solutions: SolutionArray,
    /// The owner of the second tree when it is a curve on a surface, with the
    /// interval of each node of that tree, for
    /// `SM_SO_SURFACE_PROJECTED_INTERSECT`; `None` is SMLib's null owner.
    pub surface_curve: Option<SurfaceCurveSpans<'a>>,
}

/// What the second tree's nodes say of a curve on a surface: the owner
/// (`GetOwnerObject`) and each node's Bezier span interval (`m_pData`).
pub struct SurfaceCurveSpans<'a> {
    /// The curve on a surface the tree was built for.
    pub curve: CurveOnSurface<'a>,
    /// The span interval of each node that has one.
    pub spans: Vec<(NodeId, Extent1)>,
}

impl<'a> GlobalSolver<'a> {
    /// The default constructor: no trees, an undefined tolerance.
    pub fn new(nodes: &'a SearchTrees, objects: &'a dyn SolutionObjects) -> Self {
        Self {
            nodes,
            objects,
            tops: [None; MAX_TREES],
            operation: SolverOperation::Intersect,
            category: OperationCategory::Other,
            projected: false,
            requested: SolutionRequested::All,
            num_trees: 0,
            num_variables: 0,
            tolerance: UNDEFINED,
            best_answer_so_far_sq: 0.0,
            at_distance: 0.0,
            at_distance_sq: 0.0,
            opt_vectors: Vec::new(),
            solutions: SolutionArray::new(),
            surface_curve: None,
        }
    }

    /// The constructor for one to three trees, given their top nodes.
    pub fn with_trees(
        nodes: &'a SearchTrees,
        objects: &'a dyn SolutionObjects,
        tops: &[Option<NodeId>],
        num_variables: usize,
    ) -> Self {
        let mut solver = Self::new(nodes, objects);
        for (slot, top) in solver.tops.iter_mut().zip(tops) {
            *slot = *top;
        }
        solver.num_trees = tops.len();
        solver.num_variables = num_variables;
        solver.tolerance = 0.0;
        solver
    }

    /// The constructor for one tree whose top node is `top`.
    pub fn with_one_tree(
        nodes: &'a SearchTrees,
        objects: &'a dyn SolutionObjects,
        top: Option<NodeId>,
        num_variables: usize,
    ) -> Self {
        Self::with_trees(nodes, objects, &[top], num_variables)
    }

    /// The constructor for two trees.
    pub fn with_two_trees(
        nodes: &'a SearchTrees,
        objects: &'a dyn SolutionObjects,
        top1: Option<NodeId>,
        top2: Option<NodeId>,
        num_variables: usize,
    ) -> Self {
        Self::with_trees(nodes, objects, &[top1, top2], num_variables)
    }

    /// `SetNumTrees`: the trees the solver walks; SMLib allows 1 to 3.
    pub fn set_num_trees(&mut self, num_trees: usize) {
        self.num_trees = num_trees;
    }

    /// `SetAtDistance`: sets the distance and its square together.
    pub fn set_at_distance(&mut self, at_distance: f64) {
        self.at_distance = at_distance;
        self.at_distance_sq = at_distance * at_distance;
    }

    /// `GetOptVector`.
    pub fn opt_vector(&self, index: usize) -> Result<Vector3, SolverRefusal> {
        self.opt_vectors
            .get(index)
            .copied()
            .ok_or(SolverRefusal::Failed(Fault::InvalidInput))
    }

    fn node(&self, branch: &Branch, tree: usize) -> Result<&'a TreeNode, SolverRefusal> {
        branch
            .get(tree)
            .copied()
            .flatten()
            .and_then(|id| self.nodes.get(id))
            .ok_or(SolverRefusal::Failed(Fault::Degenerate))
    }

    fn bbox(&self, branch: &Branch, tree: usize) -> Result<Extent3, SolverRefusal> {
        Ok(self.node(branch, tree)?.bbox)
    }

    fn child(&self, node: &TreeNode, which: usize) -> Result<&'a TreeNode, SolverRefusal> {
        node.children
            .and_then(|children| self.nodes.get(children[which]))
            .ok_or(SolverRefusal::Failed(Fault::Degenerate))
    }

    /// `SolveIt`: walks the trees from their top nodes and appends answers to
    /// `solutions`.
    ///
    /// Answers found before a refusal stay in `solutions`, as in SMLib.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_it<H: SolverHooks>(
        &mut self,
        hooks: &mut H,
        operation: SolverOperation,
        requested: SolutionRequested,
        tolerance: f64,
        best_answer_so_far_sq: f64,
        opt_vectors: &[Vector3],
        solutions: &mut SolutionArray,
    ) -> Result<(), SolverRefusal> {
        self.projected = operation.is_projected();
        self.category = operation.category();
        self.operation = operation;
        self.requested = requested;
        self.best_answer_so_far_sq = best_answer_so_far_sq;
        self.opt_vectors = opt_vectors.to_vec();
        self.tolerance = tolerance;

        if self.best_answer_so_far_sq == BIG_DOUBLE
            && matches!(
                operation,
                SolverOperation::Maximize | SolverOperation::DirectedMaximize
            )
        {
            self.best_answer_so_far_sq = 0.0;
        }

        self.solutions = std::mem::take(solutions);
        let mut start: Branch = [None; MAX_TREES];
        for (slot, top) in start.iter_mut().zip(self.tops).take(self.num_trees) {
            *slot = top;
        }
        let result = self.solve_branch(hooks, start);
        *solutions = std::mem::take(&mut self.solutions);
        result
    }

    /// `SolveBranch`: the depth-first search below one branch.
    ///
    /// SMLib recurses; this keeps its own stack and visits the same branches
    /// in the same order, so a deep tree cannot overflow the call stack.
    fn solve_branch<H: SolverHooks>(
        &mut self,
        hooks: &mut H,
        start: Branch,
    ) -> Result<(), SolverRefusal> {
        let mut pending = vec![start];
        while let Some(branch) = pending.pop() {
            if !hooks.branch_may_contain_answers(self, &branch)? {
                continue;
            }
            if hooks.are_ready_for_local_solve(self, &branch)? {
                let needs_more_subdivision = hooks.local_solve(self, &branch)?;
                if !needs_more_subdivision {
                    continue;
                }
            }
            let (first, second) = hooks.subdivide(self, &branch)?;
            pending.push(second);
            pending.push(first);
        }
        Ok(())
    }

    /// `SmGlobalSolver::BranchMayContainAnswers`: whether the bounding boxes,
    /// the operation and the tolerances leave room for an answer in the branch.
    #[allow(clippy::too_many_lines)]
    pub fn base_branch_may_contain_answers(&self, branch: &Branch) -> Result<bool, SolverRefusal> {
        if self.num_trees == 1 {
            return Ok(true);
        }
        let tolerance = self.tolerance;
        match self.operation {
            SolverOperation::Intersect | SolverOperation::IntersectionTest => {
                let mut i = 0;
                while i + 1 < self.num_trees {
                    let mut grown = self.bbox(branch, i)?;
                    grown.expand_absolute(tolerance);
                    for j in i + 1..self.num_trees {
                        if grown.are_disjoint(self.bbox(branch, j)?, 0.0) {
                            return Ok(false);
                        }
                    }
                    i += 1;
                }
                Ok(true)
            }
            SolverOperation::Minimize => {
                let min_sq = minimum_distance_squared(self.bbox(branch, 0)?, self.bbox(branch, 1)?);
                if min_sq > self.best_answer_so_far_sq
                    && sqrt(min_sq) - tolerance > sqrt(self.best_answer_so_far_sq)
                {
                    return Ok(false);
                }
                Ok(true)
            }
            SolverOperation::Maximize => {
                let max_sq = maximum_distance_squared(self.bbox(branch, 0)?, self.bbox(branch, 1)?);
                if max_sq < self.best_answer_so_far_sq
                    && sqrt(max_sq) + tolerance < sqrt(self.best_answer_so_far_sq)
                {
                    return Ok(false);
                }
                Ok(true)
            }
            SolverOperation::AtDistance => {
                let (box1, box2) = (self.bbox(branch, 0)?, self.bbox(branch, 1)?);
                let max_sq = maximum_distance_squared(box1, box2);
                if max_sq < self.at_distance_sq
                    && sqrt(max_sq) - tolerance < sqrt(self.at_distance_sq)
                {
                    return Ok(false);
                }
                let min_sq = minimum_distance_squared(box1, box2);
                if min_sq > self.at_distance_sq
                    && sqrt(min_sq) + tolerance > sqrt(self.at_distance_sq)
                {
                    return Ok(false);
                }
                Ok(true)
            }
            SolverOperation::ProjectedMinimize
            | SolverOperation::ProjectedMaximize
            | SolverOperation::ProjectedIntersect
            | SolverOperation::RotatedProjectedIntersect
            | SolverOperation::PerspectiveIntersect => self.projected_sphere_test(branch),
            SolverOperation::DirectedMinimize
            | SolverOperation::SignedDirectedMinimize
            | SolverOperation::DirectedMaximize => {
                let (v0, v1) = (self.opt_vector(0)?, self.opt_vector(1)?);
                let plane_normal = v1.cross(v0);
                let (center1, radius1, _) = self.bbox(branch, 0)?.sphere_bound();
                let (center2, radius2, _) = self.bbox(branch, 1)?.sphere_bound();
                let projected = point_project_to_plane(center2, center1, plane_normal);
                let sum_radius = radius1 + radius2;
                if projected.distance_squared(center2) > sum_radius * sum_radius + EFF_ZERO_SQRT {
                    return Ok(false);
                }
                let projected = point_project_to_plane(projected, center1, v0);
                match self.operation {
                    SolverOperation::SignedDirectedMinimize => {
                        let vec = projected - center1;
                        let mut dist = vec.length();
                        if vec.dot(v1) < 0.0 {
                            dist = -dist;
                        }
                        let min_proj = dist - radius1 - radius2;
                        if min_proj * min_proj.abs() > self.best_answer_so_far_sq {
                            return Ok(false);
                        }
                    }
                    SolverOperation::DirectedMinimize => {
                        let min_proj = projected.distance(center1) - radius1 - radius2;
                        if min_proj > 0.0 && min_proj * min_proj > self.best_answer_so_far_sq {
                            return Ok(false);
                        }
                    }
                    _ => {
                        let max_proj = projected.distance(center1) + radius1 + radius2;
                        if max_proj * max_proj < self.best_answer_so_far_sq {
                            return Ok(false);
                        }
                    }
                }
                Ok(true)
            }
            SolverOperation::SignedPivotMinimize => {
                let (v0, v1) = (self.opt_vector(0)?, self.opt_vector(1)?);
                let (center1, radius1, _) = self.bbox(branch, 0)?.sphere_bound();
                let (center2, radius2, _) = self.bbox(branch, 1)?.sphere_bound();
                let vec1 = v0.cross(center1 - v1).cross(v0);
                let vec2 = v0.cross(center2 - v1).cross(v0);
                let sum_radius = radius1 + radius2;
                if (vec1.length() - vec2.length()).abs() - tolerance > sum_radius {
                    return Ok(false);
                }
                Ok(true)
            }
            SolverOperation::Normalize | SolverOperation::Find => Ok(true),
            SolverOperation::Rayfire => {
                let v0 = self.opt_vector(0)?;
                let line_point = self.bbox(branch, 0)?.min;
                let (center, radius, _) = self.bbox(branch, 1)?.sphere_bound();
                let param = line_closest_point(line_point, v0, center);
                let closest = line_point + param * v0;
                let dist = closest.distance(center);
                if dist > radius + tolerance {
                    return Ok(false);
                }
                if param + radius < 0.0 {
                    return Ok(false);
                }
                let min_possible = param - radius - tolerance;
                if min_possible > 0.0 && min_possible * min_possible > self.best_answer_so_far_sq {
                    return Ok(false);
                }
                let length = 2.0 * (param.abs() + radius);
                let end = line_point + length * v0;
                let start = line_point - length * v0;
                let mut ray_box = Extent3::singleton(start);
                ray_box.add_point(end);
                ray_box.expand_absolute(tolerance);
                if ray_box.are_disjoint(self.bbox(branch, 1)?, 0.0) {
                    return Ok(false);
                }
                Ok(true)
            }
            SolverOperation::SignedDirectedMinimize3d => {
                let v0 = self.opt_vector(0)?;
                let (center1, radius1, _) = self.bbox(branch, 0)?.sphere_bound();
                let (center2, radius2, _) = self.bbox(branch, 1)?.sphere_bound();
                let param = line_closest_point(center1, v0, center2);
                let closest = center1 + param * v0;
                let dist = closest.distance(center2);
                let sum_radius = radius1 + radius2;
                if dist > sum_radius + tolerance {
                    return Ok(false);
                }
                let vec = closest - center1;
                let mut dist_to = vec.length();
                if vec.dot(v0) < 0.0 {
                    dist_to = -dist_to;
                }
                let min_dist = dist_to - radius1 - radius2 - tolerance;
                if min_dist * min_dist.abs() > self.best_answer_so_far_sq {
                    return Ok(false);
                }
                Ok(true)
            }
            SolverOperation::SurfaceProjectedIntersect => {
                let mut box1 = self.bbox(branch, 0)?;
                box1.expand_absolute(tolerance);
                let box2 = self.bbox(branch, 1)?;
                let min_distance = minimum_distance(box1, box2);
                // Overlapping boxes may hold an answer.
                if min_distance < tolerance {
                    return Ok(true);
                }
                // Otherwise offset the sampled curve along the surface normal
                // by the nearest and the farthest the boxes are apart.
                let max_distance = sqrt(maximum_distance_squared(box1, box2));
                let Some(spans) = &self.surface_curve else {
                    return Ok(false);
                };
                let second = branch[1].ok_or(Fault::Degenerate)?;
                let Some((_, interval)) = spans.spans.iter().find(|(node, _)| *node == second)
                else {
                    return Ok(false);
                };
                let towards = (box1.min + box1.max) / 2.0 - (box2.min + box2.max) / 2.0;
                let mut offsets: Option<Extent3> = None;
                const SAMPLES: usize = 10;
                for i in 0..=SAMPLES {
                    let t = interval.evaluate(i as f64 / SAMPLES as f64);
                    // SMLib ignores both statuses and goes on with the
                    // undefined marker the failed output holds.
                    let point = spans.curve.evaluate_point(t).unwrap_or_default();
                    let mut normal = spans
                        .curve
                        .evaluate_surface_normal(t, false)
                        .map(|(_, normal)| normal)
                        .unwrap_or_default();
                    if normal.dot(towards) < 0.0 {
                        normal *= -1.0;
                    }
                    for distance in [min_distance, max_distance] {
                        let tip = point + distance * normal;
                        match offsets.as_mut() {
                            Some(offsets) => offsets.add_point(tip),
                            None => offsets = Some(Extent3::singleton(tip)),
                        }
                    }
                }
                let mut offsets = offsets.ok_or(Fault::Degenerate)?;
                // Allow for the sampling.
                offsets.expand_absolute(2.0 * tolerance);
                Ok(!offsets.are_disjoint(box1, 0.0))
            }
            SolverOperation::ProjectedTangentThroughPoint => Ok(true),
            // SMLib reports "not handled yet" and then answers TRUE.
            SolverOperation::IntersectWireframe
            | SolverOperation::AngleMinimize
            | SolverOperation::SignedAngleMinimize
            | SolverOperation::ProjectedTangency
            | SolverOperation::SignedDirectedMaximize3d => Ok(true),
        }
    }

    /// The sphere test shared by the projected and perspective operations.
    fn projected_sphere_test(&self, branch: &Branch) -> Result<bool, SolverRefusal> {
        let tolerance = self.tolerance;
        let (center1, radius1, _) = self.bbox(branch, 0)?.sphere_bound();
        let (center2, radius2, _) = self.bbox(branch, 1)?.sphere_bound();
        let mut radius1 = radius1 + tolerance;
        let mut radius2 = radius2;

        let (projected1, projected2);
        match self.operation {
            SolverOperation::RotatedProjectedIntersect => {
                let (v0, v1, v2) = (
                    self.opt_vector(0)?,
                    self.opt_vector(1)?,
                    self.opt_vector(2)?,
                );
                projected1 = rotate_project_point_to_plane(center1, v0, v1, v2);
                projected2 = rotate_project_point_to_plane(center2, v0, v1, v2);
            }
            SolverOperation::PerspectiveIntersect => {
                let (v0, v1, v2) = (
                    self.opt_vector(0)?,
                    self.opt_vector(1)?,
                    self.opt_vector(2)?,
                );
                projected1 = perspective_project_point_to_plane(center1, v0, v1, v2);
                projected2 = perspective_project_point_to_plane(center2, v0, v1, v2);

                // The spheres project to ellipses; find the largest radius of
                // each projection.
                let eye_vec1 = center1 - v2;
                let mut dir1 = eye_vec1.cross(v1);
                let _ = dir1.unitize();
                let mut dir2 = dir1.cross(eye_vec1);
                let _ = dir2.unitize();
                let mut new_radius1_max_sq = 0.0;
                radius1 -= tolerance;
                for ii in 0..4 {
                    let offset = match ii {
                        0 => dir1 * radius1,
                        1 => -dir1 * radius1,
                        2 => dir2 * radius1,
                        _ => -dir2 * radius1,
                    };
                    let projected =
                        perspective_project_vector_to_plane(center1 + offset, center1, v0, v1, v2);
                    let radius_sq = projected.length_squared();
                    if radius_sq > new_radius1_max_sq {
                        new_radius1_max_sq = radius_sq;
                    }
                }

                let eye_vec2 = center2 - v2;
                let mut dir1 = eye_vec2.cross(v1);
                let _ = dir1.unitize();
                let mut dir2 = dir1.cross(eye_vec2);
                let _ = dir2.unitize();
                let mut new_radius2_max_sq = 0.0;
                // SMLib uses the first sphere's radius for the second sphere here.
                for ii in 0..4 {
                    let offset = match ii {
                        0 => dir1 * radius1,
                        1 => -dir1 * radius1,
                        2 => dir2 * radius1,
                        _ => -dir2 * radius1,
                    };
                    let projected =
                        perspective_project_vector_to_plane(center2 + offset, center1, v0, v1, v2);
                    let radius_sq = projected.length_squared();
                    if radius_sq > new_radius2_max_sq {
                        new_radius2_max_sq = radius_sq;
                    }
                }
                radius1 = sqrt(new_radius1_max_sq) + tolerance;
                radius2 = sqrt(new_radius2_max_sq);
            }
            _ => {
                let normal = self.opt_vector(0)?;
                projected2 = point_project_to_plane(center2, center1, normal);
                projected1 = center1;
            }
        }

        match self.operation {
            SolverOperation::ProjectedIntersect
            | SolverOperation::RotatedProjectedIntersect
            | SolverOperation::PerspectiveIntersect => {
                let min_proj = projected1.distance(projected2) - radius1 - radius2;
                if min_proj > 0.0 {
                    return Ok(false);
                }
            }
            SolverOperation::ProjectedMinimize => {
                let min_proj = projected1.distance(projected2) - radius1 - radius2;
                if min_proj > 0.0 && min_proj * min_proj > self.best_answer_so_far_sq {
                    return Ok(false);
                }
            }
            SolverOperation::ProjectedMaximize => {
                let max_proj = projected1.distance(projected2) + radius1 + radius2;
                if max_proj * max_proj < self.best_answer_so_far_sq {
                    return Ok(false);
                }
            }
            _ => {}
        }
        Ok(true)
    }

    /// `SmGlobalSolver::AreReadyForLocalSolve`: every node is a leaf, except
    /// that a Bezier-surface node whose children are not Bezier-surface nodes
    /// counts as ready.
    pub fn base_are_ready_for_local_solve(&self, branch: &Branch) -> Result<bool, SolverRefusal> {
        for tree in 0..self.num_trees {
            let node = self.node(branch, tree)?;
            if node.children.is_some() {
                let first = self.child(node, 0)?;
                if node.aux != AuxData::BezierSurface || first.aux == AuxData::BezierSurface {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    /// `sm_FindLargestNonLeaf`: the tree whose splittable node has the
    /// largest box diagonal, the last of equals, or `None` when none splits.
    fn find_largest_non_leaf(&self, branch: &Branch) -> Result<Option<usize>, SolverRefusal> {
        let mut found = None;
        let mut max_size = 0.0;
        for tree in 0..self.num_trees {
            let node = self.node(branch, tree)?;
            if node.children.is_none() {
                continue;
            }
            if node.aux == AuxData::BezierSurface
                && self.child(node, 0)?.aux != AuxData::BezierSurface
            {
                continue;
            }
            let size = (node.bbox.max - node.bbox.min).length_squared();
            if size >= max_size {
                max_size = size;
                found = Some(tree);
            }
        }
        Ok(found)
    }

    /// `SmGlobalSolver::FindBestSubdivisionNode`: the node to split, and
    /// whether to visit its second child first so that minimising and
    /// maximising searches meet good answers early.
    pub fn base_find_best_subdivision_node(
        &self,
        branch: &Branch,
    ) -> Result<Subdivision, SolverRefusal> {
        let mid = Vector3::new(0.5, 0.5, 0.5);
        let tree = self
            .find_largest_non_leaf(branch)?
            .ok_or(SolverRefusal::Failed(Fault::Degenerate))?;
        let keep = |swap_traversal_order| {
            Ok(Subdivision {
                tree,
                swap_traversal_order,
            })
        };

        if self.operation == SolverOperation::Rayfire {
            let second_tree = self.node(branch, 1)?;
            let sub1 = self.child(second_tree, 0)?.bbox;
            let sub2 = self.child(second_tree, 1)?.bbox;
            let mid1 = sub1.evaluate(mid);
            let mid2 = sub2.evaluate(mid);
            let line_point = self.bbox(branch, 0)?.min;
            let v0 = self.opt_vector(0)?;
            let param1 = line_closest_point(line_point, v0, mid1);
            let param2 = line_closest_point(line_point, v0, mid2);
            return keep(param2 < param1);
        }

        if self.operation == SolverOperation::SignedDirectedMinimize3d {
            let other = 1usize
                .checked_sub(tree)
                .ok_or(SolverRefusal::Failed(Fault::Degenerate))?;
            let other_box = self.bbox(branch, other)?;
            let split = self.node(branch, tree)?;
            let sub1 = self.child(split, 0)?.bbox;
            let sub2 = self.child(split, 1)?.bbox;
            let mid1 = sub1.evaluate(mid);
            let mid2 = sub2.evaluate(mid);
            let line_point = other_box.evaluate(mid);
            let v0 = self.opt_vector(0)?;
            let param1 = line_closest_point(line_point, v0, mid1);
            let param2 = line_closest_point(line_point, v0, mid2);
            return keep(param2 < param1);
        }

        if self.num_trees == 1 {
            return keep(false);
        }

        let minimizing = matches!(
            self.operation,
            SolverOperation::Minimize
                | SolverOperation::DirectedMinimize
                | SolverOperation::SignedDirectedMinimize
                | SolverOperation::ProjectedMinimize
        );
        let maximizing = matches!(
            self.operation,
            SolverOperation::Maximize
                | SolverOperation::DirectedMaximize
                | SolverOperation::ProjectedMaximize
        );
        if minimizing || maximizing {
            let other = 1usize
                .checked_sub(tree)
                .ok_or(SolverRefusal::Failed(Fault::Degenerate))?;
            let other_box = self.bbox(branch, other)?;
            let split = self.node(branch, tree)?;
            let sub1 = self.child(split, 0)?.bbox;
            let sub2 = self.child(split, 1)?.bbox;
            let mid1 = sub1.evaluate(mid);
            let mid2 = sub2.evaluate(mid);
            let mid_other = other_box.evaluate(mid);
            let v1 = (mid1 - mid_other).length_squared();
            let v2 = (mid2 - mid_other).length_squared();
            return keep(if minimizing { v2 < v1 } else { v2 > v1 });
        }
        keep(false)
    }

    /// `SmGlobalSolver::Subdivide`: two sub-branches that differ in the
    /// children of the node `find_best_subdivision_node` picks.
    pub fn base_subdivide<H: SolverHooks + ?Sized>(
        &self,
        hooks: &H,
        original: &Branch,
    ) -> Result<(Branch, Branch), SolverRefusal> {
        let best = hooks.find_best_subdivision_node(self, original)?;
        let mut first = [None; MAX_TREES];
        let mut second = [None; MAX_TREES];
        let used = self.num_trees.min(MAX_TREES);
        first[..used].copy_from_slice(&original[..used]);
        second[..used].copy_from_slice(&original[..used]);
        let split = self.node(original, best.tree)?;
        let children = split
            .children
            .ok_or(SolverRefusal::Failed(Fault::Degenerate))?;
        first[best.tree] = Some(children[0]);
        second[best.tree] = Some(children[1]);
        if best.swap_traversal_order {
            std::mem::swap(&mut first[best.tree], &mut second[best.tree]);
        }
        Ok((first, second))
    }

    /// `SmGlobalSolver::LocalSolve`: with no context the best it can do is
    /// record the leaf nodes where it stopped. Returns whether the branch
    /// needs more subdivision, which it never does.
    pub fn base_local_solve(&mut self, branch: &Branch) -> Result<bool, SolverRefusal> {
        let mut solution = Solution {
            kind: SolutionType::NodesOnly,
            ..Solution::default()
        };
        for (slot, node) in solution.nodes.iter_mut().zip(branch).take(self.num_trees) {
            *slot = *node;
        }
        self.solutions.add(solution);
        Ok(false)
    }

    /// `AddSortedSolution`: adds `solution` to the answers in its place under
    /// `key`, unless it duplicates an answer or is far from the best, then
    /// merges answers that overlap or touch and keeps only the best ones for
    /// single, minimising and maximising searches.
    ///
    /// `combine_contiguous_ranges` false keeps ranges that merely touch apart.
    #[allow(clippy::too_many_lines)]
    pub fn add_sorted_solution<H: SolverHooks + ?Sized>(
        &mut self,
        hooks: &H,
        solution: &Solution,
        key: SortKey,
        combine_contiguous_ranges: bool,
    ) {
        if solution.start.value == BIG_DOUBLE {
            return;
        }
        let value = solution.start.value;
        let this_sq = value * value;
        let mut better_than_best = false;
        let mut far_from_best = false;
        let best_so_far = sqrt(self.best_answer_so_far_sq.abs());

        if self.category == OperationCategory::Signed {
            let sign = if value >= 0.0 { 1.0 } else { -1.0 };
            better_than_best = sign * this_sq < self.best_answer_so_far_sq;
            let best_sign = if self.best_answer_so_far_sq < 0.0 {
                -1.0
            } else {
                1.0
            };
            far_from_best = value > best_sign * best_so_far + self.tolerance;
        }
        if better_than_best {
            self.best_answer_so_far_sq = this_sq;
        } else if far_from_best {
            return;
        }

        match self.operation {
            SolverOperation::Minimize
            | SolverOperation::DirectedMinimize
            | SolverOperation::ProjectedMinimize
            | SolverOperation::Rayfire => {
                better_than_best = this_sq < self.best_answer_so_far_sq;
                far_from_best = value > best_so_far + self.tolerance;
            }
            SolverOperation::Maximize
            | SolverOperation::DirectedMaximize
            | SolverOperation::ProjectedMaximize => {
                better_than_best = this_sq > self.best_answer_so_far_sq;
                far_from_best = value < best_so_far - self.tolerance;
            }
            _ => {}
        }
        if better_than_best {
            self.best_answer_so_far_sq = this_sq;
        } else if far_from_best {
            return;
        }

        let maximizing = self.operation.is_maximizing();
        let mut found_place = false;
        let mut insertion_index = self.solutions.len();
        let mut ii = 0;
        while ii < self.solutions.len() {
            let test = self.solutions[ii];
            let (pair, preferred) = hooks.identical_solutions(self.objects, solution, &test);
            if pair == SolutionPair::Identical {
                match preferred {
                    Preferred::First => self.solutions[ii] = *solution,
                    Preferred::Second => self.solutions[ii] = test,
                    Preferred::Neither => {
                        if solution.start.value.abs() < test.start.value.abs() {
                            self.solutions[ii] = *solution;
                        }
                    }
                }
                return;
            }

            match key {
                SortKey::ByValue => {
                    let by_value = if maximizing {
                        test.start.value < solution.start.value
                    } else {
                        test.start.value > solution.start.value
                    };
                    if by_value {
                        insertion_index = ii;
                        found_place = true;
                    }
                }
                SortKey::ByFirstTwoParameters => {
                    if test.start.parameters[0] > solution.start.parameters[0]
                        && test.start.parameters[1] > solution.start.parameters[1]
                    {
                        insertion_index = ii;
                        found_place = true;
                    }
                }
                SortKey::ByFirstParameter => {
                    if test.start.parameters[0] > solution.start.parameters[0] {
                        insertion_index = ii;
                        found_place = true;
                    }
                    // Single values go in front of ranges with the same start,
                    // which helps closed curves.
                    let scale = (1.0 + test.start.parameters[0]) * EFF_ZERO;
                    if (test.start.parameters[0] - solution.start.parameters[0]).abs() < scale
                        && test.kind == SolutionType::RangeOfValues
                    {
                        insertion_index = ii;
                        found_place = true;
                    }
                }
                SortKey::ByValueWithParameter | SortKey::ByAllParameters => {
                    let mut by_parameters = true;
                    if key == SortKey::ByValueWithParameter
                        && !crate::scalar::is_zero(test.start.value - solution.start.value)
                    {
                        let by_value = if maximizing {
                            test.start.value < solution.start.value
                        } else {
                            test.start.value > solution.start.value
                        };
                        if by_value {
                            insertion_index = ii;
                            found_place = true;
                        }
                        by_parameters = false;
                    }
                    if by_parameters {
                        let mut found_one = false;
                        let variables = test.variables.min(crate::solution::MAX_PARAMETERS);
                        for j in 0..variables {
                            if !found_one {
                                if test.start.parameters[j] > solution.start.parameters[j] {
                                    found_one = true;
                                    insertion_index = ii;
                                    found_place = true;
                                } else if test.start.parameters[j] < solution.start.parameters[j] {
                                    found_one = true;
                                    found_place = false;
                                }
                            }
                            let scale = (1.0 + test.start.parameters[j]) * EFF_ZERO_SQRT;
                            if (test.start.parameters[j] - solution.start.parameters[j]).abs()
                                > scale
                            {
                                break;
                            }
                            if j + 1 < variables {
                                let scale = (1.0 + test.start.parameters[j + 1]) * EFF_ZERO_SQRT;
                                if (test.start.parameters[j + 1] - solution.start.parameters[j + 1])
                                    .abs()
                                    > scale
                                {
                                    if test.start.parameters[j + 1]
                                        > solution.start.parameters[j + 1]
                                    {
                                        insertion_index = ii;
                                        found_place = true;
                                        break;
                                    }
                                    found_place = false;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            if found_place {
                break;
            }
            ii += 1;
        }

        if found_place {
            self.solutions.insert_at(insertion_index, *solution, 1);
        } else {
            self.solutions.add(*solution);
        }

        // Combine solutions that overlap or touch. They are not necessarily
        // sequential, so every pair is tried.
        let mut done = false;
        while !done && self.solutions.len() < 100 && !self.solutions.is_empty() {
            done = true;
            let mut k = 0;
            while k + 1 < self.solutions.len() {
                let mut j = k + 1;
                while j < self.solutions.len() {
                    let mut lower = self.solutions[k];
                    let mut upper = self.solutions[j];
                    let combined = combine_min_max(
                        self.operation,
                        &mut lower,
                        &mut upper,
                        combine_contiguous_ranges,
                    );
                    // SMLib's combine may reshape either input in place.
                    self.solutions[k] = lower;
                    self.solutions[j] = upper;
                    if let Some(combined) = combined {
                        self.solutions[k] = combined;
                        let _ = self.solutions.remove_at(j, 1);
                        done = false;
                        break;
                    }
                    j += 1;
                }
                if !done {
                    break;
                }
                k += 1;
            }
        }

        if self.requested == SolutionRequested::Single {
            self.solutions.set_size(1);
        }

        if matches!(
            self.operation,
            SolverOperation::Minimize
                | SolverOperation::Maximize
                | SolverOperation::Rayfire
                | SolverOperation::SignedDirectedMinimize3d
                | SolverOperation::DirectedMinimize
                | SolverOperation::DirectedMaximize
                | SolverOperation::ProjectedMinimize
                | SolverOperation::ProjectedMaximize
                | SolverOperation::SignedDirectedMinimize
                | SolverOperation::SignedPivotMinimize
        ) && !self.solutions.is_empty()
        {
            let first = self.solutions[0];
            let mut tolerance = self.tolerance;
            if self.operation == SolverOperation::Rayfire {
                let last = self.solutions[self.solutions.len() - 1];
                let first_tolerance = first.object_tolerance(self.objects, 0);
                let last_tolerance = last.object_tolerance(self.objects, 0);
                tolerance = max(self.tolerance, first_tolerance + last_tolerance);
            }
            while self.solutions.len() > 1 {
                let last = self.solutions[self.solutions.len() - 1];
                if (first.start.value - last.start.value).abs() > tolerance {
                    self.solutions.set_size(self.solutions.len() - 1);
                } else {
                    break;
                }
            }
        }
    }
}

/// `sm_CombineMinMax`: combines two solutions into one when they describe the
/// same objects and overlap, touch, or one point lies inside a range.
///
/// Returns the combined solution, or `None` when they stay apart. As in SMLib
/// it may also change `min` or `max` in place: a point inside a range whose
/// ends are far less accurate becomes the start of a second range, and the
/// first range is cut short there.
fn combine_min_max(
    operation: SolverOperation,
    min: &mut Solution,
    max: &mut Solution,
    combine_contiguous_ranges: bool,
) -> Option<Solution> {
    if min.kind == SolutionType::RangeOfValues {
        if max.kind == SolutionType::SingleValue {
            let class = min.classify_solution_end(max, true);
            if class == 3 {
                return None;
            }
            if class == -1 || class == 1 || (class == 0 && operation != SolverOperation::Intersect)
            {
                return Some(*min);
            } else if class == 0 && operation == SolverOperation::Intersect {
                let point_deviation = max.start.value;
                let min_deviation = min.start.value;
                let max_deviation = min.end.value;
                if point_deviation < min_deviation * 100.0
                    || point_deviation < max_deviation * 100.0
                {
                    max.kind = SolutionType::RangeOfValues;
                    max.end = min.end;
                    min.end = max.start;
                    return None;
                }
                return Some(*min);
            }
        } else {
            let class1 = min.classify_solution_end(max, true);
            let class2 = min.classify_solution_end(max, false);
            if class1 == 3 || class2 == 3 {
                return None;
            }
            if class1 > 1 || class2 < -1 {
                return None;
            }
            if !combine_contiguous_ranges
                && ((class1 == 1 && class2 == 2) || (class1 == -2 && class2 == -1))
            {
                return None;
            }
            if class1 == 0 && class2 == 0 {
                return Some(*min);
            }
            let mut combined = *min;
            if class1 < -1 {
                combined.start = max.start;
            }
            if class2 > 1 {
                combined.end = max.end;
            }
            return Some(combined);
        }
    }

    // The remaining case is a single value that falls into the range of `max`.
    if max.kind == SolutionType::RangeOfValues {
        let class = max.classify_solution_end(min, true);
        if class.abs() <= 1 {
            return Some(*max);
        }
        if class == -1 || class == 1 || (class == 0 && operation != SolverOperation::Intersect) {
            return Some(*max);
        } else if class == 0 && operation == SolverOperation::Intersect {
            let point_deviation = min.start.value;
            let min_deviation = max.start.value;
            let max_deviation = max.end.value;
            if point_deviation < min_deviation * 100.0 || point_deviation < max_deviation * 100.0 {
                min.kind = SolutionType::RangeOfValues;
                min.end = max.end;
                max.end = min.start;
                return None;
            }
            return Some(*max);
        }
    }
    None
}
