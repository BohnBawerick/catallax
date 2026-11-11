#![allow(missing_docs)]
//! Behaviour of the step 2 port that the SMLib gate cannot reach: object
//! kinds the oracle cannot build, SMLib's null-pointer paths, which return
//! typed refusals here, and the hooks a derived solver overrides.

use catallax_trace::Fault;
use catallax_trace::extent::{Extent1, Extent3};
use catallax_trace::global_solver::{
    AuxData, BaseHooks, Branch, GlobalSolver, SearchTrees, SolutionRequested, SolverHooks,
    SolverOperation, SolverRefusal, TreeNode,
};
use catallax_trace::solution::{
    NoObjects, NodeId, ObjectId, ObjectKind, Solution, SolutionArray, SolutionObjects,
    SolutionType, SortKey,
};
use catallax_trace::vector::Vector3;

/// One volume, one curve, one vertex and one unknown object.
struct Mixed;

impl SolutionObjects for Mixed {
    fn kind(&self, id: ObjectId) -> ObjectKind {
        match id.0 {
            0 => ObjectKind::Volume,
            1 => ObjectKind::Curve,
            2 => ObjectKind::Parameterless,
            _ => ObjectKind::Other,
        }
    }

    fn parameter_interval(&self, id: ObjectId, direction: usize) -> Extent1 {
        match (id.0, direction) {
            (0, 0) => Extent1::new(0.0, 1.0),
            (0, 1) => Extent1::new(10.0, 20.0),
            (0, _) => Extent1::new(100.0, 300.0),
            _ => Extent1::new(-1.0, 1.0),
        }
    }

    fn zone_tolerance(&self, _id: ObjectId) -> f64 {
        0.5
    }
}

fn solution(variables: usize, objects: &[Option<u32>]) -> Solution {
    let mut solution = Solution {
        kind: SolutionType::SingleValue,
        variables,
        object_count: objects.len(),
        ..Solution::default()
    };
    for (slot, object) in solution.objects.iter_mut().zip(objects) {
        *slot = object.map(ObjectId);
    }
    solution
}

fn unit_box(x: f64) -> Extent3 {
    Extent3::from_corners(Vector3::new(x, 0.0, 0.0), Vector3::new(x + 1.0, 1.0, 1.0))
}

fn leaf(trees: &mut SearchTrees, x: f64) -> NodeId {
    trees.add(TreeNode {
        bbox: unit_box(x),
        children: None,
        aux: AuxData::None,
    })
}

fn split(trees: &mut SearchTrees, x: f64, children: [NodeId; 2]) -> NodeId {
    trees.add(TreeNode {
        bbox: Extent3::from_corners(Vector3::new(x, 0.0, 0.0), Vector3::new(x + 2.0, 1.0, 1.0)),
        children: Some(children),
        aux: AuxData::None,
    })
}

#[test]
fn a_volume_has_three_parameter_directions() {
    let solution = solution(3, &[Some(0)]);
    let object = solution.object_for_parameter_index(&Mixed, 2);
    assert_eq!(object, Ok((0, 2)));
    assert_eq!(
        solution.parameter_interval(&Mixed, 2),
        Ok(Extent1::new(100.0, 300.0))
    );
    assert_eq!(
        solution.parameter_interval(&Mixed, 1),
        Ok(Extent1::new(10.0, 20.0))
    );
    // `GetIndex` counts curves, edges, faces and surfaces, not volumes.
    assert_eq!(
        solution.index_of_object(&Mixed, Some(ObjectId(9)), true),
        Err((1, 0))
    );
}

#[test]
fn a_missing_object_is_refused_where_smlib_would_read_a_null_pointer() {
    // The first object is null: SMLib returns SM_ERR before using it.
    let first_missing = solution(2, &[None, Some(1)]);
    assert_eq!(
        first_missing.object_for_parameter_index(&Mixed, 0),
        Err(Fault::Degenerate)
    );
    assert_eq!(
        first_missing.parameter_interval(&Mixed, 0),
        Err(Fault::Degenerate)
    );
    // The second object is null: SMLib would call `IsKindOf` on it.
    let second_missing = solution(2, &[Some(1), None]);
    assert_eq!(
        second_missing.parameter_interval(&Mixed, 1),
        Err(Fault::Degenerate)
    );
    // No object at all, or one with no parameters, has no interval.
    assert!(solution(1, &[]).parameter_interval(&Mixed, 0).is_err());
    assert!(
        solution(1, &[Some(2)])
            .parameter_interval(&Mixed, 0)
            .is_err()
    );
    assert!(
        solution(1, &[Some(3)])
            .parameter_interval(&Mixed, 0)
            .is_err()
    );
}

#[test]
fn a_solution_counts_at_most_five_objects_and_ten_parameters() {
    let mut wild = solution(40, &[Some(1)]);
    wild.object_count = 99;
    let other = wild;
    assert!(wild.is_same_as(&other));
    assert!(wild.is_equal(&NoObjects, &other));
    let (pair, _) = wild.compare_for_identity(&NoObjects, &other);
    assert_eq!(pair as i64, 0);
    assert_eq!(wild.object_tolerance(&Mixed, 7), 0.0);
}

#[test]
fn an_array_grown_by_set_size_holds_zeroed_solutions() {
    let mut array = SolutionArray::new();
    array.set_size(2);
    assert_eq!(array[1].kind, SolutionType::NodesOnly);
    assert_eq!(Solution::default().kind, SolutionType::Unknown);
    array.insert_at(5, solution(1, &[Some(1)]), 2);
    assert_eq!(array.len(), 7);
    assert_eq!(array[3].kind, SolutionType::NodesOnly);
    assert_eq!(array[6].kind, SolutionType::SingleValue);
    assert!(array.remove_at(6, 2).is_err());
    assert!(array.swap(0, 7).is_err());
    assert!(array.set_at(7, Solution::default()).is_err());
}

#[test]
fn surface_projected_intersect_without_a_curve_on_a_surface_answers_false_apart() {
    // SMLib answers true for overlapping boxes, and false when the boxes are
    // apart and the second tree has no owning curve on a surface. The sampled
    // branch is gated in `curve_on_surface.rs`.
    let mut trees = SearchTrees::new();
    let near = leaf(&mut trees, 0.5);
    let far = leaf(&mut trees, 9.0);
    let mut solver = GlobalSolver::new(&trees, &NoObjects);
    solver.operation = SolverOperation::SurfaceProjectedIntersect;
    solver.set_num_trees(2);
    solver.tolerance = 1e-3;
    let overlapping = [Some(near), Some(near), None];
    assert_eq!(
        solver.base_branch_may_contain_answers(&overlapping),
        Ok(true)
    );
    let apart = [Some(near), Some(far), None];
    assert_eq!(solver.base_branch_may_contain_answers(&apart), Ok(false));
}

#[test]
fn missing_nodes_children_and_vectors_are_refusals() {
    let mut trees = SearchTrees::new();
    let a = leaf(&mut trees, 0.0);
    let b = leaf(&mut trees, 0.5);
    let mut solver = GlobalSolver::new(&trees, &NoObjects);
    solver.set_num_trees(2);
    solver.tolerance = 1e-3;

    solver.operation = SolverOperation::Intersect;
    assert_eq!(
        solver.base_branch_may_contain_answers(&[Some(a), None, None]),
        Err(SolverRefusal::Failed(Fault::Degenerate))
    );

    solver.operation = SolverOperation::Rayfire;
    assert_eq!(
        solver.base_branch_may_contain_answers(&[Some(a), Some(b), None]),
        Err(SolverRefusal::Failed(Fault::InvalidInput))
    );
    solver.opt_vectors = vec![Vector3::new(1.0, 0.0, 0.0)];
    assert!(
        solver
            .base_branch_may_contain_answers(&[Some(a), Some(b), None])
            .is_ok()
    );

    // Two leaves: nothing to subdivide.
    assert_eq!(
        solver.base_find_best_subdivision_node(&[Some(a), Some(b), None]),
        Err(SolverRefusal::Failed(Fault::Degenerate))
    );

    // The ray-fire order looks at the second tree's children, which a leaf lacks.
    let c = leaf(&mut trees, 2.0);
    let d = leaf(&mut trees, 3.0);
    let root = split(&mut trees, 2.0, [c, d]);
    let mut solver = GlobalSolver::new(&trees, &NoObjects);
    solver.operation = SolverOperation::Rayfire;
    solver.set_num_trees(2);
    solver.opt_vectors = vec![Vector3::new(1.0, 0.0, 0.0)];
    assert_eq!(
        solver.base_find_best_subdivision_node(&[Some(root), Some(a), None]),
        Err(SolverRefusal::Failed(Fault::Degenerate))
    );
    assert!(
        solver
            .base_find_best_subdivision_node(&[Some(a), Some(root), None])
            .is_ok()
    );
}

/// Records the order in which leaf pairs reach the local solver.
#[derive(Default)]
struct Recorder {
    visited: Vec<Branch>,
}

impl SolverHooks for Recorder {
    fn local_solve(
        &mut self,
        _solver: &mut GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<bool, SolverRefusal> {
        self.visited.push(*branch);
        Ok(false)
    }
}

#[test]
fn traversal_visits_the_first_child_before_the_second() {
    let mut trees = SearchTrees::new();
    let l0 = leaf(&mut trees, 0.0);
    let l1 = leaf(&mut trees, 1.0);
    let a = split(&mut trees, 0.0, [l0, l1]);
    let m0 = leaf(&mut trees, 0.0);
    let m1 = leaf(&mut trees, 1.0);
    let b = split(&mut trees, 0.0, [m0, m1]);
    let mut solver = GlobalSolver::with_two_trees(&trees, &NoObjects, Some(a), Some(b), 0);
    let mut recorder = Recorder::default();
    let mut answers = SolutionArray::new();
    solver
        .solve_it(
            &mut recorder,
            SolverOperation::Intersect,
            SolutionRequested::All,
            0.0,
            0.0,
            &[],
            &mut answers,
        )
        .expect("a full search");
    // Equal boxes: the later tree splits first (`>=`), so the second tree's
    // first child pairs with both children of the first tree before its
    // second child does. Boxes that only touch still count at tolerance 0.
    let leaves: Vec<_> = recorder
        .visited
        .iter()
        .map(|branch| (branch[0], branch[1]))
        .collect();
    assert_eq!(
        leaves,
        [
            (Some(l0), Some(m0)),
            (Some(l1), Some(m0)),
            (Some(l0), Some(m1)),
            (Some(l1), Some(m1)),
        ]
    );
}

/// A derived solver: adds one sorted answer per leaf pair and stops early.
struct Collect {
    stop_after: usize,
    calls: usize,
}

impl SolverHooks for Collect {
    fn local_solve(
        &mut self,
        solver: &mut GlobalSolver<'_>,
        branch: &Branch,
    ) -> Result<bool, SolverRefusal> {
        self.calls += 1;
        let first = branch[0].expect("first node").0;
        let mut found = Solution {
            kind: SolutionType::SingleValue,
            variables: 1,
            object_count: 0,
            ..Solution::default()
        };
        found.start.value = f64::from(first) * 1e-9;
        found.start.parameters[0] = f64::from(first) / 10.0;
        solver.add_sorted_solution(&*self, &found, SortKey::ByAllParameters, true);
        if self.calls >= self.stop_after {
            return Err(SolverRefusal::Failed(Fault::Degenerate));
        }
        Ok(false)
    }
}

#[test]
fn answers_found_before_a_refusal_stay_in_the_output_array() {
    let mut trees = SearchTrees::new();
    let l0 = leaf(&mut trees, 0.0);
    let l1 = leaf(&mut trees, 0.5);
    let a = split(&mut trees, 0.0, [l0, l1]);
    let b = leaf(&mut trees, 0.25);
    let mut answers = SolutionArray::new();
    let mut solver = GlobalSolver::with_two_trees(&trees, &NoObjects, Some(a), Some(b), 1);
    let mut hooks = Collect {
        stop_after: 2,
        calls: 0,
    };
    let result = solver.solve_it(
        &mut hooks,
        SolverOperation::Intersect,
        SolutionRequested::All,
        1e-3,
        0.0,
        &[],
        &mut answers,
    );
    assert_eq!(result, Err(SolverRefusal::Failed(Fault::Degenerate)));
    assert_eq!(hooks.calls, 2);
    assert_eq!(answers.len(), 2);
    assert!(answers[0].start.parameters[0] < answers[1].start.parameters[0]);
}

#[test]
fn the_same_search_gives_the_same_answers_every_time() {
    let run = || {
        let mut trees = SearchTrees::new();
        let l0 = leaf(&mut trees, 0.0);
        let l1 = leaf(&mut trees, 0.5);
        let a = split(&mut trees, 0.0, [l0, l1]);
        let b = leaf(&mut trees, 0.25);
        let mut answers = SolutionArray::new();
        let mut solver = GlobalSolver::with_two_trees(&trees, &NoObjects, Some(a), Some(b), 1);
        solver
            .solve_it(
                &mut Collect {
                    stop_after: usize::MAX,
                    calls: 0,
                },
                SolverOperation::Intersect,
                SolutionRequested::All,
                1e-3,
                0.0,
                &[],
                &mut answers,
            )
            .expect("a full search");
        answers
    };
    assert_eq!(run(), run());
}

#[test]
fn the_default_solver_records_leaf_pairs_as_node_solutions() {
    let mut trees = SearchTrees::new();
    let l0 = leaf(&mut trees, 0.0);
    let l1 = leaf(&mut trees, 5.0);
    let a = split(&mut trees, 0.0, [l0, l1]);
    let b = leaf(&mut trees, 0.5);
    let mut answers = SolutionArray::new();
    let mut solver = GlobalSolver::with_two_trees(&trees, &NoObjects, Some(a), Some(b), 0);
    solver
        .solve_it(
            &mut BaseHooks,
            SolverOperation::Intersect,
            SolutionRequested::Nodes,
            1e-3,
            0.0,
            &[],
            &mut answers,
        )
        .expect("a full search");
    assert_eq!(answers.len(), 1);
    assert_eq!(answers[0].kind, SolutionType::NodesOnly);
    assert_eq!(answers[0].nodes[0], Some(l0));
    assert_eq!(answers[0].nodes[1], Some(b));
}
