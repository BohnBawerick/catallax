#![allow(missing_docs)]
//! Differential gate for port step 2 (`SmGlobalSolver`, `SmSolutionArray`)
//! against saved outputs of the pinned C++ SMLib. See
//! `tests/fixtures/search/README.md` for provenance and regeneration, and
//! `tests/fixtures/search/oracle.cpp` for the token layouts read here.

mod common;

use catallax_trace::extent::{Extent1, Extent3};
use catallax_trace::global_solver::{
    AuxData, BaseHooks, GlobalSolver, OperationCategory, SearchTrees, SolutionRequested,
    SolverOperation, SolverRefusal, TreeNode,
};
use catallax_trace::scalar::ZONE_TOL_3D;
use catallax_trace::solution::{
    NodeId, ObjectId, ObjectKind, Preferred, Solution, SolutionArray, SolutionEnd, SolutionObjects,
    SolutionType, SortKey,
};
use catallax_trace::vector::Vector3;
use common::{Input, Output, SUCCESS, run_gate, status_code};

/// The object table of `oracle.cpp`.
struct Table;

impl SolutionObjects for Table {
    fn kind(&self, id: ObjectId) -> ObjectKind {
        match id.0 {
            0 | 1 => ObjectKind::Curve,
            2 | 3 => ObjectKind::Surface,
            4 => ObjectKind::Parameterless,
            5 => ObjectKind::Edge,
            6 => ObjectKind::Face,
            _ => ObjectKind::Other,
        }
    }

    fn parameter_interval(&self, id: ObjectId, direction: usize) -> Extent1 {
        match (id.0, direction) {
            (0, _) => Extent1::new(0.0, 1.0),
            (1, _) => Extent1::new(-2.0, 3.5),
            (2, 0) => Extent1::new(0.0, 2.0),
            (2, _) => Extent1::new(-1.0, 1.0),
            (3, 0) => Extent1::new(0.5, 4.0),
            (3, _) => Extent1::new(2.0, 3.0),
            (5, _) => Extent1::new(0.25, 0.75),
            (6, 0) => Extent1::new(1.0, 3.0),
            (6, _) => Extent1::new(2.0, 5.0),
            _ => Extent1::default(),
        }
    }

    fn zone_tolerance(&self, id: ObjectId) -> f64 {
        match self.kind(id) {
            ObjectKind::Curve | ObjectKind::Surface | ObjectKind::Other => ZONE_TOL_3D,
            _ => 0.0,
        }
    }
}

fn index(value: i64) -> usize {
    usize::try_from(value).expect("non-negative fixture index")
}

fn optional<T>(value: i64, make: impl Fn(u32) -> T) -> Option<T> {
    u32::try_from(value).ok().map(make)
}

fn read_end(input: &mut Input) -> SolutionEnd {
    let mut end = SolutionEnd {
        value: input.d(),
        ..SolutionEnd::default()
    };
    for slot in &mut end.parameters[..4] {
        *slot = input.d();
    }
    end
}

fn read_solution(input: &mut Input) -> Solution {
    let kind = SolutionType::from_code(index(input.i())).expect("solution type");
    let variables = index(input.i());
    let object_count = index(input.i());
    let mut solution = Solution {
        kind,
        variables,
        object_count,
        ..Solution::default()
    };
    for slot in &mut solution.objects {
        *slot = optional(input.i(), ObjectId);
    }
    solution.start = read_end(input);
    solution.end = read_end(input);
    solution
}

fn object_code(id: Option<ObjectId>) -> i64 {
    id.map_or(-1, |id| i64::from(id.0))
}

fn node_code(id: Option<NodeId>) -> i64 {
    id.map_or(-1, |id| i64::from(id.0))
}

fn write_end(out: &mut Output, end: &SolutionEnd) {
    out.d(end.value);
    assert!(
        end.parameters[4..].iter().all(|p| *p == 0.0),
        "stray parameter"
    );
    for parameter in &end.parameters[..4] {
        out.d(*parameter);
    }
}

fn write_solution(out: &mut Output, solution: &Solution) {
    out.i(solution.kind as i64);
    out.i(solution.variables as i64);
    out.i(solution.object_count as i64);
    for object in solution.objects {
        out.i(object_code(object));
    }
    for node in solution.nodes {
        out.i(node_code(node));
    }
    write_end(out, &solution.start);
    write_end(out, &solution.end);
}

fn write_brief(out: &mut Output, solution: &Solution) {
    out.i(solution.kind as i64);
    out.i(solution.variables as i64);
    out.i(solution.object_count as i64);
    for node in solution.nodes {
        out.i(node_code(node));
    }
    out.d(solution.start.value);
}

fn write_array(out: &mut Output, array: &SolutionArray) {
    out.i(array.len() as i64);
    for solution in array.as_slice() {
        write_solution(out, solution);
    }
}

fn read_vectors(input: &mut Input) -> Vec<Vector3> {
    vec![input.v3(), input.v3(), input.v3()]
}

fn read_box(input: &mut Input) -> Extent3 {
    Extent3::from_corners(input.v3(), input.v3())
}

/// Reads `N` nodes and links them into an arena.
fn read_tree(input: &mut Input) -> SearchTrees {
    let count = index(input.i());
    let mut trees = SearchTrees::new();
    for _ in 0..count {
        let bbox = read_box(input);
        let aux = AuxData::from_code(index(input.i())).expect("aux data type");
        let (first, second) = (input.i(), input.i());
        let children = (first >= 0).then(|| {
            assert!(second >= 0, "a node has both children or neither");
            [
                NodeId(u32::try_from(first).expect("child index")),
                NodeId(u32::try_from(second).expect("child index")),
            ]
        });
        trees.add(TreeNode {
            bbox,
            children,
            aux,
        });
    }
    trees
}

fn read_branch(input: &mut Input) -> [Option<NodeId>; 3] {
    std::array::from_fn(|_| optional(input.i(), NodeId))
}

fn refusal_status(refusal: SolverRefusal) -> i64 {
    match refusal {
        SolverRefusal::Failed(fault) => status_code(fault),
        SolverRefusal::NotPorted(operation) => {
            panic!("the gate does not cover {operation:?}, which SMLib answers by sampling a curve")
        }
    }
}

fn status<T>(out: &mut Output, result: Result<T, SolverRefusal>) -> Option<T> {
    match result {
        Ok(value) => {
            out.i(SUCCESS);
            Some(value)
        }
        Err(refusal) => {
            out.i(refusal_status(refusal));
            None
        }
    }
}

fn operation(code: i64) -> SolverOperation {
    SolverOperation::from_code(index(code)).expect("solver operation")
}

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    match op {
        "solcmp" => {
            let (a, b) = (read_solution(input), read_solution(input));
            out.flag(a.is_same_as(&b));
            out.i(a.classify_solution_end(&b, true));
            out.i(a.classify_solution_end(&b, false));
            out.flag(a.is_equal(&Table, &b));
            let (pair, preferred) = a.compare_for_identity(&Table, &b);
            out.i(pair as i64);
            out.i(match preferred {
                Preferred::Neither => 0,
                Preferred::First => 1,
                Preferred::Second => 2,
            });
        }
        "solidx" => {
            let solution = read_solution(input);
            let (parameter, target, flag) = (index(input.i()), input.i(), input.i());
            let target = optional(target, ObjectId);
            let (found, object, first_parameter) =
                match solution.index_of_object(&Table, target, flag != 0) {
                    Ok((object, parameter)) => (SUCCESS, object, parameter),
                    Err((object, parameter)) => (
                        status_code(catallax_trace::Fault::Degenerate),
                        object,
                        parameter,
                    ),
                };
            out.i(found);
            out.i(object as i64);
            out.i(first_parameter as i64);
            match solution.object_for_parameter_index(&Table, parameter) {
                Ok((object, direction)) => {
                    out.i(SUCCESS);
                    out.i(object as i64);
                    out.i(direction as i64);
                }
                Err(fault) => {
                    out.i(status_code(fault));
                    out.i(0);
                    out.i(0);
                }
            }
            if let Some(interval) = out.status(solution.parameter_interval(&Table, parameter)) {
                out.d(interval.min);
                out.d(interval.max);
            }
            out.d(solution.object_tolerance(&Table, parameter));
        }
        "addsorted" => {
            let operation = operation(input.i());
            let requested = SolutionRequested::from_code(index(input.i())).expect("request");
            let category = OperationCategory::from_code(index(input.i())).expect("category");
            let (tolerance, best) = (input.d(), input.d());
            let variables = index(input.i());
            let key = SortKey::from_code(index(input.i())).expect("sort key");
            let combine = input.i() != 0;
            let trees = SearchTrees::new();
            let mut solver = GlobalSolver::new(&trees, &Table);
            solver.operation = operation;
            solver.requested = requested;
            solver.category = category;
            solver.tolerance = tolerance;
            solver.best_answer_so_far_sq = best;
            solver.num_variables = variables;
            for _ in 0..input.i() {
                let solution = read_solution(input);
                solver.solutions.add(solution);
            }
            for _ in 0..input.i() {
                let solution = read_solution(input);
                solver.add_sorted_solution(&BaseHooks, &solution, key, combine);
                out.i(solver.solutions.len() as i64);
                out.d(solver.best_answer_so_far_sq);
            }
            write_array(out, &solver.solutions);
        }
        "array" => {
            let mut array = SolutionArray::new();
            for _ in 0..input.i() {
                const NAMES: [&str; 13] = [
                    "add", "ins", "rem", "size", "set", "grow", "swap", "dedup", "sorted",
                    "append", "copy", "insarr", "reset",
                ];
                match NAMES[index(input.i())] {
                    "add" => out.i(array.add(read_solution(input)) as i64),
                    "ins" => {
                        let (at, count) = (index(input.i()), index(input.i()));
                        array.insert_at(at, read_solution(input), count);
                    }
                    "rem" => {
                        let (at, count) = (index(input.i()), index(input.i()));
                        array.remove_at(at, count).expect("valid removal");
                    }
                    "size" => array.set_size(index(input.i())),
                    "set" => {
                        let at = index(input.i());
                        array.set_at(at, read_solution(input)).expect("valid index");
                    }
                    "grow" => {
                        let at = index(input.i());
                        array.set_at_grow(at, read_solution(input));
                    }
                    "swap" => {
                        let (a, b) = (index(input.i()), index(input.i()));
                        array.swap(a, b).expect("valid indices");
                    }
                    "dedup" => array.remove_duplicate_solutions(),
                    "sorted" => {
                        let key = SortKey::from_code(index(input.i())).expect("sort key");
                        let solution = read_solution(input);
                        out.i(array
                            .add_sorted_solution(&Table, solution, key)
                            .map_or(99_999_999, |at| at as i64));
                    }
                    word @ ("append" | "copy" | "insarr") => {
                        let at = if word == "insarr" {
                            index(input.i())
                        } else {
                            0
                        };
                        let mut other = SolutionArray::new();
                        for _ in 0..input.i() {
                            other.add(read_solution(input));
                        }
                        match word {
                            "append" => {
                                array.append(&other);
                            }
                            "copy" => array.copy_from(&other),
                            _ => array.insert_array_at(at, &other),
                        }
                    }
                    "reset" => array.clear(),
                    word => panic!("unknown array operation {word}"),
                }
                out.i(array.len() as i64);
            }
            write_array(out, &array);
        }
        "branch" => {
            let operation = operation(input.i());
            let trees_used = index(input.i());
            let (tolerance, best, at_distance) = (input.d(), input.d(), input.d());
            let vectors = read_vectors(input);
            let mut trees = SearchTrees::new();
            let mut branch = [None; 3];
            for slot in &mut branch {
                *slot = Some(trees.add(TreeNode {
                    bbox: read_box(input),
                    children: None,
                    aux: AuxData::None,
                }));
            }
            let mut solver = GlobalSolver::new(&trees, &Table);
            solver.operation = operation;
            solver.set_num_trees(trees_used);
            solver.tolerance = tolerance;
            solver.best_answer_so_far_sq = best;
            solver.set_at_distance(at_distance);
            solver.opt_vectors = vectors;
            out.flag(
                solver
                    .base_branch_may_contain_answers(&branch)
                    .expect("a branch SMLib answers"),
            );
        }
        "tree" => {
            let operation = operation(input.i());
            let trees_used = index(input.i());
            let (tolerance, best) = (input.d(), input.d());
            let vectors = read_vectors(input);
            let trees = read_tree(input);
            let branch = read_branch(input);
            let mut solver = GlobalSolver::new(&trees, &Table);
            solver.operation = operation;
            solver.set_num_trees(trees_used);
            solver.tolerance = tolerance;
            solver.best_answer_so_far_sq = best;
            solver.opt_vectors = vectors;
            out.flag(
                solver
                    .base_are_ready_for_local_solve(&branch)
                    .expect("nodes present"),
            );
            if let Some(best) = status(out, solver.base_find_best_subdivision_node(&branch)) {
                out.i(best.tree as i64);
                out.flag(best.swap_traversal_order);
            }
            if let Some((first, second)) = status(out, solver.base_subdivide(&BaseHooks, &branch)) {
                for tree in 0..trees_used {
                    out.i(node_code(first[tree]));
                    out.i(node_code(second[tree]));
                }
            }
        }
        "solve" => {
            let operation = operation(input.i());
            let requested = SolutionRequested::from_code(index(input.i())).expect("request");
            let trees_used = index(input.i());
            let variables = index(input.i());
            let (tolerance, best) = (input.d(), input.d());
            let vectors = read_vectors(input);
            let mut solutions = SolutionArray::new();
            for _ in 0..input.i() {
                solutions.add(read_solution(input));
            }
            let trees = read_tree(input);
            let tops = read_branch(input);
            let mut solver =
                GlobalSolver::with_trees(&trees, &Table, &tops[..trees_used], variables);
            let result = solver.solve_it(
                &mut BaseHooks,
                operation,
                requested,
                tolerance,
                best,
                &vectors,
                &mut solutions,
            );
            status(out, result);
            out.i(solutions.len() as i64);
            for solution in solutions.as_slice() {
                write_brief(out, solution);
            }
        }
        _ => panic!("unknown operation {op}"),
    }
}

#[test]
fn search_matches_smlib() {
    let report = run_gate(include_str!("fixtures/search/cases.txt"), |_| 4, evaluate);
    report.assert_passed(&["regular", "tangent", "touching", "singular"]);
}
