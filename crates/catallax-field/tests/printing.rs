//! Print topology and complete boundary-coverage regressions.
use catallax_field::{
    Plan,
    printing::{self, Mesh, Options},
};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::Recipe;
use std::collections::{BTreeMap, BTreeSet};

fn mesh(steps: &str, result: &str) -> Mesh {
    let recipe = Recipe::parse(&format!(r#"{{"format":"catallax-recipe/3","name":"print regression","steps":[{steps}],"result":"{result}"}}"#)).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    printing::mesh_plan(
        &plan,
        &Options {
            tolerance: Rational::one(),
            ..Options::default()
        },
        &mut Decider::default(),
    )
    .unwrap()
}
fn topology(mesh: &Mesh) -> (usize, isize) {
    let mut edges = BTreeSet::new();
    let mut adjacency: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for &[a, b, c] in &mesh.triangles {
        for (x, y) in [(a, b), (b, c), (c, a)] {
            edges.insert(if x < y { (x, y) } else { (y, x) });
            adjacency.entry(x).or_default().push(y);
            adjacency.entry(y).or_default().push(x);
        }
    }
    let mut seen = BTreeSet::new();
    let mut components = 0;
    for &start in adjacency.keys() {
        if seen.contains(&start) {
            continue;
        }
        components += 1;
        let mut pending = vec![start];
        while let Some(v) = pending.pop() {
            if seen.insert(v) {
                pending.extend(&adjacency[&v]);
            }
        }
    }
    (
        components,
        mesh.positions.len() as isize - edges.len() as isize + mesh.triangles.len() as isize,
    )
}
#[test]
fn a_solid_cube_has_one_closed_genus_zero_surface() {
    let mesh = mesh(
        r#"{"id":"box","op":"block","size":["1","1","1"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}"#,
        "box",
    );
    assert_eq!(topology(&mesh), (1, 2));
    assert!(mesh.max_deviation <= Rational::one());
}
#[test]
fn interval_coverage_does_not_drop_a_component_between_initial_grid_corners() {
    let mesh = mesh(
        r#"{"id":"box","op":"block","size":["1","1","1"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
    {"id":"tiny","op":"sphere","radius":"0.03125","pose":{"rotation":["1","0","0","0"],"translation":["3.37","0.41","0.43"]}},
    {"id":"both","op":"boolean","kind":"union","target":"box","tool":"tiny"}"#,
        "both",
    );
    assert_eq!(topology(&mesh), (2, 4));
}
#[test]
fn a_through_hole_preserves_the_handle_in_the_printed_surface() {
    let mesh = mesh(
        r#"{"id":"box","op":"block","size":["2","2","2"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
    {"id":"bore","op":"cylinder","radius":"0.4","height":"4","pose":{"rotation":["1","0","0","0"],"translation":["1","1","-1"]}},
    {"id":"part","op":"boolean","kind":"subtract","target":"box","tool":"bore"}"#,
        "part",
    );
    assert_eq!(topology(&mesh), (1, 0));
}
