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

fn vertical_crossings(mesh: &Mesh, x: f64, y: f64) -> Vec<(f64, bool)> {
    let mut hits = Vec::new();
    for triangle in &mesh.triangles {
        let [a, b, c] = triangle.map(|i| mesh.positions[i].map(f64::from));
        let determinant = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        if determinant.abs() < f64::EPSILON {
            continue;
        }
        let u = ((x - a[0]) * (c[1] - a[1]) - (y - a[1]) * (c[0] - a[0])) / determinant;
        let v = ((b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0])) / determinant;
        if u > 0.0 && v > 0.0 && u + v < 1.0 {
            hits.push((
                a[2] + u * (b[2] - a[2]) + v * (c[2] - a[2]),
                determinant > 0.0,
            ));
        }
    }
    hits.sort_by(|a, b| a.0.total_cmp(&b.0));
    hits
}

#[test]
fn public_lattice_print_preserves_cavity_bore_and_outward_orientation() {
    let recipe = Recipe::parse(include_str!("../../../parts/P6-lattice-bracket.json")).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    let options = Options {
        tolerance: "0.5".parse().unwrap(),
        ..Options::default()
    };
    let mesh = printing::mesh_plan(&plan, &options, &mut Decider::default()).unwrap();
    assert!(mesh.max_deviation <= options.tolerance);
    assert!(mesh.cells <= options.max_cells);
    assert!(mesh.positions.len() <= options.max_vertices);
    let mut edges: BTreeMap<(usize, usize), (usize, i64)> = BTreeMap::new();
    let mut volume6 = 0.0;
    for &[a, b, c] in &mesh.triangles {
        for (x, y) in [(a, b), (b, c), (c, a)] {
            let (key, sense) = if x < y { ((x, y), 1) } else { ((y, x), -1) };
            let entry = edges.entry(key).or_default();
            entry.0 += 1;
            entry.1 += sense;
        }
        let [a, b, c] = [a, b, c].map(|i| mesh.positions[i].map(f64::from));
        volume6 += a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0]);
    }
    assert!(edges.values().all(|value| *value == (2, 0)));
    assert!(volume6 > 0.0);
    // This line misses every radius-0.3 beam and crosses the base cavity.
    let crossings = vertical_crossings(&mesh, 4.43, 4.57);
    assert_eq!(crossings.len(), 4, "{crossings:?}");
    for ((height, positive), (expected_height, expected_positive)) in
        crossings
            .iter()
            .zip([(0.0, false), (1.0, true), (2.0, false), (3.0, true)])
    {
        assert!((height - expected_height).abs() <= 0.5, "{crossings:?}");
        assert_eq!(*positive, expected_positive);
    }
    // An off-axis line strictly inside the radius-0.5 through bore.
    assert!(vertical_crossings(&mesh, 6.03, 5.07).is_empty());
}

#[test]
fn a_regularised_union_does_not_print_its_buried_contact_plane() {
    let mut mesh = mesh(
        r#"{"id":"left","op":"block","size":["1","1","1"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
    {"id":"right","op":"block","size":["1","1","1"],"pose":{"rotation":["1","0","0","0"],"translation":["1","0","0"]}},
    {"id":"both","op":"boolean","kind":"union","target":"left","tool":"right"}"#,
        "both",
    );
    assert_eq!(topology(&mesh), (1, 2));
    for point in &mut mesh.positions {
        point.swap(0, 2);
    }
    let crossings = vertical_crossings(&mesh, 0.43, 0.57);
    assert_eq!(crossings.len(), 2, "{crossings:?}");
    assert!(crossings[0].0 <= 0.5 && crossings[1].0 >= 1.5);
}

#[test]
fn printing_keeps_cell_and_vertex_refusals() {
    let recipe = Recipe::parse(r#"{"format":"catallax-recipe/3","name":"caps","steps":[
        {"id":"box","op":"block","size":["1","1","1"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}
    ],"result":"box"}"#).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    for (max_cells, max_vertices) in [(1, 1_000_000), (2_000_000, 4)] {
        let options = Options {
            tolerance: Rational::one(),
            max_cells,
            max_vertices,
        };
        assert_eq!(
            printing::mesh_plan(&plan, &options, &mut Decider::default())
                .unwrap_err()
                .code,
            "field.print-budget"
        );
    }
}
