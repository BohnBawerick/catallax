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

/// Reference serializers: one small write per field, as before chunking.
fn reference_stl(mesh: &Mesh) -> Vec<u8> {
    let mut out = vec![0u8; 80];
    let label = b"catallax field mesh; millimetres";
    out[..label.len()].copy_from_slice(label);
    out.extend((mesh.triangles.len() as u32).to_le_bytes());
    for triangle in &mesh.triangles {
        let p = triangle.map(|i| mesh.positions[i].map(f64::from));
        let a: [f64; 3] = std::array::from_fn(|i| p[1][i] - p[0][i]);
        let b: [f64; 3] = std::array::from_fn(|i| p[2][i] - p[0][i]);
        let normal: [f64; 3] = std::array::from_fn(|i| {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            a[j] * b[k] - a[k] * b[j]
        });
        let length = normal.iter().map(|n| n * n).sum::<f64>().sqrt();
        for component in normal {
            out.extend(((component / length) as f32).to_le_bytes());
        }
        for &vertex in triangle {
            for coordinate in mesh.positions[vertex] {
                out.extend(coordinate.to_le_bytes());
            }
        }
        out.extend([0u8; 2]);
    }
    out
}
fn reference_obj(mesh: &Mesh) -> Vec<u8> {
    use std::io::Write;
    let mut out = Vec::new();
    writeln!(out, "# catallax field print mesh; millimetres").unwrap();
    writeln!(out, "# proven_two_sided_deviation {}", mesh.max_deviation).unwrap();
    for p in &mesh.positions {
        writeln!(out, "v {} {} {}", p[0] as f64, p[1] as f64, p[2] as f64).unwrap();
    }
    for t in &mesh.triangles {
        writeln!(out, "f {} {} {}", t[0] + 1, t[1] + 1, t[2] + 1).unwrap();
    }
    out
}
/// Large enough that output spans several internal chunks.
fn sphere_mesh() -> Mesh {
    let recipe = Recipe::parse(r#"{"format":"catallax-recipe/3","name":"print bytes","steps":[{"id":"ball","op":"sphere","radius":"10","pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"ball"}"#).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    printing::mesh_plan(
        &plan,
        &Options {
            tolerance: Rational::from_ratio(1, 2).unwrap(),
            ..Options::default()
        },
        &mut Decider::default(),
    )
    .unwrap()
}
#[test]
fn stl_and_obj_bytes_match_the_unchunked_reference() {
    for mesh in [
        mesh(
            r#"{"id":"box","op":"block","size":["1","1","1"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}"#,
            "box",
        ),
        sphere_mesh(),
    ] {
        let (mut stl, mut obj) = (Vec::new(), Vec::new());
        mesh.write_stl(&mut stl).unwrap();
        mesh.write_obj(&mut obj).unwrap();
        assert_eq!(stl, reference_stl(&mesh));
        assert_eq!(obj, reference_obj(&mesh));
    }
    assert!(reference_stl(&sphere_mesh()).len() > 3 << 16);
}
#[test]
fn small_sinks_see_few_writes() {
    struct Counting(usize, usize);
    impl std::io::Write for Counting {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0 += 1;
            self.1 += bytes.len();
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mesh = sphere_mesh();
    let mut stl = Counting(0, 0);
    mesh.write_stl(&mut stl).unwrap();
    assert_eq!(stl.1, reference_stl(&mesh).len());
    assert!(stl.0 <= stl.1 / (1 << 16) + 2);
}
#[test]
fn a_failing_sink_reports_the_error_for_both_formats() {
    struct Failing(usize);
    impl std::io::Write for Failing {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            if self.0 < bytes.len() {
                return Err(std::io::Error::other("disk full"));
            }
            self.0 -= bytes.len();
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mesh = sphere_mesh();
    for budget in [0, 100, 70_000] {
        let stl = mesh.write_stl(&mut Failing(budget)).unwrap_err();
        assert_eq!(stl.to_string(), "disk full");
        let obj = mesh.write_obj(&mut Failing(budget)).unwrap_err();
        assert_eq!(obj.to_string(), "disk full");
    }
    mesh.write_stl(&mut Failing(usize::MAX)).unwrap();
    mesh.write_obj(&mut Failing(usize::MAX)).unwrap();
}
#[test]
fn the_artifact_path_writes_identical_bytes_to_memory() {
    let mesh = sphere_mesh();
    let dir = std::env::temp_dir().join(format!("catallax-print-bytes-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    mesh.write_stl(&mut std::fs::File::create(dir.join("print.stl")).unwrap())
        .unwrap();
    mesh.write_obj(&mut std::fs::File::create(dir.join("display.obj")).unwrap())
        .unwrap();
    assert_eq!(
        std::fs::read(dir.join("print.stl")).unwrap(),
        reference_stl(&mesh)
    );
    assert_eq!(
        std::fs::read(dir.join("display.obj")).unwrap(),
        reference_obj(&mesh)
    );
    let full = std::fs::File::open("/dev/full");
    if let Ok(full) = full {
        drop(full);
        let mut device = std::fs::OpenOptions::new()
            .write(true)
            .open("/dev/full")
            .unwrap();
        assert!(mesh.write_stl(&mut device).is_err());
        assert!(mesh.write_obj(&mut device).is_err());
    }
}
