//! Mesh Booleans and decimation (roadmap 3.4, RECIPE.md sections 4.21 and
//! 4.22): build, save, reopen and STEP.
//!
//! Set `CATALLAX_MESH_PROOF` to a directory to keep the STEP files and their
//! exact volumes for the raw OpenCASCADE check in
//! `tools/step-regressions/test_mesh_boolean.py`.
use catallax_features::{BuildError, Built, ImportFiles, evaluate_with_files};
use catallax_import::sha256_hex;
use catallax_measure::{Accuracy, volume};
use catallax_numbers::{Decider, OperationWork, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::{Value, json};

/// A closed UV sphere as ASCII STL, coordinates written to six decimals.
fn sphere(centre: [f64; 3], radius: f64, rings: usize, sectors: usize) -> Vec<u8> {
    let mut points = vec![[centre[0], centre[1], centre[2] + radius]];
    for i in 1..rings {
        let phi = std::f64::consts::PI * i as f64 / rings as f64;
        for j in 0..sectors {
            let theta = std::f64::consts::TAU * j as f64 / sectors as f64;
            points.push([
                centre[0] + radius * phi.sin() * theta.cos(),
                centre[1] + radius * phi.sin() * theta.sin(),
                centre[2] + radius * phi.cos(),
            ]);
        }
    }
    points.push([centre[0], centre[1], centre[2] - radius]);
    let south = points.len() - 1;
    let ring = |i: usize, j: usize| 1 + (i - 1) * sectors + j % sectors;
    let mut triangles = Vec::new();
    for j in 0..sectors {
        triangles.push([0, ring(1, j), ring(1, j + 1)]);
        triangles.push([south, ring(rings - 1, j + 1), ring(rings - 1, j)]);
    }
    for i in 1..rings - 1 {
        for j in 0..sectors {
            triangles.push([ring(i, j), ring(i + 1, j), ring(i + 1, j + 1)]);
            triangles.push([ring(i, j), ring(i + 1, j + 1), ring(i, j + 1)]);
        }
    }
    stl(&points, &triangles)
}

fn stl(points: &[[f64; 3]], triangles: &[[usize; 3]]) -> Vec<u8> {
    let mut s = String::from("solid part\n");
    for t in triangles {
        s.push_str("facet normal 0 0 0\nouter loop\n");
        for &v in t {
            let [x, y, z] = points[v];
            s.push_str(&format!("vertex {x:.6} {y:.6} {z:.6}\n"));
        }
        s.push_str("endloop\nendfacet\n");
    }
    s.push_str("endsolid part\n");
    s.into_bytes()
}

/// An axis box as twelve outward triangles.
fn cuboid(lo: [f64; 3], hi: [f64; 3]) -> Vec<u8> {
    let points: Vec<[f64; 3]> = (0..8)
        .map(|i| std::array::from_fn(|a| if i >> a & 1 == 0 { lo[a] } else { hi[a] }))
        .collect();
    let faces: [[usize; 4]; 6] = [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ];
    let triangles: Vec<[usize; 3]> = faces
        .iter()
        .flat_map(|f| [[f[0], f[1], f[2]], [f[0], f[2], f[3]]])
        .collect();
    stl(&points, &triangles)
}

/// Two spheres in general position: radius 10 at the origin and radius 9
/// at (7, 1, 0.5).
struct Pair {
    a: Vec<u8>,
    b: Vec<u8>,
}

impl Pair {
    fn new(rings: usize, sectors: usize) -> Self {
        Self {
            a: sphere([0.0; 3], 10.0, rings, sectors),
            b: sphere([7.0, 1.0, 0.5], 9.0, rings, sectors),
        }
    }

    fn files(&self) -> ImportFiles<'_> {
        ImportFiles::from([
            (sha256_hex(&self.a), &self.a[..]),
            (sha256_hex(&self.b), &self.b[..]),
        ])
    }

    fn import(&self, id: &str, bytes: &[u8], pose: Value) -> Value {
        json!({"id": id, "op": "import_mesh", "file": {"name": format!("{id}.stl"), "sha256": sha256_hex(bytes)},
            "pose": pose, "scale": "1"})
    }

    /// The two imports, then `steps`; the last step is the result.
    fn recipe(&self, steps: Vec<Value>) -> Recipe {
        let identity = json!({"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]});
        // The tool is placed by a pose so the Boolean reads placed triangles.
        let moved = json!({"rotation": ["1", "0", "0", "0"], "translation": ["0.5", "0", "0"]});
        let mut all = vec![
            self.import("a", &self.a, identity),
            self.import("b", &self.b, moved),
        ];
        all.extend(steps);
        let result = all.last().expect("a step")["id"].clone();
        Recipe::parse(
            &json!({"format": "catallax-recipe/5", "name": "mesh tools", "result": result, "steps": all})
                .to_string(),
        )
        .unwrap()
    }
}

fn boolean(kind: &str) -> Value {
    json!({"id": kind, "op": "mesh_boolean", "kind": kind, "target": "a", "tool": "b"})
}

/// The application's 100M allowance.
fn build(recipe: &Recipe, files: &ImportFiles) -> (Result<Built, BuildError>, OperationWork) {
    with_operation_budget(100_000_000, || {
        evaluate_with_files(recipe, files, &mut Decider::default())
    })
}

fn exact_volume(built: &Built) -> Rational {
    let v = volume(
        &built.model,
        built.result,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert_eq!(v.lo, v.hi, "a planar polyhedron has an exact volume");
    v.lo
}

fn step_bytes(built: &Built) -> Vec<u8> {
    let mut output = Vec::new();
    catallax_step::write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut Decider::default(),
        &mut output,
    )
    .unwrap();
    output
}

fn keep(name: &str, bytes: &[u8], volume: &Rational) {
    if let Some(directory) = std::env::var_os("CATALLAX_MESH_PROOF") {
        let directory = std::path::Path::new(&directory);
        std::fs::write(directory.join(format!("{name}.step")), bytes).unwrap();
        std::fs::write(
            directory.join(format!("{name}.volume")),
            format!("{:?}\n", volume.to_f64_nearest()),
        )
        .unwrap();
    }
}

/// Builds, checks the work, exports, and reopens from the saved text.
fn built(recipe: &Recipe, files: &ImportFiles) -> (Built, Vec<u8>, OperationWork) {
    let (outcome, work) = build(recipe, files);
    assert!(!work.exhausted, "{work:?}");
    let built = outcome.unwrap();
    let step = step_bytes(&built);
    let reopened = Recipe::parse(&recipe.to_text()).unwrap();
    assert_eq!(&reopened, recipe);
    let (again, rework) = build(&reopened, files);
    let again = again.unwrap();
    assert_eq!(
        step_bytes(&again),
        step,
        "a reopened recipe writes the same STEP bytes"
    );
    assert_eq!(rework.used, work.used);
    println!(
        "{}: {} faces, work {}, STEP {} bytes, sha256 {}",
        recipe.result,
        built.model.faces_of_body(built.result).len(),
        work.used,
        step.len(),
        sha256_hex(&step)
    );
    (built, step, work)
}

#[test]
fn mesh_booleans_build_save_reopen_and_export_with_exact_volumes() {
    let pair = Pair::new(8, 16);
    let files = pair.files();
    let mut volumes = std::collections::BTreeMap::new();
    for kind in ["union", "subtract", "intersect"] {
        let recipe = pair.recipe(vec![boolean(kind)]);
        let (built, step, work) = built(&recipe, &files);
        assert!(work.used <= 20_000_000, "{kind}: {work:?}");
        // The faces keep the imports' names: a Boolean makes no faces of its own.
        let steps: std::collections::BTreeSet<&str> = built
            .model
            .faces_of_body(built.result)
            .iter()
            .flat_map(|&f| built.model.face_origins(f).iter().map(|o| o.step.as_str()))
            .collect();
        assert_eq!(steps, ["a", "b"].into(), "{kind}");
        let v = exact_volume(&built);
        keep(kind, &step, &v);
        volumes.insert(kind, v);
    }
    let alone = |id: &str| {
        let recipe = pair.recipe(vec![]);
        let mut text: Value = serde_json::from_str(&recipe.to_text()).unwrap();
        text["result"] = json!(id);
        let (outcome, _) = build(&Recipe::parse(&text.to_string()).unwrap(), &files);
        exact_volume(&outcome.unwrap())
    };
    let (a, b) = (alone("a"), alone("b"));
    let (union, subtract, intersect) = (
        &volumes["union"],
        &volumes["subtract"],
        &volumes["intersect"],
    );
    // Exact, not approximate: the pieces partition both solids.
    assert_eq!(union + intersect, &a + &b);
    assert_eq!(subtract + intersect, a);
    assert!(intersect.sign() == catallax_numbers::Sign::Positive);
}

#[test]
fn a_finer_mesh_boolean_stays_inside_the_work_target() {
    let pair = Pair::new(16, 32);
    let files = pair.files();
    let (_, _, work) = built(&pair.recipe(vec![boolean("subtract")]), &files);
    assert!(work.used <= 20_000_000, "{work:?}");
}

#[test]
fn decimation_builds_saves_reopens_and_exports_and_chains_after_a_boolean() {
    let pair = Pair::new(16, 32);
    let files = pair.files();
    // 960 triangles down to at most 200.
    let decimated = pair.recipe(vec![
        json!({"id": "fewer", "op": "mesh_decimate", "body": "a", "triangles": 200}),
    ]);
    let (built_once, step, work) = built(&decimated, &files);
    assert!(work.used <= 20_000_000, "{work:?}");
    let faces = built_once.model.faces_of_body(built_once.result).len();
    assert!(faces <= 200, "{faces} faces");
    let v = exact_volume(&built_once);
    // The sphere of radius 10 holds 4188.79 mm^3; collapsing edges of a
    // convex mesh cuts material away, never all of it.
    let v64 = v.to_f64_nearest();
    assert!((3600.0..4188.8).contains(&v64), "{v64}");
    keep("decimated", &step, &v);

    // A decimated Boolean: the union, then fewer triangles.
    let chained = pair.recipe(vec![
        boolean("union"),
        json!({"id": "fewer", "op": "mesh_decimate", "body": "union", "triangles": 300}),
    ]);
    let (built_chain, step, work) = built(&chained, &files);
    assert!(work.used <= 20_000_000, "{work:?}");
    assert!(built_chain.model.faces_of_body(built_chain.result).len() <= 300);
    let v = exact_volume(&built_chain);
    keep("union-decimated", &step, &v);

    // Asking for more triangles than the mesh has keeps it as it is.
    let (whole, _, _) = built(
        &pair.recipe(vec![
            json!({"id": "same", "op": "mesh_decimate", "body": "a", "triangles": 5000}),
        ]),
        &files,
    );
    let recipe = pair.recipe(vec![]);
    let mut text: Value = serde_json::from_str(&recipe.to_text()).unwrap();
    text["result"] = json!("a");
    let (a, _) = build(&Recipe::parse(&text.to_string()).unwrap(), &files);
    assert_eq!(exact_volume(&whole), exact_volume(&a.unwrap()));
}

/// The exact volume and lump count of the union, difference and intersection
/// of the box from 0 to 4 and the box from `tool_lo` to `tool_hi`.
fn box_booleans(tool_lo: [f64; 3], tool_hi: [f64; 3]) -> [(Rational, usize); 3] {
    let pair = Pair {
        a: cuboid([0.0; 3], [4.0; 3]),
        b: cuboid(tool_lo, tool_hi),
    };
    let files = pair.files();
    ["union", "subtract", "intersect"].map(|kind| {
        let mut text: Value =
            serde_json::from_str(&pair.recipe(vec![boolean(kind)]).to_text()).unwrap();
        text["steps"][1]["pose"]["translation"] = json!(["0", "0", "0"]);
        let (built, _, _) = built(&Recipe::parse(&text.to_string()).unwrap(), &files);
        let faces = built.model.faces_of_body(built.result).len();
        let volume = if faces == 0 {
            Rational::zero()
        } else {
            exact_volume(&built)
        };
        (volume, built.model.body(built.result).lumps.len())
    })
}

#[test]
fn touching_edges_nested_and_disjoint_meshes_give_exact_volumes() {
    let q = |n: i64| Rational::from_i64(n);
    // The tool's edges at (y, z) = (1, 1) and (3, 3) cross the target's face
    // x = 4 exactly on that face's diagonal, an edge of the target's mesh.
    let [union, subtract, intersect] = box_booleans([2.0, 1.0, 1.0], [6.0, 3.0, 3.0]);
    assert_eq!(union, (q(72), 1));
    assert_eq!(subtract, (q(56), 1));
    assert_eq!(intersect, (q(8), 1));
    // A tool inside the target leaves a cavity when subtracted.
    let [union, subtract, intersect] = box_booleans([1.0; 3], [2.0; 3]);
    assert_eq!(union, (q(64), 1));
    assert_eq!(subtract, (q(63), 1));
    assert_eq!(intersect, (q(1), 1));
    // Apart: two lumps, an empty intersection.
    let [union, subtract, intersect] = box_booleans([5.0; 3], [6.0; 3]);
    assert_eq!(union, (q(65), 2));
    assert_eq!(subtract, (q(64), 1));
    assert_eq!(intersect, (q(0), 0));
}

#[test]
fn coplanar_contact_and_non_mesh_inputs_refuse_by_name() {
    // Two boxes sharing the planes y = 0, y = 2, z = 0 and z = 2.
    let a = cuboid([0.0; 3], [2.0; 3]);
    let b = cuboid([1.0, 0.0, 0.0], [3.0, 2.0, 2.0]);
    let pair = Pair { a, b };
    let files = pair.files();
    let mut text: Value =
        serde_json::from_str(&pair.recipe(vec![boolean("union")]).to_text()).unwrap();
    text["steps"][1]["pose"]["translation"] = json!(["0", "0", "0"]);
    let (outcome, _) = build(&Recipe::parse(&text.to_string()).unwrap(), &files);
    let BuildError::Geometry { refusal, .. } = outcome.unwrap_err() else {
        panic!("expected a geometric refusal");
    };
    assert_eq!(refusal.code, "mesh-boolean.coplanar", "{}", refusal.detail);

    // A block has no triangles to read.
    let mut text: Value =
        serde_json::from_str(&pair.recipe(vec![boolean("union")]).to_text()).unwrap();
    text["steps"][1] = json!({"id": "b", "op": "block", "size": ["1", "1", "1"],
        "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]}});
    let error = Recipe::parse(&text.to_string()).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("mesh_boolean reads closed triangle meshes"),
        "{error}"
    );

    // Both operations need format 5, and a positive triangle count.
    let mut text: Value =
        serde_json::from_str(&pair.recipe(vec![boolean("union")]).to_text()).unwrap();
    text["format"] = json!("catallax-recipe/4");
    let error = Recipe::parse(&text.to_string()).unwrap_err();
    assert!(error.to_string().contains("catallax-recipe/5"), "{error}");
    let zero = json!({"id": "none", "op": "mesh_decimate", "body": "a", "triangles": 0});
    let mut text: Value =
        serde_json::from_str(&pair.recipe(vec![boolean("union")]).to_text()).unwrap();
    text["steps"][2] = zero;
    text["result"] = json!("none");
    let error = Recipe::parse(&text.to_string()).unwrap_err();
    assert!(
        error.to_string().contains("positive triangle count"),
        "{error}"
    );
}
