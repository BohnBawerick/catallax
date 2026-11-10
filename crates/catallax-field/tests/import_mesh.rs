//! Closed mesh import in the field, and exact volume agreement with the B-rep,
//! RECIPE.md section 4.19 and PLAN M22.

use std::collections::BTreeMap;

use catallax_field::mass::{self, Accuracy, Options};
use catallax_field::{Location, Plan, Value, Workspace};
use catallax_import::sha256_hex;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
use serde_json::json;

type Tri = [usize; 3];

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

fn cuboid(lo: [i64; 3], hi: [i64; 3], first: usize, inward: bool) -> (Vec<[String; 3]>, Vec<Tri>) {
    let corners = (0..8)
        .map(|i| std::array::from_fn(|a| (if i >> a & 1 == 0 { lo[a] } else { hi[a] }).to_string()))
        .collect();
    let faces: [[usize; 4]; 6] = [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ];
    let triangles = faces
        .iter()
        .flat_map(|f| [[f[0], f[1], f[2]], [f[0], f[2], f[3]]])
        .map(|t: Tri| {
            let t = t.map(|v| v + first);
            if inward { [t[0], t[2], t[1]] } else { t }
        })
        .collect();
    (corners, triangles)
}

fn ascii(points: &[[String; 3]], triangles: &[Tri]) -> Vec<u8> {
    let mut s = String::from("solid part\n");
    for t in triangles {
        s.push_str("facet normal 0 0 0\nouter loop\n");
        for &v in t {
            let [x, y, z] = &points[v];
            s.push_str(&format!("vertex {x} {y} {z}\n"));
        }
        s.push_str("endloop\nendfacet\n");
    }
    s.push_str("endsolid part\n");
    s.into_bytes()
}

fn binary(points: &[[f32; 3]], triangles: &[Tri]) -> Vec<u8> {
    let mut bytes = vec![b' '; 80];
    bytes.extend_from_slice(&(triangles.len() as u32).to_le_bytes());
    for t in triangles {
        bytes.extend_from_slice(&[0u8; 12]);
        for &v in t {
            for x in points[v] {
                bytes.extend_from_slice(&x.to_le_bytes());
            }
        }
        bytes.extend_from_slice(&[0u8; 2]);
    }
    bytes
}

/// A latitude-longitude sphere of radius `r`, binary32 coordinates, outward.
fn sphere(r: f64, bands: usize, sectors: usize) -> (Vec<[f32; 3]>, Vec<Tri>) {
    let mut points = vec![[0.0, 0.0, r as f32]];
    for i in 1..bands {
        let theta = std::f64::consts::PI * i as f64 / bands as f64;
        for j in 0..sectors {
            let phi = 2.0 * std::f64::consts::PI * j as f64 / sectors as f64;
            points.push([
                (r * theta.sin() * phi.cos()) as f32,
                (r * theta.sin() * phi.sin()) as f32,
                (r * theta.cos()) as f32,
            ]);
        }
    }
    points.push([0.0, 0.0, -r as f32]);
    let ring = |i: usize, j: usize| 1 + (i - 1) * sectors + j % sectors;
    let bottom = points.len() - 1;
    let mut triangles = Vec::new();
    for j in 0..sectors {
        triangles.push([0, ring(1, j), ring(1, j + 1)]);
        triangles.push([bottom, ring(bands - 1, j + 1), ring(bands - 1, j)]);
    }
    for i in 1..bands - 1 {
        for j in 0..sectors {
            let (a, b, c, d) = (
                ring(i, j),
                ring(i, j + 1),
                ring(i + 1, j),
                ring(i + 1, j + 1),
            );
            triangles.push([a, c, d]);
            triangles.push([a, d, b]);
        }
    }
    (points, triangles)
}

fn recipe(hash: &str, scale: &str, rotation: [&str; 4], translation: [&str; 3]) -> Recipe {
    Recipe::parse(
        &json!({
            "format": "catallax-recipe/4",
            "name": "mesh import",
            "steps": [{
                "id": "scan", "op": "import_mesh",
                "file": {"name": "part.stl", "sha256": hash},
                "pose": {"rotation": rotation, "translation": translation},
                "scale": scale
            }],
            "result": "scan"
        })
        .to_string(),
    )
    .unwrap()
}

struct Both {
    plan: Plan,
    brep: catallax_features::Built,
}

fn both(bytes: &[u8], scale: &str, rotation: [&str; 4], translation: [&str; 3]) -> Both {
    let hash = sha256_hex(bytes);
    let saved = recipe(&hash, scale, rotation, translation);
    let files = BTreeMap::from([(hash, bytes)]);
    Both {
        plan: Plan::compile_with_files(&saved, &files, &mut Decider::default()).unwrap(),
        brep: catallax_features::evaluate_with_files(&saved, &files, &mut Decider::default())
            .unwrap(),
    }
}

fn field_refusal(bytes: &[u8]) -> String {
    let hash = sha256_hex(bytes);
    let files = BTreeMap::from([(hash.clone(), bytes)]);
    Plan::compile_with_files(
        &recipe(&hash, "1", ["1", "0", "0", "0"], ["0", "0", "0"]),
        &files,
        &mut Decider::default(),
    )
    .unwrap_err()
    .code
    .to_string()
}

fn field_properties(plan: &Plan) -> mass::Properties {
    mass::properties(
        plan,
        &Accuracy::default(),
        &Options::default(),
        &mut Decider::default(),
    )
    .unwrap()
}

fn brep_volume(built: &catallax_features::Built) -> Rational {
    let v = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert_eq!(v.lo, v.hi, "a planar polyhedron has an exact B-rep volume");
    v.lo
}

fn exact_field_volume(plan: &Plan) -> Rational {
    let p = field_properties(plan);
    assert_eq!(
        p.volume.lo, p.volume.hi,
        "a lone mesh has an exact field volume"
    );
    p.volume.lo
}

fn at(plan: &Plan, x: &str, y: &str, z: &str) -> catallax_field::Query {
    plan.proven_checked(
        &RPoint3::new(q(x), q(y), q(z)),
        &mut Workspace::default(),
        &mut Decider::default(),
    )
    .unwrap()
}

fn cube_with_cavity() -> Vec<u8> {
    let (mut points, mut triangles) = cuboid([0, 0, 0], [6, 6, 6], 0, false);
    let (inner, inner_triangles) = cuboid([2, 2, 2], [4, 4, 4], 8, true);
    points.extend(inner);
    triangles.extend(inner_triangles);
    ascii(&points, &triangles)
}

#[test]
fn the_field_decides_inside_by_winding_number_with_true_distance() {
    let built = both(
        &cube_with_cavity(),
        "1",
        ["1", "0", "0", "0"],
        ["0", "0", "0"],
    );
    let plan = &built.plan;
    let cases = [
        (["1", "3", "3"], Location::Inside, "-1"),
        (["3", "3", "3"], Location::Outside, "1"),
        (["-2", "3", "3"], Location::Outside, "2"),
        (["6", "3", "3"], Location::Boundary, "0"),
        (["4", "4", "4"], Location::Boundary, "0"),
        // A ray from here runs through edges and vertices of both shells.
        (["-1", "2", "2"], Location::Outside, "1"),
        (["5", "6", "3"], Location::Boundary, "0"),
    ];
    for ([x, y, z], location, distance) in cases {
        let query = at(plan, x, y, z);
        assert_eq!(query.location, location, "({x}, {y}, {z})");
        let Value::TrueDistance(value) = query.value else {
            panic!("a mesh has true distance");
        };
        let exact = q(distance).to_f64_nearest();
        assert!(
            value.lo() <= exact && exact <= value.hi(),
            "({x}, {y}, {z}): {value:?}"
        );
        assert!(value.width() < 1e-9);
    }
    // No threshold: one part in 10^30 decides.
    assert_eq!(
        at(plan, "6.000000000000000000000000000001", "3", "3").location,
        Location::Outside
    );
    assert_eq!(
        at(plan, "5.999999999999999999999999999999", "3", "3").location,
        Location::Inside
    );
    let fast = plan
        .fast([1.0, 3.0, 3.0], &mut Workspace::default())
        .unwrap();
    assert!((fast.value + 1.0).abs() < 1e-12, "{fast:?}");
    assert!(fast.true_distance_formula);
}

#[test]
fn both_kernels_give_the_exact_polyhedron_volume() {
    let built = both(
        &cube_with_cavity(),
        "1",
        ["1", "0", "0", "0"],
        ["0", "0", "0"],
    );
    let properties = field_properties(&built.plan);
    assert_eq!(properties.volume.lo, q("208"));
    assert_eq!(properties.volume.hi, q("208"));
    let centroid = properties.centroid.unwrap();
    for c in &centroid {
        assert_eq!((c.lo.clone(), c.hi.clone()), (q("3"), q("3")));
    }
    // Centroidal inertia of a 6-cube less a 2-cube: (6^5 - 2^5) / 6 on the diagonal.
    let inertia = properties.inertia.unwrap();
    assert_eq!(inertia[0][0].lo, q("7744") / q("6"));
    assert_eq!(inertia[0][1].lo, q("0"));
    assert_eq!(brep_volume(&built.brep), q("208"));
}

#[test]
fn scaled_rotated_meshes_agree_exactly_between_kernels() {
    // An octahedron with decimal vertices: volume (4/3) a^3 for half-diagonal a.
    let a = "1.25";
    let points: Vec<[String; 3]> = [
        [a, "0", "0"],
        ["-1.25", "0", "0"],
        ["0", a, "0"],
        ["0", "-1.25", "0"],
        ["0", "0", a],
        ["0", "0", "-1.25"],
    ]
    .iter()
    .map(|p| p.map(String::from))
    .collect();
    let mut triangles = Vec::new();
    for (x, y) in [(0, 2), (2, 1), (1, 3), (3, 0)] {
        triangles.push([x, y, 4]);
        triangles.push([y, x, 5]);
    }
    let built = both(
        &ascii(&points, &triangles),
        "2.5",
        ["3", "1", "-2", "5"],
        ["7", "-11", "0.125"],
    );
    let side = q(a) * q("2.5");
    let expected = q("4") / q("3") * side.clone() * side.clone() * side;
    assert_eq!(exact_field_volume(&built.plan), expected);
    assert_eq!(brep_volume(&built.brep), expected);
    let centroid = field_properties(&built.plan).centroid.unwrap();
    assert_eq!(centroid[0].lo, q("7"));
    assert_eq!(centroid[2].hi, q("0.125"));
    assert_eq!(
        at(&built.plan, "7", "-11", "0.125").location,
        Location::Inside
    );
}

#[test]
fn a_curved_binary_mesh_has_one_exact_volume_in_both_kernels() {
    let (points, triangles) = sphere(50.0, 24, 48);
    let built = both(
        &binary(&points, &triangles),
        "1",
        ["1", "0", "0", "0"],
        ["0", "0", "0"],
    );
    let field = exact_field_volume(&built.plan);
    assert_eq!(brep_volume(&built.brep), field);
    let approximate = field.to_f64_nearest();
    let ball = 4.0 / 3.0 * std::f64::consts::PI * 50.0f64.powi(3);
    assert!(
        approximate < ball && approximate > 0.98 * ball,
        "{approximate}"
    );
    assert_eq!(at(&built.plan, "0", "0", "0").location, Location::Inside);
    assert_eq!(at(&built.plan, "0", "0", "51").location, Location::Outside);
}

#[test]
fn the_field_refuses_open_crossing_and_missing_meshes() {
    let (points, mut triangles) = cuboid([0, 0, 0], [2, 2, 2], 0, false);
    let closed = triangles.clone();
    triangles.pop();
    assert_eq!(
        field_refusal(&ascii(&points, &triangles)),
        "import.mesh-open"
    );
    let (mut points, mut triangles) = (points, closed.clone());
    let (other, other_triangles) = cuboid([1, 1, 1], [3, 3, 3], 8, false);
    points.extend(other);
    triangles.extend(other_triangles);
    assert_eq!(
        field_refusal(&ascii(&points, &triangles)),
        "import.mesh-crossing"
    );
    let (points, _) = cuboid([0, 0, 0], [2, 2, 2], 0, false);
    let bytes = ascii(&points, &closed);
    let saved = recipe(
        &sha256_hex(&bytes),
        "1",
        ["1", "0", "0", "0"],
        ["0", "0", "0"],
    );
    let missing = Plan::compile_checked(&saved, &mut Decider::default()).unwrap_err();
    assert_eq!(missing.code, "field.recipe");
    assert!(
        missing.detail.contains("missing mesh source"),
        "{}",
        missing.detail
    );
}

#[test]
fn a_mesh_above_the_stated_limit_is_refused_in_both_kernels() {
    let announced = catallax_import::mesh::MAX_TRIANGLES + 1;
    let mut bytes = vec![b' '; 80];
    bytes.extend_from_slice(&(announced as u32).to_le_bytes());
    bytes.resize(84 + 50 * announced, 0);
    assert_eq!(field_refusal(&bytes), "import.mesh-too-large");
    let hash = sha256_hex(&bytes);
    let files = BTreeMap::from([(hash.clone(), &bytes[..])]);
    let refused = catallax_features::evaluate_with_files(
        &recipe(&hash, "1", ["1", "0", "0", "0"], ["0", "0", "0"]),
        &files,
        &mut Decider::default(),
    )
    .unwrap_err();
    assert!(
        refused.to_string().contains("import.mesh-too-large"),
        "{refused}"
    );
}

#[test]
fn an_empty_mesh_is_the_empty_field() {
    let built = both(
        b"solid nothing\nendsolid nothing\n",
        "1",
        ["1", "0", "0", "0"],
        ["0", "0", "0"],
    );
    assert!(built.plan.is_proven_empty());
    assert_eq!(field_properties(&built.plan).volume.hi, q("0"));
}

#[test]
fn a_mesh_in_a_field_boolean_classifies_but_has_no_mass_coverage() {
    let (points, triangles) = cuboid([0, 0, 0], [2, 2, 2], 0, false);
    let bytes = ascii(&points, &triangles);
    let hash = sha256_hex(&bytes);
    let identity = json!({"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]});
    let saved = Recipe::parse(
        &json!({
            "format": "catallax-recipe/4",
            "name": "mesh boolean",
            "steps": [
                {"id": "scan", "op": "import_mesh", "file": {"name": "part.stl", "sha256": hash},
                 "pose": identity, "scale": "1"},
                {"id": "block", "op": "block", "size": ["1", "1", "4"], "pose": identity},
                {"id": "cut", "op": "boolean", "kind": "subtract", "target": "scan", "tool": "block"}
            ],
            "result": "cut"
        })
        .to_string(),
    )
    .unwrap();
    let files = BTreeMap::from([(hash, &bytes[..])]);
    let plan = Plan::compile_with_files(&saved, &files, &mut Decider::default()).unwrap();
    assert_eq!(at(&plan, "3/2", "3/2", "1").location, Location::Inside);
    assert_eq!(at(&plan, "1/2", "1/2", "1").location, Location::Outside);
    assert_eq!(at(&plan, "3", "1", "1").location, Location::Outside);
    let refused = mass::properties(
        &plan,
        &Accuracy::default(),
        &Options::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(refused.code, "field.print-coverage");
}

#[test]
fn field_mesh_work_is_charged_to_the_callers_allowance() {
    let bytes = cube_with_cavity();
    let hash = sha256_hex(&bytes);
    let saved = recipe(&hash, "1", ["1", "0", "0", "0"], ["0", "0", "0"]);
    let files = BTreeMap::from([(hash, &bytes[..])]);
    let (refused, work) = catallax_numbers::with_operation_budget(100, || {
        Plan::compile_with_files(&saved, &files, &mut Decider::default()).map(|_| ())
    });
    assert!(work.exhausted);
    assert_eq!(refused.unwrap_err().code, "numbers.operation-budget");
    let plan = Plan::compile_with_files(&saved, &files, &mut Decider::default()).unwrap();
    let mut late = Decider::default();
    late.set_deadline(Some(std::time::Instant::now()));
    let point = RPoint3::new(q("1"), q("3"), q("3"));
    assert!(
        plan.proven_checked(&point, &mut Workspace::default(), &mut late)
            .is_err()
    );
    let (moments, work) = catallax_numbers::with_operation_budget(10, || {
        mass::properties(
            &plan,
            &Accuracy::default(),
            &Options::default(),
            &mut Decider::default(),
        )
        .map(|_| ())
    });
    assert!(work.exhausted);
    assert_eq!(moments.unwrap_err().code, "numbers.operation-budget");
}

#[test]
fn brep_mesh_work_is_charged_to_the_callers_allowance() {
    let bytes = cube_with_cavity();
    let hash = sha256_hex(&bytes);
    let saved = recipe(&hash, "1", ["1", "0", "0", "0"], ["0", "0", "0"]);
    let files = BTreeMap::from([(hash, &bytes[..])]);
    let (refused, work) = catallax_numbers::with_operation_budget(100, || {
        catallax_features::evaluate_with_files(&saved, &files, &mut Decider::default()).map(|_| ())
    });
    assert!(work.exhausted);
    assert!(
        refused
            .unwrap_err()
            .to_string()
            .contains("numbers.operation-budget")
    );
}
