//! Closed mesh import in the B-rep, RECIPE.md section 4.19 and PLAN M22.

use catallax_features::{BuildError, ImportFiles, evaluate_with_files};
use catallax_import::sha256_hex;
use catallax_measure::{Accuracy, volume};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Origin, Recipe, Role, StepId};
use serde_json::json;

type Tri = [usize; 3];

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

/// An axis box's corners and twelve outward triangles; `inward` reverses them.
fn cuboid(lo: [i64; 3], hi: [i64; 3], first: usize, inward: bool) -> (Vec<[i64; 3]>, Vec<Tri>) {
    let corners = (0..8)
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

fn ascii(points: &[[i64; 3]], triangles: &[Tri]) -> Vec<u8> {
    let mut s = String::from("solid part\n");
    for t in triangles {
        s.push_str("facet normal 0 0 0\nouter loop\n");
        for &v in t {
            let [x, y, z] = points[v];
            s.push_str(&format!("vertex {x} {y} {z}\n"));
        }
        s.push_str("endloop\nendfacet\n");
    }
    s.push_str("endsolid part\n");
    s.into_bytes()
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

fn build(bytes: &[u8]) -> Result<catallax_features::Built, BuildError> {
    build_posed(bytes, "1", ["1", "0", "0", "0"], ["0", "0", "0"])
}

fn build_posed(
    bytes: &[u8],
    scale: &str,
    rotation: [&str; 4],
    translation: [&str; 3],
) -> Result<catallax_features::Built, BuildError> {
    let hash = sha256_hex(bytes);
    let files = ImportFiles::from([(hash.clone(), bytes)]);
    evaluate_with_files(
        &recipe(&hash, scale, rotation, translation),
        &files,
        &mut Decider::default(),
    )
}

fn exact_volume(built: &catallax_features::Built) -> Rational {
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

fn refusal(bytes: &[u8]) -> String {
    match build(bytes) {
        Err(BuildError::Geometry { step, refusal }) => {
            assert_eq!(step.as_str(), "scan");
            refusal.code.to_string()
        }
        other => panic!("expected a typed geometry refusal, got {other:?}"),
    }
}

/// A square frame: a 6 by 6 block with a 2 by 2 hole through it, height `h`.
/// Its top and bottom faces each have a hole loop, and its genus is one.
fn frame(h: i64) -> (Vec<[i64; 3]>, Vec<Tri>) {
    let outer = [[0, 0], [6, 0], [6, 6], [0, 6]];
    let inner = [[2, 2], [4, 2], [4, 4], [2, 4]];
    let mut points = Vec::new();
    for z in [0, h] {
        for [x, y] in outer.iter().chain(&inner) {
            points.push([*x, *y, z]);
        }
    }
    // Indices: bottom outer 0-3, bottom inner 4-7, top outer 8-11, top inner 12-15.
    let mut triangles = Vec::new();
    for i in 0..4 {
        let j = (i + 1) % 4;
        let (o, oj, n, nj) = (i, j, 4 + i, 4 + j);
        // Top annulus faces up, bottom faces down.
        triangles.push([8 + o, 8 + oj, 8 + nj]);
        triangles.push([8 + o, 8 + nj, 8 + n]);
        triangles.push([o, nj, oj]);
        triangles.push([o, n, nj]);
        // Outer wall faces away from the centre, inner wall towards it.
        triangles.push([o, oj, 8 + oj]);
        triangles.push([o, 8 + oj, 8 + o]);
        triangles.push([n, 8 + nj, nj]);
        triangles.push([n, 8 + n, 8 + nj]);
    }
    (points, triangles)
}

#[test]
fn a_closed_mesh_builds_an_exact_solid_with_merged_planar_faces() {
    let (points, triangles) = cuboid([0, 0, 0], [3, 4, 5], 0, false);
    let built = build(&ascii(&points, &triangles)).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let counts = built.model.counts(built.result);
    assert_eq!(
        (counts.lumps, counts.shells, counts.faces, counts.edges),
        (1, 1, 6, 12)
    );
    assert_eq!(exact_volume(&built), q("60"));
    let mesh = Origin {
        step: StepId::new("scan").unwrap(),
        role: Role::Mesh,
    };
    for face in built.model.faces_of_body(built.result) {
        assert_eq!(built.model.face_origins(face), std::slice::from_ref(&mesh));
    }
}

#[test]
fn a_cavity_is_an_inner_shell_and_a_hole_is_an_inner_loop() {
    let (mut points, mut triangles) = cuboid([0, 0, 0], [6, 6, 6], 0, false);
    let (inner, inner_triangles) = cuboid([2, 2, 2], [4, 4, 4], 8, true);
    points.extend(inner);
    triangles.extend(inner_triangles);
    let built = build(&ascii(&points, &triangles)).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let counts = built.model.counts(built.result);
    assert_eq!((counts.lumps, counts.shells, counts.faces), (1, 2, 12));
    assert_eq!(exact_volume(&built), q("208"));

    let (points, triangles) = frame(3);
    let built = build(&ascii(&points, &triangles)).unwrap();
    assert!(
        built.model.check(built.result).is_empty(),
        "{:?}",
        built.model.check(built.result)
    );
    let counts = built.model.counts(built.result);
    assert_eq!((counts.faces, counts.loops), (10, 12));
    assert_eq!(exact_volume(&built), q("96"));
}

#[test]
fn scale_and_pose_place_the_exact_solid() {
    let (points, triangles) = cuboid([0, 0, 0], [1, 2, 3], 0, false);
    let bytes = ascii(&points, &triangles);
    // A quarter turn about z, then a translation, with 25.4 mm per file unit.
    let built = build_posed(&bytes, "25.4", ["1", "0", "0", "1"], ["10", "0", "-1"]).unwrap();
    assert!(built.model.check(built.result).is_empty());
    assert_eq!(
        exact_volume(&built),
        q("6") * q("25.4") * q("25.4") * q("25.4")
    );
    let bounds = catallax_measure::bounding_box(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value
    .unwrap();
    // The measured box is an outward enclosure of the exact placed corners.
    let near = |enclosure: &Rational, exact: &str| {
        let gap = (enclosure - &q(exact)).abs();
        assert!(gap < q("1e-9"), "{enclosure} is not near {exact}");
    };
    near(&bounds.min[0], "-40.8");
    near(&bounds.max[0], "10");
    near(&bounds.max[1], "25.4");
    near(&bounds.min[2], "-1");
    near(&bounds.max[2], "75.2");
}

#[test]
fn open_and_crossing_meshes_are_refused_as_impossible_with_the_reason() {
    let (points, mut triangles) = cuboid([0, 0, 0], [2, 2, 2], 0, false);
    let closed = triangles.clone();
    triangles.pop();
    assert_eq!(refusal(&ascii(&points, &triangles)), "import.mesh-open");
    let (mut points, mut triangles) = (points, closed);
    let (other, other_triangles) = cuboid([1, 1, 1], [3, 3, 3], 8, false);
    points.extend(other);
    triangles.extend(other_triangles);
    assert_eq!(refusal(&ascii(&points, &triangles)), "import.mesh-crossing");
}

#[test]
fn missing_or_changed_bytes_are_recipe_errors_and_zero_scale_is_impossible() {
    let (points, triangles) = cuboid([0, 0, 0], [2, 2, 2], 0, false);
    let bytes = ascii(&points, &triangles);
    let hash = sha256_hex(&bytes);
    let saved = recipe(&hash, "1", ["1", "0", "0", "0"], ["0", "0", "0"]);
    let missing =
        evaluate_with_files(&saved, &ImportFiles::new(), &mut Decider::default()).unwrap_err();
    assert!(matches!(&missing, BuildError::Recipe(_)), "{missing}");
    assert!(
        missing.to_string().contains("missing mesh source"),
        "{missing}"
    );
    let changed = b"solid x\nendsolid x\n";
    let files = ImportFiles::from([(hash.clone(), &changed[..])]);
    let wrong = evaluate_with_files(&saved, &files, &mut Decider::default()).unwrap_err();
    assert!(wrong.to_string().contains("sha256 mismatch"), "{wrong}");
    let zero = build_posed(&bytes, "0", ["1", "0", "0", "0"], ["0", "0", "0"]).unwrap_err();
    assert!(zero.to_string().contains("import.mesh-scale"), "{zero}");
}

#[test]
fn an_empty_mesh_builds_the_empty_solid() {
    let built = build(b"solid nothing\nendsolid nothing\n").unwrap();
    assert_eq!(built.model.counts(built.result).faces, 0);
}
