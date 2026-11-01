//! An edge-length limit on a cylindrical fillet strip at the drawing target, through
//! the real recipe build. The block has one rounded vertical edge, so its strip has
//! no spherical corner. The chord and 2-degree angle are the hidden-line drawing
//! row's, and aspect 1 is asked for because no triangle can meet it.

use std::collections::BTreeMap;

use catallax_features::evaluate;
use catallax_mesh::{Mesh, MeshControl, MeshOptions, mesh_body};
use catallax_numbers::{Decider, OperationWork, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{Pose, Surface};
use catallax_topology::{self as topology, BodyId, FaceId, Model};
use serde_json::{Value, json};

const REQUEST_WORK: u64 = 20_000_000;
const CHORD: &str = "28284271247461907/200000000000000000000";
const ANGLE: &str = "314159265/9000000000";
const HEIGHT: i64 = 10;

fn rounded_block(rotation: [&str; 4]) -> Recipe {
    let side = |segment: &str| json!({"step": "base", "role": {"side": {"segment": segment}}});
    let segments: Vec<Value> = [(8, 0), (8, 6), (0, 6), (0, 0)]
        .iter()
        .enumerate()
        .map(|(id, (x, y))| {
            json!({"kind": "line", "id": id.to_string(), "to": [x.to_string(), y.to_string()]})
        })
        .collect();
    Recipe::parse(
        &json!({
            "format": "catallax-recipe/4",
            "name": "one rounded edge",
            "steps": [
                {"id": "section", "op": "sketch",
                 "plane": {"rotation": rotation, "translation": ["0", "0", "0"]},
                 "loops": [segments]},
                {"id": "base", "op": "extrude", "profile": "section", "distance": HEIGHT.to_string()},
                {"id": "rounded", "op": "fillet", "body": "base", "radius": "1",
                 "edges": [{"between": [side("0"), side("1")]}]}
            ],
            "result": "rounded"
        })
        .to_string(),
    )
    .unwrap()
}

fn options(limit: Option<&str>) -> MeshOptions {
    MeshOptions {
        chord_tolerance: CHORD.parse().unwrap(),
        angular_deviation: Some(ANGLE.parse().unwrap()),
        max_edge_length: limit.map(|limit| limit.parse().unwrap()),
        max_aspect_ratio: Some(Rational::one()),
        ..MeshOptions::default()
    }
}

fn meshed(recipe: &Recipe, options: &MeshOptions) -> (Model, BodyId, Mesh, OperationWork) {
    let (answer, work) = with_operation_budget(REQUEST_WORK, || {
        let mut dec = Decider::default();
        let built = evaluate(recipe, &mut dec).unwrap();
        let mesh = mesh_body(&built.model, built.result, options, &mut dec);
        (built, mesh)
    });
    let (built, mesh) = answer;
    let limit = options.max_edge_length.as_ref().map(ToString::to_string);
    let mesh = mesh.unwrap_or_else(|refusal| panic!("limit {limit:?}: {refusal:?}"));
    (built.model, built.result, mesh.value, work)
}

fn squared_length(mesh: &Mesh, a: usize, b: usize) -> Rational {
    (0..3).fold(Rational::zero(), |sum, axis| {
        let delta = Rational::from_f64_exact(mesh.positions[a][axis]).unwrap()
            - Rational::from_f64_exact(mesh.positions[b][axis]).unwrap();
        sum + &delta * &delta
    })
}

/// Checks the returned mesh itself: closure, the two certificates, every edge against
/// the limit, and that only the unreachable aspect request is reported unavailable.
fn check(mesh: &Mesh, options: &MeshOptions) {
    let limit = options.max_edge_length.as_ref().unwrap();
    let squared_limit = limit * limit;
    let mut uses = BTreeMap::<(usize, usize), (usize, i32)>::new();
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let edge = uses.entry((a.min(b), a.max(b))).or_default();
            edge.0 += 1;
            edge.1 += if a < b { 1 } else { -1 };
            assert!(squared_length(mesh, a, b) <= squared_limit, "limit {limit}");
        }
    }
    assert!(uses.values().all(|edge| *edge == (2, 0)), "limit {limit}");
    for face in &mesh.faces {
        assert!(!face.triangles.is_empty());
        assert!(face.max_chord_error <= options.chord_tolerance);
        assert_eq!(face.max_angular_deviation, options.angular_deviation);
        assert!(face.max_aspect_ratio.is_some());
    }
    for edge in &mesh.edges {
        assert!(edge.max_chord_error <= options.chord_tolerance);
        for pair in edge.vertices.windows(2) {
            assert!(squared_length(mesh, pair[0], pair[1]) <= squared_limit);
            let key = (pair[0].min(pair[1]), pair[0].max(pair[1]));
            assert_eq!(uses[&key], (2, 0));
        }
    }
    assert!(
        mesh.unavailable_controls
            .iter()
            .all(|unavailable| unavailable.control == MeshControl::AspectRatio)
    );
}

fn cylinder_face(model: &Model, body: BodyId) -> FaceId {
    let cylinders: Vec<_> = model
        .faces_of_body(body)
        .into_iter()
        .filter(|&face| {
            matches!(
                model.geom().surface(model.face(face).surface),
                Surface::Cylinder(_)
            )
        })
        .collect();
    let [cylinder] = cylinders.as_slice() else {
        panic!("{} cylinders", cylinders.len());
    };
    *cylinder
}

/// The fault counted this face in hundreds of thousands of vertices, so the bound is
/// one lattice: a column for each sample of both rails, and the dyadic rows that hold
/// a cell at half the limit. Vertices come from the rails and the limit, not the mesh.
fn check_proportionate(model: &Model, body: BodyId, mesh: &Mesh, limit: i64) {
    let cylinder = cylinder_face(model, body);
    let face = mesh
        .faces
        .iter()
        .find(|face| face.face == cylinder)
        .unwrap();
    let mut vertices: Vec<_> = mesh.triangles[face.triangles.clone()]
        .iter()
        .flat_map(|triangle| triangle.vertices)
        .collect();
    vertices.sort_unstable();
    vertices.dedup();
    let rail = mesh
        .edges
        .iter()
        .map(|edge| {
            edge.vertices
                .iter()
                .filter(|vertex| vertices.binary_search(vertex).is_ok())
                .count()
        })
        .max()
        .unwrap();
    let rows = ((2 * HEIGHT / limit) as usize + 1).next_power_of_two();
    assert!(
        vertices.len() <= 2 * rail * (rows + 1),
        "limit {limit}: {} vertices from {rail} rail samples",
        vertices.len()
    );
}

#[test]
fn an_edge_limit_keeps_the_strip_a_lattice_at_the_drawing_target() {
    let recipe = rounded_block(["1", "0", "0", "0"]);
    // 1 is the reported row. 8 splits the ruling, and 10 keeps the ruling whole while
    // the first stitched diagonal exceeds it, so both ways off the strip are covered.
    for limit in [1, 8, HEIGHT] {
        let options = options(Some(&limit.to_string()));
        let (model, body, mesh, work) = meshed(&recipe, &options);
        assert!(!work.exhausted, "limit {limit}: {work:?}");
        check(&mesh, &options);
        check_proportionate(&model, body, &mesh, limit);
    }
}

#[test]
fn an_edge_limit_holds_on_a_skew_strip() {
    let recipe = rounded_block(["1", "2", "3", "4"]);
    let options = options(Some("1"));
    let (model, body, mesh, work) = meshed(&recipe, &options);
    assert!(!work.exhausted, "{work:?}");
    check(&mesh, &options);
    check_proportionate(&model, body, &mesh, 1);
}

/// A limit the stitched strip already meets changes nothing about its mesh.
#[test]
fn a_limit_the_strip_already_meets_keeps_its_mesh() {
    let recipe = rounded_block(["1", "0", "0", "0"]);
    let limited = options(Some("16"));
    let (_, _, within, work) = meshed(&recipe, &limited);
    assert!(!work.exhausted, "{work:?}");
    check(&within, &limited);
    assert_eq!(within.positions.len(), 990);
    let (_, _, unset, _) = meshed(&recipe, &options(None));
    assert_eq!(within, unset);
}

/// A closed cylinder reaches the same grid. The first limit refused with an open
/// triangulation before the face took rows. The second is the pinned plain cylinder.
#[test]
fn a_closed_cylinder_takes_rows_under_a_short_limit() {
    for (radius, height, chord, angle, limit) in [
        (1, 8, "1/100", "0.2617993877991494", "1/4"),
        (2, 5, "0.15", "0.8", "1.5"),
    ] {
        let mut model = Model::new();
        let body = topology::cylinder(&mut model, radius.into(), height.into(), &Pose::identity())
            .unwrap();
        let options = MeshOptions {
            chord_tolerance: chord.parse().unwrap(),
            angular_deviation: Some(angle.parse().unwrap()),
            max_edge_length: Some(limit.parse().unwrap()),
            max_aspect_ratio: Some(Rational::one()),
            ..MeshOptions::default()
        };
        let (mesh, work) = with_operation_budget(REQUEST_WORK, || {
            mesh_body(&model, body, &options, &mut Decider::default())
        });
        assert!(!work.exhausted, "{work:?}");
        check(&mesh.unwrap().value, &options);
    }
}
