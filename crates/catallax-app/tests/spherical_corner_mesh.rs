//! Fully filleted corners at the app's 15-degree display target, through the
//! real recipe build and both display paths. Parts, inventories and options
//! come from the geometry and the app's option conversion, never from a mesh.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Display;

use catallax_features::{Built, evaluate};
use catallax_mesh::{
    Mesh, MeshControl, MeshOptions, PreviewMesh, PreviewOptions, mesh_body, preview_mesh_body,
};
use catallax_numbers::{
    Answered, Decider, OperationWork, Rational, Refusal, with_operation_budget,
};
use catallax_recipe::Recipe;
use catallax_shapes::{Circle, Curve, Pose, RPoint3, RVector3, Rotation, Surface};
use catallax_topology::{
    self as topology, Body, BodyId, BodyKind, Coedge, Edge, Face, Loop, Lump, Model, Sense, Shell,
};
use serde_json::json;

/// One app request: construction plus the display it asked for.
const REQUEST_WORK: u64 = 20_000_000;
/// The app's 15 degrees, from its lower bound 3.14159265 on pi.
const ANGLE: &str = "314159265/1200000000";

struct Part {
    name: &'static str,
    radius: &'static str,
    outline: &'static [(i64, i64)],
    placed: bool,
    /// The app's Auto chord, 0.002 of the body diagonal, halved by its adapter.
    chord: &'static str,
    /// Planes, cylinders and spheres of the offset prism's parallel body.
    kinds: [usize; 3],
    edges: usize,
    vertices: usize,
}

const TRIANGLE: &[(i64, i64)] = &[(8, 0), (0, 6), (0, 0)];
const RECTANGLE: &[(i64, i64)] = &[(8, 0), (8, 6), (0, 6), (0, 0)];
const H_N: Part = Part {
    name: "H-N",
    radius: "1",
    outline: TRIANGLE,
    placed: false,
    chord: "158607219255619/12500000000000000",
    kinds: [5, 9, 6],
    edges: 36,
    vertices: 18,
};
const H_E: Part = Part {
    name: "H-E",
    radius: "1/2",
    outline: TRIANGLE,
    placed: false,
    chord: "2677685567799177/200000000000000000",
    kinds: [5, 9, 6],
    edges: 36,
    vertices: 18,
};
const H_B: Part = Part {
    name: "H-B",
    radius: "1",
    outline: RECTANGLE,
    placed: false,
    chord: "282842712474619/20000000000000000",
    kinds: [6, 12, 8],
    edges: 48,
    vertices: 24,
};
const H_R: Part = Part {
    name: "H-R",
    radius: "1",
    outline: TRIANGLE,
    placed: true,
    chord: "663271353888353/40000000000000000",
    kinds: [5, 9, 6],
    edges: 36,
    vertices: 18,
};

impl Part {
    fn recipe(&self) -> Recipe {
        let (rotation, translation) = if self.placed {
            (["1", "2", "3", "4"], ["7", "11", "-13"])
        } else {
            (["1", "0", "0", "0"], ["0", "0", "0"])
        };
        let segments: Vec<_> = self
            .outline
            .iter()
            .enumerate()
            .map(|(i, (x, y))| {
                json!({
                    "kind": "line", "id": i.to_string(),
                    "to": [x.to_string(), y.to_string()]
                })
            })
            .collect();
        Recipe::parse(
            &json!({
                "format": "catallax-recipe/4", "name": self.name,
                "steps": [
                    {"id": "section", "op": "sketch",
                     "plane": {"rotation": rotation, "translation": translation},
                     "loops": [segments]},
                    {"id": "base", "op": "extrude", "profile": "section", "distance": "10"},
                    {"id": "rounded", "op": "fillet", "body": "base",
                     "radius": self.radius, "edges": "all"}
                ],
                "result": "rounded"
            })
            .to_string(),
        )
        .unwrap()
    }

    fn options(&self) -> MeshOptions {
        MeshOptions {
            chord_tolerance: self.chord.parse().unwrap(),
            angular_deviation: Some(ANGLE.parse().unwrap()),
            ..MeshOptions::default()
        }
    }

    fn built(&self) -> Built {
        let (built, work) = with_operation_budget(REQUEST_WORK, || {
            evaluate(&self.recipe(), &mut Decider::default())
        });
        assert!(!work.exhausted, "{}: {work:?}", self.name);
        built.unwrap()
    }

    /// Construction and the certified display in one request scope.
    fn certified(&self, options: &MeshOptions) -> (Built, Answered<Mesh>, OperationWork) {
        let ((built, answer), work) = with_operation_budget(REQUEST_WORK, || {
            let mut dec = Decider::default();
            let built = evaluate(&self.recipe(), &mut dec).unwrap();
            let answer = mesh_body(&built.model, built.result, options, &mut dec);
            (built, answer)
        });
        (built, answer, work)
    }

    fn preview(
        &self,
        options: &PreviewOptions,
    ) -> (Built, Result<PreviewMesh, Refusal>, OperationWork) {
        let ((built, mesh), work) = with_operation_budget(REQUEST_WORK, || {
            let mut dec = Decider::default();
            let built = evaluate(&self.recipe(), &mut dec).unwrap();
            let mesh = preview_mesh_body(&built.model, built.result, options, &mut dec);
            (built, mesh)
        });
        (built, mesh, work)
    }
}

fn check_inventory(part: &Part, built: &Built) {
    assert_eq!(built.model.check(built.result), vec![], "{}", part.name);
    assert!(built.model.body(built.result).marks.is_empty());
    let mut kinds = [0; 3];
    for face in built.model.faces_of_body(built.result) {
        match built.model.geom().surface(built.model.face(face).surface) {
            Surface::Plane(_) => kinds[0] += 1,
            Surface::Cylinder(_) => kinds[1] += 1,
            Surface::Sphere(_) => kinds[2] += 1,
            other => panic!("{}: unexpected surface {other:?}", part.name),
        }
    }
    assert_eq!(kinds, part.kinds, "{}", part.name);
    assert_eq!(
        built.model.edges_of_body(built.result).len(),
        part.edges,
        "{}",
        part.name
    );
    assert_eq!(
        built.model.vertices_of_body(built.result).len(),
        part.vertices,
        "{}",
        part.name
    );
}

/// Requires a positive facet area and returns six times its signed volume.
fn check_positive_area(positions: &[[f64; 3]], vertices: [usize; 3]) -> f64 {
    let [a, b, c] = vertices.map(|v| positions[v]);
    let ab: [f64; 3] = std::array::from_fn(|i| b[i] - a[i]);
    let ac: [f64; 3] = std::array::from_fn(|i| c[i] - a[i]);
    let cross: [f64; 3] = std::array::from_fn(|i| {
        ab[(i + 1) % 3] * ac[(i + 2) % 3] - ab[(i + 2) % 3] * ac[(i + 1) % 3]
    });
    assert!(
        cross.iter().map(|x| x * x).sum::<f64>() > 0.0,
        "zero-area facet at {a:?} {b:?} {c:?}"
    );
    (0..3).map(|i| a[i] * cross[i]).sum()
}

fn check_certified(part: &Part, built: &Built, answer: Answered<Mesh>, options: &MeshOptions) {
    let name = part.name;
    let answer = answer.unwrap_or_else(|refusal| panic!("{name}: {refusal}"));
    assert!(answer.marks.is_empty(), "{name}: {:?}", answer.marks);
    let mesh = answer.value;
    let faces: BTreeSet<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .collect();
    let edges: BTreeSet<_> = built
        .model
        .edges_of_body(built.result)
        .into_iter()
        .collect();
    assert_eq!(
        mesh.faces.iter().map(|f| f.face).collect::<BTreeSet<_>>(),
        faces
    );
    assert_eq!(mesh.faces.len(), faces.len());
    assert_eq!(
        mesh.edges.iter().map(|e| e.edge).collect::<BTreeSet<_>>(),
        edges
    );
    assert_eq!(mesh.edges.len(), edges.len());
    assert!(mesh.positions.iter().flatten().all(|x| x.is_finite()));
    assert_eq!(mesh.unavailable_controls, vec![], "{name}");
    let angle = options.angular_deviation.as_ref().unwrap();
    for report in &mesh.faces {
        assert!(!report.triangles.is_empty(), "{name}: {:?}", report.face);
        assert!(
            mesh.triangles[report.triangles.clone()]
                .iter()
                .all(|triangle| triangle.face == report.face)
        );
        assert!(report.max_chord_error <= options.chord_tolerance);
        let bound = report.max_angular_deviation.as_ref().unwrap();
        assert!(bound <= angle, "{name}: {bound} above {angle}");
        let aspect = report.max_aspect_ratio.as_ref().unwrap();
        assert!(aspect.to_f64_nearest().is_finite());
    }
    let mut volume = 0.0;
    let mut incidence = BTreeMap::<_, (usize, i32, BTreeSet<_>)>::new();
    for triangle in &mesh.triangles {
        volume += check_positive_area(&mesh.positions, triangle.vertices);
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let uses = incidence.entry((a.min(b), a.max(b))).or_default();
            uses.0 += 1;
            uses.1 += if a < b { 1 } else { -1 };
            uses.2.insert(triangle.face);
        }
    }
    assert!(
        incidence.values().all(|uses| uses.0 == 2 && uses.1 == 0),
        "{name}"
    );
    assert!(volume > 0.0, "{name}: {volume}");
    for edge in &mesh.edges {
        assert!(edge.max_chord_error <= options.chord_tolerance);
        for pair in edge.vertices.windows(2) {
            // Both incident faces reference this segment by the same indices.
            let uses = &incidence[&(pair[0].min(pair[1]), pair[0].max(pair[1]))];
            assert_eq!(uses.2.len(), 2, "{name}: {:?}", edge.edge);
        }
    }
}

fn check_preview(part: &Part, built: &Built, mesh: Result<PreviewMesh, Refusal>) {
    let name = part.name;
    let mesh = mesh.unwrap_or_else(|refusal| panic!("{name}: {refusal}"));
    let faces: BTreeSet<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .collect();
    let edges: BTreeSet<_> = built
        .model
        .edges_of_body(built.result)
        .into_iter()
        .collect();
    assert_eq!(
        mesh.faces.iter().map(|f| f.face).collect::<BTreeSet<_>>(),
        faces
    );
    assert_eq!(mesh.faces.len(), faces.len());
    assert_eq!(
        mesh.edges.iter().map(|e| e.edge).collect::<BTreeSet<_>>(),
        edges
    );
    assert_eq!(mesh.edges.len(), edges.len());
    assert!(mesh.positions.iter().flatten().all(|x| x.is_finite()));
    for report in &mesh.faces {
        assert!(!report.triangles.is_empty(), "{name}: {:?}", report.face);
        for triangle in &mesh.triangles[report.triangles.clone()] {
            assert_eq!(triangle.face, report.face);
            check_positive_area(&mesh.positions, triangle.vertices);
        }
    }
    assert_eq!(
        mesh.faces.iter().map(|f| f.triangles.len()).sum::<usize>(),
        mesh.triangles.len()
    );
}

fn check_part(part: &Part) {
    let options = part.options();
    let (built, answer, work) = part.certified(&options);
    eprintln!("{} certified request work {}", part.name, work.used);
    assert!(!work.exhausted, "{}: {work:?}", part.name);
    check_inventory(part, &built);
    check_certified(part, &built, answer, &options);
    let (built, mesh, work) = part.preview(&PreviewOptions::from(&options));
    eprintln!("{} preview request work {}", part.name, work.used);
    assert!(!work.exhausted, "{}: {work:?}", part.name);
    check_preview(part, &built, mesh);
}

#[test]
fn skew_corners_display_at_the_app_target() {
    check_part(&H_N);
}

#[test]
fn small_radius_skew_corners_display_at_the_app_target() {
    check_part(&H_E);
}

#[test]
fn block_corners_display_at_the_app_target() {
    check_part(&H_B);
}

#[test]
fn placed_skew_corners_display_at_the_app_target() {
    check_part(&H_R);
}

#[test]
fn invalid_options_still_refuse_on_a_corner_part() {
    let built = H_N.built();
    let app = H_N.options();
    for options in [
        MeshOptions {
            chord_tolerance: Rational::zero(),
            ..app.clone()
        },
        MeshOptions {
            angular_deviation: Some(Rational::zero()),
            ..app.clone()
        },
        MeshOptions {
            max_vertices: 3,
            ..app.clone()
        },
    ] {
        let mut dec = Decider::default();
        let refusal = mesh_body(&built.model, built.result, &options, &mut dec).unwrap_err();
        assert_eq!(refusal.code, "mesh.invalid-options");
        let preview = PreviewOptions::from(&options);
        let refusal =
            preview_mesh_body(&built.model, built.result, &preview, &mut dec).unwrap_err();
        assert_eq!(refusal.code, "mesh.preview-invalid-options");
    }
}

#[test]
fn a_valid_vertex_cap_too_small_for_the_part_refuses_whole() {
    let built = H_N.built();
    let options = MeshOptions {
        max_vertices: 4,
        ..H_N.options()
    };
    let mut dec = Decider::default();
    let refusal = mesh_body(&built.model, built.result, &options, &mut dec).unwrap_err();
    assert_eq!(refusal.code, "mesh.resource-limit");
    let preview = PreviewOptions::from(&options);
    let refusal = preview_mesh_body(&built.model, built.result, &preview, &mut dec).unwrap_err();
    assert_eq!(refusal.code, "mesh.precision-or-budget");
}

#[test]
fn other_requested_controls_still_apply_on_corner_sectors() {
    let limit = Rational::one();
    let options = MeshOptions {
        max_edge_length: Some(limit.clone()),
        max_aspect_ratio: Some(limit.clone()),
        ..H_B.options()
    };
    let (built, answer, work) = H_B.certified(&options);
    assert!(!work.exhausted, "{work:?}");
    let mesh = answer.unwrap().value;
    let length_squared = |a: usize, b: usize| {
        (0..3).fold(Rational::zero(), |sum, i| {
            let delta = Rational::from_f64_exact(mesh.positions[a][i]).unwrap()
                - Rational::from_f64_exact(mesh.positions[b][i]).unwrap();
            &sum + &(&delta * &delta)
        })
    };
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            assert!(length_squared(a, b) <= limit);
        }
    }
    let angle = options.angular_deviation.as_ref().unwrap();
    for edge in &mesh.edges {
        assert!(edge.max_chord_error <= options.chord_tolerance);
        for pair in edge.vertices.windows(2) {
            assert!(length_squared(pair[0], pair[1]) <= limit);
        }
    }
    for face in &mesh.faces {
        assert!(face.max_chord_error <= options.chord_tolerance);
        assert!(face.max_angular_deviation.as_ref().unwrap() <= angle);
        assert!(face.max_aspect_ratio.is_some());
    }
    // A nondegenerate triangle's longest edge over shortest altitude is at
    // least 2/sqrt(3), so aspect 1 cannot be met and edge length 1 can.
    assert!(!mesh.unavailable_controls.is_empty());
    assert!(
        mesh.unavailable_controls
            .iter()
            .all(|unavailable| unavailable.control == MeshControl::AspectRatio)
    );
    assert_eq!(
        mesh.faces.len(),
        built.model.faces_of_body(built.result).len()
    );
}

#[test]
fn the_caller_owns_the_display_work() {
    let built = H_N.built();
    let options = H_N.options();
    let (answer, work) = with_operation_budget(1, || {
        mesh_body(
            &built.model,
            built.result,
            &options,
            &mut Decider::default(),
        )
    });
    assert!(answer.is_err() && work.exhausted, "{work:?}");
    let preview = PreviewOptions::from(&options);
    let (mesh, work) = with_operation_budget(1, || {
        preview_mesh_body(
            &built.model,
            built.result,
            &preview,
            &mut Decider::default(),
        )
    });
    assert!(mesh.is_err() && work.exhausted, "{work:?}");
}

#[test]
fn both_display_paths_repeat_exactly() {
    let built = H_N.built();
    let options = H_N.options();
    let certified = || {
        with_operation_budget(REQUEST_WORK, || {
            mesh_body(
                &built.model,
                built.result,
                &options,
                &mut Decider::default(),
            )
        })
    };
    let (first, first_work) = certified();
    let (second, second_work) = certified();
    let (first, second) = (first.unwrap(), second.unwrap());
    assert_eq!(first.value, second.value);
    assert_eq!(first.marks, second.marks);
    assert_eq!(first_work.used, second_work.used);
    assert_eq!(mesh_digest(&first.value), mesh_digest(&second.value));
    let visual = PreviewOptions::from(&options);
    let preview = || {
        with_operation_budget(REQUEST_WORK, || {
            preview_mesh_body(&built.model, built.result, &visual, &mut Decider::default())
        })
    };
    let (first, first_work) = preview();
    let (second, second_work) = preview();
    let (first, second) = (first.unwrap(), second.unwrap());
    assert_eq!(preview_digest(&first), preview_digest(&second));
    assert_eq!(first_work.used, second_work.used);
}

/// A stable FNV-1a over exact bits, so a pinned baseline survives toolchains.
struct Digest(u64);

impl Digest {
    fn new() -> Self {
        Self(0xcbf2_9ce4_8422_2325)
    }
    fn bytes(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.0 = (self.0 ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    fn number(&mut self, value: usize) {
        self.bytes(&(value as u64).to_le_bytes());
    }
    fn text(&mut self, value: impl Display) {
        let text = value.to_string();
        self.number(text.len());
        self.bytes(text.as_bytes());
    }
    fn positions(&mut self, positions: &[[f64; 3]]) {
        self.number(positions.len());
        for value in positions.iter().flatten() {
            self.bytes(&value.to_bits().to_le_bytes());
        }
    }
    fn indices(&mut self, indices: &[usize]) {
        self.number(indices.len());
        for &index in indices {
            self.number(index);
        }
    }
}

fn obj_digest(mesh: &Mesh) -> u64 {
    let mut bytes = Vec::new();
    mesh.write_obj(&mut bytes).unwrap();
    let mut digest = Digest::new();
    digest.bytes(&bytes);
    digest.0
}

/// Binary64 positions, indices, per-face and per-edge reports and control order.
fn mesh_digest(mesh: &Mesh) -> u64 {
    let mut digest = Digest::new();
    digest.positions(&mesh.positions);
    digest.number(mesh.triangles.len());
    for triangle in &mesh.triangles {
        digest.indices(&triangle.vertices);
        digest.number(triangle.face.index() as usize);
    }
    digest.number(mesh.faces.len());
    for face in &mesh.faces {
        digest.number(face.face.index() as usize);
        digest.number(face.triangles.start);
        digest.number(face.triangles.end);
        digest.text(&face.max_chord_error);
        for bound in [&face.max_aspect_ratio, &face.max_angular_deviation] {
            digest.text(bound.as_ref().map_or(String::new(), Rational::to_string));
        }
    }
    digest.number(mesh.edges.len());
    for edge in &mesh.edges {
        digest.number(edge.edge.index() as usize);
        digest.indices(&edge.vertices);
        digest.text(&edge.max_chord_error);
    }
    digest.number(mesh.unavailable_controls.len());
    for unavailable in &mesh.unavailable_controls {
        digest.number(unavailable.face.index() as usize);
        digest.text(format!("{:?}", unavailable.control));
        digest.text(unavailable.reason);
    }
    digest.0
}

fn preview_digest(mesh: &PreviewMesh) -> u64 {
    let mut digest = Digest::new();
    digest.positions(&mesh.positions);
    digest.number(mesh.triangles.len());
    for triangle in &mesh.triangles {
        digest.indices(&triangle.vertices);
        digest.number(triangle.face.index() as usize);
    }
    digest.number(mesh.faces.len());
    for face in &mesh.faces {
        digest.number(face.face.index() as usize);
        digest.number(face.triangles.start);
        digest.number(face.triangles.end);
    }
    digest.number(mesh.edges.len());
    for edge in &mesh.edges {
        digest.number(edge.edge.index() as usize);
        digest.indices(&edge.vertices);
    }
    digest.0
}

/// One display call, as a comparable line: its digests or its refusal.
fn certified_line(model: &Model, body: BodyId, options: &MeshOptions) -> String {
    match mesh_body(model, body, options, &mut Decider::default()) {
        Ok(answer) => format!(
            "mesh {:016x} obj {:016x}",
            mesh_digest(&answer.value),
            obj_digest(&answer.value)
        ),
        Err(refusal) => format!("refused {refusal}"),
    }
}

fn preview_line(model: &Model, body: BodyId, options: &PreviewOptions) -> String {
    match preview_mesh_body(model, body, options, &mut Decider::default()) {
        Ok(mesh) => format!("preview {:016x}", preview_digest(&mesh)),
        Err(refusal) => format!("refused {refusal}"),
    }
}

fn compare(baseline: &[(&str, &str)], actual: Vec<(String, String)>) {
    if std::env::var_os("SPHERICAL_CORNER_CAPTURE").is_some() {
        for (name, line) in &actual {
            println!("    (\"{name}\", \"{line}\"),");
        }
        return;
    }
    assert_eq!(baseline.len(), actual.len());
    for ((name, line), expected) in actual.iter().zip(baseline) {
        assert_eq!((name.as_str(), line.as_str()), *expected);
    }
}

/// Captured from main 485f172 before the corner mesh change.
const NO_ANGLE_BASELINE: &[(&str, &str)] = &[
    (
        "H-N certified",
        "mesh e0abaec724a8a4a9 obj 77c30998202840b2",
    ),
    (
        "H-N preview",
        "refused mesh.precision-or-budget: preview sampling reached its refinement limit",
    ),
    (
        "H-E certified",
        "mesh 4984bb0a62c1b7ba obj 77dd09a4c76fdb13",
    ),
    ("H-E preview", "preview 94a275338ae952ac"),
    (
        "H-B certified",
        "mesh 5fa4aba0d7a08b1a obj 6afc4d71307b0316",
    ),
    ("H-B preview", "preview 5ea3baffce21275e"),
];

#[test]
fn removing_the_angular_target_keeps_earlier_output() {
    let mut actual = Vec::new();
    for part in [H_N, H_E, H_B] {
        let built = part.built();
        let options = MeshOptions {
            angular_deviation: None,
            ..part.options()
        };
        actual.push((
            format!("{} certified", part.name),
            certified_line(&built.model, built.result, &options),
        ));
        actual.push((
            format!("{} preview", part.name),
            preview_line(&built.model, built.result, &PreviewOptions::from(&options)),
        ));
    }
    compare(NO_ANGLE_BASELINE, actual);
}

fn primitives(model: &mut Model) -> Vec<(&'static str, BodyId)> {
    let identity = Pose::identity();
    let transformed = Pose::new(
        Rotation::from_quaternion(1.into(), 2.into(), 3.into(), 4.into()).unwrap(),
        RVector3::new(7.into(), (-2).into(), 3.into()).into(),
        &mut Decider::default(),
    )
    .unwrap();
    vec![
        (
            "block",
            topology::block(model, [2.into(), 3.into(), 5.into()], &transformed).unwrap(),
        ),
        (
            "cylinder",
            topology::cylinder(model, 2.into(), 5.into(), &transformed).unwrap(),
        ),
        (
            "sphere",
            topology::sphere(model, 2.into(), &identity).unwrap(),
        ),
        (
            "cone",
            topology::cone(model, 3.into(), 0.into(), 4.into(), &identity).unwrap(),
        ),
        (
            "frustum",
            topology::cone(model, 1.into(), 3.into(), 4.into(), &transformed).unwrap(),
        ),
        (
            "torus",
            topology::torus(model, 3.into(), 1.into(), &identity).unwrap(),
        ),
    ]
}

// Two hemispheres sharing one closed circle through the definition poles.
fn sphere_split_through_its_poles(model: &mut Model) -> BodyId {
    let original = topology::sphere(model, 2.into(), &Pose::identity()).unwrap();
    let surface = model.face(model.faces_of_body(original)[0]).surface;
    let curve = model.geom_mut().add_curve(Curve::Circle(
        Circle::new(
            RPoint3::origin().into(),
            RVector3::new(1.into(), 0.into(), 0.into()).into(),
            4.into(),
            RVector3::new(0.into(), 1.into(), 0.into()).into(),
            &mut Decider::default(),
        )
        .unwrap(),
    ));
    let edge = model.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: None,
        end: None,
    });
    let mut faces = Vec::new();
    for sense in [Sense::Same, Sense::Reversed] {
        let coedge = model.add_coedge(Coedge { edge, sense });
        let boundary = model.add_loop(Loop::Coedges(vec![coedge]));
        faces.push(model.add_face(Face {
            surface,
            sense: Sense::Same,
            loops: vec![boundary],
        }));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: catallax_numbers::MarkLog::new(),
    })
}

// A sphere inside a block corner: three circular arcs, none of them great, so
// only the great-circle test separates this face from a fillet corner sector.
fn small_circle_sector() -> Built {
    let identity = json!({"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]});
    let corner = json!({"rotation": ["1", "0", "0", "0"], "translation": ["1/2", "1/2", "1/2"]});
    let recipe = Recipe::parse(
        &json!({
            "format": "catallax-recipe/4", "name": "small circle sector",
            "steps": [
                {"id": "ball", "op": "sphere", "radius": "3", "pose": identity},
                {"id": "window", "op": "block", "size": ["9", "9", "9"], "pose": corner},
                {"id": "sector", "op": "boolean", "kind": "intersect",
                 "target": "ball", "tool": "window"}
            ],
            "result": "sector"
        })
        .to_string(),
    )
    .unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let three_arcs = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .any(|face| {
            let face = built.model.face(face);
            matches!(built.model.geom().surface(face.surface), Surface::Sphere(_))
                && matches!(face.loops.as_slice(), [boundary]
                if matches!(built.model.loop_(*boundary), Loop::Coedges(c) if c.len() == 3))
        });
    assert!(three_arcs);
    built
}

/// Captured from main 485f172 before the corner mesh change.
const UNAFFECTED_BASELINE: &[(&str, &str)] = &[
    (
        "block bare certified",
        "mesh 9a86203972b02993 obj b07d820a1de38450",
    ),
    ("block bare preview", "preview 06f7a23870c2945f"),
    (
        "block angular certified",
        "mesh b654785060e5d5fd obj b07d820a1de38450",
    ),
    ("block angular preview", "preview 06f7a23870c2945f"),
    (
        "block controlled certified",
        "mesh f921bf4624d6dbad obj 10b430b261e8cc61",
    ),
    ("block controlled preview", "preview 06f7a23870c2945f"),
    (
        "cylinder bare certified",
        "mesh a2493a55b9a8804e obj 9c7d2ecdebb36563",
    ),
    ("cylinder bare preview", "preview e5a5ab841c71b3d9"),
    (
        "cylinder angular certified",
        "mesh 3120c5757e396dd2 obj 8475465eb6a0078d",
    ),
    ("cylinder angular preview", "preview dfd37d84ce6c1c28"),
    (
        "cylinder controlled certified",
        "mesh 477c65b7b30ca553 obj 1ce4ec4b84dbde74",
    ),
    ("cylinder controlled preview", "preview e5a5ab841c71b3d9"),
    (
        "sphere bare certified",
        "mesh f6417d6e522dae1a obj bff9205fa6acc784",
    ),
    ("sphere bare preview", "preview b19df3fb3b4eee00"),
    (
        "sphere angular certified",
        "mesh 9b6ecea61842b199 obj a500dd66a406757f",
    ),
    ("sphere angular preview", "preview e8f57f0b98e53ee2"),
    (
        "sphere controlled certified",
        "mesh c23499e647fdd642 obj bff9205fa6acc784",
    ),
    ("sphere controlled preview", "preview b19df3fb3b4eee00"),
    (
        "cone bare certified",
        "mesh 0362232c099c1e5c obj 26c2807a116828bf",
    ),
    ("cone bare preview", "preview cecc1b36942e8454"),
    (
        "cone angular certified",
        "mesh 9a26082a4d1264aa obj 58430e11ff5efc98",
    ),
    ("cone angular preview", "preview d52f961b90c8665a"),
    (
        "cone controlled certified",
        "mesh f7273d0e1b4bb553 obj 517942cdeea96d63",
    ),
    ("cone controlled preview", "preview cecc1b36942e8454"),
    (
        "frustum bare certified",
        "mesh 70443497bfddf0dd obj 59731c8aa9d35e22",
    ),
    ("frustum bare preview", "preview 4b41bfb24446312a"),
    (
        "frustum angular certified",
        "mesh 9a3d41c3bd5fd492 obj 95bbab73360c55bb",
    ),
    ("frustum angular preview", "preview 4c2c33a3f7fa18b8"),
    (
        "frustum controlled certified",
        "mesh 30c641a1cd032dc5 obj aba9b2885a905fdb",
    ),
    ("frustum controlled preview", "preview b41ab6d9d7ded20b"),
    (
        "torus bare certified",
        "mesh 0ac6aebe7a4abdcb obj f6008a2fab3cd0c0",
    ),
    ("torus bare preview", "preview c852684354fc5596"),
    (
        "torus angular certified",
        "mesh b4328ee0e206ca48 obj 8f0c1f752b2c7eb6",
    ),
    ("torus angular preview", "preview 902d052c0ea8370b"),
    (
        "torus controlled certified",
        "mesh b01e83073bf41951 obj f6008a2fab3cd0c0",
    ),
    ("torus controlled preview", "preview c852684354fc5596"),
    (
        "split sphere bare certified",
        "mesh af9da8c406c718a4 obj 079013b32e981ea2",
    ),
    ("split sphere bare preview", "preview 58760f41710a411c"),
    (
        "split sphere angular certified",
        "mesh 400c20b45d3887c4 obj 3c83de46c651a873",
    ),
    ("split sphere angular preview", "preview b2cfb1c990c5f5de"),
    (
        "split sphere controlled certified",
        "mesh 249fb16e9719ec5b obj 079013b32e981ea2",
    ),
    (
        "split sphere controlled preview",
        "preview 58760f41710a411c",
    ),
    (
        "small circles bare certified",
        "mesh 58358b2523d4b71d obj a727986772c34a38",
    ),
    ("small circles bare preview", "preview 086288b4c6c2da1b"),
    (
        "small circles angular certified",
        "mesh 809e26ae886a2fa5 obj 64121fbb24fff360",
    ),
    ("small circles angular preview", "preview 3b2e1a52933987e4"),
    (
        "small circles controlled certified",
        "mesh 92d9426a86e72e93 obj c68e37a8efb3ff84",
    ),
    (
        "small circles controlled preview",
        "preview 086288b4c6c2da1b",
    ),
];

#[test]
fn faces_outside_the_family_keep_their_exact_output() {
    let chord = Rational::parse("0.15").unwrap();
    let bare = MeshOptions {
        chord_tolerance: chord,
        ..MeshOptions::default()
    };
    let angular = MeshOptions {
        angular_deviation: Some(ANGLE.parse().unwrap()),
        ..bare.clone()
    };
    let controlled = MeshOptions {
        angular_deviation: Some(Rational::parse("0.8").unwrap()),
        max_edge_length: Some(Rational::parse("1.5").unwrap()),
        max_aspect_ratio: Some(Rational::one()),
        ..bare.clone()
    };
    let mut model = Model::new();
    let mut bodies = primitives(&mut model);
    bodies.push(("split sphere", sphere_split_through_its_poles(&mut model)));
    let sector = small_circle_sector();
    let mut parts = vec![(&model, bodies)];
    parts.push((&sector.model, vec![("small circles", sector.result)]));
    let mut actual = Vec::new();
    for (model, bodies) in parts {
        for (name, body) in bodies {
            for (label, options) in [
                ("bare", &bare),
                ("angular", &angular),
                ("controlled", &controlled),
            ] {
                actual.push((
                    format!("{name} {label} certified"),
                    certified_line(model, body, options),
                ));
                actual.push((
                    format!("{name} {label} preview"),
                    preview_line(model, body, &PreviewOptions::from(options)),
                ));
            }
        }
    }
    compare(UNAFFECTED_BASELINE, actual);
}
