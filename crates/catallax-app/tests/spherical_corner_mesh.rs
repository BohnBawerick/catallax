//! Fully filleted corners at the app's 15-degree display target and its 2-degree
//! drawing target, through the real recipe build and both display paths. Parts,
//! inventories and options come from the geometry and the app's option conversion,
//! never from a mesh.

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
use catallax_shapes::{Circle, Curve, PointDef, Pose, RPoint3, RVector3, Rotation, Surface};
use catallax_topology::{
    self as topology, Body, BodyId, BodyKind, Coedge, Edge, EdgeId, Face, FaceId, Loop, Lump,
    Model, Sense, Shell, Vertex,
};
use serde_json::json;

/// One app request: construction plus the display it asked for.
const REQUEST_WORK: u64 = 20_000_000;
/// The app's 15 degrees, from its lower bound 3.14159265 on pi.
const ANGLE: &str = "314159265/1200000000";
/// The hidden-line drawing view's 2 degrees, from the same bound on pi.
const DRAWING_ANGLE: &str = "314159265/9000000000";

struct Part {
    name: &'static str,
    radius: &'static str,
    outline: &'static [(i64, i64)],
    placed: bool,
    /// The app's Auto chord, 0.002 of the body diagonal, halved by its adapter.
    chord: &'static str,
    /// The drawing view's chord, 2e-5 of the body diagonal, halved by the adapter.
    drawing_chord: &'static str,
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
    drawing_chord: "5075431016179809/40000000000000000000",
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
    drawing_chord: "26776855677991773/200000000000000000000",
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
    drawing_chord: "28284271247461907/200000000000000000000",
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
    drawing_chord: "1658178384720883/10000000000000000000",
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

    fn drawing(&self) -> MeshOptions {
        MeshOptions {
            chord_tolerance: self.drawing_chord.parse().unwrap(),
            angular_deviation: Some(DRAWING_ANGLE.parse().unwrap()),
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
    check_counts(part.name, built, part.kinds, part.edges, part.vertices);
}

fn check_counts(name: &str, built: &Built, expected: [usize; 3], edges: usize, vertices: usize) {
    assert_eq!(built.model.check(built.result), vec![], "{name}");
    assert!(built.model.body(built.result).marks.is_empty());
    let mut kinds = [0; 3];
    for face in built.model.faces_of_body(built.result) {
        match built.model.geom().surface(built.model.face(face).surface) {
            Surface::Plane(_) => kinds[0] += 1,
            Surface::Cylinder(_) => kinds[1] += 1,
            Surface::Sphere(_) => kinds[2] += 1,
            other => panic!("{name}: unexpected surface {other:?}"),
        }
    }
    assert_eq!(kinds, expected, "{name}");
    assert_eq!(
        built.model.edges_of_body(built.result).len(),
        edges,
        "{name}"
    );
    assert_eq!(
        built.model.vertices_of_body(built.result).len(),
        vertices,
        "{name}"
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

fn check_certified(
    part: &Part,
    built: &Built,
    answer: Answered<Mesh>,
    options: &MeshOptions,
) -> Mesh {
    let mesh = check_certified_reports(part.name, built, answer, options);
    assert_eq!(mesh.unavailable_controls, vec![], "{}", part.name);
    mesh
}

/// Every certified report and indexed closure, whatever controls were unavailable.
fn check_certified_reports(
    name: &str,
    built: &Built,
    answer: Answered<Mesh>,
    options: &MeshOptions,
) -> Mesh {
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
    mesh
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

fn boundary_edges(model: &Model, face: FaceId) -> BTreeSet<EdgeId> {
    model
        .face(face)
        .loops
        .iter()
        .flat_map(|&boundary| match model.loop_(boundary) {
            Loop::Coedges(coedges) => coedges.as_slice(),
            Loop::Vertex(_) => &[],
        })
        .map(|&coedge| model.coedge(coedge).edge)
        .collect()
}

fn on_surface(surface: &Surface, p: [f64; 3], tolerance: f64) -> bool {
    let along = |origin: [f64; 3]| -> [f64; 3] { std::array::from_fn(|i| p[i] - origin[i]) };
    let norm = |v: [f64; 3]| v.iter().map(|x| x * x).sum::<f64>().sqrt();
    match surface {
        Surface::Sphere(s) => {
            let center = [&s.center.x, &s.center.y, &s.center.z].map(|x| x.to_f64());
            (norm(along(center)) - s.radius.to_f64()).abs() <= tolerance
        }
        Surface::Cylinder(s) => {
            let origin = [&s.origin.x, &s.origin.y, &s.origin.z].map(|x| x.to_f64());
            let axis = [&s.axis.x, &s.axis.y, &s.axis.z].map(|x| x.to_f64());
            let length = norm(axis);
            let axis = axis.map(|x| x / length);
            let d = along(origin);
            let height: f64 = (0..3).map(|i| d[i] * axis[i]).sum();
            let radial = std::array::from_fn(|i| d[i] - height * axis[i]);
            (norm(radial) - s.radius().to_f64()).abs() <= tolerance
        }
        other => panic!("unexpected surface {other:?}"),
    }
}

/// Reaches a corner sphere, the cylinder strip beside it and their shared edge
/// through retained handles, never through compacted mesh indices.
fn check_corner_and_strip(built: &Built, mesh: &Mesh) {
    let model = &built.model;
    let faces = model.faces_of_body(built.result);
    let surface = |face: FaceId| model.geom().surface(model.face(face).surface);
    let sphere = *faces
        .iter()
        .find(|&&face| matches!(surface(face), Surface::Sphere(_)))
        .unwrap();
    let arcs = boundary_edges(model, sphere);
    assert_eq!(arcs.len(), 3);
    let strip = *faces
        .iter()
        .find(|&&face| {
            matches!(surface(face), Surface::Cylinder(_))
                && !boundary_edges(model, face).is_disjoint(&arcs)
        })
        .unwrap();
    let mut shared = BTreeSet::new();
    for face in [sphere, strip] {
        let report = mesh.faces.iter().find(|report| report.face == face);
        let triangles = &mesh.triangles[report.unwrap().triangles.clone()];
        assert!(!triangles.is_empty());
        for triangle in triangles {
            assert_eq!(triangle.face, face);
            for i in 0..3 {
                let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                assert!(on_surface(surface(face), mesh.positions[a], 1e-3));
                shared.insert((face, a.min(b), a.max(b)));
            }
        }
    }
    let edge = boundary_edges(model, strip)
        .intersection(&arcs)
        .copied()
        .next()
        .unwrap();
    let polyline = mesh.edges.iter().find(|report| report.edge == edge);
    for pair in polyline.unwrap().vertices.windows(2) {
        let segment = (pair[0].min(pair[1]), pair[0].max(pair[1]));
        assert!(shared.contains(&(sphere, segment.0, segment.1)));
        assert!(shared.contains(&(strip, segment.0, segment.1)));
    }
}

/// One fresh construction-plus-display request per path, at the drawing target.
fn check_drawing(part: &Part) {
    let name = part.name;
    let options = part.drawing();
    let retained = format!("{:?}", part.built());
    let (built, answer, work) = part.certified(&options);
    eprintln!("{name} drawing certified request work {}", work.used);
    assert!(!work.exhausted, "{name}: {work:?}");
    check_inventory(part, &built);
    let mesh = check_certified(part, &built, answer, &options);
    check_corner_and_strip(&built, &mesh);
    assert!(format!("{built:?}") == retained, "{name}: display changed");
    let (built, mesh, work) = part.preview(&PreviewOptions::from(&options));
    eprintln!("{name} drawing preview request work {}", work.used);
    assert!(!work.exhausted, "{name}: {work:?}");
    check_preview(part, &built, mesh);
    assert!(format!("{built:?}") == retained, "{name}: preview changed");
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
fn skew_corners_display_at_the_drawing_target() {
    check_drawing(&H_N);
}

#[test]
fn small_radius_skew_corners_display_at_the_drawing_target() {
    check_drawing(&H_E);
}

#[test]
fn block_corners_display_at_the_drawing_target() {
    check_drawing(&H_B);
}

#[test]
fn placed_skew_corners_display_at_the_drawing_target() {
    check_drawing(&H_R);
}

#[test]
fn invalid_options_still_refuse_on_a_corner_part() {
    let built = H_N.built();
    for app in [H_N.options(), H_N.drawing()] {
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
}

#[test]
fn a_valid_vertex_cap_too_small_for_the_part_refuses_whole() {
    let built = H_N.built();
    for app in [H_N.options(), H_N.drawing()] {
        let options = MeshOptions {
            max_vertices: 4,
            ..app
        };
        let mut dec = Decider::default();
        let refusal = mesh_body(&built.model, built.result, &options, &mut dec).unwrap_err();
        assert_eq!(refusal.code, "mesh.resource-limit");
        let preview = PreviewOptions::from(&options);
        let refusal =
            preview_mesh_body(&built.model, built.result, &preview, &mut dec).unwrap_err();
        assert_eq!(refusal.code, "mesh.precision-or-budget");
    }
}

/// Exact squared length of a mesh segment, from its actual binary64 coordinates.
fn length_squared(mesh: &Mesh, a: usize, b: usize) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| {
        let delta = Rational::from_f64_exact(mesh.positions[a][i]).unwrap()
            - Rational::from_f64_exact(mesh.positions[b][i]).unwrap();
        &sum + &(&delta * &delta)
    })
}

fn check_segment_lengths(mesh: &Mesh, squared_limit: &Rational) {
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            assert!(length_squared(mesh, a, b) <= *squared_limit);
        }
    }
    for edge in &mesh.edges {
        for pair in edge.vertices.windows(2) {
            assert!(length_squared(mesh, pair[0], pair[1]) <= *squared_limit);
        }
    }
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
    check_segment_lengths(&mesh, &limit);
    let angle = options.angular_deviation.as_ref().unwrap();
    for edge in &mesh.edges {
        assert!(edge.max_chord_error <= options.chord_tolerance);
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

/// The positive octant of a ball is one admitted corner with no cylinder strip, so
/// these controls reach the drawing selector. A long strip under an edge limit at
/// this target is a separate cylinder refinement fault, not part of this family.
#[test]
fn other_requested_controls_still_apply_at_the_drawing_target() {
    let identity = json!({"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]});
    let recipe = three_arc_recipe(
        "drawing octant optional controls",
        json!([
            {"id": "ball", "op": "sphere", "radius": "3", "pose": identity},
            {"id": "window", "op": "block", "size": ["9", "9", "9"], "pose": identity},
            {"id": "result", "op": "boolean", "kind": "intersect",
             "target": "ball", "tool": "window"}
        ]),
    );
    let limit: Rational = "1/8".parse().unwrap();
    let options = MeshOptions {
        max_edge_length: Some(limit.clone()),
        max_aspect_ratio: Some(Rational::one()),
        max_vertices: 1_000_000,
        ..H_B.drawing()
    };
    let ((built, answer), work) = with_operation_budget(REQUEST_WORK, || {
        let mut dec = Decider::default();
        let built = evaluate(&recipe, &mut dec).unwrap();
        let answer = mesh_body(&built.model, built.result, &options, &mut dec);
        (built, answer)
    });
    eprintln!("octant drawing certified request work {}", work.used);
    assert!(!work.exhausted, "{work:?}");
    check_counts("octant", &built, [3, 0, 1], 6, 4);
    let body = built.model.body(built.result);
    assert_eq!(body.lumps.len(), 1);
    assert_eq!(built.model.lump(body.lumps[0]).shells.len(), 1);
    let corners: BTreeSet<_> = built
        .model
        .vertices_of_body(built.result)
        .into_iter()
        .map(|vertex| {
            let point = built.model.vertex(vertex).point;
            let p = built.model.geom().position(point);
            [&p.x, &p.y, &p.z].map(|x| x.rational_ref().unwrap().to_string())
        })
        .collect();
    let expected = [[0, 0, 0], [3, 0, 0], [0, 3, 0], [0, 0, 3]];
    assert_eq!(corners, expected.map(|p| p.map(|x| x.to_string())).into());
    let mesh = check_certified_reports("octant", &built, answer, &options);
    assert_eq!(
        mesh.faces.iter().map(|f| f.triangles.len()).sum::<usize>(),
        mesh.triangles.len()
    );
    check_segment_lengths(&mesh, &(&limit * &limit));
    // A nondegenerate triangle's aspect is at least 2/sqrt(3), so every face
    // reports this limit as unavailable and no other control joins it.
    let faces: BTreeSet<_> = mesh.faces.iter().map(|f| f.face).collect();
    assert_eq!(mesh.unavailable_controls.len(), 4);
    assert_eq!(
        mesh.unavailable_controls
            .iter()
            .map(|unavailable| unavailable.face)
            .collect::<BTreeSet<_>>(),
        faces
    );
    assert!(
        mesh.unavailable_controls
            .iter()
            .all(|unavailable| unavailable.control == MeshControl::AspectRatio)
    );
}

#[test]
fn the_caller_owns_the_display_work() {
    let built = H_N.built();
    for options in [H_N.options(), H_N.drawing()] {
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
}

#[test]
fn both_display_paths_repeat_exactly() {
    check_repeats(H_N.options());
}

#[test]
fn both_display_paths_repeat_exactly_at_the_drawing_target() {
    check_repeats(H_N.drawing());
}

fn check_repeats(options: MeshOptions) {
    let built = H_N.built();
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

/// The spherical outside of a removed octant: three arcs, but the long way round
/// the corner's cone, so the face is the complement of the convex sector.
fn sphere_less_an_octant() -> Built {
    let identity = json!({"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]});
    let recipe = three_arc_recipe(
        "sphere less an octant",
        json!([
            {"id": "ball", "op": "sphere", "radius": "3", "pose": identity},
            {"id": "window", "op": "block", "size": ["9", "9", "9"], "pose": identity},
            {"id": "result", "op": "boolean", "kind": "subtract",
             "target": "ball", "tool": "window"}
        ]),
    );
    evaluate(&recipe, &mut Decider::default()).unwrap()
}

/// A ball removed around a block vertex: a convex sector in the rays, but its
/// face normal points into the ball, so this corner is concave in the material.
fn block_less_a_ball_at_its_vertex() -> Built {
    let identity = json!({"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]});
    let recipe = three_arc_recipe(
        "block less a ball",
        json!([
            {"id": "stock", "op": "block", "size": ["9", "9", "9"], "pose": identity},
            {"id": "ball", "op": "sphere", "radius": "3", "pose": identity},
            {"id": "result", "op": "boolean", "kind": "subtract",
             "target": "stock", "tool": "ball"}
        ]),
    );
    evaluate(&recipe, &mut Decider::default()).unwrap()
}

fn three_arc_recipe(name: &str, steps: serde_json::Value) -> Recipe {
    Recipe::parse(
        &json!({"format": "catallax-recipe/4", "name": name, "steps": steps, "result": "result"})
            .to_string(),
    )
    .unwrap()
}

// Two faces of one sphere split along three great-circle arcs, the first of them
// taking the long way from x to y. Neither side is the cone of the three rays.
fn sphere_split_along_a_long_arc(model: &mut Model) -> BodyId {
    let original = topology::sphere(model, 2.into(), &Pose::identity()).unwrap();
    let surface = model.face(model.faces_of_body(original)[0]).surface;
    let axes: [[i64; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let [x, y, z] = axes.map(|axis| {
        let [a, b, c] = axis.map(|v| Rational::from_i64(2 * v));
        let point = model
            .geom_mut()
            .add_point(PointDef::Explicit(RPoint3::new(a, b, c).into()));
        model.add_vertex(Vertex { point })
    });
    let mut dec = Decider::default();
    let mut arc = |model: &mut Model, normal: [i64; 3], seam: [i64; 3], ends| {
        let [a, b, c] = normal.map(Rational::from_i64);
        let [d, e, f] = seam.map(Rational::from_i64);
        let circle = Circle::new(
            RPoint3::origin().into(),
            RVector3::new(a, b, c).into(),
            4.into(),
            RVector3::new(d, e, f).into(),
            &mut dec,
        )
        .unwrap();
        let curve = model.geom_mut().add_curve(Curve::Circle(circle));
        let (start, end) = ends;
        model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(start),
            end: Some(end),
        })
    };
    let edges = [
        arc(model, [0, 0, -1], [1, 0, 0], (x, y)),
        arc(model, [1, 0, 0], [0, 1, 0], (y, z)),
        arc(model, [0, 1, 0], [0, 0, 1], (z, x)),
    ];
    let mut faces = Vec::new();
    for sense in [Sense::Same, Sense::Reversed] {
        let mut order = edges;
        if sense == Sense::Reversed {
            order.reverse();
        }
        let coedges = order
            .map(|edge| model.add_coedge(Coedge { edge, sense }))
            .to_vec();
        let boundary = model.add_loop(Loop::Coedges(coedges));
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

/// Captured from main 96ef134 before the convex corner admission.
const NEIGHBOUR_BASELINE: &[(&str, &str)] = &[
    (
        "complement certified",
        "mesh a6eaac5422835e8c obj 4c0192c904056c49",
    ),
    ("complement preview", "preview 93f35bbe54053944"),
    (
        "concave certified",
        "mesh 23363c3ee6b53ccc obj 75b53058bc18948a",
    ),
    ("concave preview", "preview 1a4343cd6b137998"),
    (
        "long arc certified",
        "mesh 6a4f2ca239db69d4 obj a44477bcfff123e9",
    ),
    ("long arc preview", "preview 3a901db32a00ba03"),
];

/// Three finite great-circle arcs alone do not make a convex corner. A trim with
/// the same arcs and any other region keeps its earlier route and exact output.
#[test]
fn other_three_arc_trims_keep_their_exact_output() {
    let options = MeshOptions {
        chord_tolerance: Rational::parse("0.15").unwrap(),
        angular_deviation: Some(ANGLE.parse().unwrap()),
        ..MeshOptions::default()
    };
    let complement = sphere_less_an_octant();
    let concave = block_less_a_ball_at_its_vertex();
    for built in [&complement, &concave] {
        let faces = built.model.faces_of_body(built.result);
        assert!(faces.iter().any(|&face| {
            let surface = built.model.geom().surface(built.model.face(face).surface);
            matches!(surface, Surface::Sphere(_)) && boundary_edges(&built.model, face).len() == 3
        }));
    }
    let mut split = Model::new();
    let long_arc = sphere_split_along_a_long_arc(&mut split);
    let mut actual = Vec::new();
    for (name, model, body) in [
        ("complement", &complement.model, complement.result),
        ("concave", &concave.model, concave.result),
        ("long arc", &split, long_arc),
    ] {
        actual.push((
            format!("{name} certified"),
            certified_line(model, body, &options),
        ));
        actual.push((
            format!("{name} preview"),
            preview_line(model, body, &PreviewOptions::from(&options)),
        ));
    }
    compare(NEIGHBOUR_BASELINE, actual);
}
