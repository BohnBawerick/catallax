//! Unproved display selection, geometry and caller-owned allowance regressions.
use std::collections::BTreeSet;

use catallax_mesh::{PreviewOptions, preview_mesh_body};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::{Recipe, Role};

#[test]
fn preview_picks_retain_complete_face_origins_and_edge_names_after_edits() {
    let mut input: serde_json::Value =
        serde_json::from_str(include_str!("../benches/fine_mesh/countersink.json")).unwrap();
    let mut previous = None;
    for diameter in ["5", "4"] {
        input["steps"][1]["diameter"] = diameter.into();
        let recipe = Recipe::parse(&input.to_string()).unwrap();
        let body = catallax_app::build_preview(
            &recipe,
            &PreviewOptions::default(),
            &mut Decider::default(),
        )
        .unwrap();
        let mesh = &body.preview.mesh;
        let model = &body.built.model;
        assert_eq!(
            mesh.faces.iter().map(|f| f.face).collect::<BTreeSet<_>>(),
            model.faces_of_body(body.built.result).into_iter().collect()
        );
        assert_eq!(
            mesh.edges.iter().map(|e| e.edge).collect::<BTreeSet<_>>(),
            model.edges_of_body(body.built.result).into_iter().collect()
        );
        for face in &mesh.faces {
            assert_eq!(
                body.preview.names.faces[&face.face],
                model.face_origins(face.face)
            );
            assert!(!face.triangles.is_empty());
            for triangle in &mesh.triangles[face.triangles.clone()] {
                assert_eq!(triangle.face, face.face);
                assert!(
                    triangle
                        .vertices
                        .iter()
                        .all(|&v| mesh.positions[v].iter().all(|x| x.is_finite()))
                );
            }
        }
        for edge in &mesh.edges {
            assert_eq!(
                body.preview.names.edges[&edge.edge],
                model.edge_name(body.built.result, edge.edge)
            );
            assert!(edge.vertices.len() >= 2);
        }
        let origins: BTreeSet<_> = body.preview.names.faces.values().cloned().collect();
        assert!(
            origins
                .iter()
                .flatten()
                .any(|o| matches!(o.role, Role::Hole { .. }))
        );
        if let Some(previous) = previous {
            assert_eq!(origins, previous);
        }
        previous = Some(origins);
    }
}

#[test]
fn preview_does_not_fill_planar_hole_trims() {
    let recipe = Recipe::parse(include_str!("../benches/fine_mesh/countersink.json")).unwrap();
    let body =
        catallax_app::build_preview(&recipe, &PreviewOptions::default(), &mut Decider::default())
            .unwrap();
    for triangle in &body.preview.mesh.triangles {
        let positions = triangle.vertices.map(|v| body.preview.mesh.positions[v]);
        if positions.iter().all(|p| (p[2] - 5.0).abs() < 1e-10) {
            for center in [15.0, 45.0] {
                // A hole center must lie outside every top-face triangle.
                let turns: [f64; 3] = std::array::from_fn(|i| {
                    let a = positions[i];
                    let b = positions[(i + 1) % 3];
                    (b[0] - a[0]) * (20.0 - a[1]) - (b[1] - a[1]) * (center - a[0])
                });
                assert!(!(turns.iter().all(|x| *x >= 0.0) || turns.iter().all(|x| *x <= 0.0)));
            }
        }
    }
}

#[test]
fn preview_consumes_the_existing_allowance_and_refuses_without_partial_data() {
    let recipe = Recipe::parse(include_str!("../benches/fine_mesh/countersink.json")).unwrap();
    let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
    let options = PreviewOptions::default();
    let (result, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
        catallax_app::preview_body(&built, &options, &mut Decider::default())
    });
    let reference = result.unwrap();
    assert!(!work.exhausted);
    let (result, exact) = with_operation_budget(work.used, || {
        catallax_app::preview_body(&built, &options, &mut Decider::default())
    });
    assert_eq!(result.unwrap().mesh, reference.mesh);
    assert_eq!(exact.used, work.used);
    assert!(!exact.exhausted);
    let (result, short) = with_operation_budget(work.used - 1, || {
        catallax_app::preview_body(&built, &options, &mut Decider::default())
    });
    assert!(result.is_err());
    assert!(short.exhausted);
    let ((first, second), shared) = with_operation_budget(work.used, || {
        let first = catallax_app::preview_body(&built, &options, &mut Decider::default());
        let second = catallax_app::preview_body(&built, &options, &mut Decider::default());
        (first, second)
    });
    assert!(first.is_ok());
    assert!(second.is_err());
    assert!(shared.exhausted);
    let constrained = PreviewOptions {
        max_vertices: 4,
        ..options
    };
    assert!(
        preview_mesh_body(
            &built.model,
            built.result,
            &constrained,
            &mut Decider::default()
        )
        .is_err()
    );
}

#[test]
fn curved_fillet_end_trims_do_not_leave_large_flat_preview_facets() {
    use catallax_shapes::Surface;
    let recipe = Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    let options = PreviewOptions {
        chord_tolerance: "0.005".parse().unwrap(),
        angular_deviation: Some("0.1".parse().unwrap()),
        ..PreviewOptions::default()
    };
    let body = catallax_app::build_preview(&recipe, &options, &mut Decider::default()).unwrap();
    let dec = Decider::default();
    let scalar =
        |value: &catallax_numbers::Real| dec.enclose_at(value, 128).unwrap().mid().to_f64_nearest();
    for face in &body.preview.mesh.faces {
        let Surface::Cylinder(cylinder) = body
            .built
            .model
            .geom()
            .surface(body.built.model.face(face.face).surface)
        else {
            continue;
        };
        let origin = [&cylinder.origin.x, &cylinder.origin.y, &cylinder.origin.z].map(scalar);
        let mut axis = [&cylinder.axis.x, &cylinder.axis.y, &cylinder.axis.z].map(scalar);
        let norm = axis.iter().map(|x| x * x).sum::<f64>().sqrt();
        axis = axis.map(|x| x / norm);
        let radius = scalar(&cylinder.radius());
        for triangle in &body.preview.mesh.triangles[face.triangles.clone()] {
            let p = triangle.vertices.map(|i| body.preview.mesh.positions[i]);
            for weights in [
                [0.5, 0.5, 0.0],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [1.0 / 3.0; 3],
            ] {
                let d: [f64; 3] = std::array::from_fn(|i| {
                    (0..3).map(|j| weights[j] * p[j][i]).sum::<f64>() - origin[i]
                });
                let z = (0..3).map(|i| d[i] * axis[i]).sum::<f64>();
                let radial = (0..3)
                    .map(|i| (d[i] - z * axis[i]).powi(2))
                    .sum::<f64>()
                    .sqrt();
                assert!(
                    (radial - radius).abs() <= options.chord_tolerance.to_f64_nearest(),
                    "face {}: {p:?}",
                    face.face.index()
                );
            }
        }
    }
}

#[test]
fn retained_spline_preview_has_shared_rings_and_outward_caps() {
    use std::collections::BTreeMap;
    let recipe = Recipe::parse(&serde_json::json!({
        "format":"catallax-recipe/3","name":"preview retained spline",
        "steps":[
            {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[{"kind":"circle","id":"wall","centre":["0","0"],"radius":"2"}]]},
            {"id":"sweep","op":"extrude","profile":"profile","distance":"10","direction":["3","0","4"]},
            {"id":"window","op":"block","size":["8","8","4"],"pose":{"rotation":["3","0","1","0"],"translation":["-2","-4","4"]}},
            {"id":"source","op":"boolean","kind":"intersect","target":"sweep","tool":"window"},
            {"id":"offset","op":"offset","body":"source","distance":"-1/4"}
        ],"result":"offset"
    }).to_string()).unwrap();
    let built = catallax_features::evaluate_certified(&recipe, &mut Decider::default()).unwrap();
    let preview = catallax_app::preview_body(
        built.built(),
        &PreviewOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let mesh = preview.mesh;
    let mut incidence = BTreeMap::<_, (usize, i32)>::new();
    let mut volume6 = 0.0;
    for triangle in &mesh.triangles {
        let [a, b, c] = triangle.vertices.map(|i| mesh.positions[i]);
        volume6 += a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0]);
        for i in 0..3 {
            let a = triangle.vertices[i];
            let b = triangle.vertices[(i + 1) % 3];
            let entry = incidence.entry((a.min(b), a.max(b))).or_default();
            entry.0 += 1;
            entry.1 += if a < b { 1 } else { -1 };
        }
    }
    assert!(volume6 > 0.0);
    assert!(incidence.values().all(|&uses| uses == (2, 0)));
    assert_eq!(
        mesh.faces.iter().map(|f| f.face).collect::<BTreeSet<_>>(),
        built
            .built()
            .model
            .faces_of_body(built.built().result)
            .into_iter()
            .collect()
    );
    for edge in &mesh.edges {
        assert_eq!(edge.vertices.first(), edge.vertices.last());
        for pair in edge.vertices.windows(2) {
            assert_eq!(
                incidence[&(pair[0].min(pair[1]), pair[0].max(pair[1]))],
                (2, 0)
            );
        }
    }
}
