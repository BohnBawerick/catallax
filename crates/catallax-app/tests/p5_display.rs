//! The full P5 housing through the certified application route: one closed
//! mesh across every boundary, per-face bounds and stable picking names.

use std::collections::{BTreeMap, BTreeSet};

use catallax_app::{ApplicationBody, build_certified};
use catallax_mesh::{MeshControl, MeshOptions, PreviewOptions, preview_mesh_body};
use catallax_numbers::{BudgetPolicy, Decider, Rational};
use catallax_recipe::{Origin, Recipe, Role};
use serde_json::{Value, json};

fn recipe(radius: &str) -> Recipe {
    let mut recipe: Value =
        serde_json::from_str(include_str!("../../../parts/P5-flanged-housing.json")).unwrap();
    recipe["steps"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|step| step["id"] == "filleted")
        .unwrap()["radius"] = json!(radius);
    Recipe::parse(&recipe.to_string()).unwrap()
}
// The definition under the exact rotation (3/5,-4/5;4/5,3/5) about z, quaternion
// (2,0,0,1), then the translation (7,11,-13). Every frame, point and axis moves
// before evaluation, so nothing rotates the finished body.
fn placed(radius: &str) -> Recipe {
    let source = include_str!("../../../parts/P5-flanged-housing.json");
    let mut recipe: Value = serde_json::from_str(source).unwrap();
    let moved = |vector: &mut Value, by: [i64; 3]| {
        let [x, y, z] = [0, 1, 2].map(|i| Rational::parse(vector[i].as_str().unwrap()).unwrap());
        let fifths = |value: i64| Rational::from_ratio(value, 5).unwrap();
        let new = [
            fifths(3) * x.clone() - fifths(4) * y.clone() + Rational::from_i64(by[0]),
            fifths(4) * x + fifths(3) * y + Rational::from_i64(by[1]),
            z + Rational::from_i64(by[2]),
        ];
        *vector = json!(new.map(|part| part.to_string()));
    };
    let compose = |rotation: &mut Value| {
        let [w, x, y, z] =
            [0, 1, 2, 3].map(|i| rotation[i].as_str().unwrap().parse::<i64>().unwrap());
        *rotation =
            json!([2 * w - z, 2 * x - y, 2 * y + x, 2 * z + w].map(|part| part.to_string()));
    };
    for step in recipe["steps"].as_array_mut().unwrap() {
        for key in ["plane", "pose"] {
            if !step[key].is_null() {
                moved(&mut step[key]["translation"], [7, 11, -13]);
                compose(&mut step[key]["rotation"]);
            }
        }
        for path in ["/axis", "/layout/circular/axis"] {
            if let Some(axis) = step.pointer_mut(path) {
                moved(&mut axis["point"], [7, 11, -13]);
                moved(&mut axis["direction"], [0, 0, 0]);
            }
        }
        if step["id"] == "filleted" {
            step["radius"] = json!(radius);
        }
    }
    Recipe::parse(&recipe.to_string()).unwrap()
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    })
}

#[test]
fn the_mixed_body_displays_as_one_closed_mesh() {
    for (radius, recipe) in [
        ("2", recipe("2")),
        ("5/2", recipe("5/2")),
        ("placed 2", placed("2")),
        ("placed 5/2", placed("5/2")),
    ] {
        let options = MeshOptions {
            chord_tolerance: Rational::from_ratio(1, 50).unwrap(),
            angular_deviation: Some(Rational::from_ratio(1, 2).unwrap()),
            ..Default::default()
        };
        let mut dec = decider();
        let app = build_certified(&recipe, Some(&options), &mut dec).unwrap();
        let mesh = &app.display().unwrap().as_ref().unwrap().value;
        assert_eq!(mesh.faces.len(), 22);
        assert_eq!(mesh.edges.len(), 41);
        for face in &mesh.faces {
            assert!(!face.triangles.is_empty(), "{radius}: {face:?}");
            assert!(face.max_chord_error <= options.chord_tolerance);
        }
        // Every undirected chord meets exactly two triangles, in opposite
        // directions, including the rows between envelopes and their neighbours.
        let mut uses = BTreeMap::<(usize, usize), (u32, i32)>::new();
        for triangle in &mesh.triangles {
            for i in 0..3 {
                let (a, c) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                let entry = uses.entry((a.min(c), a.max(c))).or_default();
                entry.0 += 1;
                entry.1 += if a < c { 1 } else { -1 };
            }
        }
        assert!(uses.values().all(|&use_| use_ == (2, 0)));
        assert_eq!(
            mesh.unavailable_controls
                .iter()
                .filter(|control| control.control == MeshControl::AngularDeviation)
                .count(),
            4,
            "only the fitted faces lack a normal certificate"
        );
        let built = app.built().built();
        let names = app.names();
        let blends: Vec<_> = names
            .faces
            .iter()
            .filter(|(_, origins)| {
                origins
                    .iter()
                    .any(|origin: &Origin| origin.step.as_str() == "filleted")
            })
            .collect();
        assert_eq!(blends.len(), 8);
        for (_, origins) in &blends {
            assert_eq!(origins, &blends[0].1);
        }
        assert!(
            blends[0]
                .1
                .iter()
                .any(|origin| matches!(origin.role, Role::Blend { .. })),
        );
        let preview = preview_mesh_body(
            &built.model,
            built.result,
            &PreviewOptions::default(),
            &mut dec,
        )
        .unwrap();
        assert_eq!(preview.faces.len(), 22);
        assert_eq!(preview.edges.len(), 41);
    }
}

#[test]
fn the_separated_junction_displays_one_blend_region() {
    let options = MeshOptions {
        chord_tolerance: Rational::from_ratio(1, 50).unwrap(),
        angular_deviation: Some(Rational::from_ratio(1, 2).unwrap()),
        ..Default::default()
    };
    let mut dec = decider();
    let body = catallax_app::build(&recipe("999/1000"), &options, &mut dec).unwrap();
    let mesh = &body.display.value;
    assert_eq!(mesh.faces.len(), 15);
    // 38 loops of one ring each, every ring shared by two faces.
    assert_eq!(mesh.edges.len(), 19);
    for face in &mesh.faces {
        assert!(!face.triangles.is_empty(), "{face:?}");
        assert!(face.max_chord_error <= options.chord_tolerance);
    }
    let mut uses = BTreeMap::<(usize, usize), (u32, i32)>::new();
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let (a, c) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let entry = uses.entry((a.min(c), a.max(c))).or_default();
            entry.0 += 1;
            entry.1 += if a < c { 1 } else { -1 };
        }
    }
    assert!(uses.values().all(|&use_| use_ == (2, 0)));
    assert!(
        mesh.unavailable_controls.is_empty(),
        "{:?}",
        mesh.unavailable_controls
    );
    assert_eq!(body.names.faces.len(), 15);
    let blends: Vec<_> = body
        .names
        .faces
        .iter()
        .filter(|(_, origins)| {
            origins
                .iter()
                .any(|origin: &Origin| matches!(origin.role, Role::Blend { .. }))
        })
        .collect();
    assert_eq!(blends.len(), 1);
    // An exact body has no fitted region: the opt-in route claims nothing.
    let refusal = build_certified(&recipe("999/1000"), None, &mut decider()).unwrap_err();
    assert!(
        format!("{refusal:?}").contains("features.unsupported-offset-family"),
        "{refusal:?}"
    );
}

#[test]
fn the_tangent_junction_displays_one_closed_mesh() {
    let options = MeshOptions {
        chord_tolerance: Rational::from_ratio(1, 50).unwrap(),
        angular_deviation: Some(Rational::from_ratio(1, 2).unwrap()),
        ..Default::default()
    };
    let mut dec = decider();
    let body = catallax_app::build(&recipe("1"), &options, &mut dec).unwrap();
    let mesh = &body.display.value;
    assert!(
        mesh.unavailable_controls.is_empty(),
        "{:?}",
        mesh.unavailable_controls
    );
    // Every patch is a selectable name with its own certified triangles.
    assert_eq!(mesh.faces.len(), body.names.faces.len());
    for face in &mesh.faces {
        assert!(!face.triangles.is_empty(), "{face:?}");
        assert!(face.max_chord_error <= options.chord_tolerance);
        assert!(!body.names.faces[&face.face].is_empty());
    }
    for edge in &mesh.edges {
        assert!(edge.max_chord_error <= options.chord_tolerance);
    }
    let blends = body
        .names
        .faces
        .values()
        .filter(|origins| {
            origins
                .iter()
                .any(|origin: &Origin| matches!(origin.role, Role::Blend { .. }))
        })
        .count();
    assert_eq!(blends, 1);
    let mut uses = BTreeMap::<(usize, usize), (u32, i32)>::new();
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let (a, c) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let entry = uses.entry((a.min(c), a.max(c))).or_default();
            entry.0 += 1;
            entry.1 += if a < c { 1 } else { -1 };
        }
    }
    assert!(uses.values().all(|&use_| use_ == (2, 0)));
    // Each contact is one mesh vertex whose triangles close one fan through
    // the torus, both plane cusps and the cone.
    for contact in [
        [13.0, 5.0, 0.0],
        [0.0, 5.0, 13.0],
        [-13.0, 5.0, 0.0],
        [0.0, 5.0, -13.0],
    ] {
        let found: Vec<_> = mesh
            .positions
            .iter()
            .enumerate()
            .filter(|(_, p)| p.iter().zip(contact).all(|(a, b)| (a - b).abs() < 1e-9))
            .map(|(i, _)| i)
            .collect();
        let [vertex] = found[..] else {
            panic!("one mesh vertex at {contact:?}: {found:?}");
        };
        let mut link = BTreeMap::new();
        let mut around = BTreeSet::new();
        for triangle in &mesh.triangles {
            if let Some(i) = triangle.vertices.iter().position(|v| *v == vertex) {
                let next = triangle.vertices[(i + 1) % 3];
                let last = triangle.vertices[(i + 2) % 3];
                assert!(link.insert(next, last).is_none());
                around.insert(triangle.face);
            }
        }
        let first = *link.keys().next().unwrap();
        let (mut at, mut steps) = (link[&first], 1);
        while at != first {
            (at, steps) = (link[&at], steps + 1);
        }
        assert_eq!(steps, link.len(), "one fan at {contact:?}");
        assert_eq!(around.len(), 4, "{contact:?}");
    }
}

#[test]
fn the_ordinary_application_build_displays_the_same_faces() {
    let mut dec = decider();
    let options = MeshOptions {
        chord_tolerance: Rational::from_ratio(1, 50).unwrap(),
        ..Default::default()
    };
    let body: ApplicationBody = catallax_app::build(&recipe("2"), &options, &mut dec).unwrap();
    assert_eq!(body.display.value.faces.len(), 22);
}
