//! The full P5 housing through the certified application route: one closed
//! mesh across every boundary, per-face bounds and stable picking names.

use std::collections::BTreeMap;

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
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    })
}

#[test]
fn the_mixed_body_displays_as_one_closed_mesh() {
    for radius in ["2", "5/2"] {
        let options = MeshOptions {
            chord_tolerance: Rational::from_ratio(1, 50).unwrap(),
            angular_deviation: Some(Rational::from_ratio(1, 2).unwrap()),
            ..Default::default()
        };
        let mut dec = decider();
        let app = build_certified(&recipe(radius), Some(&options), &mut dec).unwrap();
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
fn the_ordinary_application_build_displays_the_same_faces() {
    let mut dec = decider();
    let options = MeshOptions {
        chord_tolerance: Rational::from_ratio(1, 50).unwrap(),
        ..Default::default()
    };
    let body: ApplicationBody = catallax_app::build(&recipe("2"), &options, &mut dec).unwrap();
    assert_eq!(body.display.value.faces.len(), 22);
}
