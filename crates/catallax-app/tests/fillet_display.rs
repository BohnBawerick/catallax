//! Shared-rail display closure through the application's exact recipe entry.

use std::collections::{BTreeMap, BTreeSet};

use catallax_mesh::MeshOptions;
use catallax_numbers::{BudgetPolicy, Decider};
use catallax_recipe::{Operation, Recipe, Role};

#[test]
fn p3_fillet_display_closes_shared_rails_before_and_after_a_radius_edit() {
    let mut recipe = Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    let mut previous_names = None;
    for radius in ["1", "3/2"] {
        if let Operation::Fillet { radius: value, .. } = &mut recipe.steps[4].op {
            *value = serde_json::from_value(serde_json::json!(radius)).unwrap();
        }
        let options = MeshOptions::default();
        let mut dec = Decider::new(BudgetPolicy {
            max_work: 1_000_000,
            ..Default::default()
        });
        let body = catallax_app::build(&recipe, &options, &mut dec).unwrap();
        let mesh = &body.display.value;
        let expected: BTreeSet<_> = body
            .built
            .model
            .faces_of_body(body.built.result)
            .into_iter()
            .collect();
        let actual: BTreeSet<_> = mesh.triangles.iter().map(|t| t.face).collect();
        assert_eq!(actual, expected);
        let mut incidence = BTreeMap::<_, (usize, i32)>::new();
        for triangle in &mesh.triangles {
            assert!(body.names.faces.contains_key(&triangle.face));
            for i in 0..3 {
                let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                let uses = incidence.entry((a.min(b), a.max(b))).or_default();
                uses.0 += 1;
                uses.1 += if a < b { 1 } else { -1 };
            }
        }
        assert!(incidence.values().all(|uses| *uses == (2, 0)));
        for edge in &mesh.edges {
            assert!(edge.max_chord_error <= options.chord_tolerance);
            for pair in edge.vertices.windows(2) {
                assert_eq!(
                    incidence.get(&(pair[0].min(pair[1]), pair[0].max(pair[1]))),
                    Some(&(2, 0))
                );
            }
        }
        for face in &mesh.faces {
            assert!(face.max_chord_error <= options.chord_tolerance);
        }
        let mut names: Vec<_> = body
            .names
            .faces
            .values()
            .flatten()
            .filter(|origin| matches!(origin.role, Role::Blend { .. }))
            .cloned()
            .collect();
        names.sort();
        assert_eq!(names.len(), 5);
        if let Some(previous) = &previous_names {
            assert_eq!(&names, previous);
        }
        previous_names = Some(names);
    }
}
