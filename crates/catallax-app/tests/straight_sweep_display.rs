//! The first user-facing sweep has bounded display at application controls.

use std::collections::BTreeMap;

use catallax_mesh::MeshOptions;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::Recipe;

#[test]
fn rectangular_oblique_sweep_is_closed_and_bounded_at_soliddraft_controls() {
    let source = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/straight-sweep.json"
    ))
    .unwrap();
    for options in [
        MeshOptions::default(),
        MeshOptions {
            chord_tolerance: Rational::from_ratio(1, 50).unwrap(),
            angular_deviation: Some(Rational::from_ratio(1, 2).unwrap()),
            ..Default::default()
        },
    ] {
        let app = catallax_app::build(&source, &options, &mut Decider::default()).unwrap();
        let mesh = &app.display.value;
        assert!(app.display.marks.is_empty());
        assert!(mesh.unavailable_controls.is_empty());
        let mut uses = BTreeMap::<_, (usize, i32)>::new();
        for triangle in &mesh.triangles {
            for i in 0..3 {
                let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                let entry = uses.entry((a.min(b), a.max(b))).or_default();
                entry.0 += 1;
                entry.1 += if a < b { 1 } else { -1 };
            }
        }
        assert!(uses.values().all(|&entry| entry == (2, 0)));
        assert_eq!(mesh.faces.len(), 6);
        for face in &mesh.faces {
            assert!(face.max_chord_error <= options.chord_tolerance);
            assert_eq!(face.max_angular_deviation, options.angular_deviation);
            assert!(!face.triangles.is_empty());
            assert_eq!(app.names.faces[&face.face].len(), 1);
        }
    }
}
