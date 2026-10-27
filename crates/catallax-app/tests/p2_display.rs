//! Public collar regression: the edited cap must retain the original display controls.
use catallax_mesh::{MeshOptions, PreviewOptions, mesh_body};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use std::collections::BTreeMap;

#[test]
fn original_and_edited_collar_keep_bounded_display_and_named_picking() {
    for distance in ["1", "0.75"] {
        let mut source: serde_json::Value =
            serde_json::from_str(include_str!("../../../parts/P2-collar.json")).unwrap();
        source["steps"][2]["distance"] = distance.into();
        let recipe = Recipe::parse(&source.to_string()).unwrap();
        let (built, work) = with_operation_budget(1_000_000, || {
            catallax_features::evaluate(&recipe, &mut Decider::default())
        });
        let built = built.unwrap();
        assert!(!work.exhausted);
        for (limit, chord, angle, edge) in [
            (1_000_000, "0.01", None, None),
            (20_000_000, "0.00036", Some("0.03490658503988659"), None),
            (20_000_000, "0.01", Some("0.2617993877991494"), Some("1.5")),
        ] {
            let options = MeshOptions {
                chord_tolerance: chord.parse().unwrap(),
                angular_deviation: angle.map(|v| v.parse().unwrap()),
                max_edge_length: edge.map(|v| v.parse().unwrap()),
                ..Default::default()
            };
            let (result, work) = with_operation_budget(limit, || {
                mesh_body(
                    &built.model,
                    built.result,
                    &options,
                    &mut Decider::default(),
                )
            });
            let mesh = result
                .unwrap_or_else(|e| {
                    panic!("distance={distance}, chord={chord}, edge={edge:?}: {e:?}")
                })
                .value;
            assert!(!work.exhausted);
            assert!(mesh.unavailable_controls.is_empty());
            assert_eq!((mesh.faces.len(), mesh.edges.len()), (6, 6));
            let mut incidence = BTreeMap::<_, (usize, i32)>::new();
            for triangle in &mesh.triangles {
                for i in 0..3 {
                    let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                    let entry = incidence.entry((a.min(b), a.max(b))).or_default();
                    entry.0 += 1;
                    entry.1 += if a < b { 1 } else { -1 };
                    if let Some(limit) = &options.max_edge_length {
                        let length = (0..3).fold(Rational::zero(), |sum, j| {
                            let d = Rational::from_f64_exact(mesh.positions[a][j]).unwrap()
                                - Rational::from_f64_exact(mesh.positions[b][j]).unwrap();
                            sum + &d * &d
                        });
                        assert!(length <= limit * limit);
                    }
                }
            }
            assert!(incidence.values().all(|v| *v == (2, 0)));
            for face in &mesh.faces {
                assert!(!face.triangles.is_empty());
                assert!(face.max_chord_error <= options.chord_tolerance);
                assert_eq!(face.max_angular_deviation, options.angular_deviation);
                assert!(!built.model.face_origins(face.face).is_empty());
            }
            for edge in &mesh.edges {
                assert!(edge.max_chord_error <= options.chord_tolerance);
                for pair in edge.vertices.windows(2) {
                    assert_eq!(
                        incidence[&(pair[0].min(pair[1]), pair[0].max(pair[1]))],
                        (2, 0)
                    );
                }
            }
            let (preview, work) = with_operation_budget(limit, || {
                catallax_app::preview_body(
                    &built,
                    &PreviewOptions::from(&options),
                    &mut Decider::default(),
                )
            });
            let preview = preview.unwrap();
            assert!(!work.exhausted);
            assert_eq!(
                (preview.names.faces.len(), preview.names.edges.len()),
                (6, 6)
            );
            assert!(preview.names.faces.values().all(|names| !names.is_empty()));
            assert!(preview.names.edges.values().all(Option::is_some));
        }
    }
}
