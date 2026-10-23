//! Canonical stock equivalence includes provenance, incidence and trimmed domains.

use catallax_features::{Built, evaluate, evaluate_canonical_stock};
use catallax_field::{FaceState, Plan, Workspace};
use catallax_numbers::Decider;
use catallax_recipe::{CanonicalStock, Operation, Recipe, Step, StepId};
use serde_json::json;

fn names(built: &Built) -> (Vec<String>, Vec<String>, Vec<String>) {
    let m = &built.model;
    let b = built.result;
    let mut faces: Vec<_> = m
        .faces_of_body(b)
        .into_iter()
        .map(|f| serde_json::to_string(m.face_origins(f)).unwrap())
        .collect();
    let mut edges: Vec<_> = m
        .edges_of_body(b)
        .into_iter()
        .map(|e| {
            serde_json::to_string(&m.edge_name(b, e).expect("physical edge must be named")).unwrap()
        })
        .collect();
    let mut vertices: Vec<_> = m
        .vertices_of_body(b)
        .into_iter()
        .map(|v| {
            serde_json::to_string(&m.vertex_name(b, v).expect("physical vertex must be named"))
                .unwrap()
        })
        .collect();
    faces.sort();
    edges.sort();
    vertices.sort();
    (faces, edges, vertices)
}

#[test]
fn canonical_features_keep_stock_face_edge_and_vertex_names_in_both_kernels() {
    for rotation in [["1", "0", "0", "0"], ["2", "1", "1", "0"]] {
        for fields in [
            json!({"op":"block", "size":["4","6","8"]}),
            json!({"op":"cylinder", "radius":"3", "height":"4"}),
            json!({"op":"sphere", "radius":"3"}),
            json!({"op":"torus", "major":"6", "minor":"2"}),
            json!({"op":"cone", "base_radius":"3", "top_radius":"1", "height":"4"}),
            json!({"op":"cone", "base_radius":"1", "top_radius":"3", "height":"4"}),
            json!({"op":"cone", "base_radius":"3", "top_radius":"0", "height":"4"}),
            json!({"op":"cone", "base_radius":"3", "top_radius":"3", "height":"4"}),
        ] {
            let mut value = fields;
            value["pose"] = json!({"rotation":rotation,"translation":["7","-2","5"]});
            let operation: Operation = serde_json::from_value(value.clone()).unwrap();
            let step = StepId::new("stock").unwrap();
            let recipe = Recipe::new(
                "stock",
                vec![Step {
                    id: step.clone(),
                    op: operation.clone(),
                }],
                step.clone(),
            );
            let direct = evaluate(&recipe, &mut Decider::default()).unwrap();
            let canonical = evaluate_canonical_stock(&operation, &step, &mut Decider::default())
                .unwrap()
                .unwrap();
            assert!(
                canonical.model.check(canonical.result).is_empty(),
                "{value}"
            );
            assert_eq!(names(&direct), names(&canonical), "{value}");
            let mut geometry = direct.model.geom().clone();
            for face in direct.model.faces_of_body(direct.result) {
                for origin in direct.model.face_origins(face) {
                    let matches = canonical.model.faces_with_origin(canonical.result, origin);
                    assert_eq!(matches.len(), 1, "{value}: {origin:?}");
                    let equivalent = geometry.add_surface(
                        canonical
                            .model
                            .geom()
                            .surface(canonical.model.face(matches[0]).surface)
                            .clone(),
                    );
                    assert!(
                        matches!(
                            geometry.same_surface(
                                direct.model.face(face).surface,
                                equivalent,
                                &mut Decider::default()
                            ),
                            catallax_numbers::Outcome::ProvenCoincident(_)
                        ),
                        "{value}: {origin:?} must name the corresponding exact carrier"
                    );
                }
            }
            let direct_field = Plan::compile(&recipe).unwrap();
            let canonical_field =
                Plan::compile_canonical_stock(&operation, &step, &mut Decider::default())
                    .unwrap()
                    .unwrap();
            let mut expected: Vec<_> = direct_field
                .faces()
                .iter()
                .map(|f| f.origin.clone())
                .collect();
            let mut actual: Vec<_> = canonical_field
                .faces()
                .iter()
                .map(|f| f.origin.clone())
                .collect();
            expected.sort();
            actual.sort();
            assert_eq!(expected, actual, "{value}");
            let mut witnessed = std::collections::BTreeSet::new();
            canonical_field
                .visit_feature_candidates(2, |point, origin| {
                    let id = canonical_field.select(origin).unwrap();
                    if canonical_field
                        .face_query(id, &point, &mut Workspace::default())
                        .unwrap()
                        .state
                        == FaceState::OnFace
                    {
                        assert_eq!(
                            direct_field
                                .face_query(
                                    direct_field.select(origin).unwrap(),
                                    &point,
                                    &mut Workspace::default()
                                )
                                .unwrap()
                                .state,
                            FaceState::OnFace,
                            "{value}: {origin:?} at {point:?}"
                        );
                        witnessed.insert(origin.clone());
                    }
                    true
                })
                .unwrap();
            assert_eq!(
                witnessed.into_iter().collect::<Vec<_>>(),
                expected,
                "every named trimmed domain needs a witness: {value}"
            );
        }
    }
}

#[test]
fn cone_correspondence_does_not_rename_an_authored_radial_segment() {
    let operation: Operation = serde_json::from_value(json!({"op":"cone", "base_radius":"3", "top_radius":"1", "height":"4", "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}})).unwrap();
    let step = StepId::new("canonical_profile").unwrap();
    let expansion = CanonicalStock::new(&operation, &step).unwrap();
    let raw = evaluate(expansion.recipe(), &mut Decider::default()).unwrap();
    let named = evaluate_canonical_stock(&operation, &step, &mut Decider::default())
        .unwrap()
        .unwrap();
    assert_ne!(names(&raw), names(&named));
    let cap = catallax_recipe::Origin {
        step,
        role: catallax_recipe::Role::Start,
    };
    assert!(
        !raw.model
            .faces_of_body(raw.result)
            .iter()
            .any(|f| raw.model.face_origins(*f).contains(&cap))
    );
    assert!(
        named
            .model
            .faces_of_body(named.result)
            .iter()
            .any(|f| named.model.face_origins(*f).contains(&cap))
    );
}
