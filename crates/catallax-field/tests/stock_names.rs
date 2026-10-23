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
                .visit_candidates(2, |point, id| {
                    let origin = &canonical_field.faces()[id.index()].origin;
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
    let expansion = CanonicalStock::new(&operation, &step).unwrap().unwrap();
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

#[test]
fn canonical_expansion_cannot_turn_impossible_stock_dimensions_into_legal_sweeps() {
    for fields in [
        json!({"op":"block","size":["-1","2","3"]}),
        json!({"op":"block","size":["1","0","3"]}),
        json!({"op":"block","size":["1","2","-3"]}),
        json!({"op":"cylinder","radius":"-1","height":"3"}),
        json!({"op":"cylinder","radius":"1","height":"0"}),
        json!({"op":"cylinder","radius":"1","height":"-3"}),
        json!({"op":"sphere","radius":"-3"}),
        json!({"op":"cone","base_radius":"0","top_radius":"1","height":"4"}),
        json!({"op":"cone","base_radius":"3","top_radius":"-1","height":"4"}),
        json!({"op":"cone","base_radius":"3","top_radius":"1","height":"0"}),
        json!({"op":"cone","base_radius":"3","top_radius":"1","height":"-4"}),
        json!({"op":"torus","major":"6","minor":"0"}),
        json!({"op":"torus","major":"1","minor":"2"}),
        json!({"op":"torus","major":"-6","minor":"2"}),
    ] {
        let mut value = fields;
        value["pose"] = json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]});
        let operation: Operation = serde_json::from_value(value.clone()).unwrap();
        let step = StepId::new("stock").unwrap();
        let recipe = Recipe::new(
            "invalid stock",
            vec![Step {
                id: step.clone(),
                op: operation.clone(),
            }],
            step.clone(),
        );
        assert!(
            evaluate(&recipe, &mut Decider::default()).is_err(),
            "{value}"
        );
        assert!(Plan::compile(&recipe).is_err(), "{value}");
        let error =
            evaluate_canonical_stock(&operation, &step, &mut Decider::default()).unwrap_err();
        assert!(
            matches!(error, catallax_features::BuildError::Geometry { refusal, .. }
            if refusal.code == "recipe.impossible-stock"),
            "{value}"
        );
        assert_eq!(
            Plan::compile_canonical_stock(&operation, &step, &mut Decider::default())
                .unwrap_err()
                .code,
            "recipe.impossible-stock",
            "{value}"
        );
    }
}

#[test]
fn pattern_and_nested_mirror_copies_resolve_canonical_stock_names() {
    use catallax_features::{copy_origin, mirror_placement, pattern_placement};
    use catallax_recipe::{EdgeName, Origin, PatternLayout, Plane, Point3, Vector3};
    use std::collections::BTreeSet;

    let seed = StepId::new("seed").unwrap();
    let row = StepId::new("row").unwrap();
    let image = StepId::new("image").unwrap();
    let layout = PatternLayout::Linear {
        vector: Vector3([30.into(), 0.into(), 0.into()]),
        count: 3,
    };
    let plane = Plane {
        point: Point3([0.into(), 0.into(), 0.into()]),
        normal: Vector3([1.into(), 0.into(), 0.into()]),
    };
    for mut value in [
        json!({"op":"block","size":["4","6","8"]}),
        json!({"op":"cylinder","radius":"3","height":"4"}),
        json!({"op":"cone","base_radius":"3","top_radius":"1","height":"4"}),
        json!({"op":"sphere","radius":"3"}),
        json!({"op":"torus","major":"6","minor":"2"}),
    ] {
        value["pose"] = json!({"rotation":["1","0","0","0"],"translation":["12","0","0"]});
        let operation: Operation = serde_json::from_value(value.clone()).unwrap();
        let canonical = evaluate_canonical_stock(&operation, &seed, &mut Decider::default())
            .unwrap()
            .unwrap();
        let source_field =
            Plan::compile_canonical_stock(&operation, &seed, &mut Decider::default())
                .unwrap()
                .unwrap();
        for expanded in [false, true] {
            // Cone cap correspondence belongs to the canonical entry, not to
            // generic authored radial segments. Its copied names are still
            // checked against that entry's result below.
            if expanded && matches!(operation, Operation::Cone { .. }) {
                continue;
            }
            for mirrored in [false, true] {
                let mut recipe = if expanded {
                    CanonicalStock::new(&operation, &seed)
                        .unwrap()
                        .unwrap()
                        .recipe()
                        .clone()
                } else {
                    Recipe::new(
                        "copied stock",
                        vec![Step {
                            id: seed.clone(),
                            op: operation.clone(),
                        }],
                        seed.clone(),
                    )
                };
                recipe.steps.push(Step {
                    id: row.clone(),
                    op: Operation::Pattern {
                        body: seed.clone(),
                        layout: layout.clone(),
                    },
                });
                recipe.result = row.clone();
                if mirrored {
                    recipe.steps.push(Step {
                        id: image.clone(),
                        op: Operation::Mirror {
                            body: row.clone(),
                            plane: plane.clone(),
                            keep_original: false,
                        },
                    });
                    recipe.result = image.clone();
                }
                let built = evaluate(&recipe, &mut Decider::default()).unwrap();
                let field = Plan::compile(&recipe).unwrap();
                for instance in 0..3 {
                    let copied = |origin: &Origin| {
                        let origin = copy_origin(&row, instance, origin);
                        if mirrored {
                            copy_origin(&image, 1, &origin)
                        } else {
                            origin
                        }
                    };
                    let mut map =
                        pattern_placement(&layout, instance, &mut Decider::default()).unwrap();
                    if mirrored {
                        map = map.then(&mirror_placement(&plane).unwrap());
                    }
                    let mut expected = BTreeSet::new();
                    for face in canonical.model.faces_of_body(canonical.result) {
                        for origin in canonical.model.face_origins(face) {
                            let name = copied(origin);
                            assert_eq!(
                                built.model.faces_with_origin(built.result, &name).len(),
                                1,
                                "{value}: {name:?}, expanded={expanded}, mirrored={mirrored}"
                            );
                            expected.insert(name);
                        }
                    }
                    for edge in canonical.model.edges_of_body(canonical.result) {
                        let name = canonical.model.edge_name(canonical.result, edge).unwrap();
                        let [a, b] = name.between();
                        let name = EdgeName::new(copied(a), copied(b)).unwrap();
                        assert_eq!(
                            built.model.edges_with_name(built.result, &name).len(),
                            1,
                            "{value}: {name:?}, expanded={expanded}, mirrored={mirrored}"
                        );
                    }
                    let mut witnessed = BTreeSet::new();
                    source_field
                        .visit_candidates(1, |point, id| {
                            let origin = &source_field.faces()[id.index()].origin;
                            if source_field
                                .face_query(id, &point, &mut Workspace::default())
                                .unwrap()
                                .state
                                == FaceState::OnFace
                            {
                                let name = copied(origin);
                                let placed = map.apply_rational_point(&point).unwrap();
                                assert_eq!(
                                    field
                                        .face_query(
                                            field.select(&name).unwrap(),
                                            &placed,
                                            &mut Workspace::default()
                                        )
                                        .unwrap()
                                        .state,
                                    FaceState::OnFace,
                                    "{value}: {name:?}, expanded={expanded}, mirrored={mirrored}"
                                );
                                if mirrored {
                                    assert_eq!(
                                        field
                                            .face_query(
                                                field.select(origin).unwrap(),
                                                &placed,
                                                &mut Workspace::default()
                                            )
                                            .unwrap()
                                            .state,
                                        FaceState::NotOnFace
                                    );
                                }
                                witnessed.insert(name);
                            }
                            true
                        })
                        .unwrap();
                    assert_eq!(
                        witnessed, expected,
                        "{value}: every copied origin needs a boundary witness"
                    );
                }
                if expanded {
                    // The final feature-sampling pass must not publish aliases
                    // belonging only to the seed or an unkept mirror input.
                    let mut sampled = BTreeSet::new();
                    field
                        .visit_candidates(1, |_, id| {
                            let origin = &field.faces()[id.index()].origin;
                            assert_eq!(
                                built.model.faces_with_origin(built.result, origin).len(),
                                1,
                                "{value}: sampled an unused source alias {origin:?}"
                            );
                            sampled.insert(origin.clone());
                            true
                        })
                        .unwrap();
                    let expected: BTreeSet<_> = built
                        .model
                        .faces_of_body(built.result)
                        .into_iter()
                        .flat_map(|face| built.model.face_origins(face).iter().cloned())
                        .collect();
                    assert_eq!(
                        sampled, expected,
                        "{value}: sample every named result face, not its unused aliases"
                    );
                }
            }
        }
    }
}
