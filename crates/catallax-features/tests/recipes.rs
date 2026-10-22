//! End-to-end recipe construction, through-all bounds and boolean provenance.

use catallax_features::{Built, evaluate};
use catallax_numbers::{Ball, Decider, Outcome, Rational};
use catallax_recipe::{ElementId, Origin, Recipe, Role, StepId, VertexName};
use serde_json::{Value, json};

fn rectangle(id: &str, x: i64, y: i64, width: i64, height: i64, z: i64) -> Value {
    json!({"id": id, "op": "sketch", "plane": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", z.to_string()]}, "loops": [[
        {"kind":"line", "id":"bottom", "to":[(x+width).to_string(), y.to_string()]},
        {"kind":"line", "id":"right", "to":[(x+width).to_string(), (y+height).to_string()]},
        {"kind":"line", "id":"top", "to":[x.to_string(), (y+height).to_string()]},
        {"kind":"line", "id":"left", "to":[x.to_string(), y.to_string()]}
    ]]})
}
fn build(steps: Vec<Value>, result: &str) -> Built {
    let text = json!({"format": catallax_recipe::FORMAT, "name": "feature", "steps": steps, "result":result}).to_string();
    let recipe = Recipe::parse(&text).unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(
        built.model.check(built.result).is_empty(),
        "{:?}",
        built.model.check(built.result)
    );
    built
}
fn mass(built: &Built) -> catallax_measure::MassProperties {
    let answer = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(answer.marks.is_empty(), "{:?}", answer.marks);
    answer.value
}
fn export(built: &Built) -> Vec<u8> {
    let mut bytes = Vec::new();
    let answer = catallax_step::write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    assert!(answer.marks.is_empty(), "{:?}", answer.marks);
    bytes
}

#[test]
fn merged_coplanar_faces_keep_both_creating_steps() {
    let built = build(
        vec![
            rectangle("a_profile", 0, 0, 4, 2, 0),
            json!({"id":"a", "op":"extrude", "profile":"a_profile", "distance":"3"}),
            rectangle("b_profile", 4, 0, 2, 2, 0),
            json!({"id":"b", "op":"extrude", "profile":"b_profile", "distance":"3"}),
            json!({"id":"joined", "op":"boolean", "kind":"union", "target":"a", "tool":"b"}),
        ],
        "joined",
    );
    let mass = mass(&built);
    assert_eq!(mass.volume.lo, Rational::from_i64(36));
    assert_eq!(mass.volume.hi, Rational::from_i64(36));
    assert_eq!(built.model.counts(built.result).faces, 6);
    let top = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .find(|f| {
            built
                .model
                .face_origins(*f)
                .iter()
                .any(|o| o.step.as_str() == "a" && o.role == Role::End)
        })
        .unwrap();
    let origins = built.model.face_origins(top);
    assert_eq!(origins.len(), 2);
    assert!(
        origins
            .iter()
            .any(|o| o.step.as_str() == "b" && o.role == Role::End)
    );
    assert!(origins.windows(2).all(|w| w[0].key() < w[1].key()));
    assert_eq!(
        built
            .model
            .counts(built.bodies[&catallax_recipe::StepId::new("a").unwrap()])
            .faces,
        6
    );
}

#[test]
fn through_all_cut_uses_target_extents_and_keeps_only_side_origins() {
    for z in [-100, 0, 100] {
        let built = build(
            vec![
                rectangle("base_profile", 0, 0, 10, 10, z),
                json!({"id":"base", "op":"extrude", "profile":"base_profile", "distance":"4"}),
                json!({"id":"drill", "op":"sketch", "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}, "loops":[[{"kind":"circle", "id":"bore", "centre":["5","5"], "radius":"1"}]]}),
                json!({"id":"cut", "op":"extrude_through", "target":"base", "profile":"drill"}),
            ],
            "cut",
        );
        let mass = mass(&built);
        let volume = Ball::from_rational(&Rational::from_i64(400), 192)
            - Ball::from_rational(&Rational::from_i64(4), 192) * Ball::pi(192);
        assert!(mass.volume.lo <= volume.lo() && mass.volume.hi >= volume.hi());
        let area = Ball::from_rational(&Rational::from_i64(360), 192)
            + Ball::from_rational(&Rational::from_i64(6), 192) * Ball::pi(192);
        assert!(mass.area.lo <= area.lo() && mass.area.hi >= area.hi());
        for (value, expected) in mass.centroid.unwrap().iter().zip([5, 5, z + 2]) {
            assert!(value.contains(&Rational::from_i64(expected)));
        }
        let cut_origins: Vec<_> = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .flat_map(|f| built.model.face_origins(f))
            .filter(|o| o.step.as_str() == "cut")
            .collect();
        assert_eq!(cut_origins.len(), 1);
        assert!(
            matches!(&cut_origins[0].role, Role::Side { segment, path: None } if segment.as_str() == "bore")
        );
    }
}

#[test]
fn repeated_recipe_builds_have_identical_step_bytes_and_origin_sets() {
    let steps = vec![
        rectangle("profile", 0, 0, 4, 2, 0),
        json!({"id":"body", "op":"extrude", "profile":"profile", "distance":"-3"}),
    ];
    let a = build(steps.clone(), "body");
    let expected_step = export(&a);
    let expected_mass = mass(&a);
    let names = |b: &Built| {
        b.model
            .faces_of_body(b.result)
            .into_iter()
            .map(|f| b.model.face_origins(f).to_vec())
            .collect::<Vec<_>>()
    };
    let expected_names = names(&a);
    for count in [1, 2, 4] {
        std::thread::scope(|scope| {
            let handles: Vec<_> = (0..count)
                .map(|_| {
                    let steps = &steps;
                    scope.spawn(move || {
                        let built = build(steps.clone(), "body");
                        (export(&built), mass(&built), names(&built))
                    })
                })
                .collect();
            for handle in handles {
                let (step, properties, origins) = handle.join().unwrap();
                assert_eq!(step, expected_step);
                assert_eq!(properties, expected_mass);
                assert_eq!(origins, expected_names);
            }
        });
    }
}

#[test]
fn subtracting_a_blind_slot_splits_faces_without_renaming_them() {
    let built = build(
        vec![
            rectangle("base_profile", 0, 0, 10, 4, 0),
            json!({"id":"base", "op":"extrude", "profile":"base_profile", "distance":"2"}),
            rectangle("tool_profile", 4, -1, 2, 6, 1),
            json!({"id":"tool", "op":"extrude", "profile":"tool_profile", "distance":"2"}),
            json!({"id":"cut", "op":"boolean", "kind":"subtract", "target":"base", "tool":"tool"}),
        ],
        "cut",
    );
    let properties = mass(&built);
    assert!(properties.volume.contains(&Rational::from_i64(72)));
    assert!(properties.area.contains(&Rational::from_i64(140)));
    assert!(properties.volume.width() < Rational::parse("1e-9").unwrap());
    assert!(properties.area.width() < Rational::parse("1e-9").unwrap());
    let tops: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|f| {
            built
                .model
                .face_origins(*f)
                .iter()
                .any(|o| o.step.as_str() == "base" && o.role == Role::End)
        })
        .collect();
    assert_eq!(tops.len(), 2);
    assert_eq!(
        built.model.face_origins(tops[0]),
        built.model.face_origins(tops[1])
    );
    assert!(
        built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .any(|f| built
                .model
                .face_origins(f)
                .iter()
                .any(|o| o.step.as_str() == "tool" && o.role == Role::Start))
    );
}

#[test]
fn consecutive_collinear_segments_merge_their_side_origins() {
    let mut sketch = rectangle("profile", 0, 0, 4, 2, 0);
    sketch["loops"][0].as_array_mut().unwrap().insert(
        0,
        json!({"kind":"line", "id":"bottom_left", "to":["2","0"]}),
    );
    let built = build(
        vec![
            sketch,
            json!({"id":"body", "op":"extrude", "profile":"profile", "distance":"3"}),
        ],
        "body",
    );
    assert_eq!(built.model.counts(built.result).faces, 6);
    assert_eq!(mass(&built).volume.lo, Rational::from_i64(24));
    let merged = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .find(|f| built.model.face_origins(*f).len() == 2)
        .unwrap();
    let ids: Vec<_> = built
        .model
        .face_origins(merged)
        .iter()
        .map(|o| match &o.role {
            Role::Side { segment, .. } => segment.as_str(),
            _ => panic!("merged cap"),
        })
        .collect();
    assert_eq!(ids, ["bottom", "bottom_left"]);
}

#[test]
fn coincident_boolean_faces_keep_all_origins_without_double_geometry() {
    for kind in ["union", "intersect"] {
        let built = build(
            vec![
                rectangle("profile", 0, 0, 4, 2, 0),
                json!({"id":"a", "op":"extrude", "profile":"profile", "distance":"3"}),
                json!({"id":"b", "op":"extrude", "profile":"profile", "distance":"3"}),
                json!({"id":"result", "op":"boolean", "kind":kind, "target":"a", "tool":"b"}),
            ],
            "result",
        );
        assert!(mass(&built).volume.contains(&Rational::from_i64(24)));
        assert_eq!(built.model.counts(built.result).faces, 6);
        for face in built.model.faces_of_body(built.result) {
            let origins = built.model.face_origins(face);
            assert_eq!(origins.len(), 2, "{kind}: {origins:?}");
            assert_eq!(origins[0].step.as_str(), "a");
            assert_eq!(origins[1].step.as_str(), "b");
            assert_eq!(origins[0].role, origins[1].role);
        }
    }
}

#[test]
fn existing_stock_body_intersects_a_named_extrusion() {
    let built = build(
        vec![
            json!({"id":"stock", "op":"block", "size":["4","2","3"], "pose":{"rotation":["1","0","0","0"], "translation":["0","0","0"]}}),
            rectangle("profile", 2, 0, 4, 2, 0),
            json!({"id":"feature", "op":"extrude", "profile":"profile", "distance":"3"}),
            json!({"id":"result", "op":"boolean", "kind":"intersect", "target":"stock", "tool":"feature"}),
        ],
        "result",
    );
    let properties = mass(&built);
    assert!(properties.volume.contains(&Rational::from_i64(12)));
    assert!(properties.area.contains(&Rational::from_i64(32)));
    for (value, expected) in properties.centroid.unwrap().iter().zip([
        Rational::from_i64(3),
        Rational::from_i64(1),
        Rational::from_ratio(3, 2).unwrap(),
    ]) {
        assert!(value.contains(&expected));
    }
    for face in built.model.faces_of_body(built.result) {
        assert!(!built.model.face_origins(face).is_empty());
    }
}

#[test]
fn cone_apex_names_survive_singular_loops_and_split_edge_fans() {
    let origin = |step: &str, role| Origin {
        step: StepId::new(step).unwrap(),
        role,
    };
    let side = |name: &str| Role::Side {
        segment: ElementId::new(name).unwrap(),
        path: None,
    };
    for (tool, tool_side, split_apex) in [
        (
            json!({"id":"tool", "op":"cylinder", "radius":"1", "height":"4", "pose":{"rotation":["1","0","0","0"], "translation":["0","0","-1"]}}),
            "wall",
            false,
        ),
        (
            json!({"id":"tool", "op":"block", "size":["4","6","4"], "pose":{"rotation":["1","0","0","0"], "translation":["0","-3","-1"]}}),
            "x_min",
            true,
        ),
    ] {
        let built = build(
            vec![
                json!({"id":"cone", "op":"cone", "base_radius":"2", "top_radius":"0", "height":"2", "pose":{"rotation":["1","0","0","0"], "translation":["0","0","0"]}}),
                tool,
                json!({"id":"result", "op":"boolean", "kind":"intersect", "target":"cone", "tool":"tool"}),
            ],
            "result",
        );
        let source = built.bodies[&StepId::new("cone").unwrap()];
        let source_apex = built.model.vertices_of_body(source)[0];
        let point = built.model.vertex(source_apex).point;
        let mut dec = Decider::default();
        let apex = built
            .model
            .vertices_of_body(built.result)
            .into_iter()
            .find(|vertex| {
                matches!(
                    built.model.geom().same_point(
                        point,
                        built.model.vertex(*vertex).point,
                        &mut dec,
                    ),
                    Outcome::ProvenCoincident(_)
                )
            })
            .unwrap();
        let wall = origin("cone", side("wall"));
        let cut = origin("tool", side(tool_side));
        let mut incident = vec![wall.clone()];
        if split_apex {
            incident.push(cut.clone());
        }
        assert_eq!(
            built.model.vertex_name(built.result, apex),
            Some(VertexName::new(incident).unwrap())
        );
        let actual = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .flat_map(|face| {
                let origins = built.model.face_origins(face);
                assert_eq!(origins.len(), 1);
                origins.iter().cloned()
            })
            .collect::<std::collections::BTreeSet<_>>();
        let expected = [origin("cone", Role::Start), wall, cut]
            .into_iter()
            .collect();
        assert_eq!(actual, expected);
        let volume = Ball::from_rational(&Rational::parse("4/3").unwrap(), 192) * Ball::pi(192);
        let properties = mass(&built);
        assert!(properties.volume.lo <= volume.lo() && properties.volume.hi >= volume.hi());
    }
}
