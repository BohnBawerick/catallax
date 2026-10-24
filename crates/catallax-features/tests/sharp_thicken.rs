//! Sharp selected joins, straight open rims and independently compiled field signs.

use catallax_features::evaluate;
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Decider, Outcome, Rational, with_operation_budget};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

fn recipe(stock: Value, faces: Vec<Value>, side: &str) -> Recipe {
    Recipe::parse(
        &json!({"format":"catallax-recipe/4","name":"sharp thicken","steps":[
            stock,
            {"id":"result","op":"thicken","body":"stock","faces":faces,"thickness":"1/4","side":side}
        ],"result":"result"})
        .to_string(),
    )
    .unwrap()
}
fn face(role: Value) -> Value {
    json!({"face":{"step":"stock","role":role}})
}
fn stock(kind: &str) -> Value {
    let mut value = match kind {
        "block" => json!({"op":"block","size":["10","8","6"]}),
        "cylinder" => json!({"op":"cylinder","radius":"2","height":"4"}),
        "cone" => json!({"op":"cone","base_radius":"4","top_radius":"1","height":"4"}),
        _ => unreachable!(),
    };
    value["id"] = json!("stock");
    value["pose"] = json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]});
    value
}
fn check(recipe: &Recipe, points: &[([&str; 3], Location)], round: bool) {
    let mut built = evaluate(recipe, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let origins: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|f| built.model.face_origins(f))
        .collect();
    assert_eq!(
        origins.iter().any(|o| matches!(o.role, Role::Round { .. })),
        round
    );
    assert!(origins.iter().any(|o| matches!(o.role, Role::Rim { .. })));
    let plan = Plan::compile(recipe).unwrap();
    let mut workspace = Workspace::default();
    for (point, expected) in points {
        let [x, y, z] = point.map(|s| Rational::parse(s).unwrap());
        let point3 = RPoint3::new(x, y, z);
        let brep = match catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &point3,
            &mut Decider::default(),
        ) {
            Outcome::Proven(catallax_boolean::Containment::Inside) => Location::Inside,
            Outcome::Proven(catallax_boolean::Containment::Outside) => Location::Outside,
            Outcome::ProvenCoincident(_) => Location::Boundary,
            other => panic!("B-rep classification at {point:?}: {other:?}"),
        };
        assert_eq!(brep, *expected, "B-rep at {point:?}");
        let result = plan
            .proven_checked(&point3, &mut workspace, &mut Decider::default())
            .unwrap();
        assert_eq!(result.location, *expected, "at {point:?}: {result:?}");
    }
}

#[test]
fn planar_front_fills_the_selected_edge_without_rounding_open_boundaries() {
    let recipe = recipe(
        stock("block"),
        vec![
            face(json!("end")),
            face(json!({"side":{"segment":"x_max"}})),
        ],
        "front",
    );
    check(
        &recipe,
        &[
            (["203/20", "4", "123/20"], Location::Inside),
            (["51/5", "4", "31/5"], Location::Outside),
            (["203/20", "-1/100", "123/20"], Location::Outside),
            (["5", "4", "61/10"], Location::Inside),
            (["5", "4", "3"], Location::Outside),
        ],
        true,
    );
}

#[test]
fn planar_back_regularises_the_overlap_and_keeps_the_far_face_names() {
    let recipe = recipe(
        stock("block"),
        vec![
            face(json!("end")),
            face(json!({"side":{"segment":"x_max"}})),
        ],
        "back",
    );
    check(
        &recipe,
        &[
            (["99/10", "4", "59/10"], Location::Inside),
            (["5", "4", "59/10"], Location::Inside),
            (["99/10", "4", "3"], Location::Inside),
            (["5", "4", "3"], Location::Outside),
            (["99/10", "-1/100", "59/10"], Location::Outside),
        ],
        false,
    );
}

#[test]
fn circular_front_joins_close_only_the_selected_cap() {
    for (kind, radius) in [("cylinder", "43/20"), ("cone", "23/20")] {
        let recipe = recipe(
            stock(kind),
            vec![face(json!("end")), face(json!({"side":{"segment":"wall"}}))],
            "front",
        );
        check(
            &recipe,
            &[
                ([radius, "0", "83/20"], Location::Inside),
                (["0", "0", "41/10"], Location::Inside),
                (["0", "0", "2"], Location::Outside),
                (["4", "0", "-1/100"], Location::Outside),
            ],
            true,
        );
    }
}

#[test]
fn circular_back_joins_regularise_overlapping_normal_images() {
    for kind in ["cylinder", "cone"] {
        let recipe = recipe(
            stock(kind),
            vec![face(json!("end")), face(json!({"side":{"segment":"wall"}}))],
            "back",
        );
        check(
            &recipe,
            &[
                (["0", "0", "39/10"], Location::Inside),
                (["0", "0", "2"], Location::Outside),
                (["0", "0", "41/10"], Location::Outside),
            ],
            false,
        );
    }
}

#[test]
fn selected_interior_vertex_gets_a_spherical_sector_not_a_cubic_mitre() {
    let recipe = recipe(
        stock("block"),
        vec![
            face(json!("end")),
            face(json!({"side":{"segment":"x_max"}})),
            face(json!({"side":{"segment":"y_max"}})),
        ],
        "front",
    );
    check(
        &recipe,
        &[
            (["101/10", "81/10", "61/10"], Location::Inside),
            (["51/5", "41/5", "31/5"], Location::Outside),
        ],
        true,
    );
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(built.model.faces_of_body(built.result).iter().any(|f| {
        built
            .model
            .face_origins(*f)
            .iter()
            .any(|o| matches!(o.role, Role::RoundCorner { .. }))
    }));
}

#[test]
fn boundary_sector_fans_keep_both_incident_rim_origins() {
    let recipe = recipe(
        stock("block"),
        vec![
            face(json!("end")),
            face(json!({"side":{"segment":"x_max"}})),
        ],
        "front",
    );
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let fans: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|f| {
            let origins = built.model.face_origins(*f);
            origins.len() == 2 && origins.iter().all(|o| matches!(o.role, Role::Rim { .. }))
        })
        .collect();
    assert_eq!(fans.len(), 2);
    for fan in fans {
        assert!(matches!(
            built.model.geom().surface(built.model.face(fan).surface),
            catallax_shapes::Surface::Plane(_)
        ));
    }
}

fn symmetric_base(text: &str, p: RPoint3) {
    let recipe = Recipe::parse(text).unwrap();
    let (built, work) =
        with_operation_budget(100_000_000, || evaluate(&recipe, &mut Decider::default()));
    let built = built.unwrap();
    assert!(!work.exhausted);
    assert!(built.model.check(built.result).is_empty());
    assert!(built.model.faces_of_body(built.result).iter().all(|f| {
        built
            .model
            .face_origins(*f)
            .iter()
            .all(|o| o.step.as_str() != "stock")
    }));
    let plan = Plan::compile(&recipe).unwrap();
    assert_eq!(
        plan.proven_checked(&p, &mut Workspace::default(), &mut Decider::default())
            .unwrap()
            .location,
        Location::Inside
    );
}

#[test]
fn planar_two_sided_union_removes_the_shared_base() {
    symmetric_base(
        include_str!("fixtures/sharp-thicken/planar-symmetric.json"),
        RPoint3::new(5.into(), 4.into(), 6.into()),
    );
}

#[test]
fn cylindrical_two_sided_union_removes_the_shared_base() {
    symmetric_base(
        include_str!("fixtures/sharp-thicken/cylinder-symmetric.json"),
        RPoint3::new(0.into(), 0.into(), 4.into()),
    );
}

#[test]
fn conical_two_sided_union_removes_the_shared_base() {
    symmetric_base(
        include_str!("fixtures/sharp-thicken/cone-symmetric.json"),
        RPoint3::new(0.into(), 0.into(), 4.into()),
    );
}

#[test]
fn conical_tips_and_horn_projection_carriers_still_refuse_cleanly() {
    for (top, thickness, code) in [
        ("0", "1/4", "distance.cone-domain"),
        ("1", "1", "offset.spindle-carrier"),
    ] {
        let mut cone = stock("cone");
        cone["top_radius"] = json!(top);
        let mut text = json!({"format":"catallax-recipe/4","name":"unproved sharp sector","steps":[cone,
            {"id":"result","op":"thicken","body":"stock","faces":[face(json!("end")),face(json!({"side":{"segment":"wall"}}))],"thickness":thickness,"side":"front"}
        ],"result":"result"});
        if top == "0" {
            text["steps"][1]["faces"][0] = face(json!("start"));
        }
        let recipe = Recipe::parse(&text.to_string()).unwrap();
        let catallax_features::BuildError::Geometry { refusal, .. } =
            evaluate(&recipe, &mut Decider::default()).unwrap_err()
        else {
            panic!("expected typed geometry refusal");
        };
        assert_eq!(refusal.code, code);
        assert_eq!(Plan::compile(&recipe).unwrap_err().code, code);
    }
}

#[test]
fn coaxial_cone_cylinder_selection_uses_its_shared_circular_sector() {
    let recipe = Recipe::parse(include_str!(
        "fixtures/sharp-thicken/cone-cylinder-front.json"
    ))
    .unwrap();
    check(
        &recipe,
        &[
            (["269/100", "0", "97/50"], Location::Inside),
            (["11/4", "0", "9/5"], Location::Outside),
            (["0", "0", "2"], Location::Outside),
        ],
        true,
    );
}

#[test]
fn version_three_keeps_the_earlier_sharp_selection_refusal() {
    for kind in ["block", "cylinder", "cone"] {
        let roles = if kind == "block" {
            vec![
                face(json!("end")),
                face(json!({"side":{"segment":"x_max"}})),
            ]
        } else {
            vec![face(json!("end")), face(json!({"side":{"segment":"wall"}}))]
        };
        let mut recipe = recipe(stock(kind), roles, "front");
        recipe.format = catallax_recipe::FORMAT.into();
        let catallax_features::BuildError::Geometry { refusal, .. } =
            evaluate(&recipe, &mut Decider::default()).unwrap_err()
        else {
            panic!("expected legacy geometry refusal");
        };
        assert_eq!(refusal.code, "distance.unsupported-envelope");
        assert_eq!(Plan::compile(&recipe).unwrap_err().code, refusal.code);
    }
}

#[test]
fn closed_selected_skins_retain_the_input_void() {
    for kind in ["block", "cylinder", "cone"] {
        let recipe = recipe(stock(kind), vec![json!({"face":{"step":"stock"}})], "front");
        let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
        assert!(built.model.check(built.result).is_empty());
        let plan = Plan::compile(&recipe).unwrap();
        let top = if kind == "block" { "61/10" } else { "41/10" };
        let centre = if kind == "block" {
            ["5", "4", "3"]
        } else {
            ["0", "0", "2"]
        };
        for (p, expected) in [
            (centre, Location::Outside),
            (["0", "0", top], Location::Inside),
        ] {
            let [x, y, z] = p.map(|s| Rational::parse(s).unwrap());
            let point = RPoint3::new(x, y, z);
            let brep = catallax_boolean::classify_point(
                &mut built.model,
                built.result,
                &point,
                &mut Decider::default(),
            );
            assert!(
                matches!(
                    (&brep, expected),
                    (
                        Outcome::Proven(catallax_boolean::Containment::Inside),
                        Location::Inside
                    ) | (
                        Outcome::Proven(catallax_boolean::Containment::Outside),
                        Location::Outside
                    )
                ),
                "{kind} at {p:?}: {brep:?}"
            );
            assert_eq!(
                plan.proven_checked(&point, &mut Workspace::default(), &mut Decider::default())
                    .unwrap()
                    .location,
                expected,
                "{kind} at {p:?}"
            );
        }
    }
}
