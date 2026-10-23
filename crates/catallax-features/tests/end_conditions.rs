//! Versioned one-sided tools and their real start cap.
use catallax_features::evaluate;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Recipe, Role};
use serde_json::json;

#[test]
fn one_sided_cuts_preserve_half_the_material_and_name_the_sketch_cap_start() {
    for side in ["forward", "backward"] {
        let recipe=Recipe::parse(&json!({"format":"catallax-recipe/4","name":"half cut","steps":[
            {"id":"stock","op":"block","size":["10","10","10"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","-5"]}},
            {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[
                {"kind":"line","id":"a","to":["6","4"]},{"kind":"line","id":"b","to":["6","6"]},{"kind":"line","id":"c","to":["4","6"]},{"kind":"line","id":"d","to":["4","4"]}]]},
            {"id":"cut","op":"extrude_through","target":"stock","profile":"profile","side":side}],"result":"cut"}).to_string()).unwrap();
        let built = evaluate(&recipe, &mut Decider::default()).unwrap();
        assert!(built.model.check(built.result).is_empty());
        let volume = catallax_measure::volume(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        let expected = Rational::from_i64(980);
        assert!(volume.value.lo <= expected && expected <= volume.value.hi);
        let caps: Vec<_> = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .filter(|f| {
                built
                    .model
                    .face_origins(*f)
                    .iter()
                    .any(|o| o.step.as_str() == "cut" && o.role == Role::Start)
            })
            .collect();
        assert_eq!(caps.len(), 1);
        assert!(
            !built
                .model
                .faces_of_body(built.result)
                .iter()
                .flat_map(|f| built.model.face_origins(*f))
                .any(|o| o.step.as_str() == "cut" && o.role == Role::End)
        );
    }
}

fn capped_box(op: &str, profile_z: &str, limit: Option<serde_json::Value>) -> Recipe {
    let mut feature = json!({"id":"cap","op":op,"target":"stock","profile":"profile"});
    if let Some(limit) = limit {
        feature["limit"] = limit;
    }
    Recipe::parse(&json!({"format":"catallax-recipe/4","name":"finite cap","steps":[
        {"id":"stock","op":"block","size":["10","10","10"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","5"]}},
        {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0",profile_z]},"loops":[[
            {"kind":"line","id":"a","to":["4","2"]},{"kind":"line","id":"b","to":["4","4"]},{"kind":"line","id":"c","to":["2","4"]},{"kind":"line","id":"d","to":["2","2"]}
        ]]},feature],"result":"cap"}).to_string()).unwrap()
}

#[test]
fn next_uses_first_entry_or_exit_while_to_resolves_the_selected_trimmed_cap() {
    for (op, z, limit, volume) in [
        ("extrude_next", "0", None, 20),
        ("extrude_next", "10", None, 20),
        ("extrude_next", "5", None, 40),
        (
            "extrude_to",
            "0",
            Some(json!({"face":{"step":"stock","role":"end"}})),
            60,
        ),
    ] {
        let built = evaluate(&capped_box(op, z, limit), &mut Decider::default()).unwrap();
        assert!(built.model.check(built.result).is_empty());
        let measured = catallax_measure::volume(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        let volume = Rational::from_i64(volume);
        assert!(
            measured.value.lo <= volume && volume <= measured.value.hi,
            "{op} at {z}: {:?}",
            measured.value
        );
        assert!(
            built
                .model
                .faces_of_body(built.result)
                .iter()
                .flat_map(|f| built.model.face_origins(*f))
                .all(|o| o.step.as_str() == "cap")
        );
    }
}

#[test]
fn closed_trimmed_hits_create_the_interior_curtain_without_target_origins() {
    let recipe = Recipe::parse(include_str!("fixtures/stepped-limit.json")).unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(
        built.model.check(built.result).is_empty(),
        "{:?}",
        built.model.check(built.result)
    );
    let measured = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let volume = Rational::from_ratio(15, 2).unwrap();
    assert!(
        measured.value.lo <= volume && volume <= measured.value.hi,
        "{:?}",
        measured.value
    );
    let ends: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|f| {
            built
                .model
                .face_origins(*f)
                .iter()
                .any(|o| o.role == Role::End)
        })
        .collect();
    assert_eq!(
        ends.len(),
        3,
        "two limiting levels and one discontinuity curtain"
    );
    assert!(
        built
            .model
            .faces_of_body(built.result)
            .iter()
            .flat_map(|f| built.model.face_origins(*f))
            .all(|o| o.step.as_str() == "cap")
    );
}

#[test]
fn next_requires_a_material_crossing_even_on_an_isolated_interior_profile_line() {
    let mut value: serde_json::Value =
        serde_json::from_str(include_str!("fixtures/stepped-limit.json")).unwrap();
    let cap = value["steps"].as_array_mut().unwrap().last_mut().unwrap();
    cap["op"] = json!("extrude_next");
    cap.as_object_mut().unwrap().remove("limit");
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let catallax_features::BuildError::Geometry { refusal, .. } =
        evaluate(&recipe, &mut Decider::default()).unwrap_err()
    else {
        panic!("valid recipe must reach geometric end-condition validation")
    };
    assert_eq!(refusal.code, "extrude.no-finite-cap");
}

fn stepped_value() -> serde_json::Value {
    serde_json::from_str(include_str!("fixtures/stepped-limit.json")).unwrap()
}
fn near_profile(value: &mut serde_json::Value) {
    value["steps"][2]["plane"] = json!({"rotation":["1","1","0","0"],"translation":["0","6","0"]});
    value["steps"][2]["loops"] = json!([[
        {"kind":"line","id":"a","to":["4/5","26/5"]},{"kind":"line","id":"b","to":["4/5","29/5"]},
        {"kind":"line","id":"c","to":["1/5","29/5"]},{"kind":"line","id":"d","to":["1/5","26/5"]}
    ]]);
}

#[test]
fn near_uses_closed_trimmed_distance_and_touching_uses_shared_edges() {
    for reference in [
        json!({"face":{"step":"obstacles","role":"end"},"near":["1/2","2","11/2"]}),
        json!({"face":{"step":"obstacles","role":"end"},"touching":[{"step":"obstacles","role":{"side":{"segment":"near-bottom"}}}]}),
    ] {
        let mut value = stepped_value();
        near_profile(&mut value);
        value["steps"][3]["limit"] = reference;
        let built = evaluate(
            &Recipe::parse(&value.to_string()).unwrap(),
            &mut Decider::default(),
        )
        .unwrap();
        let measured = catallax_measure::volume(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        let expected = Rational::from_ratio(36, 25).unwrap();
        assert!(measured.value.lo <= expected && expected <= measured.value.hi);
    }
    let mut value = stepped_value();
    near_profile(&mut value);
    value["steps"][3]["limit"] =
        json!({"face":{"step":"obstacles","role":"end"},"near":["1","2","8"]});
    let error = evaluate(
        &Recipe::parse(&value.to_string()).unwrap(),
        &mut Decider::default(),
    )
    .unwrap_err();
    let catallax_features::BuildError::Geometry { refusal, .. } = error else {
        panic!("near tie is discovered against geometry")
    };
    assert_eq!(refusal.code, "recipe.ambiguous-reference");
}

#[test]
fn circular_profiles_keep_exact_semicircular_caps_and_their_curtain() {
    let mut value = stepped_value();
    value["steps"][2]["loops"] =
        json!([[{"kind":"circle","id":"rim","centre":["1","3"],"radius":"1/2"}]]);
    let built = evaluate(
        &Recipe::parse(&value.to_string()).unwrap(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(built.model.check(built.result).is_empty());
    let measured = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let expected = catallax_numbers::Ball::pi(128)
        * catallax_numbers::Ball::from_rational(&Rational::from_ratio(15, 8).unwrap(), 128);
    assert!(measured.value.lo <= expected.hi() && expected.lo() <= measured.value.hi);
}

#[test]
fn selected_support_extensions_and_zero_hits_are_not_finite_trimmed_caps() {
    for (z, translate, role) in [("0", "20", "end"), ("5", "0", "start")] {
        let recipe = capped_box(
            "extrude_to",
            z,
            Some(json!({"face":{"step":"stock","role":role}})),
        );
        let mut value: serde_json::Value = serde_json::from_str(&recipe.to_text()).unwrap();
        value["steps"][1]["plane"]["translation"][0] = json!(translate);
        let error = evaluate(
            &Recipe::parse(&value.to_string()).unwrap(),
            &mut Decider::default(),
        )
        .unwrap_err();
        let catallax_features::BuildError::Geometry { refusal, .. } = error else {
            panic!("valid references require a geometric cap proof")
        };
        assert_eq!(refusal.code, "extrude.no-finite-cap");
    }
}

#[test]
fn next_with_a_negative_authored_direction_preserves_the_start_plane() {
    let recipe = capped_box("extrude_next", "20", None);
    let mut value: serde_json::Value = serde_json::from_str(&recipe.to_text()).unwrap();
    value["steps"][2]["direction"] = json!(["0", "0", "-2"]);
    let built = evaluate(
        &Recipe::parse(&value.to_string()).unwrap(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(built.model.check(built.result).is_empty());
    let measured = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let expected = Rational::from_i64(20);
    assert!(measured.value.lo <= expected && expected <= measured.value.hi);
}

#[test]
fn tilted_trimmed_caps_keep_the_authored_sweep_direction() {
    let recipe = capped_box(
        "extrude_to",
        "0",
        Some(json!({"face":{"step":"stock","role":"start"}})),
    );
    let mut value: serde_json::Value = serde_json::from_str(&recipe.to_text()).unwrap();
    value["steps"][0]["pose"]["rotation"] = json!(["2", "0", "1", "0"]);
    value["steps"][1]["loops"] = json!([[
        {"kind":"line","id":"a","to":["3/2","2"]},{"kind":"line","id":"b","to":["3/2","3"]},
        {"kind":"line","id":"c","to":["1/2","3"]},{"kind":"line","id":"d","to":["1/2","2"]}
    ]]);
    let built = evaluate(
        &Recipe::parse(&value.to_string()).unwrap(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(built.model.check(built.result).is_empty());
    let measured = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let expected = Rational::from_ratio(11, 3).unwrap();
    assert!(measured.value.lo <= expected && expected <= measured.value.hi);
}

#[test]
fn empty_profile_needs_only_structural_reference_validation() {
    let recipe=Recipe::parse(&json!({"format":"catallax-recipe/4","name":"empty cap","steps":[
        {"id":"stock","op":"sphere","radius":"2","pose":{"rotation":["1","0","0","0"],"translation":["0","0","5"]}},
        {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[]},
        {"id":"cap","op":"extrude_to","target":"stock","profile":"profile","limit":{"face":{"step":"stock","role":{"side":{"segment":"wall"}}}}}
    ],"result":"cap"}).to_string()).unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert_eq!(
        built.model.body(built.result).kind,
        catallax_topology::BodyKind::Empty
    );
    let mut value: serde_json::Value = serde_json::from_str(&recipe.to_text()).unwrap();
    value["steps"][2]["limit"]["face"]["step"] = json!("missing");
    assert_eq!(
        Recipe::parse(&value.to_string()).unwrap_err().code(),
        "recipe.unknown-origin"
    );
}

#[test]
fn a_zero_height_ridge_separates_touching_lumps_after_regularisation() {
    let recipe=Recipe::parse(&json!({"format":"catallax-recipe/4","name":"zero ridge","steps":[
        {"id":"roof-profile","op":"sketch","plane":{"rotation":["1","1","0","0"],"translation":["0","4","0"]},"loops":[[
            {"kind":"line","id":"left","to":["0","0"]},{"kind":"line","id":"right","to":["2","2"]},{"kind":"line","id":"top","to":["-2","2"]}
        ]]},
        {"id":"roof","op":"extrude","profile":"roof-profile","distance":"2"},
        {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[
            {"kind":"line","id":"a","to":["1","5/2"]},{"kind":"line","id":"b","to":["1","7/2"]},{"kind":"line","id":"c","to":["-1","7/2"]},{"kind":"line","id":"d","to":["-1","5/2"]}
        ]]},
        {"id":"cap","op":"extrude_next","target":"roof","profile":"profile"}
    ],"result":"cap"}).to_string()).unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    assert_eq!(built.model.body(built.result).lumps.len(), 2);
    let measured = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let expected = Rational::one();
    assert!(measured.value.lo <= expected && expected <= measured.value.hi);
}
