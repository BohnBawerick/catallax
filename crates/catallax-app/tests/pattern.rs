//! Saved instances are ordinary replayable operations with independent names.
use catallax_app::pattern::{FeatureInstance, expand_feature};
use catallax_features::evaluate;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Recipe, Role, StepId};
use catallax_shapes::{Isometry, RPoint3, RVector3};
use serde_json::json;
use std::collections::BTreeMap;
fn id(s: &str) -> StepId {
    StepId::new(s).unwrap()
}
fn source() -> Recipe {
    Recipe::parse(&json!({"format":catallax_recipe::FORMAT,"name":"saved-cuts","steps":[
        {"id":"base","op":"block","size":["12","4","3"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
        {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["1","1","0"]},"loops":[[
            {"kind":"line","id":"a","to":["1","0"]},{"kind":"line","id":"b","to":["1","1"]},{"kind":"line","id":"c","to":["0","1"]},{"kind":"line","id":"d","to":["0","0"]}
        ]]},
        {"id":"cut","op":"extrude_through","target":"base","profile":"profile"}
    ],"result":"cut"}).to_string()).unwrap()
}
fn instance(k: i64, previous: &str) -> FeatureInstance {
    FeatureInstance {
        placement: Isometry::translation(RVector3::new(
            Rational::from_i64(3 * k),
            0.into(),
            0.into(),
        )),
        steps: BTreeMap::from([
            (id("profile"), id(&format!("profile{k}"))),
            (id("cut"), id(&format!("cut{k}"))),
        ]),
        inputs: BTreeMap::from([(id("base"), id(previous))]),
        origins: BTreeMap::new(),
    }
}
#[test]
fn sequential_cuts_round_trip_with_stable_ids_and_independent_origins() {
    let r = source();
    let graph = [id("profile"), id("cut")];
    let one = expand_feature(
        &r,
        &graph,
        &id("cut"),
        &[instance(1, "cut")],
        &mut Decider::default(),
    )
    .unwrap();
    let two = expand_feature(
        &r,
        &graph,
        &id("cut"),
        &[instance(1, "cut"), instance(2, "cut1")],
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(one.steps, two.steps[..one.steps.len()]);
    let saved = Recipe::parse(&two.to_text()).unwrap();
    let built = evaluate(&saved, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let mass = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert!(mass.volume.contains(&Rational::from_i64(135)));
    assert!(mass.volume.width() <= Rational::parse("1e-9").unwrap());
    let origins = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|f| built.model.face_origins(f))
        .collect::<Vec<_>>();
    for name in ["cut", "cut1", "cut2"] {
        assert!(origins.iter().any(|o| o.step == id(name)));
    }
    assert!(origins.iter().all(|o| !matches!(o.role, Role::Copy { .. })));
}
#[test]
fn reflected_feature_keeps_its_original_extrusion_direction() {
    let mut r = source();
    r.steps[2].op = catallax_recipe::Operation::Extrude {
        profile: id("profile"),
        distance: 3.into(),
        direction: None,
    };
    let mut copy = instance(1, "base");
    copy.placement = Isometry::reflection(
        &RPoint3::origin(),
        &RVector3::new(0.into(), 0.into(), 1.into()),
    )
    .unwrap();
    let expanded = expand_feature(
        &r,
        &[id("profile"), id("cut")],
        &id("cut"),
        &[copy],
        &mut Decider::default(),
    )
    .unwrap();
    let built = evaluate(&expanded, &mut Decider::default()).unwrap();
    let mass = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert_eq!(mass.volume.lo, Rational::from_i64(3));
    assert!(mass.centroid.unwrap()[2].contains(&Rational::parse("-3/2").unwrap()));
}
#[test]
fn expansion_rejects_missing_bindings_and_algebraic_parameter_rounding_atomically() {
    let r = source();
    let saved = r.to_text();
    let mut copy = instance(1, "cut");
    copy.inputs.clear();
    assert_eq!(
        expand_feature(
            &r,
            &[id("profile"), id("cut")],
            &id("cut"),
            &[copy],
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "pattern.expansion"
    );
    let mut copy = instance(1, "cut");
    copy.placement = Isometry::around_axis(
        &RPoint3::origin(),
        &RVector3::new(0.into(), 0.into(), 1.into()),
        &Rational::from_i64(45),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(
        expand_feature(
            &r,
            &[id("profile"), id("cut")],
            &id("cut"),
            &[copy],
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "pattern.expansion-exact-parameter"
    );
    assert_eq!(r.to_text(), saved);
}

#[test]
fn a_rational_composite_of_algebraic_rotations_saves_without_rounding() {
    let quarter = Isometry::around_axis(
        &RPoint3::origin(),
        &RVector3::new(0.into(), 0.into(), 1.into()),
        &Rational::from_i64(45),
        &mut Decider::default(),
    )
    .unwrap();
    let mut copy = instance(1, "base");
    copy.placement = quarter.then(&quarter);
    let expanded = expand_feature(
        &source(),
        &[id("profile"), id("cut")],
        &id("cut"),
        &[copy],
        &mut Decider::default(),
    )
    .unwrap();
    let built = evaluate(&expanded, &mut Decider::default()).unwrap();
    let mass = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert!(mass.volume.contains(&Rational::from_i64(144)));
}

#[test]
fn reflected_first_hit_features_preserve_direction_and_limit_bindings() {
    for operation in ["extrude_to", "extrude_next"] {
        let mut value: serde_json::Value = serde_json::from_str(&source().to_text()).unwrap();
        value["format"] = json!("catallax-recipe/4");
        value["steps"][2] = json!({"id":"cut","op":operation,"target":"base","profile":"profile"});
        if operation == "extrude_to" {
            value["steps"][2]["limit"] = json!({"face":{"step":"base","role":"end"}});
        }
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let copy = FeatureInstance {
            placement: Isometry::reflection(
                &RPoint3::origin(),
                &RVector3::new(0.into(), 0.into(), 1.into()),
            )
            .unwrap(),
            steps: BTreeMap::from([
                (id("base"), id("base1")),
                (id("profile"), id("profile1")),
                (id("cut"), id("cut1")),
            ]),
            inputs: BTreeMap::new(),
            origins: BTreeMap::new(),
        };
        let expanded = expand_feature(
            &recipe,
            &[id("base"), id("profile"), id("cut")],
            &id("cut"),
            &[copy],
            &mut Decider::default(),
        )
        .unwrap();
        let saved = Recipe::parse(&expanded.to_text()).unwrap();
        let built = evaluate(&saved, &mut Decider::default()).unwrap();
        let mass = catallax_measure::mass_properties(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap()
        .value;
        assert!(mass.volume.contains(&Rational::from_i64(3)));
        assert!(mass.centroid.unwrap()[2].contains(&Rational::parse("-3/2").unwrap()));
    }
}
