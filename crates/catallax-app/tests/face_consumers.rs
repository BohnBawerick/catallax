//! Face-mounted features remain usable by patterns and first-hit operations.
use catallax_app::pattern::{FeatureInstance, expand_feature};
use catallax_boolean::{Containment, classify_point};
use catallax_numbers::{Decider, Outcome, Rational};
use catallax_recipe::{Recipe, StepId};
use catallax_shapes::{Isometry, RPoint3, RVector3};
use serde_json::json;
use std::collections::BTreeMap;

fn id(name: &str) -> StepId {
    StepId::new(name).unwrap()
}

fn seed(boss: bool) -> Recipe {
    let points = if boss {
        [[-10, -10], [10, -10], [10, 10], [-10, 10]]
    } else {
        [[5, -1], [7, -1], [7, 1], [5, 1]]
    };
    Recipe::parse(
        &json!({"format":"catallax-recipe/5","name":"face copies","steps":[
        {"id":"base","op":"block","size":["20","20","2"],
            "pose":{"rotation":["1","0","0","0"],"translation":["-10","-10","0"]}},
        {"id":"profile","op":"sketch",
            "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
            "loops":[(0..4).map(|i|json!({"kind":"line","id":format!("s{i}"),
                "to":points[(i+1)%4].map(|n|n.to_string())})).collect::<Vec<_>>() ]},
        {"id":"tool","op":"extrude_on_face","body":"base","profile":"profile",
            "face":{"step":"base","role":"end"},"normal_direction":["0","0","1"],
            "x_direction":["1","0","0"],"offset":"0","distance":if boss {"2"} else {"-2"}},
        {"id":"feature","op":"boolean","kind":if boss {"union"} else {"subtract"},
            "target":"base","tool":"tool"}
    ],"result":"feature"})
        .to_string(),
    )
    .unwrap()
}

fn copy(index: u32, previous: &StepId, placement: Isometry) -> FeatureInstance {
    FeatureInstance {
        placement,
        steps: ["profile", "tool", "feature"]
            .into_iter()
            .map(|name| (id(name), id(&format!("{name}{index}"))))
            .collect(),
        inputs: BTreeMap::from([(id("base"), previous.clone())]),
        origins: BTreeMap::from([(id("base"), id("base"))]),
    }
}

fn volume(recipe: &Recipe, expected: i64) {
    let built = catallax_features::evaluate(recipe, &mut Decider::default()).unwrap();
    let mass = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert!(
        mass.volume.contains(&Rational::from_i64(expected)),
        "{mass:?}"
    );
}

fn location(recipe: &Recipe, at: [i64; 3], inside: bool) {
    let point = RPoint3::new(at[0].into(), at[1].into(), at[2].into());
    let mut dec = Decider::default();
    let mut built = catallax_features::evaluate(recipe, &mut dec).unwrap();
    assert!(
        matches!(classify_point(&mut built.model, built.result, &point, &mut dec),
        Outcome::Proven(found) if found == if inside {Containment::Inside} else {Containment::Outside})
    );
    let plan = catallax_field::Plan::compile_checked(recipe, &mut dec).unwrap();
    let found = plan
        .proven_checked(&point, &mut Default::default(), &mut dec)
        .unwrap();
    assert_eq!(
        found.location,
        if inside {
            catallax_field::Location::Inside
        } else {
            catallax_field::Location::Outside
        }
    );
}

#[test]
fn circular_copies_move_the_tool_without_moving_its_local_profile() {
    let mut recipe = seed(false);
    for index in 1..4 {
        let map = Isometry::around_axis(
            &RPoint3::origin(),
            &RVector3::new(0.into(), 0.into(), 1.into()),
            &Rational::from_i64(i64::from(index) * 90),
            &mut Decider::default(),
        )
        .unwrap();
        let instance = copy(index, &recipe.result, map);
        recipe = expand_feature(
            &recipe,
            &[id("profile"), id("tool"), id("feature")],
            &id("feature"),
            &[instance],
            &mut Decider::default(),
        )
        .unwrap();
    }
    let recipe = Recipe::parse(&recipe.to_text()).unwrap();
    volume(&recipe, 768);
    for point in [[6, 0, 1], [0, 6, 1], [-6, 0, 1], [0, -6, 1]] {
        location(&recipe, point, false);
    }
    location(&recipe, [0, 0, 1], true);
}

#[test]
fn a_boss_copy_retains_its_support_after_the_seed_consumes_that_face() {
    let recipe = seed(true);
    let instance = copy(
        1,
        &recipe.result,
        Isometry::translation(RVector3::new(0.into(), 0.into(), 2.into())),
    );
    let recipe = expand_feature(
        &recipe,
        &[id("profile"), id("tool"), id("feature")],
        &id("feature"),
        &[instance],
        &mut Decider::default(),
    )
    .unwrap();
    volume(&recipe, 2400);
    location(&recipe, [0, 0, 5], true);
}

#[test]
fn copying_the_support_or_reflecting_the_face_feature_refuses() {
    let recipe = seed(false);
    let mut instance = copy(
        1,
        &recipe.result,
        Isometry::translation(RVector3::new(1.into(), 0.into(), 0.into())),
    );
    instance.steps.insert(id("base"), id("base1"));
    instance.inputs.clear();
    instance.origins.clear();
    let error = expand_feature(
        &recipe,
        &[id("base"), id("profile"), id("tool"), id("feature")],
        &id("feature"),
        &[instance],
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code, "pattern.expansion");
    let reflected = copy(
        1,
        &recipe.result,
        Isometry::reflection(
            &RPoint3::origin(),
            &RVector3::new(1.into(), 0.into(), 0.into()),
        )
        .unwrap(),
    );
    let error = expand_feature(
        &recipe,
        &[id("profile"), id("tool"), id("feature")],
        &id("feature"),
        &[reflected],
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code, "pattern.expansion");
}

#[test]
fn up_to_next_accepts_a_target_built_from_face_mounted_bosses() {
    let mut input: serde_json::Value =
        serde_json::from_str(include_str!("fixtures/face-extrusion-next.json")).unwrap();
    let tool = Recipe::parse(&input.to_string()).unwrap();
    volume(&tool, 180);
    input["steps"]
        .as_array_mut()
        .unwrap()
        .push(json!({"id":"cut","op":"boolean",
        "kind":"subtract","target":"op_98","tool":"op_117"}));
    input["result"] = json!("cut");
    let recipe = Recipe::parse(&input.to_string()).unwrap();
    volume(&recipe, 35820);
    location(&recipe, [52, 32, 37], false);
    location(&recipe, [52, 32, 2], true);
}
