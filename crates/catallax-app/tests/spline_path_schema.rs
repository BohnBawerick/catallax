//! Saved path instances preserve rational controls and oriented corner axes.

use std::collections::BTreeMap;

use catallax_app::pattern::{FeatureInstance, expand_feature};
use catallax_numbers::Decider;
use catallax_recipe::{Operation, PathSegment, Recipe, StepId};
use catallax_shapes::{Isometry, RPoint3, RVector3};
use serde_json::json;

fn id(value: &str) -> StepId {
    StepId::new(value).unwrap()
}

fn source() -> Recipe {
    Recipe::parse(&json!({"format":"catallax-recipe/4","name":"spline definition","steps":[
        {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
         "loops":[[{"kind":"circle","id":"rim","centre":["0","0"],"radius":"1"}]]},
        {"id":"path","op":"path","start":["-1","0","0"],"segments":[
            {"kind":"line","id":"lead","to":["0","0","0"]},
            {"kind":"spline","id":"rail","degree":2,"control_points":[["0","0","0"],["-1","0","0"],["-2","0","2"]],
             "weights":["1","2","1"],"knots":["0","0","0","1","1","1"],"corner_axis":["0","0","1"]}
        ]},
        {"id":"swept","op":"sweep","profile":"profile","path":"path","orientation":"fixed"}
    ],"result":"swept"}).to_string()).unwrap()
}

#[test]
fn reflected_translated_path_controls_and_half_turn_axis_save_exactly() {
    let recipe = source();
    let original = recipe.to_text();
    let placement = Isometry::reflection(
        &RPoint3::origin(),
        &RVector3::new(1.into(), 0.into(), 0.into()),
    )
    .unwrap()
    .then(&Isometry::translation(RVector3::new(
        3.into(),
        5.into(),
        7.into(),
    )));
    let graph = [id("profile"), id("path"), id("swept")];
    let instance = FeatureInstance {
        placement,
        steps: BTreeMap::from([
            (id("profile"), id("profile-copy")),
            (id("path"), id("path-copy")),
            (id("swept"), id("swept-copy")),
        ]),
        inputs: BTreeMap::new(),
        origins: BTreeMap::new(),
    };
    let expanded = expand_feature(
        &recipe,
        &graph,
        &id("swept"),
        &[instance],
        &mut Decider::default(),
    )
    .unwrap();
    let saved = Recipe::parse(&expanded.to_text()).unwrap();
    let Operation::Path { start, segments } = &saved.step(&id("path-copy")).unwrap().op else {
        panic!()
    };
    assert_eq!(start.0, [4, 5, 7].map(Into::into));
    let PathSegment::Spline {
        control_points,
        weights,
        knots,
        corner_axis,
        ..
    } = &segments[1]
    else {
        panic!()
    };
    assert_eq!(
        control_points
            .iter()
            .map(|p| p.0.clone())
            .collect::<Vec<_>>(),
        vec![
            [3, 5, 7].map(Into::into),
            [4, 5, 7].map(Into::into),
            [5, 5, 9].map(Into::into)
        ]
    );
    assert_eq!(weights, &[1.into(), 2.into(), 1.into()]);
    assert_eq!(
        knots,
        &[0.into(), 0.into(), 0.into(), 1.into(), 1.into(), 1.into()]
    );
    assert_eq!(corner_axis.as_ref().unwrap().0, [0, 0, -1].map(Into::into));
    let Operation::Sweep { profile, path, .. } = &saved.step(&id("swept-copy")).unwrap().op else {
        panic!()
    };
    assert_eq!(profile, &id("profile-copy"));
    assert_eq!(path, &id("path-copy"));
    assert_eq!(recipe.to_text(), original);
}

#[test]
fn spline_paths_are_not_replaced_with_sampled_line_sweeps() {
    let recipe = source();
    let error =
        catallax_app::build(&recipe, &Default::default(), &mut Decider::default()).unwrap_err();
    assert_eq!(error.code(), "features.unsupported-path");
    assert!(
        matches!(error, catallax_app::Error::Build(catallax_features::BuildError::Geometry { step, .. }) if step == id("swept"))
    );
}
