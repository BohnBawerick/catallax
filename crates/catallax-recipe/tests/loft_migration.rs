//! Rational-chart correspondence is not a relabelling of legacy arc-length lofts.

use catallax_recipe::{FORMAT, FORMAT_V4, FORMAT_V5, LoftParameterization, Recipe, RecipeError};
use serde_json::{Value, json};

fn pose(z: &str) -> Value {
    json!({"rotation":["1","0","0","0"],"translation":["0","0",z]})
}

fn loft() -> Value {
    json!({"format":FORMAT_V4,"name":"loft schema","steps":[
        {"id":"lower","op":"sketch","plane":pose("0"),"loops":[[
            {"kind":"line","id":"a","to":["2","0"]},
            {"kind":"line","id":"b","to":["2","2"]},
            {"kind":"line","id":"c","to":["0","2"]},
            {"kind":"line","id":"d","to":["0","0"]}
        ],[{"kind":"circle","id":"hole","centre":["1","1"],"radius":"1/4"}]]},
        {"id":"upper","op":"sketch","plane":pose("3"),"loops":[[
            {"kind":"line","id":"x","to":["3","0"]},
            {"kind":"line","id":"y","to":["0","3"]},
            {"kind":"line","id":"z","to":["0","0"]}
        ],[{"kind":"circle","id":"bore","centre":["1","1"],"radius":"1/4"}]]},
        {"id":"solid","op":"loft","shape":"ruled","sections":[
            {"profile":"lower","start":["hole","b"]},
            {"profile":"upper","start":["bore","z"]}
        ]}
    ],"result":"solid"})
}

fn parse(value: &Value) -> Result<Recipe, RecipeError> {
    Recipe::parse(&value.to_string())
}

#[test]
fn rational_lofts_accept_unequal_segment_counts_and_keep_authored_loop_order() {
    for format in [FORMAT_V4, FORMAT_V5] {
        for shape in ["ruled", "smooth"] {
            let mut value = loft();
            value["format"] = json!(format);
            value["steps"][2]["shape"] = json!(shape);
            let recipe = parse(&value).unwrap();
            assert_eq!(
                recipe.loft_parameterization().unwrap(),
                LoftParameterization::RationalCharts
            );
            let restored = Recipe::parse(&recipe.to_text()).unwrap();
            assert_eq!(restored.steps[2], recipe.steps[2]);
            assert_eq!(
                serde_json::to_value(&restored.steps[2]).unwrap()["sections"][0]["start"],
                json!(["hole", "b"])
            );
        }
    }
}

#[test]
fn loop_counts_must_agree_between_sections() {
    let mut value = loft();
    value["steps"][1]["loops"].as_array_mut().unwrap().pop();
    value["steps"][2]["sections"][1]["start"] = json!(["x"]);
    let error = parse(&value).unwrap_err();
    assert!(matches!(
        error,
        RecipeError::LoftLoops {
            expected: 2,
            found: 1,
            ..
        }
    ));
    assert_eq!(error.code(), "recipe.loft-loops");
}

#[test]
fn a_missing_loop_start_is_not_a_loop_count_error() {
    let mut value = loft();
    value["steps"][2]["sections"][1]["start"] = json!(["x"]);
    assert_eq!(parse(&value).unwrap_err().code(), "recipe.loft-starts");
}

#[test]
fn repeated_start_in_one_loop_does_not_cover_another_loop() {
    let mut value = loft();
    value["steps"][2]["sections"][1]["start"] = json!(["x", "z"]);
    assert_eq!(parse(&value).unwrap_err().code(), "recipe.loft-starts");
}

#[test]
fn a_start_from_another_section_is_not_a_correspondence_guess() {
    let mut value = loft();
    value["steps"][2]["sections"][1]["start"] = json!(["hole", "z"]);
    assert_eq!(parse(&value).unwrap_err().code(), "recipe.unknown-element");
}

#[test]
fn a_loft_requires_at_least_two_sections() {
    let mut value = loft();
    value["steps"][2]["sections"].as_array_mut().unwrap().pop();
    assert_eq!(parse(&value).unwrap_err().code(), "recipe.loft-sections");
}

#[test]
fn legacy_loft_reads_and_writes_without_changing_its_law() {
    let mut value = loft();
    value["format"] = json!(FORMAT);
    let recipe = parse(&value).unwrap();
    assert_eq!(
        recipe.loft_parameterization().unwrap(),
        LoftParameterization::ArcLength
    );
    let text = recipe.to_text();
    let restored = Recipe::parse(&text).unwrap();
    assert_eq!(restored.format, FORMAT);
    assert_eq!(
        restored.loft_parameterization().unwrap(),
        LoftParameterization::ArcLength
    );
    assert_eq!(restored.to_text(), text);
}

#[test]
fn legacy_loft_cannot_be_automatically_upgraded_even_with_only_lines() {
    let mut value = loft();
    value["format"] = json!(FORMAT);
    for section in [0, 1] {
        value["steps"][section]["loops"]
            .as_array_mut()
            .unwrap()
            .pop();
        value["steps"][2]["sections"][section]["start"]
            .as_array_mut()
            .unwrap()
            .remove(0);
    }
    let recipe = parse(&value).unwrap();
    let original = recipe.to_text();
    for target in [FORMAT_V4, FORMAT_V5] {
        let error = recipe.upgrade_to(target).unwrap_err();
        assert!(
            matches!(&error, RecipeError::MigrationReview { step: Some(step), .. } if step.as_str() == "solid")
        );
        assert_eq!(error.code(), "recipe.migration-review");
    }
    assert_eq!(recipe.to_text(), original);
}

#[test]
fn reviewed_rational_loft_upgrade_keeps_stations_starts_and_shape() {
    let recipe = parse(&loft()).unwrap();
    let upgraded = recipe.upgrade_to(FORMAT_V5).unwrap();
    assert_eq!(upgraded.steps, recipe.steps);
    assert_eq!(upgraded.result, recipe.result);
    assert_eq!(
        upgraded.loft_parameterization().unwrap(),
        LoftParameterization::RationalCharts
    );
    assert_eq!(Recipe::parse(&upgraded.to_text()).unwrap(), upgraded);
}

fn block_with(op: Value) -> Recipe {
    Recipe::parse(
        &json!({"format":FORMAT,"name":"migration review","steps":[
        {"id":"base","op":"block","size":["4","4","4"],"pose":pose("0")}, op
    ],"result":"edit"})
        .to_string(),
    )
    .unwrap()
}

#[test]
fn version_dependent_feature_defaults_also_require_review() {
    for op in [
        json!({"id":"edit","op":"fillet","body":"base","edges":"all","radius":"1"}),
        json!({"id":"edit","op":"chamfer","body":"base","edges":"all","distance":"1"}),
        json!({"id":"edit","op":"thicken","body":"base","faces":[{"face":{"step":"base","role":"end"}}],"thickness":"1","side":"front"}),
        json!({"id":"edit","op":"draft","body":"base","faces":[{"face":{"step":"base","role":"end"}}],"pull":["0","0","1"],"neutral":{"point":["0","0","0"],"normal":["0","0","1"]},"angle":"5"}),
    ] {
        let recipe = block_with(op);
        assert!(
            matches!(recipe.upgrade_to(FORMAT_V4), Err(RecipeError::MigrationReview { step: Some(step), .. }) if step.as_str() == "edit")
        );
    }
}

#[test]
fn legacy_sweep_transport_requires_review_in_both_orientation_modes() {
    let mut recipe = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/straight-sweep.json"
    ))
    .unwrap();
    recipe.format = FORMAT.into();
    for orientation in [
        catallax_recipe::SweepOrientation::Fixed,
        catallax_recipe::SweepOrientation::RotationMinimising,
    ] {
        let catallax_recipe::Operation::Sweep {
            orientation: mode, ..
        } = &mut recipe.steps[2].op
        else {
            panic!()
        };
        *mode = orientation;
        assert_eq!(
            recipe.upgrade_to(FORMAT_V4).unwrap_err().code(),
            "recipe.migration-review"
        );
    }
}

#[test]
fn unchanged_operation_graph_can_upgrade_without_rewriting_names() {
    let recipe = block_with(
        json!({"id":"edit","op":"mirror","body":"base","keep_original":true,"plane":{"point":["6","0","0"],"normal":["1","0","0"]}}),
    );
    let original = recipe.to_text();
    let upgraded = recipe.upgrade_to(FORMAT_V4).unwrap();
    assert_eq!(upgraded.format, FORMAT_V4);
    assert_eq!(upgraded.steps, recipe.steps);
    assert_eq!(recipe.to_text(), original);
    assert_eq!(upgraded.upgrade_to(FORMAT_V4).unwrap(), upgraded);
    assert_eq!(Recipe::parse(&upgraded.to_text()).unwrap(), upgraded);
}

#[test]
fn migration_rejects_downgrades_unknown_formats_and_invalid_input() {
    let recipe = parse(&loft()).unwrap();
    assert_eq!(
        recipe.upgrade_to(FORMAT).unwrap_err().code(),
        "recipe.migration-review"
    );
    assert_eq!(
        recipe.upgrade_to("catallax-recipe/6").unwrap_err().code(),
        "recipe.format"
    );
    let mut invalid = recipe;
    invalid.format = "catallax-recipe/6".into();
    assert_eq!(
        invalid.loft_parameterization().unwrap_err().code(),
        "recipe.format"
    );
    assert_eq!(
        invalid.upgrade_to(FORMAT_V5).unwrap_err().code(),
        "recipe.format"
    );
}
