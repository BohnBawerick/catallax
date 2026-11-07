//! Explicit edge orientation, legal zero ends and strict radius-law input.

use catallax_recipe::{FORMAT, FORMAT_V4, FORMAT_V5, Operation, RadiusLaw, Recipe, RecipeError};
use serde_json::{Value, json};

fn source(start: &str, end: &str) -> Value {
    json!({"format":FORMAT_V5,"name":"variable edge","result":"round","steps":[
        {"id":"b","op":"block","size":["20","20","10"],
         "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
        {"id":"round","op":"fillet_variable","body":"b","edges":[{
            "edge":{"between":[{"step":"b","role":"end"},
                {"step":"b","role":{"side":{"segment":"x_max"}}}]},
            "from":{"faces":[{"step":"b","role":"end"},
                {"step":"b","role":{"side":{"segment":"x_max"}}},
                {"step":"b","role":{"side":{"segment":"y_min"}}}]},
            "radius":{"kind":"linear","start":start,"end":end}
        }]}
    ]})
}

#[test]
fn endpoint_zero_is_legal_but_negative_and_identically_zero_laws_are_not() {
    for (start, end, valid) in [
        ("0", "2", true),
        ("2", "0", true),
        ("1/3", "2/3", true),
        ("2", "2", true),
        ("0", "0", false),
        ("-1/1000000000", "2", false),
        ("2", "-1/1000000000", false),
    ] {
        assert_eq!(
            Recipe::parse(&source(start, end).to_string()).is_ok(),
            valid
        );
    }
}

#[test]
fn orientation_survives_canonical_vertex_sorting_without_reversing_the_law() {
    let mut value = source("0", "3/2");
    value["steps"][1]["edges"][0]["from"]["faces"]
        .as_array_mut()
        .unwrap()
        .reverse();
    let saved = Recipe::parse(&value.to_string()).unwrap().to_text();
    let recipe = Recipe::parse(&saved).unwrap();
    let Operation::FilletVariable { edges, .. } = &recipe.steps[1].op else {
        panic!("expected a variable fillet");
    };
    let RadiusLaw::Linear { start, end } = &edges[0].radius else {
        panic!("expected a linear law");
    };
    assert_eq!(start.value(), &catallax_numbers::Rational::zero());
    assert_eq!(
        end.value(),
        &catallax_numbers::Rational::from_ratio(3, 2).unwrap()
    );
    assert!(edges[0].from.faces().iter().any(|face| {
        serde_json::to_value(face).unwrap()
            == json!({"step":"b","role":{"side":{"segment":"y_min"}}})
    }));
}

#[test]
fn laws_require_v5_and_an_explicit_orientation() {
    for format in [FORMAT, FORMAT_V4] {
        let mut value = source("1", "2");
        value["format"] = json!(format);
        assert!(matches!(
            Recipe::parse(&value.to_string()),
            Err(RecipeError::Parse { .. })
        ));
    }
    let mut value = source("1", "2");
    value["steps"][1]["edges"][0]
        .as_object_mut()
        .unwrap()
        .remove("from");
    assert!(Recipe::parse(&value.to_string()).is_err());
}

#[test]
fn empty_duplicate_unknown_and_inexact_laws_are_rejected() {
    let base = source("1", "2");
    let mut invalid = Vec::new();
    let mut empty = base.clone();
    empty["steps"][1]["edges"] = json!([]);
    invalid.push(empty);
    let mut duplicate = base.clone();
    let law = duplicate["steps"][1]["edges"][0].clone();
    duplicate["steps"][1]["edges"]
        .as_array_mut()
        .unwrap()
        .push(law);
    invalid.push(duplicate);
    for radius in [
        json!({"kind":"cubic","start":"1","end":"2"}),
        json!({"kind":"linear","start":1,"end":"2"}),
        json!({"kind":"linear","start":"1","end":"2","samples":[]}),
    ] {
        let mut value = base.clone();
        value["steps"][1]["edges"][0]["radius"] = radius;
        invalid.push(value);
    }
    for value in invalid {
        assert!(matches!(
            Recipe::parse(&value.to_string()),
            Err(RecipeError::Parse { .. })
        ));
    }
}

#[test]
fn endpoint_and_edge_references_cannot_bypass_origin_validation() {
    for field in ["edge", "from"] {
        let mut value = source("1", "2");
        let key = if field == "edge" { "between" } else { "faces" };
        value["steps"][1]["edges"][0][field][key][0]["step"] = json!("missing");
        assert!(matches!(
            Recipe::parse(&value.to_string()),
            Err(RecipeError::UnknownOrigin { .. })
        ));
    }
    let mut value = source("1", "2");
    value["steps"][1]["edges"][0]["from"]["faces"][0]["step"] = json!("round");
    assert!(matches!(
        Recipe::parse(&value.to_string()),
        Err(RecipeError::UnknownOrigin { .. })
    ));
}
