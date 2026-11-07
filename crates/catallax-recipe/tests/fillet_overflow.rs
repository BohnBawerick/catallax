use catallax_recipe::{FilletOverflow, Operation, Recipe, RecipeError};
use serde_json::{Value, json};

fn face(role: Value) -> Value {
    json!({"step": "b", "role": role})
}

fn recipe() -> Value {
    let end = face(json!("end"));
    let side = face(json!({"side": {"segment": "x_max"}}));
    let edge = json!({"between": [end, side]});
    json!({
        "format": "catallax-recipe/4", "name": "setback", "result": "round",
        "steps": [
            {"id": "b", "op": "block", "size": ["20", "20", "10"],
             "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]}},
            {"id": "round", "op": "fillet", "body": "b", "radius": "1",
             "edges": [edge], "overflow": {"keep_neighbours": [
                 {"edge": edge, "at": {"faces": [end, side, face(json!({"side": {"segment": "y_min"}}))]},
                  "setback": "2"}
             ]}}
        ]
    })
}

fn parse(value: &Value) -> Result<Recipe, RecipeError> {
    Recipe::parse(&value.to_string())
}

#[test]
fn station_policy_retains_exact_distance_and_canonical_endpoint() {
    let mut value = recipe();
    value["steps"][1]["overflow"]["keep_neighbours"][0]["setback"] = json!("3/2");
    let parsed = parse(&value).unwrap();
    let Operation::Fillet {
        overflow: Some(FilletOverflow::KeepNeighbours(limits)),
        ..
    } = &parsed.steps[1].op
    else {
        panic!("missing authored station")
    };
    assert_eq!(
        limits[0].setback,
        catallax_recipe::Number::parse("3/2").unwrap()
    );
    assert_eq!(
        limits[0].at,
        serde_json::from_value(value["steps"][1]["overflow"]["keep_neighbours"][0]["at"].clone())
            .unwrap()
    );
    assert_eq!(Recipe::parse(&parsed.to_text()).unwrap(), parsed);
    for version in ["catallax-recipe/4", "catallax-recipe/5"] {
        value["format"] = json!(version);
        value["steps"][1]["overflow"] = json!("roll_over");
        assert!(parse(&value).is_ok());
    }
    value["format"] = json!("catallax-recipe/3");
    assert!(matches!(parse(&value), Err(RecipeError::Parse { .. })));
    value["steps"][1]
        .as_object_mut()
        .unwrap()
        .remove("overflow");
    assert!(parse(&value).is_ok());
}

#[test]
fn station_references_are_checked_before_geometry() {
    let mut value = recipe();
    value["steps"][1]["overflow"]["keep_neighbours"][0]["at"]["faces"][0]["step"] = json!("round");
    assert!(matches!(
        parse(&value),
        Err(RecipeError::UnknownOrigin { .. })
    ));
    let mut value = recipe();
    value["steps"][1]["overflow"]["keep_neighbours"][0]["edge"]["between"][0]["role"] =
        json!("beams");
    assert!(matches!(parse(&value), Err(RecipeError::WrongRole { .. })));
}

#[test]
fn station_syntax_never_accepts_nonpositive_or_inexact_distances() {
    for invalid in [json!("0"), json!("-1"), json!(2), Value::Null] {
        let mut value = recipe();
        value["steps"][1]["overflow"]["keep_neighbours"][0]["setback"] = invalid;
        assert!(matches!(parse(&value), Err(RecipeError::Parse { .. })));
    }
    let mut value = recipe();
    value["steps"][1]["overflow"]["keep_neighbours"][0]["at"] = json!({"faces": []});
    assert!(matches!(parse(&value), Err(RecipeError::Parse { .. })));
    let mut value = recipe();
    value["steps"][1]["overflow"]["keep_neighbours"][0]["distance"] = json!("2");
    assert!(matches!(parse(&value), Err(RecipeError::Parse { .. })));
}
