//! Versioned rational paths, exact joins and complete structural invariants.

use catallax_recipe::{FORMAT, FORMAT_V4, FORMAT_V5, Operation, PathSegment, Recipe};
use serde_json::{Value, json};

fn spline() -> Value {
    json!({"kind":"spline","id":"rail","degree":2,
        "control_points":[["0","0","0"],["1","0","1"],["2","0","2"],["3","0","3"],["4","0","4"]],
        "weights":["1","2","3/2","2","1"],
        "knots":["0","0","0","1/2","1/2","1","1","1"]})
}

fn source(segment: Value) -> Value {
    json!({"format":FORMAT_V4,"name":"rational path","steps":[
        {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
         "loops":[[{"kind":"circle","id":"rim","centre":["0","0"],"radius":"1"}]]},
        {"id":"path","op":"path","start":["0","0","0"],"segments":[segment]},
        {"id":"swept","op":"sweep","profile":"profile","path":"path","orientation":"fixed"}
    ],"result":"swept"})
}

fn parse(value: &Value) -> Result<Recipe, catallax_recipe::RecipeError> {
    Recipe::parse(&value.to_string())
}

#[test]
fn weighted_repeated_knot_paths_keep_their_exact_definition_in_both_modes() {
    for format in [FORMAT_V4, FORMAT_V5] {
        for orientation in ["fixed", "rotation_minimising"] {
            let mut value = source(spline());
            value["format"] = json!(format);
            value["steps"][2]["orientation"] = json!(orientation);
            let recipe = parse(&value).unwrap();
            let restored = Recipe::parse(&recipe.to_text()).unwrap();
            assert_eq!(restored, recipe);
            let Operation::Path { segments, .. } = &restored.steps[1].op else {
                panic!()
            };
            let PathSegment::Spline { weights, knots, .. } = &segments[0] else {
                panic!()
            };
            assert_eq!(weights[2].canonical(), "1.5");
            assert_eq!(knots[3], knots[4]);
            assert_eq!(segments[0].id().as_str(), "rail");
        }
    }
}

#[test]
fn version_three_rejects_spline_paths() {
    let mut value = source(spline());
    value["format"] = json!(FORMAT);
    assert_eq!(parse(&value).unwrap_err().code(), "recipe.parse");
}

macro_rules! invalid_spline {
    ($name:ident, $field:literal, $value:expr, $code:literal) => {
        #[test]
        fn $name() {
            let mut segment = spline();
            segment[$field] = json!($value);
            assert_eq!(parse(&source(segment)).unwrap_err().code(), $code);
        }
    };
}

invalid_spline!(zero_degree, "degree", 0, "recipe.path-spline");
invalid_spline!(degree_at_control_count, "degree", 5, "recipe.path-spline");
invalid_spline!(degree_overflow, "degree", usize::MAX, "recipe.path-spline");
invalid_spline!(negative_degree, "degree", -1, "recipe.parse");
invalid_spline!(
    empty_controls,
    "control_points",
    json!([]),
    "recipe.path-spline"
);
invalid_spline!(wrong_weight_count, "weights", ["1"], "recipe.path-spline");
invalid_spline!(
    zero_weight,
    "weights",
    ["1", "1", "0", "1", "1"],
    "recipe.path-spline"
);
invalid_spline!(
    negative_weight,
    "weights",
    ["1", "1", "-1", "1", "1"],
    "recipe.path-spline"
);
invalid_spline!(
    numeric_weight,
    "weights",
    json!(["1", "1", 1, "1", "1"]),
    "recipe.parse"
);
invalid_spline!(
    wrong_knot_count,
    "knots",
    ["0", "0", "0", "1", "1", "1"],
    "recipe.path-spline"
);
invalid_spline!(
    decreasing_knots,
    "knots",
    ["0", "0", "0", "3/4", "1/2", "1", "1", "1"],
    "recipe.path-spline"
);
invalid_spline!(
    unclamped_start,
    "knots",
    ["0", "0", "1/4", "1/2", "1/2", "1", "1", "1"],
    "recipe.path-spline"
);
invalid_spline!(
    unclamped_end,
    "knots",
    ["0", "0", "0", "1/2", "1/2", "3/4", "1", "1"],
    "recipe.path-spline"
);
invalid_spline!(
    extra_start_knot,
    "knots",
    ["0", "0", "0", "0", "1/2", "1", "1", "1"],
    "recipe.path-spline"
);
invalid_spline!(
    extra_end_knot,
    "knots",
    ["0", "0", "0", "1/2", "1", "1", "1", "1"],
    "recipe.path-spline"
);
invalid_spline!(
    outside_domain,
    "knots",
    ["-1", "-1", "-1", "1/2", "1/2", "1", "1", "1"],
    "recipe.path-spline"
);
invalid_spline!(unknown_spline_field, "to", ["4", "0", "4"], "recipe.parse");

#[test]
fn interior_knot_multiplicity_cannot_exceed_the_degree() {
    let mut segment = spline();
    segment["control_points"]
        .as_array_mut()
        .unwrap()
        .push(json!(["5", "0", "5"]));
    segment["weights"].as_array_mut().unwrap().push(json!("1"));
    segment["knots"] = json!(["0", "0", "0", "1/2", "1/2", "1/2", "1", "1", "1"]);
    assert_eq!(
        parse(&source(segment)).unwrap_err().code(),
        "recipe.path-spline"
    );
}

#[test]
fn joins_follow_the_exact_last_control_across_mixed_segment_kinds() {
    let mut value = source(spline());
    value["steps"][1]["segments"] = json!([
        {"kind":"line","id":"lead","to":["0","0","0"]},
        spline(),
        {"kind":"arc","id":"turn","centre":["4","1","4"],"to":["5","1","4"],"axis":["0","0","1"]},
        {"kind":"spline","id":"tail","degree":1,"control_points":[["5","1","4"],["6","1","4"]],"weights":["1","1"],"knots":["0","0","1","1"]}
    ]);
    value["steps"][1]["start"] = json!(["-1", "0", "0"]);
    let recipe = parse(&value).unwrap();
    assert_eq!(Recipe::parse(&recipe.to_text()).unwrap(), recipe);
    value["steps"][1]["segments"][3]["control_points"][0][0] = json!("5.00000000000000000001");
    assert_eq!(parse(&value).unwrap_err().code(), "recipe.path-spline");
}

#[test]
fn first_control_must_equal_path_start_without_a_join_tolerance() {
    let mut value = source(spline());
    value["steps"][1]["start"][0] = json!("1e-30");
    assert_eq!(parse(&value).unwrap_err().code(), "recipe.path-spline");
}

#[test]
fn repeated_controls_and_stationary_points_are_not_schema_errors() {
    let mut segment = spline();
    segment["control_points"][1] = segment["control_points"][0].clone();
    let recipe = parse(&source(segment)).unwrap();
    assert_eq!(Recipe::parse(&recipe.to_text()).unwrap(), recipe);
}

#[test]
fn spline_corner_axes_are_retained_and_ids_stay_unique() {
    let mut segment = spline();
    segment["corner_axis"] = json!(["0", "1", "0"]);
    let mut value = source(segment);
    let recipe = parse(&value).unwrap();
    let Operation::Path { segments, .. } = &recipe.steps[1].op else {
        panic!()
    };
    assert_eq!(segments[0].corner_axis().unwrap().0[1].canonical(), "1");
    value["steps"][1]["segments"]
        .as_array_mut()
        .unwrap()
        .push(json!({"kind":"line","id":"rail","to":["5","0","4"]}));
    assert_eq!(
        parse(&value).unwrap_err().code(),
        "recipe.duplicate-element"
    );
}
