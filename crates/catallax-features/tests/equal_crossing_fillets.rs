//! A collapsed equal-cylinder pole is an unfinished corner, not a fitted face.

use catallax_boolean::Containment;
use catallax_features::{
    BuildError, evaluate,
    fillet::{UNSUPPORTED_CORNER, radius_limit},
};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

fn recipe(kind: &str, radius: &str, second_radius: &str, gap: &str) -> Value {
    json!({"format":"catallax-recipe/3","name":"equal-radius crossing", "steps":[
        {"id":"a","op":"cylinder","radius":"8","height":"32",
         "pose":{"rotation":["1","0","1","0"],"translation":["-16","0","0"]}},
        {"id":"b","op":"cylinder","radius":second_radius,"height":"32",
         "pose":{"rotation":["1","0","0","0"],"translation":["0",gap,"-16"]}},
        {"id":"joined","op":"boolean","kind":kind,"target":"a","tool":"b"},
        {"id":"round","op":"fillet","body":"joined","radius":radius,"edges":[
            {"between":[{"step":"a","role":{"side":{"segment":"wall"}}},
                        {"step":"b","role":{"side":{"segment":"wall"}}}]}]}
    ],"result":"round"})
}
fn parsed(value: &Value) -> Recipe {
    Recipe::parse(&value.to_string()).unwrap()
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 100_000_000,
        ..Default::default()
    })
}
fn refusal(value: &Value) -> catallax_numbers::Refusal {
    match evaluate(&parsed(value), &mut decider()) {
        Err(BuildError::Geometry { step, refusal }) => {
            assert_eq!(step.as_str(), "round");
            refusal
        }
        other => panic!("expected a fillet refusal: {other:?}"),
    }
}

#[test]
fn equal_crossing_poles_report_an_unsupported_corner() {
    let mut tilted = recipe("union", "3/2", "8", "0");
    tilted["steps"][1]["pose"] = json!({"rotation":["3","0","1","0"],
        "translation":["-48/5","0","-64/5"]});
    for value in [
        recipe("union", "1", "8", "0"),
        recipe("intersect", "1", "8", "0"),
        tilted,
    ] {
        let error = refusal(&value);
        assert_eq!(error.code, UNSUPPORTED_CORNER);
        assert_eq!(radius_limit(&error), None);
    }
}

#[test]
fn different_ball_sides_and_clipped_arcs_keep_their_contact_refusals() {
    let mut clipped = recipe("union", "1", "8", "0");
    let steps = clipped["steps"].as_array_mut().unwrap();
    let mut fillet = steps.pop().unwrap();
    steps.extend([
        json!({"id":"window","op":"block","size":["40","40","2"],
            "pose":{"rotation":["1","0","0","0"],"translation":["-20","-20","2"]}}),
        json!({"id":"trimmed","op":"boolean","kind":"intersect","target":"joined","tool":"window"}),
    ]);
    fillet["body"] = json!("trimmed");
    steps.push(fillet);
    for value in [recipe("subtract", "1", "8", "0"), clipped] {
        assert_eq!(refusal(&value).code, "fillet.unsupported-contact");
    }
}

#[test]
fn equal_radius_offset_axes_and_unequal_crossing_axes_still_build() {
    for value in [
        recipe("union", "1", "8", "1"),
        recipe("union", "1", "4", "0"),
    ] {
        let mut built = evaluate(&parsed(&value), &mut decider()).unwrap();
        assert!(built.model.check(built.result).is_empty());
        for (point, expected) in [
            (["12", "0", "0"], Containment::Inside),
            (["0", "1", "12"], Containment::Inside),
            (["12", "12", "12"], Containment::Outside),
        ] {
            let [x, y, z] = point.map(|v| Rational::parse(v).unwrap());
            let p = RPoint3::new(x, y, z);
            let outcome = catallax_boolean::classify_point(
                &mut built.model,
                built.result,
                &p,
                &mut decider(),
            );
            assert!(
                matches!(outcome, Outcome::Proven(found) if found == expected),
                "{p:?}: {outcome:?}"
            );
        }
    }
}
