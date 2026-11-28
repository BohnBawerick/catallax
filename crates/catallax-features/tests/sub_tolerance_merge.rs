//! The captured SolidDraft bracket's base and stand have parallel right walls
//! 3.5164312091753774e-36 mm apart. The union merges that gap once, in the
//! model, instead of leaving a sliver face for every export to repair.
use catallax_boolean::SUB_TOLERANCE_MERGE;
use catallax_features::evaluate;
use catallax_numbers::{BoundKind, Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;

const CAPTURE: &str = include_str!("../../catallax-cli/tests/fixtures/soliddraft-slot-cut.json");

fn prefix(result: &str) -> Recipe {
    let question: serde_json::Value = serde_json::from_str(CAPTURE).unwrap();
    let mut recipe = Recipe::parse(&question["problem"].to_string()).unwrap();
    let end = recipe
        .steps
        .iter()
        .position(|s| s.id.as_str() == result)
        .unwrap();
    recipe.steps.truncate(end + 1);
    recipe.result = recipe.steps[end].id.clone();
    recipe
}

#[test]
fn the_captured_bracket_union_merges_its_sub_tolerance_gap() {
    let recipe = prefix("op_61");
    let (built, work) =
        with_operation_budget(100_000_000, || evaluate(&recipe, &mut Decider::default()));
    assert!(!work.exhausted, "{work:?}");
    let built = built.unwrap();
    let body = built.model.body(built.result);
    // One sliver fewer; the stand's tilted walls stay for the exporter.
    assert_eq!(built.model.counts(built.result).faces, 10);
    assert!(built.model.check(built.result).is_empty());
    let marks: Vec<_> = body.marks.iter().collect();
    let [mark] = marks.as_slice() else {
        panic!("{marks:?}");
    };
    assert_eq!(mark.code, SUB_TOLERANCE_MERGE);
    assert_eq!(mark.kind, BoundKind::Distance);
    let gap = Rational::parse("3.5164312091753774e-36").unwrap();
    // The corners slide along the stand's nearly level edges, so the bound is
    // the gap and a sliver more.
    assert!(mark.bound >= gap, "{}", mark.bound);
    assert!(mark.bound < &gap + &gap, "{}", mark.bound);
}
