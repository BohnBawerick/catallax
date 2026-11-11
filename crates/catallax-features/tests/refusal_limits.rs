//! R12 preserves regularised geometry while qualifying contact diagnostics.

use catallax_features::{BuildError, evaluate};
use catallax_numbers::{BudgetPolicy, Decider, Rational};
use catallax_recipe::Recipe;
use serde_json::{Value, json};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    })
}

fn case(id: &str) -> Value {
    include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|row| row["id"] == id)
        .unwrap()["recipe"]
        .clone()
}

#[test]
fn the_ten_mm_cube_keeps_its_terminal_roll_over_at_and_above_ten() {
    for (id, requested) in [
        ("r12-cube-at-contact", None),
        ("r12-cube-past-contact", Some("10.001")),
    ] {
        let recipe = Recipe::parse(&case(id).to_string()).unwrap();
        let built = evaluate(&recipe, &mut decider()).unwrap();
        assert!(built.model.check(built.result).is_empty());
        let [definition] = built.fillets.as_slice() else {
            panic!("expected one exact strip");
        };
        assert_eq!(
            definition.radius.as_rational(),
            Some(Rational::from_i64(10))
        );
        assert_eq!(
            definition
                .requested_radius
                .as_ref()
                .and_then(|r| r.as_rational()),
            requested.map(|r| Rational::parse(r).unwrap())
        );
    }
}

#[test]
fn opposing_shell_walls_consuming_the_cavity_still_succeed() {
    let mut value = case("r12-cube-at-contact");
    value["steps"].as_array_mut().unwrap().truncate(1);
    value["steps"].as_array_mut().unwrap().push(json!({
        "id":"shell", "op":"shell", "body":"b", "open":[], "thickness":"6"
    }));
    value["result"] = json!("shell");
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let built = evaluate(&recipe, &mut decider()).unwrap();
    assert_eq!(
        built.result,
        built.bodies[&catallax_recipe::StepId::new("b").unwrap()]
    );
    let mass = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut decider(),
    )
    .unwrap()
    .value;
    assert_eq!(mass.volume.lo, Rational::from_i64(1000));
    assert_eq!(mass.volume.hi, Rational::from_i64(1000));
}

#[test]
fn an_uncovered_end_reports_an_unknown_maximum_without_inventing_a_bound() {
    let recipe = Recipe::parse(&case("partial-oblique").to_string()).unwrap();
    let BuildError::Geometry { refusal, .. } = evaluate(&recipe, &mut decider()).unwrap_err()
    else {
        panic!("expected geometric coverage refusal");
    };
    assert_eq!(refusal.code, "fillet.unsupported-end");
    assert!(refusal.detail.contains("maximum unknown"), "{refusal:?}");
    assert_eq!(catallax_features::fillet::radius_limit(&refusal), None);
}

#[test]
fn an_uncovered_shell_does_not_treat_a_gap_as_a_thickness_failure() {
    let mut value = case("box-single");
    value["steps"].as_array_mut().unwrap().push(json!({
        "id":"shell", "op":"shell", "body":"round", "open":[], "thickness":"1/4"
    }));
    value["result"] = json!("shell");
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let BuildError::Geometry { refusal, .. } = evaluate(&recipe, &mut decider()).unwrap_err()
    else {
        panic!("expected geometric coverage refusal");
    };
    assert!(refusal.code.starts_with("distance."), "{refusal:?}");
    assert!(refusal.detail.contains("maximum unknown"), "{refusal:?}");
    assert!(!refusal.detail.contains("thickness bound"), "{refusal:?}");
}
