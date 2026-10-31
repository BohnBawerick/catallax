//! Independent acceptance for an early refusal of the buried capsule trim.
use catallax_app::{ConstructionOptions, Error, construction_preview};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;

#[test]
fn slotted_bracket_refuses_with_a_typed_reason_before_half_the_reserve() {
    let source = include_str!("fixtures/construction-slot.json");
    for depth in ["11", "14", "20"] {
        let mut value: serde_json::Value = serde_json::from_str(source).unwrap();
        value["steps"][1]["distance"] = depth.into();
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let options = ConstructionOptions {
            tolerance: "0.16158836963098594".parse().unwrap(),
            ..ConstructionOptions::default()
        };
        let (answer, work) = with_operation_budget(2_000_000, || {
            construction_preview(&recipe, &options, &mut Decider::default())
        });
        assert!(
            matches!(&answer, Err(Error::ConstructionPreview(r)) if r.code == "field.display-trim-obstruction"),
            "depth {depth}: expected the structural typed refusal, got {answer:?}; {work:?}"
        );
        assert!(!work.exhausted, "depth {depth}: {work:?}");
        assert!(work.used <= 1_000_000, "depth {depth}: {work:?}");
    }
}
