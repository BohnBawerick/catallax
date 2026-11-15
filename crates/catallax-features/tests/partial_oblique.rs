//! Open oblique rims share their end section with the vertical wall.

use catallax_features::evaluate;
use catallax_numbers::{BudgetPolicy, Decider};
use catallax_recipe::Recipe;

#[test]
fn the_public_partial_oblique_rim_builds_and_exports() {
    let recipe = Recipe::parse(include_str!(
        "../../../evidence/phase1/P03/inputs/partial-oblique.json"
    ))
    .unwrap();
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    });
    let built = evaluate(&recipe, &mut dec).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let mut bytes = Vec::new();
    catallax_step::write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut dec,
        &mut bytes,
    )
    .unwrap();
    assert!(!bytes.is_empty());
}
