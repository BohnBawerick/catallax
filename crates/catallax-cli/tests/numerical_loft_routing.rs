//! A top-stock route needs overlap, not only a wide projection.

use catallax_boolean::numerical::loft::stock_overlaps_loft_top;
use catallax_numbers::{Decider, Rational, current_operation_work, with_operation_budget};
use catallax_recipe::{Recipe, StepId};
use serde_json::{Value, json};

#[test]
fn top_stock_routing_keeps_disjoint_and_boundary_cases_on_the_ordinary_path() {
    for (name, size, at, expected) in [
        ("original", ["10", "10", "1"], ["-3", "5/2", "1/2"], true),
        ("beyond-cap", ["10", "10", "1"], ["-3", "5/2", "3"], false),
        ("above-roof", ["10", "10", "1"], ["-3", "10", "1/2"], false),
        ("below-floor", ["10", "1", "1"], ["-3", "-2", "1/2"], false),
        ("on-cap", ["10", "10", "1"], ["-3", "5/2", "0"], false),
        ("side-stock", ["2", "10", "1"], ["1", "1", "1/2"], false),
    ] {
        let (routing_work, work) = with_operation_budget(20_000_000, || {
            let mut input: Value = serde_json::from_str(include_str!(
                "../../catallax-features/tests/fixtures/m10-smooth-loft.json"
            ))
            .unwrap();
            input["steps"].as_array_mut().unwrap().push(json!({
                "id":"box","op":"block","size":size,
                "pose":{"rotation":["1","0","0","0"],"translation":at}
            }));
            let recipe = Recipe::parse(&input.to_string()).unwrap();
            let mut dec = Decider::default();
            let built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
            let loft = built.bodies[&StepId::new("loft").unwrap()];
            let stock = built.bodies[&StepId::new("box").unwrap()];
            let start = current_operation_work().used;
            let selected = stock_overlaps_loft_top(
                &built.model,
                loft,
                stock,
                &Rational::parse("1/10000000").unwrap(),
                &mut dec,
            )
            .unwrap();
            assert_eq!(selected, expected, "{name}");
            current_operation_work().used - start
        });
        assert!(!work.exhausted, "{name}");
        println!("{name}: {routing_work} routing work, {} total", work.used);
    }
}
