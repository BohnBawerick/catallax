//! Numerical continuation of the actual ruled-loft crossing and its physical hole.

use catallax_numbers::{Decider, NumericalCheck, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_topology::numerical::{self, NumericalOptions};
use serde_json::{Value, json};

#[test]
fn the_loft_crossing_keeps_its_curved_edges_and_two_loop_face_numerically() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut input: Value =
            serde_json::from_str(include_str!("fixtures/loft-spline-prism.json")).unwrap();
        let steps = input["steps"].as_array_mut().unwrap();
        steps.push(json!({"id":"box","op":"block","size":["10","10","1"],
            "pose":{"rotation":["1","0","0","0"],"translation":["-3","1999999/1000000","1"]}}));
        steps.push(
            json!({"id":"cut","op":"boolean","kind":"subtract","target":"loft","tool":"box"}),
        );
        input["result"] = json!("cut");
        let recipe = Recipe::parse(&input.to_string()).unwrap();
        let mut dec = Decider::default();
        let built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
        let definition = numerical::definition_from_model(
            &built.model,
            built.result,
            "loft crossing".into(),
            &mut dec,
        )
        .unwrap();
        assert!(definition.faces.iter().any(|face| face.loops.len() == 2));
        assert!(
            definition
                .edges
                .iter()
                .any(|edge| matches!(edge.curve, catallax_shapes::Curve::BSpline(_)))
        );
        let count = built.model.counts(built.result);
        let result = numerical::build(
            definition,
            NumericalOptions {
                tolerance_mm: Rational::parse("1/10000000").unwrap(),
            },
            &mut dec,
        )
        .unwrap();
        assert_eq!(result.counts(), count);
        assert!(result.evidence().checks().iter().any(|check| matches!(check, NumericalCheck::Passed {name,..} if name=="trim.sampled-simple-oriented-parameter-loop")));
    });
    assert!(!work.exhausted, "{} work units", work.used);
    eprintln!(
        "loft crossing and numerical admission: {} work units",
        work.used
    );
}
