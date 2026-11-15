//! The actual loft crossing retains its native curves and physical patch hole.

use catallax_numbers::{
    Decider, NumericalCheck, Rational, current_operation_work, with_operation_budget,
};
use catallax_recipe::Recipe;
use catallax_step::numerical::{self as step, NumericalStepOptions};
use catallax_topology::numerical::{self, NumericalOptions};
use serde_json::{Value, json};

#[test]
fn the_loft_crossing_reaches_approximate_step_with_its_hole() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut input: Value = serde_json::from_str(include_str!(
            "../../catallax-features/tests/fixtures/loft-spline-prism.json"
        ))
        .unwrap();
        let steps = input["steps"].as_array_mut().unwrap();
        steps.push(json!({"id":"box","op":"block","size":["10","10","1"],"pose":{"rotation":["1","0","0","0"],"translation":["-3","1999999/1000000","1"]}}));
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
        let body = numerical::build(
            definition,
            NumericalOptions {
                tolerance_mm: Rational::parse("1/10000000").unwrap(),
            },
            &mut dec,
        )
        .unwrap();
        let before = current_operation_work().used;
        let file = step::build(&body, &NumericalStepOptions::default(), &mut dec).unwrap();
        let text = std::str::from_utf8(file.bytes()).unwrap();
        assert!(text.contains("numerical approximate"));
        assert!(text.contains("complete_body_bound=unestablished"));
        assert_eq!(text.matches("FACE_OUTER_BOUND(").count(), 7);
        assert_eq!(text.matches("=FACE_BOUND(").count(), 1);
        assert_eq!(text.matches("=EDGE_CURVE(").count(), 12);
        assert!(text.contains("RATIONAL_B_SPLINE_SURFACE("));
        assert!(text.contains("RATIONAL_B_SPLINE_CURVE("));
        assert!(text.contains("=PCURVE("));
        assert_eq!(file.report().counts, [7, 12, 8]);
        assert_eq!(file.report().source_definition.faces[0].loops.len(), 2);
        assert_eq!(
            file.report().source_evidence.tolerance_mm(),
            body.evidence().tolerance_mm()
        );
        assert!(file.report().spline_rounding.is_some());
        assert!(file.report().emitted_evidence.checks().iter().any(|check| matches!(check,NumericalCheck::Passed {name,..} if name=="trim.sampled-hole-contained-disjoint")));
        let output = std::env::temp_dir().join(format!(
            "catallax-numerical-curved-{}.step",
            std::process::id()
        ));
        std::fs::write(&output, file.bytes()).unwrap();
        let reread = std::fs::read(&output).unwrap();
        std::fs::remove_file(&output).unwrap();
        file.check_written(&reread, &mut dec).unwrap();
        eprintln!(
            "loft admission={before} approximate STEP={} work units",
            current_operation_work().used - before
        );
        let coarse = NumericalStepOptions {
            significant_digits: 1,
            max_scalar_rounding_error_mm: Rational::one(),
            ..Default::default()
        };
        assert!(step::build(&body, &coarse, &mut dec).is_err());
    });
    assert!(!work.exhausted, "{} work units", work.used);
}
