//! Full mass queries for mixed fillet corners and their radius and placement edits.

use catallax_features::evaluate;
use catallax_measure::{Accuracy, Enclosure, mass_properties};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::{Value, json};

fn source() -> Value {
    include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|case| case["id"] == "l-mixed-all")
        .unwrap()["recipe"]
        .clone()
}

fn close(actual: &Enclosure, expected: &Value, tolerance: &str) {
    let expected = Rational::parse(&expected.to_string()).unwrap();
    let tolerance = Rational::parse(tolerance).unwrap();
    assert!(
        actual.lo >= &expected - &tolerance && actual.hi <= &expected + &tolerance,
        "{actual:?} differs from {expected} by more than {tolerance}"
    );
}

#[test]
fn mixed_corner_mass_matches_occt_under_the_work_target() {
    // Generated with OCCT 7.9.3 by reference_mixed_fillet_mass.py.
    let references: Vec<Value> =
        serde_json::from_str(include_str!("fixtures/mixed-fillet-mass-occt.json")).unwrap();
    let accuracy = Accuracy {
        absolute: Rational::parse("1e-7").unwrap(),
        relative: Rational::zero(),
    };
    for reference in references {
        let mut source = source();
        source["steps"][3]["radius"] = json!(reference["radius"].to_string());
        if reference["posed"].as_bool().unwrap() {
            for i in [0, 1] {
                source["steps"][i]["pose"] = json!({
                    "rotation": ["1", "2", "3", "4"],
                    "translation": ["7/3", "-11/7", "13/5"]
                });
            }
        }
        let recipe = Recipe::parse(&source.to_string()).unwrap();
        let (built, work) = with_operation_budget(20_000_000, || {
            evaluate(&recipe, &mut Decider::default()).unwrap()
        });
        assert!(!work.exhausted);
        assert!(built.model.check(built.result).is_empty());
        let (answer, work) = with_operation_budget(20_000_000, || {
            mass_properties(
                &built.model,
                built.result,
                &accuracy,
                &mut Decider::default(),
            )
        });
        eprintln!(
            "radius={} posed={} mass={work:?}",
            reference["radius"], reference["posed"]
        );
        assert!(!work.exhausted);
        let answer = answer.unwrap().value;
        close(&answer.volume, &reference["volume"], "1e-6");
        close(&answer.area, &reference["area"], "1e-6");
        for (actual, expected) in answer
            .centroid
            .as_ref()
            .unwrap()
            .iter()
            .zip(reference["centroid"].as_array().unwrap())
        {
            close(actual, expected, "1e-7");
        }
        let tensor = answer.inertia.as_ref().unwrap();
        for (i, row) in tensor.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                close(value, &reference["inertia"][i][j], "1e-5");
                assert!(value.width() <= accuracy.absolute);
            }
        }
    }
}
