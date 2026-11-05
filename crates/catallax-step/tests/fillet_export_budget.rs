//! Fillet exports must either keep their certificates or refuse without output.
use catallax_features::evaluate;
use catallax_numbers::{BudgetPolicy, Decider, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_step::{StepOptions, write_step};
use serde_json::Value;

fn recipe(name: &str) -> Recipe {
    let case = include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|case| case["id"] == name)
        .unwrap();
    Recipe::parse(&case["recipe"].to_string()).unwrap()
}

#[test]
fn local_chart_budget_refuses_instead_of_writing_an_invalid_legacy_sphere() {
    std::thread::Builder::new()
        .stack_size(8 * 1024 * 1024)
        .spawn(|| {
            let ((), work) = with_operation_budget(1_000_000, || {
                let mut dec = Decider::new(BudgetPolicy {
                    max_work: 1_000_000,
                    ..Default::default()
                });
                let built = evaluate(&recipe("sphere-dome"), &mut dec).unwrap();
                assert!(built.model.check(built.result).is_empty());
                let mut out = b"untouched".to_vec();
                let refusal = write_step(
                    &built.model,
                    &[built.result],
                    &StepOptions::default(),
                    &mut dec,
                    &mut out,
                )
                .unwrap_err();
                assert_eq!(refusal.code, "step.budget-exhausted");
                assert_eq!(out, b"untouched");
            });
            assert!(
                !work.exhausted,
                "the failure must be local, not request-wide"
            );
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn rounded_top_and_uprights_export_within_one_million_work() {
    std::thread::Builder::new()
        .stack_size(8 * 1024 * 1024)
        .spawn(|| {
            let ((), work) = with_operation_budget(1_000_000, || {
                let mut dec = Decider::default();
                let built = evaluate(&recipe("roll-over-mitres"), &mut dec).unwrap();
                assert!(built.model.check(built.result).is_empty());
                let options = StepOptions::default();
                let mut out = Vec::new();
                let answer =
                    write_step(&built.model, &[built.result], &options, &mut dec, &mut out)
                        .unwrap();
                assert_eq!(answer.value.solids, 1);
                assert_eq!(answer.value.faces, 18);
                assert!(answer.value.max_deviation <= options.curve_tolerance);
                assert!(answer.marks.is_empty());
                assert!(!out.is_empty());
            });
            assert!(!work.exhausted, "{} work used", work.used);
        })
        .unwrap()
        .join()
        .unwrap();
}
