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
fn a_rounded_spindle_band_needs_certified_axis_clearance() {
    std::thread::Builder::new()
        .stack_size(8 * 1024 * 1024)
        .spawn(|| {
            let mut case = include_str!("../../../parts/fillets/cases.jsonl")
                .lines()
                .map(|line| serde_json::from_str::<Value>(line).unwrap())
                .find(|case| case["id"] == "narrow-cone")
                .unwrap();
            case["recipe"]["steps"][0]["height"] = "4".into();
            case["recipe"]["steps"][1]["radius"] = "1999999999999999999/1000000000000000000".into();
            let recipe = Recipe::parse(&case["recipe"].to_string()).unwrap();
            let ((), work) = with_operation_budget(1_000_000, || {
                let mut dec = Decider::default();
                let built = evaluate(&recipe, &mut dec).unwrap();
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
                assert!(refusal.detail.contains("radial clearance"));
                assert_eq!(out, b"untouched");
            });
            assert!(!work.exhausted);
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn a_coarse_spindle_band_needs_certified_rounding_clearance() {
    std::thread::Builder::new()
        .stack_size(8 * 1024 * 1024)
        .spawn(|| {
            let mut case = include_str!("../../../parts/fillets/cases.jsonl")
                .lines()
                .map(|line| serde_json::from_str::<Value>(line).unwrap())
                .find(|case| case["id"] == "narrow-cone")
                .unwrap();
            case["recipe"]["steps"][0]["base_radius"] = "151/25".into();
            case["recipe"]["steps"][0]["height"] = "21/10".into();
            let recipe = Recipe::parse(&case["recipe"].to_string()).unwrap();
            let mut dec = Decider::default();
            let built = evaluate(&recipe, &mut dec).unwrap();
            let options = StepOptions {
                significant_digits: 3,
                curve_tolerance: 100.into(),
                ..StepOptions::default()
            };
            let mut out = b"untouched".to_vec();
            let refusal = write_step(&built.model, &[built.result], &options, &mut dec, &mut out)
                .unwrap_err();
            assert_eq!(refusal.code, "step.budget-exhausted");
            assert!(refusal.detail.contains("torus sheet band"), "{refusal:?}");
            assert_eq!(out, b"untouched");
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn fillet_exports_stay_certified_within_one_million_work() {
    for (name, faces) in [
        ("roll-over-mitres", 18),
        ("torus-slice", 4),
        ("unequal-ordered", 10),
        ("acute-5deg", 6),
        ("acute-halfdeg", 6),
        ("narrow-cone", 4),
        ("sphere-dome", 5),
    ] {
        std::thread::Builder::new()
            .stack_size(8 * 1024 * 1024)
            .spawn(move || {
                let ((), work) = with_operation_budget(1_000_000, || {
                    let mut dec = Decider::default();
                    let built = evaluate(&recipe(name), &mut dec).unwrap();
                    assert!(built.model.check(built.result).is_empty());
                    let options = StepOptions::default();
                    let mut out = Vec::new();
                    let answer =
                        write_step(&built.model, &[built.result], &options, &mut dec, &mut out)
                            .unwrap_or_else(|refusal| panic!("{name}: {refusal:?}"));
                    assert_eq!(answer.value.solids, 1);
                    assert_eq!(answer.value.faces, faces);
                    assert!(answer.value.max_deviation <= options.curve_tolerance);
                    assert!(answer.marks.is_empty());
                    assert!(!out.is_empty());
                });
                assert!(!work.exhausted, "{name}: {} work used", work.used);
            })
            .unwrap()
            .join()
            .unwrap();
    }
}

#[test]
fn rounded_corners_vertices_and_tangent_chains_export_within_one_million_work() {
    let mut failures = Vec::new();
    for name in [
        "box-all",
        "l-mixed-all",
        "cylinder-collapse",
        "cylinder-past-collapse",
        "tangent-boss",
        "tangent-pocket",
        "vertex-4-one",
    ] {
        let (answer, work) = std::thread::Builder::new()
            .stack_size(8 * 1024 * 1024)
            .spawn(move || {
                with_operation_budget(1_000_000, || {
                    let mut dec = Decider::default();
                    let built = evaluate(&recipe(name), &mut dec).map_err(|e| format!("{e:?}"))?;
                    assert!(built.model.check(built.result).is_empty());
                    let options = StepOptions::default();
                    let mut out = Vec::new();
                    let answer =
                        write_step(&built.model, &[built.result], &options, &mut dec, &mut out)
                            .map_err(|e| format!("{e:?}"))?;
                    assert_eq!(answer.value.solids, 1);
                    assert!(answer.value.max_deviation <= options.curve_tolerance);
                    assert!(answer.marks.is_empty());
                    assert!(!out.is_empty());
                    Ok::<_, String>(())
                })
            })
            .unwrap()
            .join()
            .unwrap();
        if answer.is_err() || work.exhausted {
            failures.push(format!("{name}: {answer:?}, {work:?}"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
