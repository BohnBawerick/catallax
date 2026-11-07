//! A rounded ring where a cylinder passes through a sphere writes a closed
//! fitted band whose rails keep their exact sphere and cylinder walls.
use catallax_features::evaluate;
use catallax_numbers::{BudgetPolicy, Decider, Rational};
use catallax_recipe::Recipe;
use catallax_step::{CertifiedSplineStepOptions, write_rolling_junction_step};
use serde_json::Value;

fn recipe() -> Recipe {
    let case = include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|case| case["id"] == "sphere-cylinder")
        .unwrap();
    Recipe::parse(&case["recipe"].to_string()).unwrap()
}

#[test]
fn both_rings_write_within_the_file_tolerance() {
    std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(|| {
            let mut dec = Decider::new(BudgetPolicy {
                max_work: 200_000_000,
                ..Default::default()
            });
            let built = evaluate(&recipe(), &mut dec).unwrap();
            assert!(built.model.check(built.result).is_empty());
            let mut out = Vec::new();
            let options = CertifiedSplineStepOptions {
                product_name: "sphere-cylinder".into(),
                file_tolerance: Rational::from_ratio(1, 10_000_000).unwrap(),
                decimal_places: 15,
            };
            let report = write_rolling_junction_step(
                &built.model,
                built.result,
                &options,
                &mut dec,
                &mut out,
            )
            .unwrap()
            .value;
            assert!(report.max_deviation <= options.file_tolerance);
            let text = String::from_utf8(out).unwrap();
            if let Ok(path) = std::env::var("CATALLAX_SPHERE_RING_STEP") {
                std::fs::write(path, &text).unwrap();
            }
            // The sphere, a wall and a cap on each side, and a rational band
            // for each ring.
            assert_eq!(text.matches("ADVANCED_FACE(").count(), 7, "{report:?}");
            assert_eq!(text.matches("SPHERICAL_SURFACE(").count(), 1);
            assert_eq!(text.matches("CYLINDRICAL_SURFACE(").count(), 2);
            assert_eq!(text.matches("RATIONAL_B_SPLINE_SURFACE(").count(), 2);
        })
        .unwrap()
        .join()
        .unwrap();
}
