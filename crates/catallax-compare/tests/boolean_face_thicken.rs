//! Whole-result sampled comparison is separate from the exact volume/point oracles.

use catallax_compare::{Completion, Controller, Verdict};
use catallax_numbers::with_operation_budget;
use catallax_recipe::Recipe;

#[test]
fn retained_polygonal_face_thickens_have_no_comparison_findings() {
    for text in [
        include_str!("../../catallax-features/tests/fixtures/boolean-face-thicken/l-front.json"),
        include_str!("../../catallax-features/tests/fixtures/boolean-face-thicken/l-back.json"),
        include_str!("../../catallax-features/tests/fixtures/boolean-face-thicken/l-tilted.json"),
        include_str!(
            "../../catallax-features/tests/fixtures/boolean-face-thicken/l-merged-origins.json"
        ),
        include_str!("../../catallax-features/tests/fixtures/boolean-face-thicken/hole-front.json"),
        include_str!("../../catallax-features/tests/fixtures/boolean-face-thicken/hole-back.json"),
        include_str!("../../catallax-features/tests/fixtures/boolean-face-thicken/symmetric.json"),
        include_str!(
            "../../catallax-features/tests/fixtures/boolean-face-thicken/reversed-wall.json"
        ),
        include_str!("../../catallax-features/tests/fixtures/boolean-face-thicken/split-near.json"),
    ] {
        let recipe = Recipe::parse(text).unwrap();
        let controller = Controller::default();
        let revision = controller.begin();
        let (completion, work) =
            with_operation_budget(100_000_000, || revision.check(&recipe, &Default::default()));
        assert_eq!(completion, Completion::Published);
        assert!(!work.exhausted);
        let report = controller.latest().unwrap();
        assert_eq!(
            report.verdict,
            Verdict::Agree,
            "{}: {:?}; {:?}",
            recipe.name,
            report.findings,
            report.not_comparable
        );
        assert!(report.findings.is_empty());
        assert!(
            report.not_comparable.is_empty(),
            "{:?}",
            report.not_comparable
        );
        assert!(
            report
                .sampled_distance_upper
                .as_ref()
                .is_some_and(|u| u <= &catallax_compare::Options::default().accuracy)
        );
    }
}
