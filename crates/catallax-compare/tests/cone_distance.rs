//! Conical distance features must agree on origins, gaps and material direction.

use catallax_compare::{Completion, Controller, Verdict};
use catallax_numbers::with_operation_budget;
use catallax_recipe::Recipe;

#[test]
fn conical_offsets_shells_and_sweeps_agree_without_direction_or_origin_findings() {
    for text in [
        include_str!("../../catallax-features/tests/fixtures/cone-distance/frustum-inward.json"),
        include_str!("../../catallax-features/tests/fixtures/cone-distance/thicken-front.json"),
        include_str!("../../catallax-features/tests/fixtures/cone-distance/thicken-back.json"),
        include_str!("../../catallax-features/tests/fixtures/cone-distance/shell-open-wide.json"),
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
            "{}: {:?}",
            recipe.name,
            report.findings
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
