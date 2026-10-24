//! Both kernels must agree on sharp projection sectors and unrounded open rims.

use catallax_compare::{Completion, Controller, Verdict};
use catallax_numbers::with_operation_budget;
use catallax_recipe::Recipe;

#[test]
fn sharp_thickening_preserves_origins_distances_and_material_direction() {
    for text in [
        include_str!("../../catallax-features/tests/fixtures/sharp-thicken/planar-front.json"),
        include_str!("../../catallax-features/tests/fixtures/sharp-thicken/planar-back.json"),
        include_str!("../../catallax-features/tests/fixtures/sharp-thicken/planar-symmetric.json"),
        include_str!("../../catallax-features/tests/fixtures/sharp-thicken/cylinder-front.json"),
        include_str!("../../catallax-features/tests/fixtures/sharp-thicken/cylinder-back.json"),
        include_str!(
            "../../catallax-features/tests/fixtures/sharp-thicken/cylinder-symmetric.json"
        ),
        include_str!("../../catallax-features/tests/fixtures/sharp-thicken/cone-front.json"),
        include_str!("../../catallax-features/tests/fixtures/sharp-thicken/cone-back.json"),
        include_str!("../../catallax-features/tests/fixtures/sharp-thicken/cone-symmetric.json"),
        include_str!(
            "../../catallax-features/tests/fixtures/sharp-thicken/cone-cylinder-front.json"
        ),
    ] {
        let recipe = Recipe::parse(text).unwrap();
        let controller = Controller::default();
        let revision = controller.begin();
        let (completion, work) =
            with_operation_budget(100_000_000, || revision.check(&recipe, &Default::default()));
        assert_eq!(completion, Completion::Published);
        assert!(!work.exhausted, "{}", recipe.name);
        let report = controller.latest().unwrap();
        assert_eq!(
            report.verdict,
            Verdict::Agree,
            "{}: {:?}",
            recipe.name,
            report.findings
        );
        assert!(
            report.findings.is_empty(),
            "{}: {:?}",
            recipe.name,
            report.findings
        );
        assert!(
            report.not_comparable.is_empty(),
            "{}: {:?}",
            recipe.name,
            report.not_comparable
        );
        assert!(
            report
                .sampled_distance_upper
                .as_ref()
                .is_some_and(|d| d <= &catallax_compare::Options::default().accuracy)
        );
    }
}
