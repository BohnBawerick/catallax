//! Certified trimmed distances must resolve the real default comparison.

use catallax_compare::{Completion, Controller, Verdict};
use catallax_numbers::with_operation_budget;
use catallax_recipe::Recipe;

#[test]
fn p2_compares_in_both_directions_with_the_default_allowance() {
    let recipe = Recipe::parse(include_str!("../../../parts/P2-collar.json")).unwrap();
    let controller = Controller::default();
    let revision = controller.begin();
    let (completion, work) =
        with_operation_budget(100_000_000, || revision.check(&recipe, &Default::default()));
    assert_eq!(completion, Completion::Published);
    assert!(!work.exhausted);
    let report = controller.latest().unwrap();
    assert_eq!(report.verdict, Verdict::Agree, "{:?}", report.findings);
    assert!(report.findings.is_empty());
    assert!(report.not_comparable.is_empty());
    assert!(report.forward_vertices > 0 && report.reverse_samples > 0);
    assert!(report.sampled_distance_upper.as_ref().is_some_and(|upper| {
        upper <= &catallax_compare::Options::default().accuracy
    }));
}
