//! End-condition dispatch reaches the unchanged two-kernel comparison.
use catallax_compare::{Completion, Controller, Options, Verdict};
use catallax_recipe::Recipe;

#[test]
fn curved_first_hit_runs_both_comparison_passes() {
    let recipe = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/curved-limits/sphere-entry-next.json"
    ))
    .unwrap();
    let controller = Controller::default();
    assert_eq!(
        controller.begin().check(&recipe, &Options::default()),
        Completion::Published
    );
    let report = controller.latest().unwrap();
    assert!(report.not_comparable.is_empty());
    assert_eq!(report.verdict, Verdict::Agree, "{report:?}");
    assert!(report.forward_vertices > 0);
    assert!(report.reverse_samples > 0);
    assert!(report.refusal.is_none());
}
