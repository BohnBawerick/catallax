//! Independent bore-shell distance, source-name and material-direction comparison.

use catallax_compare::{Completion, Controller, Options, Verdict};
use catallax_numbers::with_operation_budget;
use catallax_recipe::Recipe;

#[test]
fn bored_shell_and_thickness_edit_agree_on_distances_names_and_directions() {
    let fixture = include_str!("../../catallax-features/tests/fixtures/bored-shell.json");
    let mut source: serde_json::Value = serde_json::from_str(fixture).unwrap();
    for thickness in ["2", "3/2"] {
        source["steps"][2]["thickness"] = serde_json::json!(thickness);
        let recipe = Recipe::parse(&source.to_string()).unwrap();
        let controller = Controller::default();
        let revision = controller.begin();
        let (completion, work) =
            with_operation_budget(100_000_000, || revision.check(&recipe, &Options::default()));
        assert_eq!(completion, Completion::Published);
        assert!(!work.exhausted, "{work:?}");
        let report = controller.latest().unwrap();
        assert_eq!(
            report.verdict,
            Verdict::Agree,
            "{thickness}: {:?}; {:?}",
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
                .is_some_and(|upper| upper <= &Options::default().accuracy)
        );
    }
}
