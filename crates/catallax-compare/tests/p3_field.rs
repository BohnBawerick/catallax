//! Full public P3 comparison. Scoreboard and verification allowances are distinct.
use std::time::Instant;

use catallax_compare::{Completion, Controller, Options, Report, Verdict};
use catallax_numbers::{OperationWork, with_operation_budget};
use catallax_recipe::Recipe;

fn check(limit: u64) -> (Report, OperationWork) {
    let recipe = Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    let started = Instant::now();
    let (report, work) = with_operation_budget(limit, || {
        let controller = Controller::default();
        assert_eq!(
            controller.begin().check(&recipe, &Options::default()),
            Completion::Published
        );
        controller.latest().unwrap()
    });
    eprintln!(
        "P3 comparison allowance={limit}, work={work:?}, elapsed={:?}, forward={}, reverse={}, findings={}",
        started.elapsed(),
        report.forward_vertices,
        report.reverse_samples,
        report.findings.len()
    );
    (report, work)
}

#[test]
fn p3_verification_comparison_has_two_way_coverage_without_unresolved_observations() {
    // The verification request default, not a replacement scoreboard allowance.
    let (report, work) = check(100_000_000);
    assert!(!work.exhausted);
    assert!(report.refusal.is_none(), "{report:#?}");
    assert!(report.not_comparable.is_empty(), "{report:#?}");
    assert!(report.forward_vertices > 0);
    assert!(report.reverse_samples > 0);
    assert!(report.findings.is_empty(), "{report:#?}");
    assert_eq!(report.verdict, Verdict::Agree);
    assert!(
        report
            .sampled_distance_upper
            .as_ref()
            .is_some_and(|bound| bound <= &Options::default().accuracy)
    );
}

#[test]
fn p3_scoreboard_comparison_resource_exhaustion_remains_a_refusal() {
    let (report, work) = check(1_000_000);
    assert!(work.exhausted);
    assert_eq!(report.refusal.unwrap().code, "numbers.operation-budget");
    assert_eq!(report.verdict, Verdict::Unresolved);
    assert!(!report.findings.is_empty());
    assert!(report.sampled_distance_upper.is_none());
}
