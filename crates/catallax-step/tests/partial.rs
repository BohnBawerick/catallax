//! Public partial-periodic export and unresolved-rounding boundaries.
mod common;
use catallax_numbers::with_operation_budget;
use common::*;

#[test]
fn cylinder_sphere_trim_preserves_the_existing_request_allowance() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut model = Model::new();
        let mut dec = Decider::default();
        let cylinder =
            cylinder(&mut model, n(1), r("11/4"), &translated("0", "0", "-11/8")).unwrap();
        let sphere = sphere(&mut model, r("7/4"), &translated("1/8", "0", "0")).unwrap();
        let body = catallax_boolean::boolean(
            &mut model,
            cylinder,
            sphere,
            catallax_boolean::BoolOp::Intersect,
            &mut dec,
        )
        .unwrap()
        .body;
        let mut bytes = Vec::new();
        let answer = write_step(
            &model,
            &[body],
            &StepOptions::default(),
            &mut dec,
            &mut bytes,
        )
        .unwrap();
        let report = answer.value;
        assert!(answer.marks.is_empty(), "{:?}", answer.marks);
        assert_eq!(
            (report.solids, report.faces, report.edges, report.vertices),
            (1, 5, 9, 6)
        );
        assert_eq!((report.added_edges, report.added_vertices), (3, 2));
        assert!(report.max_deviation <= r("1e-7"));
        let original = model.counts(body);
        assert_eq!(
            (original.faces, original.edges, original.vertices),
            (5, 6, 4)
        );
        assert!(model.check(body).is_empty());
    });
    assert!(work.used <= 100_000_000);
}

#[test]
fn unresolved_sliver_keeps_legacy_topology_instead_of_periodic_chart_coverage() {
    let (mut model, a, b) = boolean_pair(13);
    let mut dec = Decider::default();
    let body = catallax_boolean::boolean(
        &mut model,
        a,
        b,
        catallax_boolean::BoolOp::Intersect,
        &mut dec,
    )
    .unwrap()
    .body;
    let mut bytes = Vec::new();
    let report = write_step(
        &model,
        &[body],
        &StepOptions::default(),
        &mut dec,
        &mut bytes,
    )
    .unwrap()
    .value;
    assert_eq!((report.faces, report.edges, report.vertices), (4, 6, 4));
    assert_eq!((report.added_edges, report.added_vertices), (0, 0));
    assert!(report.regularizations.is_empty());
}
