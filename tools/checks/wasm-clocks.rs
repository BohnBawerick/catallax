//! Standalone wasm smoke module compiled by ci/wasm-clocks, without native dev dependencies.
use catallax_numbers::{
    AheadSite, Decider, Helpers, OperationSpan, Real, RunAhead, ValuesAhead,
    record_operation_path, set_run_ahead_loops, set_run_ahead_wait, set_threads, threads,
    with_operation_budget, with_operation_log,
};

#[unsafe(no_mangle)]
pub extern "C" fn operation_log() {
    let decide = || {
        let mut dec = Decider::default();
        let span = OperationSpan::begin("brep", &"check", "check", &dec);
        record_operation_path("browser.check");
        let result = dec.sign(&Real::from_i64(2).sqrt());
        span.finish(&dec, None, None);
        format!("{result:?} {:?}", dec.stats())
    };
    let (plain, plain_work) = with_operation_budget(100_000, decide);
    for _ in 0..2 {
        let ((logged, entries), work) = with_operation_budget(100_000, || with_operation_log(decide));
        assert_eq!(logged, plain);
        assert_eq!(work.used, plain_work.used);
        assert!(work.used > 0 && !work.exhausted);
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].path, Some("browser.check"));
        assert_eq!(entries[0].work_after - entries[0].work_before, work.used);
        assert_eq!(entries[0].wall_ms, None);
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn run_ahead() {
    set_threads(std::num::NonZeroUsize::new(4));
    set_run_ahead_loops(true);
    set_run_ahead_wait(true);
    assert_eq!(threads(), 1);
    let mut dec = Decider::default();
    assert!(Helpers::new(&dec).is_none());
    static SITE: AheadSite = AheadSite::new();
    let mut ahead = RunAhead::new(&SITE, &dec, 3);
    let mut values = ValuesAhead::<usize>::new(&SITE, &dec, 3);
    for index in 0..3 {
        ahead.at(index, || panic!("browser requested worker body"));
        assert!(values.take(index, &mut dec, || panic!("browser requested worker value")).is_none());
    }
    set_run_ahead_loops(false);
    set_run_ahead_wait(false);
    set_threads(None);
}

#[unsafe(no_mangle)]
pub extern "C" fn display() {
    use catallax_field::{Plan, display::{Diagnostics, Options, mesh_plan, mesh_plan_reported}};
    use catallax_recipe::{Operation, Pose, Recipe, Step, StepId};
    let id = StepId::new("box").unwrap();
    let plan = Plan::compile(&Recipe::new("browser", vec![Step {
        id: id.clone(),
        op: Operation::Block { size: [2.into(), 3.into(), 4.into()], pose: Pose::identity() },
    }], id)).unwrap();
    let options = Options::default();
    let (plain, plain_work) = with_operation_budget(100_000_000, || {
        mesh_plan(&plan, &options, &mut Decider::default()).unwrap()
    });
    let mut diagnostics = Diagnostics::default();
    let (reported, work) = with_operation_budget(100_000_000, || {
        mesh_plan_reported(&plan, &options, &mut Decider::default(), &mut diagnostics).unwrap()
    });
    let mut plain_obj = Vec::new();
    plain.write_obj(&mut plain_obj).unwrap();
    let mut reported_obj = Vec::new();
    reported.write_obj(&mut reported_obj).unwrap();
    assert!(!reported.triangles.is_empty());
    assert_eq!(plain_obj, reported_obj);
    assert_eq!(plain_work.used, work.used);
    assert!(!work.exhausted);
    assert_eq!(diagnostics.total.work, work.used);
    assert_eq!(diagnostics.total.elapsed, None);
    assert_eq!(diagnostics.common.elapsed, None);
    for support in diagnostics.supports {
        for cost in [support.total, support.classification, support.queries, support.emission] {
            assert_eq!(cost.elapsed, None);
        }
    }
}
