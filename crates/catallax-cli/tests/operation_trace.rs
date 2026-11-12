//! Operation outcomes keep their established dump when trace metadata is added.
use catallax_numbers::{
    Decider, OperationSpan, record_operation_path, with_operation_budget, with_operation_log,
};

#[test]
fn path_metadata_preserves_the_exact_operation_outcome_dump() {
    for path in [None, Some("distance.shell.sphere-radius")] {
        let (((), mut entries), work) = with_operation_budget(1, || {
            with_operation_log(|| {
                let dec = Decider::default();
                let span = OperationSpan::begin("brep", &"op_467", "sketch", &dec);
                if let Some(path) = path {
                    record_operation_path(path);
                }
                span.finish(&dec, None, None);
            })
        });
        assert_eq!(work.used, 0);
        assert!(!work.exhausted);
        assert_eq!(entries.len(), 1);
        let entry = &mut entries[0];
        assert_eq!(entry.path, path);
        entry.wall_ms = 0.0;
        assert_eq!(
            format!("{entry:#?}\n"),
            include_str!("fixtures/operation-entry-debug.txt")
        );
    }
}
