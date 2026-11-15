//! Check the actual source and decimal native joins through the public writer.
#[allow(dead_code)]
#[path = "../examples/numerical_uv_join.rs"]
mod fixture;

use catallax_numbers::Decider;
use catallax_topology::numerical::{NumericalBody, NumericalHistory};
use std::{
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "catallax-native-join-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
#[test]
fn a_real_source_join_cannot_hide_behind_declared_vertices() {
    let out = Scratch::new();
    let error = fixture::run(&out.0, "source-gap").unwrap_err();
    assert!(error.starts_with("source admission:"), "{error}");
    assert!(error.contains("trim.native-join-position"), "{error}");
    assert!(!out.0.join("answer.step").exists());
}
#[test]
fn decimal_rounding_must_recheck_the_actual_join() {
    let out = Scratch::new();
    let error = fixture::run(&out.0, "rounded-gap").unwrap_err();
    assert!(error.starts_with("decimal admission:"), "{error}");
    assert!(error.contains("trim.native-join-position"), "{error}");
    assert!(!out.0.join("answer.step").exists());
}
#[test]
fn matching_joins_and_current_tagged_handovers_remain_usable() {
    let out = Scratch::new();
    let result = fixture::run(&out.0, "valid").unwrap();
    assert_eq!(result["status"], "numerical");
    let mut dec = Decider::default();
    let body: NumericalBody = catallax_numbers::transfer::decode(
        &std::fs::read(out.0.join("body.handover")).unwrap(),
        &mut dec,
    )
    .unwrap();
    assert_eq!(body.native_trims().len(), 2);
    assert_eq!(body.history().len(), 1);
    let history: NumericalHistory = catallax_numbers::transfer::decode(
        &std::fs::read(out.0.join("history.handover")).unwrap(),
        &mut dec,
    )
    .unwrap();
    assert_eq!(history.native_trims.len(), 2);
}
