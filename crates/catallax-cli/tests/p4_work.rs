//! Complete P4 export keeps its proven bytes as its work cost falls.
use catallax_cli::verification::{Request, answer};
use sha2::{Digest, Sha256};

#[test]
fn shelled_tray_preserves_step_with_less_work() {
    let out = std::env::temp_dir().join(format!("catallax-p4-work-{}", std::process::id()));
    let response = answer(&Request {
        format: "catallax-request/1".into(),
        source: "recipe".into(),
        problem: serde_json::from_str(include_str!("../../../parts/P4-shelled-tray.json")).unwrap(),
        out: out.clone(),
        max_work: 100_000_000,
    })
    .unwrap();
    let bytes = std::fs::read(out.join("answer.step")).unwrap();
    std::fs::remove_dir_all(out).unwrap();
    assert_eq!(response["manifest"]["status"], "answered");
    assert_eq!(response["manifest"]["exact"], true);
    assert!(!response["work"]["exhausted"].as_bool().unwrap());
    assert!(response["work"]["used"].as_u64().unwrap() <= 1_366_808);
    assert_eq!(
        format!("{:x}", Sha256::digest(&bytes)),
        "f9a8b9a15a576a4989690f75e1bd9f697e052299b80b33c1bfb4fdc4634c7381"
    );
}
