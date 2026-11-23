//! Fitted sweep bounds through the public verification request.

use catallax_cli::verification::{Request, answer};
use catallax_numbers::Rational;
use serde_json::{Value, json};

#[test]
fn compares_a_fitted_quarter_helix_at_the_existing_accuracy() {
    let out = std::env::temp_dir().join(format!("catallax-fitted-compare-{}", std::process::id()));
    let problem = serde_json::from_str(include_str!("fixtures/fitted-quarter-helix.json")).unwrap();
    let mut request = Request {
        format: "catallax-request/1".into(),
        source: "compare".into(),
        problem,
        out: out.clone(),
        max_work: 20_000_000,
    };
    let response = answer(&request).unwrap();
    let manifest = &response["manifest"];
    let comparison = &manifest["comparison"];
    assert_eq!(manifest["status"], "answered", "{manifest}");
    assert_eq!(manifest["exact"], false);
    assert_eq!(response["work"]["exhausted"], false);
    assert_eq!(comparison["verdict"], "agree");
    assert_eq!(comparison["proof"], "fitted-sweep-distance");
    assert_eq!(comparison["faces"], 10);
    assert!(comparison["display_vertices"].as_u64().unwrap() > 0);
    let read = |key: &str| Rational::parse(comparison[key].as_str().unwrap()).unwrap();
    let source = read("source_upper_mm");
    let display = read("display_upper_mm");
    let composed = read("composed_upper_mm");
    assert!(source <= read("modelling_tolerance_mm"));
    assert!(display <= Rational::parse("0.005").unwrap());
    assert!(composed >= source && composed >= display);
    assert!(composed <= Rational::parse("0.01").unwrap());

    request.problem["tolerance_mm"] = json!("0");
    let invalid: Value = answer(&request).unwrap();
    assert_eq!(invalid["manifest"]["status"], "invalid-question");
    assert_eq!(invalid["manifest"]["refusal"]["code"], "cli.fitted-sweep");
    std::fs::remove_dir_all(out).unwrap();
}
