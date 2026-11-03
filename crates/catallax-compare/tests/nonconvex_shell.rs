//! Independent planar-shell distance, source-name and material-direction comparison.

use catallax_compare::{Completion, Controller, Options, Verdict};
use catallax_numbers::with_operation_budget;
use catallax_recipe::Recipe;
use serde_json::{Value, json};

fn document() -> Value {
    serde_json::from_str(include_str!(
        "../../catallax-features/tests/fixtures/nonconvex-shell.json"
    ))
    .unwrap()
}
fn agree(source: &Value, case: &str) {
    let recipe = Recipe::parse(&source.to_string()).unwrap();
    let controller = Controller::default();
    let revision = controller.begin();
    let options = Options::default();
    let (completion, work) =
        with_operation_budget(100_000_000, || revision.check(&recipe, &options));
    assert_eq!(completion, Completion::Published, "{case}: {work:?}");
    let report = controller.latest().unwrap();
    assert!(
        !work.exhausted,
        "{case}: {work:?}; {:?}; {:?}",
        report.findings, report.not_comparable
    );
    assert_eq!(
        report.verdict,
        Verdict::Agree,
        "{case}: {:?}; {:?}",
        report.findings,
        report.not_comparable
    );
    assert!(report.findings.is_empty(), "{case}: {:?}", report.findings);
    assert!(
        report.not_comparable.is_empty(),
        "{case}: {:?}",
        report.not_comparable
    );
    assert!(
        report
            .sampled_distance_upper
            .as_ref()
            .is_some_and(|upper| upper <= &options.accuracy),
        "{case}: {:?}",
        report.sampled_distance_upper
    );
}
fn chain(points: &[[i64; 2]], prefix: &str) -> Value {
    json!(
        points
            .iter()
            .enumerate()
            .map(|(i, _)| json!({
                "kind":"line", "id":format!("{prefix}-{i}"),
                "to":points[(i + 1) % points.len()].map(|n| n.to_string())
            }))
            .collect::<Vec<_>>()
    )
}
fn open_except(source: &mut Value, segment: &str) {
    let mut open = vec![
        json!({"face":{"step":"stock","role":"start"}}),
        json!({"face":{"step":"stock","role":"end"}}),
    ];
    for edge in source["steps"][0]["loops"]
        .as_array()
        .unwrap()
        .iter()
        .flat_map(|cycle| cycle.as_array().unwrap())
    {
        if edge["id"] != segment {
            open.push(json!({"face":{"step":"stock","role":{"side":{"segment":edge["id"]}}}}));
        }
    }
    source["steps"][2]["open"] = json!(open);
}
fn holed_document() -> Value {
    let mut source = document();
    source["steps"][0]["loops"] = json!([
        chain(&[[0, 0], [10, 0], [10, 10], [0, 10]], "outer"),
        chain(&[[3, 3], [3, 7], [7, 7], [7, 3]], "hole")
    ]);
    source["steps"][2]["open"] = json!([
        {"face":{"step":"stock","role":"start"}},
        {"face":{"step":"stock","role":"end"}}
    ]);
    source
}

#[test]
fn reentrant_tray_agrees_independently() {
    agree(&document(), "original reentrant tray");
}

#[test]
fn both_sided_kept_wall_agrees_independently() {
    let mut source = document();
    open_except(&mut source, "notch-bottom");
    agree(&source, "both-sided finite kept wall");
}

#[test]
fn cross_gap_opening_agrees_independently() {
    let mut source = document();
    source["steps"][0]["loops"] = json!([chain(
        &[
            [0, 0],
            [8, 0],
            [8, 8],
            [5, 8],
            [5, 2],
            [3, 2],
            [3, 8],
            [0, 8]
        ],
        "u",
    )]);
    open_except(&mut source, "u-5");
    source["steps"][2]["thickness"] = json!("3");
    agree(&source, "material across the separated-arm gap");
}

#[test]
fn holed_shell_and_consumption_agree_independently() {
    let mut source = holed_document();
    for thickness in ["3/2", "2"] {
        source["steps"][2]["thickness"] = json!(thickness);
        agree(&source, thickness);
    }
}

#[test]
fn thickness_edit_retains_independent_names_and_directions() {
    let mut source = holed_document();
    for thickness in ["1", "5/4"] {
        source["steps"][2]["thickness"] = json!(thickness);
        agree(&source, thickness);
    }
}
