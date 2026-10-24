//! Public exchange checks actual rebuilt names, not display indices.
use catallax_cli::verification::{Request, answer};
use serde_json::{Value, json};

fn block() -> Value {
    json!({"format":"catallax-recipe/3","name":"edit","steps":[{"id":"stock","op":"block","size":["4","5","6"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"stock"})
}
fn check(before: Value, after: Value) -> Value {
    check_with_work(before, after, 1_000_000)
}
fn check_with_work(before: Value, after: Value, max_work: u64) -> Value {
    let dir = std::env::temp_dir().join(format!(
        "catallax-name-edit-{}-{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    let result = answer(&Request {
        format: "catallax-request/1".into(),
        source: "name-edit".into(),
        problem: json!({"format":"catallax-name-edit/1","recipe":before,"edited_recipe":after}),
        out: dir.clone(),
        max_work,
    })
    .unwrap();
    std::fs::remove_dir_all(dir).unwrap();
    result
}
#[test]
fn dimension_edit_rebuilds_and_resolves_all_faces_and_edges() {
    let before = block();
    let mut after = before.clone();
    after["steps"][0]["size"][2] = json!("7");
    let report = check(before, after);
    assert_eq!(
        report["manifest"]["name_edit"]["status"], "passed",
        "{report}"
    );
    assert_eq!(report["manifest"]["name_edit"]["brep"]["faces_checked"], 6);
    assert_eq!(report["manifest"]["name_edit"]["brep"]["edges_checked"], 12);
}
#[test]
fn deleted_profile_identity_fails_even_when_geometry_is_unchanged() {
    let mut before = block();
    before["steps"] = json!([
        {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[{"kind":"line","id":"bottom","to":["4","0"]},{"kind":"line","id":"right","to":["4","5"]},{"kind":"line","id":"top","to":["0","5"]},{"kind":"line","id":"left","to":["0","0"]}]]},
        {"id":"stock","op":"extrude","profile":"profile","distance":"6"}
    ]);
    let mut after = before.clone();
    after["steps"][0]["loops"][0][0]["id"] = json!("replacement");
    let report = check(before, after);
    assert_eq!(
        report["manifest"]["name_edit"]["status"], "failed",
        "{report}"
    );
    assert!(
        report["manifest"]["name_edit"]["brep"]["missing_faces"]
            .as_array()
            .unwrap()
            .iter()
            .any(|f| f["name"].to_string().contains("bottom"))
    );
    assert!(
        report["manifest"]["name_edit"]["brep"]["missing_edges"]
            .as_array()
            .unwrap()
            .iter()
            .any(|f| f["name"].to_string().contains("bottom"))
    );
}
#[test]
fn unfinished_feature_reports_wait_with_step_and_operation() {
    let mut before = block();
    before["steps"]
        .as_array_mut()
        .unwrap()
        .push(json!({"id":"tray","op":"shell","body":"stock","thickness":"1","open":[]}));
    before["result"] = json!("tray");
    let mut after = before.clone();
    after["steps"][1]["thickness"] = json!("0.75");
    let report = check(before, after);
    let edit = &report["manifest"]["name_edit"];
    assert_eq!(edit["status"], "waiting", "{report}");
    assert_eq!(edit["waiting"][0]["refusal"]["step"], "tray");
    assert_eq!(edit["waiting"][0]["refusal"]["operation"], "shell");
}
#[test]
fn no_edit_cannot_pass_the_edit_gate() {
    let before = block();
    let mut after = before.clone();
    for name in ["edit", "renamed"] {
        after["name"] = json!(name);
        let report = check(before.clone(), after.clone());
        assert_eq!(report["manifest"]["name_edit"]["code"], "names.no-edit");
    }
}

#[test]
fn seam_height_edit_keeps_names_when_merged_caps_separate() {
    let before: Value = serde_json::from_str(include_str!("../../../parts/P0/seam.json")).unwrap();
    let mut after = before.clone();
    after["steps"][0]["height"] = json!("12");
    let report = check_with_work(before, after, 100_000_000);
    assert_eq!(
        report["manifest"]["name_edit"]["status"], "passed",
        "{report}"
    );
}

#[test]
fn collar_chamfer_parameter_edit_keeps_faces_and_edges_in_both_kernels() {
    let before: Value =
        serde_json::from_str(include_str!("../../../parts/P2-collar.json")).unwrap();
    let mut after = before.clone();
    after["steps"][2]["distance"] = json!("0.75");
    let report = check_with_work(before, after, 100_000_000);
    assert_eq!(
        report["manifest"]["name_edit"]["status"], "passed",
        "{report}"
    );
}
