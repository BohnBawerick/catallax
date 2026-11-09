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

#[test]
fn complete_bracket_and_tray_parameter_edits_keep_construction_names() {
    for (fixture, step, parameter, value) in [
        (
            include_str!("../../../parts/P1-bracket.json"),
            "boss",
            "radius",
            "9",
        ),
        (
            include_str!("../../../parts/P4-shelled-tray.json"),
            "tray",
            "thickness",
            "1.5",
        ),
    ] {
        let before: Value = serde_json::from_str(fixture).unwrap();
        let mut after = before.clone();
        let changed = after["steps"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|candidate| candidate["id"] == step)
            .unwrap();
        changed[parameter] = json!(value);
        let report = check_with_work(before, after, 100_000_000);
        let names = &report["manifest"]["name_edit"];
        assert_eq!(names["brep"]["status"], "passed", "{step}: {report}");
        for kind in ["faces", "edges"] {
            assert_eq!(
                names["brep"]["summary"][kind]["lost"], 0,
                "{step}: {report}"
            );
            assert_eq!(
                names["brep"]["summary"][kind]["renamed"], 0,
                "{step}: {report}"
            );
        }
    }
}

#[test]
fn field_feature_gap_does_not_hide_complete_brep_edit_survival() {
    let before = include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|case| case["id"] == "explicit-setbacks")
        .unwrap()["recipe"]
        .clone();
    let mut after = before.clone();
    after["steps"][1]["radius"] = json!("3/4");
    let report = check_with_work(before, after, 100_000_000);
    let names = &report["manifest"]["name_edit"];
    assert_eq!(names["status"], "passed", "{report}");
    assert_eq!(names["brep"]["status"], "passed", "{report}");
    assert_eq!(names["field"]["status"], "waiting", "{report}");
    assert_eq!(names["waiting"][0]["phase"], "before", "{report}");
    assert_eq!(names["waiting"][1]["phase"], "after", "{report}");
    assert!(
        names["waiting"]
            .as_array()
            .unwrap()
            .iter()
            .all(|wait| wait["kernel"] == "field"
                && wait["refusal"]["code"] == "fillet.unsupported-overflow")
    );
}
