//! Exercise the persistent public command through stdin and output artifacts.
use serde_json::{Value, json};
use std::{
    fs,
    io::Write,
    process::{Command, Stdio},
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);

fn question() -> Value {
    json!({"format":"exact-brep-question/2","id":"stream-test",
        "operands":{"a":{"shape":"block","size":["1","2","3"]}},
        "steps":[],"result":"a","limits":{"deadline_ms":0}})
}

#[test]
fn stream_ignores_wall_limit_and_resets_operation_budget() {
    let dir = std::env::temp_dir().join(format!(
        "catallax-stream-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&dir).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
        .arg("stream")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let mut input = child.stdin.take().unwrap();
    for (i, limit) in [100_000_000u64, 0, 100_000_000].into_iter().enumerate() {
        writeln!(
            input,
            "{}",
            json!({"format":"catallax-request/1","source":"question",
            "problem":question(),"max_work":limit,"out":dir.join(i.to_string())})
        )
        .unwrap();
    }
    drop(input);
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    let rows: Vec<Value> = String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(rows.len(), 3);
    assert_eq!(rows[0]["manifest"]["status"], "answered");
    assert_eq!(
        rows[1]["manifest"]["refusal"]["code"],
        "numbers.operation-budget"
    );
    assert_eq!(rows[1]["operation_log"]["operations"][0]["id"], "a");
    assert_eq!(
        rows[1]["operation_log"]["operations"][0]["refusal"]["code"],
        "numbers.operation-budget"
    );
    assert!(!dir.join("1/answer.step").exists());
    assert_eq!(rows[0]["work"], rows[2]["work"]);
    let mut a = rows[0]["manifest"].clone();
    let mut b = rows[2]["manifest"].clone();
    a.as_object_mut().unwrap().remove("stats");
    b.as_object_mut().unwrap().remove("stats");
    assert_eq!(a, b);
    assert_eq!(
        fs::read(dir.join("0/answer.step")).unwrap(),
        fs::read(dir.join("2/answer.step")).unwrap()
    );
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn operation_scope_restores_after_panic_and_is_thread_local() {
    use catallax_numbers::{Decider, with_operation_budget};
    let _ = std::panic::catch_unwind(|| with_operation_budget(0, || panic!("test")));
    assert!(Decider::default().check_deadline().is_ok());
    let handles: Vec<_> = [0, 5]
        .into_iter()
        .map(|limit| {
            std::thread::spawn(move || {
                let (result, work) =
                    with_operation_budget(limit, || Decider::default().check_deadline());
                (result.is_ok(), work.used, work.exhausted)
            })
        })
        .collect();
    let results: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();
    assert_eq!(results, [(false, 0, true), (true, 1, false)]);
}

#[test]
fn request_does_not_inherit_callers_numeric_thread_state() {
    use catallax_cli::verification::{Request, answer};
    use catallax_numbers::with_operation_budget;
    let dir = std::env::temp_dir().join(format!(
        "catallax-isolated-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let request = Request {
        format: "catallax-request/1".into(),
        source: "question".into(),
        problem: question(),
        out: dir.clone(),
        max_work: 100_000_000,
    };
    let cold = answer(&request).unwrap();
    let (nested, outer) = with_operation_budget(0, || answer(&request).unwrap());
    assert_eq!(nested["manifest"]["status"], "answered");
    assert_eq!(nested["work"], cold["work"]);
    assert_eq!(outer.used, 0);
    assert!(!outer.exhausted);
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn exported_recipe_replays_geometry_work_and_operation_outcomes() {
    use catallax_cli::verification::{Request, answer};
    let root = std::env::temp_dir().join(format!(
        "catallax-export-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let recipe = json!({"format":"catallax-recipe/3","name":"exported",
        "steps":[{"id":"base","op":"block","size":["4","5","6"],
        "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],
        "result":"base"});
    for unsupported in [false, true] {
        let mut problem = recipe.clone();
        if unsupported {
            problem["steps"].as_array_mut().unwrap().extend([
                json!({"id":"offset","op":"offset","body":"base","distance":"1"}),
                json!({"id":"later","op":"offset","body":"offset","distance":"1"}),
            ]);
            problem["result"] = json!("later");
        }
        let directory = root.join(unsupported.to_string());
        let original = answer(&Request {
            format: "catallax-request/1".into(),
            source: "recipe".into(),
            problem,
            out: directory.clone(),
            max_work: 12_345_678,
        })
        .unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
            .args([
                "replay",
                directory.join("question.json").to_str().unwrap(),
                "--out",
                directory.join("replay").to_str().unwrap(),
            ])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let replay: Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(original["manifest"], replay["manifest"]);
        assert_eq!(original["work"], replay["work"]);
        let mut logs = [
            original["operation_log"].clone(),
            replay["operation_log"].clone(),
        ];
        for log in &mut logs {
            for operation in log["operations"].as_array_mut().unwrap() {
                operation.as_object_mut().unwrap().remove("wall_ms");
            }
        }
        assert_eq!(logs[0], logs[1]);
        let operations = logs[0]["operations"].as_array().unwrap();
        assert_eq!(operations[0]["id"], "base");
        assert_eq!(operations[0]["status"], "answered");
        if unsupported {
            assert_eq!(
                original["manifest"]["refusal"]["code"],
                "features.unsupported-operation"
            );
            assert_eq!(operations[1]["id"], "offset");
            assert_eq!(
                operations[1]["refusal"]["code"],
                "features.unsupported-operation"
            );
            assert!(!operations.iter().any(|entry| entry["id"] == "later"));
            assert!(!directory.join("answer.step").exists());
        } else {
            assert_eq!(original["manifest"]["status"], "answered");
            assert_eq!(
                fs::read(directory.join("answer.step")).unwrap(),
                fs::read(directory.join("replay/answer.step")).unwrap()
            );
        }
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn construction_refusal_retains_boolean_context_after_measurement_split() {
    use catallax_cli::verification::{Request, answer};
    let dir = std::env::temp_dir().join(format!(
        "catallax-construction-refusal-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let response = answer(&Request {
        format: "catallax-request/1".into(),
        source: "question".into(),
        problem: json!({
            "format":"exact-brep-question/2","id":"crossing-cones",
            "operands":{
                "a":{"shape":"cone","base_radius":"5","top_radius":"2","height":"10"},
                "b":{"shape":"cone","base_radius":"5","top_radius":"2","height":"10",
                    "pose":{"rotation":["2","1","1","0"],"translation":["2","1","1"]}}
            },
            "steps":[{"id":"joined","op":"unite","args":["a","b"]}],
            "result":"joined","measure":{"properties":[]}
        }),
        out: dir.clone(),
        max_work: 100_000_000,
    })
    .unwrap();
    assert_eq!(response["manifest"]["status"], "refused");
    assert_eq!(
        response["manifest"]["refusal"]["code"],
        "intersect.traced-pair-required"
    );
    assert_eq!(response["manifest"]["refusal"]["step"], "joined");
    assert_eq!(response["manifest"]["refusal"]["stage"], "intersect");
    let operations = response["operation_log"]["operations"].as_array().unwrap();
    let failed = operations.last().unwrap();
    assert_eq!(failed["id"], "joined");
    assert_eq!(failed["refusal"]["code"], "intersect.traced-pair-required");
    assert!(
        !operations
            .iter()
            .any(|entry| entry["operation"] == "measure")
    );
    assert!(!dir.join("answer.step").exists());
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn declared_face_checks_survive_replay_and_report_radius_edit_origins() {
    use catallax_cli::verification::{Request, answer};
    let dir = std::env::temp_dir().join(format!(
        "catallax-face-check-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let problem = json!({
        "format":"catallax-face-check/1",
        "recipe":{"format":"catallax-recipe/3","name":"tube","steps":[{
            "id":"tube","op":"torus","major":"4","minor":"1.0002",
            "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}
        }],"result":"tube"},
        "definitions":[{"face":{"step":"tube"},"definition":{
            "kind":"fillet","radius":"1",
            "path":{"kind":"circle","center":["0","0","0"],"normal":["0","0","1"],
                "radius":"4","reference":["1","0","0"]}
        }}],
        "edit":{"step":"tube","parameter":"minor","value":"1.1"}
    });
    let response = answer(&Request {
        format: "catallax-request/1".into(),
        source: "recipe".into(),
        problem,
        out: dir.clone(),
        max_work: 100_000_000,
    })
    .unwrap();
    assert_eq!(response["manifest"]["status"], "answered", "{response}");
    let report = &response["manifest"]["face_measurements"];
    let face = &report["faces"][0];
    let lower =
        catallax_numbers::Rational::parse(face["deviation_mm"]["lower"].as_str().unwrap()).unwrap();
    let upper =
        catallax_numbers::Rational::parse(face["deviation_mm"]["upper"].as_str().unwrap()).unwrap();
    let truth = catallax_numbers::Rational::parse("0.0002").unwrap();
    assert!(lower <= truth && truth <= upper);
    assert_eq!(face["local_tolerance_mm"], face["deviation_mm"]["upper"]);
    assert_eq!(
        report["edit"]["resolved"][0]["faces"][0]["origins"],
        face["origins"]
    );
    assert_eq!(report["edit"]["resolved"][0]["status"], "resolved");
    let operations = response["operation_log"]["operations"].as_array().unwrap();
    let measured = operations
        .iter()
        .find(|o| o["operation"] == "face_measure")
        .unwrap();
    assert_eq!(measured["local_tolerances"][0], *face);
    let exported: catallax_cli::verification::ExportedQuestion =
        serde_json::from_slice(&fs::read(dir.join("question.json")).unwrap()).unwrap();
    let replay = answer(&exported.request(dir.join("replay")).unwrap()).unwrap();
    assert_eq!(
        replay["manifest"]["face_measurements"]["faces"],
        report["faces"]
    );
    assert_eq!(replay["work"], response["work"]);
    assert_eq!(
        fs::read(dir.join("answer.step")).unwrap(),
        fs::read(dir.join("replay/answer.step")).unwrap()
    );
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn unsupported_field_feature_removes_previous_display_without_brep_fallback() {
    use catallax_cli::verification::{Request, answer};
    let dir = std::env::temp_dir().join(format!(
        "catallax-field-display-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&dir).unwrap();
    fs::write(dir.join("display.obj"), "stale display").unwrap();
    fs::write(dir.join("answer.step"), "stale STEP").unwrap();
    let response = answer(&Request {
        format: "catallax-request/1".into(),
        source: "field".into(),
        problem: serde_json::from_str(include_str!("../../../parts/P2-collar.json")).unwrap(),
        out: dir.clone(),
        max_work: 100_000_000,
    })
    .unwrap();
    assert_eq!(
        response["manifest"]["refusal"]["code"],
        "field.unsupported-operation"
    );
    assert!(!dir.join("display.obj").exists());
    assert!(!dir.join("answer.step").exists());
    assert!(
        response["operation_log"]["operations"]
            .as_array()
            .unwrap()
            .iter()
            .all(|operation| operation["kernel"] == "field")
    );
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn curved_compare_preserves_work_refusals_and_reports_finite_agreement() {
    use catallax_cli::verification::{Request, answer};
    let root = std::env::temp_dir().join(format!(
        "catallax-distance-budget-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let problem = json!({"format":"catallax-recipe/3","name":"distance-budget","steps":[
        {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[
            {"kind":"line","id":"bottom","to":["2/5","0"]},
            {"kind":"line","id":"right","to":["2/5","2/5"]},
            {"kind":"line","id":"top","to":["1/5","2/5"]},
            {"kind":"line","id":"left","to":["1/5","0"]}
        ]]},
        {"id":"feature","op":"revolve","profile":"profile","axis":{"point":["0","0","0"],"direction":["0","1","0"]},"angle":"270"}
    ],"result":"feature"});
    for limit in [0, 1_000_000, 30_000_000, 100_000_000] {
        let out = root.join(limit.to_string());
        fs::create_dir_all(&out).unwrap();
        fs::write(out.join("answer.step"), "stale STEP").unwrap();
        let response = answer(&Request {
            format: "catallax-request/1".into(),
            source: "compare".into(),
            problem: problem.clone(),
            out: out.clone(),
            max_work: limit,
        })
        .unwrap();
        assert!(
            response["work"]["used"].as_u64().unwrap() <= limit,
            "{response}"
        );
        if response["work"]["exhausted"] == true {
            assert_eq!(response["manifest"]["status"], "refused", "{response}");
            assert_eq!(
                response["manifest"]["refusal"]["code"], "numbers.operation-budget",
                "{response}"
            );
        } else {
            assert_eq!(
                response["manifest"]["comparison"]["verdict"], "agree",
                "{response}"
            );
            let upper = catallax_numbers::Rational::parse(
                response["manifest"]["comparison"]["sampled_distance_upper_mm"]
                    .as_str()
                    .unwrap(),
            )
            .unwrap();
            assert!(
                upper <= catallax_numbers::Rational::parse("0.01").unwrap(),
                "{response}"
            );
        }
        if limit == 100_000_000 {
            assert_eq!(response["manifest"]["status"], "answered", "{response}");
        }
        assert!(!out.join("answer.step").exists());
    }
    fs::remove_dir_all(root).unwrap();
}
