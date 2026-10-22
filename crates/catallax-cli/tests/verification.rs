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
