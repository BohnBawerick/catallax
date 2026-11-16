//! Real stream requests keep numerical STEP separate from exact answers.

use serde_json::{Value, json};
use std::{
    fs,
    io::Write,
    process::{Command, Stdio},
};

fn recipe() -> Value {
    let mut recipe: Value = serde_json::from_str(include_str!(
        "../../catallax-features/tests/fixtures/loft-spline-prism.json"
    ))
    .unwrap();
    let steps = recipe["steps"].as_array_mut().unwrap();
    steps.push(json!({"id":"box","op":"block","size":["2","4","1"],
        "pose":{"rotation":["1","0","0","0"],"translation":["1","1","1"]}}));
    steps.push(json!({"id":"cut","op":"boolean","kind":"subtract","target":"loft","tool":"box"}));
    recipe["result"] = json!("cut");
    recipe
}

#[test]
fn numerical_step_requests_keep_tolerance_evidence_and_refused_artifacts_separate() {
    let dir = std::env::temp_dir().join(format!("catallax-numerical-step-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let problem = json!({"format":"catallax-numerical-step/1","recipe":recipe(),
        "tolerance_mm":"1/10000000"});
    let mut extra = problem.clone();
    extra["precision"] = json!(17);
    let mut invalid: Vec<_> = ["0", "-1", "garbage"]
        .into_iter()
        .map(|tolerance| {
            let mut value = problem.clone();
            value["tolerance_mm"] = json!(tolerance);
            ("recipe", value, 20_000_000)
        })
        .collect();
    invalid.extend([
        ("recipe", extra, 20_000_000),
        ("compare", problem.clone(), 20_000_000),
        ("recipe", problem.clone(), 0),
        (
            "recipe",
            json!({"format":"catallax-step/1","recipe":recipe()}),
            20_000_000,
        ),
        ("recipe", problem.clone(), 20_000_000),
    ]);
    let requests: Vec<_> = std::iter::once(("recipe", problem, 20_000_000))
        .chain(invalid)
        .collect();
    let mut child = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
        .arg("stream")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let mut input = child.stdin.take().unwrap();
    let input_dir = dir.clone();
    let writer = std::thread::spawn(move || {
        for (i, (source, problem, limit)) in requests.into_iter().enumerate() {
            let out = input_dir.join(i.to_string());
            fs::create_dir_all(&out).unwrap();
            fs::write(out.join("answer.step"), b"STALE").unwrap();
            writeln!(
                input,
                "{}",
                json!({"format":"catallax-request/1","source":source,
            "problem":problem,"max_work":limit,"out":out})
            )
            .unwrap();
        }
        drop(input);
    });
    let output = child.wait_with_output().unwrap();
    writer.join().unwrap();
    assert!(output.status.success());
    let rows: Vec<Value> = String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(rows.len(), 9);
    for i in [0, 8] {
        let row = &rows[i];
        let manifest = &row["manifest"];
        assert_eq!(manifest["format"], "catallax-numerical-manifest/1");
        assert_eq!(manifest["status"], "answered");
        assert_eq!(manifest["proof"], "numerical");
        assert_eq!(manifest["exact"], false);
        assert_eq!(manifest["category"], "approximate");
        assert_eq!(manifest["result"]["faces"], 9);
        for scope in ["source", "emitted"] {
            assert_eq!(manifest["numerical"][scope]["tolerance_mm"], "1/10000000");
            assert!(
                !manifest["numerical"][scope]["checks"]
                    .as_array()
                    .unwrap()
                    .is_empty()
            );
        }
        assert_eq!(manifest["step"]["complete_body_bound"], "unestablished");
        assert_eq!(manifest["step"]["seams"], "required");
        let work = row["work"]["used"].as_u64().unwrap();
        assert_eq!(
            manifest["work_stages"]
                .as_object()
                .unwrap()
                .values()
                .map(|v| v.as_u64().unwrap())
                .sum::<u64>(),
            work
        );
        assert!(work < 20_000_000);
        let file = fs::read(dir.join(i.to_string()).join("answer.step")).unwrap();
        assert!(
            file.windows(b"numerical approximate".len())
                .any(|s| s == b"numerical approximate")
        );
        assert!(catallax_step::check_seams(&file).is_ok());
        let saved: Value = serde_json::from_slice(
            &fs::read(dir.join(i.to_string()).join("manifest.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(saved, *manifest);
    }
    for (i, row) in rows.iter().enumerate().take(8).skip(1) {
        assert!(!dir.join(i.to_string()).join("answer.step").exists());
        assert!(row["manifest"]["step"].is_null());
        assert!(row["manifest"]["numerical"].is_null());
        assert_eq!(row["manifest"]["exact"], false);
        if i <= 5 {
            assert_eq!(row["manifest"]["status"], "invalid-question");
            assert_eq!(
                row["manifest"]["refusal"]["code"],
                "cli.numerical-step-wrapper"
            );
            assert_eq!(row["work"]["used"], 0);
        }
    }
    assert_eq!(
        rows[6]["manifest"]["format"],
        "catallax-numerical-manifest/1"
    );
    assert_eq!(
        rows[6]["manifest"]["refusal"]["code"],
        "numbers.operation-budget"
    );
    assert_eq!(
        rows[7]["manifest"]["refusal"]["code"],
        "numbers.geometric-gap-unproven"
    );
    assert_eq!(rows[0]["manifest"], rows[8]["manifest"]);
    assert_eq!(rows[0]["work"], rows[8]["work"]);
    assert_eq!(
        fs::read(dir.join("0/answer.step")).unwrap(),
        fs::read(dir.join("8/answer.step")).unwrap()
    );
    let mut portable: Value =
        serde_json::from_slice(&fs::read(dir.join("0/question.json")).unwrap()).unwrap();
    portable["max_work"] = json!(rows[0]["work"]["used"].as_u64().unwrap() - 1);
    let question = dir.join("late-budget.json");
    fs::write(&question, serde_json::to_vec(&portable).unwrap()).unwrap();
    let late = dir.join("late-budget");
    let output = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
        .arg("replay")
        .arg(&question)
        .arg("--out")
        .arg(&late)
        .output()
        .unwrap();
    assert!(output.status.success());
    let row: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(row["manifest"]["format"], "catallax-numerical-manifest/1");
    assert_eq!(
        row["manifest"]["refusal"]["code"],
        "numbers.operation-budget"
    );
    assert_eq!(
        row["operation_log"]["operations"]
            .as_array()
            .unwrap()
            .last()
            .unwrap()["operation"],
        "numerical-step-reread"
    );
    assert!(!late.join("answer.step").exists());
    let broken = dir.join("broken-metadata");
    fs::create_dir_all(broken.join("operations.json")).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
        .arg("replay")
        .arg(dir.join("0/question.json"))
        .arg("--out")
        .arg(&broken)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(!broken.join("answer.step").exists());
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn skew_collar_numerical_request_exports_the_actual_hole_with_step_trace_spans() {
    let proof = std::env::var_os("CATALLAX_SKEW_CLI_PROOF");
    let dir = proof.as_ref().map_or_else(
        || std::env::temp_dir().join(format!("catallax-skew-hole-{}", std::process::id())),
        std::path::PathBuf::from,
    );
    fs::create_dir_all(&dir).unwrap();
    let mut recipe: Value =
        serde_json::from_str(include_str!("../../../parts/P2-collar.json")).unwrap();
    recipe["steps"].as_array_mut().unwrap().push(json!({
        "id":"hole","op":"hole","target":"collar","diameter":"2","depth":"through","tip":"flat",
        "at":[{"id":"radial-0","point":["21/2","12","0"],"direction":["1","-1","1"]}]
    }));
    recipe["result"] = json!("hole");
    let request = json!({"format":"catallax-request/1","source":"recipe",
        "problem":{"format":"catallax-numerical-step/1","recipe":recipe,"tolerance_mm":"1e-7"},
        "max_work":20_000_000,"out":dir});
    let mut child = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
        .arg("stream")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    writeln!(child.stdin.take().unwrap(), "{request}").unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    let row: Value = serde_json::from_slice(&output.stdout).unwrap();
    let manifest = &row["manifest"];
    assert_eq!(manifest["status"], "answered", "{row}");
    assert_eq!(manifest["proof"], "numerical");
    assert_eq!(manifest["exact"], false);
    assert_eq!(manifest["category"], "approximate");
    assert_eq!(manifest["result"]["faces"], 23);
    assert_eq!(manifest["step"]["tolerance_mm"], "1/10000000");
    assert!(row["work"]["used"].as_u64().unwrap() < 20_000_000);
    assert_eq!(
        manifest["work_stages"]
            .as_object()
            .unwrap()
            .values()
            .map(|v| v.as_u64().unwrap())
            .sum::<u64>(),
        row["work"]["used"].as_u64().unwrap()
    );
    let operations = row["operation_log"]["operations"].as_array().unwrap();
    for required in ["hole", "numerical-step", "numerical-step-reread"] {
        assert!(
            operations
                .iter()
                .any(|operation| operation["operation"] == required),
            "{required}"
        );
    }
    let file = fs::read(dir.join("answer.step")).unwrap();
    assert!(catallax_step::check_seams(&file).is_ok());
    fs::write(
        dir.join("test-response.json"),
        serde_json::to_vec_pretty(&row).unwrap(),
    )
    .unwrap();
    if proof.is_none() {
        fs::remove_dir_all(dir).unwrap();
    }
}
