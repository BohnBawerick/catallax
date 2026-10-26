//! Explicit IGES verification requests through the public stream and replay commands.
use serde_json::{Value, json};
use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "catallax-iges-{name}-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&dir).unwrap();
    dir
}
fn block(size: [&str; 3]) -> Value {
    json!({"format":"catallax-recipe/3","name":"box",
        "steps":[{"id":"box","op":"block","size":size,
        "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"box"})
}
fn bored() -> Value {
    serde_json::from_str(r#"{
      "format":"catallax-recipe/3","name":"bored plate",
      "steps":[{"id":"profile","op":"sketch",
      "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
      "loops":[[{"kind":"line","id":"a","to":["10","0"]},
                {"kind":"line","id":"b","to":["10","8"]},
                {"kind":"line","id":"c","to":["0","8"]},
                {"kind":"line","id":"d","to":["0","0"]}],
               [{"kind":"circle","id":"bore","centre":["5","4"],"radius":"2"}]]},
               {"id":"plate","op":"extrude","profile":"profile","distance":"2"}],"result":"plate"}"#)
    .unwrap()
}
fn collar() -> Value {
    serde_json::from_str(include_str!("../../../parts/P2-collar.json")).unwrap()
}
fn sphere() -> Value {
    json!({"format":"catallax-recipe/3","name":"ball",
        "steps":[{"id":"ball","op":"sphere","radius":"2",
        "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"ball"})
}
fn iges(recipe: Value) -> Value {
    json!({"format":"catallax-iges/1","recipe":recipe})
}
fn request(problem: Value, out: &Path, max_work: u64) -> Value {
    json!({"format":"catallax-request/1","source":"iges","problem":problem,
        "out":out,"max_work":max_work})
}
/// One worker process answering its requests in order.
fn stream(requests: &[Value]) -> Vec<Value> {
    let mut child = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
        .arg("stream")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let mut stdin = child.stdin.take().unwrap();
    for request in requests {
        writeln!(stdin, "{request}").unwrap();
    }
    drop(stdin);
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect()
}
fn answer(problem: Value, max_work: u64) -> (Value, PathBuf) {
    let dir = scratch("one");
    let response = stream(&[request(problem, &dir, max_work)]).remove(0);
    (response, dir)
}
fn is_igs(dir: &Path) -> bool {
    dir.join("answer.igs").exists()
}
fn strip_timing(mut response: Value) -> Value {
    response.as_object_mut().unwrap().remove("wall_ms");
    for operation in response["operation_log"]["operations"]
        .as_array_mut()
        .unwrap()
    {
        operation.as_object_mut().unwrap().remove("wall_ms");
    }
    response
}

#[test]
fn supported_parts_publish_the_actual_report_and_artifact_identity() {
    // Faces, solids and the bore's cylinder pin the topology the writer reports.
    for (name, recipe, faces) in [("box", block(["2", "3", "5"]), 6), ("bored", bored(), 7)] {
        let (response, dir) = answer(iges(recipe), 100_000_000);
        let manifest = &response["manifest"];
        assert_eq!(manifest["status"], "answered", "{name}: {manifest}");
        assert!(manifest["result"].is_null() && manifest["properties"].is_null());
        assert!(manifest["step"].is_null());
        let report = &manifest["iges"];
        assert_eq!(report["file"], "answer.igs");
        assert_eq!(report["units"], "mm");
        assert_eq!(report["solids"], 1);
        assert_eq!(report["faces"], faces, "{name}");
        let bytes = fs::read(dir.join("answer.igs")).unwrap();
        assert_eq!(report["bytes"], bytes.len());
        assert_eq!(report["sha256"].as_str().unwrap().len(), 64);
        assert!(String::from_utf8(bytes).unwrap().contains("S      1G"));
        assert!(!dir.join("answer.step").exists() && !dir.join("display.obj").exists());
        assert_eq!(response["work"]["exhausted"], false);
    }
}

#[test]
fn the_collar_exports_and_the_manifest_matches_the_library_report() {
    let (response, dir) = answer(iges(collar()), 100_000_000);
    let manifest = &response["manifest"];
    assert_eq!(manifest["status"], "answered", "{manifest}");
    let recipe = catallax_recipe::Recipe::parse(&collar().to_string()).unwrap();
    let mut dec = catallax_numbers::Decider::default();
    let built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
    let mut direct = Vec::new();
    let report = catallax_iges::write_iges(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut dec,
        &mut direct,
    )
    .unwrap();
    assert_eq!(fs::read(dir.join("answer.igs")).unwrap(), direct);
    let published = &manifest["iges"];
    assert_eq!(published["faces"], report.faces);
    assert_eq!(published["edges"], report.edges);
    assert_eq!(published["vertices"], report.vertices);
    assert_eq!(published["added_edges"], report.added_edges);
    assert_eq!(published["added_vertices"], report.added_vertices);
    assert_eq!(
        published["max_deviation_mm"],
        report.max_deviation.to_string()
    );
    assert_eq!(published["scale"], report.scale.to_string());
    assert_eq!(
        published["coordinate_bound_mm"],
        report.coordinate_bound.to_string()
    );
}

#[test]
fn refusals_are_typed_and_remove_the_previous_artifact() {
    let dir = scratch("transition");
    let ok = stream(&[request(iges(block(["2", "3", "5"])), &dir, 100_000_000)]).remove(0);
    assert_eq!(ok["manifest"]["status"], "answered");
    assert!(is_igs(&dir));
    let tight = json!({"format":"catallax-iges/1","recipe":block(["1/3","2/7","5/11"]),
        "tolerance":"1/1000000000000000000000000000000000000"});
    let cases = [
        (iges(sphere()), 100_000_000, "iges.unsupported-surface"),
        (tight, 100_000_000, "iges.tolerance-exceeded"),
        (
            json!({"format":"catallax-iges/1","recipe":block(["2","3","5"]),"scale_digits":0}),
            100_000_000,
            "iges.invalid-options",
        ),
        (
            json!({"format":"catallax-iges/1","recipe":block(["2","3","5"]),"tolerance":"-1"}),
            100_000_000,
            "iges.invalid-options",
        ),
        (iges(collar()), 0, "numbers.operation-budget"),
        (iges(bored()), 5_000, "numbers.operation-budget"),
    ];
    for (problem, work, code) in cases {
        // Each case follows a success, so a surviving file would be stale.
        let ok = stream(&[request(iges(block(["2", "3", "5"])), &dir, 100_000_000)]).remove(0);
        assert_eq!(ok["manifest"]["status"], "answered");
        assert!(is_igs(&dir));
        let response = stream(&[request(problem, &dir, work)]).remove(0);
        let manifest = &response["manifest"];
        assert_eq!(manifest["status"], "refused", "{code}: {manifest}");
        assert_eq!(manifest["refusal"]["code"], code, "{manifest}");
        assert!(manifest["iges"].is_null());
        assert!(!is_igs(&dir), "{code} left a stale artifact");
    }
}

#[test]
fn malformed_requests_are_invalid_questions() {
    for problem in [
        json!({"format":"catallax-iges/2","recipe":block(["1","1","1"])}),
        json!({"format":"catallax-iges/1","recipe":block(["1","1","1"]),"extra":1}),
        json!({"format":"catallax-iges/1","recipe":block(["1","1","1"]),"units":"in"}),
        json!({"format":"catallax-iges/1","recipe":block(["1","1","1"]),"scale_digits":-1}),
        json!({"format":"catallax-iges/1","recipe":{"format":"nope"}}),
    ] {
        let (response, dir) = answer(problem, 100_000_000);
        assert_eq!(response["manifest"]["status"], "invalid-question");
        assert!(!is_igs(&dir));
    }
}

#[test]
fn a_process_io_failure_is_not_a_geometric_refusal() {
    // A directory where the previous artifact should be makes cleanup fail.
    let dir = scratch("io");
    fs::create_dir_all(dir.join("answer.igs")).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
        .arg("stream")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    writeln!(
        child.stdin.take().unwrap(),
        "{}",
        request(iges(block(["2", "3", "5"])), &dir, 100_000_000)
    )
    .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty(), "no response claims an outcome");
    assert!(dir.join("answer.igs").is_dir());
}

#[test]
fn replay_and_worker_counts_reproduce_outcomes_work_and_bytes() {
    let problems = [
        iges(block(["2", "3", "5"])),
        iges(bored()),
        iges(collar()),
        iges(sphere()),
    ];
    let baseline_dir = scratch("base");
    let requests: Vec<_> = problems
        .iter()
        .enumerate()
        .map(|(i, p)| request(p.clone(), &baseline_dir.join(i.to_string()), 100_000_000))
        .collect();
    let baseline: Vec<_> = stream(&requests).into_iter().map(strip_timing).collect();
    for workers in [1usize, 2, 4] {
        let results: Vec<_> = std::thread::scope(|scope| {
            (0..workers)
                .map(|w| {
                    let root = scratch(&format!("w{workers}-{w}"));
                    let requests: Vec<_> = problems
                        .iter()
                        .enumerate()
                        .map(|(i, p)| request(p.clone(), &root.join(i.to_string()), 100_000_000))
                        .collect();
                    scope.spawn(move || (root, stream(&requests)))
                })
                .collect::<Vec<_>>()
                .into_iter()
                .map(|handle| handle.join().unwrap())
                .collect()
        });
        for (root, responses) in results {
            for (i, response) in responses.into_iter().enumerate() {
                assert_eq!(
                    strip_timing(response),
                    baseline[i],
                    "{workers} workers, {i}"
                );
                let a = baseline_dir.join(i.to_string()).join("answer.igs");
                let b = root.join(i.to_string()).join("answer.igs");
                assert_eq!(a.exists(), b.exists());
                if a.exists() {
                    assert_eq!(fs::read(a).unwrap(), fs::read(b).unwrap());
                }
            }
        }
    }
    // The portable question replays the same outcome, work and bytes.
    for (i, expected) in baseline.iter().enumerate() {
        let original = baseline_dir.join(i.to_string());
        let replay = original.join("replay");
        let output = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
            .args([
                "replay",
                original.join("question.json").to_str().unwrap(),
                "--out",
                replay.to_str().unwrap(),
            ])
            .output()
            .unwrap();
        assert!(output.status.success());
        let replayed = strip_timing(serde_json::from_slice(&output.stdout).unwrap());
        assert_eq!(replayed["manifest"], expected["manifest"]);
        assert_eq!(replayed["work"], expected["work"]);
        assert_eq!(replayed["operation_log"], expected["operation_log"]);
        assert_eq!(
            original.join("answer.igs").exists(),
            replay.join("answer.igs").exists()
        );
        if replay.join("answer.igs").exists() {
            assert_eq!(
                fs::read(original.join("answer.igs")).unwrap(),
                fs::read(replay.join("answer.igs")).unwrap()
            );
        }
    }
}

#[test]
fn existing_recipe_requests_keep_their_artifacts_and_gain_no_iges_file() {
    let dir = scratch("recipe");
    let response = stream(&[json!({"format":"catallax-request/1","source":"recipe",
        "problem":block(["2","3","5"]),"out":dir,"max_work":100_000_000})])
    .remove(0);
    assert_eq!(response["manifest"]["status"], "answered");
    assert!(dir.join("answer.step").exists() && !is_igs(&dir));
    assert!(response["manifest"]["iges"].is_null());
}

#[test]
fn about_advertises_the_source_and_problem_format() {
    let output = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
        .arg("about")
        .output()
        .unwrap();
    let about: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(
        about["sources"]
            .as_array()
            .unwrap()
            .contains(&json!("iges"))
    );
    assert_eq!(about["iges"], "catallax-iges/1");
}
