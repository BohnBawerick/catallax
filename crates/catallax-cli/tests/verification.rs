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
    let directory = root.join("recipe");
    let original = answer(&Request {
        format: "catallax-request/1".into(),
        source: "recipe".into(),
        problem: recipe,
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
    assert_eq!(original["manifest"]["status"], "answered");
    assert_eq!(
        fs::read(directory.join("answer.step")).unwrap(),
        fs::read(directory.join("replay/answer.step")).unwrap()
    );
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

fn certified_problem() -> Value {
    serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../.claude/skills/verify-catallax/references/examples/certified-recipe.json"
    )))
    .unwrap()
}

#[test]
fn certified_output_refusal_preserves_the_complete_valid_body() {
    use catallax_cli::verification::{Request, answer};
    let dir = std::env::temp_dir().join(format!(
        "catallax-certified-output-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let mut problem = certified_problem();
    problem["outputs"] = json!({"step":{"file_tolerance":"1e-12","decimal_places":0}});
    problem["properties"] = Value::Null;
    problem["field_points"] = json!([]);
    problem["edit"] = Value::Null;
    let response = answer(&Request {
        format: "catallax-request/1".into(),
        source: "recipe".into(),
        problem,
        out: dir.clone(),
        max_work: 100_000_000,
    })
    .unwrap();
    let manifest = &response["manifest"];
    assert_eq!(manifest["status"], "answered", "{response}");
    assert_eq!(manifest["result"]["kind"], "solid");
    assert_eq!(manifest["result"]["faces_by_kind"]["bspline"], 1);
    assert_eq!(manifest["result"]["ring_edges"], 2);
    assert_eq!(manifest["step"]["status"], "refused");
    assert_eq!(
        manifest["step"]["refusal"]["code"],
        "step.tolerance-not-reached"
    );
    assert!(!dir.join("answer.step").exists());
    assert!(!response["work"]["exhausted"].as_bool().unwrap());
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn certified_budget_refusal_retains_its_version_and_publishes_no_partial_body() {
    use catallax_cli::verification::{Request, answer};
    let dir = std::env::temp_dir().join(format!(
        "catallax-certified-budget-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let response = answer(&Request {
        format: "catallax-request/1".into(),
        source: "recipe".into(),
        problem: certified_problem(),
        out: dir.clone(),
        max_work: 0,
    })
    .unwrap();
    let manifest = &response["manifest"];
    assert_eq!(manifest["format"], "catallax-manifest/4");
    assert_eq!(manifest["status"], "refused");
    assert_eq!(manifest["refusal"]["code"], "numbers.operation-budget");
    assert!(manifest["result"].is_null() && manifest["certificate"].is_null());
    assert!(!dir.join("answer.step").exists());
    assert!(!dir.join("display.obj").exists());
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

#[test]
fn comparison_query_exhaustion_is_a_budget_refusal_not_an_internal_error() {
    use catallax_cli::verification::{Request, answer};
    let dir = std::env::temp_dir().join(format!(
        "catallax-compare-query-budget-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let mut problem: Value =
        serde_json::from_str(include_str!("../../../parts/P5-flanged-housing.json")).unwrap();
    problem["steps"].as_array_mut().unwrap().truncate(2);
    problem["steps"].as_array_mut().unwrap().push(json!({
        "id":"drilled","op":"hole","target":"housing",
        "at":[{"id":"bolt","point":["17","5","0"],"direction":["0","-1","0"]}],
        "diameter":"4","depth":"through","tip":"flat",
        "head":{"countersink":{"diameter":"8","angle":"90"}}
    }));
    problem["result"] = json!("drilled");
    let response = answer(&Request {
        format: "catallax-request/1".into(),
        source: "compare".into(),
        problem,
        out: dir.clone(),
        max_work: 100_000_000,
    })
    .unwrap();
    if response["work"]["exhausted"] == true {
        assert_eq!(response["manifest"]["status"], "refused", "{response}");
        assert_eq!(
            response["manifest"]["refusal"]["code"],
            "numbers.operation-budget"
        );
        assert!(!dir.join("answer.step").exists());
    } else {
        assert_eq!(response["manifest"]["status"], "answered", "{response}");
        assert!(matches!(
            response["manifest"]["comparison"]["verdict"].as_str(),
            Some("agree" | "unresolved")
        ));
    }
    fs::remove_dir_all(dir).unwrap();
}

/// Standard algebraic drill angles retain exact construction intent within the
/// existing whole-request allowance, including summary and STEP export.
#[test]
fn standard_hole_angles_build_at_the_existing_request_allowance() {
    use catallax_cli::verification::{Request, answer};
    for (tip, head) in [
        (json!({"cone":"118"}), None),
        (
            json!("flat"),
            Some(json!({"countersink":{"diameter":"8","angle":"82"}})),
        ),
        (
            json!("flat"),
            Some(json!({"countersink":{"diameter":"8","angle":"100"}})),
        ),
    ] {
        let dir = std::env::temp_dir().join(format!(
            "catallax-standard-hole-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let mut hole = json!({"id":"drill","op":"hole","target":"stock","at":[{"id":"centre","point":["0","0","0"],"direction":["0","0","1"]}],"diameter":"2","depth":{"blind":"4"},"tip":tip});
        if let Some(head) = head {
            hole["head"] = head;
        }
        let problem = json!({"format":"catallax-recipe/3","name":"standard hole angle coverage","steps":[
            {"id":"stock","op":"block","size":["10","10","10"],"pose":{"rotation":["1","0","0","0"],"translation":["-5","-5","0"]}},hole
        ],"result":"drill"});
        let response = answer(&Request {
            format: "catallax-request/1".into(),
            source: "recipe".into(),
            problem,
            out: dir.clone(),
            max_work: 100_000_000,
        })
        .unwrap();
        assert_eq!(response["manifest"]["status"], "answered", "{response}");
        assert_eq!(response["manifest"]["exact"], true, "{response}");
        assert_eq!(response["manifest"]["marks"], json!([]));
        assert_eq!(response["work"]["exhausted"], false, "{response}");
        assert!(dir.join("answer.step").is_file());
        assert!(fs::metadata(dir.join("answer.step")).unwrap().len() > 0);
        fs::remove_dir_all(dir).unwrap();
    }
}

#[test]
fn planar_drafts_export_checked_step_at_the_application_allowance() {
    use catallax_cli::verification::{Request, answer};
    use catallax_measure::{Accuracy, face_centroid};
    use catallax_numbers::{Decider, Rational};
    use catallax_recipe::{ElementId, Origin, Recipe, Role, StepId};
    let side = |segment: &str| Role::Side {
        segment: ElementId::new(segment).unwrap(),
        path: None,
    };
    for (angle, fixture) in [
        (3.0_f64, include_str!("fixtures/planar-draft/draft-3.json")),
        (5.0, include_str!("fixtures/planar-draft/draft-5.json")),
    ] {
        let problem: Value = serde_json::from_str(fixture).unwrap();
        let dir = std::env::temp_dir().join(format!(
            "catallax-planar-draft-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let request = |out: std::path::PathBuf, max_work| Request {
            format: "catallax-request/1".into(),
            source: "recipe".into(),
            problem: problem.clone(),
            out,
            max_work,
        };
        let response = answer(&request(dir.join("application"), 20_000_000)).unwrap();
        let manifest = &response["manifest"];
        assert_eq!(manifest["status"], "answered", "{response}");
        assert_eq!(manifest["format"], "catallax-manifest/3");
        assert_eq!(manifest["exact"], true);
        assert_eq!(manifest["marks"], json!([]));
        assert_eq!(response["work"]["exhausted"], false, "{response}");
        assert_eq!(response["work"]["limit"], 20_000_000);
        assert_eq!(manifest["result"]["kind"], "solid");
        assert_eq!(manifest["result"]["lumps"], 1);
        assert_eq!(manifest["result"]["shells"], 1);
        assert_eq!(manifest["result"]["faces_by_kind"], json!({"plane":6}));
        assert_eq!(manifest["result"]["edges_by_kind"], json!({"line":12}));
        assert_eq!(manifest["result"]["vertices"], 8);
        assert_eq!(manifest["step"]["approximated_curves"], 0);
        let deviation =
            Rational::parse(manifest["step"]["max_deviation"].as_str().unwrap()).unwrap();
        assert!(!deviation.is_zero() && deviation <= Rational::from_ratio(1, 10_000_000).unwrap());
        let bytes = fs::read(dir.join("application/answer.step")).unwrap();
        let text = String::from_utf8_lossy(&bytes);
        assert_eq!(text.matches("=ADVANCED_FACE(").count(), 6);
        assert_eq!(text.matches("=PLANE(").count(), 6);
        assert!(!text.contains("B_SPLINE"));

        let replay = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
            .args([
                "replay",
                dir.join("application/question.json").to_str().unwrap(),
                "--out",
                dir.join("fresh").to_str().unwrap(),
            ])
            .output()
            .unwrap();
        assert!(replay.status.success());
        let replay: Value = serde_json::from_slice(&replay.stdout).unwrap();
        assert_eq!(replay["work"]["limit"], 20_000_000);
        assert_eq!(replay["manifest"], response["manifest"]);
        assert_eq!(fs::read(dir.join("fresh/answer.step")).unwrap(), bytes);
        let wider = answer(&request(dir.join("replay"), 100_000_000)).unwrap();
        assert_eq!(wider["manifest"], response["manifest"]);
        assert_eq!(fs::read(dir.join("replay/answer.step")).unwrap(), bytes);

        // The manifest has no origins. Each one must name the right support.
        let recipe = Recipe::parse(fixture).unwrap();
        let mut dec = Decider::default();
        let built =
            catallax_app::build_incremental_without_display(&recipe, None, &mut dec).unwrap();
        let t = angle.to_radians().tan();
        let expected = [
            (Role::Start, [0.0, 0.0, 0.0]),
            (Role::End, [0.0, 0.0, 20.0]),
            (side("right"), [20.0 - 10.0 * t, 0.0, 10.0]),
            (side("left"), [10.0 * t - 20.0, 0.0, 10.0]),
            (
                side("top"),
                [0.0, 15.0, (10.0 - 40.0 * t / 6.0) / (1.0 - t / 2.0)],
            ),
            (
                side("bottom"),
                [0.0, -15.0, (10.0 - 40.0 * t / 6.0) / (1.0 - t / 2.0)],
            ),
        ];
        assert_eq!(built.names.faces.len(), expected.len());
        for (role, centroid) in expected {
            let origin = [Origin {
                step: StepId::new("base").unwrap(),
                role,
            }];
            let (face, _) = built
                .names
                .faces
                .iter()
                .find(|(_, origins)| origins[..] == origin)
                .unwrap_or_else(|| panic!("no face for {origin:?}"));
            let model = &built.checkpoint.built().model;
            let found = face_centroid(model, *face, &Accuracy::default(), &mut dec).unwrap();
            for (bounds, expected) in found.value.iter().zip(centroid) {
                assert!(
                    (bounds.lo.to_f64_nearest() - expected).abs() < 1e-9,
                    "{origin:?}"
                );
            }
        }
        fs::remove_dir_all(dir).unwrap();
    }
}

#[test]
fn placed_preserved_bore_draft_answers_at_the_application_allowance() {
    use catallax_cli::verification::{Request, answer};
    let problem: Value = serde_json::from_str(include_str!(
        "../../catallax-features/tests/fixtures/draft-preserved-bore/K-R.json"
    ))
    .unwrap();
    let dir = std::env::temp_dir().join(format!(
        "catallax-preserved-bore-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let response = answer(&Request {
        format: "catallax-request/1".into(),
        source: "recipe".into(),
        problem,
        out: dir.clone(),
        max_work: 20_000_000,
    })
    .unwrap();
    let manifest = &response["manifest"];
    assert_eq!(manifest["status"], "answered", "{response}");
    assert_eq!(manifest["format"], "catallax-manifest/3");
    assert_eq!(manifest["exact"], true);
    assert_eq!(manifest["marks"], json!([]));
    assert_eq!(response["work"]["exhausted"], false, "{response}");
    assert_eq!(
        manifest["result"]["faces_by_kind"],
        json!({"plane":6,"cylinder":1})
    );
    assert_eq!(
        manifest["result"]["edges_by_kind"],
        json!({"line":12,"circle":2})
    );
    assert_eq!(manifest["step"]["approximated_curves"], 0);
    assert!(!fs::read(dir.join("answer.step")).unwrap().is_empty());
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn captured_slot_cut_answers_with_the_default_work_allowance() {
    use catallax_cli::verification::{Request, answer};
    let captured: Value =
        serde_json::from_str(include_str!("fixtures/soliddraft-slot-cut.json")).unwrap();
    let directory = std::env::temp_dir().join(format!(
        "catallax-slot-cut-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let response = answer(&Request {
        format: "catallax-request/1".into(),
        source: "recipe".into(),
        problem: captured["problem"].clone(),
        out: directory.clone(),
        max_work: 100_000_000,
    })
    .unwrap();
    assert_eq!(
        response["manifest"]["status"], "answered",
        "{}",
        response["manifest"]
    );
    assert_eq!(response["manifest"]["result"]["kind"], "solid");
    assert!(directory.join("answer.step").is_file());
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn refused_print_request_removes_the_previous_print_and_display_artifacts() {
    use catallax_cli::verification::{Request, answer};
    let dir = std::env::temp_dir().join(format!(
        "catallax-print-transition-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let recipe: Value = serde_json::from_str(include_str!(
        "../../../.claude/skills/verify-catallax/references/examples/recipe.json"
    ))
    .unwrap();
    let mut request = Request {
        format: "catallax-request/1".into(),
        source: "print".into(),
        problem: json!({"format":"catallax-print/1","recipe":recipe,"tolerance":"1"}),
        out: dir.clone(),
        max_work: 100_000_000,
    };
    let first = answer(&request).unwrap();
    assert_eq!(first["manifest"]["status"], "answered", "{first}");
    assert!(dir.join("print.stl").exists());
    assert!(dir.join("display.obj").exists());
    assert!(!dir.join("answer.step").exists());
    request.max_work = 0;
    let refused = answer(&request).unwrap();
    assert_eq!(
        refused["manifest"]["refusal"]["code"],
        "numbers.operation-budget"
    );
    assert!(!dir.join("print.stl").exists());
    assert!(!dir.join("display.obj").exists());
    assert!(!dir.join("answer.step").exists());
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn the_field_adapter_shows_a_lattice_as_an_unproved_preview() {
    use catallax_cli::verification::{Request, answer};
    let root = std::env::temp_dir().join(format!(
        "catallax-lattice-preview-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let pose = json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]});
    let recipe = json!({"format":"catallax-recipe/3","name":"preview",
        "steps":[
            {"id":"region","op":"block","size":["4","4","4"],"pose":pose},
            {"id":"beams","op":"lattice","region":"region","cell":"cubic","size":"2",
             "frame":pose,"radius":{"constant":"0.3"}}],
        "result":"beams"});
    let request = |problem: Value, name: &str| Request {
        format: "catallax-request/1".into(),
        source: "field".into(),
        problem,
        out: root.join(name),
        max_work: 100_000_000,
    };
    let lattice = answer(&request(recipe, "lattice")).unwrap();
    let manifest = &lattice["manifest"];
    assert_eq!(manifest["status"], "answered");
    assert_eq!(manifest["exact"], false);
    let display = &manifest["face_measurements"]["display"];
    assert_eq!(display["proof"], "unproved");
    assert!(display["triangles"].as_u64().unwrap() > 0);
    assert!(root.join("lattice/display.obj").exists());
    // The certified request is untouched: a plain block is still exact and carries no proof label.
    let block = json!({"format":"catallax-recipe/3","name":"block",
        "steps":[{"id":"region","op":"block","size":["4","4","4"],"pose":pose}],"result":"region"});
    let plain = answer(&request(block, "block")).unwrap();
    assert_eq!(plain["manifest"]["exact"], true);
    assert!(plain["manifest"]["face_measurements"]["display"]["proof"].is_null());
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn a_fitted_recipe_answers_with_its_certificate_through_both_routes() {
    let recipe: Value =
        serde_json::from_str(include_str!("../../../parts/P5-flanged-housing.json")).unwrap();
    let wrapped = json!({"format":"catallax-certified-recipe/1","recipe":recipe,"outputs":{}});
    let dir = std::env::temp_dir().join(format!(
        "catallax-fitted-{}-{}",
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
    for (i, problem) in [recipe, wrapped].iter().enumerate() {
        writeln!(
            input,
            "{}",
            json!({"format":"catallax-request/1","source":"recipe","problem":problem,
                "max_work":20_000_000,"out":dir.join(i.to_string())})
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
    assert_eq!(rows.len(), 2);
    for row in &rows {
        let manifest = &row["manifest"];
        assert_eq!(manifest["format"], "catallax-manifest/4");
        assert_eq!(manifest["status"], "answered");
        assert_eq!(manifest["exact"], false);
        let regions = manifest["certificate"]["regions"].as_array().unwrap();
        assert_eq!(regions.len(), 22);
        assert_eq!(
            regions
                .iter()
                .filter(|region| !region["exception"].is_null())
                .count(),
            4
        );
    }
    assert_eq!(rows[0]["manifest"]["step"]["status"], "answered");
    assert!(dir.join("0/answer.step").exists());
    assert!(!dir.join("1/answer.step").exists());
    assert_eq!(
        rows[0]["manifest"]["certificate"]["identity"]["target"],
        rows[1]["manifest"]["certificate"]["identity"]["target"]
    );
    fs::remove_dir_all(&dir).unwrap();
}
