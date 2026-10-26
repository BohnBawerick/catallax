//! End-to-end checks of the kernel-side exchange, without judge dependencies.

use std::{
    fs,
    path::PathBuf,
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
};

use catallax_numbers::Rational;
use serde_json::{Value, json};

static NEXT: AtomicU64 = AtomicU64::new(0);

struct Sandbox(PathBuf);

impl Sandbox {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "exact-brep-cli-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }

    fn run(&self, args: &[&str]) -> Output {
        Command::new(env!("CARGO_BIN_EXE_exact-brep"))
            .args(args)
            .env_clear()
            .current_dir(&self.0)
            .output()
            .unwrap()
    }

    fn answer_text(&self, text: &str) -> (Value, String) {
        fs::write(self.0.join("question.json"), text).unwrap();
        let output = self.run(&["answer", "question.json", "--out", "out"]);
        assert_eq!(output.status.code(), Some(0), "{output:?}");
        assert!(output.stdout.is_empty(), "{output:?}");
        assert!(output.stderr.is_empty(), "{output:?}");
        let text = fs::read_to_string(self.0.join("out/manifest.json")).unwrap();
        let manifest: Value = serde_json::from_str(&text).unwrap();
        envelope(&manifest);
        (manifest, text)
    }

    fn answer(&self, question: &Value) -> Value {
        self.answer_text(&question.to_string()).0
    }
}

impl Drop for Sandbox {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}

fn rational(value: &Value) -> Rational {
    let text = value.as_str().unwrap();
    assert!(
        !text.contains('/'),
        "manifest number is not decimal: {text}"
    );
    Rational::parse(text).unwrap()
}

fn enclosure(value: &Value, expected_lo: &str, expected_hi: &str, relative: &str, absolute: &str) {
    keys(value, &["lo", "hi"]);
    let lo = rational(&value["lo"]);
    let hi = rational(&value["hi"]);
    assert!(lo <= Rational::parse(expected_lo).unwrap(), "{value}");
    assert!(hi >= Rational::parse(expected_hi).unwrap(), "{value}");
    let magnitude = if lo.sign() == hi.sign() {
        lo.abs().min(hi.abs())
    } else {
        Rational::zero()
    };
    let allowed =
        (Rational::parse(relative).unwrap() * magnitude).max(Rational::parse(absolute).unwrap());
    assert!(hi - lo <= allowed, "enclosure too wide: {value}");
}

fn keys(value: &Value, expected: &[&str]) {
    let mut actual: Vec<_> = value
        .as_object()
        .unwrap()
        .keys()
        .map(String::as_str)
        .collect();
    actual.sort_unstable();
    let mut expected = expected.to_vec();
    expected.sort_unstable();
    assert_eq!(actual, expected);
}

fn envelope(manifest: &Value) {
    keys(
        manifest,
        &[
            "format",
            "question_id",
            "kernel",
            "status",
            "exact",
            "result",
            "properties",
            "step",
            "marks",
            "refusal",
            "stats",
        ],
    );
    assert_eq!(manifest["format"], "exact-brep-manifest/2");
    keys(&manifest["kernel"], &["name", "version", "commit"]);
    assert_eq!(manifest["kernel"]["name"], "exact-brep");
    assert_eq!(manifest["kernel"]["version"], env!("CARGO_PKG_VERSION"));
    let commit = manifest["kernel"]["commit"].as_str().unwrap();
    assert_eq!(commit.len(), 40);
    assert!(commit.bytes().all(|c| c.is_ascii_hexdigit()));
    let stats = &manifest["stats"];
    keys(
        stats,
        &[
            "wall_ms",
            "decisions",
            "filter",
            "refined",
            "exact",
            "coincident_exact",
            "coincident_separation",
            "coincident_fact",
            "not_exact",
            "refused",
            "max_bits",
            "work",
            "per_step",
        ],
    );
    for (key, value) in stats.as_object().unwrap() {
        if key != "per_step" {
            assert!(value.is_u64(), "{key}: {value}");
        }
    }
    let outcomes: u64 = [
        "filter",
        "refined",
        "exact",
        "coincident_exact",
        "coincident_separation",
        "coincident_fact",
        "not_exact",
        "refused",
    ]
    .iter()
    .map(|name| stats[name].as_u64().unwrap())
    .sum();
    assert_eq!(stats["decisions"].as_u64().unwrap(), outcomes);
    for step in stats["per_step"].as_array().unwrap() {
        keys(step, &["id", "wall_ms", "decisions", "not_exact"]);
        assert!(step["id"].is_string());
        for key in ["wall_ms", "decisions", "not_exact"] {
            assert!(step[key].is_u64());
        }
    }
    assert_eq!(
        manifest["exact"],
        manifest["status"] == "answered" && manifest["marks"] == json!([])
    );
    for mark in manifest["marks"].as_array().unwrap() {
        keys(mark, &["code", "kind", "bound", "step", "detail"]);
        assert!(["distance", "angle", "value"].contains(&mark["kind"].as_str().unwrap()));
        assert!(rational(&mark["bound"]) >= Rational::zero());
        assert!(mark["step"].is_null() || mark["step"].is_string());
        assert!(mark["detail"].is_string());
    }
    if manifest["status"] != "answered" {
        assert!(manifest["result"].is_null());
        assert!(manifest["properties"].is_null());
        assert!(manifest["step"].is_null());
        keys(&manifest["refusal"], &["code", "step", "stage", "detail"]);
        assert!(manifest["refusal"]["detail"].is_string());
    } else {
        assert!(manifest["refusal"].is_null());
    }
}

fn question(shape: Value) -> Value {
    json!({
        "format": "exact-brep-question/2", "id": "opaque-id",
        "operands": {"a": shape}, "steps": [], "result": "a"
    })
}

fn complete_example() -> Value {
    json!({
        "format": "exact-brep-question/2", "id": "q-7f3a9c01",
        "operands": {
            "a": {"shape": "block", "size": ["10", "10", "10"]},
            "b": {"shape": "cylinder", "radius": "3", "height": "12",
                "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["5", "5", "-1"]}}
        },
        "steps": [{"id": "s1", "op": "subtract", "args": ["a", "b"]}],
        "result": "s1",
        "measure": {"properties": ["volume", "area", "centroid"], "relative": "1e-15", "absolute": "1e-30"}
    })
}

#[test]
fn documented_complete_example_and_exact_answer() {
    let sandbox = Sandbox::new();
    let input = complete_example();
    let (manifest, text) = sandbox.answer_text(&input.to_string());
    assert_eq!(manifest["question_id"], "q-7f3a9c01");
    assert_eq!(manifest["status"], "answered");
    assert_eq!(manifest["exact"], true);
    assert_eq!(
        manifest["result"],
        json!({
            "kind": "solid", "lumps": 1, "shells": 1, "faces": 7, "edges": 14,
            "ring_edges": 2, "vertices": 8, "faces_by_kind": {"plane": 6, "cylinder": 1},
            "edges_by_kind": {"circle": 2, "line": 12}
        })
    );
    let properties = &manifest["properties"];
    keys(properties, &["volume", "area", "centroid"]);
    enclosure(
        &properties["volume"],
        "717.2566611769186085383",
        "717.2566611769186085384",
        "1e-15",
        "1e-30",
    );
    enclosure(
        &properties["area"],
        "731.9468914507713160154",
        "731.9468914507713160155",
        "1e-15",
        "1e-30",
    );
    assert_eq!(properties["centroid"].as_array().unwrap().len(), 3);
    for bound in properties["centroid"].as_array().unwrap() {
        enclosure(bound, "5", "5", "1e-15", "1e-30");
    }
    let step_report = &manifest["step"];
    assert_eq!(step_report["added_vertices"], 2);
    assert_eq!(step_report["added_edges"], 1);
    let deviation = Rational::parse(step_report["max_deviation"].as_str().unwrap()).unwrap();
    assert!(deviation >= Rational::zero());
    assert!(deviation <= Rational::parse("1e-7").unwrap());
    assert_eq!(manifest["marks"], json!([]));
    assert_eq!(manifest["refusal"], Value::Null);
    assert_eq!(manifest["stats"]["per_step"][0]["id"], "s1");
    assert_eq!(manifest["stats"]["per_step"].as_array().unwrap().len(), 1);
    let mut previous = 0;
    for key in [
        "format",
        "question_id",
        "kernel",
        "status",
        "exact",
        "result",
        "properties",
        "step",
        "marks",
        "refusal",
        "stats",
    ] {
        let position = text.find(&format!("\"{key}\":")).unwrap();
        assert!(position >= previous);
        previous = position;
    }
    let step = fs::read_to_string(sandbox.0.join("out/answer.step")).unwrap();
    assert!(step.starts_with("ISO-10303-21;"));
    assert!(step.contains("AP242"));
    assert!(step.contains("q-7f3a9c01"));
    assert!(step.contains("1970-01-01T00:00:00"));
    assert_eq!(step.matches("MANIFOLD_SOLID_BREP(").count(), 1);
    assert_eq!(step.matches("ADVANCED_FACE(").count(), 7);
    let mut rerun = sandbox.answer(&input);
    let mut first = manifest;
    rerun.as_object_mut().unwrap().remove("stats");
    first.as_object_mut().unwrap().remove("stats");
    assert_eq!(rerun, first);
    assert_eq!(
        fs::read_to_string(sandbox.0.join("out/answer.step")).unwrap(),
        step
    );
}

#[test]
fn three_step_chain_covers_all_operations() {
    let sandbox = Sandbox::new();
    let input = json!({
        "format": "exact-brep-question/2", "id": "chain",
        "operands": {
            "a": {"shape": "block", "size": ["4", "4", "4"]},
            "b": {"shape": "block", "size": ["4", "4", "4"], "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["2", "0", "0"]}},
            "c": {"shape": "block", "size": ["4", "4", "4"], "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["1", "0", "0"]}}
        },
        "steps": [
            {"id": "s1", "op": "unite", "args": ["a", "b"]},
            {"id": "s2", "op": "subtract", "args": ["s1", "c"]},
            {"id": "s3", "op": "intersect", "args": ["s2", "a"]}
        ], "result": "s3", "export": {"step": false}
    });
    let manifest = sandbox.answer(&input);
    assert_eq!(manifest["status"], "answered", "{manifest}");
    assert_eq!(
        manifest["result"],
        json!({
            "kind": "solid", "lumps": 1, "shells": 1, "faces": 6, "edges": 12,
            "ring_edges": 0, "vertices": 8, "faces_by_kind": {"plane": 6}, "edges_by_kind": {"line": 12}
        })
    );
    assert_eq!(
        manifest["properties"]["volume"],
        json!({"lo": "16", "hi": "16"})
    );
    assert_eq!(
        manifest["properties"]["area"],
        json!({"lo": "48", "hi": "48"})
    );
    assert_eq!(
        manifest["properties"]["centroid"],
        json!([
            {"lo": "0.5", "hi": "0.5"}, {"lo": "2", "hi": "2"}, {"lo": "2", "hi": "2"}
        ])
    );
    let ids: Vec<_> = manifest["stats"]["per_step"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| s["id"].as_str().unwrap())
        .collect();
    assert_eq!(ids, ["s1", "s2", "s3"]);
    assert!(manifest["step"].is_null());
    assert!(!sandbox.0.join("out/answer.step").exists());
}

#[test]
fn rational_quaternion_rotates_before_translation() {
    let sandbox = Sandbox::new();
    let mut input = question(json!({
        "shape": "block", "size": ["2", "4", "6"],
        "pose": {"rotation": ["2/3", "0", "0", "2/3"], "translation": ["10", "20", "30"]}
    }));
    input["export"] = json!({"step": false});
    input["measure"] = json!({"properties": ["centroid"]});
    let manifest = sandbox.answer(&input);
    assert_eq!(manifest["status"], "answered");
    assert_eq!(
        manifest["properties"],
        json!({"centroid": [
            {"lo": "8", "hi": "8"}, {"lo": "21", "hi": "21"}, {"lo": "33", "hi": "33"}
        ]})
    );
}

#[test]
fn sphere_topology_excludes_step_only_entities() {
    let sandbox = Sandbox::new();
    let mut input = question(json!({"shape": "sphere", "radius": "3"}));
    input["measure"] = json!({"properties": ["area"], "relative": "1e-15", "absolute": "1e-30"});
    input["export"] = json!({"step": true, "curve_tolerance": "1e-9", "significant_digits": 24});
    let manifest = sandbox.answer(&input);
    assert_eq!(manifest["status"], "answered");
    assert_eq!(
        manifest["result"],
        json!({
            "kind": "solid", "lumps": 1, "shells": 1, "faces": 1, "edges": 0,
            "ring_edges": 0, "vertices": 0, "faces_by_kind": {"sphere": 1}, "edges_by_kind": {}
        })
    );
    keys(&manifest["properties"], &["area"]);
    enclosure(
        &manifest["properties"]["area"],
        "113.0973355292325565846551617980621",
        "113.0973355292325565846551617980622",
        "1e-15",
        "1e-30",
    );
    assert_eq!(
        manifest["step"],
        json!({
            "file": "answer.step", "max_deviation": "0", "approximated_curves": 0,
            "added_vertices": 1, "added_edges": 0
        })
    );
}

#[test]
fn empty_answer_removes_stale_step() {
    let sandbox = Sandbox::new();
    let mut input = question(json!({"shape": "block", "size": ["2", "3", "4"]}));
    assert_eq!(sandbox.answer(&input)["status"], "answered");
    assert!(sandbox.0.join("out/answer.step").exists());
    input["steps"] = json!([{"id": "s1", "op": "subtract", "args": ["a", "a"]}]);
    input["result"] = json!("s1");
    let manifest = sandbox.answer(&input);
    assert_eq!(manifest["status"], "answered");
    assert_eq!(manifest["exact"], true);
    assert_eq!(
        manifest["result"],
        json!({
            "kind": "empty", "lumps": 0, "shells": 0, "faces": 0, "edges": 0,
            "ring_edges": 0, "vertices": 0, "faces_by_kind": {}, "edges_by_kind": {}
        })
    );
    assert_eq!(
        manifest["properties"],
        json!({
            "volume": {"lo": "0", "hi": "0"}, "area": {"lo": "0", "hi": "0"}, "centroid": null
        })
    );
    assert!(manifest["step"].is_null());
    assert!(!sandbox.0.join("out/answer.step").exists());
    assert_eq!(manifest["stats"]["per_step"][0]["decisions"], 0);
    assert_eq!(manifest["stats"]["per_step"][0]["not_exact"], 0);
}

#[test]
fn identical_bodies_and_aliased_names_follow_set_rules() {
    let sandbox = Sandbox::new();
    let mut input = question(json!({"shape": "block", "size": ["2", "3", "4"]}));
    let original = sandbox.answer(&input);
    for op in ["unite", "intersect"] {
        input["steps"] = json!([{"id": "s1", "op": op, "args": ["a", "a"]}]);
        input["result"] = json!("s1");
        let manifest = sandbox.answer(&input);
        assert_eq!(manifest["status"], "answered", "{manifest}");
        assert_eq!(manifest["result"], original["result"]);
        assert_eq!(manifest["properties"], original["properties"]);
        assert_eq!(manifest["marks"], original["marks"]);
        assert_eq!(manifest["stats"]["per_step"][0]["decisions"], 0);
        assert_eq!(manifest["stats"]["per_step"][0]["not_exact"], 0);
    }
    input["steps"] = json!([
        {"id": "s1", "op": "unite", "args": ["a", "a"]},
        {"id": "s2", "op": "intersect", "args": ["s1", "a"]},
        {"id": "s3", "op": "subtract", "args": ["s2", "a"]}
    ]);
    input["result"] = json!("s3");
    let manifest = sandbox.answer(&input);
    assert_eq!(manifest["status"], "answered", "{manifest}");
    assert_eq!(manifest["result"]["kind"], "empty");
    assert_eq!(
        manifest["properties"],
        json!({
            "volume": {"lo": "0", "hi": "0"},
            "area": {"lo": "0", "hi": "0"},
            "centroid": null
        })
    );
    assert_eq!(manifest["marks"], original["marks"]);
    assert!(manifest["step"].is_null());
    assert!(!sandbox.0.join("out/answer.step").exists());
    let steps = manifest["stats"]["per_step"].as_array().unwrap();
    assert_eq!(steps.len(), 3);
    for (step, id) in steps.iter().zip(["s1", "s2", "s3"]) {
        assert_eq!(step["id"], id);
        assert_eq!(step["decisions"], 0);
        assert_eq!(step["not_exact"], 0);
    }
}

#[test]
fn deadline_breach_is_a_normal_refusal() {
    let sandbox = Sandbox::new();
    let mut input = complete_example();
    input["limits"] = json!({"deadline_ms": 0});
    let manifest = sandbox.answer(&input);
    assert_eq!(manifest["status"], "refused");
    assert_eq!(manifest["refusal"]["code"], "numbers.deadline");
    assert!(!sandbox.0.join("out/answer.step").exists());
}

#[test]
fn boolean_deadlines_abandon_work_and_publish_only_a_refusal() {
    use std::time::{Duration, Instant};

    let pose = |translation| json!({"rotation": ["1", "0", "0", "0"], "translation": translation});
    let cylinder = json!({
        "shape": "cylinder", "radius": "2", "height": "8",
        "pose": pose(["0", "0", "-4"])
    });
    let cases = [
        ("3", ["1", "0.5", "0.3"]),
        ("3", ["1", "0", "0"]),
        ("2", ["1", "0", "0"]),
        ("3", ["4.9999999999", "0", "0"]),
        ("3", ["0", "0", "0"]),
    ];
    for (radius, translation) in cases {
        for op in ["unite", "subtract", "intersect"] {
            let operand = json!({
                "shape": "sphere", "radius": radius, "pose": pose(translation)
            });
            assert_deadline_refusal(&cylinder, &operand, op);
        }
    }
    let crossing = json!({
        "shape": "cylinder", "radius": "2", "height": "8",
        "pose": {"rotation": ["1", "0", "1", "0"], "translation": ["-4", "0", "0"]}
    });
    for op in ["unite", "subtract", "intersect"] {
        assert_deadline_refusal(&cylinder, &crossing, op);
    }

    fn assert_deadline_refusal(a: &Value, b: &Value, op: &str) {
        let sandbox = Sandbox::new();
        let question = json!({
            "format": "exact-brep-question/2", "id": "deadline-regression",
            "operands": {"a": a, "b": b},
            "steps": [{"id": "s1", "op": op, "args": ["a", "b"]}],
            "result": "s1", "limits": {"deadline_ms": 200}
        });
        fs::write(sandbox.0.join("question.json"), question.to_string()).unwrap();
        fs::create_dir(sandbox.0.join("out")).unwrap();
        fs::write(sandbox.0.join("out/answer.step"), "stale answer").unwrap();
        let start = Instant::now();
        let mut child = Command::new(env!("CARGO_BIN_EXE_exact-brep"))
            .args(["answer", "question.json", "--out", "out"])
            .env_clear()
            // The dev-only numbers feature freezes the deadline until the
            // boolean has entered numeric work, then expires it on check eight.
            .env("CATALLAX_TEST_DEADLINE_CHECKS", "8")
            .current_dir(&sandbox.0)
            .spawn()
            .unwrap();
        // Hang guard only, not a deadline latency or machine-speed assertion.
        let limit = Duration::from_secs(60);
        let status = loop {
            if let Some(status) = child.try_wait().unwrap() {
                break status;
            }
            if start.elapsed() > limit {
                child.kill().unwrap();
                child.wait().unwrap();
                panic!("answer did not exit after deterministic deadline expiry: {question}");
            }
            std::thread::sleep(Duration::from_millis(5));
        };
        assert_eq!(status.code(), Some(0), "{question}");
        let manifest: Value =
            serde_json::from_str(&fs::read_to_string(sandbox.0.join("out/manifest.json")).unwrap())
                .unwrap();
        envelope(&manifest);
        assert_eq!(manifest["status"], "refused", "{question}");
        assert_eq!(
            manifest["refusal"]["code"], "numbers.deadline",
            "{question}"
        );
        assert_eq!(manifest["exact"], false);
        assert!(manifest["result"].is_null());
        assert!(manifest["properties"].is_null());
        assert!(manifest["step"].is_null());
        assert!(!sandbox.0.join("out/answer.step").exists());
    }
}

#[test]
fn malformed_questions_are_invalid_not_crashes() {
    let sandbox = Sandbox::new();
    let base = question(json!({"shape": "block", "size": ["1", "2", "3"]}));
    let mutations = [
        ("/format", json!("exact-brep-question/1")),
        ("/id", json!("bad id")),
        ("/id", json!("x".repeat(201))),
        ("/operands", json!({})),
        (
            "/operands",
            json!({"UPPER": {"shape": "sphere", "radius": "1"}}),
        ),
        ("/operands/a/size", json!([1, "2", "3"])),
        ("/operands/a/size", json!(["1/0", "2", "3"])),
        ("/operands/a/size", json!(["0x10", "2", "3"])),
        ("/operands/a/size", json!(["0", "2", "3"])),
        ("/operands/a/size", json!(["1", "2"])),
        ("/operands/a", json!({"shape": "sphere", "radius": "-1"})),
        (
            "/operands/a",
            json!({"shape": "cone", "base_radius": "1", "top_radius": "1", "height": "2"}),
        ),
        (
            "/operands/a",
            json!({"shape": "torus", "major": "1", "minor": "2"}),
        ),
        (
            "/operands/a",
            json!({"shape": "sphere", "radius": "1", "pose": {"rotation": ["0", "0", "0", "0"], "translation": ["0", "0", "0"]}}),
        ),
        (
            "/operands/a",
            json!({"shape": "sphere", "radius": "1", "pose": {"rotation": [["1", "0", "0"], ["0", "1", "0"], ["0", "0", "1"]], "translation": ["0", "0", "0"]}}),
        ),
        (
            "/operands/a",
            json!({"shape": "sphere", "radius": "1", "pose": null}),
        ),
        ("/result", json!("missing")),
        (
            "/steps",
            json!([{"id": "s1", "op": "unite", "args": ["a", "later"]}]),
        ),
        (
            "/steps",
            json!([{"id": "a", "op": "unite", "args": ["a", "a"]}]),
        ),
        (
            "/steps",
            json!([{"id": "s1", "op": "xor", "args": ["a", "a"]}]),
        ),
    ];
    let mut invalid: Vec<String> = mutations
        .into_iter()
        .map(|(pointer, value)| {
            let mut q = base.clone();
            *q.pointer_mut(pointer).unwrap() = value;
            q.to_string()
        })
        .collect();
    for (key, value) in [
        ("unknown", json!(true)),
        ("measure", json!({"properties": ["mass"]})),
        ("measure", json!({"relative": 0.1})),
        ("measure", json!({"absolute": "-1"})),
        ("measure", json!({"relative": null})),
        ("measure", Value::Null),
        ("export", json!({"significant_digits": 5})),
        ("export", json!({"significant_digits": 41})),
        ("export", json!({"curve_tolerance": "0"})),
        ("export", json!({"step": "true"})),
        ("export", json!({"unknown": true})),
        ("limits", json!({"deadline_ms": -1})),
        ("limits", json!({"deadline_ms": 1.5})),
        ("limits", json!({"deadline_ms": 1, "unknown": true})),
    ] {
        let mut q = base.clone();
        q[key] = value;
        invalid.push(q.to_string());
    }
    invalid.extend([
        "{".into(),
        "{\"id\":\"already-read\",\"operands\":".into(),
        base.to_string()
            .replace("\"size\":", "\"extra\":true,\"size\":"),
        base.to_string()
            .replace("\"id\":", "\"id\":\"duplicate\",\"id\":"),
        base.to_string().replace(
            "\"a\":",
            "\"a\":{\"shape\":\"sphere\",\"radius\":\"1\"},\"a\":",
        ),
    ]);
    for text in invalid {
        fs::write(sandbox.0.join("question.json"), &text).unwrap();
        let check = sandbox.run(&["check", "question.json"]);
        assert_eq!(check.status.code(), Some(65), "{text}: {check:?}");
        assert!(check.stdout.is_empty());
        assert!(!check.stderr.is_empty());
        let manifest = sandbox.answer_text(&text).0;
        assert_eq!(manifest["status"], "invalid-question", "{text}: {manifest}");
        assert_eq!(manifest["refusal"]["code"], "cli.invalid-question");
        assert!(manifest["refusal"]["step"].is_null());
        assert!(manifest["refusal"]["stage"].is_null());
        if text.starts_with("{\"id\":\"already-read\"") {
            assert_eq!(manifest["question_id"], "already-read");
        }
    }
}

#[test]
fn usage_io_and_about_contracts() {
    let sandbox = Sandbox::new();
    let about = sandbox.run(&["about"]);
    assert!(about.status.success());
    assert!(about.stderr.is_empty());
    let about: Value = serde_json::from_slice(&about.stdout).unwrap();
    keys(
        &about,
        &[
            "name",
            "version",
            "commit",
            "question_formats",
            "manifest_format",
            "shapes",
            "operations",
            "properties",
        ],
    );
    assert_eq!(about["name"], "exact-brep");
    assert_eq!(about["question_formats"], json!(["exact-brep-question/2"]));
    assert_eq!(about["manifest_format"], "exact-brep-manifest/2");
    assert_eq!(
        about["shapes"],
        json!(["block", "cylinder", "sphere", "cone", "torus"])
    );
    assert_eq!(
        about["operations"],
        json!(["unite", "subtract", "intersect"])
    );
    for args in [
        vec![],
        vec!["answer"],
        vec!["about", "extra"],
        vec!["answer", "q", "--wrong", "out"],
    ] {
        assert_eq!(sandbox.run(&args).status.code(), Some(64));
        assert!(!sandbox.0.join("out/manifest.json").exists());
    }
    assert_eq!(sandbox.run(&["check", "missing"]).status.code(), Some(65));
    assert_eq!(
        sandbox
            .run(&["answer", "missing", "--out", "out"])
            .status
            .code(),
        Some(0)
    );
    let manifest: Value =
        serde_json::from_str(&fs::read_to_string(sandbox.0.join("out/manifest.json")).unwrap())
            .unwrap();
    assert_eq!(manifest["status"], "invalid-question");
    assert!(manifest["question_id"].is_null());
    fs::write(sandbox.0.join("not-a-directory"), "file").unwrap();
    assert_eq!(
        sandbox
            .run(&["answer", "missing", "--out", "not-a-directory/out"])
            .status
            .code(),
        Some(74)
    );
    fs::create_dir_all(sandbox.0.join("blocked/manifest.json")).unwrap();
    assert_eq!(
        sandbox
            .run(&["answer", "missing", "--out", "blocked"])
            .status
            .code(),
        Some(74)
    );
}

#[test]
fn precision_limits_and_decimal_rounding_are_marked() {
    let sandbox = Sandbox::new();
    let mut input = question(json!({"shape": "block", "size": ["1/3", "1", "1"]}));
    input["measure"] = json!({"properties": ["volume"], "relative": "0", "absolute": "0"});
    input["export"] = json!({"step": false});
    let manifest = sandbox.answer(&input);
    assert_eq!(manifest["status"], "answered");
    assert_eq!(manifest["exact"], false);
    assert_eq!(manifest["marks"][0]["code"], "measure.accuracy-not-reached");
    let lo = rational(&manifest["properties"]["volume"]["lo"]);
    let hi = rational(&manifest["properties"]["volume"]["hi"]);
    let third = Rational::parse("1/3").unwrap();
    assert!(lo < third && third < hi);
    assert!(rational(&manifest["marks"][0]["bound"]) >= hi - lo);
    input["measure"] = json!({"properties": ["volume"], "relative": "1e-90", "absolute": "1e-100"});
    let precise = sandbox.answer(&input);
    assert_eq!(precise["status"], "answered");
    assert_eq!(precise["exact"], true);
    let lo = rational(&precise["properties"]["volume"]["lo"]);
    let hi = rational(&precise["properties"]["volume"]["hi"]);
    assert!(lo < third && third < hi);
    assert!(&hi - &lo <= Rational::parse("1e-90").unwrap() * lo);
    input["operands"]["a"]["size"][0] = json!("1000000000000000000000000000000/3");
    input["measure"] = json!({"properties": []});
    input["export"] = json!({"step": true, "significant_digits": 6, "curve_tolerance": "1e-7"});
    let rounded = sandbox.answer(&input);
    assert_eq!(rounded["status"], "answered", "{rounded}");
    assert_eq!(rounded["exact"], false);
    assert!(
        rounded["marks"]
            .as_array()
            .unwrap()
            .iter()
            .any(|mark| mark["code"] == "step.tolerance-not-reached")
    );
    assert!(rational(&rounded["step"]["max_deviation"]) > Rational::parse("1e-7").unwrap());
    assert_eq!(rounded["properties"], json!({}));
}
