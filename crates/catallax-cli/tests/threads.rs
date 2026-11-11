//! The kernel's thread count must change nothing a caller can see except
//! elapsed time: the manifest, the operation log, the work count, a refusal
//! when the allowance runs out, and the bytes of every file are those of one
//! thread.
use serde_json::{Value, json};
use std::{
    fs,
    io::Write,
    path::Path,
    process::{Command, Stdio},
};

// Small enough to run ten times, wide enough to pass through the Boolean
// stages, a fillet-free STEP export with curved pcurves, and a refusal.
const RECIPES: [(&str, &str); 6] = [
    (
        "line-tangent",
        include_str!("../../../parts/P0/line-tangent.json"),
    ),
    (
        "near-miss",
        include_str!("../../../parts/P0/near-miss.json"),
    ),
    (
        "equal-cylinder-sphere",
        include_str!("../../../parts/P0/equal-cylinder-sphere.json"),
    ),
    (
        "general-cone",
        include_str!("../../../parts/P0/general-cone.json"),
    ),
    ("collar", include_str!("../../../parts/P2-collar.json")),
    ("tray", include_str!("../../../parts/P4-shelled-tray.json")),
];

// A waiting run answers the recipes before the tray. It keeps a record of
// every decision and proves each match by walking the whole question, which
// takes seconds on the small recipes and many minutes on the tray.
const BEFORE_THE_TRAY: usize = RECIPES.len() - 1;

// The full allowance, then allowances that run out inside the Boolean and
// inside the export of the larger parts.
const ALLOWANCES: [u64; 5] = [100_000_000, 1_000_000, 100_000, 20_000, 3_000];

fn requests(dir: &Path, recipes: usize) -> Vec<(String, Value)> {
    let mut requests = Vec::new();
    for (name, recipe) in &RECIPES[..recipes] {
        let problem: Value = serde_json::from_str(recipe).unwrap();
        for allowance in ALLOWANCES {
            let id = format!("{name}-{allowance}");
            requests.push((
                id.clone(),
                json!({"format":"catallax-request/1","source":"recipe","problem":problem,
                    "max_work":allowance,"out":dir.join(&id)}),
            ));
        }
    }
    requests
}

// Elapsed time is the one thing allowed to differ.
fn without_clocks(value: &mut Value) {
    match value {
        Value::Object(map) => {
            map.remove("wall_ms");
            map.values_mut().for_each(without_clocks);
        }
        Value::Array(values) => values.iter_mut().for_each(without_clocks),
        _ => {}
    }
}

fn files(dir: &Path) -> Vec<(String, Vec<u8>)> {
    let mut found = Vec::new();
    let Ok(entries) = fs::read_dir(dir) else {
        return found;
    };
    for entry in entries {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_string_lossy().into_owned();
        // Manifests repeat the response; the response is compared without clocks.
        if path.is_file() && !name.ends_with(".json") {
            found.push((name, fs::read(&path).unwrap()));
        }
    }
    found.sort();
    found
}

type Run = Vec<(String, Value, Vec<(String, Vec<u8>)>)>;

// `mode` is a `CATALLAX_RUN_AHEAD` value. Without one, only the STEP pcurve
// fit uses the other threads. `loops` also hands the Boolean's loops to
// workers; `wait` makes every loop wait for them, so every record is adopted.
fn run(label: &str, threads: Option<usize>, mode: Option<&str>) -> Run {
    let recipes = if mode == Some("wait") {
        BEFORE_THE_TRAY
    } else {
        RECIPES.len()
    };
    let dir = std::env::temp_dir().join(format!("catallax-threads-{}-{label}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let mut command = Command::new(env!("CARGO_BIN_EXE_catallax-verify"));
    command.env_remove("CATALLAX_THREADS");
    command.env_remove("CATALLAX_RUN_AHEAD");
    if let Some(threads) = threads {
        command.args(["--threads", &threads.to_string()]);
    }
    if let Some(mode) = mode {
        command.env("CATALLAX_RUN_AHEAD", mode);
    }
    let mut child = command
        .arg("stream")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let requests = requests(&dir, recipes);
    let mut input = child.stdin.take().unwrap();
    for (_, request) in &requests {
        writeln!(input, "{request}").unwrap();
    }
    drop(input);
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success(), "{label}");
    let lines = String::from_utf8(output.stdout).unwrap();
    let rows: Vec<Value> = lines
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(rows.len(), requests.len(), "{label}");
    let run = requests
        .into_iter()
        .zip(rows)
        .map(|((id, _), mut row)| {
            without_clocks(&mut row);
            let files = files(&dir.join(&id));
            (id, row, files)
        })
        .collect();
    fs::remove_dir_all(dir).unwrap();
    run
}

#[test]
fn thread_count_changes_only_elapsed_time() {
    let one = run("one", Some(1), None);
    let answered = one
        .iter()
        .filter(|(_, row, _)| row["manifest"]["status"] == "answered")
        .count();
    let out_of_work = one
        .iter()
        .filter(|(_, row, _)| row["work"]["exhausted"] == true)
        .count();
    assert!(answered >= 6, "answers to compare: {answered}");
    assert!(
        out_of_work >= 6,
        "exhausted allowances to compare: {out_of_work}"
    );
    for (label, threads, mode) in [
        ("two", Some(2), None),
        ("four", Some(4), None),
        ("eight", Some(8), None),
        ("default", None, None),
        ("two-waiting", Some(2), Some("wait")),
        ("four-waiting", Some(4), Some("wait")),
        ("two-loops", Some(2), Some("loops")),
        ("eight-loops", Some(8), Some("loops")),
        ("default-loops", None, Some("loops")),
    ] {
        let other = run(label, threads, mode);
        assert!(other.len() >= BEFORE_THE_TRAY * ALLOWANCES.len(), "{label}");
        for (expected, actual) in one.iter().zip(&other) {
            assert_eq!(expected.0, actual.0, "{label}");
            assert_eq!(expected.1, actual.1, "{label}: response of {}", expected.0);
            assert!(expected.2 == actual.2, "{label}: files of {}", expected.0);
        }
    }
}

#[test]
fn a_thread_count_must_be_positive() {
    for count in ["0", "many"] {
        let output = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
            .args(["--threads", count, "about"])
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(64), "--threads {count}");
    }
    let output = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
        .args(["--threads", "3", "about"])
        .output()
        .unwrap();
    assert!(output.status.success());
}
