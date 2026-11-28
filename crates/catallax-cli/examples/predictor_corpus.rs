//! Development-only runner of the predictor's frozen corpus. Never linked into
//! an application entry.
//!
//! Checkpoint 2 decisions, section 4.5, state the predictor's done-when lines
//! on a frozen corpus. `tools/predictor/corpus.json` is that corpus and
//! `docs/archive/milestones/3.0B-predictor.md` says how it was chosen. `run` builds
//! every row through the application's ordinary build, the one that tries the
//! exact builders and then the fitted construction, exports the result, and
//! writes one line for each row and thread count. `check` reads those lines
//! and tests the done-when lines.
use catallax_app::preview::{display_description, numerical_display_description};
use catallax_app::{
    DEFAULT_BUILD_WORK, EvaluatedApplicationOutcome, TARGET_BUILD_WORK,
    build_incremental_with_files_without_display,
};
use catallax_features::{ImportFiles, predictor};
use catallax_mesh::float::DisplayOptions;
use catallax_numbers::{Decider, OperationWork, with_operation_budget, with_operation_log};
use catallax_recipe::Recipe;
use catallax_step::{StepOptions, write_step};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    env, fs,
    io::Write,
    num::NonZeroUsize,
    path::{Path, PathBuf},
    time::Instant,
};

const USAGE: &str = "usage:
  predictor_corpus run --corpus FILE --out FILE [--root DIR] [--sources DIR --selection FILE]
                       [--threads 1,2,4,24] [--repeat N] [--rows ID,ID] [--label TEXT]
                       [--describe | --display | --require-preview] [--max-work N]
  predictor_corpus check --corpus FILE FILE... [--base FILE]
  predictor_corpus census-preview --sources DIR --selection FILE [--max-work N]

run   builds and exports every row at every thread count, N times each, and
      appends one JSON line for each to --out. Rows that import a private
      source need --sources and --selection, as tools/checkpoint/checkpoint.py
      takes them; without them those rows are skipped and named.
      --describe records unproved display descriptions in place of STEP export.
      --display records the preview and certified mesh in place of STEP export.
      --require-preview checks that every built row has an unproved preview.
      census-preview applies the same preview check to a private STEP selection,
      reporting aggregate counts only.
check tests the predictor's done-when lines on the lines of the runs. It
      needs every row of --corpus at 1, 2, 4 and 24 threads on two machines,
      a machine being a --label, and fails on any that is missing. --base is
      a run of the same corpus on the revision the predictor is compared
      with, and it too must hold every row.";

/// The thread counts of the done-when line: "path choice and work are
/// identical at 1, 2, 4 and 24 threads and on two machines" (Checkpoint 2
/// decisions, section 4.5).
const THREADS: [u64; 4] = [1, 2, 4, 24];

/// The machines of the same line.
const MACHINES: usize = 2;

/// A row's time is read only where its best time reaches this many
/// milliseconds. Below it the clock reads process and allocation overhead,
/// which is not arithmetic and not what the band line is about.
const BAND_FLOOR_MS: f64 = 10.0;

/// The band the decisions allow between the fastest and the slowest row, in
/// wall time per million work units.
const BAND: f64 = 10.0;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("run") => run(&args[1..]),
        Some("check") => check(&args[1..]),
        Some("census-preview") => census_preview(&args[1..]),
        _ => {
            eprintln!("{USAGE}");
            std::process::exit(64);
        }
    }
}

fn flag<'a>(args: &'a [String], name: &str) -> Option<&'a str> {
    args.iter()
        .position(|arg| arg == name)
        .and_then(|at| args.get(at + 1))
        .map(String::as_str)
}

fn sha256(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

/// A private source a row imports: its SHA-256 and its bytes.
type Source = (String, Vec<u8>);

/// A row ready to run: its recipe and the private source it imports, if any.
type Row = (Value, Option<Source>);

/// A row's recipe and the private source it imports, if any.
fn recipe_of(
    row: &Value,
    root: &Path,
    sources: Option<&BTreeMap<String, PathBuf>>,
) -> Result<Option<Row>, Box<dyn std::error::Error>> {
    let read = |path: &str| -> Result<Value, Box<dyn std::error::Error>> {
        Ok(serde_json::from_slice(&fs::read(root.join(path))?)?)
    };
    let mut source = None;
    let mut recipe = if let Some(path) = row["recipe"].as_str() {
        read(path)?
    } else if row["inline"].is_object() {
        row["inline"].clone()
    } else if let Some(path) = row["cases"].as_str() {
        let case = row["case"].as_str().ok_or("a cases row names its case")?;
        let text = fs::read_to_string(root.join(path))?;
        let mut found = None;
        for line in text.lines().filter(|line| !line.trim().is_empty()) {
            let value: Value = serde_json::from_str(line)?;
            if value["id"] == case {
                found = Some(value["recipe"].clone());
            }
        }
        found.ok_or_else(|| format!("{path} has no case {case}"))?
    } else if let Some(model) = row["import"]["model"].as_str() {
        let Some(path) = sources.and_then(|sources| sources.get(model)) else {
            return Ok(None);
        };
        let bytes = fs::read(path)?;
        let hash = sha256(&bytes);
        if row["import"]["sha256"] != hash.as_str() {
            return Err(format!("{model}: source bytes differ from the corpus").into());
        }
        source = Some((hash.clone(), bytes));
        import_recipe(&hash, &format!("ABC {model}"))
    } else {
        return Err("a row has a recipe, an inline recipe, a case or an import".into());
    };
    if let Some(steps) = row["append"].as_array() {
        recipe["steps"]
            .as_array_mut()
            .ok_or("a recipe has steps")?
            .extend(steps.iter().cloned());
        recipe["result"] = row["result"].clone();
    }
    if let Some(format) = row["format"].as_str() {
        recipe["format"] = json!(format);
    }
    Ok(Some((recipe, source)))
}

fn import_recipe(hash: &str, name: &str) -> Value {
    json!({"format": "catallax-recipe/3", "name": name,
        "steps": [{"id": "imported", "op": "import_step",
            "file": {"name": "part.step", "sha256": hash},
            "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]}}],
        "result": "imported"})
}

#[derive(Clone, Copy)]
enum Report {
    Step,
    Description,
    Preview,
    Display,
}

fn record_mesh<T: std::fmt::Debug>(mesh: impl FnOnce() -> Result<T, String>) -> Value {
    let started = Instant::now();
    let (answer, work) = with_operation_budget(u64::MAX, mesh);
    let mut record = work_json(work);
    record["ms"] = json!(started.elapsed().as_secs_f64() * 1000.0);
    match answer {
        Ok(mesh) if !work.exhausted => {
            record["sha256"] = json!(sha256(format!("{mesh:?}").as_bytes()));
        }
        Ok(_) => record["refusal"] = json!("numbers.operation-budget"),
        Err(error) => record["refusal"] = json!(error),
    }
    record
}

fn record_display(outcome: &EvaluatedApplicationOutcome, certified: bool) -> Value {
    let options = catallax_mesh::PreviewOptions::default();
    match outcome {
        EvaluatedApplicationOutcome::Exact(body) => {
            let built = body.checkpoint.built();
            let preview = record_mesh(|| {
                catallax_app::preview_body(built, &options, &mut Decider::default())
                    .map(|display| display.mesh)
                    .map_err(|error| error.to_string())
            });
            let mesh = certified.then(|| {
                let options = catallax_mesh::MeshOptions::default();
                let mut bounds = Value::Null;
                let mut record = record_mesh(|| {
                    let mesh = catallax_mesh::mesh_body(
                        &built.model, built.result, &options, &mut Decider::default(),
                    ).map_err(|error| error.to_string())?.value;
                    bounds = json!({
                        "vertices": mesh.positions.len(),
                        "triangles": mesh.triangles.len(),
                        "face_chord": mesh.faces.iter().map(|face| &face.max_chord_error).max().map(ToString::to_string),
                        "edge_chord": mesh.edges.iter().map(|edge| &edge.max_chord_error).max().map(ToString::to_string),
                    });
                    Ok(mesh)
                });
                record["bounds"] = bounds;
                record
            });
            json!({"preview": preview, "certified": mesh})
        }
        EvaluatedApplicationOutcome::Numerical(body) => {
            let preview = record_mesh(|| {
                catallax_mesh::preview_numerical_body(
                    body.body(),
                    &options,
                    &mut Decider::default(),
                )
                .map_err(|error| error.to_string())
            });
            json!({"preview": preview})
        }
    }
}

fn missing_preview(line: &Value) -> bool {
    matches!(line["build"]["status"].as_str(), Some("exact" | "fitted"))
        && !line["display"]["preview"]["sha256"].is_string()
}

fn census_preview(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let sources = Path::new(flag(args, "--sources").ok_or("census-preview needs --sources")?);
    let selection = flag(args, "--selection").ok_or("census-preview needs --selection")?;
    let members: Vec<(String, String)> = serde_json::from_slice(&fs::read(selection)?)?;
    let max_work = flag(args, "--max-work").map_or(Ok(DEFAULT_BUILD_WORK), str::parse)?;
    let mut built = 0;
    let mut missing = 0;
    let mut unreadable = 0;
    for (_, member) in &members {
        let Ok(bytes) = fs::read(sources.join(member)) else {
            unreadable += 1;
            continue;
        };
        let hash = sha256(&bytes);
        let line = measure(
            &import_recipe(&hash, "STEP import"),
            Some(&(hash, bytes)),
            max_work,
            Report::Preview,
        );
        if matches!(line["build"]["status"].as_str(), Some("exact" | "fitted")) {
            built += 1;
            missing += usize::from(missing_preview(&line));
        }
    }
    println!(
        "{}",
        json!({"files": members.len(), "built": built, "missing_preview": missing, "unreadable": unreadable})
    );
    if missing != 0 || unreadable != 0 {
        return Err("census did not preview every built file".into());
    }
    Ok(())
}

fn work_json(work: OperationWork) -> Value {
    json!({"used": work.used, "refunded": work.refunded, "uncharged": work.uncharged,
        "exhausted": work.exhausted})
}

/// Builds and exports one row on a thread of its own, as a request runs.
fn measure(recipe: &Value, source: Option<&Source>, max_work: u64, report: Report) -> Value {
    let text = recipe.to_string();
    std::thread::scope(|scope| {
        std::thread::Builder::new()
            .name("catallax-request".into())
            .stack_size(8 * 1024 * 1024)
            .spawn_scoped(scope, || measure_here(&text, source, max_work, report))
            .expect("a request thread")
            .join()
            .unwrap_or_else(|_| json!({"build": {"status": "panicked"}}))
    })
}

fn measure_here(text: &str, source: Option<&Source>, max_work: u64, report: Report) -> Value {
    let recipe = match Recipe::parse(text) {
        Ok(recipe) => recipe,
        Err(error) => {
            return json!({"build": {"status": "invalid", "detail": format!("{error:?}")}});
        }
    };
    let mut files = ImportFiles::new();
    if let Some((hash, bytes)) = source {
        files.insert(hash.clone(), bytes.as_slice());
    }
    let started = Instant::now();
    let (((outcome, work), entries), choices) = predictor::with_log(|| {
        with_operation_log(|| {
            with_operation_budget(max_work, || {
                build_incremental_with_files_without_display(
                    &recipe,
                    None,
                    &files,
                    &mut Decider::default(),
                )
            })
        })
    });
    let build_ms = started.elapsed().as_secs_f64() * 1000.0;
    let operations: Vec<Value> = entries
        .iter()
        .map(|entry| {
            json!({"id": entry.id, "operation": entry.operation,
                "used": entry.work_after - entry.work_before,
                "refusal": entry.refusal.as_ref().map(|refusal| refusal.code)})
        })
        .collect();
    let choices: Vec<Value> = choices
        .iter()
        .map(|choice| {
            let look = choice.look.map(|look| {
                json!({"degree": look.degree, "limbs": look.limbs, "score": look.score(),
                    "tier": look.tier().name(), "cap": look.tier().cap(), "read": look.read})
            });
            json!({"step": choice.step, "look": look, "path": choice.path.name(),
                "attempt_work": choice.attempt_work, "operation_work": choice.operation_work})
        })
        .collect();
    let (status, refusal) = match &outcome {
        _ if work.exhausted => ("refused", Some("numbers.operation-budget".to_owned())),
        Ok(EvaluatedApplicationOutcome::Exact(_)) => ("exact", None),
        Ok(EvaluatedApplicationOutcome::Numerical(_)) => ("fitted", None),
        Err(catallax_app::Error::Build(catallax_features::BuildError::Geometry {
            refusal,
            ..
        })) => ("refused", Some(refusal.code.to_owned())),
        Err(error) => ("refused", Some(error.code().to_owned())),
    };
    let mut build = work_json(work);
    build["status"] = json!(status);
    build["refusal"] = json!(refusal);
    build["ms"] = json!(build_ms);
    if !matches!(report, Report::Step) {
        let display = match outcome {
            Ok(body) if !work.exhausted => match report {
                Report::Description => describe_outcome(&body),
                Report::Preview | Report::Display => {
                    record_display(&body, matches!(report, Report::Display))
                }
                Report::Step => unreachable!(),
            },
            _ => Value::Null,
        };
        return json!({"build": build, "display": display, "operations": operations, "choices": choices});
    }
    let export = if work.exhausted {
        Value::Null
    } else {
        match &outcome {
            Ok(EvaluatedApplicationOutcome::Exact(body)) => {
                let built = body.checkpoint.built();
                let options = StepOptions {
                    product_name: "Part1".into(),
                    ..StepOptions::default()
                };
                let mut bytes = Vec::new();
                let started = Instant::now();
                let (answer, work) = with_operation_budget(max_work, || {
                    write_step(
                        &built.model,
                        &[built.result],
                        &options,
                        &mut Decider::default(),
                        &mut bytes,
                    )
                });
                exported(
                    work,
                    started,
                    answer.as_ref().err().map(|refusal| refusal.code),
                    &bytes,
                )
            }
            Ok(EvaluatedApplicationOutcome::Numerical(body)) => {
                let started = Instant::now();
                let (step, work) = with_operation_budget(max_work, || {
                    catallax_step::numerical::build(
                        body.body(),
                        &Default::default(),
                        &mut Decider::default(),
                    )
                });
                let bytes = step.as_ref().map_or(&[][..], |step| step.bytes());
                exported(
                    work,
                    started,
                    step.as_ref().err().map(|refusal| refusal.code),
                    bytes,
                )
            }
            Err(_) => Value::Null,
        }
    };
    json!({"build": build, "export": export, "operations": operations, "choices": choices})
}

fn describe_outcome(outcome: &EvaluatedApplicationOutcome) -> Value {
    fn record<T: std::fmt::Debug>(read: impl FnOnce() -> Result<T, catallax_app::Error>) -> Value {
        let started = Instant::now();
        let (answer, work) = with_operation_budget(0, read);
        let mut report = work_json(work);
        report["ms"] = json!(started.elapsed().as_secs_f64() * 1000.0);
        match answer {
            Ok(description) => {
                // Rust's shortest round-trip representation preserves each
                // binary64 value. This fingerprint is for corpus comparisons,
                // not a wire encoding of the display API.
                report["sha256"] = json!(sha256(format!("{description:?}").as_bytes()));
            }
            Err(error) => report["refusal"] = json!(format!("{error:?}")),
        }
        report
    }
    let options = DisplayOptions::default();
    match outcome {
        EvaluatedApplicationOutcome::Exact(body) => {
            record(|| display_description(body.checkpoint.built(), &options))
        }
        EvaluatedApplicationOutcome::Numerical(body) => {
            record(|| numerical_display_description(body, &options))
        }
    }
}

fn exported(work: OperationWork, started: Instant, refusal: Option<&str>, bytes: &[u8]) -> Value {
    let ms = started.elapsed().as_secs_f64() * 1000.0;
    let mut export = work_json(work);
    let refusal = if work.exhausted {
        Some("numbers.operation-budget")
    } else {
        refusal
    };
    export["status"] = json!(if refusal.is_some() {
        "refused"
    } else {
        "written"
    });
    export["refusal"] = json!(refusal);
    export["sha256"] = json!(refusal.is_none().then(|| sha256(bytes)));
    export["ms"] = json!(ms);
    export
}

fn run(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let corpus: Value = serde_json::from_slice(&fs::read(
        flag(args, "--corpus").ok_or("run needs --corpus")?,
    )?)?;
    if corpus["format"] != "catallax-predictor-corpus/1" {
        return Err("expected catallax-predictor-corpus/1".into());
    }
    let out = flag(args, "--out").ok_or("run needs --out")?;
    let root = flag(args, "--root").map_or_else(
        || Path::new(env!("CARGO_MANIFEST_DIR")).join("../.."),
        PathBuf::from,
    );
    let threads: Vec<usize> = flag(args, "--threads")
        .unwrap_or("1,2,4,24")
        .split(',')
        .map(str::parse)
        .collect::<Result<_, _>>()?;
    let repeat: usize = flag(args, "--repeat").map_or(Ok(1), str::parse)?;
    let only: Option<Vec<&str>> = flag(args, "--rows").map(|rows| rows.split(',').collect());
    let label = flag(args, "--label").unwrap_or("");
    let require_preview = args.iter().any(|arg| arg == "--require-preview");
    let report = if require_preview {
        Report::Preview
    } else if args.iter().any(|arg| arg == "--display") {
        Report::Display
    } else if args.iter().any(|arg| arg == "--describe") {
        Report::Description
    } else {
        Report::Step
    };
    let max_work = flag(args, "--max-work")
        .map(str::parse::<u64>)
        .transpose()?;
    let sources = match (flag(args, "--sources"), flag(args, "--selection")) {
        (Some(sources), Some(selection)) => {
            let members: Vec<(String, String)> = serde_json::from_slice(&fs::read(selection)?)?;
            Some(
                members
                    .into_iter()
                    .map(|(model, member)| (model, Path::new(sources).join(member)))
                    .collect::<BTreeMap<_, _>>(),
            )
        }
        _ => None,
    };
    let mut file = fs::OpenOptions::new().create(true).append(true).open(out)?;
    let mut missing = 0;
    for row in corpus["rows"].as_array().ok_or("a corpus has rows")? {
        let id = row["id"].as_str().ok_or("a row has an id")?;
        if only.as_ref().is_some_and(|only| !only.contains(&id)) {
            continue;
        }
        let Some((recipe, source)) = recipe_of(row, &root, sources.as_ref())? else {
            eprintln!("{id}: skipped, its private source was not given");
            continue;
        };
        let max_work = max_work
            .or_else(|| row["max_work"].as_u64())
            .unwrap_or(DEFAULT_BUILD_WORK);
        for &count in &threads {
            catallax_numbers::set_threads(NonZeroUsize::new(count));
            for pass in 0..repeat {
                let mut line = measure(&recipe, source.as_ref(), max_work, report);
                if require_preview && missing_preview(&line) {
                    missing += 1;
                    eprintln!(
                        "{id}: missing preview: {}",
                        line["display"]["preview"]["refusal"]
                    );
                }
                line["id"] = json!(id);
                line["group"] = row["group"].clone();
                line["threads"] = json!(count);
                line["pass"] = json!(pass);
                line["max_work"] = json!(max_work);
                line["label"] = json!(label);
                writeln!(file, "{line}")?;
                file.flush()?;
                eprintln!(
                    "{id} threads={count} pass={pass} {} used={} {:.0} ms",
                    line["build"]["status"].as_str().unwrap_or("?"),
                    line["build"]["used"],
                    line["build"]["ms"].as_f64().unwrap_or(0.0),
                );
            }
        }
    }
    catallax_numbers::set_threads(None);
    if missing != 0 {
        return Err(format!("{missing} built corpus rows have no preview").into());
    }
    Ok(())
}

/// What must be the same in every line of one row: everything but the clock,
/// the thread count and the labels.
fn answer_of(line: &Value) -> Value {
    let mut answer = line.clone();
    for key in ["threads", "pass", "label"] {
        answer.as_object_mut().expect("a line").remove(key);
    }
    for part in ["build", "export", "display"] {
        if let Some(part) = answer[part].as_object_mut() {
            part.remove("ms");
        }
    }
    answer
}

fn units(part: &Value) -> u64 {
    ["used", "refunded", "uncharged"]
        .iter()
        .map(|key| part[*key].as_u64().unwrap_or(0))
        .sum()
}

fn lines(path: &str) -> Result<Vec<Value>, Box<dyn std::error::Error>> {
    fs::read_to_string(path)?
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| Ok(serde_json::from_str(line)?))
        .collect()
}

fn check(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let base = flag(args, "--base");
    let corpus: Value = serde_json::from_slice(&fs::read(
        flag(args, "--corpus").ok_or("check needs --corpus")?,
    )?)?;
    let expected: Vec<&str> = corpus["rows"]
        .as_array()
        .ok_or("a corpus has rows")?
        .iter()
        .map(|row| row["id"].as_str().ok_or("a row has an id"))
        .collect::<Result<_, _>>()?;
    let mut all = Vec::new();
    let mut skip = false;
    for arg in args {
        if skip {
            skip = false;
        } else if arg == "--base" || arg == "--corpus" {
            skip = true;
        } else {
            all.extend(lines(arg)?);
        }
    }
    let mut rows: BTreeMap<String, Vec<&Value>> = BTreeMap::new();
    for line in &all {
        rows.entry(line["id"].as_str().unwrap_or("?").to_owned())
            .or_default()
            .push(line);
    }
    let mut failures = Vec::new();

    // 0. Coverage: every row of the corpus, on each machine, at each thread
    // count. A row, a machine or a thread count that was not run is a
    // failure, and so is a row the corpus does not have.
    let machines: std::collections::BTreeSet<&str> = all
        .iter()
        .map(|line| line["label"].as_str().unwrap_or(""))
        .collect();
    if machines.len() != MACHINES {
        failures.push(format!(
            "the runs are of {} machines ({machines:?}), and the line asks for {MACHINES}",
            machines.len()
        ));
    }
    let mut missing = 0;
    for id in &expected {
        for machine in &machines {
            for threads in THREADS {
                let ran = rows.get(*id).is_some_and(|lines| {
                    lines.iter().any(|line| {
                        line["label"].as_str().unwrap_or("") == *machine
                            && line["threads"] == threads
                    })
                });
                if !ran {
                    missing += 1;
                    // The first few name what is missing; the count says how much.
                    if missing <= 12 {
                        failures.push(format!(
                            "{id}: no line at {threads} threads on machine {machine:?}"
                        ));
                    }
                }
            }
        }
    }
    if missing > 12 {
        failures.push(format!("{missing} lines are missing in all"));
    }
    for id in rows.keys() {
        if !expected.contains(&id.as_str()) {
            failures.push(format!("{id}: not a row of the corpus"));
        }
    }
    println!(
        "coverage: {} rows of the corpus, machines {machines:?}, threads {THREADS:?}, {missing} lines missing",
        expected.len()
    );

    // 1. Path choice and work are the same in every line of a row.
    let mut runs = 0;
    for (id, lines) in &rows {
        let first = answer_of(lines[0]);
        runs = runs.max(lines.len());
        for line in &lines[1..] {
            if answer_of(line) != first {
                failures.push(format!(
                    "{id}: threads {} label {:?} differs from threads {} label {:?}",
                    line["threads"], line["label"], lines[0]["threads"], lines[0]["label"]
                ));
            }
        }
    }
    println!(
        "same path and work in every line: {} rows, up to {runs} lines a row",
        rows.len()
    );

    // 2. The quick look, against the work of its operation. 3. Refunded work
    // of a capped attempt, against its cap. An attempt that ended before the
    // high tier's cap was never looked at: it read nothing, and what it gave
    // back is under that cap.
    let (mut attempts, mut looks) = (0, 0);
    let (mut worst_share, mut worst_row) = (0.0_f64, String::new());
    let mut paths: BTreeMap<String, usize> = BTreeMap::new();
    for (id, lines) in &rows {
        let line = lines[0];
        for choice in line["choices"].as_array().into_iter().flatten() {
            attempts += 1;
            let path = choice["path"].as_str().unwrap_or("?");
            let tier = match choice["look"]["tier"].as_str() {
                Some(tier) => format!("{tier} tier"),
                None => "not looked at".into(),
            };
            *paths.entry(format!("{tier}, {path}")).or_default() += 1;
            let attempt = choice["attempt_work"].as_u64().unwrap_or(0);
            let cap = choice["look"]["cap"]
                .as_u64()
                .unwrap_or(predictor::HIGH_COMPLEXITY_CAP);
            if path != "exact" && attempt > cap {
                failures.push(format!(
                    "{id} step {}: refunded {attempt} past its cap of {cap}",
                    choice["step"]
                ));
            }
            let Some(read) = choice["look"]["read"].as_u64() else {
                continue;
            };
            looks += 1;
            let operation = choice["operation_work"].as_u64().unwrap_or(0);
            let share = read as f64 / operation.max(1) as f64;
            if share > worst_share {
                worst_share = share;
                worst_row = id.clone();
            }
            if read.saturating_mul(100) > operation {
                failures.push(format!(
                    "{id} step {}: the quick look read {read} against {operation} units of work",
                    choice["step"]
                ));
            }
        }
    }
    println!(
        "quick look: {attempts} attempts, {looks} looked at, largest share {:.4}% ({worst_row})",
        worst_share * 100.0
    );
    for (path, count) in &paths {
        println!("  {path}: {count}");
    }

    // 4. Wall time per million units, over rows the clock can read. The
    // clock is read at one thread, where the time is the arithmetic's, and
    // on one machine at a time: two machines' clocks are not one scale.
    for machine in &machines {
        let mut rates: Vec<(f64, String)> = Vec::new();
        for (id, lines) in &rows {
            for part in ["build", "export"] {
                let best = lines
                    .iter()
                    .filter(|line| line["threads"] == 1)
                    .filter(|line| line["label"].as_str().unwrap_or("") == *machine)
                    .filter_map(|line| line[part]["ms"].as_f64())
                    .fold(f64::INFINITY, f64::min);
                let counted = units(&lines[0][part]);
                if best.is_finite() && best >= BAND_FLOOR_MS && counted > 0 {
                    rates.push((best / (counted as f64 / 1e6), format!("{id} {part}")));
                }
            }
        }
        rates.sort_by(|a, b| a.0.total_cmp(&b.0));
        if let (Some(low), Some(high)) = (rates.first(), rates.last()) {
            let band = high.0 / low.0;
            println!(
                "wall time per million units on {machine:?}, {} measurements of {BAND_FLOOR_MS} ms or more: {:.1} ({}) to {:.1} ({}), band {band:.1}x",
                rates.len(),
                low.0,
                low.1,
                high.0,
                high.1
            );
            // The ends of the band, for whoever reads a failure.
            let ends = rates.iter().take(5).chain(rates.iter().rev().take(5).rev());
            for (rate, row) in ends {
                println!("  {rate:8.1} ms per million  {row}");
            }
            if band > BAND {
                failures.push(format!(
                    "the band on {machine:?} is {band:.1}x, past {BAND}x"
                ));
            }
        }
    }

    // 5. What built exact within the target on the base revision still does.
    if let Some(base) = base {
        let (mut compared, mut exact) = (0, 0);
        let base = lines(base)?;
        for id in &expected {
            if !base.iter().any(|line| line["id"] == *id) {
                failures.push(format!("{id}: not in the base run"));
            }
        }
        for line in base {
            let id = line["id"].as_str().unwrap_or("?");
            let Some(head) = rows.get(id).map(|lines| lines[0]) else {
                failures.push(format!("{id}: in the base run and not in these runs"));
                continue;
            };
            compared += 1;
            let was = &line["build"];
            if was["status"] != "exact"
                || was["used"].as_u64().unwrap_or(u64::MAX) > TARGET_BUILD_WORK
            {
                continue;
            }
            exact += 1;
            let now = &head["build"];
            if now["status"] != "exact"
                || now["used"] != was["used"]
                || head["export"]["sha256"] != line["export"]["sha256"]
                || head["export"]["used"] != line["export"]["used"]
            {
                failures.push(format!(
                    "{id}: exact at {} on the base, now {} at {}",
                    was["used"], now["status"], now["used"]
                ));
            }
        }
        println!("against the base: {compared} rows compared, {exact} exact within the target");
    }

    if failures.is_empty() {
        println!("every line holds");
        Ok(())
    } else {
        for failure in &failures {
            println!("FAILED {failure}");
        }
        Err(format!("{} failures", failures.len()).into())
    }
}
