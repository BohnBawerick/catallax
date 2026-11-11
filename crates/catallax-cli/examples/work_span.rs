//! Development-only work and span measurement. Never linked into an application entry.
//! Answers one verification request exactly as `catallax-verify replay` does and
//! reports its work, its critical-path work (span) on an ideal machine with
//! unlimited workers, and where the span lies.
use catallax_cli::verification::{DEFAULT_WORK, ExportedQuestion, answer};
use catallax_numbers::work_span::{MODELS, Model, WorkSpanReport, with_work_span};
use serde_json::{Value, json};
use std::{collections::BTreeMap, env, fs};

const USAGE: &str = "usage: work_span QUESTION_OR_RECIPE --out DIR [--max-work UNITS]\n\
QUESTION_OR_RECIPE is a catallax-question/1 file, such as a driver case's question.json,\n\
or a bare recipe answered with source `recipe` and the default allowance.";

const TOP: usize = 25;

// Arena kinds are full type paths; keep the type name.
fn short(kind: &str) -> String {
    match kind.rsplit_once("::") {
        Some((_, name)) => format!("arena {name}"),
        None => kind.to_owned(),
    }
}

fn ratio(work: u64, span: u64) -> f64 {
    if span == 0 {
        1.0
    } else {
        work as f64 / span as f64
    }
}

// The build is the recipe step loop; everything after it in the request
// (summary, STEP export) is listed separately.
const BUILD: &str = "request > recipe steps";

fn model_name(model: Option<Model>) -> &'static str {
    match model {
        Some(Model::Steps) => "steps",
        Some(Model::Independent) => "independent",
        Some(Model::Speculative) => "speculative",
        Some(Model::Unshared) => "unshared",
        Some(Model::Ordered) => "ordered",
        None => "unknown",
    }
}

// Boolean stages inside the build, then what the request does after it.
const STAGES: [&str; 4] = ["Intersect", "Imprint", "Classify", "Sew"];

// Work and span by Boolean stage, summed over every Boolean of the build,
// and for the measurement and export that follow it. Work is inclusive of
// nested regions; a region's span is its own time open on its strand.
fn stages(report: &WorkSpanReport) -> Value {
    let inclusive = |prefix: &str| -> u64 {
        report
            .regions
            .iter()
            .filter(|row| row.region == prefix || row.region.starts_with(&format!("{prefix} > ")))
            .map(|row| row.work)
            .sum()
    };
    let mut rows = serde_json::Map::new();
    let mut add = |name: &str, regions: Vec<&str>| {
        let work: u64 = regions.iter().map(|region| inclusive(region)).sum();
        let mut span = [0u64; MODELS.len()];
        for region in &regions {
            if let Some(row) = report.regions.iter().find(|row| row.region == *region) {
                for (total, value) in span.iter_mut().zip(row.span) {
                    *total += value;
                }
            }
        }
        rows.insert(name.to_owned(), json!({"work": work, "span": span}));
    };
    for stage in STAGES {
        let regions = report
            .regions
            .iter()
            .filter(|row| {
                row.region.starts_with(BUILD) && row.region.ends_with(&format!(" > {stage}"))
            })
            .map(|row| row.region.as_str())
            .collect();
        add(stage, regions);
    }
    let after = |kind: &str| -> Vec<&str> {
        report
            .regions
            .iter()
            .filter(|row| {
                row.region
                    .strip_prefix("request > ")
                    .is_some_and(|rest| !rest.contains(" > ") && rest.starts_with(kind))
            })
            .map(|row| row.region.as_str())
            .collect()
    };
    add("measurement", after("brep summarize"));
    add("export", after("brep step"));
    let direct = report
        .regions
        .iter()
        .find(|row| row.region == "request")
        .map_or(0, |row| row.work);
    rows.insert(
        "request, outside any operation".into(),
        json!({"work": direct}),
    );
    Value::Object(rows)
}

fn summary(report: &WorkSpanReport) -> Value {
    let build_work: u64 = report
        .regions
        .iter()
        .filter(|row| row.region.starts_with(BUILD))
        .map(|row| row.work)
        .sum();
    let build_loop = report.loops.iter().find(|row| row.region == BUILD);
    let spans: Vec<_> = report
        .spans
        .iter()
        .enumerate()
        .map(|(index, row)| {
            let mut by_region = BTreeMap::<&str, u64>::new();
            for segment in &row.path {
                *by_region.entry(segment.region.as_str()).or_default() += segment.work;
            }
            let mut path: Vec<_> = by_region.into_iter().collect();
            path.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));
            let build_span = build_loop.map(|row| row.span[index]);
            json!({
                "model": model_name(row.model),
                "span": row.span,
                "parallelism": ratio(report.work, row.span),
                "build_span": build_span,
                "build_parallelism": build_span.map(|span| ratio(build_work, span)),
                "path": path.iter().take(TOP).map(|(region, work)| json!({
                    "region": region, "work": work, "share_of_span": ratio(*work, row.span),
                })).collect::<Vec<_>>(),
            })
        })
        .collect();
    let mut regions: Vec<_> = report.regions.iter().collect();
    regions.sort_by(|a, b| b.work.cmp(&a.work).then_with(|| a.region.cmp(&b.region)));
    let mut loops: Vec<_> = report.loops.iter().collect();
    loops.sort_by(|a, b| b.work.cmp(&a.work).then_with(|| a.region.cmp(&b.region)));
    let mut dependencies: Vec<_> = report.dependencies.iter().collect();
    dependencies.sort_by(|a, b| b.delay.cmp(&a.delay).then_with(|| a.kind.cmp(&b.kind)));
    json!({
        "work": report.work,
        "build": build_loop.map(|row| json!({
            "work": build_work, "steps": row.tasks, "largest_step": row.largest_task,
        })),
        "spans": spans,
        "largest_decision": report.decisions.first().map(|row| json!({
            "label": row.label, "caller": row.caller, "region": row.region, "work": row.largest,
            "work_over_largest": ratio(report.work, row.largest),
        })),
        "stages": stages(report),
        "regions": regions.iter().take(TOP).map(|row| json!({
            "region": row.region, "work": row.work, "span": row.span,
        })).collect::<Vec<_>>(),
        "loops": loops.iter().take(TOP).map(|row| json!({
            "region": row.region, "search": row.search, "runs": row.runs, "tasks": row.tasks,
            "work": row.work, "span": row.span,
            "parallelism": row.span.map(|span| ratio(row.work, span)),
            "largest_task": row.largest_task,
        })).collect::<Vec<_>>(),
        "dependencies": dependencies.iter().map(|row| json!({
            "kind": short(&row.kind), "shared": row.shared, "waits": row.waits,
            "delay": row.delay, "sibling_waits": row.sibling_waits,
        })).collect::<Vec<_>>(),
        "decisions": report.decisions.iter().take(10).map(|row| json!({
            "label": row.label, "caller": row.caller, "region": row.region,
            "calls": row.calls, "work": row.work, "largest": row.largest,
        })).collect::<Vec<_>>(),
        "tracked": {"reads": report.reads, "writes": report.writes},
    })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args().skip(1).collect();
    let (input, out, max_work) = match args.as_slice() {
        [input, flag, out] if flag == "--out" => (input, out, None),
        [input, flag, out, work_flag, work] if flag == "--out" && work_flag == "--max-work" => {
            (input, out, Some(work.parse::<u64>()?))
        }
        _ => {
            eprintln!("{USAGE}");
            std::process::exit(64);
        }
    };
    let problem: Value = serde_json::from_slice(&fs::read(input)?)?;
    let question = if problem["format"] == "catallax-question/1" {
        serde_json::from_value(problem)?
    } else {
        ExportedQuestion {
            format: "catallax-question/1".into(),
            source: "recipe".into(),
            problem,
            max_work: DEFAULT_WORK,
        }
    };
    let mut request = question.request(out.into())?;
    if let Some(max_work) = max_work {
        request.max_work = max_work;
    }
    let (response, report) = with_work_span(|| answer(&request));
    let response = response?;
    let used = response["work"]["used"].as_u64().unwrap_or_default();
    let paths_match = report
        .spans
        .iter()
        .all(|row| row.path.iter().map(|segment| segment.work).sum::<u64>() == row.span);
    println!(
        "{}",
        json!({
            "input": input,
            "max_work": request.max_work,
            "status": response["manifest"]["status"],
            "refusal": response["manifest"]["refusal"]["code"],
            "exhausted": response["work"]["exhausted"],
            "work_used": used,
            "observed_work_matches": used == report.work,
            "paths_match_spans": paths_match,
            "summary": summary(&report),
        })
    );
    Ok(())
}
