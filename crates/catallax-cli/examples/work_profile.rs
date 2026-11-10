//! Development-only request work attribution. Never linked into an application entry.
//! Answers one verification request exactly as `catallax-verify replay` does and
//! reports which decision kinds spent its deterministic work.
use catallax_cli::verification::{DEFAULT_WORK, ExportedQuestion, answer};
use catallax_numbers::{WorkProfileRow, with_work_profile};
use serde_json::{Value, json};
use std::{collections::BTreeMap, env, fs};

const USAGE: &str = "usage: work_profile QUESTION_OR_RECIPE --out DIR [--max-work UNITS]\n\
QUESTION_OR_RECIPE is a catallax-question/1 file, such as a driver case's question.json,\n\
or a bare recipe answered with source `recipe` and the default allowance.";

const TOP_CALLERS: usize = 20;

// Coarse classes for the stage 2 cause inventory. Exact-zero signs and
// enclosures that never separated zero bound the work a 1 nm decision could
// replace; every other class remains under any tolerance rule.
fn class(label: &str) -> &'static str {
    match label {
        "sign.zero.rational"
        | "sign.zero.algebraic"
        | "sign.zero.separation"
        | "sign.zero.fact" => "exact-zero sign",
        "sign.filter" | "sign.refined" | "sign.exact-nonzero" => "nonzero sign",
        "sign.not-exact" => "not-exact sign",
        "sign.refused" => "refused sign",
        "real_roots" | "certify_root" => "root isolation",
        "real_poly_eliminant"
        | "real_poly_resultant"
        | "real_poly_tower_split"
        | "tower_root_system" => "elimination",
        "enclose.contains-zero" => "unseparated enclosure",
        "enclose.separated"
        | "enclose.refused"
        | "enclose_at"
        | "enclose_all_at"
        | "enclose_budgeted_at" => "enclosure",
        "rational_identity" | "rational_value" | "angle_degrees" => "exact reduction",
        "check_deadline" => "checkpoint",
        _ => "other",
    }
}

fn summary(rows: &[WorkProfileRow], used: u64) -> Value {
    let mut classes = BTreeMap::<&str, (u64, u64)>::new();
    let mut sections = BTreeMap::<&str, BTreeMap<&str, u64>>::new();
    for row in rows {
        let class = class(row.label);
        let entry = classes.entry(class).or_default();
        entry.0 += row.calls;
        entry.1 += row.work;
        *sections
            .entry(row.section.as_str())
            .or_default()
            .entry(class)
            .or_default() += row.work;
    }
    let attributed: u64 = rows.iter().map(|row| row.work).sum();
    let mut callers = BTreeMap::<(&str, &str), u64>::new();
    for row in rows {
        *callers
            .entry((row.caller.as_str(), class(row.label)))
            .or_default() += row.work;
    }
    let mut callers: Vec<_> = callers.into_iter().collect();
    callers.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let zero: Vec<_> = rows
        .iter()
        .filter(|row| class(row.label) == "exact-zero sign")
        .collect();
    let zero_work: u64 = zero.iter().map(|row| row.work).sum();
    let first: u64 = zero.iter().map(|row| row.first_enclosure_work).sum();
    let unseparated = classes
        .get("unseparated enclosure")
        .map_or(0, |(_, work)| *work);
    let share = |work: u64| {
        if used == 0 {
            0.0
        } else {
            work as f64 / used as f64
        }
    };
    json!({
        "work_used": used,
        "attributed": attributed,
        "unattributed": used.saturating_sub(attributed),
        "classes": classes.iter().map(|(class, (calls, work))| {
            (class.to_string(), json!({"calls":calls,"work":work,"share":share(*work)}))
        }).collect::<serde_json::Map<_, _>>(),
        "exact_zero": {
            "calls": zero.iter().map(|row| row.calls).sum::<u64>(),
            "work": zero_work,
            "share": share(zero_work),
            "first_enclosure_work": first,
            "beyond_first_enclosure": zero_work.saturating_sub(first),
            "beyond_share": share(zero_work.saturating_sub(first)),
        },
        "tolerance_ceiling": {
            "work": zero_work + unseparated,
            "share": share(zero_work + unseparated),
        },
        "sections": sections,
        "top_callers": callers.iter().take(TOP_CALLERS).map(|((caller, class), work)| {
            json!({"caller":caller,"class":class,"work":work,"share":share(*work)})
        }).collect::<Vec<_>>(),
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
    let (response, rows) = with_work_profile(|| answer(&request));
    let response = response?;
    let used = response["work"]["used"].as_u64().unwrap_or_default();
    let rows_json: Vec<_> = rows
        .iter()
        .map(|row| {
            json!({"section":row.section,"label":row.label,"caller":row.caller,
                "calls":row.calls,"work":row.work,
                "first_enclosure_work":row.first_enclosure_work})
        })
        .collect();
    println!(
        "{}",
        json!({
            "input": input,
            "max_work": request.max_work,
            "status": response["manifest"]["status"],
            "refusal": response["manifest"]["refusal"]["code"],
            "exhausted": response["work"]["exhausted"],
            "summary": summary(&rows, used),
            "rows": rows_json,
        })
    );
    Ok(())
}
