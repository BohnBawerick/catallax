//! Development-only work attribution for the eager application build.
//! Runs `catallax_app::build` on one recipe under a caller-owned allowance and
//! reports which decider calls spent its deterministic work. `work_profile`
//! measures the verification path instead; the two entries charge differently.
use catallax_app::DEFAULT_BUILD_WORK;
use catallax_mesh::MeshOptions;
use catallax_numbers::{
    Decider, Rational, WorkProfileRow, with_operation_budget, with_operation_log, with_work_profile,
};
use catallax_recipe::Recipe;
use serde_json::{Value, json};
use std::{collections::BTreeMap, env, fs};

const USAGE: &str = "usage: app_work_profile RECIPE [--max-work UNITS] [--chord TOLERANCE] \
[--exact-only]\n\
RECIPE is a catallax-recipe file. The default allowance is the application's own.\n\
--exact-only stops after the exact body and its topology check, without display.";

const TOP_CALLERS: usize = 30;

fn rows_by<'a>(
    rows: &'a [WorkProfileRow],
    key: impl Fn(&'a WorkProfileRow) -> String,
    used: u64,
) -> Vec<Value> {
    let mut grouped = BTreeMap::<String, (u64, u64)>::new();
    for row in rows {
        let entry = grouped.entry(key(row)).or_default();
        entry.0 += row.calls;
        entry.1 += row.work;
    }
    let mut grouped: Vec<_> = grouped.into_iter().collect();
    grouped.sort_by(|a, b| b.1.1.cmp(&a.1.1).then_with(|| a.0.cmp(&b.0)));
    grouped
        .into_iter()
        .map(|(key, (calls, work))| {
            let share = if used == 0 {
                0.0
            } else {
                work as f64 / used as f64
            };
            json!({"key":key,"calls":calls,"work":work,"share":share})
        })
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let Some(input) = args.next() else {
        eprintln!("{USAGE}");
        std::process::exit(64);
    };
    let mut max_work = DEFAULT_BUILD_WORK;
    let mut chord = "1e-2".to_owned();
    let mut exact_only = false;
    while let Some(flag) = args.next() {
        match (flag.as_str(), args.next()) {
            ("--exact-only", next) => {
                exact_only = true;
                if next.is_some() {
                    eprintln!("{USAGE}");
                    std::process::exit(64);
                }
            }
            ("--max-work", Some(value)) => max_work = value.parse()?,
            ("--chord", Some(value)) => chord = value,
            _ => {
                eprintln!("{USAGE}");
                std::process::exit(64);
            }
        }
    }
    let recipe = Recipe::parse(&fs::read_to_string(&input)?).map_err(|e| format!("{e:?}"))?;
    let options = MeshOptions {
        chord_tolerance: Rational::parse(&chord).map_err(|e| format!("{e:?}"))?,
        max_vertices: 200_000,
        ..Default::default()
    };
    let (((outcome, work), entries), rows) = with_work_profile(|| {
        with_operation_log(|| {
            with_operation_budget(max_work, || {
                let mut decider = Decider::default();
                if exact_only {
                    catallax_app::build_incremental_without_display(&recipe, None, &mut decider)
                        .map(|body| {
                            let built = body.checkpoint.built();
                            json!({"faces":built.model.faces_of_body(built.result).len()})
                        })
                } else {
                    catallax_app::build(&recipe, &options, &mut decider).map(|body| {
                        json!({
                            "faces":body.names.faces.len(),
                            "edges":body.names.edges.len(),
                            "vertices":body.names.vertices.len(),
                            "display_vertices":body.display.value.positions.len(),
                            "display_triangles":body.display.value.triangles.len(),
                        })
                    })
                }
            })
        })
    });
    let attributed: u64 = rows.iter().map(|row| row.work).sum();
    let status = match &outcome {
        _ if work.exhausted => json!({"status":"refused","code":"application.work-limit"}),
        Ok(counts) => json!({"status":"built","counts":counts}),
        Err(error) => json!({"status":"refused","code":error.code(),"detail":error.to_string()}),
    };
    let rows_json: Vec<_> = rows
        .iter()
        .map(|row| {
            json!({"section":row.section,"label":row.label,"caller":row.caller,
                "calls":row.calls,"work":row.work})
        })
        .collect();
    println!(
        "{}",
        json!({
            "input": input,
            "entry": if exact_only { "build_incremental_without_display" } else { "build" },
            "max_work": max_work,
            "chord_tolerance": chord,
            "outcome": status,
            "work_used": work.used,
            "exhausted": work.exhausted,
            "attributed": attributed,
            "unattributed": work.used.saturating_sub(attributed),
            "operations": entries.iter().map(|entry| json!({
                "id":entry.id,"operation":entry.operation,
                "work":entry.work_after - entry.work_before,
                "refusal":entry.refusal.as_ref().map(|refusal| refusal.code),
            })).collect::<Vec<_>>(),
            "sections": rows_by(&rows, |row| row.section.clone(), work.used),
            "labels": rows_by(&rows, |row| row.label.to_owned(), work.used),
            "callers": rows_by(&rows, |row| format!("{} {}", row.caller, row.label), work.used)
                .into_iter().take(TOP_CALLERS).collect::<Vec<_>>(),
            "rows": rows_json,
        })
    );
    Ok(())
}
