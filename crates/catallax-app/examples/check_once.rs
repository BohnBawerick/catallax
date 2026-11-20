//! Check once (roadmap 3.0A item 5): the work and the time of the numerical
//! STEP export of the owner's pocket rim at R2, and of reopening the written
//! file. Export and reopening are measured apart.
//!
//! Usage: check_once [--runs N] [--write OUT.step] [REOPEN.step ...]
//!
//! One JSON line per measurement. Work is deterministic; the clock is not,
//! so each measurement runs N times (3 by default) and prints every run with
//! the machine's load average beside it. Reopening is measured twice: under
//! a work allowance large enough to count all of it, and as the application
//! opens a file, with no allowance of its own.
use std::time::Instant;

use catallax_app::{EvaluatedApplicationOutcome, build_incremental_with_files_without_display};
use catallax_import::{
    sha256_hex,
    tolerant::opened::{BodyView, open_step},
};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::{Recipe, StepId};
use serde_json::{Value, json};

const OWNER: &str = include_str!("../../catallax-features/tests/fixtures/owner-curved-pocket.json");

const WHOLE_RIM: [(&str, &str); 6] = [
    ("top", "near"),
    ("top", "left"),
    ("top", "far"),
    ("right-round", "near"),
    ("right-round", "right"),
    ("right-round", "far"),
];

/// An allowance no measured case reaches, so the count is complete.
const COUNT_ALL: u64 = 1_000_000_000_000;

fn rim(radius: &str) -> Recipe {
    let mut value: Value = serde_json::from_str(OWNER).unwrap();
    let edges: Vec<Value> = WHOLE_RIM
        .iter()
        .map(|(outer, wall)| {
            json!({"between": [
                {"step": "blank", "role": {"side": {"segment": outer}}},
                {"step": "pocket-tool", "role": {"side": {"segment": wall}}}
            ]})
        })
        .collect();
    value["steps"].as_array_mut().unwrap().push(json!({
        "id": "round", "op": "fillet", "body": "pocket", "radius": radius, "edges": edges
    }));
    value["result"] = json!("round");
    Recipe::parse(&value.to_string()).unwrap()
}

fn load() -> String {
    std::fs::read_to_string("/proc/loadavg")
        .map(|text| text.split(' ').take(3).collect::<Vec<_>>().join(" "))
        .unwrap_or_default()
}

fn reopen(label: &str, bytes: &[u8], run: usize, counted: bool) {
    let sha = sha256_hex(bytes);
    let step = StepId::new("import").unwrap();
    let before = load();
    let start = Instant::now();
    let open = || open_step(bytes, &sha, &step, Default::default(), &mut Decider::default());
    let (opened, work) = if counted {
        let (opened, work) = with_operation_budget(COUNT_ALL, open);
        (opened, Some(work))
    } else {
        (open(), None)
    };
    let wall_ms = start.elapsed().as_secs_f64() * 1e3;
    let opened = match opened {
        Ok(opened) => opened,
        Err(error) => {
            println!(
                "{}",
                json!({"measure":"reopen","file":label,"run":run,"counted":counted,
                    "refused":format!("{error:?}"),"wall_ms":wall_ms})
            );
            return;
        }
    };
    let bodies: Vec<Value> = (0..opened.bodies().len())
        .map(|index| {
            let view = match opened.body(index) {
                Some(BodyView::Certified { .. }) => "certified",
                Some(BodyView::Numerical(_)) => "numerical",
                Some(BodyView::Tolerant { .. }) => "tolerant",
                Some(BodyView::Surfaces { .. }) => "surfaces",
                None => "none",
            };
            json!({"tier": opened.bodies()[index].tier(), "view": view,
                "refusal": opened.bodies()[index].refusal().map(|e| format!("{e:?}")),
                "upgrade_refusal": opened.bodies()[index].upgrade_refusal().map(|e| format!("{e:?}"))})
        })
        .collect();
    println!(
        "{}",
        json!({"measure":"reopen","file":label,"run":run,"counted":counted,
            "bytes":bytes.len(),"sha256":sha,
            "work":work.map(|w| w.used),"exhausted":work.map(|w| w.exhausted),
            "read_work":opened.read_work_used(),
            "file_admission_work":opened.file_admission_work_used(),
            "native_admission_work":opened.native_admission_work_used(),
            "checked_upgrade_work":opened.checked_upgrade_work_used(),
            "bodies":bodies,"wall_ms":wall_ms,"load_before":before,"load_after":load()})
    );
}

fn main() {
    let mut runs = 3usize;
    let mut write = None;
    let mut files = Vec::new();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--runs" => runs = args.next().and_then(|n| n.parse().ok()).expect("--runs N"),
            "--write" => write = Some(args.next().expect("--write OUT.step")),
            _ => files.push(arg),
        }
    }
    let recipe = rim("2");
    let start = Instant::now();
    let (outcome, build) = with_operation_budget(100_000_000, || {
        build_incremental_with_files_without_display(
            &recipe,
            None,
            &Default::default(),
            &mut Decider::default(),
        )
    });
    let EvaluatedApplicationOutcome::Numerical(result) = outcome.unwrap() else {
        panic!("the exact builders refuse this rim");
    };
    println!(
        "{}",
        json!({"measure":"build","part":"owner pocket rim R2","work":build.used,
            "exhausted":build.exhausted,"wall_ms":start.elapsed().as_secs_f64() * 1e3,
            "faces":result.body().definition().faces.len()})
    );
    let options = catallax_step::numerical::NumericalStepOptions {
        product_name: "Owner pocket R2".into(),
        ..Default::default()
    };
    let mut written = Vec::new();
    for run in 0..runs {
        let before = load();
        let start = Instant::now();
        let (step, export) = with_operation_budget(100_000_000, || {
            catallax_step::numerical::build(result.body(), &options, &mut Decider::default())
        });
        let wall_ms = start.elapsed().as_secs_f64() * 1e3;
        let step = step.unwrap();
        println!(
            "{}",
            json!({"measure":"export","part":"owner pocket rim R2","run":run,
                "work":export.used,"exhausted":export.exhausted,"wall_ms":wall_ms,
                "bytes":step.bytes().len(),"sha256":sha256_hex(step.bytes()),
                "load_before":before,"load_after":load()})
        );
        written = step.bytes().to_vec();
    }
    if let Some(path) = write {
        std::fs::write(&path, &written).unwrap();
    }
    let mut sources = vec![("written by this build".to_owned(), written)];
    for path in files {
        let bytes = std::fs::read(&path).unwrap();
        sources.push((path, bytes));
    }
    for (label, bytes) in &sources {
        for run in 0..runs {
            reopen(label, bytes, run, true);
        }
        for run in 0..runs {
            reopen(label, bytes, run, false);
        }
    }
}
