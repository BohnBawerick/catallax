//! Non-gating measurement of SolidDraft's rounded STEP handoff.
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_step::{StepOptions, write_step};
use serde_json::{Value, json};
use std::{env, fs, path::PathBuf, time::Instant};

fn main() {
    let out = PathBuf::from(env::args_os().nth(1).expect("output directory"));
    fs::create_dir_all(&out).unwrap();
    let mut source: Value =
        serde_json::from_str(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    source["steps"][3]["radius"] = json!("2.5");
    let recipe = Recipe::parse(&source.to_string()).unwrap();
    let options = StepOptions {
        product_name: "Part1".into(),
        ..StepOptions::default()
    };
    let ((answer, report), work) = with_operation_budget(20_000_000, || {
        let body =
            catallax_app::build_incremental_without_display(&recipe, None, &mut Decider::default())
                .unwrap();
        let built = body.checkpoint.built();
        let mut bytes = Vec::new();
        eprintln!("starting rounded STEP writer");
        let start = Instant::now();
        let answer = write_step(
            &built.model,
            &[built.result],
            &options,
            &mut Decider::default(),
            &mut bytes,
        )
        .unwrap();
        let elapsed = start.elapsed();
        assert!(answer.marks.is_empty(), "{:#?}", answer.marks);
        assert!(answer.value.max_deviation <= options.curve_tolerance);
        fs::write(out.join("answer.step"), bytes).unwrap();
        (
            format!("{answer:#?}"),
            json!({"writer_ms": elapsed.as_secs_f64() * 1000.0}),
        )
    });
    assert!(!work.exhausted);
    fs::write(out.join("answer.txt"), answer).unwrap();
    fs::write(out.join("work.txt"), work.used.to_string()).unwrap();
    println!("{report}");
}
