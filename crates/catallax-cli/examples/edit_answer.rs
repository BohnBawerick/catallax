//! Exact answer and certified-display capture for the hosted edit comparison.
use catallax_numbers::{Decider, with_operation_budget, with_operation_log};
use catallax_recipe::Recipe;
use catallax_step::{StepOptions, write_step};
use serde_json::{Value, json};
use std::{env, fs, path::PathBuf, time::Instant};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = PathBuf::from(env::args_os().nth(1).ok_or("native recipe path")?);
    let out = PathBuf::from(env::args_os().nth(2).ok_or("output directory")?);
    fs::create_dir_all(&out)?;
    let source: Value = serde_json::from_slice(&fs::read(input)?)?;
    let recipe = Recipe::parse(&source.to_string())?;
    let ((answer, mut operations), work) = with_operation_budget(20_000_000, || {
        with_operation_log(|| -> Result<Value, Box<dyn std::error::Error>> {
            let mut dec = Decider::default();
            let start = Instant::now();
            let body = catallax_app::build_incremental_without_display(&recipe, None, &mut dec)?;
            let build_ms = start.elapsed().as_secs_f64() * 1000.0;
            let start = Instant::now();
            let mesh = body.display(&Default::default(), &mut dec)?;
            let display_ms = start.elapsed().as_secs_f64() * 1000.0;
            let built = body.checkpoint.built();
            let mut obj = Vec::new();
            mesh.value.write_obj(&mut obj)?;
            fs::write(out.join("display.obj"), obj)?;
            fs::write(out.join("mesh.txt"), format!("{mesh:#?}\n"))?;
            fs::write(out.join("names.txt"), format!("{:#?}\n", body.names))?;
            let mut bytes = Vec::new();
            let step = write_step(
                &built.model,
                &[built.result],
                &StepOptions::default(),
                &mut dec,
                &mut bytes,
            )?;
            fs::write(out.join("answer.step"), bytes)?;
            fs::write(out.join("answer.txt"), format!("{step:#?}\n"))?;
            fs::write(out.join("stats.txt"), format!("{:#?}\n", dec.stats()))?;
            Ok(json!({"build_ms":build_ms,"display_ms":display_ms}))
        })
    });
    let report = answer?;
    assert!(!work.exhausted);
    fs::write(out.join("work.txt"), format!("{}\n", work.used))?;
    let timings: Vec<_> = operations
        .iter()
        .map(|operation| json!({"id":operation.id,"operation":operation.operation,"wall_ms":operation.wall_ms}))
        .collect();
    // Compare operation outcomes, not dispatch counts or cumulative work prefixes.
    // Actual total work and final decision counters remain in their own files.
    for operation in &mut operations {
        operation.wall_ms = 0.0;
        operation.work_before = 0;
        operation.work_after = 0;
        operation.before = Default::default();
        if let Some(after) = &mut operation.after {
            *after = Default::default();
        }
    }
    fs::write(out.join("operations.txt"), format!("{operations:#?}\n"))?;
    println!("{}", json!({"timing":report,"operations":timings}));
    Ok(())
}
