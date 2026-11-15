//! Measures body adoption separately from construction on a saved recipe/question.

use catallax_app::handover::{read_body, write_body};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;
use std::{hint::black_box, time::Instant};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let text = match std::env::args().nth(1) {
        Some(path) => std::fs::read_to_string(path)?,
        None => include_str!("../../catallax-cli/tests/fixtures/soliddraft-gui-bracket.json")
            .to_string(),
    };
    let document: serde_json::Value = serde_json::from_str(&text)?;
    let problem = document.get("problem").unwrap_or(&document);
    let recipe = Recipe::parse(&problem.to_string())?;
    let source = catallax_features::evaluate(&recipe, &mut Decider::default())?;
    let export_start = Instant::now();
    let (bytes, export_work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
        write_body(&source, &mut Decider::default())
    });
    let bytes = bytes?;
    let export_us = export_start.elapsed().as_micros();
    let mut imports_us = Vec::new();
    let mut import_work = Vec::new();
    for _ in 0..7 {
        let start = Instant::now();
        let (body, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
            read_body(&bytes, &mut Decider::default())
        });
        let body = body?;
        imports_us.push(start.elapsed().as_micros());
        import_work.push(work.used);
        black_box(body.model.counts(body.result));
    }
    println!(
        "{}",
        serde_json::json!({
            "part": recipe.name, "bytes": bytes.len(), "imports_us": imports_us,
            "import_work": import_work,
            "export_us": export_us, "export_work": export_work.used,
        })
    );
    Ok(())
}
