//! Report exact construction, native embedding and numerical admission per recipe.
//! Usage: cargo run -p catallax-features --example embedding_admission -- RECIPE...

use catallax_numbers::{
    BudgetPolicy, Decider, Rational, current_operation_work, with_operation_budget,
};
use catallax_recipe::Recipe;
use catallax_shapes::Curve;
use catallax_topology::numerical::{NumericalOptions, build_with_native_trims};
use serde_json::json;
use std::{error::Error, fs, process::ExitCode};

fn run(path: &str) -> Result<bool, Box<dyn Error>> {
    let recipe = Recipe::parse(&fs::read_to_string(path)?)?;
    let (passed, work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::new(BudgetPolicy {
            max_work: 100_000_000,
            ..Default::default()
        });
        let built = match catallax_features::evaluate(&recipe, &mut dec) {
            Ok(built) => built,
            Err(error) => {
                println!("{}", json!({"recipe": path, "stage": "construction", "error": error.to_string()}));
                return false;
            }
        };
        let construction_work = current_operation_work().used;
        let source_edges = built.model.edges_of_body(built.result);
        let tolerance = Rational::parse("0.0000001").unwrap();
        let (definition, native) = match catallax_features::numerical_embedding::definition_from_exact_model(
            &built.model, built.result, recipe.name.clone(), &tolerance, &mut dec,
        ) {
            Ok(definition) => definition,
            Err(error) => {
                println!("{}", json!({"recipe": path, "stage": "embedding", "error": error.to_string(), "construction_work": construction_work, "embedding_work": current_operation_work().used - construction_work}));
                return false;
            }
        };
        let embedding_end = current_operation_work().used;
        let edge_rows: Vec<_> = definition.edges.iter().zip(&source_edges).map(|(edge, source)| {
            let family = match built.model.geom().curve(built.model.edge(*source).curve) {
                Curve::Line(_) => "line",
                Curve::Circle(_) => "circle",
                Curve::BSpline(_) => "spline",
                _ => "traced",
            };
            json!({"edge": edge.key, "source_curve": family, "native_trims": native.iter().filter(|trim| trim.edge == edge.key).count()})
        }).collect();
        let body = match build_with_native_trims(
            definition, native, NumericalOptions { tolerance_mm: tolerance }, &mut dec,
        ) {
            Ok(body) => body,
            Err(error) => {
                println!("{}", json!({"recipe": path, "stage": "admission", "error": error.to_string(), "construction_work": construction_work, "embedding_work": embedding_end - construction_work, "admission_work": current_operation_work().used - embedding_end}));
                return false;
            }
        };
        println!("{}", json!({
            "recipe": path, "stage": "admitted", "tolerance_mm": "0.0000001",
            "faces": body.definition().faces.len(), "edges": body.definition().edges.len(),
            "construction_work": construction_work, "embedding_work": embedding_end - construction_work,
            "admission_work": current_operation_work().used - embedding_end,
        }));
        for mut row in edge_rows {
            row["recipe"] = json!(path);
            row["stage"] = json!("admitted");
            println!("{row}");
        }
        true
    });
    println!("{}", json!({"recipe": path, "work": work.used, "exhausted": work.exhausted}));
    Ok(passed && !work.exhausted)
}

fn main() -> ExitCode {
    let paths: Vec<_> = std::env::args().skip(1).collect();
    if paths.is_empty() {
        eprintln!("usage: embedding_admission RECIPE...");
        return ExitCode::FAILURE;
    }
    let mut passed = true;
    for path in paths {
        match run(&path) {
            Ok(result) => passed &= result,
            Err(error) => {
                eprintln!("{path}: {error}");
                passed = false;
            }
        }
    }
    if passed { ExitCode::SUCCESS } else { ExitCode::FAILURE }
}
