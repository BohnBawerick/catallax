//! Development-only sampled observer shared by all three B-rep producers.
use catallax_field::{FaceState, Plan, Workspace};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{RPoint3, distance::Options};
use serde_json::{Value, json};
use std::{env, fs};

fn observe(input: &Value) -> Result<Value, Box<dyn std::error::Error>> {
    let recipe = Recipe::parse(&serde_json::to_string(&input["recipe"])?)?;
    let mut dec = Decider::default();
    let plan = Plan::compile_checked(&recipe, &mut dec)?;
    if input["action"] == "samples" {
        let mut points = Vec::new();
        let mut unresolved = 0usize;
        let mut workspace = Workspace::default();
        let mut failure = None;
        let completed = plan.visit_candidates(4, |point, face| {
            match plan.face_query_checked(face, &point, &mut workspace, &mut dec) {
                Ok(query) => match query.state {
                    FaceState::OnFace => points.push(json!([
                        point.x.to_string(), point.y.to_string(), point.z.to_string()
                    ])),
                    FaceState::NotOnFace => (),
                    _ => unresolved += 1,
                },
                Err(error) => {
                    failure = Some(error.to_string());
                    return false;
                }
            }
            true
        })?;
        return Ok(json!({"status":"answered","points":points,"completed":completed,
            "unresolved_candidates":unresolved,"failure":failure,
            "scope":"trim-checked definition-native rational candidates, four divisions; finite sampling, not complete boundary coverage"}));
    }
    let options = Options {
        width: Rational::parse(input["distance_width_mm"].as_str().ok_or("width required")?)?,
        max_cells: 256,
    };
    let prepared = plan.prepare_distances(&mut dec)?;
    let mut distances = Vec::new();
    for coordinates in input["points"].as_array().ok_or("points required")? {
        let values = coordinates.as_array().ok_or("point array required")?;
        if values.len() != 3 { return Err("three coordinates required".into()); }
        let coordinate = |index: usize| -> Result<Rational, Box<dyn std::error::Error>> {
            Ok(Rational::parse(values[index].as_str().ok_or("exact coordinate required")?)?)
        };
        let point = RPoint3::new(coordinate(0)?, coordinate(1)?, coordinate(2)?);
        let distance = prepared.point_to_boundary(&point.to_real(), &options, &mut dec)?;
        distances.push(json!({"lower_mm":distance.lower.to_string(),
            "upper_mm":distance.upper.map(|q|q.to_string()),"converged":distance.converged}));
    }
    Ok(json!({"status":"answered","distances":distances,
        "scope":"certified point-to-definition distance intervals at supplied samples; not a global distance certificate"}))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args_os().skip(1).collect();
    if args.len() != 2 { return Err("usage: kernel_compare INPUT.json OUTPUT.json".into()); }
    let input: Value = serde_json::from_slice(&fs::read(&args[0])?)?;
    let ((result, work), elapsed) = {
        let start = std::time::Instant::now();
        let result = with_operation_budget(100_000_000, || observe(&input));
        (result, start.elapsed())
    };
    let mut report = match result {
        Ok(value) if !work.exhausted => value,
        Ok(_) => json!({"status":"unresolved","detail":"observer work allowance exhausted"}),
        Err(error) => json!({"status":"unresolved","detail":error.to_string()}),
    };
    report["work"] = json!(work.used);
    report["wall_ms"] = json!(elapsed.as_secs_f64()*1000.0);
    fs::write(&args[1], serde_json::to_vec_pretty(&report)?)?;
    Ok(())
}
