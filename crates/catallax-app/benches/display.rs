//! Display timing and complete mesh capture, outside the application build path.
use std::{env, fs, io::BufWriter, path::PathBuf, time::Instant};

use catallax_features::evaluate;
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::Decider;
use catallax_recipe::Recipe;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args: Vec<_> = env::args_os().skip(1).collect();
    let incremental = args.first().is_some_and(|arg| arg == "--incremental");
    if incremental {
        args.remove(0);
    }
    if args.len() != 3 {
        return Err(
            "usage: display [--incremental] INPUT.json OUTPUT_DIRECTORY REPETITIONS".into(),
        );
    }
    let input = fs::read_to_string(&args[0])?;
    let recipe =
        if serde_json::from_str::<serde_json::Value>(&input)?["format"] == "soliddraft-recipe/1" {
            catallax_app::soliddraft::convert(&input)?.recipe
        } else {
            Recipe::parse(&input)?
        };
    let out = PathBuf::from(&args[1]);
    fs::create_dir_all(&out)?;
    let repetitions: usize = args[2].to_str().ok_or("invalid repetitions")?.parse()?;
    if repetitions == 0 {
        return Err("repetitions must be positive".into());
    }
    let options = MeshOptions::default();
    let mut rows = Vec::new();
    if incremental {
        return incremental_edits(&recipe, &options, &out, repetitions);
    }
    for iteration in 0..repetitions {
        let mut dec = Decider::default();
        let start = Instant::now();
        let built = evaluate(&recipe, &mut dec)?;
        let faults = built.model.check(built.result);
        if !faults.is_empty() {
            return Err(format!("invalid body: {faults:?}").into());
        }
        let build_ms = start.elapsed().as_secs_f64() * 1000.0;
        let start = Instant::now();
        let mesh = mesh_body(&built.model, built.result, &options, &mut dec)?;
        let display_ms = start.elapsed().as_secs_f64() * 1000.0;
        rows.push(serde_json::json!({
            "iteration": iteration, "build_ms": build_ms, "display_ms": display_ms,
            "vertices": mesh.value.positions.len(), "triangles": mesh.value.triangles.len(),
            "chord_tolerance_mm": options.chord_tolerance.to_string(),
            "max_vertices": options.max_vertices,
        }));
        if iteration + 1 == repetitions {
            mesh.value.write_obj(&mut BufWriter::new(fs::File::create(
                out.join("display.obj"),
            )?))?;
            fs::write(out.join("mesh.txt"), format!("{mesh:?}\n"))?;
        }
    }
    println!("{}", serde_json::to_string(&rows)?);
    fs::write(out.join("samples.json"), serde_json::to_vec_pretty(&rows)?)?;
    Ok(())
}

fn incremental_edits(
    recipe: &Recipe,
    options: &MeshOptions,
    out: &std::path::Path,
    repetitions: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut rows = Vec::new();
    for iteration in 0..repetitions {
        let mut checkpoint = None;
        for (index, step) in recipe.steps.iter().enumerate() {
            if step.op.kind() != catallax_recipe::ValueKind::Body {
                continue;
            }
            let mut edited = recipe.clone();
            edited.steps.truncate(index + 1);
            edited.result = step.id.clone();
            let start = Instant::now();
            let body = catallax_app::build_incremental(
                &edited,
                checkpoint,
                options,
                &mut Decider::default(),
            )?;
            let elapsed = start.elapsed().as_secs_f64() * 1000.0;
            let built = body.checkpoint.built();
            assert_eq!(
                body.names.faces.len(),
                built.model.counts(built.result).faces
            );
            assert!(
                body.display
                    .value
                    .triangles
                    .iter()
                    .all(|triangle| body.names.faces.contains_key(&triangle.face))
            );
            rows.push(serde_json::json!({
                "iteration":iteration,"step":step.id.as_str(),"edit_ms":elapsed,
                "faces":body.names.faces.len(),"triangles":body.display.value.triangles.len(),
            }));
            if iteration + 1 == repetitions && index + 1 == recipe.steps.len() {
                body.display
                    .value
                    .write_obj(&mut BufWriter::new(fs::File::create(
                        out.join("display.obj"),
                    )?))?;
            }
            checkpoint = Some(body.checkpoint);
        }
    }
    println!("{}", serde_json::to_string(&rows)?);
    fs::write(out.join("samples.json"), serde_json::to_vec_pretty(&rows)?)?;
    Ok(())
}
