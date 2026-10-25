//! Display timing and complete mesh capture, outside the application build path.
use std::{env, fs, io::BufWriter, path::PathBuf, time::Instant};

use catallax_features::evaluate;
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::Decider;
use catallax_recipe::Recipe;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args: Vec<_> = env::args_os().skip(1).collect();
    let previous_input = if args.first().is_some_and(|arg| arg == "--edit-from") {
        if args.len() < 2 {
            return Err("--edit-from requires a native recipe path".into());
        }
        args.remove(0);
        Some(args.remove(0))
    } else {
        None
    };
    let incremental = args.first().is_some_and(|arg| arg == "--incremental");
    let application = args.first().is_some_and(|arg| arg == "--application");
    if incremental || application {
        args.remove(0);
    }
    if !(3..=5).contains(&args.len()) {
        return Err(
            "usage: display [--incremental|--application|--edit-from OLD_NATIVE.json] INPUT.json OUTPUT_DIRECTORY REPETITIONS [CHORD_MM|drawing [ANGLE_RADIANS]]".into(),
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
    let drawing = args.get(3).is_some_and(|value| value == "drawing");
    if (incremental || application || previous_input.is_some()) && drawing {
        return Err("drawing chord mode requires independent display measurements".into());
    }
    let options = MeshOptions {
        chord_tolerance: args
            .get(3)
            .filter(|_| !drawing)
            .map(|value| value.to_string_lossy().parse())
            .transpose()?
            .unwrap_or_else(|| MeshOptions::default().chord_tolerance),
        angular_deviation: args
            .get(4)
            .map(|value| value.to_string_lossy().parse())
            .transpose()?,
        ..MeshOptions::default()
    };
    if let Some(previous_input) = previous_input {
        let previous = Recipe::parse(&fs::read_to_string(previous_input)?)?;
        return edit_rebuilds(&previous, &recipe, &options, &out, repetitions);
    }
    let mut rows = Vec::new();
    if incremental {
        return incremental_edits(&recipe, &options, &out, repetitions);
    }
    if application {
        return application_builds(&recipe, &options, &out, repetitions);
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
        let mut options = options.clone();
        if drawing {
            let bounds = catallax_measure::bounding_box(
                &built.model,
                built.result,
                &catallax_measure::Accuracy {
                    relative: "1e-8".parse()?,
                    absolute: "1e-10".parse()?,
                },
                &mut dec,
            )?
            .value
            .ok_or("empty drawing body")?;
            let diagonal = (0..3)
                .map(|i| (&bounds.max[i] - &bounds.min[i]).to_f64_nearest().powi(2))
                .sum::<f64>()
                .sqrt()
                .max(1.0);
            options.chord_tolerance = catallax_numbers::Rational::from_f64_exact(diagonal * 1e-5)
                .ok_or("non-finite drawing diagonal")?;
        }
        // The unproved preview is timed on its own and excluded from the
        // certified wait, which keeps its historical definition.
        let preview_start = Instant::now();
        let preview = catallax_app::preview_body(&built, &Default::default(), &mut dec)?;
        let preview_ms = preview_start.elapsed().as_secs_f64() * 1000.0;
        let mesh_start = Instant::now();
        let mesh = mesh_body(&built.model, built.result, &options, &mut dec)?;
        let display_ms = mesh_start.elapsed().as_secs_f64() * 1000.0;
        let requested_wait_ms = start.elapsed().as_secs_f64() * 1000.0 - preview_ms;
        rows.push(serde_json::json!({
            "iteration": iteration, "build_ms": build_ms, "display_ms": display_ms,
            "requested_wait_ms": requested_wait_ms,
            "preview_ms": preview_ms, "preview_wait_ms": build_ms + preview_ms,
            "preview_state": "unproved", "preview_target_ms": 1000,
            "preview_triangles": preview.mesh.triangles.len(),
            "vertices": mesh.value.positions.len(), "triangles": mesh.value.triangles.len(),
            "chord_tolerance_mm": options.chord_tolerance.to_string(),
            "angular_deviation_radians": options.angular_deviation.as_ref().map(ToString::to_string),
            "boundary_segments": mesh.value.edges.iter().map(|edge| edge.vertices.len() - 1).sum::<usize>(),
            "faces": mesh.value.faces.iter().map(|face| serde_json::json!({
                "face": face.face.index(), "triangles": face.triangles.len(),
                "chord_bound_mm": face.max_chord_error.to_f64_nearest(),
                "angular_bound_radians": face.max_angular_deviation.as_ref().map(|v| v.to_f64_nearest()),
            })).collect::<Vec<_>>(),
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

fn application_builds(
    recipe: &Recipe,
    options: &MeshOptions,
    out: &std::path::Path,
    repetitions: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut rows = Vec::new();
    for iteration in 0..repetitions {
        let start = Instant::now();
        let body = catallax_app::build(recipe, options, &mut Decider::default())?;
        let app_ms = start.elapsed().as_secs_f64() * 1000.0;
        rows.push(serde_json::json!({
            "iteration":iteration,"app_ms":app_ms,
            "faces":body.names.faces.len(),
            "triangles":body.display.value.triangles.len(),
            "chord_tolerance_mm":options.chord_tolerance.to_string(),
            "max_work":catallax_app::DEFAULT_BUILD_WORK,
        }));
        if iteration + 1 == repetitions {
            body.display
                .value
                .write_obj(&mut BufWriter::new(fs::File::create(
                    out.join("display.obj"),
                )?))?;
            fs::write(out.join("mesh.txt"), format!("{:?}\n", body.display))?;
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

fn edit_rebuilds(
    previous: &Recipe,
    recipe: &Recipe,
    options: &MeshOptions,
    out: &std::path::Path,
    repetitions: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut rows = Vec::new();
    for iteration in 0..repetitions {
        let before =
            catallax_app::build_incremental(previous, None, options, &mut Decider::default())?;
        let mut before_obj = Vec::new();
        before.display.value.write_obj(&mut before_obj)?;
        let start = Instant::now();
        let body = catallax_app::build_incremental(
            recipe,
            Some(before.checkpoint.clone()),
            options,
            &mut Decider::default(),
        )?;
        let edit_ms = start.elapsed().as_secs_f64() * 1000.0;
        let mut edited_obj = Vec::new();
        body.display.value.write_obj(&mut edited_obj)?;
        let fresh = catallax_app::build(recipe, options, &mut Decider::default())?;
        let mut fresh_obj = Vec::new();
        fresh.display.value.write_obj(&mut fresh_obj)?;
        assert_eq!(
            edited_obj, fresh_obj,
            "edited display differs from fresh build"
        );
        assert_eq!(body.names.faces.len(), fresh.names.faces.len());
        let retained = catallax_app::build_incremental(
            previous,
            Some(before.checkpoint),
            options,
            &mut Decider::default(),
        )?;
        let mut retained_obj = Vec::new();
        retained.display.value.write_obj(&mut retained_obj)?;
        assert_eq!(before_obj, retained_obj, "edit changed retained display");
        rows.push(serde_json::json!({
            "iteration":iteration,"edit_ms":edit_ms,
            "faces":body.names.faces.len(),
            "triangles":body.display.value.triangles.len(),
            "chord_tolerance_mm":options.chord_tolerance.to_string(),
            "max_work":catallax_app::DEFAULT_BUILD_WORK,
        }));
        if iteration + 1 == repetitions {
            fs::write(out.join("display.obj"), edited_obj)?;
        }
    }
    println!("{}", serde_json::to_string(&rows)?);
    fs::write(out.join("samples.json"), serde_json::to_vec_pretty(&rows)?)?;
    Ok(())
}
