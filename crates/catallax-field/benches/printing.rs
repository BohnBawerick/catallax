//! Explicit print timing and byte capture, separate from ordinary construction.
use std::{env, fs, path::PathBuf, time::Instant};

use catallax_field::{Plan, printing};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::{Number, Recipe};
use serde_json::json;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args_os().skip(1).collect();
    if args.len() != 3 {
        return Err("usage: printing REQUEST.json OUTPUT_DIRECTORY REPETITIONS".into());
    }
    let request: serde_json::Value = serde_json::from_str(&fs::read_to_string(&args[0])?)?;
    if request["format"] != "catallax-print/1" {
        return Err("expected catallax-print/1".into());
    }
    let recipe = Recipe::parse(&request["recipe"].to_string())?;
    let tolerance: Number = serde_json::from_value(request["tolerance"].clone())?;
    let defaults = printing::Options::default();
    let options = printing::Options {
        tolerance: tolerance.value().clone(),
        max_cells: request["max_cells"]
            .as_u64()
            .map(|v| v as usize)
            .unwrap_or(defaults.max_cells),
        max_vertices: request["max_vertices"]
            .as_u64()
            .map(|v| v as usize)
            .unwrap_or(defaults.max_vertices),
    };
    let out = PathBuf::from(&args[1]);
    fs::create_dir_all(&out)?;
    let repetitions: usize = args[2].to_str().ok_or("invalid repetitions")?.parse()?;
    if repetitions == 0 {
        return Err("repetitions must be positive".into());
    }
    let mut rows = Vec::new();
    for iteration in 0..repetitions {
        let start = Instant::now();
        let (result, work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::default();
            let plan = Plan::compile_checked(&recipe, &mut dec)?;
            let build_ms = start.elapsed().as_secs_f64() * 1000.0;
            let mesh_start = Instant::now();
            let mesh = printing::mesh_plan(&plan, &options, &mut dec)?;
            Ok::<_, catallax_numbers::Refusal>((
                mesh,
                build_ms,
                mesh_start.elapsed().as_secs_f64() * 1000.0,
            ))
        });
        let (mesh, build_ms, mesh_ms) = result?;
        let stl_start = Instant::now();
        mesh.write_stl(&mut fs::File::create(out.join("print.stl"))?)?;
        let stl_ms = stl_start.elapsed().as_secs_f64() * 1000.0;
        let obj_start = Instant::now();
        mesh.write_obj(&mut fs::File::create(out.join("display.obj"))?)?;
        let obj_ms = obj_start.elapsed().as_secs_f64() * 1000.0;
        let total_ms = start.elapsed().as_secs_f64() * 1000.0;
        // Serializer-only timings into memory, outside the complete-request total.
        let (mut stl, mut obj) = (Vec::new(), Vec::new());
        let memory_start = Instant::now();
        mesh.write_stl(&mut stl)?;
        let stl_memory_ms = memory_start.elapsed().as_secs_f64() * 1000.0;
        let memory_start = Instant::now();
        mesh.write_obj(&mut obj)?;
        let obj_memory_ms = memory_start.elapsed().as_secs_f64() * 1000.0;
        if stl != fs::read(out.join("print.stl"))? || obj != fs::read(out.join("display.obj"))? {
            return Err("memory and file output differ".into());
        }
        rows.push(json!({
            "iteration": iteration, "build_ms": build_ms, "mesh_ms": mesh_ms,
            "stl_ms": stl_ms, "obj_ms": obj_ms,
            "stl_memory_ms": stl_memory_ms, "obj_memory_ms": obj_memory_ms,
            "total_ms": total_ms,
            "target_ms": 1000, "work": work.used, "exhausted": work.exhausted,
            "max_work": 100_000_000, "max_cells": options.max_cells,
            "max_vertices": options.max_vertices, "cells": mesh.cells,
            "vertices": mesh.positions.len(), "triangles": mesh.triangles.len(),
            "tolerance_mm": options.tolerance.to_string(),
            "deviation_bound_mm": mesh.max_deviation.to_string(),
        }));
    }
    fs::write(out.join("samples.json"), serde_json::to_vec_pretty(&rows)?)?;
    println!("{}", serde_json::to_string(&rows)?);
    Ok(())
}
