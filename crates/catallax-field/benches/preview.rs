//! Shared-machine diagnostics and exact display payloads, without speed gates.

use catallax_field::{
    Plan,
    display::{Options, PreviewOptions, mesh_plan, preview_plan},
};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::json;
use std::{io::Write, path::Path, time::Instant};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 5 {
        return Err(
            "usage: preview RECIPE OUTPUT_DIR certified|preview TOLERANCE_DENOMINATOR".into(),
        );
    }
    let recipe = Recipe::parse(&std::fs::read_to_string(&args[1])?)?;
    let out = Path::new(&args[2]);
    std::fs::create_dir_all(out)?;
    for name in ["mesh.bin", "display.obj", "report.json"] {
        match std::fs::remove_file(out.join(name)) {
            Ok(()) => (),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => (),
            Err(error) => return Err(error.into()),
        }
    }
    let options = Options {
        tolerance: Rational::from_ratio(1, args[4].parse()?).unwrap(),
        ..Options::default()
    };
    let mut dec = Decider::default();
    if let Ok(seconds) = std::env::var("CATALLAX_PREVIEW_DEADLINE_SECONDS") {
        dec.set_deadline(Some(
            Instant::now() + std::time::Duration::from_secs(seconds.parse()?),
        ));
    }
    let mut obj = None;
    let started = Instant::now();
    let (result, work) = with_operation_budget(100_000_000, || {
        let plan = Plan::compile_checked(&recipe, &mut dec)?;
        if args[3] == "certified" {
            mesh_plan(&plan, &options, &mut dec).map(|m| {
                let mut bytes = Vec::new();
                m.write_obj(&mut bytes).expect("write to memory");
                obj = Some(bytes);
                (m.positions, m.triangles, m.cells, 0)
            })
        } else if args[3] == "preview" {
            preview_plan(&plan, &PreviewOptions::from(&options), &mut dec)
                .map(|m| (m.positions, m.triangles, m.cells, m.unresolved_cells))
        } else {
            panic!("expected certified or preview");
        }
    });
    let mut row = json!({"environment":"shared machine; informational only; unloaded machine not guaranteed",
        "mode":args[3], "tolerance":options.tolerance.to_string(),
        "proof":if args[3] == "preview" { "unproved" } else { "certified" },
        "wall_ms":started.elapsed().as_secs_f64()*1000.0, "work":work.used,
        "exhausted":work.exhausted});
    if let Some(bytes) = obj {
        std::fs::write(out.join("display.obj"), bytes)?;
    }
    match result {
        Ok((positions, triangles, cells, gaps)) => {
            row["vertices"] = json!(positions.len());
            row["triangles"] = json!(triangles.len());
            row["cells"] = json!(cells);
            row["unresolved_cells"] = json!(gaps);
            let mut payload = std::fs::File::create(out.join("mesh.bin"))?;
            payload.write_all(&(positions.len() as u64).to_le_bytes())?;
            for point in positions {
                for coordinate in point {
                    payload.write_all(&coordinate.to_bits().to_le_bytes())?;
                }
            }
            for triangle in triangles {
                for index in triangle.vertices.into_iter().chain([triangle.face.index()]) {
                    payload.write_all(&(index as u64).to_le_bytes())?;
                }
            }
        }
        Err(refusal) => {
            row["refusal"] = json!({"code":refusal.code,"detail":refusal.detail});
        }
    }
    std::fs::write(out.join("report.json"), serde_json::to_vec_pretty(&row)?)?;
    println!("{row}");
    Ok(())
}
