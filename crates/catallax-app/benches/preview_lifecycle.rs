//! Development-only complete preview wait and deferred certified-display timing.
use catallax_app::lifecycle::{Controller, Status, WorkAllocation};
use catallax_mesh::{MeshOptions, PreviewOptions};
use catallax_numbers::Decider;
use catallax_recipe::Recipe;
use serde_json::json;
use std::{env, fs, time::Instant};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args().skip(1).filter(|arg| arg != "--bench").collect();
    if args.is_empty() || args.len() > 2 {
        return Err("usage: preview_lifecycle OUTPUT.json [REPETITIONS]".into());
    }
    let repetitions: usize = args.get(1).map_or(Ok(3), |v| v.parse())?;
    if repetitions == 0 {
        return Err("repetitions must be positive".into());
    }
    let fixtures = [
        (
            "fillet",
            include_str!("fine_mesh/fillet.json"),
            "0.0005462600113499065",
        ),
        (
            "chamfer-P2",
            include_str!("../../../parts/P2-collar.json"),
            "0.00036",
        ),
        ("chamfer", include_str!("fine_mesh/chamfer.json"), "0.00036"),
        (
            "countersink",
            include_str!("fine_mesh/countersink.json"),
            "0.000722841614740048",
        ),
        (
            "P1",
            include_str!("../../../parts/P1-bracket.json"),
            "0.0016638509548634457",
        ),
        (
            "P3",
            include_str!("../../../parts/P3-filleted-boss.json"),
            "0.0005462600113499065",
        ),
    ];
    let mut rows = Vec::new();
    for (name, input, chord) in fixtures {
        let recipe = Recipe::parse(input)?;
        let options = MeshOptions {
            chord_tolerance: chord.parse()?,
            angular_deviation: Some("0.03490658503988659".parse()?),
            ..MeshOptions::default()
        };
        for iteration in 0..=repetitions {
            let controller = Controller::default();
            // Explicit equal split of one 20M request, not 20M per phase.
            let allocation = WorkAllocation::application(10_000_000).unwrap();
            let start = Instant::now();
            let mut preview_decider = Decider::default();
            let result = controller.begin().preview(
                &recipe,
                &PreviewOptions::from(&options),
                &options,
                allocation,
                &mut preview_decider,
            );
            let preview_ms = start.elapsed().as_secs_f64() * 1000.0;
            let mut row = json!({"fixture":name,"iteration":iteration,"warmup":iteration == 0,
                "total_allowance":catallax_app::DEFAULT_BUILD_WORK,"preview_cap":allocation.preview,
                "proof_cap":allocation.proof,"preview_wait_ms":preview_ms});
            match result {
                Ok(Some(job)) => {
                    let preview = controller.latest().unwrap();
                    assert!(matches!(preview.status(), Status::Unproved));
                    row["preview_work"] = json!(preview.work().preview.used);
                    let display = preview.preview().brep().unwrap();
                    row["preview_vertices"] = json!(display.mesh.positions.len());
                    row["preview_triangles"] = json!(display.mesh.triangles.len());
                    // No artificial delay. Real callers can postpone this job.
                    let proof_start = Instant::now();
                    let completion = job.prove(&mut Decider::default());
                    assert!(completion.publish());
                    row["proof_phase_ms"] = json!(proof_start.elapsed().as_secs_f64() * 1000.0);
                    row["proof_completion_ms"] = json!(start.elapsed().as_secs_f64() * 1000.0);
                    let completed = controller.latest().unwrap();
                    let work = completed.work().proof.unwrap();
                    row["proof_work"] = json!(work.used);
                    row["proof_exhausted"] = json!(work.exhausted);
                    match completed.status() {
                        Status::Confirmed(body) => {
                            row["status"] = json!("confirmed");
                            row["certified_triangles"] = json!(body.display.value.triangles.len());
                        }
                        Status::Failed(error) => {
                            row["status"] = json!("failed");
                            row["proof_refusal"] = json!(error.code());
                        }
                        Status::Unproved => unreachable!(),
                    }
                }
                Err(error) => {
                    row["status"] = json!("preview-refused");
                    row["preview_refusal"] = json!(error.code());
                    row["detail"] = json!(error.to_string());
                }
                Ok(None) => unreachable!(),
            }
            println!("{row}");
            rows.push(row);
        }
    }
    fs::write(&args[0], serde_json::to_vec_pretty(&rows)?)?;
    Ok(())
}
