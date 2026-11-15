//! Replay an opt-in numerical recipe and publish approximate STEP after its checks.
use catallax_numbers::{Decider, with_operation_budget};
use catallax_step::numerical::{self, NumericalStepOptions};
use serde_json::json;
use std::{error::Error, path::Path};

fn run(input: &Path, output: &Path) -> Result<(), Box<dyn Error>> {
    let saved = std::fs::read_to_string(input)?;
    let (result, work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let built = catallax_features::replay_numerical(&saved, &mut dec)?;
        let artifact = numerical::build(built.body(), &NumericalStepOptions::default(), &mut dec)?;
        Ok::<_, catallax_numbers::Refusal>((built, artifact))
    });
    let (built, artifact) = result?;
    if work.exhausted {
        return Err("shared replay/export allowance exhausted".into());
    }
    std::fs::write(output, artifact.bytes())?;
    artifact.check_written(&std::fs::read(output)?, &mut Decider::default())?;
    println!(
        "{}",
        json!({"status":"numerical","exact":false,"step":"approximate",
        "tolerance_mm":built.body().evidence().tolerance_mm().to_string(),
        "checks":built.body().evidence().checks().len(),
        "repairs":built.body().evidence().repairs().len(),
        "history":built.body().history().iter().map(|stage| json!({"tolerance_mm":stage.options.tolerance_mm.to_string(),"checks":stage.evidence.checks().len(),"repairs":stage.evidence.repairs().len()})).collect::<Vec<_>>(),
        "origins":built.face_origins(),"work":work.used})
    );
    Ok(())
}
fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args_os().collect();
    if args.len() != 3 {
        return Err("usage: numerical_replay INPUT.json OUTPUT.step".into());
    }
    let output = Path::new(&args[2]);
    if let Err(error) = run(Path::new(&args[1]), output) {
        match std::fs::remove_file(output) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(e.into()),
        }
        eprintln!(
            "{}",
            json!({"status":"refused","exact":false,"step":null,"detail":error.to_string()})
        );
        return Err(error);
    }
    Ok(())
}
