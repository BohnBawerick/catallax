//! Read, reconstruct and reopen tolerance-bearing source topology.
//! Usage: cargo run --release -p catallax-import --example tolerant_topology --
//! SOURCE.step EXPECTED_SHA256 RECORD.json
//! This checks storage fidelity, not solid validity or repaired geometry.
use catallax_import::{read_step, tolerant};
use catallax_numbers::{Decider, with_operation_budget};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 3 {
        return Err("expected SOURCE.step EXPECTED_SHA256 RECORD.json".into());
    }
    let bytes = std::fs::read(&args[0])?;
    let (result, work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let file = tolerant::read_step(
            &bytes,
            &args[1],
            tolerant::ToleranceOptions::default(),
            &mut dec,
        )?;
        let original = read_step(&bytes)?;
        assert_eq!(
            file.definition(&mut dec)?,
            original,
            "source geometry/topology changed"
        );
        let record = file.record();
        assert_eq!(tolerant::reopen(&bytes, &record, &mut dec)?, file);
        Ok::<_, tolerant::Error>((file, record))
    });
    let (file, record) = result.map_err(|e| format!("{e:?}"))?;
    if work.exhausted {
        return Err("resource_limit".into());
    }
    std::fs::write(&args[2], serde_json::to_vec_pretty(&record)?)?;
    println!(
        "{}",
        serde_json::json!({
            "source_sha256": file.source_sha256(), "definition_unchanged": true,
            "reopened": true, "work_used": work.used,
            "solid_validity_checked": false,
            "bodies": file.bodies().iter().map(|b| serde_json::json!({
                "source": b.name(), "declared_solid": b.declared_solid(),
                "shells": b.shells().len(), "faces": b.shells().iter().map(|s| s.faces.len()).sum::<usize>(),
                "edges": b.edges().len(), "vertices": b.vertices().len(),
                "maximum_tolerance_mm": b.maximum_tolerance_mm().to_string(),
                "tolerance_changes": b.tolerance_changes().len(),
            })).collect::<Vec<_>>()
        })
    );
    Ok(())
}
