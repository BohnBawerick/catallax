//! Compare two actual STEP imports, retaining separate charged stage costs.

use std::collections::BTreeMap;

use catallax_import::{ImportRefusal, apply, correspondence, sha256_hex};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_topology::Model;
use serde_json::json;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let paths: Vec<_> = std::env::args_os().skip(1).collect();
    if paths.len() != 2 {
        return Err("expected old STEP path and replacement STEP path".into());
    }
    let sources = [std::fs::read(&paths[0])?, std::fs::read(&paths[1])?];
    let mut imports = Vec::new();
    let mut costs = BTreeMap::new();
    for (stage, source) in ["old_import", "new_import"].into_iter().zip(&sources) {
        let (result, work) = with_operation_budget(20_000_000, || {
            apply(
                &mut Model::new(),
                source,
                "source.step",
                &Default::default(),
                &[],
                None,
            )
        });
        costs.insert(stage, work.used);
        let result = if work.exhausted {
            Err(ImportRefusal::ResourceLimit)
        } else {
            result
        };
        match result {
            Ok(imported) => imports.push(imported),
            Err(refusal) => {
                println!(
                    "{}",
                    json!({"status":"refused", "stage":stage, "refusal":refusal,
                        "source_sha256":sha256_hex(source), "work":costs,
                        "stage_work_limit":20_000_000})
                );
                return Ok(());
            }
        }
    }
    let (report, work) = with_operation_budget(20_000_000, || {
        correspondence::compare(&imports[0], &imports[1], &mut Decider::default())
    });
    costs.insert("comparison", work.used);
    println!(
        "{}",
        match report {
            Ok(report) => json!({"status":"compared", "correspondence":report,
                "work":costs, "stage_work_limit":20_000_000}),
            Err(refusal) => json!({"status":"refused", "stage":"comparison",
                "refusal":{"code":refusal.code,"detail":refusal.detail},
                "work":costs, "stage_work_limit":20_000_000}),
        }
    );
    Ok(())
}
