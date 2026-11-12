//! Observe the public recorded-repair entry under the existing import allowance.

use catallax_import::{ImportOptions, ImportRefusal, apply, sha256_hex};
use catallax_numbers::with_operation_budget;
use catallax_topology::Model;
use serde_json::{Value, json};

const IMPORT_WORK: u64 = 20_000_000;

fn evaluate(bytes: &[u8]) -> Value {
    let mut model = Model::new();
    let (result, work) = with_operation_budget(IMPORT_WORK, || {
        apply(
            &mut model,
            bytes,
            "source.step",
            &ImportOptions::default(),
            &[],
            None,
        )
    });
    let result = if work.exhausted {
        Err(ImportRefusal::ResourceLimit)
    } else {
        result
    };
    let mut row = match result {
        Ok(imported) => {
            let valid = imported
                .bodies
                .iter()
                .all(|(_, body)| model.check(*body).is_empty());
            json!({"outcome": "admitted", "bodies": imported.bodies.len(),
                "left_open": imported.left_open.len(), "topology_valid": valid,
                "record_format": imported.record.format,
                "changes": imported.record.changes.len(),
                "source_incidence": imported.record.source_incidence,
                "face_tolerances": imported.record.face_tolerances,
                "derivative_sha256": imported.derivative.as_deref().map(sha256_hex)})
        }
        Err(error) => json!({"outcome": "refused", "refusal": error}),
    };
    row["source_sha256"] = json!(sha256_hex(bytes));
    row["work_used"] = json!(work.used);
    row["work_limit"] = json!(IMPORT_WORK);
    row
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args_os().skip(1);
    let input = args.next().ok_or("expected one STEP input")?;
    if args.next().is_some() {
        return Err("expected one STEP input".into());
    }
    println!("{}", evaluate(&std::fs::read(input)?));
    Ok(())
}

#[cfg(test)]
#[path = "../tests/common/mod.rs"]
mod common;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clean_and_repaired_sources_use_the_recorded_entry() {
        let clean = common::box_step(common::BoxFixture::new(["1", "1", "1"]));
        let row = evaluate(clean.as_bytes());
        assert_eq!(row["outcome"], "admitted");
        assert_eq!(row["record_format"], "catallax-import-repair/2");
        assert_eq!(row["source_incidence"], json!([]));
        let source =
            common::box_with_placement(5, Some(["0.000001", "0", "1"]), Some(["1", "0", "0"]));
        let row = evaluate(source.as_bytes());
        assert_eq!(row["outcome"], "admitted");
        assert_eq!(row["record_format"], "catallax-import-repair/3");
        assert_eq!(row["source_incidence"].as_array().unwrap().len(), 1);
        assert_eq!(row["topology_valid"], true);
        assert_eq!(row["work_limit"], IMPORT_WORK);
        assert!(row["derivative_sha256"].is_string());
    }
}
