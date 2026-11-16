//! Observe import and export separately for a frozen set of census inputs.
//!
//! Admission policy is unchanged. Each import and the aggregate export of all
//! its bodies have the existing 20M kernel allowance, with separate work rows.
//! Observer I/O, metadata hashing and external OCCT work are not kernel work.

use std::io::{self, Write};
use std::path::Path;

use catallax_import::{numerical::StepImport, sha256_hex};
use catallax_numbers::{Decider, Refusal, with_operation_budget};
use catallax_recipe::StepId;
use catallax_step::numerical::{NumericalStepOptions, build as numerical_step};
use catallax_topology::numerical::{DEFAULT_WORK, NumericalOptions};
use serde_json::{Value, json};

fn import_refusal(error: catallax_import::numerical::ImportError) -> Value {
    match error {
        catallax_import::numerical::ImportError::Source(error) => {
            serde_json::to_value(error).unwrap()
        }
        catallax_import::numerical::ImportError::Numerical(error) => {
            json!({"code":error.code,"detail":error.detail})
        }
    }
}

fn export_one(imported: &StepImport, index: usize, dec: &mut Decider) -> Result<Vec<u8>, Refusal> {
    match imported {
        StepImport::Certified(imported) => {
            let mut bytes = Vec::new();
            catallax_step::write_step(
                imported.model(),
                &[imported.bodies()[index]],
                &catallax_step::StepOptions::default(),
                dec,
                &mut bytes,
            )?;
            Ok(bytes)
        }
        StepImport::Numerical(imported) => {
            let file = numerical_step(
                &imported.bodies()[index].1,
                &NumericalStepOptions::default(),
                dec,
            )?;
            file.check_written(file.bytes(), dec)?;
            Ok(file.bytes().to_vec())
        }
    }
}

fn run(bytes: &[u8], expected: &str, out: &Path, mut checkpoint: impl FnMut(&Value)) -> Value {
    let actual = sha256_hex(bytes);
    let (result, import_work) = with_operation_budget(DEFAULT_WORK, || {
        catallax_import::numerical::import_step(
            bytes,
            expected,
            &StepId::new("export-validity-census").unwrap(),
            &NumericalOptions::default(),
            &mut Decider::default(),
        )
    });
    let imported = match result {
        Ok(imported) if !import_work.exhausted => imported,
        result => {
            return json!({"stage":"import","outcome":"refused","source_sha256":actual,
                "work_used":import_work.used,"work_limit":DEFAULT_WORK,
                "refusal":if import_work.exhausted { json!({"code":"resource_limit"}) }
                    else { import_refusal(result.unwrap_err()) }});
        }
    };
    let (category, bodies) = match &imported {
        StepImport::Certified(imported) => ("certified", imported.bodies().len()),
        StepImport::Numerical(imported) => ("numerical approximate", imported.bodies().len()),
    };
    // Observation is outside kernel work. Bind numerical source face keys and
    // kinds to their generating definition order; external checks must map the
    // actual emitted entities rather than assume this is OCCT traversal order.
    let numerical_definitions = match &imported {
        StepImport::Numerical(imported) => Some(imported.bodies().iter().enumerate().map(|(index, (key, body))| {
            json!({"body":index,"source_body":key,"faces":body.definition().faces.iter().enumerate().map(|(index,face)| {
                json!({"definition_index":index,"source_key":face.key,"surface_kind":format!("{:?}",face.surface.kind())})
            }).collect::<Vec<_>>()})
        }).collect::<Vec<_>>()),
        StepImport::Certified(_) => None,
    };
    let admission = json!({"stage":"import","outcome":"admitted","source_sha256":actual,
        "import_category":category,"bodies":bodies,"work_used":import_work.used,"work_limit":DEFAULT_WORK,
        "numerical_definitions":numerical_definitions});
    // An export timeout must not erase the independently completed import row.
    checkpoint(&admission);
    let mut files = Vec::new();
    let mut failed_body = None;
    let (result, export_work) = with_operation_budget(DEFAULT_WORK, || {
        let mut dec = Decider::default();
        for index in 0..bodies {
            failed_body = Some(index);
            let bytes = export_one(&imported, index, &mut dec)?;
            let name = format!("body-{index}.step");
            let path = out.join(&name);
            std::fs::write(&path, &bytes)
                .map_err(|e| Refusal::new("observer.write", e.to_string()))?;
            let reread =
                std::fs::read(&path).map_err(|e| Refusal::new("observer.reread", e.to_string()))?;
            if reread != bytes {
                return Err(Refusal::new(
                    "observer.byte-identity",
                    "written STEP bytes changed",
                ));
            }
            files.push(
                json!({"body":index,"file":name,"bytes":bytes.len(),"sha256":sha256_hex(&reread)}),
            );
        }
        failed_body = None;
        Ok(())
    });
    let refusal = if export_work.exhausted {
        Some(json!({"code":"resource_limit"}))
    } else {
        result
            .err()
            .map(|e| json!({"code":e.code,"detail":e.detail}))
    };
    json!({"stage":"export","outcome":if refusal.is_none() {"complete"} else {"refused"},
        "source_sha256":actual,"import_category":category,"bodies":bodies,"files":files,
        "import":admission,"work_used":export_work.used,"work_limit":DEFAULT_WORK,
        "failed_body":failed_body,"refusal":refusal})
}

fn emit(value: &Value) {
    println!("{value}");
    io::stdout().flush().unwrap();
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args_os().skip(1).collect::<Vec<_>>();
    if args.len() != 3 {
        return Err("expected source STEP, pinned SHA-256 and new capture directory".into());
    }
    let expected = args[1].to_str().ok_or("SHA-256 must be UTF-8")?;
    if expected.len() != 64 || !expected.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err("expected a 64-digit source SHA-256".into());
    }
    let bytes = std::fs::read(&args[0])?;
    let out = Path::new(&args[2]);
    std::fs::create_dir(out)?;
    let response = run(&bytes, expected, out, emit);
    std::fs::write(out.join("response.json"), format!("{response}\n"))?;
    emit(&response);
    Ok(())
}

#[cfg(test)]
#[path = "../tests/common/mod.rs"]
mod common;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn both_source_bodies_are_exported_after_an_import_checkpoint() {
        let a = common::BoxFixture::new(["1", "2", "3"]);
        let mut b = common::BoxFixture::new(["2", "3", "4"]);
        b.origin[0] = common::r("10");
        let bytes = common::boxes_step(&[a, b]).into_bytes();
        let out =
            std::env::temp_dir().join(format!("catallax-export-census-{}", std::process::id()));
        std::fs::create_dir(&out).unwrap();
        let mut checkpoints = Vec::new();
        let response = run(&bytes, &sha256_hex(&bytes), &out, |row| {
            checkpoints.push(row.clone())
        });
        assert_eq!(checkpoints.len(), 1);
        assert_eq!(checkpoints[0]["outcome"], "admitted");
        assert_eq!(checkpoints[0]["bodies"], 2);
        assert_eq!(response["outcome"], "complete");
        let files = response["files"].as_array().unwrap();
        assert_eq!(files.len(), 2);
        for entry in files {
            let bytes = std::fs::read(out.join(entry["file"].as_str().unwrap())).unwrap();
            assert_eq!(sha256_hex(&bytes), entry["sha256"]);
            assert_eq!(
                std::str::from_utf8(&bytes)
                    .unwrap()
                    .matches("MANIFOLD_SOLID_BREP(")
                    .count(),
                1
            );
        }
        if let Ok(path) = std::env::var("CATALLAX_EXPORT_VALIDITY_EVIDENCE") {
            std::fs::write(Path::new(&path).join("two-boxes.step"), &bytes).unwrap();
        }
        std::fs::remove_dir_all(out).unwrap();
    }

    #[test]
    fn a_changed_source_hash_produces_no_admission_checkpoint_or_export() {
        let bytes = common::box_step(common::BoxFixture::new(["1", "2", "3"]));
        let response = run(
            bytes.as_bytes(),
            &"0".repeat(64),
            Path::new("unused"),
            |_| panic!("wrong source admitted"),
        );
        assert_eq!(response["stage"], "import");
        assert_eq!(response["outcome"], "refused");
        assert!(response.get("files").is_none());
    }
}
