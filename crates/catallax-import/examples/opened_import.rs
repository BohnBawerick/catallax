//! Frozen-census observer for the tier consumer, not another admission policy.
//! Usage: opened_import SOURCE.step SHA256 NEW_CAPTURE_DIR [--retain-controls]
//! One 20M open query and one aggregate 20M export query. Observer metadata,
//! serialization and disk I/O are outside numerical work; process time is bounded.

use std::{
    io::{self, Write},
    path::Path,
};

use catallax_import::{
    numerical::StepImport,
    sha256_hex,
    tolerant::{
        self,
        opened::{self, BodyView, OpenedStep},
    },
};
use catallax_numbers::{Decider, Refusal, with_operation_budget};
use catallax_recipe::StepId;
use catallax_topology::{BodyKind, numerical::DEFAULT_WORK};
use serde_json::{Value, json};

fn error(error: &tolerant::Error) -> Value {
    match error {
        tolerant::Error::Source(error) => serde_json::to_value(error).unwrap(),
        tolerant::Error::Numerical(error) => json!({"code":error.code,"detail":error.detail}),
    }
}

fn emit(value: &Value) {
    println!("{value}");
    io::stdout().flush().unwrap();
}

fn bodies(opened: &OpenedStep) -> Vec<Value> {
    opened.bodies().iter().enumerate().map(|(index, outcome)| {
        let source = &opened.source().bodies()[index];
        let (proof, counts, exportable) = match opened.body(index).unwrap() {
            BodyView::Certified { model, body } => {
                let c = model.counts(body);
                ("certified", json!({"faces":c.faces,"edges":c.edges,"vertices":c.vertices}), true)
            }
            BodyView::Numerical(body) => {
                let c = body.counts();
                ("numerical", json!({"faces":c.faces,"edges":c.edges,"vertices":c.vertices}), true)
            }
            BodyView::Surfaces { .. } => ("retained_source_surfaces", Value::Null, false),
        };
        json!({"body":index,"source_body":source.name(),"declared_solid":source.declared_solid(),
            "tier":outcome.tier(),"approximate":outcome.approximate(),"proof":proof,
            "checked_brep_counts":counts,"exportable":exportable,
            "source_entity_inventory":{"faces":source.shells().iter().map(|s|s.faces.len()).sum::<usize>(),"edges":source.edges().len(),"vertices":source.vertices().len()},
            "refusal":outcome.refusal().map(error),
            "healing_refusal":outcome.healing_failure().map(|f|json!({"refusal":error(&f.error),"stage":f.failed_stage,"body":f.body,
                "stages":f.stages,"preparation_work_used":f.preparation_work_used,"publication_work_used":f.publication_work_used})),
            "healing_issues":outcome.healing().map(|h|&h.report().issues)})
    }).collect()
}

fn diagnostics(opened: &OpenedStep) -> Vec<Value> {
    opened.bodies().iter().enumerate().map(|(index, outcome)| {
        let admission = match opened.admission(index) {
            Some(StepImport::Numerical(a)) => json!({"source_record":a.record(),"source_interpretation":a.source_admission().map(|s|json!({
                "method":s.method,"strict_refusal":s.strict_refusal,"tolerance_mm":s.tolerance_mm.to_string(),
                "linear_spline_images":s.linear_spline_images,
                "circular_spline_images":s.circular_spline_images,
                "topology_repairs":s.topology_repairs,"intersection_checks":s.intersection_checks}))}),
            _ => Value::Null,
        };
        json!({"body":index,"admission":admission,
            "healing":outcome.healing().map(|h|json!({"report":h.report(),"tolerances":h.tolerances()}))})
    }).collect()
}

fn export(
    opened: &OpenedStep,
    index: usize,
    dec: &mut Decider,
) -> Result<Option<(Vec<u8>, &'static str)>, Refusal> {
    match opened.body(index).unwrap() {
        BodyView::Certified { model, body } => {
            let mut bytes = Vec::new();
            let mut options = catallax_step::StepOptions::default();
            if opened.bodies()[index].approximate() {
                options.product_name = "catallax approximate source repair".into();
            }
            catallax_step::write_step(model, &[body], &options, dec, &mut bytes)?;
            Ok(Some((
                bytes,
                if model.body(body).kind == BodyKind::Solid {
                    "solid"
                } else {
                    "surface"
                },
            )))
        }
        BodyView::Numerical(body) => {
            let file = catallax_step::numerical::build(body, &Default::default(), dec)?;
            file.check_written(file.bytes(), dec)?;
            Ok(Some((file.bytes().to_vec(), "solid")))
        }
        // No exact/numerical B-rep was admitted, so there is no checked STEP
        // authority. Retained input bytes are separate from a body export.
        BodyView::Surfaces { .. } => Ok(None),
    }
}

fn run(
    bytes: &[u8],
    expected: &str,
    out: &Path,
    retain: bool,
) -> Result<Value, Box<dyn std::error::Error>> {
    let actual = sha256_hex(bytes);
    let (result, work) = with_operation_budget(DEFAULT_WORK, || {
        opened::open_step(
            bytes,
            expected,
            &StepId::new("opened-import-census").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
    });
    let opened = match result {
        Ok(opened) => opened, // An optional admission refusal may retain surfaces.
        Err(refusal) => {
            return Ok(
                json!({"stage":"open","outcome":"refused","source_sha256":actual,
            "work_used":work.used,"work_exhausted":work.exhausted,"work_limit":DEFAULT_WORK,"refusal":error(&refusal)}),
            );
        }
    };
    let body_info = bodies(&opened);
    let checkpoint = json!({"stage":"open","outcome":"opened","source_sha256":actual,
        "bodies":body_info,"work_used":work.used,"work_exhausted":work.exhausted,"work_limit":DEFAULT_WORK,
        "read_work_used":opened.read_work_used(),"file_admission_work_used":opened.file_admission_work_used(),
        "first_refusal":opened.first_refusal().map(error)});
    emit(&checkpoint);
    std::fs::write(out.join("opened.json"), serde_json::to_vec(&checkpoint)?)?;
    if retain {
        std::fs::write(out.join("source.step"), bytes)?;
        std::fs::write(
            out.join("controls.json"),
            serde_json::to_vec(&opened.saved_controls())?,
        )?;
    }
    let mut files = Vec::new();
    let mut failed_body = None;
    let (exported, export_work) = with_operation_budget(DEFAULT_WORK, || {
        let mut dec = Decider::default();
        for index in 0..opened.bodies().len() {
            failed_body = Some(index);
            dec.charge_work(1)?;
            if let Some((bytes, kind)) = export(&opened, index, &mut dec)? {
                let file = format!("body-{index}.step");
                let path = out.join(&file);
                std::fs::write(&path, &bytes)
                    .map_err(|e| Refusal::new("observer.write", e.to_string()))?;
                let saved = std::fs::read(path)
                    .map_err(|e| Refusal::new("observer.reread", e.to_string()))?;
                if saved != bytes {
                    return Err(Refusal::new(
                        "observer.byte-identity",
                        "written bytes changed",
                    ));
                }
                files.push(json!({"body":index,"kind":kind,"file":file,"bytes":saved.len(),"sha256":sha256_hex(&saved)}));
            }
        }
        failed_body = None;
        Ok::<_, Refusal>(())
    });
    let refusal = exported
        .err()
        .map(|r| json!({"code":r.code,"detail":r.detail}));
    Ok(
        json!({"stage":"export","outcome":if refusal.is_none()&&!export_work.exhausted {"complete"} else {"refused"},
        "source_sha256":actual,"open":checkpoint,"files":files,"failed_body":failed_body,"refusal":refusal,
        "work_used":export_work.used,"work_exhausted":export_work.exhausted,"work_limit":DEFAULT_WORK,
        "diagnostics":diagnostics(&opened),"complete_source_body_bound_established":false}),
    )
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if !(args.len() == 3 || (args.len() == 4 && args[3] == "--retain-controls")) {
        return Err("expected SOURCE.step SHA256 NEW_CAPTURE_DIR [--retain-controls]".into());
    }
    let expected = args[1].to_str().ok_or("SHA-256 must be text")?;
    if expected.len() != 64 || !expected.bytes().all(|c| c.is_ascii_hexdigit()) {
        return Err("expected SHA-256".into());
    }
    let bytes = std::fs::read(&args[0])?;
    let out = Path::new(&args[2]);
    std::fs::create_dir(out)?;
    let response = run(&bytes, expected, out, args.len() == 4)?;
    std::fs::write(out.join("response.json"), serde_json::to_vec(&response)?)?;
    emit(
        &json!({"stage":response["stage"],"outcome":response["outcome"],"source_sha256":response["source_sha256"]}),
    );
    Ok(())
}
