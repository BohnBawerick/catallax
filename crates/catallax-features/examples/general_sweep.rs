//! Build, parameter-edit, save and reopen a numerical sweep before exporting STEP.

use catallax_features::general_sweep::{self, SweepRecipe, Swept};
use catallax_numbers::{Decider, with_operation_budget, with_operation_log};
use catallax_step::numerical::{self, NumericalStepOptions};
use serde::Deserialize;
use serde_json::json;
use std::{error::Error, fs, path::Path};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Input {
    format: String,
    definition: SweepRecipe,
}
fn read(path: &str) -> Result<SweepRecipe, Box<dyn Error>> {
    let input: Input = serde_json::from_slice(&fs::read(path)?)?;
    if input.format != "catallax-general-sweep/1" {
        return Err("expected catallax-general-sweep/1".into());
    }
    Ok(input.definition)
}
fn export(swept: &Swept, path: &Path, dec: &mut Decider) -> Result<Vec<u8>, Box<dyn Error>> {
    let artifact = numerical::build(swept.body(), &NumericalStepOptions::default(), dec)?;
    fs::write(path, artifact.bytes())?;
    let bytes = fs::read(path)?;
    artifact.check_written(&bytes, dec)?;
    Ok(bytes)
}
fn run() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 4 {
        return Err("usage: general_sweep INITIAL.json EDITED.json OUTPUT_DIRECTORY".into());
    }
    let mut dec = Decider::default();
    let initial_definition=read(&args[1])?;
    let initial=request("build",||Ok(general_sweep::build(&initial_definition,&mut dec)?))?;
    eprintln!("initial numerical body admitted");
    let out = Path::new(&args[3]);
    fs::create_dir_all(out)?;
    request("initial-export",||export(&initial,&out.join("initial.step"),&mut dec))?;
    eprintln!("initial STEP checked");
    let edited_definition=read(&args[2])?;
    let edited=request("edit",||Ok(initial.edit(edited_definition,&mut dec)?))?;
    eprintln!("edited numerical body admitted");
    let bytes=request("edited-export",||export(&edited,&out.join("edited.step"),&mut dec))?;
    let saved = json!({"format":"catallax-general-sweep/1","definition":edited.recipe()});
    fs::write(out.join("saved.json"), serde_json::to_vec_pretty(&saved)?)?;
    let reopened_definition=read(out.join("saved.json").to_str().ok_or("non-UTF8 output path")?)?;
    let reopened=request("reopen",||Ok(general_sweep::build(&reopened_definition,&mut dec)?))?;
    eprintln!("reopened numerical body admitted");
    let replay=request("reopened-export",||export(&reopened,&out.join("reopened.step"),&mut dec))?;
    if replay != bytes {
        return Err("reopened STEP differs from edited STEP".into());
    }
    let counts = edited.body().counts();
    let report = json!({"proof":"numerical","exact":false,"tolerance_mm":edited.recipe().tolerance_mm,
        "fitting":edited.checks(),"faces":counts.faces,"edges":counts.edges,"vertices":counts.vertices,
        "history":edited.body().history().len(),"reopened_step_identical":true,
        "complete_body_bound":"unestablished"});
    fs::write(out.join("report.json"), serde_json::to_vec_pretty(&report)?)?;
    println!("{report}");
    Ok(())
}
fn request<T>(id:&str,run:impl FnOnce()->Result<T,Box<dyn Error>>)->Result<T,Box<dyn Error>> {
    let ((result,work),operations)=with_operation_log(||with_operation_budget(1_000_000_000,run));
    eprintln!("{}",json!({"request":id,"work":work.used,"work_limit":1000000000,"exhausted":work.exhausted}));
    for operation in operations {
        eprintln!(
            "{}",
            json!({"operation":operation.operation,"work":operation.work_after-operation.work_before,
            "refusal":operation.refusal.as_ref().map(|r|r.code)})
        );
    }
    if work.exhausted {
        return Err(format!("{id} exhausted its explicit 1G request allowance").into());
    }
    result
}

fn main()->Result<(),Box<dyn Error>> {run()}
