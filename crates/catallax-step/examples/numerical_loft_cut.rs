//! Reproduce two retained loft cuts and their numerical STEP evidence.

use catallax_numbers::{
    Decider, NumericalCheck, NumericalEvidence, Rational, current_operation_work,
    with_operation_budget,
};
use catallax_recipe::Recipe;
use catallax_step::numerical::{self as step, NumericalStepOptions};
use catallax_topology::{
    Sense,
    numerical::{self, NumericalOptions},
};
use serde_json::{Value, json};
use std::{path::Path, process::ExitCode};

const WORK: u64 = 20_000_000;

fn evidence(source: &NumericalEvidence) -> Value {
    let checks: Vec<_> = source.checks().iter().map(|check| match check {
        NumericalCheck::Passed {name,scope} => json!({"kind":"passed","name":name,"scope":scope}),
        NumericalCheck::Samples {name,scope,samples} => json!({"kind":"samples","name":name,"scope":scope,"samples":samples.iter().map(|sample| json!({"parameters":sample.parameters.iter().map(ToString::to_string).collect::<Vec<_>>(),"distance_upper_mm":sample.distance_upper_mm.to_string()})).collect::<Vec<_>>()}),
    }).collect();
    json!({"tolerance_mm":source.tolerance_mm().to_string(),"checks":checks,"repair_count":source.repairs().len()})
}

fn run(out: &Path, case: &str, control: &str) -> Result<Value, String> {
    let (name, size, at) = match case {
        "tiny" => (
            "loft crossing",
            ["10", "10", "1"],
            ["-3", "1999999/1000000", "1"],
        ),
        "m10-008" => ("m7-cut", ["2", "4", "1"], ["1", "1", "1"]),
        _ => return Err("unknown part".into()),
    };
    let mut input: Value = serde_json::from_str(include_str!(
        "../../catallax-features/tests/fixtures/loft-spline-prism.json"
    ))
    .map_err(|e| e.to_string())?;
    input["name"] = json!(name);
    let steps = input["steps"].as_array_mut().unwrap();
    steps.push(json!({"id":"box","op":"block","size":size,"pose":{"rotation":["1","0","0","0"],"translation":at}}));
    steps.push(json!({"id":"cut","op":"boolean","kind":"subtract","target":"loft","tool":"box"}));
    input["result"] = json!("cut");
    std::fs::write(
        out.join("recipe.json"),
        serde_json::to_vec_pretty(&input).unwrap(),
    )
    .map_err(|e| e.to_string())?;
    let recipe = Recipe::parse(&input.to_string()).map_err(|e| format!("{e:?}"))?;
    let mut dec = Decider::default();
    let built = catallax_features::evaluate(&recipe, &mut dec).map_err(|e| format!("{e:?}"))?;
    let recipe_work = current_operation_work().used;
    let mut definition =
        numerical::definition_from_model(&built.model, built.result, name.into(), &mut dec)
            .map_err(|e| format!("{e:?}"))?;
    if control != "source" {
        let face = definition
            .faces
            .iter_mut()
            .find(|face| face.loops.len() == 2)
            .ok_or("part has no physical hole")?;
        match control {
            "reverse-hole" => {
                let inner = &mut face.loops[1];
                inner.reverse();
                for trim in inner {
                    trim.uv.reverse();
                    trim.sense = if trim.sense == Sense::Same {
                        Sense::Reversed
                    } else {
                        Sense::Same
                    };
                }
            }
            "missing-hole" => {
                face.loops.pop();
            }
            _ => return Err("unknown control".into()),
        }
    }
    let body = numerical::build(
        definition,
        NumericalOptions {
            tolerance_mm: Rational::parse("1/10000000").unwrap(),
        },
        &mut dec,
    )
    .map_err(|e| format!("{e:?}"))?;
    let admitted_work = current_operation_work().used;
    let file = step::build(&body, &NumericalStepOptions::default(), &mut dec)
        .map_err(|e| format!("{e:?}"))?;
    std::fs::write(out.join("answer.step"), file.bytes()).map_err(|e| e.to_string())?;
    let reread = std::fs::read(out.join("answer.step")).map_err(|e| e.to_string())?;
    file.check_written(&reread, &mut dec)
        .map_err(|e| format!("{e:?}"))?;
    let report = file.report();
    let rounding = report.spline_rounding.as_ref().unwrap();
    Ok(json!({
        "status":"numerical","category":"approximate","case":case,"control":control,
        "counts":report.counts,"two_loop_faces":report.source_definition.faces.iter().filter(|face| face.loops.len()==2).count(),
        "bytes":file.bytes().len(),"source_evidence":evidence(&report.source_evidence),"emitted_evidence":evidence(&report.emitted_evidence),
        "file_checks":report.checks,"max_scalar_rounding_error_mm":report.max_scalar_rounding_error_mm.to_string(),
        "max_direction_ratio_rounding_error":report.max_direction_ratio_rounding_error.to_string(),
        "max_native_parameter_error":rounding.max_parameter_error.to_string(),"max_weight_error":rounding.max_weight_error.to_string(),
        "complete_body_bound":"unestablished",
        "work_stages":{"recipe":recipe_work,"numerical_admission":admitted_work-recipe_work,"step_and_reread":current_operation_work().used-admitted_work}
    }))
}

fn main() -> ExitCode {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if !(2..=3).contains(&args.len())
        || !matches!(args[1].as_str(), "tiny" | "m10-008")
        || args
            .get(2)
            .is_some_and(|s| !matches!(s.as_str(), "source" | "reverse-hole" | "missing-hole"))
    {
        eprintln!(
            "usage: numerical_loft_cut OUTPUT_DIRECTORY tiny|m10-008 [source|reverse-hole|missing-hole]"
        );
        return ExitCode::from(2);
    }
    let out = Path::new(&args[0]);
    let control = args.get(2).map_or("source", String::as_str);
    if let Err(error) = std::fs::create_dir_all(out) {
        eprintln!("{error}");
        return ExitCode::FAILURE;
    }
    if let Err(error) = std::fs::remove_file(out.join("answer.step"))
        && error.kind() != std::io::ErrorKind::NotFound
    {
        eprintln!("{error}");
        return ExitCode::FAILURE;
    }
    let (result, work) = with_operation_budget(WORK, || run(out, &args[1], control));
    let mut response = match result {
        Ok(value) if !work.exhausted => value,
        other => {
            let _ = std::fs::remove_file(out.join("answer.step"));
            let reason = match other {
                Err(reason) => reason,
                Ok(_) => "operation work allowance exhausted".into(),
            };
            json!({"status":"refused","case":args[1],"control":control,"reason":reason})
        }
    };
    response["work"] =
        json!({"used":work.used,"limit":WORK,"target":WORK,"exhausted":work.exhausted});
    if let Err(error) = std::fs::write(
        out.join("response.json"),
        serde_json::to_vec_pretty(&response).unwrap(),
    ) {
        eprintln!("{error}");
        return ExitCode::FAILURE;
    }
    println!(
        "{}",
        json!({"status":response["status"],"case":args[1],"control":control,"work":response["work"]})
    );
    if response["status"] == "numerical" {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}
