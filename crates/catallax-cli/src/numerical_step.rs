//! Explicit numerical recipe admission and approximate STEP delivery.

use std::{fs, io, path::Path};

use catallax_numbers::{Decider, OperationSpan, Rational, Sign, current_operation_work};
use catallax_step::numerical::{self as step, NumericalStepOptions};
use catallax_topology::numerical::{self, NumericalOptions};
use serde::Deserialize;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

pub(crate) const FORMAT: &str = "catallax-numerical-step/1";
const MANIFEST: &str = "catallax-numerical-manifest/1";

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Problem {
    format: String,
    recipe: Value,
    tolerance_mm: String,
}

pub(crate) fn is_problem(value: &Value) -> bool {
    value["format"] == FORMAT
}

pub(crate) fn refused(code: &str, detail: impl ToString) -> Value {
    json!({"format":MANIFEST,"status":"refused","proof":"refused","exact":false,
        "category":"approximate","result":null,"numerical":null,"step":null,"marks":[],
        "refusal":{"code":code,"detail":detail.to_string(),"step":null,"stage":null}})
}

pub(crate) fn evaluate(value: &Value, source: &str, out: &Path) -> io::Result<Value> {
    let invalid = |detail: String| {
        let mut manifest = refused("cli.numerical-step-wrapper", detail);
        manifest["status"] = json!("invalid-question");
        Ok(manifest)
    };
    if source != "recipe" {
        return invalid("numerical STEP requires the recipe source".into());
    }
    let problem: Problem = match serde_json::from_value(value.clone()) {
        Ok(problem) => problem,
        Err(error) => return invalid(error.to_string()),
    };
    if problem.format != FORMAT {
        return invalid(format!("expected {FORMAT}"));
    }
    let tolerance = match Rational::parse(&problem.tolerance_mm) {
        Ok(tolerance) if tolerance.sign() == Sign::Positive => tolerance,
        _ => return invalid("tolerance_mm must be a positive rational string".into()),
    };
    let imported = match crate::import_recipe::parse(&problem.recipe) {
        Ok(imported) => imported,
        Err(error) => return Ok(refused(error.code, error.detail)),
    };
    let empty_files = catallax_features::ImportFiles::new();
    let (input, files) = imported
        .as_ref()
        .map_or((&problem.recipe, &empty_files), |(recipe, files)| {
            (*recipe, files)
        });
    let recipe = match catallax_recipe::Recipe::parse(&input.to_string()) {
        Ok(recipe) => recipe,
        Err(error) => {
            let mut manifest = refused("cli.invalid-recipe", error);
            manifest["status"] = json!("invalid-question");
            return Ok(manifest);
        }
    };
    let mut dec = Decider::default();
    let admitted = match admit(&recipe, files, &tolerance, &mut dec) {
        Ok(admitted) => admitted,
        Err(manifest) => return Ok(manifest),
    };
    let body = admitted.body;
    let recipe_work = admitted.recipe_work;
    let admitted_work = current_operation_work().used;
    let options = NumericalStepOptions {
        max_scalar_rounding_error_mm: tolerance
            .clone()
            .min(NumericalStepOptions::default().max_scalar_rounding_error_mm),
        ..Default::default()
    };
    let span = OperationSpan::begin("brep", &recipe.result, "numerical-step", &dec);
    let artifact = step::build(&body, &options, &mut dec);
    span.finish(&dec, artifact.as_ref().err(), None);
    let artifact = match artifact {
        Ok(artifact) => artifact,
        Err(error) => return Ok(refused(error.code, error.detail)),
    };
    let path = out.join("answer.step");
    let reread = match fs::write(&path, artifact.bytes()).and_then(|()| fs::read(&path)) {
        Ok(bytes) => bytes,
        Err(error) => {
            let _ = fs::remove_file(&path);
            return Err(error);
        }
    };
    let span = OperationSpan::begin("brep", &recipe.result, "numerical-step-reread", &dec);
    let checked = artifact.check_written(&reread, &mut dec);
    span.finish(&dec, checked.as_ref().err(), None);
    if let Err(error) = checked {
        return Ok(refused(error.code, error.detail));
    }
    let report = artifact.report();
    let spline = report.spline_rounding.as_ref().map(|rounding| {
        json!({
            "max_native_parameter_error":rounding.max_parameter_error.to_string(),
            "max_weight_error":rounding.max_weight_error.to_string()
        })
    });
    let counts = body.counts();
    Ok(
        json!({"format":MANIFEST,"status":"answered","proof":"numerical","exact":false,
            "category":"approximate","marks":[],"refusal":null,
            "result":{"lumps":counts.lumps,"shells":counts.shells,"faces":counts.faces,"edges":counts.edges,"vertices":counts.vertices},
            "numerical":{"source":crate::numerical_evidence_summary(&report.source_evidence),
                "emitted":crate::numerical_evidence_summary(&report.emitted_evidence),
                "complete_body_bound":"unestablished"},
            "step":{"file":"answer.step","category":"approximate","seams":"required",
                "sha256":format!("{:x}",Sha256::digest(&reread)),"bytes":reread.len(),
                "checks":report.checks,"tolerance_mm":tolerance.to_string(),
                "complete_body_bound":"unestablished",
                "scalar_rounding_limit_mm":report.scalar_rounding_limit_mm.to_string(),
                "max_scalar_rounding_error_mm":report.max_scalar_rounding_error_mm.to_string(),
                "max_direction_ratio_rounding_error":report.max_direction_ratio_rounding_error.to_string(),
                "spline_rounding":spline},
            "work_stages":{"recipe":recipe_work,"numerical_admission":admitted_work-recipe_work,
                "step_and_reread":current_operation_work().used-admitted_work}
        }),
    )
}

/// Numerical-first STEP for an exact recipe body whose geometry has no exact
/// STEP form (ROADMAP "Numerical first"). The tolerance is the exact writer's
/// requested curve tolerance, fixed before any embedding work. `None` keeps
/// the exact refusal when this embedding has no form for a carrier of the
/// body. A refusal leaves no STEP file.
pub(crate) fn embedded_exact(
    model: &catallax_topology::Model,
    body: catallax_topology::BodyId,
    recipe: &catallax_recipe::Recipe,
    out: &Path,
    dec: &mut Decider,
) -> io::Result<Option<Result<Value, catallax_numbers::Refusal>>> {
    if !step::embeds(model, body) {
        return Ok(None);
    }
    let tolerance = catallax_step::StepOptions::default().curve_tolerance;
    let options = NumericalStepOptions {
        max_scalar_rounding_error_mm: tolerance
            .clone()
            .min(NumericalStepOptions::default().max_scalar_rounding_error_mm),
        ..Default::default()
    };
    let span = OperationSpan::begin("brep", &recipe.result, "numerical-step", dec);
    let built = step::build_from_exact_model(
        model,
        body,
        recipe.name.clone(),
        NumericalOptions {
            tolerance_mm: tolerance.clone(),
        },
        &options,
        dec,
    );
    span.finish(dec, built.as_ref().err(), None);
    let artifact = match built {
        Ok((_, artifact)) => artifact,
        Err(error) => return Ok(Some(Err(error))),
    };
    let path = out.join("answer.step");
    let reread = match fs::write(&path, artifact.bytes()).and_then(|()| fs::read(&path)) {
        Ok(bytes) => bytes,
        Err(error) => {
            let _ = fs::remove_file(&path);
            return Err(error);
        }
    };
    let span = OperationSpan::begin("brep", &recipe.result, "numerical-step-reread", dec);
    let checked = artifact.check_written(&reread, dec);
    span.finish(dec, checked.as_ref().err(), None);
    if let Err(error) = checked {
        let _ = fs::remove_file(&path);
        return Ok(Some(Err(error)));
    }
    let report = artifact.report();
    let spline = report.spline_rounding.as_ref().map(|rounding| {
        json!({
            "max_native_parameter_error":rounding.max_parameter_error.to_string(),
            "max_weight_error":rounding.max_weight_error.to_string()
        })
    });
    Ok(Some(Ok(json!({
        "file":"answer.step","approximate":true,"category":"approximate","proof":"numerical",
        "seams":"required","sha256":format!("{:x}",Sha256::digest(&reread)),"bytes":reread.len(),
        "checks":report.checks,"tolerance_mm":tolerance.to_string(),
        "complete_body_bound":"unestablished",
        "numerical":{"source":crate::numerical_evidence_summary(&report.source_evidence),
            "emitted":crate::numerical_evidence_summary(&report.emitted_evidence)},
        "native_trims":report.source_native_trims.len(),
        "scalar_rounding_limit_mm":report.scalar_rounding_limit_mm.to_string(),
        "max_scalar_rounding_error_mm":report.max_scalar_rounding_error_mm.to_string(),
        "max_direction_ratio_rounding_error":report.max_direction_ratio_rounding_error.to_string(),
        "spline_rounding":spline
    }))))
}

struct Admitted {
    body: numerical::NumericalBody,
    recipe_work: u64,
}
fn admit(
    recipe: &catallax_recipe::Recipe,
    files: &catallax_features::ImportFiles<'_>,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<Admitted, Value> {
    if files.is_empty() {
        match catallax_features::evaluate_numerical_hole(
            recipe,
            NumericalOptions {
                tolerance_mm: tolerance.clone(),
            },
            dec,
        ) {
            Ok(Some(built)) => {
                let assembly = built.work_stages()[1];
                let recipe_work = current_operation_work().used - assembly;
                return Ok(Admitted {
                    body: built.into_hole().into_body(),
                    recipe_work,
                });
            }
            Ok(None) => {}
            Err(catallax_features::BuildError::Geometry { step, refusal }) => {
                let mut manifest = refused(refusal.code, refusal.detail);
                manifest["refusal"]["step"] = json!(step.to_string());
                return Err(manifest);
            }
            Err(error) => return Err(refused("cli.invalid-recipe", error)),
        }
    }
    let built = match catallax_features::evaluate_with_files(recipe, files, dec) {
        Ok(built) => built,
        Err(catallax_features::BuildError::Geometry { step, refusal }) => {
            let mut manifest = refused(refusal.code, refusal.detail);
            manifest["refusal"]["step"] = json!(step.to_string());
            return Err(manifest);
        }
        Err(error) => return Err(refused("cli.invalid-recipe", error)),
    };
    if !built.stand_ins.is_empty() {
        return Err(refused(
            "field.stand-in-export",
            "numerical STEP requires retained boundary geometry",
        ));
    }
    let recipe_work = current_operation_work().used;
    let span = OperationSpan::begin("brep", &recipe.result, "numerical-admission", dec);
    let body =
        numerical::definition_from_model(&built.model, built.result, recipe.name.clone(), dec)
            .and_then(|definition| {
                numerical::build(
                    definition,
                    NumericalOptions {
                        tolerance_mm: tolerance.clone(),
                    },
                    dec,
                )
            });
    span.finish(dec, body.as_ref().err(), None);
    let body = match body {
        Ok(body) => body,
        Err(error) => return Err(refused(error.code, error.detail)),
    };
    Ok(Admitted { body, recipe_work })
}
