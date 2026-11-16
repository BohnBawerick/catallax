//! Numerical top-stock Boolean delivery using the shared STEP embedding.

use std::{io, path::Path};

use catallax_boolean::numerical::loft::{loft_box_selection, stock_spans_loft_width};
use catallax_numbers::{BoundKind, Decider, Mark, OperationSpan, Refusal};
use catallax_recipe::{BooleanKind, Operation, Recipe};
use catallax_step::numerical::{self as step, NumericalStepOptions};
use serde_json::{Value, json};

fn failure(error: Refusal, stage: &str) -> Value {
    json!({"format":"catallax-manifest/5","status":"refused","exact":false,
        "result":null,"step":null,"marks":[],
        "refusal":{"code":error.code,"detail":error.detail,"step":null,"stage":stage}})
}

fn measured<T>(
    stage: &str,
    recipe: &Recipe,
    dec: &mut Decider,
    run: impl FnOnce(&mut Decider) -> Result<T, Refusal>,
) -> Result<T, Refusal> {
    let span = OperationSpan::begin("brep", &recipe.result, stage, dec);
    let result = run(dec);
    span.finish(dec, result.as_ref().err(), None);
    result
}

/// Handles an unchanged loft and a stock that spans its width at the final step.
pub(crate) fn evaluate(
    recipe: &Recipe,
    files: &catallax_features::ImportFiles<'_>,
    out: &Path,
    dec: &mut Decider,
) -> io::Result<Option<Value>> {
    let Some(last) = recipe.steps.last() else {
        return Ok(None);
    };
    let Operation::Boolean { kind, target, tool } = &last.op else {
        return Ok(None);
    };
    if last.id != recipe.result || !recipe.materials.is_empty() {
        return Ok(None);
    }
    let original_loft = recipe
        .steps
        .iter()
        .any(|step| step.id == *target && matches!(step.op, Operation::Loft { .. }));
    let original_stock = recipe
        .steps
        .iter()
        .any(|step| step.id == *tool && matches!(step.op, Operation::Block { .. }));
    if !original_loft || !original_stock {
        return Ok(None);
    }
    let mut prefix = recipe.clone();
    prefix.steps.pop();
    prefix.result = target.clone();
    let built = match catallax_features::evaluate_with_files(&prefix, files, dec) {
        Ok(built) => built,
        Err(catallax_features::BuildError::Geometry { step, refusal }) => {
            let mut value = failure(refusal, "source-build");
            value["refusal"]["step"] = json!(step.to_string());
            return Ok(Some(value));
        }
        Err(error) => {
            return Ok(Some(failure(
                Refusal::new("cli.invalid-recipe", error.to_string()),
                "source-build",
            )));
        }
    };
    let loft = built.bodies[target];
    let stock = built.bodies[tool];
    match measured("numerical-loft-routing", recipe, dec, |dec| {
        stock_spans_loft_width(&built.model, loft, stock, dec)
    }) {
        Ok(true) => {}
        Ok(false) => return Ok(None),
        Err(error) if error.code == "boolean.numerical-loft" => return Ok(None),
        Err(error) => return Ok(Some(failure(error, "width-routing"))),
    }
    let tolerance = catallax_step::StepOptions::default().curve_tolerance;
    let operation = match kind {
        BooleanKind::Union => catallax_boolean::BoolOp::Unite,
        BooleanKind::Subtract => catallax_boolean::BoolOp::Subtract,
        BooleanKind::Intersect => catallax_boolean::BoolOp::Intersect,
    };
    let span = OperationSpan::begin("brep", &recipe.result, "numerical-boolean", dec);
    let selected = loft_box_selection(&built.model, loft, stock, operation, tolerance.clone(), dec);
    span.finish(dec, selected.as_ref().err(), None);
    let selected = match selected {
        Ok(selected) => selected,
        Err(error) => return Ok(Some(failure(error, "numerical-boolean"))),
    };
    let options = NumericalStepOptions {
        max_scalar_rounding_error_mm: tolerance
            .clone()
            .min(NumericalStepOptions::default().max_scalar_rounding_error_mm),
        ..Default::default()
    };
    let mut artifact = None;
    let mut failed_stage = "numerical-step-definition";
    let span = OperationSpan::begin("brep", &recipe.result, "numerical-step", dec);
    let body = selected.admit_with(dec, |model, collection, admission, dec| {
        let (definition, native_trims) =
            measured("numerical-step-definition", recipe, dec, |dec| {
                step::definition_from_exact_model(
                    model,
                    collection,
                    recipe.name.clone(),
                    &admission.tolerance_mm,
                    dec,
                )
            })?;
        failed_stage = "numerical-native-admission";
        let body = measured("numerical-native-admission", recipe, dec, |dec| {
            catallax_topology::numerical::build_with_native_trims(
                definition,
                native_trims,
                admission,
                dec,
            )
        })?;
        failed_stage = "numerical-step-payload";
        let file = measured("numerical-step-payload", recipe, dec, |dec| {
            step::build(&body, &options, dec)
        })?;
        artifact = Some(file);
        Ok(body)
    });
    span.finish(dec, body.as_ref().err(), None);
    let body = match body {
        Ok(body) => body,
        Err(error) => return Ok(Some(failure(error, failed_stage))),
    };
    let step = match crate::numerical_step::publish_embedded(
        &artifact.expect("successful numerical embedding supplies its artifact"),
        recipe,
        out,
        dec,
    )? {
        Ok(step) => step,
        Err(error) => return Ok(Some(failure(error, "numerical-step-reread"))),
    };
    let counts = body.counts();
    let classes: Vec<_> = selected.classifications.iter().map(|class| json!({
        "face":class.face.index(),"from_loft":class.from_loft,"inside":class.inside,
        "witness_mm":class.witness_bounds.each_ref().map(|bounds| json!({"lo":bounds.lo().to_string(),"hi":bounds.hi().to_string()})),
        "evidence":class.evidence.iter().map(crate::numerical_evidence_summary).collect::<Vec<_>>()
    })).collect();
    let mark = Mark::new(
        "step.numerical-approximate",
        BoundKind::Value,
        tolerance.clone(),
        "numerically selected material and checked STEP embedding; complete body bound unestablished",
    );
    Ok(Some(
        json!({"format":"catallax-manifest/5","status":"answered","exact":false,
            "proof":"numerical","category":"approximate","refusal":null,
            "result":{"lumps":counts.lumps,"shells":counts.shells,"faces":counts.faces,"loops":counts.loops,"edges":counts.edges,"vertices":counts.vertices},
            "step":step,"marks":[crate::mark_entry(&mark,None)],
            "numerical":{"tolerance_mm":tolerance.to_string(),"material_classifications":classes,
                "work_stages":{
                    "source_admission":selected.work_stages[0],
                    "intersection":selected.work_stages[1],
                    "imprint":selected.work_stages[2],
                    "material_classification":selected.work_stages[3],
                    "face_collection":selected.work_stages[4]},"admission":crate::numerical_evidence_summary(body.evidence())}
        }),
    ))
}
