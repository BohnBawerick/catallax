//! Numerical top-stock Boolean delivery using the shared STEP embedding.

use std::{io, path::Path};

use catallax_boolean::numerical::loft::{loft_box_selection, stock_overlaps_loft_top};
use catallax_numbers::{BoundKind, Decider, Mark, OperationSpan, Rational, Refusal};
use catallax_recipe::{BooleanKind, Operation, Recipe};
use catallax_step::numerical as step;
use catallax_topology::numerical::{self, NumericalOptions, tolerance};
use serde_json::{Value, json};

use crate::numerical_step::{embedding_step_options, exact_floor};

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

/// Handles an unchanged loft and an overlapping top stock at the final step.
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
    // Rule T2 for own fitted geometry (Checkpoint 2 decisions, section 4.3;
    // owner's choice A of 20 November 2026): the traced cut and its
    // embedding fit to max(1e-6 mm, 10 eps), and the body admits under
    // max(1e-5 mm, 100 eps). The floor is rule T1's in the frame of the
    // input bodies' exact vertices.
    let eps = match exact_floor(&built.model, [loft, stock], dec) {
        Ok(eps) => eps,
        Err(error) => return Ok(Some(failure(error, "width-routing"))),
    };
    let fitted = tolerance::own_fitted(&eps, None);
    let tolerance = fitted.target_mm;
    match measured("numerical-loft-routing", recipe, dec, |dec| {
        stock_overlaps_loft_top(&built.model, loft, stock, &tolerance, dec)
    }) {
        Ok(true) => {}
        Ok(false) => return Ok(None),
        Err(error) if error.code == "boolean.numerical-loft" => return Ok(None),
        Err(error) => return Ok(Some(failure(error, "width-routing"))),
    }
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
    let mut fit_target = tolerance.clone();
    let mut artifact = None;
    let mut failed_stage = "numerical-step-definition";
    let span = OperationSpan::begin("brep", &recipe.result, "numerical-step", dec);
    let cap = NumericalOptions {
        tolerance_mm: fitted.cap_mm,
        vertex_tolerance_mm: None,
    };
    let body = selected.admit_with(dec, |model, collection, _, dec| {
        let embed = |target: &Rational, dec: &mut Decider| {
            measured("numerical-step-definition", recipe, dec, |dec| {
                catallax_features::numerical_embedding::definition_from_exact_model(
                    model,
                    collection,
                    recipe.name.clone(),
                    target,
                    dec,
                )
            })
        };
        let (mut definition, mut native_trims) = embed(&tolerance, dec)?;
        failed_stage = "numerical-step-feature-search";
        fit_target = measured(failed_stage, recipe, dec, |dec| {
            tolerance::fitted_target_for_definition(&definition, &native_trims, &eps, dec)
        })?;
        if fit_target < tolerance {
            failed_stage = "numerical-step-definition";
            (definition, native_trims) = embed(&fit_target, dec)?;
        }
        failed_stage = "numerical-native-admission";
        let body = measured(failed_stage, recipe, dec, |dec| {
            numerical::build_with_native_trims(definition, native_trims, cap, dec)
        })?;
        failed_stage = "numerical-step-payload";
        let file = measured("numerical-step-payload", recipe, dec, |dec| {
            step::build(&body, &embedding_step_options(fit_target.clone()), dec)
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
        fit_target.clone(),
        "numerically selected material and checked STEP embedding; complete body bound unestablished",
    );
    Ok(Some(
        json!({"format":"catallax-manifest/5","status":"answered","exact":false,
            "proof":"numerical","category":"approximate","refusal":null,
            "result":{"lumps":counts.lumps,"shells":counts.shells,"faces":counts.faces,"loops":counts.loops,"edges":counts.edges,"vertices":counts.vertices},
            "step":step,"marks":[crate::mark_entry(&mark,None)],
            "numerical":{"tolerance_mm":fit_target.to_string(),"material_classifications":classes,
                "work_stages":{
                    "source_admission":selected.work_stages[0],
                    "intersection":selected.work_stages[1],
                    "imprint":selected.work_stages[2],
                    "material_classification":selected.work_stages[3],
                    "face_collection":selected.work_stages[4]},"admission":crate::numerical_evidence_summary(body.evidence())}
        }),
    ))
}
