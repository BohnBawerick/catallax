//! Numerical fillet retention without an exact checkpoint or model conversion.

use super::*;
use catallax_recipe::{EdgeSelection, FilletCorner, FilletOverflow, Origin, Step};
use catallax_topology::numerical::{NumericalBody, NumericalOptions};

/// An exact continuation checkpoint or a separately admitted numerical result.
#[derive(Clone, Debug)]
pub enum EvaluationOutcome {
    /// Every retained body and continuation state remains on the exact path.
    Exact(EvaluationCheckpoint),
    /// Scoped numerical checks, never an exact continuation checkpoint.
    Numerical(Arc<NumericalFilletEvaluation>),
}

/// Immutable authored controls, numerical geometry and canonical face origins.
#[derive(Clone, Debug)]
pub struct NumericalFilletEvaluation {
    recipe: Arc<Recipe>,
    fillet: crate::fillet::numerical::NumericalFilletBuilt,
    origins: BTreeMap<String, Vec<Origin>>,
}
impl NumericalFilletEvaluation {
    /// Numerical geometry, native trims, declared tolerance and edit history.
    pub fn body(&self) -> &NumericalBody {
        self.fillet.body()
    }
    /// Authored controls retained for replay and parameter edits.
    pub fn recipe(&self) -> &Recipe {
        &self.recipe
    }
    /// Canonical origins keyed by the numerical definition's stable face keys.
    pub fn face_origins(&self) -> &BTreeMap<String, Vec<Origin>> {
        &self.origins
    }
    /// Consumes the result into its admitted geometry, authored controls and origins.
    pub fn into_parts(self) -> (NumericalBody, Arc<Recipe>, BTreeMap<String, Vec<Origin>>) {
        (self.fillet.into_body(), self.recipe, self.origins)
    }
}

#[derive(Default)]
pub(super) struct Context {
    bodies: BTreeMap<StepId, Arc<NumericalFilletEvaluation>>,
}
impl Context {
    pub(super) fn exact_inputs(&self, step: &Step) -> Result<(), Refusal> {
        if self.bodies.is_empty() {
            return Ok(());
        }
        if let Some((input, _)) = step
            .op
            .inputs()
            .into_iter()
            .find(|(id, _)| self.bodies.contains_key(*id))
        {
            return Err(Refusal::new(
                "features.numerical-continuation",
                format!(
                    "operation {} reads numerical body {input}; its native numerical continuation has no admitted construction for this operation",
                    step.op.name()
                ),
            ));
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn fillet(
        &mut self,
        model: &Model,
        target: &BodyValue,
        selection: &EdgeSelection,
        radius: &Rational,
        overflow: Option<&FilletOverflow>,
        corner: Option<FilletCorner>,
        earlier: &[crate::fillet::FilletDefinition],
        step: &StepId,
        recipe: &Recipe,
        exact_refusal: Refusal,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        if exact_refusal.code.starts_with("numbers.")
            || radius.sign() != catallax_numbers::Sign::Positive
        {
            return Err(exact_refusal);
        }
        let mut attempt = || {
            if overflow.is_some() || corner.is_some_and(|value| value != FilletCorner::Cliff) {
                return Err(Refusal::new(
                    "features.numerical-fillet-closure",
                    "the numerical spline assembly cannot honour explicit neighbour limits or a non-cliff corner closure",
                ));
            }
            let selected = crate::fillet_build::select(
                model,
                target.body,
                selection,
                recipe.format != catallax_recipe::FORMAT,
                earlier,
                dec,
            )?;
            let [edge] = selected.as_slice() else {
                return Err(Refusal::new(
                    "features.numerical-fillet-network",
                    format!(
                        "the retained spline assembly needs one outer-boundary edge, but the selection resolves to {} edges",
                        selected.len()
                    ),
                ));
            };
            let name = model.edge_name(target.body, *edge).ok_or_else(|| {
                Refusal::new(
                    "features.numerical-fillet-name",
                    "the selected edge has no canonical recipe name",
                )
            })?;
            let faults = model.check(target.body);
            if !faults.is_empty() {
                return Err(Refusal::new(
                    "features.numerical-fillet-source",
                    format!("source topology is invalid: {faults:?}"),
                ));
            }
            let tolerance = Rational::from_ratio(1, 10_000_000).unwrap();
            let (definition, trims) = crate::numerical_embedding::definition_from_exact_model(
                model,
                target.body,
                recipe.name.clone(),
                &tolerance,
                dec,
            )?;
            let source = catallax_topology::numerical::build_with_native_trims(
                definition,
                trims,
                NumericalOptions {
                    tolerance_mm: tolerance.clone(),
                },
                dec,
            )?;
            let edge_key = format!("{edge:?}");
            let fillet = crate::fillet::numerical::build(
                source,
                &edge_key,
                step.as_str(),
                radius.clone(),
                tolerance,
                dec,
            )?;
            let mut source_origins: BTreeMap<_, _> = model
                .faces_of_body(target.body)
                .into_iter()
                .map(|face| (format!("{face:?}"), model.face_origins(face).to_vec()))
                .collect();
            source_origins.insert(
                format!("{step}/blend/{edge_key}"),
                vec![Origin {
                    step: step.clone(),
                    role: Role::Blend { edge: name },
                }],
            );
            let origins = fillet
                .body()
                .definition()
                .faces
                .iter()
                .map(|face| {
                    source_origins
                        .remove(&face.key)
                        .map(|origins| (face.key.clone(), origins))
                        .ok_or_else(|| {
                            Refusal::new(
                                "features.numerical-fillet-name",
                                format!("output face {} lost its canonical origins", face.key),
                            )
                        })
                })
                .collect::<Result<_, _>>()?;
            Ok(NumericalFilletEvaluation {
                recipe: Arc::new(recipe.clone()),
                fillet,
                origins,
            })
        };
        let result = attempt().map_err(|mut refusal: Refusal| {
            refusal.detail = format!(
                "exact path refused {}: {}; numerical continuation refused {}: {}",
                exact_refusal.code, exact_refusal.detail, refusal.code, refusal.detail
            );
            refusal
        })?;
        self.bodies.insert(step.clone(), Arc::new(result));
        Ok(())
    }
}

/// Evaluates a matching exact prefix once, then tries exact fillets before
/// numerical continuation in the same active request allowance. Numerical
/// results never enter an exact checkpoint. Files retain their ordinary checks.
pub fn evaluate_incremental_outcome(
    recipe: &Recipe,
    prefix: Option<EvaluationCheckpoint>,
    files: &ImportFiles<'_>,
    dec: &mut Decider,
) -> Result<EvaluationOutcome, BuildError> {
    let (state, start) = resume_state(recipe, prefix, dec);
    let mut context = Context::default();
    let state = evaluate_impl(
        recipe,
        files,
        None,
        None,
        state,
        start,
        Some(&mut context),
        dec,
    )?;
    if let Some(result) = context.bodies.remove(&recipe.result) {
        return Ok(EvaluationOutcome::Numerical(result));
    }
    if !context.bodies.is_empty() {
        return Err(BuildError::Geometry {
            step: recipe.result.clone(),
            refusal: Refusal::new(
                "features.numerical-mixed-checkpoint",
                "an exact result with numerical intermediate branches cannot retain an exact-only checkpoint",
            ),
        });
    }
    Ok(EvaluationOutcome::Exact(checkpoint(state, recipe, dec)))
}
