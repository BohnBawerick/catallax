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
    body: NumericalBody,
    origins: BTreeMap<String, Vec<Origin>>,
}
impl NumericalFilletEvaluation {
    /// Numerical geometry, native trims, declared tolerance and edit history.
    pub fn body(&self) -> &NumericalBody {
        &self.body
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
        (self.body, self.recipe, self.origins)
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
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let result = (|| {
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
            if selected.is_empty() {
                return Err(Refusal::new(
                    "features.numerical-fillet-selection",
                    "the selection contains no edges to round",
                ));
            }
            let faults = model.check(target.body);
            if !faults.is_empty() {
                return Err(Refusal::new(
                    "features.numerical-fillet-source",
                    format!("source topology is invalid: {faults:?}"),
                ));
            }
            let tolerance = Rational::from_ratio(1, 10_000_000).unwrap();
            let chains = crate::fillet::numerical::rolling::construct(
                model,
                target.body,
                &selected,
                radius,
                &tolerance,
                dec,
            )?;
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
            let mut source_origins: BTreeMap<_, _> = model
                .faces_of_body(target.body)
                .into_iter()
                .map(|face| (format!("{face:?}"), model.face_origins(face).to_vec()))
                .collect();
            for strip in chains.iter().flat_map(|chain| &chain.strips) {
                let edge = strip.source.edge;
                let name = model.edge_name(target.body, edge).ok_or_else(|| {
                    Refusal::new(
                        "features.numerical-fillet-name",
                        format!("rolled edge {edge:?} has no canonical recipe name"),
                    )
                })?;
                source_origins.insert(
                    format!("{step}/blend/{edge:?}"),
                    vec![Origin {
                        step: step.clone(),
                        role: Role::Blend { edge: name },
                    }],
                );
            }
            let body = if chains.len() == 1 && chains[0].strips.len() == 1 && !chains[0].closed {
                crate::fillet::numerical::body::cap_strip(
                    &source,
                    &chains[0].strips[0],
                    step.as_str(),
                    &tolerance,
                    dec,
                )?
            } else {
                let corners = crate::fillet::numerical::vertex::construct(
                    model, &chains, radius, &tolerance, dec,
                )?;
                for corner in &corners {
                    let vertex =
                        model
                            .vertex_name(target.body, corner.source)
                            .ok_or_else(|| {
                                Refusal::new(
                                    "features.numerical-fillet-name",
                                    format!(
                                        "corner {:?} has no canonical recipe name",
                                        corner.source
                                    ),
                                )
                            })?;
                    source_origins.insert(
                        format!("{step}/corner/{:?}", corner.source),
                        vec![Origin {
                            step: step.clone(),
                            role: Role::Corner { vertex },
                        }],
                    );
                }
                crate::fillet::numerical::vertex_body::build(
                    &source,
                    model,
                    &chains,
                    &corners,
                    step.as_str(),
                    &tolerance,
                    dec,
                )?
            };
            let origins = body
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
                body,
                origins,
            })
        })()?;
        self.bodies.insert(step.clone(), Arc::new(result));
        Ok(())
    }
}

/// Evaluates the matching exact prefix and constructs fillets numerically first.
/// No speculative exact producer spends or gates this construction allowance.
/// Explicit exact-only evaluation remains a separate caller-scheduled request;
/// it never converts or replaces this retained numerical geometry.
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
