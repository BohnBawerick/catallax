//! Ordered feature evaluation, with bounds derived from the recipe definitions.

use crate::{Profile, revolve};
use catallax_boolean::{BoolOp, boolean};
use catallax_measure::BoundingBox;
use catallax_numbers::{Ball, BudgetPolicy, Decider, Rational, Real, Refusal};
use catallax_recipe::{BooleanKind, Operation, Recipe, RecipeError, Role, StepId, ThroughSide};
use catallax_shapes::{Isometry, RPoint3, RVector3};
use catallax_topology::{BodyId, BodyKind, Model};
use im::OrdMap;
use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
    sync::Arc,
};

#[path = "evaluate_numerical_fillet.rs"]
mod numerical_fillet;
#[path = "evaluate_step.rs"]
mod step;
pub use numerical_fillet::{
    EvaluationOutcome, NumericalFilletEvaluation, evaluate_incremental_outcome,
};

/// Immutable STEP, IGES and mesh source bytes, keyed by the SHA-256 saved in the recipe.
/// File names are labels, not paths the evaluator opens.
pub type ImportFiles<'a> = BTreeMap<String, &'a [u8]>;

/// A complete build. Intermediate input bodies remain immutable and queryable.
#[derive(Clone, Debug, catallax_numbers::Transfer)]
pub struct Built {
    /// Geometry, topology and canonical face origin sets.
    pub model: Model,
    /// The requested result.
    pub result: BodyId,
    /// Every body-producing step, in canonical step-id order.
    pub bodies: BTreeMap<StepId, BodyId>,
    /// Generating definitions of every elementary fillet and named copy in build order.
    pub fillets: Vec<crate::fillet::FilletDefinition>,
    /// Canonical removed-face history from feature operations, not selectable names.
    pub consumed: Vec<crate::ConsumedFace>,
    /// Generating definitions of every elementary chamfer and named copy.
    pub chamfers: Vec<crate::chamfer::ChamferDefinition>,
    /// Field-only steps the B-rep carries as their region envelope, in build order.
    /// The field kernel supplies the matching mass bodies; see [`crate::stand_in`].
    pub stand_ins: Vec<StandIn>,
    pub(crate) definition: Arc<Recipe>,
    import_records: BTreeMap<StepId, catallax_import::RepairRecord>,
    import_solids: BTreeMap<StepId, Vec<crate::import_edit::ImportedSolid>>,
}

/// A field-only step whose B-rep body is the envelope of its region.
#[derive(Clone, PartialEq, Eq, Debug, catallax_numbers::Transfer)]
pub struct StandIn {
    /// The field-only step.
    pub step: StepId,
    /// The region step whose body stands in for it.
    pub region: StepId,
}

impl Built {
    /// Source-bound repair evidence for an immutable import step. Distances are
    /// in source coordinates; rigid placement preserves them. This is not a
    /// certificate for subsequent booleans, display meshes or exported files.
    pub fn import_record(&self, step: &StepId) -> Option<&catallax_import::RepairRecord> {
        self.import_records.get(step)
    }

    /// Source solid entity keys and their placed lumps, in source order. These
    /// identify the immutable import step, not the result of later features.
    pub fn import_solids(&self, step: &StepId) -> Option<&[crate::import_edit::ImportedSolid]> {
        self.import_solids.get(step).map(Vec::as_slice)
    }

    /// The immutable recipe definition behind the analytic reductions and trims.
    /// Origins refer to this snapshot, not a subsequently edited input recipe.
    pub fn recipe(&self) -> &Recipe {
        &self.definition
    }
}

/// A structural recipe error or the first geometric failure in step order.
#[derive(Clone, Debug)]
pub enum BuildError {
    /// The recipe cannot be evaluated because its structure is invalid.
    Recipe(RecipeError),
    /// A geometric operation failed; no partial build is returned.
    Geometry {
        /// The first failing step.
        step: StepId,
        /// The stable refusal and its explanation.
        refusal: Refusal,
    },
}
impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Recipe(e) => e.fmt(f),
            Self::Geometry { step, refusal } => {
                write!(f, "step {step}: {}: {}", refusal.code, refusal.detail)
            }
        }
    }
}
impl std::error::Error for BuildError {}

/// An immutable evaluated prefix. Cloning shares its state in constant time.
///
/// Keep one per body for undo, branches and concurrent queries. Continuing a
/// clone copies only modified store paths, leaving saved models unchanged.
#[derive(Clone, Debug)]
pub struct EvaluationCheckpoint {
    built: Arc<Built>,
    profiles: OrdMap<StepId, Profile>,
    paths: OrdMap<StepId, crate::sweep::Path>,
    bodies: OrdMap<StepId, BodyValue>,
    meshes: OrdMap<StepId, crate::mesh_ops::MeshValue>,
    policy: BudgetPolicy,
}

impl EvaluationCheckpoint {
    /// The exact body, intermediate bodies and definitions at this checkpoint.
    pub fn built(&self) -> &Built {
        &self.built
    }

    /// Takes the build, sharing model stores if other checkpoints still retain it.
    pub fn into_built(self) -> Built {
        Arc::unwrap_or_clone(self.built)
    }

    fn into_state(self) -> EvaluationState {
        let built = Arc::unwrap_or_clone(self.built);
        EvaluationState {
            model: built.model,
            profiles: self.profiles,
            paths: self.paths,
            bodies: self.bodies,
            meshes: self.meshes,
            fillets: built.fillets,
            consumed: built.consumed,
            chamfers: built.chamfers,
            stand_ins: built.stand_ins,
            import_records: built.import_records,
            import_solids: built.import_solids,
        }
    }
}

#[derive(Default)]
struct EvaluationState {
    model: Model,
    profiles: OrdMap<StepId, Profile>,
    paths: OrdMap<StepId, crate::sweep::Path>,
    bodies: OrdMap<StepId, BodyValue>,
    /// The exact triangles of every mesh step, RECIPE.md sections 4.19 to 4.22.
    meshes: OrdMap<StepId, crate::mesh_ops::MeshValue>,
    fillets: Vec<crate::fillet::FilletDefinition>,
    consumed: Vec<crate::ConsumedFace>,
    chamfers: Vec<crate::chamfer::ChamferDefinition>,
    stand_ins: Vec<StandIn>,
    import_records: BTreeMap<StepId, catallax_import::RepairRecord>,
    import_solids: BTreeMap<StepId, Vec<crate::import_edit::ImportedSolid>>,
}

impl EvaluationState {
    fn finish(self, definition: Arc<Recipe>) -> Built {
        Built {
            result: self.bodies[&definition.result].body,
            model: self.model,
            bodies: self
                .bodies
                .into_iter()
                .map(|(id, value)| (id, value.body))
                .collect(),
            fillets: self.fillets,
            consumed: self.consumed,
            chamfers: self.chamfers,
            stand_ins: self.stand_ins,
            import_records: self.import_records,
            import_solids: self.import_solids,
            definition,
        }
    }
}
#[derive(Clone, Debug)]
struct BodyValue {
    body: BodyId,
    // A certified bounding sphere about model-space zero, not a fitted box.
    radius: Rational,
}

fn norm_upper(p: &RPoint3) -> Rational {
    p.x.abs() + p.y.abs() + p.z.abs()
}
pub(crate) fn profile_radius(profile: &Profile) -> Rational {
    let zero = RVector3::zero();
    profile
        .loops
        .iter()
        .flat_map(|b| &b.segments)
        .map(|s| {
            if let Some(controls) = &s.bezier {
                controls
                    .iter()
                    .map(|p| norm_upper(&profile.world(p, &zero)))
                    .max()
                    .expect("Bezier controls")
            } else if let Some((c, radius_sq, _)) = &s.circle {
                let radius = Ball::from_rational(radius_sq, 64)
                    .sqrt()
                    .expect("positive circle radius")
                    .hi();
                norm_upper(&profile.world(c, &zero))
                    + radius
                        * profile
                            .ellipse_frame
                            .as_ref()
                            .map_or_else(Rational::one, |f| {
                                let a = f.axis[0].abs() + f.axis[1].abs();
                                a * Rational::one().max(f.ratio.clone())
                            })
            } else {
                norm_upper(&profile.world(&s.start, &zero))
                    .max(norm_upper(&profile.world(&s.end, &zero)))
            }
        })
        .max()
        .unwrap_or_else(Rational::zero)
}

/// Builds stock solids, sketches, extrusions, analytic line/arc sweeps, revolutions,
/// elementary fillets, cuts, booleans, patterns and mirrors.
/// Other recipe operations retain their unfinished-operation refusal, not a partial body.
pub fn evaluate(recipe: &Recipe, dec: &mut Decider) -> Result<Built, BuildError> {
    evaluate_impl(
        recipe,
        &ImportFiles::new(),
        None,
        None,
        EvaluationState::default(),
        0,
        None,
        dec,
    )
    .map(|state| state.finish(Arc::new(recipe.clone())))
}

/// A separate numerical Hole result with its original recipe and named origins.
#[derive(Clone, Debug)]
pub struct NumericalHoleBuilt {
    recipe: Recipe,
    hole: crate::hole::NumericalHole,
    origins: BTreeMap<String, Vec<catallax_recipe::Origin>>,
    prefix_work: u64,
    assembly_work: u64,
}
impl NumericalHoleBuilt {
    /// Original authored controls, retained for parameter edits and replay.
    pub fn recipe(&self) -> &Recipe {
        &self.recipe
    }
    /// Numerical geometry, native curves, tolerance and edit history.
    pub fn hole(&self) -> &crate::hole::NumericalHole {
        &self.hole
    }
    /// Moves the numerical result into a numerical delivery path.
    pub fn into_hole(self) -> crate::hole::NumericalHole {
        self.hole
    }
    /// Original stock origins and compiled Hole-band origins by output face key.
    pub fn face_origins(&self) -> &BTreeMap<String, Vec<catallax_recipe::Origin>> {
        &self.origins
    }
    /// Actual work for the exact input prefix, then numerical assembly/admission.
    pub fn work_stages(&self) -> [u64; 2] {
        [self.prefix_work, self.assembly_work]
    }
    /// Re-evaluates edited controls from their immutable source definitions.
    pub fn edit(
        &self,
        recipe: &Recipe,
        options: catallax_topology::numerical::NumericalOptions,
        dec: &mut Decider,
    ) -> Result<Option<Self>, BuildError> {
        numerical_hole_request(recipe, Some(self), options, dec)
    }
}

/// Opts into a numerical continuation for one terminal through Hole placement.
/// Other input families return `None`. The exact evaluator and its model type
/// remain separate. Numerical STEP and numerical handover consume this result.
pub fn evaluate_numerical_hole(
    recipe: &Recipe,
    options: catallax_topology::numerical::NumericalOptions,
    dec: &mut Decider,
) -> Result<Option<NumericalHoleBuilt>, BuildError> {
    numerical_hole_request(recipe, None, options, dec)
}
/// A standard recipe build or a separate numerical Hole result.
#[derive(Clone, Debug)]
pub enum HoleEvaluation {
    /// Result of the ordinary feature evaluator, retaining its marks and checks.
    Standard(Box<Built>),
    /// Checked numerical Hole result with its authored controls and origins.
    Numerical(Box<NumericalHoleBuilt>),
}

/// Tries the numerical Hole continuation and resumes the ordinary evaluator
/// from its retained prefix when that stock family is unsupported. One request
/// allowance covers both the prefix and the selected suffix.
pub fn evaluate_with_numerical_hole(
    recipe: &Recipe,
    options: catallax_topology::numerical::NumericalOptions,
    dec: &mut Decider,
) -> Result<HoleEvaluation, BuildError> {
    hole_request(recipe, || {
        match numerical_hole_impl(recipe, None, options, dec)? {
            HoleAttempt::Numerical(built) => Ok(HoleEvaluation::Numerical(built)),
            HoleAttempt::Standard { state, start } => evaluate_impl(
                recipe,
                &ImportFiles::new(),
                None,
                None,
                *state,
                start,
                None,
                dec,
            )
            .map(|state| {
                HoleEvaluation::Standard(Box::new(state.finish(Arc::new(recipe.clone()))))
            }),
        }
    })
}

enum HoleAttempt {
    Numerical(Box<NumericalHoleBuilt>),
    Standard {
        state: Box<EvaluationState>,
        start: usize,
    },
}
impl HoleAttempt {
    fn standard() -> Self {
        Self::Standard {
            state: Box::default(),
            start: 0,
        }
    }
    fn into_numerical(self) -> Option<NumericalHoleBuilt> {
        match self {
            Self::Numerical(built) => Some(*built),
            Self::Standard { .. } => None,
        }
    }
}
fn numerical_hole_request(
    recipe: &Recipe,
    previous: Option<&NumericalHoleBuilt>,
    options: catallax_topology::numerical::NumericalOptions,
    dec: &mut Decider,
) -> Result<Option<NumericalHoleBuilt>, BuildError> {
    hole_request(recipe, || {
        numerical_hole_impl(recipe, previous, options, dec).map(HoleAttempt::into_numerical)
    })
}
fn hole_request<T>(
    recipe: &Recipe,
    run: impl FnOnce() -> Result<T, BuildError>,
) -> Result<T, BuildError> {
    if catallax_numbers::operation_budget_active() {
        return run();
    }
    let (result, work) =
        catallax_numbers::with_operation_budget(catallax_topology::numerical::DEFAULT_WORK, run);
    if work.exhausted {
        Err(BuildError::Geometry {
            step: recipe.result.clone(),
            refusal: Refusal::new(
                "numbers.operation-budget",
                "numerical Hole recipe exhausted its operation allowance",
            ),
        })
    } else {
        result
    }
}
fn numerical_hole_impl(
    recipe: &Recipe,
    previous: Option<&NumericalHoleBuilt>,
    options: catallax_topology::numerical::NumericalOptions,
    dec: &mut Decider,
) -> Result<HoleAttempt, BuildError> {
    recipe.check().map_err(BuildError::Recipe)?;
    let Some(step) = recipe.steps.last().filter(|step| step.id == recipe.result) else {
        return Ok(HoleAttempt::standard());
    };
    let Operation::Hole {
        target,
        at,
        diameter,
        depth,
        tip,
        head,
    } = &step.op
    else {
        return Ok(HoleAttempt::standard());
    };
    if at.len() != 1
        || !matches!(depth, catallax_recipe::HoleDepth::Through)
        || !matches!(tip, catallax_recipe::HoleTip::Flat)
        || head.is_some()
        || !recipe.materials.is_empty()
    {
        return Ok(HoleAttempt::standard());
    }
    let mut prefix = recipe.clone();
    prefix.steps.pop();
    prefix.result = target.clone();
    let before = catallax_numbers::current_operation_work().used;
    let state = evaluate_impl(
        &prefix,
        &ImportFiles::new(),
        None,
        None,
        EvaluationState::default(),
        0,
        None,
        dec,
    )?;
    let after_prefix = catallax_numbers::current_operation_work().used;
    let source = &state.bodies[target];
    let geometry = |refusal| BuildError::Geometry {
        step: step.id.clone(),
        refusal,
    };
    let operation = catallax_numbers::OperationSpan::begin("brep", &step.id, "hole", dec);
    let result = (|| -> Result<Option<NumericalHoleBuilt>, Refusal> {
        if !state.stand_ins.is_empty() {
            return Err(Refusal::new(
                "field.stand-in-export",
                "numerical Hole requires retained boundary geometry",
            ));
        }
        crate::hole::validate(diameter.value(), depth, tip, head.as_ref())?;
        let span = crate::hole::HoleTool::new(
            &at[0],
            diameter.value(),
            depth,
            tip,
            head.as_ref(),
            &source.radius,
            dec,
        )?
        .cylinder_span()
        .ok_or_else(|| {
            Refusal::new(
                "hole.numerical-cylinder",
                "compiled Hole has no constant-radius span",
            )
        })?;
        let key = format!("{}/{}", step.id, at[0].id);
        let next = if let Some(previous) = previous {
            previous
                .hole
                .rebuild(&state.model, source.body, &key, &span, options, dec)
        } else {
            crate::hole::build_numerical_hole(&state.model, source.body, &key, &span, options, dec)
        }?;
        let Some(hole) = next else { return Ok(None) };
        let mut origins = BTreeMap::new();
        for (key, source) in hole.source_faces() {
            origins.insert(key.clone(), state.model.face_origins(*source).to_vec());
        }
        let mut wall_origins: Vec<_> = span
            .bands
            .iter()
            .map(|band| catallax_recipe::Origin {
                step: step.id.clone(),
                role: Role::Hole {
                    at: at[0].id.clone(),
                    part: band.part,
                },
            })
            .collect();
        wall_origins.sort();
        wall_origins.dedup();
        for key in hole.tool_faces() {
            origins.insert(key.clone(), wall_origins.clone());
        }
        Ok(Some(NumericalHoleBuilt {
            recipe: recipe.clone(),
            hole,
            origins,
            prefix_work: after_prefix - before,
            assembly_work: catallax_numbers::current_operation_work().used - after_prefix,
        }))
    })();
    operation.finish(dec, result.as_ref().err(), None);
    match result.map_err(geometry)? {
        Some(built) => Ok(HoleAttempt::Numerical(Box::new(built))),
        None => Ok(HoleAttempt::Standard {
            state: Box::new(state),
            start: prefix.steps.len(),
        }),
    }
}

/// Builds a recipe with caller-supplied immutable STEP and mesh sources.
pub fn evaluate_with_files(
    recipe: &Recipe,
    files: &ImportFiles<'_>,
    dec: &mut Decider,
) -> Result<Built, BuildError> {
    evaluate_impl(
        recipe,
        files,
        None,
        None,
        EvaluationState::default(),
        0,
        None,
        dec,
    )
    .map(|state| state.finish(Arc::new(recipe.clone())))
}

/// Evaluates only the suffix after a matching checkpoint.
///
/// Prefix equality includes step order, IDs and every operation parameter.
/// A changed/removed prefix step, format or numeric policy discards the checkpoint
/// and rebuilds from the start. Recipe metadata and result selection may change.
/// Source-byte verification is charged on each call, including prefix sources.
/// Only suffix geometry work is charged anew; cached prefix proofs were charged
/// when that checkpoint was built. Certified builds are not
/// accepted as checkpoints. Failures return no partial checkpoint.
pub fn evaluate_incremental(
    recipe: &Recipe,
    prefix: Option<EvaluationCheckpoint>,
    dec: &mut Decider,
) -> Result<EvaluationCheckpoint, BuildError> {
    evaluate_incremental_with_files(recipe, prefix, &ImportFiles::new(), dec)
}

/// Evaluates a suffix with immutable STEP and mesh sources, including sources in the prefix.
/// Missing or changed bytes are recipe errors even when geometry could be reused.
pub fn evaluate_incremental_with_files(
    recipe: &Recipe,
    prefix: Option<EvaluationCheckpoint>,
    files: &ImportFiles<'_>,
    dec: &mut Decider,
) -> Result<EvaluationCheckpoint, BuildError> {
    let (state, start) = resume_state(recipe, prefix, dec);
    let state = evaluate_impl(recipe, files, None, None, state, start, None, dec)?;
    Ok(checkpoint(state, recipe, dec))
}

fn resume_state(
    recipe: &Recipe,
    prefix: Option<EvaluationCheckpoint>,
    dec: &Decider,
) -> (EvaluationState, usize) {
    match prefix {
        Some(prefix)
            if prefix.policy == *dec.policy()
                && prefix.built.recipe().format == recipe.format
                && recipe.steps.starts_with(&prefix.built.recipe().steps) =>
        {
            let start = prefix.built.recipe().steps.len();
            (prefix.into_state(), start)
        }
        _ => (EvaluationState::default(), 0),
    }
}

fn checkpoint(mut state: EvaluationState, recipe: &Recipe, dec: &Decider) -> EvaluationCheckpoint {
    let profiles = std::mem::take(&mut state.profiles);
    let paths = std::mem::take(&mut state.paths);
    let bodies = state.bodies.clone();
    let meshes = state.meshes.clone();
    EvaluationCheckpoint {
        built: Arc::new(state.finish(Arc::new(recipe.clone()))),
        profiles,
        paths,
        bodies,
        meshes,
        policy: *dec.policy(),
    }
}

/// Opts into the recorded 1 nm modelling policy and the first complete bounded
/// offset family. Legacy recipe evaluation remains exact-only and unchanged.
pub fn evaluate_certified(
    recipe: &Recipe,
    dec: &mut Decider,
) -> Result<crate::CertifiedBuilt, BuildError> {
    let mut certification = crate::certified::CertificationContext::new(recipe);
    let state = evaluate_impl(
        recipe,
        &ImportFiles::new(),
        Some(&mut certification),
        None,
        EvaluationState::default(),
        0,
        None,
        dec,
    )?;
    let built = state.finish(
        certification
            .snapshot()
            .unwrap_or_else(|| Arc::new(recipe.clone())),
    );
    certification
        .finish(built, dec)
        .map_err(|refusal| BuildError::Geometry {
            step: recipe.result.clone(),
            refusal,
        })
}

#[allow(clippy::too_many_arguments)]
fn evaluate_impl(
    recipe: &Recipe,
    files: &ImportFiles<'_>,
    certification: Option<&mut crate::certified::CertificationContext<'_>>,
    numerical_sweep: Option<&mut crate::general_sweep::Context>,
    state: EvaluationState,
    start: usize,
    numerical: Option<&mut numerical_fillet::Context>,
    dec: &mut Decider,
) -> Result<EvaluationState, BuildError> {
    recipe.check().map_err(BuildError::Recipe)?;
    catallax_numbers::reuse_operation_enclosures();

    let mut verified_sources = BTreeSet::new();
    for step in &recipe.steps {
        let (file, kind) = match &step.op {
            Operation::ImportStep { file, .. } => (file, "STEP"),
            Operation::ImportIges { file, .. } => (file, "IGES"),
            Operation::ImportMesh { file, .. } => (file, "mesh"),
            _ => continue,
        };
        if !verified_sources.insert(file.sha256.as_str()) {
            continue;
        }
        let source = files.get(&file.sha256).ok_or_else(|| {
            BuildError::Recipe(RecipeError::Parse {
                message: format!(
                    "step {}: missing {kind} source {} with sha256 {}",
                    step.id, file.name, file.sha256
                ),
            })
        })?;
        dec.charge_work(source.len() as u64)
            .map_err(|refusal| BuildError::Geometry {
                step: step.id.clone(),
                refusal,
            })?;
        let actual = catallax_import::sha256_hex(source);
        if actual != file.sha256 {
            return Err(BuildError::Recipe(RecipeError::Parse {
                message: format!(
                    "step {}: {kind} source {} sha256 mismatch: expected {}, actual {actual}",
                    step.id, file.name, file.sha256
                ),
            }));
        }
    }
    let EvaluationState {
        mut model,
        mut profiles,
        mut paths,
        mut bodies,
        mut meshes,
        mut fillets,
        mut consumed,
        mut chamfers,
        mut stand_ins,
        mut import_records,
        mut import_solids,
    } = state;
    // Each step is a task. It waits for the steps it names as inputs and, when
    // it reads or extends the fillet, draft and chamfer registries carried from
    // step to step, for the last step that changed them. Everything else a step
    // reads comes through the model, which tracks its own reads.
    #[cfg(feature = "work-span")]
    let tasks = catallax_numbers::work_span::steps("recipe steps");
    // Step outputs and registries of this evaluation, by step position.
    #[cfg(feature = "work-span")]
    let (outputs, position) = (
        catallax_numbers::work_span::Instance::new(),
        |id: &catallax_recipe::StepId| {
            recipe
                .steps
                .iter()
                .position(|step| &step.id == id)
                .map_or(u64::MAX, |index| index as u64)
        },
    );
    let mut modes = step::Modes {
        certification,
        numerical_sweep,
        numerical,
    };
    for step in &recipe.steps[start..] {
        let declared = step::Effects::of(&step.op);
        #[cfg(feature = "work-span")]
        let _task = tasks.task();
        #[cfg(feature = "work-span")]
        let _output = {
            use catallax_numbers::work_span;
            for (input, _) in step.op.inputs() {
                work_span::read_at("recipe step output", outputs.id(), position(input));
            }
            let registry = declared.feature_registries();
            if registry {
                work_span::read_at("feature registries", outputs.id(), u64::MAX);
            }
            (
                work_span::write_on_drop_at("recipe step output", outputs.id(), position(&step.id)),
                registry.then(|| {
                    work_span::write_on_drop_at("feature registries", outputs.id(), u64::MAX)
                }),
            )
        };
        let inputs = step::Inputs {
            profiles: &profiles,
            paths: &paths,
            bodies: &bodies,
            meshes: &meshes,
            records: step::Published::new(&fillets, &chamfers),
        };
        let step::Outcome {
            value,
            changes,
            effects,
        } = step::run(recipe, files, &mut modes, step, &inputs, &mut model, dec).map_err(
            |refusal| BuildError::Geometry {
                step: step.id.clone(),
                refusal,
            },
        )?;
        debug_assert!(
            declared.covers(effects),
            "step {} touched {effects:?}, beyond its declared {declared:?}",
            step.id
        );
        // Publish in step order: the records first, then the step's value.
        fillets.extend(changes.fillets);
        consumed.extend(changes.consumed);
        chamfers.extend(changes.chamfers);
        stand_ins.extend(changes.stand_ins);
        let id = step.id.clone();
        match value {
            step::Value::Nothing | step::Value::Numerical => {}
            step::Value::Profile(profile) => {
                profiles.insert(id, profile);
            }
            step::Value::Path(path) => {
                paths.insert(id, path);
            }
            step::Value::Body(body) => {
                bodies.insert(id, body);
            }
            step::Value::Mesh(body, mesh) => {
                meshes.insert(id.clone(), mesh);
                bodies.insert(id, body);
            }
            step::Value::Import(body, record, solids) => {
                import_solids.insert(id.clone(), solids);
                import_records.insert(id.clone(), record);
                bodies.insert(id, body);
            }
        }
    }
    Ok(EvaluationState {
        model,
        profiles,
        paths,
        bodies,
        meshes,
        fillets,
        consumed,
        chamfers,
        stand_ins,
        import_records,
        import_solids,
    })
}

pub(crate) fn import_refusal(error: catallax_import::ImportRefusal) -> Refusal {
    use catallax_import::ImportRefusal;
    let code = match &error {
        ImportRefusal::PlanarRepairNotVerified { .. } => "import.planar-repair-not-verified",
        ImportRefusal::Malformed { .. } => "import.malformed",
        ImportRefusal::Unsupported { .. } => "import.unsupported",
        ImportRefusal::UnknownUnit { .. } => "import.unknown-unit",
        ImportRefusal::HashMismatch { .. } => "import.hash-mismatch",
        ImportRefusal::NotSolid { .. } => "import.not-solid",
        ImportRefusal::IncidenceGap { .. } => "import.incidence-gap",
        ImportRefusal::InvalidDecision { .. } => "import.invalid-decision",
        ImportRefusal::ValidationFailed { .. } => "import.validation-failed",
        ImportRefusal::SelfIntersecting { .. } => "import.self-intersecting",
        ImportRefusal::SelfIntersectionNotVerified { .. } => {
            "import.self-intersection-not-verified"
        }
        ImportRefusal::ResourceLimit => "import.resource-limit",
    };
    // An undecided pair keeps the decision it stopped on, so a work or
    // precision limit is not reported as unsupported geometry.
    let undecided = match &error {
        ImportRefusal::SelfIntersectionNotVerified {
            refusal: Some(refusal),
            ..
        } => refusal.undecided().cloned(),
        _ => None,
    };
    Refusal::new(code, format!("{error:?}")).keeping_any(undecided)
}

pub(crate) fn numerical_sweep_impl(
    recipe: &Recipe,
    context: &mut crate::general_sweep::Context,
    dec: &mut Decider,
) -> Result<(Model, BodyId), BuildError> {
    let state = evaluate_impl(
        recipe,
        &ImportFiles::new(),
        None,
        Some(context),
        EvaluationState::default(),
        0,
        None,
        dec,
    )?;
    let body = state.bodies[&recipe.result].body;
    Ok((state.model, body))
}

/// Evaluate the opt-in numerical source and its edits atomically. Every stage shares
/// the active request allowance, or the unchanged 20M default. No exact checkpoint
/// or bounded-answer conversion is made.
pub fn evaluate_numerical(
    recipe: &crate::numerical::NumericalRecipe,
    dec: &mut Decider,
) -> Result<crate::numerical::NumericalBuilt, Refusal> {
    numerical_request(|| crate::numerical::evaluate(recipe, dec))
}

/// Replay a saved numerical definition. Parsing and all build/edit stages share
/// one request allowance. Saved geometry or evidence is never trusted as checked.
pub fn replay_numerical(
    input: &str,
    dec: &mut Decider,
) -> Result<crate::numerical::NumericalBuilt, Refusal> {
    numerical_request(|| {
        let recipe = crate::numerical::NumericalRecipe::parse(input, dec)?;
        crate::numerical::evaluate(&recipe, dec)
    })
}

fn numerical_request<T>(run: impl FnOnce() -> Result<T, Refusal>) -> Result<T, Refusal> {
    if catallax_numbers::operation_budget_active() {
        return run();
    }
    let (result, work) =
        catallax_numbers::with_operation_budget(catallax_topology::numerical::DEFAULT_WORK, run);
    if work.exhausted {
        Err(Refusal::new(
            "numbers.operation-budget",
            "numerical edits exhausted their shared work allowance",
        ))
    } else {
        result
    }
}

#[cfg(test)]
mod import_refusal_tests {
    use super::*;
    use catallax_numbers::{NotExact, Rational};

    // An undecided self-intersection pair keeps the decision it stopped on in
    // the kernel's refusal; an uncovered pair keeps none.
    #[test]
    fn an_undecided_import_pair_keeps_its_decision() {
        let decision = NotExact {
            bound: Rational::from_ratio(1, 1 << 20).unwrap(),
            code: "numbers.precision-limit",
            detail: "order along an arc".into(),
        };
        let pair = |refusal: Refusal| {
            import_refusal(
                catallax_import::ImportRefusal::SelfIntersectionNotVerified {
                    solid: "#1".into(),
                    faces: ["#2".into(), "#3".into()],
                    reason: refusal.to_string(),
                    refusal: Some(Box::new(refusal)),
                },
            )
        };
        let refusal = pair(Refusal::from(decision.clone()));
        assert_eq!(refusal.code, "import.self-intersection-not-verified");
        assert_eq!(refusal.undecided(), Some(&decision));
        let uncovered = pair(Refusal::new("import.pair-undecided", "an unbounded face"));
        assert!(uncovered.undecided().is_none());
    }
}
