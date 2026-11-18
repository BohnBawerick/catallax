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

fn source_fillet_indices(
    model: &Model,
    body: BodyId,
    fillets: &[crate::fillet::FilletDefinition],
) -> Vec<usize> {
    if fillets.is_empty() {
        return Vec::new();
    }
    let faces = model.faces_of_body(body);
    fillets
        .iter()
        .enumerate()
        .filter_map(|(index, definition)| {
            faces
                .iter()
                .any(|face| model.face_origins(*face).contains(&definition.origin))
                .then_some(index)
        })
        .collect()
}

fn copy_fillets(
    fillets: &mut Vec<crate::fillet::FilletDefinition>,
    indices: &[usize],
    map: &Isometry,
    step: &StepId,
    instance: u32,
) -> Result<(), Refusal> {
    for &index in indices {
        let mut definition = fillets[index].clone();
        definition.transform(map)?;
        definition.origin = crate::copy_origin(step, instance, &definition.origin);
        fillets.push(definition);
    }
    Ok(())
}

fn vector(v: &catallax_recipe::Vector3) -> RVector3 {
    RVector3::new(
        v.0[0].value().clone(),
        v.0[1].value().clone(),
        v.0[2].value().clone(),
    )
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

fn transformed_bounds(
    bounds: &BoundingBox,
    map: &Isometry,
    dec: &mut Decider,
) -> Result<BoundingBox, Refusal> {
    let local = std::array::from_fn::<_, 3, _>(|i| {
        Ball::from_rational(&bounds.min[i], 128).hull(&Ball::from_rational(&bounds.max[i], 128))
    });
    let mut enclose = |value: &Real| {
        if let Some(value) = value.rational_ref() {
            Ok(Ball::from_rational(value, 128))
        } else {
            dec.enclose_budgeted_at(value, 128)
        }
    };
    let translation = [&map.translation.x, &map.translation.y, &map.translation.z];
    let mut min = std::array::from_fn(|_| Rational::zero());
    let mut max = min.clone();
    for i in 0..3 {
        let mut coordinate = enclose(translation[i])?;
        for (j, local) in local.iter().enumerate() {
            coordinate = &coordinate + &(enclose(&map.matrix()[i][j])? * local.clone());
        }
        min[i] = coordinate.lo();
        max[i] = coordinate.hi();
    }
    Ok(BoundingBox { min, max })
}

fn separated(a: &BoundingBox, b: &BoundingBox) -> bool {
    (0..3).any(|i| a.max[i] < b.min[i] || b.max[i] < a.min[i])
}

fn require_material(model: &Model, bodies: &[BodyId]) -> Result<(), Refusal> {
    if bodies
        .iter()
        .any(|body| model.body(*body).kind == BodyKind::Surface)
    {
        return Err(Refusal::new(
            "boolean.not-solid",
            "a surface body has no solid material",
        ));
    }
    Ok(())
}

// Strictly separated enclosing boxes prove that regularised union is just the
// disjoint collection of the existing immutable lumps. No support intersection
// or replacement geometry is needed, and contact never takes this path.
fn unite_separated(model: &mut Model, a: BodyId, b: BodyId) -> BodyId {
    let mut result = model.body(a).clone();
    result.lumps.extend_from_slice(&model.body(b).lumps);
    result.marks.extend(&model.body(b).marks);
    result.kind = if result.lumps.is_empty() {
        BodyKind::Empty
    } else {
        BodyKind::Solid
    };
    model.add_body(result)
}

fn through(
    model: &mut Model,
    profile: &Profile,
    step: &StepId,
    target: &BodyValue,
    direction: Option<&catallax_recipe::Vector3>,
    side: Option<ThroughSide>,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let normal = crate::profile::rational_vector(
        &profile.pose,
        &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
    );
    let direction = direction.map(vector).unwrap_or_else(|| normal.clone());
    if direction.is_zero() {
        return Err(Refusal::new(
            "extrude.zero-direction",
            "extrusion direction must not be zero",
        ));
    }
    let dot = direction.dot(&normal);
    if dot.is_zero() {
        return Err(Refusal::new(
            "extrude.zero-thickness",
            "extrusion direction lies in the sketch plane",
        ));
    }
    if model.body(target.body).kind == BodyKind::Empty || profile.is_empty() {
        return Ok(target.body);
    }
    let origin = crate::profile::rational_point(&profile.pose, &RPoint3::origin());
    let normal_l1 = normal.x.abs() + normal.y.abs() + normal.z.abs();
    // For every target point p, |(p-o).n| <= (radius + |o|_1)|n|_1.
    // The strict extra unit leaves both construction caps outside the target.
    // Rescaling direction changes t inversely, not the resulting through-all set.
    let extent =
        ((&target.radius + &norm_upper(&origin)) * normal_l1) / dot.abs() + Rational::one();
    let (low, high, caps) = match side.unwrap_or(ThroughSide::Both) {
        ThroughSide::Both => (-&extent, extent.clone(), [None, None]),
        ThroughSide::Forward => (Rational::zero(), extent.clone(), [Some(Role::Start), None]),
        ThroughSide::Backward => (-&extent, Rational::zero(), [None, Some(Role::Start)]),
    };
    let mut shifted = profile.clone();
    shifted.pose.translation = shifted
        .pose
        .translation
        .plus(&direction.scaled(&low).to_real());
    let tool = crate::extrude::extrude_vector(
        model,
        &shifted,
        step,
        &direction.scaled(&(high - low)),
        caps,
        None,
        dec,
    )?;
    Ok(boolean(model, target.body, tool, BoolOp::Subtract, dec)
        .map_err(|e| e.refusal)?
        .body)
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
    EvaluationCheckpoint {
        built: Arc::new(state.finish(Arc::new(recipe.clone()))),
        profiles,
        paths,
        bodies,
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
    mut certification: Option<&mut crate::certified::CertificationContext<'_>>,
    mut numerical_sweep: Option<&mut crate::general_sweep::Context>,
    state: EvaluationState,
    start: usize,
    mut numerical: Option<&mut numerical_fillet::Context>,
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
    for step in &recipe.steps[start..] {
        #[cfg(feature = "work-span")]
        let _task = tasks.task();
        #[cfg(feature = "work-span")]
        let _output = {
            use catallax_numbers::work_span;
            for (input, _) in step.op.inputs() {
                work_span::read_at("recipe step output", outputs.id(), position(input));
            }
            let registry = matches!(
                step.op,
                Operation::Draft { .. }
                    | Operation::Fillet { .. }
                    | Operation::FilletVariable { .. }
                    | Operation::Chamfer { .. }
                    | Operation::Pattern { .. }
                    | Operation::Mirror { .. }
            );
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
        let span = catallax_numbers::OperationSpan::begin("brep", &step.id, step.op.name(), dec);
        let mut build = |model: &mut Model,
                         dec: &mut Decider|
         -> Result<Option<BodyValue>, Refusal> {
            if let Some(context) = numerical.as_deref() {
                context.exact_inputs(step)?;
            }
            if let Some((body, radius)) = crate::stock::build(model, &step.op, &step.id, dec)? {
                return Ok(Some(BodyValue { body, radius }));
            }
            Ok(Some(match &step.op {
                Operation::Sketch { .. } => return Ok(None),
                Operation::ImportStep {
                    source_solid,
                    file,
                    pose,
                    decisions,
                }
                | Operation::ImportIges {
                    source_solid,
                    file,
                    pose,
                    decisions,
                } => {
                    let format = if matches!(step.op, Operation::ImportIges { .. }) {
                        catallax_import::SourceFormat::Iges
                    } else {
                        catallax_import::SourceFormat::Step
                    };
                    let decisions = crate::import_edit::import_decisions(&file.sha256, decisions)?;
                    let mut imported = catallax_import::apply_source(
                        model,
                        format,
                        files[&file.sha256],
                        &file.name,
                        &Default::default(),
                        &decisions,
                        Some(&step.id),
                    )
                    .map_err(import_refusal)?;
                    if let Some(entity) = source_solid {
                        imported.bodies.retain(|(key, _)| key == entity);
                        imported.left_open.retain(|key| key == entity);
                        if imported.bodies.is_empty() && imported.left_open.is_empty() {
                            return Err(Refusal::new(
                                "import.source-solid",
                                format!("{entity} is not a source solid in {}", file.sha256),
                            ));
                        }
                    }
                    if !imported.left_open.is_empty() {
                        return Err(Refusal::new(
                            "import.left-open",
                            format!(
                                "{} not imported until the gap is resolved. The original file was never changed",
                                imported.left_open.join(", ")
                            ),
                        ));
                    }
                    let raw = crate::import_edit::collect_solids(model, &imported.bodies, dec)?;
                    let body = if pose == &catallax_recipe::Pose::identity() {
                        raw
                    } else {
                        let placement = crate::profile::placement(pose, dec)?;
                        catallax_topology::transformed_body(
                            model,
                            raw,
                            &Isometry::from(&placement),
                            Clone::clone,
                        )
                    };
                    let bounds =
                        catallax_measure::bounding_box(model, body, &Default::default(), dec)?
                            .value
                            .ok_or_else(|| {
                                Refusal::new("import.not-solid", "the imported body has no bounds")
                            })?;
                    let radius = (0..3).fold(Rational::zero(), |sum, i| {
                        sum + bounds.min[i].abs().max(bounds.max[i].abs())
                    });
                    let mut first_lump = 0;
                    let membership = imported
                        .bodies
                        .iter()
                        .map(|(entity, raw)| {
                            let count = model.body(*raw).lumps.len();
                            let lumps =
                                model.body(body).lumps[first_lump..first_lump + count].to_vec();
                            first_lump += count;
                            crate::import_edit::ImportedSolid {
                                entity: entity.clone(),
                                lumps,
                            }
                        })
                        .collect();
                    import_solids.insert(step.id.clone(), membership);
                    import_records.insert(step.id.clone(), imported.record);
                    BodyValue { body, radius }
                }
                Operation::ImportMesh { file, pose, scale } => {
                    let mesh = catallax_import::mesh::ClosedMesh::read(
                        files[&file.sha256],
                        scale.value(),
                        dec,
                    )?;
                    let placement = crate::profile::placement(pose, dec)?;
                    let (body, radius) = crate::mesh_import::build(
                        model,
                        &mesh,
                        &Isometry::from(&placement),
                        &step.id,
                        dec,
                    )?;
                    BodyValue { body, radius }
                }
                Operation::Extrude {
                    profile,
                    distance,
                    direction,
                } => {
                    let profile = &profiles[profile];
                    let direction = direction.as_ref().map(vector);
                    let body = crate::extrude::extrude_with_decider(
                        model,
                        profile,
                        &step.id,
                        distance.value(),
                        direction.as_ref(),
                        dec,
                    )?;
                    BodyValue {
                        body,
                        radius: profile_radius(profile) + distance.value().abs(),
                    }
                }
                Operation::ExtrudeTapered {
                    profile,
                    distance,
                    angle,
                } => {
                    let (body, radius) = crate::tapered_extrude::with_radius(
                        model,
                        &profiles[profile],
                        &step.id,
                        distance.value(),
                        angle.value(),
                        dec,
                    )?;
                    BodyValue { body, radius }
                }
                Operation::ExtrudeOnFace {
                    profile,
                    body,
                    face,
                    normal_direction,
                    x_direction,
                    offset,
                    distance,
                    placement,
                } => {
                    let profile = &profiles[profile];
                    crate::check_face_profile(&profile.pose, dec)?;
                    let support = &bodies[body];
                    let planes = model
                        .faces_of_body(support.body)
                        .into_iter()
                        .filter(|id| {
                            model
                                .face_origins(*id)
                                .iter()
                                .any(|origin| face.matches(origin))
                        })
                        .map(|id| {
                            let native = model.face(id);
                            let catallax_shapes::Surface::Plane(plane) =
                                model.geom().surface(native.surface)
                            else {
                                return Err(Refusal::new(
                                    "extrude.nonplanar-face",
                                    "the sketch support is not planar",
                                ));
                            };
                            Ok(plane.clone())
                        })
                        .collect::<Result<Vec<_>, Refusal>>()?;
                    let map = crate::face_placement(
                        &planes,
                        &vector(normal_direction),
                        &vector(x_direction),
                        offset.value(),
                        placement,
                        dec,
                    )?;
                    let map = Isometry::from(&map);
                    let local = crate::extrude::extrude_with_decider(
                        model,
                        profile,
                        &step.id,
                        distance.value(),
                        None,
                        dec,
                    )?;
                    let body =
                        catallax_topology::transformed_body(model, local, &map, Clone::clone);
                    BodyValue {
                        body,
                        radius: support.radius.clone()
                            + profile_radius(profile)
                            + offset.value().abs()
                            + distance.value().abs()
                            + placement
                                .translation
                                .0
                                .iter()
                                .fold(Rational::zero(), |sum, n| sum + n.value().abs()),
                    }
                }
                Operation::Loft { sections, shape } => {
                    let (body, radius) = crate::loft::build(
                        model, recipe, sections, *shape, &step.id, &profiles, dec,
                    )?;
                    BodyValue { body, radius }
                }
                Operation::Sweep {
                    profile,
                    path,
                    orientation,
                    guide,
                    twist,
                } => {
                    let helical = recipe.steps.iter().find(|s| &s.id == path)
                        .is_some_and(|s| matches!(&s.op,Operation::Path {segments,..}
                            if segments.iter().any(|s|matches!(s,catallax_recipe::PathSegment::Helix {..}))));
                    if helical
                        || guide.is_some()
                        || twist.as_ref().is_some_and(|n| !n.value().is_zero())
                    {
                        let context=numerical_sweep.as_deref_mut().ok_or_else(||Refusal::new(
                            "sweep.requires-numerical-evaluation",
                            "helical sweeps require numerical admission under a stated tolerance"))?;
                        let (body, radius) = context.sweep(model, recipe, step, dec)?;
                        BodyValue { body, radius }
                    } else {
                        let profile = &profiles[profile];
                        let path = &paths[path];
                        let body = path.build(model, profile, &step.id, *orientation, dec)?;
                        BodyValue {
                            body,
                            radius: profile_radius(profile) + path.extent(),
                        }
                    }
                }
                Operation::Revolve {
                    profile,
                    axis,
                    angle,
                } => {
                    let profile = &profiles[profile];
                    let body = revolve(model, profile, &step.id, axis, angle.value(), dec)?;
                    let origin = RPoint3::new(
                        axis.point.0[0].value().clone(),
                        axis.point.0[1].value().clone(),
                        axis.point.0[2].value().clone(),
                    );
                    BodyValue {
                        body,
                        radius: profile_radius(profile)
                            + Rational::from_i64(2) * norm_upper(&origin),
                    }
                }
                Operation::ExtrudeThrough {
                    target,
                    profile,
                    direction,
                    side,
                } => {
                    let target = &bodies[target];
                    let body = through(
                        model,
                        &profiles[profile],
                        &step.id,
                        target,
                        direction.as_ref(),
                        *side,
                        dec,
                    )?;
                    BodyValue {
                        body,
                        radius: target.radius.clone(),
                    }
                }
                Operation::ExtrudeTo {
                    target,
                    profile,
                    direction,
                    limit,
                } => {
                    profiles[profile].require_circular("extrude-to.ellipse-profile")?;
                    if !profiles[profile].is_empty()
                        && matches!(
                            recipe.step(target).map(|s| &s.op),
                            Some(Operation::Torus { .. })
                        )
                    {
                        let cap = crate::end_condition::toroidal::AnnularCap::compile(
                            recipe,
                            target,
                            profile,
                            direction.as_ref(),
                            Some(limit),
                            dec,
                        )?;
                        BodyValue {
                            body: cap.build(
                                model,
                                bodies[target].body,
                                &profiles[profile],
                                &step.id,
                                dec,
                            )?,
                            radius: cap.radius().clone(),
                        }
                    } else if !profiles[profile].is_empty()
                        && matches!(
                            recipe.step(target).map(|s| &s.op),
                            Some(
                                Operation::Sphere { .. }
                                    | Operation::Cylinder { .. }
                                    | Operation::Cone { .. }
                            )
                        )
                    {
                        let cap = crate::end_condition::curved::AuthoredCap::compile(
                            recipe,
                            target,
                            profile,
                            direction.as_ref(),
                            Some(limit),
                            dec,
                        )?;
                        BodyValue {
                            body: cap.build(model, &step.id, dec)?,
                            radius: cap.radius().clone(),
                        }
                    } else {
                        let cap = crate::end_condition::Capped::compile(
                            recipe,
                            target,
                            profile,
                            direction.as_ref(),
                            Some(limit),
                            dec,
                        )?;
                        BodyValue {
                            body: cap.build(model, &step.id, dec)?,
                            radius: cap.radius,
                        }
                    }
                }
                Operation::ExtrudeNext {
                    target,
                    profile,
                    direction,
                } => {
                    profiles[profile].require_circular("extrude-next.ellipse-profile")?;
                    if !profiles[profile].is_empty()
                        && matches!(
                            recipe.step(target).map(|s| &s.op),
                            Some(Operation::Torus { .. })
                        )
                    {
                        let cap = crate::end_condition::toroidal::AnnularCap::compile(
                            recipe,
                            target,
                            profile,
                            direction.as_ref(),
                            None,
                            dec,
                        )?;
                        BodyValue {
                            body: cap.build(
                                model,
                                bodies[target].body,
                                &profiles[profile],
                                &step.id,
                                dec,
                            )?,
                            radius: cap.radius().clone(),
                        }
                    } else if !profiles[profile].is_empty()
                        && matches!(
                            recipe.step(target).map(|s| &s.op),
                            Some(
                                Operation::Sphere { .. }
                                    | Operation::Cylinder { .. }
                                    | Operation::Cone { .. }
                            )
                        )
                    {
                        let cap = crate::end_condition::curved::AuthoredCap::compile(
                            recipe,
                            target,
                            profile,
                            direction.as_ref(),
                            None,
                            dec,
                        )?;
                        BodyValue {
                            body: cap.build(model, &step.id, dec)?,
                            radius: cap.radius().clone(),
                        }
                    } else {
                        let cap = crate::end_condition::Capped::compile(
                            recipe,
                            target,
                            profile,
                            direction.as_ref(),
                            None,
                            dec,
                        )?;
                        BodyValue {
                            body: cap.build(model, &step.id, dec)?,
                            radius: cap.radius,
                        }
                    }
                }
                Operation::Hole {
                    target,
                    at,
                    diameter,
                    depth,
                    tip,
                    head,
                } => {
                    crate::hole::validate(diameter.value(), depth, tip, head.as_ref())?;
                    let target = &bodies[target];
                    let mut body = target.body;
                    let tools = at
                        .iter()
                        .map(|placement| {
                            crate::hole::HoleTool::new(
                                placement,
                                diameter.value(),
                                depth,
                                tip,
                                head.as_ref(),
                                &target.radius,
                                dec,
                            )
                        })
                        .collect::<Result<Vec<_>, _>>()?;
                    if model.body(body).kind != BodyKind::Empty {
                        let mut disjoint = true;
                        for i in 0..tools.len() {
                            for other in &tools[..i] {
                                if !tools[i].disjoint(other, dec)? {
                                    disjoint = false;
                                    break;
                                }
                            }
                            if !disjoint {
                                break;
                            }
                        }
                        if disjoint && !tools.is_empty() {
                            let mut lumps = Vec::new();
                            for (tool, placement) in tools.iter().zip(at) {
                                let tool = tool.build(model, &step.id, &placement.id);
                                lumps.extend_from_slice(&model.body(tool).lumps);
                            }
                            let tool = model.add_body(catallax_topology::Body {
                                kind: BodyKind::Solid,
                                lumps,
                                wires: Vec::new(),
                                marks: Default::default(),
                            });
                            body = boolean(model, body, tool, BoolOp::Subtract, dec)
                                .map_err(|e| e.refusal)?
                                .body;
                        } else {
                            for (tool, placement) in tools.iter().zip(at) {
                                let tool = tool.build(model, &step.id, &placement.id);
                                body = boolean(model, body, tool, BoolOp::Subtract, dec)
                                    .map_err(|e| e.refusal)?
                                    .body;
                            }
                        }
                    }
                    BodyValue {
                        body,
                        radius: target.radius.clone(),
                    }
                }
                Operation::Boolean { kind, target, tool } => {
                    let a = &bodies[target];
                    let b = &bodies[tool];
                    require_material(model, &[a.body, b.body])?;
                    if let Some(context) = certification.as_deref_mut()
                        && (context.has_certificate(target) || context.has_certificate(tool))
                    {
                        let (body, radius) = context.boolean(
                            model,
                            &step.id,
                            kind,
                            [(target, a.body), (tool, b.body)],
                            dec,
                        )?;
                        return Ok(Some(BodyValue { body, radius }));
                    }
                    let (op, radius) = match kind {
                        BooleanKind::Union => {
                            (BoolOp::Unite, a.radius.clone().max(b.radius.clone()))
                        }
                        BooleanKind::Subtract => (BoolOp::Subtract, a.radius.clone()),
                        BooleanKind::Intersect => {
                            (BoolOp::Intersect, a.radius.clone().min(b.radius.clone()))
                        }
                    };
                    // Recipes allow repeated operands. Apply exact set identities
                    // before the two-distinct-operands intersection graph.
                    let body = if a.body == b.body {
                        if *kind == BooleanKind::Subtract {
                            // An empty result has no fitted face.
                            model.add_body(catallax_topology::Body {
                                kind: BodyKind::Empty,
                                lumps: Vec::new(),
                                wires: Vec::new(),
                                marks: catallax_topology::without_fit_mark(
                                    &model.body(a.body).marks,
                                ),
                            })
                        } else {
                            a.body
                        }
                    } else {
                        boolean(model, a.body, b.body, op, dec)
                            .map_err(|e| e.refusal)?
                            .body
                    };
                    BodyValue { body, radius }
                }
                Operation::Draft {
                    body,
                    faces,
                    pull,
                    neutral,
                    angle,
                    keep_regions,
                } => {
                    let definition = crate::draft::DraftDefinition::new(
                        vector(pull).to_real(),
                        crate::draft::neutral_plane(neutral),
                        angle.value(),
                        dec,
                    )?;
                    let target = &bodies[body];
                    let (body, history, radius) = crate::draft_build::build(
                        model,
                        target.body,
                        faces,
                        &definition,
                        keep_regions.as_deref(),
                        &step.id,
                        dec,
                    )?;
                    consumed.extend(history);
                    BodyValue { body, radius }
                }
                Operation::Fillet {
                    body,
                    edges,
                    radius,
                    overflow,
                    corner,
                } => {
                    let target = &bodies[body];
                    let mut exact = |model: &mut Model, dec: &mut Decider| {
                        if let Some(catallax_recipe::FilletOverflow::KeepNeighbours(limits)) =
                            overflow
                        {
                            crate::fillet_build::fillet_setback::build(
                                model,
                                target.body,
                                edges,
                                &radius.value().clone().into(),
                                limits,
                                &step.id,
                                dec,
                            )
                        } else {
                            crate::fillet_build::build(
                                model,
                                target.body,
                                edges,
                                &radius.value().clone().into(),
                                recipe.format != catallax_recipe::FORMAT,
                                corner.unwrap_or(catallax_recipe::FilletCorner::Cliff),
                                &step.id,
                                &fillets,
                                &mut consumed,
                                dec,
                            )
                        }
                    };
                    let (body, definitions) = match numerical.as_deref_mut() {
                        None => exact(model, dec)?,
                        // The exact builders are an optional proof of the
                        // same fillet: a body they finish is the exact result,
                        // charged as it always was. A refusal of theirs is no
                        // verdict on the fillet and returns what it charged,
                        // so the rolling-ball construction keeps the whole
                        // allowance it found.
                        Some(context) => {
                            match catallax_numbers::refund_on_refusal(|| exact(model, dec)) {
                                Ok(result) => result,
                                Err(refusal) => {
                                    context.fillet(
                                        model,
                                        target,
                                        edges,
                                        radius.value(),
                                        overflow.as_ref(),
                                        *corner,
                                        &fillets,
                                        &step.id,
                                        recipe,
                                        refusal,
                                        dec,
                                    )?;
                                    return Ok(None);
                                }
                            }
                        }
                    };
                    fillets.extend(definitions);
                    BodyValue {
                        body,
                        radius: &target.radius + &(radius.value() * &Rational::from_i64(2)),
                    }
                }
                Operation::FilletVariable { body, edges } => {
                    let target = &bodies[body];
                    let mut maximum = Rational::zero();
                    for law in edges {
                        let catallax_recipe::RadiusLaw::Linear { start, end } = &law.radius else {
                            return Err(Refusal::new(
                                "fillet.unsupported-variable-contact",
                                "the radius law has no construction definition",
                            ));
                        };
                        maximum = maximum.max(start.value().clone()).max(end.value().clone());
                    }
                    let (body, definitions) = crate::fillet_build::build_variable(
                        model,
                        target.body,
                        edges,
                        &step.id,
                        &fillets,
                        &mut consumed,
                        dec,
                    )?;
                    fillets.extend(definitions);
                    BodyValue {
                        body,
                        radius: &target.radius + &(&maximum * &Rational::from_i64(2)),
                    }
                }
                Operation::Chamfer {
                    body,
                    edges,
                    distance,
                    other_distance,
                } => {
                    let target = &bodies[body];
                    let other = other_distance.as_ref().unwrap_or(distance);
                    let distances = [
                        distance.value().clone().into(),
                        other.value().clone().into(),
                    ];
                    let (body, definitions) = crate::chamfer_build::build(
                        model,
                        target.body,
                        edges,
                        &distances,
                        other_distance.is_some(),
                        &step.id,
                        dec,
                    )?;
                    chamfers.extend(definitions);
                    BodyValue {
                        body,
                        radius: &(&target.radius + distance.value()) + other.value(),
                    }
                }
                Operation::Pattern { body, layout } => {
                    // Validate even count-one and empty patterns before set identities.
                    crate::pattern_placement(layout, 0, dec)?;
                    let source = &bodies[body];
                    if layout.count() > 1 {
                        require_material(model, &[source.body])?;
                    }
                    let mut body = source.body;
                    let mut radius = source.radius.clone();
                    let source_chamfers: Vec<_> = chamfers
                        .iter()
                        .filter(|d| {
                            model
                                .faces_of_body(source.body)
                                .iter()
                                .any(|f| model.face_origins(*f).contains(&d.origin))
                        })
                        .cloned()
                        .collect();
                    let source_fillets = if layout.count() > 1 {
                        source_fillet_indices(model, source.body, &fillets)
                    } else {
                        Vec::new()
                    };
                    let source_bounds = if layout.count() > 1 {
                        catallax_measure::bounding_box(
                            model,
                            source.body,
                            &Default::default(),
                            dec,
                        )?
                        .value
                    } else {
                        None
                    };
                    let mut placed_bounds = source_bounds.iter().cloned().collect::<Vec<_>>();
                    for instance in 1..layout.count() {
                        let map = crate::pattern_placement(layout, instance, dec)?;
                        let shift = [&map.translation.x, &map.translation.y, &map.translation.z]
                            .into_iter()
                            .try_fold(Rational::zero(), |sum, v| {
                                let b = dec.enclose_at(v, 128)?;
                                Ok::<_, Refusal>(sum + b.lo().abs().max(b.hi().abs()))
                            })?;
                        radius = radius.max(&source.radius + &shift);
                        let copy = catallax_topology::transformed_body(
                            model,
                            source.body,
                            &map,
                            |origin| crate::copy_origin(&step.id, instance, origin),
                        );
                        copy_fillets(&mut fillets, &source_fillets, &map, &step.id, instance)?;
                        for definition in &source_chamfers {
                            let mut definition = definition.clone();
                            definition.transform(&map)?;
                            definition.origin =
                                crate::copy_origin(&step.id, instance, &definition.origin);
                            chamfers.push(definition);
                        }
                        let copy_bounds = source_bounds
                            .as_ref()
                            .map(|bounds| transformed_bounds(bounds, &map, dec))
                            .transpose()?;
                        body = if copy_bounds.as_ref().is_some_and(|copy| {
                            placed_bounds.iter().all(|prior| separated(prior, copy))
                        }) {
                            unite_separated(model, body, copy)
                        } else {
                            boolean(model, body, copy, BoolOp::Unite, dec)
                                .map_err(|e| e.refusal)?
                                .body
                        };
                        if let Some(bounds) = copy_bounds {
                            placed_bounds.push(bounds);
                        }
                    }
                    BodyValue { body, radius }
                }
                Operation::Mirror {
                    body,
                    plane,
                    keep_original,
                } => {
                    let map = crate::mirror_placement(plane)?;
                    let source = &bodies[body];
                    if *keep_original {
                        require_material(model, &[source.body])?;
                    }
                    let source_chamfers: Vec<_> = chamfers
                        .iter()
                        .filter(|d| {
                            model
                                .faces_of_body(source.body)
                                .iter()
                                .any(|f| model.face_origins(*f).contains(&d.origin))
                        })
                        .cloned()
                        .collect();
                    let source_fillets = source_fillet_indices(model, source.body, &fillets);
                    let source_bounds = if *keep_original {
                        catallax_measure::bounding_box(
                            model,
                            source.body,
                            &Default::default(),
                            dec,
                        )?
                        .value
                    } else {
                        None
                    };
                    let shift = [&map.translation.x, &map.translation.y, &map.translation.z]
                        .into_iter()
                        .try_fold(Rational::zero(), |sum, v| {
                            let b = dec.enclose_at(v, 128)?;
                            Ok::<_, Refusal>(sum + b.lo().abs().max(b.hi().abs()))
                        })?;
                    let copy =
                        catallax_topology::transformed_body(model, source.body, &map, |origin| {
                            crate::copy_origin(&step.id, 1, origin)
                        });
                    copy_fillets(&mut fillets, &source_fillets, &map, &step.id, 1)?;
                    for mut definition in source_chamfers {
                        definition.transform(&map)?;
                        definition.origin = crate::copy_origin(&step.id, 1, &definition.origin);
                        chamfers.push(definition);
                    }
                    let body = if *keep_original {
                        let copy_bounds = source_bounds
                            .as_ref()
                            .map(|bounds| transformed_bounds(bounds, &map, dec))
                            .transpose()?;
                        if source_bounds
                            .as_ref()
                            .zip(copy_bounds.as_ref())
                            .is_some_and(|(a, b)| separated(a, b))
                        {
                            unite_separated(model, source.body, copy)
                        } else {
                            boolean(model, source.body, copy, BoolOp::Unite, dec)
                                .map_err(|e| e.refusal)?
                                .body
                        }
                    } else {
                        copy
                    };
                    BodyValue {
                        body,
                        radius: &source.radius + &shift,
                    }
                }
                Operation::Lattice {
                    region,
                    cell,
                    size,
                    radius,
                    ..
                } => {
                    // The lattice is field-only. The B-rep carries the PRD section 7
                    // stand-in: the region is the outer envelope, and the field kernel
                    // adds bodies matching mass, centre of gravity and inertia.
                    crate::lattice::validate(cell, size.value(), radius)?;
                    let source = &bodies[region];
                    stand_ins.push(StandIn {
                        step: step.id.clone(),
                        region: region.clone(),
                    });
                    BodyValue {
                        body: source.body,
                        radius: source.radius.clone(),
                    }
                }
                Operation::Gyroid {
                    region,
                    period,
                    level,
                    ..
                } => {
                    // The gyroid is field-only too, carried by the same region envelope.
                    crate::lattice::validate_gyroid(period.value(), level)?;
                    let source = &bodies[region];
                    stand_ins.push(StandIn {
                        step: step.id.clone(),
                        region: region.clone(),
                    });
                    BodyValue {
                        body: source.body,
                        radius: source.radius.clone(),
                    }
                }
                Operation::Offset { body, distance } if certification.is_some() => {
                    let input = &bodies[body];
                    let (body, radius) = certification
                        .as_deref_mut()
                        .expect("certified branch")
                        .offset(
                            model,
                            &profiles,
                            &step.id,
                            (body, input.body, &input.radius),
                            distance.value(),
                            dec,
                        )?;
                    BodyValue { body, radius }
                }
                Operation::Offset { body, distance } => {
                    let input = &bodies[body];
                    let result = crate::distance_build::offset(
                        model,
                        input.body,
                        &distance.value().clone().into(),
                        &step.id,
                        dec,
                    )?;
                    BodyValue {
                        body: crate::distance_build::inherit_marks(model, input.body, result),
                        radius: &input.radius + &distance.value().abs(),
                    }
                }
                Operation::Shell {
                    body,
                    open,
                    thickness,
                } => {
                    let input = &bodies[body];
                    let result = crate::distance_build::shell(
                        model,
                        input.body,
                        open,
                        &thickness.value().clone().into(),
                        &step.id,
                        dec,
                    )?;
                    BodyValue {
                        body: crate::distance_build::inherit_marks(model, input.body, result),
                        radius: input.radius.clone(),
                    }
                }
                Operation::Thicken {
                    body,
                    faces,
                    thickness,
                    side,
                } => {
                    let input = &bodies[body];
                    let result = crate::distance_build::thicken(
                        model,
                        input.body,
                        faces,
                        &thickness.value().clone().into(),
                        *side,
                        &step.id,
                        recipe.format != catallax_recipe::FORMAT,
                        dec,
                    )?;
                    BodyValue {
                        body: crate::distance_build::inherit_marks(model, input.body, result),
                        radius: &input.radius + &thickness.value().abs(),
                    }
                }
                _ => {
                    return Err(Refusal::new(
                        "features.unsupported-operation",
                        format!(
                            "B-rep recipe operation {} is not implemented",
                            step.op.name()
                        ),
                    ));
                }
            }))
        };
        let result = if let Operation::Sketch { plane, loops } = &step.op {
            Profile::new(plane, loops, dec).map(|profile| {
                profiles.insert(step.id.clone(), profile);
                None
            })
        } else if let Operation::Path { start, segments } = &step.op {
            if numerical_sweep.is_some()
                && crate::general_sweep::path_needs_numerical(recipe, &step.id)
            {
                Ok(None)
            } else {
                crate::sweep::Path::new(start, segments).map(|path| {
                    paths.insert(step.id.clone(), path);
                    None
                })
            }
        } else if let Operation::ScanImport { .. } = &step.op {
            // Reference data, RECIPE.md section 4.19.1: no body, and no body
            // operation may read it, so the B-rep has nothing to build.
            Ok(None)
        } else {
            build(&mut model, dec).map(|body| {
                body.map(|body| {
                    let id = body.body;
                    bodies.insert(step.id.clone(), body);
                    id
                })
            })
        };
        span.finish(
            dec,
            result.as_ref().err(),
            result
                .as_ref()
                .ok()
                .and_then(|body| body.map(|id| &model.body(id).marks)),
        );
        result.map_err(|refusal| {
            let limit_feature = matches!(
                step.op,
                Operation::Fillet { .. }
                    | Operation::FilletVariable { .. }
                    | Operation::Shell { .. }
                    | Operation::Thicken { .. }
            );
            let contact_refusal = refusal.code.starts_with("fillet.")
                || refusal.code.starts_with("shapes.rolling-")
                || refusal.code.starts_with("distance.")
                || refusal.code.starts_with("offset.");
            let refusal = if limit_feature && contact_refusal && refusal.undecided().is_none() {
                crate::refusal_limit::unknown(refusal)
            } else {
                refusal
            };
            BuildError::Geometry {
                step: step.id.clone(),
                refusal,
            }
        })?;
    }
    Ok(EvaluationState {
        model,
        profiles,
        paths,
        bodies,
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
