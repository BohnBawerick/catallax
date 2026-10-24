//! Ordered feature evaluation, with bounds derived from the recipe definitions.

use crate::{Profile, extrude, revolve};
use catallax_boolean::{BoolOp, boolean};
use catallax_measure::BoundingBox;
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal};
use catallax_recipe::{BooleanKind, Operation, Recipe, RecipeError, Role, StepId, ThroughSide};
use catallax_shapes::{Isometry, RPoint3, RVector3};
use catallax_topology::{BodyId, BodyKind, Model};
use std::{collections::BTreeMap, fmt, sync::Arc};

/// A complete build. Intermediate input bodies remain immutable and queryable.
#[derive(Clone, Debug)]
pub struct Built {
    /// Geometry, topology and canonical face origin sets.
    pub model: Model,
    /// The requested result.
    pub result: BodyId,
    /// Every body-producing step, in canonical step-id order.
    pub bodies: BTreeMap<StepId, BodyId>,
    /// Generating definitions of every elementary fillet and named copy in build order.
    pub fillets: Vec<crate::fillet::FilletDefinition>,
    /// Canonical removed-face history from draft operations, not selectable names.
    pub consumed: Vec<crate::ConsumedFace>,
    /// Generating definitions of every elementary chamfer and named copy.
    pub chamfers: Vec<crate::chamfer::ChamferDefinition>,
    pub(crate) definition: Arc<Recipe>,
}

impl Built {
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
            if let Some((c, radius_sq, _)) = &s.circle {
                let radius = Ball::from_rational(radius_sq, 64)
                    .sqrt()
                    .expect("positive circle radius")
                    .hi();
                norm_upper(&profile.world(c, &zero)) + radius
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
    )?;
    Ok(boolean(model, target.body, tool, BoolOp::Subtract, dec)
        .map_err(|e| e.refusal)?
        .body)
}

/// Builds stock solids, sketches, extrusions, revolutions, elementary fillets,
/// cuts, booleans, patterns and mirrors.
/// Other recipe operations retain their unfinished-operation refusal, not a partial body.
pub fn evaluate(recipe: &Recipe, dec: &mut Decider) -> Result<Built, BuildError> {
    evaluate_impl(recipe, None, dec)
}

/// Opts into the recorded 1 nm modelling policy and the first complete bounded
/// offset family. Legacy recipe evaluation remains exact-only and unchanged.
pub fn evaluate_certified(
    recipe: &Recipe,
    dec: &mut Decider,
) -> Result<crate::CertifiedBuilt, BuildError> {
    let mut certification = crate::certified::CertificationContext::new(recipe);
    let built = evaluate_impl(recipe, Some(&mut certification), dec)?;
    certification
        .finish(built, dec)
        .map_err(|refusal| BuildError::Geometry {
            step: recipe.result.clone(),
            refusal,
        })
}

fn evaluate_impl(
    recipe: &Recipe,
    mut certification: Option<&mut crate::certified::CertificationContext<'_>>,
    dec: &mut Decider,
) -> Result<Built, BuildError> {
    recipe.check().map_err(BuildError::Recipe)?;
    let mut model = Model::new();
    let mut profiles = BTreeMap::new();
    let mut bodies: BTreeMap<StepId, BodyValue> = BTreeMap::new();
    let mut fillets = Vec::new();
    let mut consumed = Vec::new();
    let mut chamfers = Vec::new();
    for step in &recipe.steps {
        let span = catallax_numbers::OperationSpan::begin("brep", &step.id, step.op.name(), dec);
        let mut build = |model: &mut Model,
                         dec: &mut Decider|
         -> Result<Option<BodyValue>, Refusal> {
            if let Some((body, radius)) = crate::stock::build(model, &step.op, &step.id, dec)? {
                return Ok(Some(BodyValue { body, radius }));
            }
            Ok(Some(match &step.op {
                Operation::Sketch { .. } => return Ok(None),
                Operation::Extrude {
                    profile,
                    distance,
                    direction,
                } => {
                    let profile = &profiles[profile];
                    let direction = direction.as_ref().map(vector);
                    let body = extrude(
                        model,
                        profile,
                        &step.id,
                        distance.value(),
                        direction.as_ref(),
                    )?;
                    BodyValue {
                        body,
                        radius: profile_radius(profile) + distance.value().abs(),
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
                    if !profiles[profile].is_empty()
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
                    if !profiles[profile].is_empty()
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
                            model.add_body(catallax_topology::Body {
                                kind: BodyKind::Empty,
                                lumps: Vec::new(),
                                marks: model.body(a.body).marks.clone(),
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
                } => {
                    let target = &bodies[body];
                    let (body, definitions) = crate::fillet_build::build(
                        model,
                        target.body,
                        edges,
                        &radius.value().clone().into(),
                        &step.id,
                        dec,
                    )?;
                    fillets.extend(definitions);
                    BodyValue {
                        body,
                        radius: &target.radius + &(radius.value() * &Rational::from_i64(2)),
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
                    let result = crate::distance_build::offset(model, input.body, &distance.value().clone().into(), &step.id, dec)?;
                    BodyValue {
                        body: crate::distance_build::inherit_marks(model, input.body, result),
                        radius: &input.radius + &distance.value().abs(),
                    }
                }
                Operation::Shell { body, open, thickness } => {
                    let input = &bodies[body];
                    let result = crate::distance_build::shell(model, input.body, open, &thickness.value().clone().into(), &step.id, dec)?;
                    BodyValue {
                        body: crate::distance_build::inherit_marks(model, input.body, result),
                        radius: input.radius.clone(),
                    }
                }
                Operation::Thicken { body, faces, thickness, side } => {
                    let input = &bodies[body];
                    let result = crate::distance_build::thicken(model, input.body, faces, &thickness.value().clone().into(), *side, &step.id, dec)?;
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
        } else {
            build(&mut model, dec).map(|body| {
                let body = body.expect("body-producing operation");
                let id = body.body;
                bodies.insert(step.id.clone(), body);
                Some(id)
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
        result.map_err(|refusal| BuildError::Geometry {
            step: step.id.clone(),
            refusal,
        })?;
    }
    let result = bodies[&recipe.result].body;
    Ok(Built {
        model,
        result,
        bodies: bodies
            .into_iter()
            .map(|(id, value)| (id, value.body))
            .collect(),
        fillets,
        chamfers,
        definition: certification
            .and_then(|context| context.snapshot())
            .unwrap_or_else(|| Arc::new(recipe.clone())),
        consumed,
    })
}
