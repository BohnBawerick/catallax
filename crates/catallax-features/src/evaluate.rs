//! Ordered feature evaluation, with bounds derived from the recipe definitions.

use crate::{Profile, extrude, revolve};
use catallax_boolean::{BoolOp, boolean};
use catallax_numbers::{Ball, Decider, Rational, Refusal};
use catallax_recipe::{BooleanKind, Operation, Recipe, RecipeError, StepId};
use catallax_shapes::{RPoint3, RVector3};
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
    definition: Arc<Recipe>,
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
fn profile_radius(profile: &Profile) -> Rational {
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

fn through(
    model: &mut Model,
    profile: &Profile,
    step: &StepId,
    target: &BodyValue,
    direction: Option<&catallax_recipe::Vector3>,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let normal = profile.pose.apply_vector(&RVector3::new(
        Rational::zero(),
        Rational::zero(),
        Rational::one(),
    ));
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
    let origin = profile.pose.apply_point(&RPoint3::origin());
    let normal_l1 = normal.x.abs() + normal.y.abs() + normal.z.abs();
    // For every target point p, |(p-o).n| <= (radius + |o|_1)|n|_1.
    // The strict extra unit leaves both construction caps outside the target.
    // Rescaling direction changes t inversely, not the resulting through-all set.
    let extent =
        ((&target.radius + &norm_upper(&origin)) * normal_l1) / dot.abs() + Rational::one();
    let offset = direction.scaled(&extent);
    let mut shifted = profile.clone();
    shifted.pose.translation = shifted.pose.translation.minus(&offset);
    let tool = crate::extrude::extrude_vector(
        model,
        &shifted,
        step,
        &offset.scaled(&Rational::from_i64(2)),
    )?;
    Ok(boolean(model, target.body, tool, BoolOp::Subtract, dec)
        .map_err(|e| e.refusal)?
        .body)
}

/// Builds sketches, extrusions, revolutions, through-all cuts and explicit booleans.
/// Other recipe operations retain their unfinished-operation refusal, not a partial body.
pub fn evaluate(recipe: &Recipe, dec: &mut Decider) -> Result<Built, BuildError> {
    recipe.check().map_err(BuildError::Recipe)?;
    let mut model = Model::new();
    let mut profiles = BTreeMap::new();
    let mut bodies: BTreeMap<StepId, BodyValue> = BTreeMap::new();
    for step in &recipe.steps {
        let build = |model: &mut Model, dec: &mut Decider| -> Result<Option<BodyValue>, Refusal> {
            if let Some((body, radius)) = crate::stock::build(model, &step.op, &step.id)? {
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
                } => {
                    let target = &bodies[target];
                    let body = through(
                        model,
                        &profiles[profile],
                        &step.id,
                        target,
                        direction.as_ref(),
                        dec,
                    )?;
                    BodyValue {
                        body,
                        radius: target.radius.clone(),
                    }
                }
                Operation::Boolean { kind, target, tool } => {
                    let a = &bodies[target];
                    let b = &bodies[tool];
                    let (op, radius) = match kind {
                        BooleanKind::Union => {
                            (BoolOp::Unite, a.radius.clone().max(b.radius.clone()))
                        }
                        BooleanKind::Subtract => (BoolOp::Subtract, a.radius.clone()),
                        BooleanKind::Intersect => {
                            (BoolOp::Intersect, a.radius.clone().min(b.radius.clone()))
                        }
                    };
                    BodyValue {
                        body: boolean(model, a.body, b.body, op, dec)
                            .map_err(|e| e.refusal)?
                            .body,
                        radius,
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
        if let Operation::Sketch { plane, loops } = &step.op {
            let profile =
                Profile::new(plane, loops, dec).map_err(|refusal| BuildError::Geometry {
                    step: step.id.clone(),
                    refusal,
                })?;
            profiles.insert(step.id.clone(), profile);
        } else {
            let body = build(&mut model, dec)
                .map_err(|refusal| BuildError::Geometry {
                    step: step.id.clone(),
                    refusal,
                })?
                .expect("body-producing operation");
            bodies.insert(step.id.clone(), body);
        }
    }
    let result = bodies[&recipe.result].body;
    Ok(Built {
        model,
        result,
        bodies: bodies
            .into_iter()
            .map(|(id, value)| (id, value.body))
            .collect(),
        definition: Arc::new(recipe.clone()),
    })
}
