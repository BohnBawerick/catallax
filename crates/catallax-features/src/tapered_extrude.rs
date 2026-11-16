//! Tapered extrusion through the existing exact draft support arrangement.

use crate::{Profile, draft::DraftDefinition};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{FaceRef, Role, StepId};
use catallax_shapes::{Plane, RVector3};
use catallax_topology::{BodyId, Model};
use std::collections::BTreeSet;

/// Extrudes normal to the sketch, then drafts every side about the starting plane.
/// Positive angles narrow the material along the signed extrusion direction.
/// The support arrangement keeps exact algebraic trigonometry and source origins.
/// This uses the same sharp-corner semantics as an Extrude followed by Draft.
pub fn extrude_tapered(
    model: &mut Model,
    profile: &Profile,
    step: &StepId,
    distance: &Rational,
    angle_degrees: &Rational,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let normal = profile
        .pose
        .rotation
        .apply(&RVector3::new(Rational::zero(), Rational::zero(), Rational::one()).to_real());
    let definition = DraftDefinition::new(
        normal.scaled(&Real::from_i64(if distance.sign() == Sign::Negative {
            -1
        } else {
            1
        })),
        Plane {
            origin: profile
                .pose
                .apply_point(&catallax_shapes::RPoint3::origin().to_real()),
            normal,
            ref_dir: profile.pose.rotation.apply(
                &RVector3::new(Rational::one(), Rational::zero(), Rational::zero()).to_real(),
            ),
        },
        angle_degrees,
        dec,
    )?;
    let body = crate::extrude::extrude_with_decider(model, profile, step, distance, None, dec)?;
    if angle_degrees.is_zero() || profile.is_empty() {
        return Ok(body);
    }
    let origins: BTreeSet<_> = model
        .faces_of_body(body)
        .into_iter()
        .flat_map(|face| model.face_origins(face).iter())
        .filter(|origin| matches!(origin.role, Role::Side { .. }))
        .cloned()
        .collect();
    let selected: Vec<_> = origins
        .into_iter()
        .map(|origin| FaceRef {
            face: Some(origin.into()),
            touching: Vec::new(),
            near: None,
        })
        .collect();
    let (body, _, _) =
        crate::draft_build::build(model, body, &selected, &definition, None, step, dec)?;
    Ok(body)
}
