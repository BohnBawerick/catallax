//! Existing exact stock builders as immutable, named recipe inputs.

use catallax_numbers::{Decider, MarkLog, Rational, Refusal};
use catallax_recipe::{ElementId, Operation, Origin, Role, StepId};
use catallax_topology::{Body, BodyId, Lump, Model, Shell};

fn side(name: &str) -> Role {
    Role::Side {
        segment: ElementId::new(name).expect("fixed stock segment id"),
        path: None,
    }
}

pub(crate) fn build(
    model: &mut Model,
    op: &Operation,
    step: &StepId,
    dec: &mut Decider,
) -> Result<Option<(BodyId, Rational)>, Refusal> {
    let (raw, roles, local_bound, placement) = match op {
        Operation::Block { size, pose } => {
            let sizes = size.each_ref().map(|n| n.value().clone());
            let bound = sizes.iter().fold(Rational::zero(), |a, b| a + b.clone());
            (
                catallax_topology::block(model, sizes, &crate::profile::placement(pose, dec)?)?,
                vec![
                    side("x_min"),
                    side("x_max"),
                    side("y_min"),
                    side("y_max"),
                    Role::Start,
                    Role::End,
                ],
                bound,
                pose,
            )
        }
        Operation::Cylinder {
            radius,
            height,
            pose,
        } => (
            catallax_topology::cylinder(
                model,
                radius.value().clone(),
                height.value().clone(),
                &crate::profile::placement(pose, dec)?,
            )?,
            vec![Role::Start, Role::End, side("wall")],
            radius.value() + height.value(),
            pose,
        ),
        Operation::Sphere { radius, pose } => (
            catallax_topology::sphere(
                model,
                radius.value().clone(),
                &crate::profile::placement(pose, dec)?,
            )?,
            vec![side("wall")],
            radius.value().clone(),
            pose,
        ),
        Operation::Cone {
            base_radius,
            top_radius,
            height,
            pose,
        } => {
            let placement = crate::profile::placement(pose, dec)?;
            let raw = if base_radius == top_radius {
                catallax_topology::cylinder(
                    model,
                    base_radius.value().clone(),
                    height.value().clone(),
                    &placement,
                )?
            } else {
                catallax_topology::cone(
                    model,
                    base_radius.value().clone(),
                    top_radius.value().clone(),
                    height.value().clone(),
                    &placement,
                )?
            };
            let mut roles = vec![Role::Start];
            if !top_radius.value().is_zero() {
                roles.push(Role::End);
            }
            roles.push(side("wall"));
            (
                raw,
                roles,
                base_radius.value().clone().max(top_radius.value().clone())
                    + height.value().clone(),
                pose,
            )
        }
        Operation::Torus { major, minor, pose } => (
            catallax_topology::torus(
                model,
                major.value().clone(),
                minor.value().clone(),
                &crate::profile::placement(pose, dec)?,
            )?,
            vec![side("wall")],
            major.value() + minor.value(),
            pose,
        ),
        _ => return Ok(None),
    };
    let old_faces = model.faces_of_body(raw);
    assert_eq!(old_faces.len(), roles.len(), "stock builder face count");
    let faces = old_faces
        .into_iter()
        .zip(roles)
        .map(|(old, role)| {
            let face = model.face(old).clone();
            model.add_named_face(
                face,
                vec![Origin {
                    step: step.clone(),
                    role,
                }],
            )
        })
        .collect();
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: model.body(raw).kind,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    let bound = local_bound
        + placement
            .translation
            .0
            .iter()
            .fold(Rational::zero(), |a, n| a + n.value().abs());
    Ok(Some((body, bound)))
}
