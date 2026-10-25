//! Direct convex conical offsets with circular edge and spherical tip joins.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{EdgeName, Origin, Role, StepId, VertexName};
use catallax_shapes::{Circle, Curve, Sphere, Surface, Torus};
use catallax_topology::{BodyId, Model};

use crate::{
    cone_skin::{ConeBoundary, ConeFace, ConeSkin, reduced},
    profile::sign,
};

fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.cone-domain", detail)
}

pub(crate) fn build(
    model: &mut Model,
    skin: &ConeSkin,
    distance: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    if skin
        .faces
        .iter()
        .any(|f| matches!(f.surface, Surface::Cylinder(_)))
    {
        return Err(bad(
            "outward cylindrical neighbours need their circular-edge join arrangement",
        ));
    }
    let (wall, cone) = skin
        .faces
        .iter()
        .enumerate()
        .find_map(|(i, f)| match &f.surface {
            Surface::Cone(c) => Some((i, c)),
            _ => None,
        })
        .ok_or_else(|| bad("an outward conical cell needs one conical wall"))?;
    let Some(surface) = skin.faces[wall].surface.parallel_support(distance, dec)? else {
        unreachable!("a cone support does not collapse")
    };
    let mut faces = vec![ConeFace {
        surface,
        origins: skin.faces[wall].origins.clone(),
    }];
    let mut boundaries = Vec::new();
    let factor_sq = &Real::one() + &cone.slope.square();
    let factor = reduced(factor_sq.sqrt(), dec)?;
    let radial = &(distance / &factor_sq) * &factor;
    let axial = -(&radial * &cone.slope);
    if skin.boundaries.len() != if skin.apex.is_some() { 1 } else { 2 } {
        return Err(bad(
            "outward conical joins need a complete circular cap at each finite end",
        ));
    }
    if let Some((apex, incident)) = &skin.apex {
        let origin = skin.faces[*incident]
            .origins
            .first()
            .ok_or_else(|| bad("a conical tip needs a named incident face"))?;
        let vertex = VertexName::new(vec![origin.clone()]).map_err(|e| bad(&e.to_string()))?;
        let tip = faces.len();
        faces.push(ConeFace {
            surface: Surface::Sphere(Sphere {
                center: apex.clone(),
                radius: distance.clone(),
                axis: skin.axis.clone(),
                ref_dir: cone.ref_dir.clone(),
            }),
            origins: vec![Origin {
                step: step.clone(),
                role: Role::RoundCorner { vertex },
            }],
        });
        // The normal at a conical generator is (radial - slope * axis) / factor.
        // Its image of the apex is the exact sphere/cone contact circle.
        boundaries.push(ConeBoundary {
            curve: Curve::Circle(Circle {
                center: apex.plus(&skin.axis.scaled(&axial)),
                normal: skin.axis.clone(),
                radius_sq: reduced(&distance.square() / &factor_sq, dec)?,
                ref_dir: cone.ref_dir.clone(),
            }),
            faces: [tip, 0],
        });
    }
    for boundary in &skin.boundaries {
        let lower = boundary.faces[0] != wall;
        let Curve::Circle(circle) = &boundary.curve else {
            return Err(bad(
                "elliptic conical edges need a certified non-toroidal distance carrier",
            ));
        };
        if sign(&circle.normal.cross(&skin.axis).norm_squared(), dec)? != Sign::Zero {
            return Err(bad(
                "an outward conical edge needs a perpendicular circular section",
            ));
        }
        let cap = boundary
            .faces
            .iter()
            .copied()
            .find(|i| *i != wall)
            .ok_or_else(|| bad("a conical edge has no cap neighbour"))?;
        let Surface::Plane(plane) = &skin.faces[cap].surface else {
            return Err(bad("a circular conical edge needs its planar neighbour"));
        };
        if sign(&plane.normal.cross(&skin.axis).norm_squared(), dec)? != Sign::Zero {
            return Err(bad(
                "an oblique conical cap needs a non-toroidal edge-distance carrier",
            ));
        }
        let radius = reduced(circle.radius_sq.sqrt(), dec)?;
        if sign(&(&radius - distance), dec)? != Sign::Positive {
            return Err(Refusal::new(
                "offset.spindle-carrier",
                "the circular-edge distance tube reaches a horn or spindle torus",
            ));
        }
        let a = skin.faces[boundary.faces[0]]
            .origins
            .first()
            .ok_or_else(|| bad("a conical edge needs named incident faces"))?;
        let b = skin.faces[boundary.faces[1]]
            .origins
            .first()
            .ok_or_else(|| bad("a conical edge needs named incident faces"))?;
        let edge = EdgeName::new(a.clone(), b.clone()).map_err(|e| bad(&e.to_string()))?;
        let round = faces.len();
        faces.push(ConeFace {
            surface: Surface::Torus(Torus {
                center: circle.center.clone(),
                axis: skin.axis.clone(),
                major: radius.clone(),
                minor: distance.clone(),
                ref_dir: circle.ref_dir.clone(),
            }),
            origins: vec![Origin {
                step: step.clone(),
                role: Role::Round { edge },
            }],
        });
        let Some(surface) = skin.faces[cap].surface.parallel_support(distance, dec)? else {
            unreachable!("a plane support does not collapse")
        };
        let cap_face = faces.len();
        faces.push(ConeFace {
            surface,
            origins: skin.faces[cap].origins.clone(),
        });
        let side = Curve::Circle(Circle {
            center: circle.center.plus(&skin.axis.scaled(&axial)),
            normal: skin.axis.clone(),
            radius_sq: (&radius + &radial).square(),
            ref_dir: circle.ref_dir.clone(),
        });
        let end = Curve::Circle(Circle {
            center: circle.center.plus(
                &skin
                    .axis
                    .scaled(&(distance * &Real::from_i64(if lower { -1 } else { 1 }))),
            ),
            normal: skin.axis.clone(),
            radius_sq: circle.radius_sq.clone(),
            ref_dir: circle.ref_dir.clone(),
        });
        boundaries.push(ConeBoundary {
            curve: side,
            faces: if lower { [round, 0] } else { [0, round] },
        });
        boundaries.push(ConeBoundary {
            curve: end,
            faces: if lower {
                [cap_face, round]
            } else {
                [round, cap_face]
            },
        });
    }
    Ok(crate::cone_build::sections(
        model,
        &faces,
        &boundaries,
        None,
    ))
}
