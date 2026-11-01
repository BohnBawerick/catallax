//! Distance shell of a convex planar envelope pierced by separated through bores.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{FaceRef, Origin, Role, StepId};
use catallax_shapes::Surface;
use catallax_topology::{BodyId, Model};

use crate::{
    cylinder_skin::CylinderSkin,
    draft::bored::{axis_point, certify_bores},
    polyhedron::{ConvexCell, PlaneFace},
    profile::sign,
};

use super::bad;

fn domain(error: Refusal) -> Refusal {
    if error.code.starts_with("draft.") {
        bad(&error.detail)
    } else {
        error
    }
}

pub(super) fn shell(
    model: &mut Model,
    body: BodyId,
    open: &[FaceRef],
    thickness: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<Option<BodyId>, Refusal> {
    let faces = model.faces_of_body(body);
    let Some((skin, mut bores)) =
        crate::draft_build::bored_distance_skin(model, &faces, dec).map_err(domain)?
    else {
        return Ok(None);
    };
    // The filled skin omits cap holes. It cannot answer a nearest-face tie or
    // adjacency to their circular rims, so do not resolve those on filled trims.
    for reference in open {
        let mut candidates = 0;
        let mut curved = false;
        for face in &faces {
            if model
                .face_origins(*face)
                .iter()
                .any(|o| reference.face.matches(o))
            {
                candidates += 1;
                curved |= !matches!(
                    model.geom().surface(model.face(*face).surface),
                    Surface::Plane(_)
                );
            }
        }
        if !reference.touching.is_empty() || (reference.near.is_some() && candidates > 1) {
            return Err(bad(
                "bored shells need unambiguous cap references without adjacency filters",
            ));
        }
        if curved {
            return Err(bad(
                "opening a bore wall needs its finite edge-distance envelopes",
            ));
        }
    }
    let cell = ConvexCell::new(
        skin.patches
            .into_iter()
            .map(|patch| PlaneFace {
                plane: patch.plane,
                origins: patch.origins,
            })
            .collect(),
        dec,
    )?;
    let selected = cell.select(open, dec)?;
    for (face, selected) in cell.faces.iter().zip(&selected) {
        if *selected
            && sign(
                &face
                    .plane
                    .normal
                    .cross(&bores[0].cylinder.axis)
                    .norm_squared(),
                dec,
            )? != Sign::Zero
        {
            return Err(bad("bored shells open only their perpendicular end caps"));
        }
    }
    let cavity = cell.cavity(&selected, thickness, step, dec)?;
    if cavity.is_empty() {
        return Ok(Some(body));
    }
    // Inside the original material, projection to a perpendicular cap cannot
    // enter a bore. Kept outer-face witnesses therefore remain on physical faces;
    // each bore contributes the radial distance to its complete finite wall.
    for bore in &mut bores {
        bore.cylinder.radius_sq = (&bore.cylinder.radius_sq.sqrt() + thickness).square();
        bore.origins = bore
            .origins
            .iter()
            .map(|of| Origin {
                step: step.clone(),
                role: Role::Inner {
                    of: Box::new(of.clone()),
                },
            })
            .collect();
    }
    let cylinders: Vec<_> = bores.iter().map(|bore| &bore.cylinder).collect();
    // Strict clearance keeps every new bore a complete cylinder between two caps.
    // Touching envelopes require a different arrangement, not clipped radii.
    let caps =
        certify_bores(cavity.faces.iter().map(|f| &f.plane), &cylinders, dec).map_err(domain)?;
    let mut hollow = cavity.build(model);
    for bore in bores {
        let mut ends = caps.map(|i| axis_point(&bore.cylinder, &cavity.faces[i].plane));
        if sign(&ends[1].minus(&ends[0]).dot(&bore.cylinder.axis), dec)? == Sign::Negative {
            ends.swap(0, 1);
        }
        let tool = super::cylinder(
            model,
            &CylinderSkin {
                cylinder: bore.cylinder,
                ends,
                origins: [Vec::new(), Vec::new(), bore.origins],
            },
        );
        hollow =
            catallax_boolean::boolean(model, hollow, tool, catallax_boolean::BoolOp::Subtract, dec)
                .map_err(|e| e.refusal)?
                .body;
    }
    Ok(Some(
        catallax_boolean::boolean(model, body, hollow, catallax_boolean::BoolOp::Subtract, dec)
            .map_err(|e| e.refusal)?
            .body,
    ))
}
