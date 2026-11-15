//! Distance shell of a convex planar envelope pierced by separated through bores.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{FaceRef, Origin, Role, StepId};
use catallax_shapes::{Surface, Vector3};
use catallax_topology::{BodyId, Model};

use crate::{
    cylinder_skin::CylinderSkin,
    draft::bored::{axis_point, certify_bores},
    polyhedron::{ConvexCell, PlaneFace},
    profile::sign,
};

use super::bad;

#[path = "bored_shell/prismatic.rs"]
mod prismatic;

fn domain(error: Refusal) -> Refusal {
    if error.code.starts_with("draft.") {
        error.with_code("distance.unsupported-envelope")
    } else {
        error
    }
}

// Preserve the established cap-only family. A side opening needs a separate
// prismatic proof; perpendicular caps alone do not certify finite side distances.
fn certify_openings(
    cell: &ConvexCell,
    selected: &[bool],
    axis: &Vector3,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let mut side_open = false;
    for (face, selected) in cell.faces.iter().zip(selected) {
        if *selected && sign(&face.plane.normal.cross(axis).norm_squared(), dec)? != Sign::Zero {
            side_open = true;
        }
    }
    if !side_open {
        return Ok(());
    }
    if selected.iter().filter(|v| **v).count() != 1 {
        return Err(bad(
            "a bored prismatic shell supports one planar side opening",
        ));
    }
    for face in &cell.faces {
        let normal = &face.plane.normal;
        if sign(&normal.cross(axis).norm_squared(), dec)? != Sign::Zero
            && sign(&normal.dot(axis), dec)? != Sign::Zero
        {
            return Err(bad(
                "a side-open bored shell needs every side plane parallel to its bore axes",
            ));
        }
    }
    Ok(())
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
    let inventory = catallax_numbers::OperationSpan::begin(
        "brep",
        &format!("{step}.opening-inventory"),
        "distance.shell.opening-inventory",
        dec,
    );
    let Some(source) =
        crate::draft_build::bored_distance_skin(model, &faces, dec).map_err(domain)?
    else {
        inventory.finish(dec, None, None);
        return Ok(None);
    };
    inventory.finish(dec, None, None);
    let crate::draft_build::BoredDistanceSource {
        skin,
        mut bores,
        openings,
    } = source;
    // The filled skin omits cap holes. It cannot answer a nearest-face tie or
    // adjacency to their circular rims, so do not resolve those on filled trims.
    for reference in open {
        let mut candidates = 0;
        let mut curved = false;
        for face in &faces {
            if reference.matches_origins(model.face_origins(*face)) {
                candidates += 1;
                curved |= !matches!(
                    model.geom().surface(model.face(*face).surface),
                    Surface::Plane(_)
                );
            }
        }
        if !reference.touching.is_empty() || (reference.near.is_some() && candidates > 1) {
            return Err(bad(
                "bored shells need unambiguous planar references without adjacency filters",
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
    certify_openings(&cell, &selected, &bores[0].cylinder.axis, dec)?;
    let cavity = cell.cavity(&selected, thickness, step, dec)?;
    if cavity.is_empty() {
        return Ok(Some(body));
    }
    // Projection to a perpendicular cap preserves the radial coordinate and
    // cannot enter a bore. For a side opening, certify_openings also proves all
    // side planes parallel to the bore axes: their complete supports miss every
    // source bore. Together with cavity's finite-trim projection proof, minimum
    // kept-plane witnesses remain on physical faces. Radial bore witnesses stay
    // between the perpendicular caps and inside the strictly cleared envelope.
    for bore in &mut bores {
        let radius = &bore.cylinder.radius_sq.sqrt() + thickness;
        // Keep exact rational carriers on the analytic integration path instead
        // of making consumers repeatedly prove a rational square-root identity.
        let radius = dec
            .rational_value(&radius)?
            .map(Real::from)
            .unwrap_or(radius);
        bore.cylinder.radius_sq = radius.square();
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
    let direct = if openings.is_empty() {
        None
    } else {
        let span = catallax_numbers::OperationSpan::begin(
            "brep",
            &format!("{step}.prismatic-cavity"),
            "distance.shell.prismatic-cavity",
            dec,
        );
        let marks = model.body(body).marks.clone();
        let result = prismatic::build(
            model,
            &cavity,
            prismatic::Holes {
                bores: &bores,
                polygons: &openings,
                caps,
                marks,
            },
            thickness,
            step,
            dec,
        )?;
        span.finish(dec, None, None);
        result
    };
    let hollow = if let Some(hollow) = direct {
        hollow
    } else {
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
            hollow = catallax_boolean::boolean(
                model,
                hollow,
                tool,
                catallax_boolean::BoolOp::Subtract,
                dec,
            )
            .map_err(|e| e.refusal)?
            .body;
        }
        // Within the original cap interval, distance to a prismatic void equals
        // transverse distance to its side walls. Its closed convex offset supplies
        // planar strips and circular corner sectors. Rounded end patches lie beyond
        // that interval and cannot replace a physical opening rim.
        for (i, opening) in openings.into_iter().enumerate() {
            let offset = catallax_numbers::OperationSpan::begin(
                "brep",
                &format!("{step}.opening-offset.{i}"),
                "distance.shell.opening-offset",
                dec,
            );
            let tool = crate::rounded_polyhedron::build(
                model,
                &opening,
                thickness,
                step,
                model.body(body).marks.clone(),
                dec,
            )?;
            offset.finish(dec, None, None);
            model.map_face_origins(tool, |origin| match origin.role {
                Role::Round { .. } | Role::RoundCorner { .. } => origin.clone(),
                _ => Origin {
                    step: step.clone(),
                    role: Role::Inner {
                        of: Box::new(origin.clone()),
                    },
                },
            });
            let cut = catallax_numbers::OperationSpan::begin(
                "brep",
                &format!("{step}.opening-cavity.{i}"),
                "distance.shell.opening-cavity",
                dec,
            );
            hollow = catallax_boolean::boolean(
                model,
                hollow,
                tool,
                catallax_boolean::BoolOp::Subtract,
                dec,
            )
            .map_err(|error| error.refusal)?
            .body;
            cut.finish(dec, None, None);
        }
        hollow
    };
    let cut = catallax_numbers::OperationSpan::begin(
        "brep",
        &format!("{step}.material"),
        "distance.shell.material",
        dec,
    );
    let result =
        catallax_boolean::boolean(model, body, hollow, catallax_boolean::BoolOp::Subtract, dec)
            .map_err(|e| e.refusal)?
            .body;
    cut.finish(dec, None, None);
    Ok(Some(result))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;
    use catallax_recipe::{Operation, Recipe};
    use serde_json::{Value, json};

    #[test]
    fn side_open_refusals_add_nothing_to_the_input_model() {
        let document: Value =
            serde_json::from_str(include_str!("../tests/fixtures/bored-side-shell.json")).unwrap();
        // Expanded r=5 touches/crosses the open x=0 plane, then the kept
        // inset x=38 plane. The r=3 source bore remains strictly internal.
        for (x, thickness, code) in [
            ("5", "2", "distance.unsupported-envelope"),
            ("4", "2", "distance.unsupported-envelope"),
            ("33", "2", "distance.unsupported-envelope"),
            ("34", "2", "distance.unsupported-envelope"),
            ("pair-touch", "2", "distance.unsupported-envelope"),
            ("pair-cross", "2", "distance.unsupported-envelope"),
            ("20", "0", "distance.non-positive-thickness"),
            ("20", "-1", "distance.non-positive-thickness"),
        ] {
            let mut source = document.clone();
            if x.starts_with("pair-") {
                let right = if x == "pair-touch" { "25" } else { "24" };
                source["steps"][1]["at"] = json!([
                    {"id":"left","point":["15","15","12"],"direction":["0","0","-1"]},
                    {"id":"right","point":[right,"15","12"],"direction":["0","0","-1"]}
                ]);
            } else {
                source["steps"][1]["at"][0]["point"][0] = json!(x);
            }
            source["steps"][2]["thickness"] = json!(thickness);
            let mut recipe = Recipe::parse(&source.to_string()).unwrap();
            let step = recipe.steps.pop().unwrap();
            recipe.result = recipe.steps.last().unwrap().id.clone();
            let (mut built, work) = with_operation_budget(20_000_000, || {
                crate::evaluate(&recipe, &mut Decider::default()).unwrap()
            });
            assert!(!work.exhausted);
            let Operation::Shell {
                open, thickness, ..
            } = step.op
            else {
                unreachable!()
            };
            built.model.journal_begin();
            let error = crate::distance_build::shell(
                &mut built.model,
                built.result,
                &open,
                &thickness.value().clone().into(),
                &step.id,
                &mut Decider::default(),
            )
            .unwrap_err();
            assert_eq!(error.code, code, "x={x}, t={thickness:?}: {error:?}");
            assert!(
                built.model.journal_end().is_empty(),
                "a refusal must publish no model writes"
            );
            assert!(built.model.check(built.result).is_empty());
        }
    }
}
