//! Complete circular-trim validation before convex rotational distance assembly.

use catallax_numbers::{Decider, Refusal};
use catallax_shapes::{Curve, PointDef, Surface};
use catallax_topology::{BodyId, Loop, Model, Sense};

use crate::circular_skin::{CircularSkin, bad};

pub(crate) fn skin(
    model: &mut Model,
    body: BodyId,
    dec: &mut Decider,
) -> Result<Option<CircularSkin>, Refusal> {
    let source = model.faces_of_body(body);
    if !source.iter().any(|f| {
        matches!(
            model.geom().surface(model.face(*f).surface),
            Surface::Cylinder(_)
        )
    }) {
        return Ok(None);
    }
    let surfaces: Vec<_> = source
        .iter()
        .map(|f| {
            (
                model.geom().surface(model.face(*f).surface).clone(),
                model.face_origins(*f).to_vec(),
            )
        })
        .collect();
    let Some(skin) = CircularSkin::new(&surfaces, dec)? else {
        return Ok(None);
    };
    if model.body(body).lumps.len() != 1
        || model.lump(model.body(body).lumps[0]).shells.len() != 1
        || source.len() != skin.faces.len()
    {
        return Err(bad(
            "the convex circular supports do not match one complete input shell",
        ));
    }
    let mut matched = vec![false; source.len()];
    for (i, face) in skin.faces.iter().enumerate() {
        let mut expected = Vec::new();
        for boundary in &skin.boundaries {
            if boundary.faces.contains(&i) {
                let Curve::Circle(circle) = &boundary.curve else {
                    unreachable!()
                };
                let scale = crate::cone_skin::reduced(
                    (&circle.radius_sq / &circle.ref_dir.norm_squared()).sqrt(),
                    dec,
                )?;
                let radial = circle.ref_dir.scaled(&scale);
                let point = model
                    .geom_mut()
                    .add_point(PointDef::Explicit(circle.center.plus(&radial)));
                let tangent = circle.normal.cross(&radial);
                expected.push((
                    model.geom_mut().add_curve(boundary.curve.clone()),
                    point,
                    tangent,
                    boundary.faces[1] == i,
                ));
            }
        }
        let mut found = false;
        for (j, source) in source.iter().enumerate() {
            let f = model.face(*source);
            if matched[j]
                || f.sense != Sense::Same
                || !model
                    .geom()
                    .surface(f.surface)
                    .same_definition(&face.surface)
                || f.loops.len() != expected.len()
            {
                continue;
            }
            let mut curves = Vec::new();
            for lp in &f.loops {
                let Loop::Coedges(cs) = model.loop_(*lp) else {
                    return Err(bad("a circular skin has an unsupported singular trim"));
                };
                let [coedge] = cs.as_slice() else {
                    return Err(bad(
                        "partial circular joins need their finite edge-pipe trim arrangement",
                    ));
                };
                let coedge = model.coedge(*coedge);
                let edge = model.edge(coedge.edge);
                if edge.start.is_some() || edge.end.is_some() {
                    return Err(bad(
                        "circular joins require complete unbroken section curves",
                    ));
                }
                curves.push((edge.curve, coedge.sense == edge.sense));
            }
            let mut used = vec![false; curves.len()];
            for (curve, point, tangent, positive) in &expected {
                let mut same = false;
                for (k, (actual, actual_positive)) in curves.iter().enumerate() {
                    if !used[k]
                        && catallax_intersect::same_curve(model.geom(), *curve, *actual, dec)
                            .is_coincident()
                    {
                        let [actual_tangent, _] =
                            catallax_intersect::tangent_at(model.geom(), *actual, *point, dec)?;
                        let orientation = crate::profile::sign(&actual_tangent.dot(tangent), dec)?;
                        if orientation != catallax_numbers::Sign::Zero
                            && *positive
                                == (*actual_positive
                                    == (orientation == catallax_numbers::Sign::Positive))
                        {
                            used[k] = true;
                            same = true;
                            break;
                        }
                    }
                }
                if !same {
                    break;
                }
            }
            if used.iter().all(|v| *v) {
                matched[j] = true;
                found = true;
                break;
            }
        }
        if !found {
            return Err(bad(
                "a source face does not have the complete active circular trim domain",
            ));
        }
    }
    Ok(Some(skin))
}

pub(crate) fn outward(
    model: &mut Model,
    skin: &CircularSkin,
    distance: &catallax_numbers::Real,
    step: &catallax_recipe::StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let (faces, boundaries) = skin.outward(distance, step, dec)?;
    Ok(crate::cone_build::sections(
        model,
        &faces,
        &boundaries,
        None,
    ))
}
