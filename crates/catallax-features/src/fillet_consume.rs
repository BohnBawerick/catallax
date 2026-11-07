//! Roll-over affected regions of straight plane/plane fillets whose contacts
//! reach the far boundary of a neighbouring face.
//!
//! Section 4.5.1 keeps a contact on its source face until its foot reaches a
//! trimmed boundary. When every foot stays in its closed face for the whole
//! radius homotopy, the affected region is the exact tangent-triangle prism
//! along the edge, cut by the planes of the two end walls. A face whose whole
//! trimmed set lies in that region is consumed. The regularised Boolean then
//! gives the definition's result directly, including a face reduced to the
//! contact line at the limit radius.
//!
//! Only this stratum is covered: one selected straight edge between planes,
//! ending at two three-face vertices whose third faces are planes normal to
//! the edge. Every other face must be strictly clear of the region's hull or
//! touch it only along a contact rail. A foot that leaves its face before the
//! full radius needs neighbour continuation and keeps refusing.

use std::collections::BTreeSet;

use catallax_numbers::{Decider, Refusal, Sign};
use catallax_shapes::{Curve, Point3, Surface};
use catallax_topology::{BodyId, FaceId, Model};

use super::{Run, Uses, fillet_setback};
use crate::{ConsumedFace, fillet::hull_clear, profile::sign};

/// The plane end wall at one endpoint of the run's edge: the only face at
/// that vertex besides the edge's own two, normal to the edge.
fn end_wall(
    model: &Model,
    run: &Run,
    uses: &Uses,
    at: catallax_topology::VertexId,
    dec: &mut Decider,
) -> Result<Option<FaceId>, Refusal> {
    let faces: BTreeSet<_> = uses
        .iter()
        .filter(|(edge, _)| {
            let e = model.edge(**edge);
            e.start == Some(at) || e.end == Some(at)
        })
        .flat_map(|(_, incident)| incident.iter().map(|(f, _)| *f))
        .filter(|f| !run.faces.contains(f))
        .collect();
    let [wall] = faces.into_iter().collect::<Vec<_>>()[..] else {
        return Ok(None);
    };
    let (Surface::Plane(plane), Curve::Line(edge)) = (
        model.geom().surface(model.face(wall).surface),
        &run.definition.edge,
    ) else {
        return Ok(None);
    };
    let normal = plane.normal.cross(&edge.dir);
    if sign(&normal.norm_squared(), dec)? != Sign::Zero {
        return Ok(None);
    }
    Ok(Some(wall))
}

/// Whether a face's carrier meets the closed hull at most along one contact
/// rail: a plane with every hull point on one closed side, whose points on
/// it all belong to the same rail.
fn touches_rail(hull: &[Point3; 6], surface: &Surface, dec: &mut Decider) -> Result<bool, Refusal> {
    let Surface::Plane(plane) = surface else {
        return Ok(false);
    };
    let mut sides = Vec::with_capacity(6);
    for p in hull {
        sides.push(sign(&p.minus(&plane.origin).dot(&plane.normal), dec)?);
    }
    if sides.contains(&Sign::Positive) && sides.contains(&Sign::Negative) {
        return Ok(false);
    }
    let on: Vec<_> = (0..6).filter(|i| sides[*i] == Sign::Zero).collect();
    // Hull points 1 and 4 lie on the first rail, 2 and 5 on the second.
    Ok(!on.is_empty() && (on.iter().all(|i| i % 3 == 1) || on.iter().all(|i| i % 3 == 2)))
}

/// Builds the run's roll-over region and its result, or `None` when the run
/// is outside this stratum. Consumed input faces are reported with the blend
/// origin.
pub(super) fn build(
    model: &mut Model,
    body: BodyId,
    run: &Run,
    uses: &Uses,
    dec: &mut Decider,
) -> Result<Option<(BodyId, Vec<ConsumedFace>)>, Refusal> {
    let definition = &run.definition;
    let ([Surface::Plane(_), Surface::Plane(_)], Curve::Line(_), [Curve::Line(_), Curve::Line(_)]) =
        (&definition.supports, &definition.edge, &definition.contacts)
    else {
        return Ok(None);
    };
    let input = model.edge(run.edge).clone();
    let (Some(start), Some(end)) = (input.start, input.end) else {
        return Ok(None);
    };
    if start == end {
        return Ok(None);
    }
    let mut walls = Vec::with_capacity(2);
    for at in [start, end] {
        let Some(wall) = end_wall(model, run, uses, at, dec)? else {
            return Ok(None);
        };
        walls.push(wall);
    }
    if walls[0] == walls[1] {
        return Ok(None);
    }
    let stations = [start, end].map(|v| model.geom().position(model.vertex(v).point));
    let along = stations[1].minus(&stations[0]);
    if sign(&along.norm_squared(), dec)? != Sign::Positive {
        return Ok(None);
    }
    let Some(hull) = definition.cut_hull([(&stations[0], &along), (&stations[1], &along)], dec)?
    else {
        return Ok(None);
    };
    // Any boundary of the run's faces inside a contact band belongs to some
    // other face whose carrier meets the hull away from the far rails, so
    // these proofs keep every foot in its closed face for every radius up to
    // the requested one.
    for other in model.faces_of_body(body) {
        if run.faces.contains(&other) || walls.contains(&other) {
            continue;
        }
        let surface = model.geom().surface(model.face(other).surface).clone();
        if !hull_clear(&hull, &surface, dec)? && !touches_rail(&hull, &surface, dec)? {
            return Ok(None);
        }
    }
    for section in [0, 3] {
        for (offset, face) in run.faces.iter().enumerate() {
            if !catallax_boolean::face_contains(model, *face, &hull[section + offset + 1], dec)? {
                return Ok(None);
            }
        }
    }
    let turn = sign(
        &hull[1]
            .minus(&hull[0])
            .cross(&hull[2].minus(&hull[0]))
            .dot(&along),
        dec,
    )?;
    if turn == Sign::Zero {
        return Ok(None);
    }
    let caps = [walls[0], walls[1]].map(|wall| model.face_origins(wall).to_vec());
    let region = fillet_setback::region(
        model,
        run.faces,
        definition,
        &hull,
        [&stations[0], &stations[1]],
        &along,
        turn,
        caps,
    )?;
    let operation = if definition.convex {
        catallax_boolean::BoolOp::Subtract
    } else {
        catallax_boolean::BoolOp::Unite
    };
    let result = catallax_boolean::boolean(model, body, region, operation, dec)
        .map_err(|failure| failure.refusal)?
        .body;
    let surviving: BTreeSet<_> = model
        .faces_of_body(result)
        .into_iter()
        .flat_map(|face| model.face_origins(face).to_vec())
        .collect();
    let mut consumed: Vec<_> = model
        .faces_of_body(body)
        .into_iter()
        .map(|face| model.face_origins(face).to_vec())
        .filter(|source| source.iter().all(|origin| !surviving.contains(origin)))
        .map(|source| ConsumedFace {
            step: definition.origin.step.clone(),
            source,
            by: vec![definition.origin.clone()],
        })
        .collect();
    consumed.sort_by(|a, b| a.source.cmp(&b.source));
    consumed.dedup_by(|a, b| a.source == b.source);
    Ok(Some((result, consumed)))
}
