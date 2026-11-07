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

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{Curve, Point3, PointDef, Surface};
use catallax_topology::{BodyId, FaceId, Loop, Model};

use super::{Run, Uses, convex, fillet_setback, on_trimmed};
use crate::{
    ConsumedFace,
    fillet::{FilletDefinition, hull_clear},
    profile::sign,
};

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

/// The contact rail along which a face's carrier meets the closed hull, if
/// it meets it only there: a plane with every hull point on one closed side,
/// whose points on it all belong to that rail. Rail 0 lies on the run's first
/// face, rail 1 on its second.
fn touches_rail(
    hull: &[Point3; 6],
    surface: &Surface,
    dec: &mut Decider,
) -> Result<Option<usize>, Refusal> {
    let Surface::Plane(plane) = surface else {
        return Ok(None);
    };
    let mut sides = Vec::with_capacity(6);
    for p in hull {
        sides.push(sign(&p.minus(&plane.origin).dot(&plane.normal), dec)?);
    }
    if sides.contains(&Sign::Positive) && sides.contains(&Sign::Negative) {
        return Ok(None);
    }
    let on: Vec<_> = (0..6).filter(|i| sides[*i] == Sign::Zero).collect();
    // Hull points 1 and 4 lie on the first rail, 2 and 5 on the second.
    for rail in 0..2 {
        if !on.is_empty() && on.iter().all(|i| i % 3 == rail + 1) {
            return Ok(Some(rail));
        }
    }
    Ok(None)
}

/// The smallest fraction of the requested radius at which a foot of the
/// run reaches a straight boundary edge of its face parallel to the selected
/// edge, when that comes before the full radius.
fn event_fraction(model: &Model, run: &Run, dec: &mut Decider) -> Result<Option<Real>, Refusal> {
    let Curve::Line(edge) = &run.definition.edge else {
        return Ok(None);
    };
    let mut least: Option<Real> = None;
    for (face, rail) in run.faces.iter().zip(&run.definition.contacts) {
        let Curve::Line(rail) = rail else {
            return Ok(None);
        };
        let offset = rail.origin.minus(&edge.origin);
        let offset = offset.minus(
            &edge
                .dir
                .scaled(&(offset.dot(&edge.dir) / edge.dir.norm_squared())),
        );
        let length = offset.norm_squared();
        if sign(&length, dec)? != Sign::Positive {
            return Ok(None);
        }
        for lp in &model.face(*face).loops {
            let Loop::Coedges(coedges) = model.loop_(*lp) else {
                continue;
            };
            for coedge in coedges {
                let id = model.coedge(*coedge).edge;
                if id == run.edge {
                    continue;
                }
                let Curve::Line(line) = model.geom().curve(model.edge(id).curve) else {
                    continue;
                };
                if sign(&line.dir.cross(&edge.dir).norm_squared(), dec)? != Sign::Zero {
                    continue;
                }
                let fraction = &line.origin.minus(&edge.origin).dot(&offset) / &length;
                if sign(&fraction, dec)? != Sign::Positive {
                    continue;
                }
                least = Some(match least {
                    Some(least) if sign(&(&least - &fraction), dec)? != Sign::Positive => least,
                    _ => fraction,
                });
            }
        }
    }
    Ok(match least {
        Some(fraction) if sign(&(&fraction - &Real::one()), dec)? == Sign::Negative => {
            Some(fraction)
        }
        _ => None,
    })
}

/// Builds the region at the first boundary event of a run whose foot leaves
/// its face before the requested radius, or `None` outside that stratum.
///
/// At the event every foot that reaches its boundary lies along a whole
/// straight edge shared with a plane face, at a convex, non-tangent edge.
/// Replacing the face by that neighbour cannot continue the centre/foot
/// limit: the foot normal of the old face is not normal to the new one, and
/// a ball inside the material cannot rest on a convex edge. With no
/// continuation, section 4.5.1 keeps the union of the regions up to the
/// event, so the result uses the event radius and records the request.
pub(super) fn terminal(
    model: &mut Model,
    body: BodyId,
    run: &Run,
    uses: &Uses,
    dec: &mut Decider,
) -> Result<Option<(BodyId, Vec<ConsumedFace>, FilletDefinition)>, Refusal> {
    let requested = &run.definition;
    let ([Surface::Plane(_), Surface::Plane(_)], Curve::Line(_)) =
        (&requested.supports, &requested.edge)
    else {
        return Ok(None);
    };
    let Some(fraction) = event_fraction(model, run, dec)? else {
        return Ok(None);
    };
    let mut definition = FilletDefinition::new(
        requested.origin.clone(),
        requested.supports.clone(),
        requested.same_sense,
        requested.edge.clone(),
        &requested.radius * &fraction,
        requested.convex,
        dec,
    )?;
    definition.requested_radius = Some(requested.radius.clone());
    let Some(built) = build_with(model, body, run, &definition, true, uses, dec)? else {
        return Ok(None);
    };
    Ok(Some((built.0, built.1, definition)))
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
    build_with(model, body, run, &run.definition, false, uses, dec)
}

/// Builds the region of `definition` along the run's edge. At a `terminal`
/// event, every rail that ends on a neighbour must lie along a whole convex,
/// non-tangent edge shared with it, and at least one must.
#[allow(clippy::too_many_arguments)]
fn build_with(
    model: &mut Model,
    body: BodyId,
    run: &Run,
    definition: &FilletDefinition,
    terminal: bool,
    uses: &Uses,
    dec: &mut Decider,
) -> Result<Option<(BodyId, Vec<ConsumedFace>)>, Refusal> {
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
    let mut touching = Vec::new();
    for other in model.faces_of_body(body) {
        if run.faces.contains(&other) || walls.contains(&other) {
            continue;
        }
        let surface = model.geom().surface(model.face(other).surface).clone();
        if hull_clear(&hull, &surface, dec)? {
            continue;
        }
        let Some(rail) = touches_rail(&hull, &surface, dec)? else {
            return Ok(None);
        };
        touching.push((other, rail));
    }
    if terminal {
        if touching.is_empty() {
            return Ok(None);
        }
        for (other, rail) in &touching {
            if !ends_on(model, run, &hull, *other, *rail, uses, dec)? {
                return Ok(None);
            }
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

/// Whether the run's rail `rail` at `hull` lies along one whole straight
/// edge between its face and `other`, a plane meeting it at a convex,
/// non-tangent edge.
fn ends_on(
    model: &mut Model,
    run: &Run,
    hull: &[Point3; 6],
    other: FaceId,
    rail: usize,
    uses: &Uses,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let face = run.faces[rail];
    let (Surface::Plane(own), Surface::Plane(next)) = (
        model.geom().surface(model.face(face).surface),
        model.geom().surface(model.face(other).surface),
    ) else {
        return Ok(false);
    };
    if sign(&own.normal.cross(&next.normal).norm_squared(), dec)? == Sign::Zero {
        return Ok(false);
    }
    let shared: Vec<_> = uses
        .iter()
        .filter(|(_, incident)| {
            incident.iter().any(|(f, _)| *f == face) && incident.iter().any(|(f, _)| *f == other)
        })
        .map(|(edge, incident)| (*edge, incident.clone()))
        .collect();
    let ends = [&hull[rail + 1], &hull[rail + 4]]
        .map(|p| model.geom_mut().add_point(PointDef::Explicit(p.clone())));
    for (edge, incident) in shared {
        if !matches!(model.geom().curve(model.edge(edge).curve), Curve::Line(_)) {
            continue;
        }
        let mut covers = true;
        for end in ends {
            covers &= on_trimmed(model, edge, end, dec)?;
        }
        if covers && convex(model, edge, &incident, dec)? {
            return Ok(true);
        }
    }
    Ok(false)
}
