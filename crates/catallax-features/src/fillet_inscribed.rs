//! Two coaxial rings whose rolling balls meet at the inscribed radius.
//!
//! A shell of two faces of revolution about one axis, joined along two
//! circles, has a meridian section bounded by two arcs or an arc and a
//! segment. At the radius of the section's inscribed circle the two offsets
//! of the supports touch instead of crossing twice, and both rings' ball
//! centres reach the same circle. Each ring keeps the half of that tube on
//! its own side, the one-sided limit of its radius homotopy, so the two
//! halves close a whole ring torus and both input faces are consumed. The
//! halves share the two contact circles as their boundary. Nothing is fitted.
//! Concave rings round a void the same way, filling it but for the torus.
//! Any other configuration at that radius keeps the regular construction's
//! refusal.

use std::collections::BTreeSet;

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{Origin, Role, StepId};
use catallax_shapes::{Circle, Curve, InscribedTorus, RollingTorus, Surface};
use catallax_topology::{BodyId, Edge, EdgeId, FaceId, Model, Sense};

use super::{Uses, convex, exact_ring_circle};
use crate::{
    ConsumedFace,
    edge_blend::{self, Replacement},
    fillet::{FilletDefinition, ball_side},
    profile::sign,
};

/// The rounded body, its strip definitions and the consumed faces.
pub(super) type Rounded = (BodyId, Vec<FilletDefinition>, Vec<ConsumedFace>);

struct Half {
    edge: EdgeId,
    faces: [FaceId; 2],
    circle: Circle,
    same_sense: [bool; 2],
    supports: [Surface; 2],
    convex: bool,
    band: InscribedTorus,
}

fn zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(value, dec)? == Sign::Zero)
}

// The same circle as a point set: centre, radius and plane.
fn same_circle(a: &Circle, b: &Circle, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(zero(&a.center.minus(&b.center).norm_squared(), dec)?
        && zero(&(&a.radius_sq - &b.radius_sq), dec)?
        && zero(&a.normal.cross(&b.normal).norm_squared(), dec)?)
}

// The ring's half tube when the radius is its exact inscribed limit.
fn half(
    model: &Model,
    edge: EdgeId,
    uses: &Uses,
    radius: &Real,
    dec: &mut Decider,
) -> Result<Option<Half>, Refusal> {
    let input = model.edge(edge);
    if input.start.is_some() || input.end.is_some() {
        return Ok(None);
    }
    let incident = &uses[&edge];
    let [(a, _), (b, _)] = incident.as_slice() else {
        return Ok(None);
    };
    let faces = [*a, *b];
    let supports = faces.map(|f| model.geom().surface(model.face(f).surface).clone());
    let circle = match model.geom().curve(input.curve) {
        Curve::Circle(circle) => circle.clone(),
        Curve::Intersection(_) => match exact_ring_circle(model, edge, &supports, dec)? {
            Some(circle) => circle,
            None => return Ok(None),
        },
        _ => return Ok(None),
    };
    let convex = convex(model, edge, incident, dec)?;
    let same_sense = faces.map(|f| model.face(f).sense == Sense::Same);
    let sides = same_sense.map(|same| ball_side(same, convex));
    let band =
        match RollingTorus::inscribed([&supports[0], &supports[1]], sides, &circle, radius, dec) {
            Ok(Some(band)) => band,
            Ok(None) => return Ok(None),
            Err(refusal) if refusal.code.starts_with("shapes.rolling-torus") => return Ok(None),
            Err(refusal) => return Err(refusal),
        };
    Ok(Some(Half {
        edge,
        faces,
        circle,
        same_sense,
        supports,
        convex,
        band,
    }))
}

// The edges `face` uses.
fn edges_of(uses: &Uses, face: FaceId) -> BTreeSet<EdgeId> {
    uses.iter()
        .filter(|(_, incident)| incident.iter().any(|(f, _)| *f == face))
        .map(|(edge, _)| *edge)
        .collect()
}

// The way `face` runs round its ring `edge`, as a multiple of the ring
// circle's normal: +1 when it turns about that normal.
fn turning(model: &Model, uses: &Uses, edge: EdgeId, face: FaceId) -> Result<i64, Refusal> {
    let uses_edge = uses[&edge]
        .iter()
        .filter(|(f, _)| *f == face)
        .map(|(_, sense)| *sense)
        .collect::<Vec<_>>();
    let [coedge] = uses_edge[..] else {
        return Err(super::bad("an inscribed ring face uses its ring once"));
    };
    Ok(
        if (model.edge(edge).sense == Sense::Same) == (coedge == Sense::Same) {
            1
        } else {
            -1
        },
    )
}

// The two rings of one two-face shell, with the second turned to list the
// faces in the first's order, when they meet on one centre circle from
// opposite sides. `None` otherwise.
fn pair(
    first: Half,
    mut second: Half,
    uses: &Uses,
    dec: &mut Decider,
) -> Result<Option<[Half; 2]>, Refusal> {
    if second.faces == [first.faces[1], first.faces[0]] {
        second.faces.swap(0, 1);
        second.same_sense.swap(0, 1);
        second.supports.swap(0, 1);
        second.band.contacts.swap(0, 1);
    }
    // Both rings remove material or both add it.
    if second.faces != first.faces
        || first.faces[0] == first.faces[1]
        || first.convex != second.convex
    {
        return Ok(None);
    }
    let rings = BTreeSet::from([first.edge, second.edge]);
    for face in first.faces {
        if edges_of(uses, face) != rings {
            return Ok(None);
        }
    }
    let (a, b) = (&first.band, &second.band);
    if !zero(&a.torus.center.minus(&b.torus.center).norm_squared(), dec)?
        || !zero(&(&a.torus.major - &b.torus.major), dec)?
        || !zero(&a.torus.axis.cross(&b.torus.axis).norm_squared(), dec)?
        || same_circle(&a.middle, &b.middle, dec)?
    {
        return Ok(None);
    }
    for i in 0..2 {
        if !same_circle(&a.contacts[i], &b.contacts[i], dec)? {
            return Ok(None);
        }
    }
    Ok(Some([first, second]))
}

// Proves that nothing else in the body lies in the region a pair rounds:
// the lens bounded by its shell. A convex pair removes material inside the
// lens, so its shell must be the only shell of its lump. No void then lies
// inside it, and every other lump is disjoint material outside it. A concave
// pair adds material inside a void lens, so the body must have that one lump.
// Its outer shell and other voids lie outside the void, and no island lump
// can sit inside it. Either proof is structural; no clearance is sampled.
fn alone(model: &Model, body: BodyId, face: FaceId, convex: bool) -> Result<(), Refusal> {
    let lumps = &model.body(body).lumps;
    for lump in lumps {
        let shells = &model.lump(*lump).shells;
        let Some(at) = shells
            .iter()
            .position(|shell| model.shell(*shell).faces.contains(&face))
        else {
            continue;
        };
        let clear = if convex {
            at == 0 && shells.len() == 1
        } else {
            at > 0 && lumps.len() == 1
        };
        return if clear {
            Ok(())
        } else {
            Err(Refusal::new(
                "fillet.unsupported-overflow",
                if convex {
                    "a void inside the inscribed shell needs its blend region trimmed against it"
                } else {
                    "another lump beside the inscribed void needs the added region trimmed against it"
                },
            ))
        };
    }
    Err(super::bad("an inscribed shell is not in the body"))
}

/// The whole tori of selected ring pairs whose balls meet at the inscribed
/// radius, with their input faces consumed. A convex pair removes the rest of
/// its shell's material; a concave pair, round a void, fills the rest of the
/// void. `None`, leaving the regular construction to answer or refuse, unless
/// every selected edge is one of the two rings of a two-face shell whose
/// rings are both selected and reach that radius on one centre circle with
/// opposite halves. Refuses `fillet.unsupported-overflow` when another shell
/// or lump could lie in the region a pair rounds.
pub(super) fn build(
    model: &mut Model,
    body: BodyId,
    selected: &[EdgeId],
    uses: &Uses,
    radius: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<Option<Rounded>, Refusal> {
    if selected.is_empty() || !selected.len().is_multiple_of(2) {
        return Ok(None);
    }
    let mut halves = Vec::with_capacity(selected.len());
    for edge in selected {
        // An edge this construction cannot read is the regular builder's to
        // explain; only an exhausted allowance stops here.
        match half(model, *edge, uses, radius, dec) {
            Ok(Some(half)) => halves.push(half),
            Ok(None) => return Ok(None),
            Err(refusal) if refusal.code == "numbers.operation-budget" => return Err(refusal),
            Err(_) => return Ok(None),
        }
    }
    let mut pairs = Vec::with_capacity(halves.len() / 2);
    while let Some(first) = halves.pop() {
        let faces = BTreeSet::from(first.faces);
        let Some(at) = halves
            .iter()
            .position(|half| BTreeSet::from(half.faces) == faces)
        else {
            return Ok(None);
        };
        let second = halves.swap_remove(at);
        let Some(pair) = pair(first, second, uses, dec)? else {
            return Ok(None);
        };
        alone(model, body, pair[0].faces[0], pair[0].convex)?;
        pairs.push(pair);
    }
    let mut edit = Replacement::default();
    let mut definitions = Vec::with_capacity(selected.len());
    let mut consumed: Vec<ConsumedFace> = Vec::new();
    for halves in &pairs {
        let band = &halves[0].band;
        let torus = Surface::Torus(band.torus.clone());
        let rails = band.contacts.clone().map(|contact| {
            let curve = model.geom_mut().add_curve(Curve::Circle(contact));
            model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: None,
                end: None,
            })
        });
        let mut by = Vec::with_capacity(2);
        for half in halves {
            let name = model.edge_name(body, half.edge).ok_or_else(|| {
                super::bad("a selected edge needs an unambiguous origin-pair name")
            })?;
            let origin = Origin {
                step: step.clone(),
                role: Role::Blend { edge: name },
            };
            let mut paths = Vec::with_capacity(2);
            for (i, face) in half.faces.iter().enumerate() {
                // The band lies across each rail from where the consumed face
                // would lie, so it runs the rail the other way round.
                let along = half
                    .circle
                    .normal
                    .scaled(&Real::from_i64(-turning(model, uses, half.edge, *face)?));
                let sense = match sign(&band.contacts[i].normal.dot(&along), dec)? {
                    Sign::Positive => Sense::Same,
                    Sign::Negative => Sense::Reversed,
                    Sign::Zero => {
                        return Err(super::bad("an inscribed rail leaves the ring's plane"));
                    }
                };
                paths.push(vec![(rails[i], sense)]);
            }
            // The torus normal points away from the tube's centre circle:
            // out of a convex pair's material, into a concave pair's.
            let face = crate::extrude::face(
                model,
                torus.clone(),
                if half.convex {
                    Sense::Same
                } else {
                    Sense::Reversed
                },
                paths,
                None,
                vec![origin.clone()],
            );
            edit.faces.push((half.faces[0], face));
            by.push(origin.clone());
            definitions.push(FilletDefinition {
                origin,
                supports: half.supports.clone(),
                same_sense: half.same_sense,
                edge: Curve::Circle(half.circle.clone()),
                radius: radius.clone(),
                convex: half.convex,
                surface: torus.clone(),
                contacts: half.band.contacts.clone().map(Curve::Circle),
            });
        }
        for face in halves[0].faces {
            edit.consumed.insert(face);
            let mut source = model.face_origins(face).to_vec();
            source.sort();
            source.dedup();
            consumed.push(ConsumedFace {
                step: step.clone(),
                source,
                by: by.clone(),
            });
        }
    }
    let result = edge_blend::apply(model, body, &edit)?;
    // One record per source origin set, by every blend that removed one.
    consumed.sort_by(|a, b| a.source.cmp(&b.source));
    let mut history: Vec<ConsumedFace> = Vec::with_capacity(consumed.len());
    for record in consumed {
        match history.last_mut() {
            Some(last) if last.source == record.source => last.by.extend(record.by),
            _ => history.push(record),
        }
    }
    for record in &mut history {
        record.by.sort();
        record.by.dedup();
    }
    Ok(Some((result, definitions, history)))
}
