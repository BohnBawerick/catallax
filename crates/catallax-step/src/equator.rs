//! Export-only cuts of torus bands across their outer equator. A written
//! torus starts its tube angle on that equator, so a band across it has no
//! regular chart. Cut there, it leaves two bands that each start on the
//! equator: the upper one on the torus itself, and the lower one on the same
//! torus written with its axis reversed, whose tube angle turns the other way.
//! Under the seam contract a band that ends on the equator from below, such as
//! the lower round of an offset cylinder, is written on that reversed torus
//! whole: its tube angle then starts on the equator like the upper band's.

use std::collections::BTreeMap;

use catallax_numbers::{Decider, Fact, FactKind, FactOrigin, Real, Refusal, Sign};
use catallax_shapes::{Circle, Curve, GeomKey, Surface, Torus};
use catallax_topology::{BodyId, Coedge, Edge, Face, FaceId, Loop, Model, Sense, Shell};

use crate::{band, trim};

pub(crate) struct Split {
    pub model: Model,
    pub bodies: Vec<BodyId>,
    /// Source face to its upper half.
    pub faces: BTreeMap<FaceId, FaceId>,
    /// Export-only equator edges added.
    pub rings: usize,
}

/// How a torus band meets its outer equator.
enum Equator {
    /// The band runs up the tube from its low ring, at this loop index,
    /// across the equator to its other ring.
    Across(usize),
    /// The band's high ring is on the equator itself, which ends the band.
    Ends,
}

// How a torus band bounded by two parallels meets the outer equator, or
// `None` when it meets it nowhere.
fn meeting(model: &Model, face: FaceId, dec: &mut Decider) -> Result<Option<Equator>, Refusal> {
    if !matches!(
        model.geom().surface(model.face(face).surface),
        Surface::Torus(_)
    ) {
        return Ok(None);
    }
    let Some((_, rings, low)) = band::tube_rings(model, face, dec)? else {
        return Ok(None);
    };
    let [first, second] = [&rings[0], &rings[1]].map(|ring| ring.tube.as_ref().unwrap());
    // The high ring on the equator itself ends the band there.
    if trim::sign(&second.height, dec)? == Sign::Zero
        && trim::sign(&second.across, dec)? == Sign::Positive
    {
        return Ok(Some(Equator::Ends));
    }
    Ok(match band::tube_below(second, first, dec)? {
        Some(true) => Some(Equator::Across(low)),
        _ => None,
    })
}

/// Cuts every torus band of `bodies` that runs across its outer equator, or
/// `None` when no band does. With `ending`, a band that ends on the equator
/// is moved onto the reversed torus, uncut.
pub(crate) fn split(
    model: &Model,
    bodies: &[BodyId],
    ending: bool,
    dec: &mut Decider,
) -> Result<Option<Split>, Refusal> {
    let mut crossing = Vec::new();
    let mut ends = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    for &body in bodies {
        for face in model.faces_of_body(body) {
            if !seen.insert(face) {
                continue;
            }
            match meeting(model, face, dec)? {
                Some(Equator::Across(low)) => crossing.push((face, low)),
                Some(Equator::Ends) if ending => ends.push(face),
                _ => {}
            }
        }
    }
    if crossing.is_empty() && ends.is_empty() {
        return Ok(None);
    }
    let mut split = model.clone();
    let mut uppers = BTreeMap::new();
    let mut lowers = BTreeMap::new();
    let count = crossing.len();
    for (face, low) in crossing {
        let source = model.face(face).clone();
        let Surface::Torus(torus) = model.geom().surface(source.surface).clone() else {
            unreachable!("an equator cut on a non-toroidal face")
        };
        let reversed = split.geom_mut().add_surface(Surface::Torus(Torus {
            axis: torus.axis.scaled(&Real::from_i64(-1)),
            ..torus.clone()
        }));
        let equator = split.geom_mut().add_curve(Curve::Circle(Circle {
            center: torus.center.clone(),
            normal: torus.axis.clone(),
            radius_sq: (&torus.major + &torus.minor).square(),
            ref_dir: torus.ref_dir.clone(),
        }));
        for surface in [source.surface, reversed] {
            split.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a: GeomKey::Curve(equator),
                b: GeomKey::Surface(surface),
                at: None,
                origin: FactOrigin::Construction("step.equator-cut"),
            });
        }
        let low_ring = source.loops[low];
        if let Loop::Coedges(coedges) = model.loop_(low_ring) {
            for coedge in coedges {
                let curve = model.edge(model.coedge(*coedge).edge).curve;
                split.geom_mut().record_fact(Fact {
                    kind: FactKind::Incident,
                    a: GeomKey::Curve(curve),
                    b: GeomKey::Surface(reversed),
                    at: None,
                    origin: FactOrigin::Construction("step.equator-cut"),
                });
            }
        }
        let edge = split.add_edge(Edge {
            curve: equator,
            sense: Sense::Same,
            start: None,
            end: None,
        });
        // The equator turns counterclockwise about the axis, and the outward
        // normal there points away from it, so a face in the torus's sense
        // keeps its upper half to the left of the equator's own direction.
        let upper = source.sense;
        let lower = match upper {
            Sense::Same => Sense::Reversed,
            Sense::Reversed => Sense::Same,
        };
        let mut cut = |sense| {
            let coedge = split.add_coedge(Coedge { edge, sense });
            split.add_loop(Loop::Coedges(vec![coedge]))
        };
        let (upper_cut, lower_cut) = (cut(upper), cut(lower));
        uppers.insert(
            face,
            split.add_face(Face {
                surface: source.surface,
                sense: source.sense,
                loops: vec![source.loops[1 - low], upper_cut],
            }),
        );
        lowers.insert(
            face,
            split.add_face(Face {
                surface: reversed,
                sense: source.sense,
                loops: vec![low_ring, lower_cut],
            }),
        );
    }
    for face in ends {
        let source = model.face(face).clone();
        let Surface::Torus(torus) = model.geom().surface(source.surface).clone() else {
            unreachable!("an equator reversal on a non-toroidal face")
        };
        let reversed = split.geom_mut().add_surface(Surface::Torus(Torus {
            axis: torus.axis.scaled(&Real::from_i64(-1)),
            ..torus
        }));
        for loop_id in &source.loops {
            let Loop::Coedges(coedges) = model.loop_(*loop_id) else {
                continue;
            };
            for coedge in coedges {
                split.geom_mut().record_fact(Fact {
                    kind: FactKind::Incident,
                    a: GeomKey::Curve(model.edge(model.coedge(*coedge).edge).curve),
                    b: GeomKey::Surface(reversed),
                    at: None,
                    origin: FactOrigin::Construction("step.equator-reversal"),
                });
            }
        }
        uppers.insert(
            face,
            split.add_face(Face {
                surface: reversed,
                sense: source.sense,
                loops: source.loops,
            }),
        );
    }
    let mut replacements = Vec::with_capacity(bodies.len());
    for &body in bodies {
        let mut body = model.body(body).clone();
        for id in &mut body.lumps {
            let mut lump = model.lump(*id).clone();
            for id in &mut lump.shells {
                let mut shell = Vec::new();
                for face in &model.shell(*id).faces {
                    shell.push(uppers.get(face).copied().unwrap_or(*face));
                    shell.extend(lowers.get(face).copied());
                }
                *id = split.add_shell(Shell { faces: shell });
            }
            *id = split.add_lump(lump);
        }
        replacements.push(split.add_body(body));
    }
    Ok(Some(Split {
        model: split,
        bodies: replacements,
        faces: uppers,
        rings: count,
    }))
}
