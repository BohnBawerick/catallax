//! Remove a redundant subdivision of one straight boundary segment.
use std::collections::BTreeSet;

use catallax_numbers::{Decider, Rational, Refusal};
use catallax_topology::{Coedge, EdgeId, FaceId, Model, Shell, ShellId, VertexId};

use super::ribbon::{position, straight, written_length_squared};
use super::{Prepared, coedges, cross, dot, ends, store_face, sub, zero};
use crate::{StepOptions, StepRegularization};

pub(super) struct Short {
    edge: EdgeId,
    other: EdgeId,
    lost: VertexId,
    kept: VertexId,
    faces: [FaceId; 2],
}

pub(super) fn find(
    model: &Model,
    shell: ShellId,
    options: &StepOptions,
    dec: &mut Decider,
) -> Result<Option<Short>, Refusal> {
    let floor = Rational::from_ratio(1, 10_000_000).unwrap();
    let floor_sq = &floor * &floor;
    let mut seen = BTreeSet::new();
    for &face in &model.shell(shell).faces {
        for wire in coedges(model, face).into_iter().flatten() {
            for c in wire {
                if !seen.insert(c.edge)
                    || written_length_squared(model, &c, options.significant_digits.max(1))
                        .is_none_or(|length| length > floor_sq)
                {
                    continue;
                }
                let Some([a, b]) = ends(model, &c) else {
                    continue;
                };
                for (lost, kept) in [(a, b), (b, a)] {
                    if let Some(found) = certify(model, shell, c.edge, lost, kept, dec)? {
                        return Ok(Some(found));
                    }
                }
            }
        }
    }
    Ok(None)
}

fn certify(
    model: &Model,
    shell: ShellId,
    edge: EdgeId,
    lost: VertexId,
    kept: VertexId,
    dec: &mut Decider,
) -> Result<Option<Short>, Refusal> {
    if lost == kept {
        return Ok(None);
    }
    let mut incident = BTreeSet::new();
    let mut faces = Vec::new();
    for &face in &model.shell(shell).faces {
        let Some(wires) = coedges(model, face) else {
            return Ok(None);
        };
        for wire in wires {
            let mut touches = false;
            for c in &wire {
                if ends(model, c).is_some_and(|pair| pair.contains(&lost)) {
                    incident.insert(c.edge);
                    touches = true;
                }
            }
            if touches {
                // Removing the subdivision must leave an ordinary closed wire.
                if wire.len() < 4 || wire.iter().filter(|c| c.edge == edge).count() != 1 {
                    return Ok(None);
                }
                faces.push(face);
            }
        }
    }
    incident.remove(&edge);
    let incident: Vec<_> = incident.into_iter().collect();
    let ([other], [fa, fb]) = (incident.as_slice(), faces.as_slice()) else {
        return Ok(None);
    };
    let other = *other;
    let e = model.edge(other);
    // Reusing the carrier preserves its endpoint incidence facts.
    if fa == fb || model.edge(edge).curve != e.curve || !straight(model, edge, dec)? {
        return Ok(None);
    }
    let far = if e.start == Some(lost) {
        e.end
    } else {
        e.start
    };
    let (Some(p), Some(q), Some(r)) = (
        position(model, kept),
        position(model, lost),
        far.and_then(|v| position(model, v)),
    ) else {
        return Ok(None);
    };
    let (before, after) = (sub(&q, &p), sub(&r, &q));
    if zero(&before)
        || zero(&after)
        || !zero(&cross(&before, &after))
        || dot(&before, &after) <= Rational::zero()
    {
        return Ok(None);
    }
    Ok(Some(Short {
        edge,
        other,
        lost,
        kept,
        faces: [*fa, *fb],
    }))
}

pub(super) fn apply(prepared: &mut Prepared, shell: ShellId, short: Short) -> ShellId {
    let mut edge = prepared.model.edge(short.other).clone();
    for end in [&mut edge.start, &mut edge.end] {
        if *end == Some(short.lost) {
            *end = Some(short.kept);
        }
    }
    let replacement = prepared.model.add_edge(edge);
    let error = [short.edge, short.other]
        .iter()
        .filter_map(|e| prepared.edge_errors.get(e))
        .cloned()
        .max()
        .unwrap_or_default();
    prepared.edge_errors.insert(replacement, error);
    let mut faces = Vec::new();
    for old in prepared.model.shell(shell).faces.clone() {
        if !short.faces.contains(&old) {
            faces.push(old);
            continue;
        }
        let wires = coedges(&prepared.model, old)
            .unwrap()
            .into_iter()
            .map(|wire| {
                wire.into_iter()
                    .filter(|c| c.edge != short.edge)
                    .map(|c| {
                        if c.edge == short.other {
                            Coedge {
                                edge: replacement,
                                ..c
                            }
                        } else {
                            c
                        }
                    })
                    .collect()
            })
            .collect();
        let new = store_face(&mut prepared.model, old, wires);
        if let Some(error) = prepared.face_errors.get(&old).cloned() {
            prepared.face_errors.insert(new, error);
        }
        faces.push(new);
    }
    prepared.records.push(StepRegularization {
        face: short.faces[0],
        edges: [short.edge, short.other],
        merged_faces: short.faces,
        max_deviation: Rational::zero(),
        removed_faces: 0,
        removed_edges: 1,
    });
    prepared.model.add_shell(Shell { faces })
}
