//! Export-only removal of a planar wedge too narrow for a STEP reader to keep.
//! A planar face's wire turns back across one straight edge at or below the
//! 1e-7 mm floor, between two long straight edges that border walls facing the
//! same way. The smaller wall moves onto the other's plane, closing that end,
//! and the short edge and the two long edges leave the file. One straight edge
//! joins their far ends unless those are already one vertex. This is the
//! one-ended case of the `ribbon` module's step, which closes two ends.

use std::collections::BTreeSet;

use catallax_numbers::{Decider, Rational, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{Coedge, Edge, EdgeId, FaceId, Model, Sense, Shell, ShellId, VertexId};

use super::ribbon::{
    Move, line, move_edges, move_wall, plane, position, rewrite, straight, written_length_squared,
};
use super::{Prepared, after, coedges, cross, dot, ends, neighbour, store_face, sub, zero};
use crate::{StepOptions, StepRegularization};

pub(super) struct Wedge {
    face: FaceId,
    wall: FaceId,
    target: FaceId,
    /// The face's three edges in wire order: long, short, long.
    wire: [Coedge; 3],
    /// Whether the first edge borders the moved wall.
    moved_first: bool,
    /// The far ends of the two long edges, in wire order, when they differ.
    rest: Option<[VertexId; 2]>,
    moving: Move,
}

impl Wedge {
    fn moved(&self) -> &Coedge {
        &self.wire[if self.moved_first { 0 } else { 2 }]
    }
    fn kept(&self) -> &Coedge {
        &self.wire[if self.moved_first { 2 } else { 0 }]
    }
}

/// Wires and positions in the wire at which a straight short edge sits between
/// two edges at least ten times longer. Reading written lengths spends no
/// decision work, so a file without one costs and writes exactly what it did.
fn short_turns(model: &Model, face: FaceId, digits: u32) -> Vec<(usize, usize)> {
    if !matches!(
        model.geom().surface(model.face(face).surface),
        Surface::Plane(_)
    ) {
        return Vec::new();
    }
    let Some(wires) = coedges(model, face) else {
        return Vec::new();
    };
    let floor = Rational::from_ratio(1, 10_000_000).unwrap();
    let floor = &floor * &floor;
    let hundred = Rational::from_i64(100);
    let mut found = Vec::new();
    for (w, wire) in wires.iter().enumerate() {
        let n = wire.len();
        if n < 3 {
            continue;
        }
        for at in 0..n {
            let Some(short) = written_length_squared(model, &wire[at], digits) else {
                continue;
            };
            if short > floor {
                continue;
            }
            let sides = [(at + n - 1) % n, (at + 1) % n]
                .map(|i| written_length_squared(model, &wire[i], digits));
            let [Some(a), Some(b)] = sides else {
                continue;
            };
            let long = a.min(b);
            if long > floor && &short * &hundred <= long {
                found.push((w, at));
            }
        }
    }
    found
}

/// Certifies moving the wall beyond `first` (or `last`) onto the plane of the
/// wall beyond the other long edge. Every check is exact: a failure keeps the
/// original topology and its original bytes.
#[allow(clippy::too_many_arguments)]
fn certify(
    model: &Model,
    shell: ShellId,
    face: FaceId,
    wire: [&Coedge; 3],
    moved_first: bool,
    limit: &Rational,
    dec: &mut Decider,
) -> Result<Option<Wedge>, Refusal> {
    // A vertex loop has no wire to rewrite, so such a shell keeps its topology.
    if model
        .shell(shell)
        .faces
        .iter()
        .any(|f| coedges(model, *f).is_none())
    {
        return Ok(None);
    }
    let [first, short, last] = wire;
    let neighbours = wire.map(|c| neighbour(model, shell, face, c.edge));
    let [Some(before), Some(side), Some(beyond)] = neighbours else {
        return Ok(None);
    };
    let (wall, target) = if moved_first {
        (before, beyond)
    } else {
        (beyond, before)
    };
    if wall == target || side == wall || side == target {
        return Ok(None);
    }
    let [Some([a, p]), Some([s0, s1]), Some([q, b])] = wire.map(|c| ends(model, c)) else {
        return Ok(None);
    };
    if p != s0 || s1 != q {
        return Ok(None);
    }
    let (Some(target_plane), Some(wall_plane)) = (plane(model, target), plane(model, wall)) else {
        return Ok(None);
    };
    let sides = dot(&target_plane.1, &wall_plane.1);
    if sides.is_zero()
        || (sides > Rational::zero()) != (model.face(target).sense == model.face(wall).sense)
    {
        return Ok(None);
    }
    for c in wire {
        if !straight(model, c.edge, dec)? {
            return Ok(None);
        }
    }
    // The moved wall's vertex on the short edge lands on its other end.
    let (moved, closed) = if moved_first { (first, q) } else { (last, p) };
    let onto = [(short.edge, closed)];
    let Some(moving) = move_wall(model, shell, wall, target, moved.edge, &onto, limit, dec)? else {
        return Ok(None);
    };
    // The two long edges must now fold back along one line from the closed end.
    let moved_to = |v: VertexId| moving.after.get(&v).cloned().or_else(|| position(model, v));
    let (Some(at), Some(far_a), Some(far_b)) = (position(model, closed), moved_to(a), moved_to(b))
    else {
        return Ok(None);
    };
    let (da, db) = (sub(&far_a, &at), sub(&far_b, &at));
    if zero(&da) || zero(&db) || !zero(&cross(&da, &db)) || dot(&da, &db) <= Rational::zero() {
        return Ok(None);
    }
    let rest = if a == b {
        None
    } else if far_a == far_b {
        return Ok(None);
    } else {
        Some([a, b])
    };
    Ok(Some(Wedge {
        face,
        wall,
        target,
        wire: [first.clone(), short.clone(), last.clone()],
        moved_first,
        rest,
        moving,
    }))
}

/// The first wedge that can close, moving the smaller-area wall of its two.
pub(super) fn find(
    model: &Model,
    shell: ShellId,
    options: &StepOptions,
    limit: &Rational,
    dec: &mut Decider,
) -> Result<Option<Wedge>, Refusal> {
    let digits = options.significant_digits.max(1);
    for &face in &model.shell(shell).faces {
        for (w, at) in short_turns(model, face, digits) {
            let wire = &coedges(model, face).unwrap()[w];
            let n = wire.len();
            let turn = [&wire[(at + n - 1) % n], &wire[at], &wire[(at + 1) % n]];
            let mut best: Option<Wedge> = None;
            for moved_first in [true, false] {
                if let Some(found) = certify(model, shell, face, turn, moved_first, limit, dec)?
                    && best
                        .as_ref()
                        .is_none_or(|b| found.moving.area < b.moving.area)
                {
                    best = Some(found);
                }
            }
            if best.is_some() {
                return Ok(best);
            }
        }
    }
    Ok(None)
}

pub(super) fn apply(
    prepared: &mut Prepared,
    shell: ShellId,
    wedge: Wedge,
    dec: &mut Decider,
) -> Result<ShellId, Refusal> {
    let removed: BTreeSet<EdgeId> = wedge.wire.iter().map(|c| c.edge).collect();
    let (vertices, edges) = move_edges(prepared, shell, wedge.wall, &wedge.moving, &removed, dec)?;
    let bound = wedge.moving.bound.clone();
    let model = &mut prepared.model;
    let rest = match wedge.rest {
        Some(ends) => {
            let [start, end] = ends.map(|v| vertices.get(&v).copied().unwrap_or(v));
            let [a, b] = [start, end].map(|v| position(model, v).unwrap());
            let curve = line(model, &a, &b, dec)?;
            let edge = model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: Some(start),
                end: Some(end),
            });
            let error = [&wedge.wire[0], &wedge.wire[2]]
                .map(|c| {
                    prepared
                        .edge_errors
                        .get(&c.edge)
                        .cloned()
                        .unwrap_or_default()
                })
                .into_iter()
                .fold(Rational::zero(), Rational::max);
            prepared.edge_errors.insert(edge, error + bound.clone());
            Some(edge)
        }
        None => None,
    };
    let model = &mut prepared.model;
    let forward = rest.map(|edge| Coedge {
        edge,
        sense: Sense::Same,
    });
    let backward = rest.map(|edge| Coedge {
        edge,
        sense: Sense::Reversed,
    });
    let wall_wire = coedges(model, wedge.wall).unwrap().remove(0);
    let old_faces = model.shell(shell).faces.clone();
    let mut faces = Vec::with_capacity(old_faces.len());
    let mut removed_faces = 1;
    for &face in &old_faces {
        if face == wedge.wall {
            continue;
        }
        let mut wires = coedges(model, face).unwrap();
        if face == wedge.face {
            // The turn closes: its far ends join, or meet already.
            let short = wedge.wire[1].edge;
            for wire in &mut wires {
                if let Some(at) = wire.iter().position(|c| c.edge == short) {
                    let n = wire.len();
                    wire.rotate_left((at + n - 1) % n);
                    wire.splice(0..3, forward.clone());
                    break;
                }
            }
            wires.retain(|wire| !wire.is_empty());
            if wires.is_empty() {
                removed_faces += 1;
                continue;
            }
        } else if face == wedge.target {
            // The kept wall takes the moved wall's other edges, and the edge
            // joining the far ends where the kept edge reached past them.
            let mut path = after(&wall_wire, wedge.moved().edge);
            if wedge.moved_first {
                path.splice(0..0, backward.clone());
            } else {
                path.extend(backward.clone());
            }
            for wire in &mut wires {
                if let Some(at) = wire.iter().position(|c| c.edge == wedge.kept().edge) {
                    wire.splice(at..=at, path.iter().cloned());
                    break;
                }
            }
        }
        let new_wires = rewrite(wires.clone(), &removed, &edges);
        if new_wires == wires && face != wedge.face {
            faces.push(face);
            continue;
        }
        let stored = store_face(model, face, new_wires);
        let own = prepared.face_errors.get(&face).cloned().unwrap_or_default();
        let wall = (face == wedge.target)
            .then(|| prepared.face_errors.get(&wedge.wall).cloned())
            .flatten()
            .unwrap_or_default();
        prepared
            .face_errors
            .insert(stored, own.max(wall) + bound.clone());
        faces.push(stored);
    }
    prepared.records.push(StepRegularization {
        face: wedge.face,
        edges: [wedge.moved().edge, wedge.kept().edge],
        merged_faces: [wedge.wall, wedge.target],
        max_deviation: bound,
        removed_faces,
        removed_edges: removed.len() - usize::from(rest.is_some()),
    });
    Ok(model.add_shell(Shell { faces }))
}
