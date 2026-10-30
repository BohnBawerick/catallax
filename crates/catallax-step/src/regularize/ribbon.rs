//! Export-only removal of a planar step too narrow for a STEP reader to keep.
//! OpenCASCADE cannot translate an edge of at most 1e-7 mm, so a ribbon face
//! with two such end edges is unreadable whatever the file declares. The wall
//! on one side of the step moves onto the other wall's plane, and the ribbon's
//! two end edges, its two long edges and the moved wall all leave the file.

use std::collections::{BTreeMap, BTreeSet};

use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::{Curve, CurveId, Line, Point3, PointDef, Surface, Vector3};
use catallax_topology::{
    Coedge, Edge, EdgeId, FaceId, Model, Sense, Shell, ShellId, Vertex, VertexId,
};

use super::{P, Prepared, Wires, after, coedges, cross, dot, ends, neighbour, point, positions};
use super::{rounded, store_face, sub};
use crate::{StepOptions, StepRegularization};

/// Faces one collapse takes out of the file: the ribbon and the moved wall.
const REMOVED_FACES: usize = 2;
/// Edges one collapse takes out of the file: the ribbon's whole wire.
const REMOVED_EDGES: usize = 4;

/// Where a vertex of the moved wall goes.
pub(super) enum Slide {
    /// Onto the kept end of a short edge that closes.
    Onto(VertexId),
    /// Along its one edge outside the wall, to the kept wall's plane.
    To(P),
}

/// A certified move of one planar wall onto another wall's plane.
pub(super) struct Move {
    /// Every wall vertex that moves. Any other stays where it is.
    pub(super) slides: BTreeMap<VertexId, Slide>,
    /// Every wall vertex's position after the move.
    pub(super) after: BTreeMap<VertexId, P>,
    /// Largest L1 vertex displacement.
    pub(super) bound: Rational,
    /// Squared doubled area of the wall, to move the smaller of two walls.
    pub(super) area: Rational,
}

pub(super) struct Ribbon {
    face: FaceId,
    wall: FaceId,
    target: FaceId,
    moved: Coedge,
    kept: Coedge,
    end_edges: [EdgeId; 2],
    moving: Move,
}

pub(super) fn plane(model: &Model, face: FaceId) -> Option<(P, P)> {
    let Surface::Plane(s) = model.geom().surface(model.face(face).surface) else {
        return None;
    };
    let normal = Point3::new(s.normal.x.clone(), s.normal.y.clone(), s.normal.z.clone());
    Some((point(&s.origin)?, point(&normal)?))
}

pub(super) fn written_length_squared(model: &Model, c: &Coedge, digits: u32) -> Option<Rational> {
    let [a, b] = positions(model, c)?.map(|p| rounded(&p, digits));
    let d = sub(&b, &a);
    Some(dot(&d, &d))
}

pub(super) fn straight(model: &Model, edge: EdgeId, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(matches!(
        curve_form(model.geom(), model.edge(edge).curve, dec)?,
        CurveForm::Line { .. }
    ))
}

pub(super) fn position(model: &Model, vertex: VertexId) -> Option<P> {
    point(&model.geom().position(model.vertex(vertex).point))
}

fn on_plane(p: &P, plane: &(P, P)) -> bool {
    dot(&sub(p, &plane.0), &plane.1).is_zero()
}

fn distance_l1(a: &P, b: &P) -> Rational {
    sub(a, b).iter().fold(Rational::zero(), |s, x| s + x.abs())
}

fn vertex_edges(model: &Model, shell: ShellId) -> BTreeMap<VertexId, BTreeSet<EdgeId>> {
    let mut found: BTreeMap<VertexId, BTreeSet<EdgeId>> = BTreeMap::new();
    for &face in &model.shell(shell).faces {
        for wire in coedges(model, face).into_iter().flatten() {
            for c in wire {
                let edge = model.edge(c.edge);
                for vertex in [edge.start, edge.end].into_iter().flatten() {
                    found.entry(vertex).or_default().insert(c.edge);
                }
            }
        }
    }
    found
}

/// The ribbon's wire, rotated to start at an end edge, when opposite edges
/// are short and long. Reading written lengths spends no decision work, so a
/// file without such a face costs and writes exactly what it did.
fn short_ended_wire(model: &Model, face: FaceId, digits: u32) -> Option<[Coedge; 4]> {
    let f = model.face(face);
    if !matches!(model.geom().surface(f.surface), Surface::Plane(_)) {
        return None;
    }
    let wires = coedges(model, face)?;
    let [wire] = wires.as_slice() else {
        return None;
    };
    let wire: [Coedge; 4] = wire.clone().try_into().ok()?;
    let floor = Rational::from_ratio(1, 10_000_000).unwrap();
    let floor = &floor * &floor;
    let lengths = wire
        .each_ref()
        .map(|c| written_length_squared(model, c, digits));
    let lengths: [Rational; 4] = lengths
        .into_iter()
        .collect::<Option<Vec<_>>>()?
        .try_into()
        .ok()?;
    // The width must also be under a tenth of the length, the feature-size
    // rule of the tolerance policy. A uniformly tiny body has no such face.
    let hundred = Rational::from_i64(100);
    for at in 0..2 {
        let (e0, e1) = (&lengths[at], &lengths[at + 2]);
        let (l0, l1) = (&lengths[at + 1], &lengths[(at + 3) % 4]);
        let long = l0.clone().min(l1.clone());
        if *e0 <= floor
            && *e1 <= floor
            && long > floor
            && &e0.clone().max(e1.clone()) * &hundred <= long
        {
            let mut wire = wire;
            wire.rotate_left(at);
            return Some(wire);
        }
    }
    None
}

/// Certifies moving `moved`'s wall onto `kept`'s. Every check is exact: a
/// failure keeps the original topology and its original bytes.
fn certify(
    model: &Model,
    shell: ShellId,
    face: FaceId,
    [kept, end_a, moved, end_b]: [&Coedge; 4],
    limit: &Rational,
    dec: &mut Decider,
) -> Result<Option<Ribbon>, Refusal> {
    // A vertex loop has no wire to rewrite, so such a shell keeps its topology.
    if model
        .shell(shell)
        .faces
        .iter()
        .any(|f| coedges(model, *f).is_none())
    {
        return Ok(None);
    }
    let neighbours = [kept, end_a, moved, end_b].map(|c| neighbour(model, shell, face, c.edge));
    let [Some(target), Some(side_a), Some(wall), Some(side_b)] = neighbours else {
        return Ok(None);
    };
    if wall == target || [wall, target].iter().any(|f| *f == side_a || *f == side_b) {
        return Ok(None);
    }
    let ends_of = [kept, end_a, moved, end_b].map(|c| ends(model, c));
    let [
        Some([a0, a1]),
        Some([b0, b1]),
        Some([c0, c1]),
        Some([d0, d1]),
    ] = ends_of
    else {
        return Ok(None);
    };
    if a1 != b0 || b1 != c0 || c1 != d0 || d1 != a0 {
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
    for c in [kept, end_a, moved, end_b] {
        if !straight(model, c.edge, dec)? {
            return Ok(None);
        }
    }
    let onto = [(end_a.edge, b0), (end_b.edge, d1)];
    let Some(moving) = move_wall(model, shell, wall, target, moved.edge, &onto, limit, dec)? else {
        return Ok(None);
    };
    Ok(Some(Ribbon {
        face,
        wall,
        target,
        moved: moved.clone(),
        kept: kept.clone(),
        end_edges: [end_a.edge, end_b.edge],
        moving,
    }))
}

/// Certifies moving `wall` onto `target`'s plane, every vertex either onto
/// the vertex `onto` names for its outside edge or along that edge. Every
/// check is exact: a failure keeps the original topology and its original bytes.
#[allow(clippy::too_many_arguments)]
pub(super) fn move_wall(
    model: &Model,
    shell: ShellId,
    wall: FaceId,
    target: FaceId,
    moved: EdgeId,
    onto: &[(EdgeId, VertexId)],
    limit: &Rational,
    dec: &mut Decider,
) -> Result<Option<Move>, Refusal> {
    let Some(target_plane) = plane(model, target) else {
        return Ok(None);
    };
    let wall_wires = coedges(model, wall);
    let (Some([wall_wire]), Some(_)) = (
        wall_wires
            .as_deref()
            .and_then(|w| <&[_; 1]>::try_from(w).ok()),
        coedges(model, target),
    ) else {
        return Ok(None);
    };
    let wall_edges: BTreeSet<EdgeId> = wall_wire.iter().map(|c| c.edge).collect();
    let mut wall_vertices = BTreeSet::new();
    for c in wall_wire {
        let Some(pair) = ends(model, c) else {
            return Ok(None);
        };
        wall_vertices.extend(pair);
        if !straight(model, c.edge, dec)? {
            return Ok(None);
        }
    }
    let at_vertex = vertex_edges(model, shell);
    let mut slides = BTreeMap::new();
    let mut after_move = BTreeMap::new();
    let mut bound = Rational::zero();
    for &vertex in &wall_vertices {
        let others: Vec<EdgeId> = at_vertex[&vertex]
            .iter()
            .copied()
            .filter(|e| !wall_edges.contains(e))
            .collect();
        let Some(from) = position(model, vertex) else {
            return Ok(None);
        };
        let [other] = others.as_slice() else {
            // A vertex where more faces meet can only stay where it is.
            if on_plane(&from, &target_plane) {
                after_move.insert(vertex, from);
                continue;
            }
            return Ok(None);
        };
        let onto = onto.iter().find(|(edge, _)| edge == other).map(|(_, v)| *v);
        let (slide, to) = if let Some(onto) = onto {
            let Some(to) = position(model, onto) else {
                return Ok(None);
            };
            (Slide::Onto(onto), to)
        } else {
            let (Some(far), true) = (
                far_end(model, *other, vertex),
                straight(model, *other, dec)?,
            ) else {
                return Ok(None);
            };
            let Some(base) = position(model, far) else {
                return Ok(None);
            };
            if wall_vertices.contains(&far) {
                return Ok(None);
            }
            let along = sub(&from, &base);
            let across = dot(&along, &target_plane.1);
            if across.is_zero() {
                return Ok(None);
            }
            let t = dot(&sub(&target_plane.0, &base), &target_plane.1) / across;
            if t <= Rational::zero() {
                return Ok(None);
            }
            let to: P = std::array::from_fn(|i| &base[i] + &(&t * &along[i]));
            (Slide::To(to.clone()), to)
        };
        if !on_plane(&to, &target_plane) {
            return Ok(None);
        }
        bound = bound.max(distance_l1(&from, &to));
        slides.insert(vertex, slide);
        after_move.insert(vertex, to);
    }
    // A planar wall moved vertex by vertex maps piecewise linearly, so no point
    // of it, the ribbon or a neighbour's trim moves further than a vertex does.
    if bound > *limit {
        return Ok(None);
    }
    // Each new wall edge must still lie on the face beyond it.
    let mut sum = [Rational::zero(), Rational::zero(), Rational::zero()];
    for c in wall_wire {
        let [a, b] = ends(model, c).unwrap();
        let (to_a, to_b) = (&after_move[&a], &after_move[&b]);
        let part = cross(&position(model, a).unwrap(), &position(model, b).unwrap());
        sum = std::array::from_fn(|i| &sum[i] + &part[i]);
        if c.edge == moved {
            continue;
        }
        let Some(face) = neighbour(model, shell, wall, c.edge) else {
            return Ok(None);
        };
        let Some(beyond) = plane(model, face).filter(|_| face != target) else {
            return Ok(None);
        };
        // A vertex landing on one the face beyond already has stays as true to
        // that face as it was.
        let lies = |v: VertexId, to: &P| {
            on_plane(to, &beyond)
                || matches!(slides.get(&v), Some(Slide::Onto(o)) if has_vertex(model, face, *o))
        };
        if to_a == to_b || !lies(a, to_a) || !lies(b, to_b) {
            return Ok(None);
        }
    }
    Ok(Some(Move {
        slides,
        after: after_move,
        bound,
        area: dot(&sum, &sum),
    }))
}

fn has_vertex(model: &Model, face: FaceId, vertex: VertexId) -> bool {
    coedges(model, face)
        .into_iter()
        .flatten()
        .flatten()
        .any(|c| ends(model, &c).is_some_and(|pair| pair.contains(&vertex)))
}

fn far_end(model: &Model, edge: EdgeId, vertex: VertexId) -> Option<VertexId> {
    let e = model.edge(edge);
    let (start, end) = (e.start?, e.end?);
    (start == vertex)
        .then_some(end)
        .or((end == vertex).then_some(start))
}

/// The smallest-area wall that can move, so the file's volume changes least.
pub(super) fn find(
    model: &Model,
    shell: ShellId,
    options: &StepOptions,
    limit: &Rational,
    dec: &mut Decider,
) -> Result<Option<Ribbon>, Refusal> {
    let digits = options.significant_digits.max(1);
    for &face in &model.shell(shell).faces {
        let Some([end_a, long_a, end_b, long_b]) = short_ended_wire(model, face, digits) else {
            continue;
        };
        let mut best: Option<Ribbon> = None;
        for wire in [
            [&long_a, &end_b, &long_b, &end_a],
            [&long_b, &end_a, &long_a, &end_b],
        ] {
            if let Some(found) = certify(model, shell, face, wire, limit, dec)?
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
    Ok(None)
}

fn rational_point(p: &P) -> Point3 {
    Point3::new(
        Real::from(p[0].clone()),
        Real::from(p[1].clone()),
        Real::from(p[2].clone()),
    )
}

/// Old-to-new vertices and edges of one move.
pub(super) type Copies = (BTreeMap<VertexId, VertexId>, BTreeMap<EdgeId, EdgeId>);

/// Creates the moved wall's vertices and edges, and every edge whose ends
/// moved, carrying `moving.bound` into their errors. Returns the old-to-new
/// vertex and edge maps; `removed` edges leave the file and get no copy.
pub(super) fn move_edges(
    prepared: &mut Prepared,
    shell: ShellId,
    wall: FaceId,
    moving: &Move,
    removed: &BTreeSet<EdgeId>,
    dec: &mut Decider,
) -> Result<Copies, Refusal> {
    let model = &mut prepared.model;
    let vertices: BTreeMap<VertexId, VertexId> = moving
        .slides
        .iter()
        .map(|(old, slide)| {
            let new = match slide {
                Slide::Onto(v) => *v,
                Slide::To(p) => {
                    let point = model
                        .geom_mut()
                        .add_point(PointDef::Explicit(rational_point(p)));
                    model.add_vertex(Vertex { point })
                }
            };
            (*old, new)
        })
        .collect();
    let mapped = |v: Option<VertexId>| v.map(|v| vertices.get(&v).copied().unwrap_or(v));
    let wall_wire = coedges(model, wall).unwrap().remove(0);
    let wall_edges: BTreeSet<EdgeId> = wall_wire.iter().map(|c| c.edge).collect();
    let mut edges = BTreeMap::new();
    let old_faces = model.shell(shell).faces.clone();
    for &face in &old_faces {
        for c in coedges(model, face).unwrap().into_iter().flatten() {
            if removed.contains(&c.edge) || edges.contains_key(&c.edge) {
                continue;
            }
            let source = model.edge(c.edge).clone();
            let (start, end) = (mapped(source.start), mapped(source.end));
            if wall_edges.contains(&c.edge) {
                let (a, b) = (
                    position(model, start.unwrap()),
                    position(model, end.unwrap()),
                );
                let (a, b) = (a.unwrap(), b.unwrap());
                let curve = line(model, &a, &b, dec)?;
                let edge = model.add_edge(Edge {
                    curve,
                    sense: Sense::Same,
                    start,
                    end,
                });
                edges.insert(c.edge, edge);
            } else if (start, end) != (source.start, source.end) {
                edges.insert(
                    c.edge,
                    model.add_edge(Edge {
                        start,
                        end,
                        ..source
                    }),
                );
            }
        }
    }
    for (old, new) in &edges {
        let error =
            prepared.edge_errors.get(old).cloned().unwrap_or_default() + moving.bound.clone();
        prepared.edge_errors.insert(*new, error);
    }
    Ok((vertices, edges))
}

/// Adds the straight carrier from `a` to `b`.
pub(super) fn line(model: &mut Model, a: &P, b: &P, dec: &mut Decider) -> Result<CurveId, Refusal> {
    let d = sub(b, a);
    let line = Line::new(
        rational_point(a),
        Vector3::new(
            Real::from(d[0].clone()),
            Real::from(d[1].clone()),
            Real::from(d[2].clone()),
        ),
        dec,
    )?;
    Ok(model.geom_mut().add_curve(Curve::Line(line)))
}

/// Replaces wire edges by their copies and drops the removed ones.
pub(super) fn rewrite(
    wires: Wires,
    removed: &BTreeSet<EdgeId>,
    edges: &BTreeMap<EdgeId, EdgeId>,
) -> Wires {
    wires
        .into_iter()
        .map(|wire| {
            wire.into_iter()
                .filter(|c| !removed.contains(&c.edge))
                .map(|c| Coedge {
                    edge: edges.get(&c.edge).copied().unwrap_or(c.edge),
                    sense: c.sense,
                })
                .collect()
        })
        .collect()
}

pub(super) fn apply(
    prepared: &mut Prepared,
    shell: ShellId,
    ribbon: Ribbon,
    dec: &mut Decider,
) -> Result<ShellId, Refusal> {
    let removed: BTreeSet<EdgeId> = ribbon
        .end_edges
        .into_iter()
        .chain([ribbon.moved.edge, ribbon.kept.edge])
        .collect();
    let (_, edges) = move_edges(prepared, shell, ribbon.wall, &ribbon.moving, &removed, dec)?;
    let model = &mut prepared.model;
    let wall_wire = coedges(model, ribbon.wall).unwrap().remove(0);
    let old_faces = model.shell(shell).faces.clone();
    let mut faces = Vec::with_capacity(old_faces.len() - 1);
    for &face in &old_faces {
        if face == ribbon.face || face == ribbon.wall {
            continue;
        }
        let mut wires = coedges(model, face).unwrap();
        if face == ribbon.target {
            let path = after(&wall_wire, ribbon.moved.edge);
            for wire in &mut wires {
                if let Some(at) = wire.iter().position(|c| c.edge == ribbon.kept.edge) {
                    wire.splice(at..=at, path.iter().cloned());
                    break;
                }
            }
        }
        let new_wires = rewrite(wires.clone(), &removed, &edges);
        if new_wires == wires {
            faces.push(face);
            continue;
        }
        let stored = store_face(model, face, new_wires);
        let own = prepared.face_errors.get(&face).cloned().unwrap_or_default();
        let wall = (face == ribbon.target)
            .then(|| prepared.face_errors.get(&ribbon.wall).cloned())
            .flatten()
            .unwrap_or_default();
        prepared
            .face_errors
            .insert(stored, own.max(wall) + ribbon.moving.bound.clone());
        faces.push(stored);
    }
    prepared.records.push(StepRegularization {
        face: ribbon.face,
        edges: [ribbon.moved.edge, ribbon.kept.edge],
        merged_faces: [ribbon.wall, ribbon.target],
        max_deviation: ribbon.moving.bound,
        removed_faces: REMOVED_FACES,
        removed_edges: REMOVED_EDGES,
    });
    Ok(model.add_shell(Shell { faces }))
}
