//! Sew selected face uses, split contact fans, then remove internal
//! boundaries and degree-two vertices before writing the result topology.

use crate::Containment;
use crate::geometry::{
    coincident, oriented, product, record_incident, refusal, reverse, sign, unit,
};
use crate::points::curve_sample;
use crate::stages::{Imprinted, Selection};
use catallax_intersect::{Component, SurfacePair, intersect_surfaces, same_curve, tangent_at};
use catallax_numbers::{Answer, Answered, Decider, MarkLog, Outcome, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Curve, CurveId, GeomKey, PointId, RVector3, SurfaceId};
use catallax_topology::{
    Body, BodyId, BodyKind, Coedge, Edge, Face, FaceId, Loop, Lump, Model, Sense, Shell, Vertex,
};
use std::collections::{BTreeMap, BTreeSet};

struct Groups {
    parents: Vec<usize>,
}
impl Groups {
    fn new(n: usize) -> Self {
        Self {
            parents: (0..n).collect(),
        }
    }
    fn root(&mut self, mut i: usize) -> usize {
        while self.parents[i] != i {
            self.parents[i] = self.parents[self.parents[i]];
            i = self.parents[i];
        }
        i
    }
    fn join(&mut self, a: usize, b: usize) {
        let a = self.root(a);
        let b = self.root(b);
        if a != b {
            let (lo, hi) = (a.min(b), a.max(b));
            self.parents[hi] = lo;
        }
    }
}
#[derive(Clone)]
struct FaceData {
    surface: SurfaceId,
    sense: Sense,
    original: usize,
    loops: Vec<Vec<(usize, Sense)>>,
    vertices: Vec<catallax_topology::VertexId>,
    origins: Vec<Origin>,
}
#[derive(Clone)]
struct Use {
    edge: catallax_topology::EdgeId,
    face: usize,
    sense: Sense,
    next: usize,
}
#[derive(Clone)]
struct EdgeData {
    curve: CurveId,
    sense: Sense,
    from: Option<usize>,
    to: Option<usize>,
    faces: [usize; 2],
    active: bool,
}
fn use_ends(edge: &EdgeData, sense: Sense) -> (Option<usize>, Option<usize>) {
    if sense == Sense::Same {
        (edge.from, edge.to)
    } else {
        (edge.to, edge.from)
    }
}

type FaceUses = (
    Vec<FaceData>,
    Vec<Use>,
    BTreeMap<catallax_topology::EdgeId, Vec<usize>>,
);
type MergedFaces = (Vec<FaceData>, Vec<bool>, Vec<usize>, Groups);

fn read(model: &Model, selection: &Selection) -> FaceUses {
    let mut faces = Vec::new();
    let mut uses = Vec::new();
    let mut by_edge: BTreeMap<_, Vec<_>> = BTreeMap::new();
    for &(face_id, selected_sense) in &selection.faces {
        let face = model.face(face_id);
        let face_index = faces.len();
        let mut loops = Vec::new();
        let mut vertices = Vec::new();
        for &l in &face.loops {
            let cs = match model.loop_(l) {
                Loop::Coedges(cs) => cs,
                Loop::Vertex(vertex) => {
                    vertices.push(*vertex);
                    continue;
                }
            };
            let mut coedges: Vec<_> = cs.iter().map(|&c| model.coedge(c).clone()).collect();
            if selected_sense == Sense::Reversed {
                coedges.reverse();
                for c in &mut coedges {
                    c.sense = reverse(c.sense);
                }
            }
            let start = uses.len();
            let count = coedges.len();
            let mut loop_uses = Vec::new();
            for (i, c) in coedges.into_iter().enumerate() {
                let index = uses.len();
                by_edge.entry(c.edge).or_default().push(index);
                uses.push(Use {
                    edge: c.edge,
                    face: face_index,
                    sense: c.sense,
                    next: start + (i + 1) % count,
                });
                loop_uses.push((index, Sense::Same));
            }
            loops.push(loop_uses);
        }
        let mut origins = model.face_origins(face_id).to_vec();
        if let Some(other) = selection.coincident.get(&face_id) {
            origins.extend(model.face_origins(*other).iter().cloned());
        }
        faces.push(FaceData {
            surface: face.surface,
            sense: product(face.sense, selected_sense),
            original: face_index,
            loops,
            vertices,
            origins,
        });
    }
    (faces, uses, by_edge)
}

fn pair(
    model: &mut Model,
    faces: &[FaceData],
    uses: &[Use],
    by_edge: &BTreeMap<catallax_topology::EdgeId, Vec<usize>>,
    dec: &mut Decider,
) -> Result<Vec<usize>, Refusal> {
    let mut twins = vec![usize::MAX; uses.len()];
    for (&edge_id, list) in by_edge {
        if list.len() == 2 {
            if uses[list[0]].sense == uses[list[1]].sense {
                return Err(refusal(
                    "boolean.merge-ambiguous",
                    "kept edge uses do not have opposite senses",
                ));
            }
            twins[list[0]] = list[1];
            twins[list[1]] = list[0];
            continue;
        }
        if list.len() < 2 || list.len() % 2 != 0 {
            return Err(refusal(
                "boolean.merge-ambiguous",
                "kept boundary does not close in edge pairs",
            ));
        }
        let edge = model.edge(edge_id).clone();
        let (a, b) = crate::geometry::natural_ends(model, edge_id);
        let point = curve_sample(model, edge.curve, a, b, dec)?;
        let tangent = oriented(
            &tangent_at(model.geom(), edge.curve, point, dec)?[1],
            edge.sense,
        );
        let mut germs = Vec::new();
        for &index in list {
            let face = &faces[uses[index].face];
            let normal = oriented(&model.geom().normal(face.surface, point), face.sense);
            let inward = normal.cross(&oriented(&tangent, uses[index].sense));
            germs.push(crate::germ::surface_germ(
                model,
                face.surface,
                point,
                inward,
                index,
                dec,
            )?);
        }
        crate::germ::sort(model, &mut germs, &tangent, dec, "boolean.merge-ambiguous")?;
        for (i, germ) in germs.iter().enumerate() {
            let index = germ.key;
            if uses[index].sense != Sense::Same {
                continue;
            }
            let other = germs[(i + germs.len() - 1) % germs.len()].key;
            if uses[other].sense != Sense::Reversed || twins[other] != usize::MAX {
                return Err(refusal(
                    "boolean.self-touching-shell",
                    "radial kept faces do not bound disjoint material wedges",
                ));
            }
            twins[index] = other;
            twins[other] = index;
        }
    }
    if twins.contains(&usize::MAX) {
        return Err(refusal(
            "boolean.merge-ambiguous",
            "a kept coedge has no sewn partner",
        ));
    }
    Ok(twins)
}

fn merge_faces(
    model: &mut Model,
    faces: &[FaceData],
    uses: &[Use],
    twins: &[usize],
    dec: &mut Decider,
) -> Result<MergedFaces, Refusal> {
    let mut components = Groups::new(faces.len());
    let mut merged = Groups::new(faces.len());
    for (i, u) in uses.iter().enumerate() {
        components.join(u.face, uses[twins[i]].face);
    }
    for (i, u) in uses.iter().enumerate() {
        let other = twins[i];
        if i > other {
            continue;
        }
        let a = &faces[u.face];
        let b = &faces[uses[other].face];
        if a.original == b.original {
            continue;
        }
        if coincident(
            model.geom().same_surface(a.surface, b.surface, dec),
            "boolean.merge-ambiguous",
        )? {
            let edge = model.edge(u.edge).clone();
            let (from, to) = crate::geometry::natural_ends(model, u.edge);
            let point = curve_sample(model, edge.curve, from, to, dec)?;
            let agree = match sign(
                &model
                    .geom()
                    .normal(a.surface, point)
                    .dot(&model.geom().normal(b.surface, point)),
                dec,
                "boolean.merge-ambiguous",
            )? {
                Sign::Positive => true,
                Sign::Negative => false,
                Sign::Zero => {
                    return Err(refusal(
                        "boolean.merge-ambiguous",
                        "coincident surface normals vanished",
                    ));
                }
            };
            if agree == (a.sense == b.sense) {
                merged.join(u.face, uses[other].face);
            }
        }
    }
    let removed: Vec<_> = uses
        .iter()
        .enumerate()
        .map(|(i, u)| merged.root(u.face) == merged.root(uses[twins[i]].face))
        .collect();
    let mut next = Vec::with_capacity(uses.len());
    for u in uses {
        let mut cursor = u.next;
        let mut count = 0;
        while removed[cursor] {
            cursor = uses[twins[cursor]].next;
            count += 1;
            if count > uses.len() {
                break;
            }
        }
        next.push(cursor);
    }
    let mut result = Vec::new();
    let mut mapping = BTreeMap::new();
    for (i, face) in faces.iter().enumerate() {
        let root = merged.root(i);
        if root == i {
            mapping.insert(root, result.len());
            result.push(FaceData {
                surface: face.surface,
                sense: face.sense,
                original: i,
                loops: Vec::new(),
                vertices: Vec::new(),
                origins: Vec::new(),
            });
        }
    }
    for (i, face) in faces.iter().enumerate() {
        let target = mapping[&merged.root(i)];
        result[target].origins.extend(face.origins.iter().cloned());
        for &vertex in &face.vertices {
            let point = model.vertex(vertex).point;
            let mut present = false;
            for &old in &result[target].vertices {
                if coincident(
                    model.geom().same_point(point, model.vertex(old).point, dec),
                    "boolean.merge-ambiguous",
                )? {
                    present = true;
                    break;
                }
            }
            if !present {
                result[target].vertices.push(vertex);
            }
        }
    }
    let mut visited = BTreeSet::new();
    for seed in 0..uses.len() {
        if removed[seed] || visited.contains(&seed) {
            continue;
        }
        let face = mapping[&merged.root(uses[seed].face)];
        let mut loop_uses = Vec::new();
        let mut current = seed;
        loop {
            if !visited.insert(current) {
                if current == seed {
                    break;
                }
                return Err(refusal(
                    "boolean.merge-ambiguous",
                    "merged face walks do not form disjoint cycles",
                ));
            }
            if removed[current] {
                return Err(refusal(
                    "boolean.merge-ambiguous",
                    "internal-edge removal left a broken face walk",
                ));
            }
            loop_uses.push((current, Sense::Same));
            current = next[current];
        }
        result[face].loops.push(loop_uses);
    }
    Ok((result, removed, next, components))
}

fn singular(
    model: &mut Model,
    edge: &EdgeData,
    point: PointId,
    faces: &[FaceData],
    nodes: &BTreeSet<PointId>,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for &node in nodes {
        if coincident(
            model.geom().same_point(node, point, dec),
            "boolean.merge-ambiguous",
        )? {
            return Ok(true);
        }
    }
    let supports = match model.geom().curve(edge.curve) {
        Curve::Intersection(c) => c.surfaces,
        _ => [faces[edge.faces[0]].surface, faces[edge.faces[1]].surface],
    };
    let pair = intersect_surfaces(model.geom_mut(), supports[0], supports[1], dec)?;
    if let SurfacePair::Meets(parts) = pair.value {
        for part in parts {
            if let Component::Curve { singular, .. } = part {
                for node in singular {
                    if coincident(
                        model.geom().same_point(node, point, dec),
                        "boolean.merge-ambiguous",
                    )? {
                        return Ok(true);
                    }
                }
            }
        }
    }
    Ok(false)
}

fn maximal_edges(
    model: &mut Model,
    faces: &mut [FaceData],
    edges: &mut Vec<EdgeData>,
    points: &BTreeMap<usize, PointId>,
    nodes: &BTreeSet<PointId>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    loop {
        let mut incident: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (i, e) in edges.iter().enumerate() {
            if e.active {
                for v in [e.from, e.to].into_iter().flatten() {
                    incident.entry(v).or_default().push(i);
                }
            }
        }
        let mut changed = false;
        for (vertex, at) in incident {
            if at.len() != 2 {
                continue;
            }
            let first = at[0];
            let second = at[1];
            if first == second {
                if !singular(model, &edges[first], points[&vertex], faces, nodes, dec)? {
                    edges[first].from = None;
                    edges[first].to = None;
                    changed = true;
                    break;
                }
                continue;
            }
            if edges[first].faces != edges[second].faces {
                continue;
            }
            if !coincident(
                same_curve(model.geom(), edges[first].curve, edges[second].curve, dec),
                "boolean.merge-ambiguous",
            )? {
                continue;
            }
            if singular(model, &edges[first], points[&vertex], faces, nodes, dec)? {
                continue;
            }
            let a = edges[first].clone();
            let b = edges[second].clone();
            let from = if a.from == Some(vertex) { a.to } else { a.from };
            let to = if b.from == Some(vertex) { b.to } else { b.from };
            let sense = if a.to == Some(vertex) {
                a.sense
            } else {
                reverse(a.sense)
            };
            let merged = edges.len();
            edges.push(EdgeData {
                curve: a.curve,
                sense,
                from,
                to,
                faces: a.faces,
                active: true,
            });
            let mut replacements = 0;
            for face in faces.iter_mut() {
                for loop_uses in &mut face.loops {
                    let count = loop_uses.len();
                    let position = (0..count).find(|&i| {
                        let (x, s) = loop_uses[i];
                        let (y, t) = loop_uses[(i + 1) % count];
                        ((x == first && y == second) || (x == second && y == first))
                            && use_ends(&edges[x], s).1 == Some(vertex)
                            && use_ends(&edges[y], t).0 == Some(vertex)
                    });
                    if let Some(i) = position {
                        loop_uses.rotate_left(i);
                        // The merged direction follows `first` into this
                        // vertex. Endpoints cannot distinguish the two walks
                        // when this merge closes a curve.
                        let travel = if loop_uses[0].0 == first {
                            Sense::Same
                        } else {
                            Sense::Reversed
                        };
                        loop_uses.splice(0..2, [(merged, travel)]);
                        replacements += 1;
                    }
                }
            }
            assert_eq!(replacements, 2, "degree-two edge merge lost a face walk");
            edges[first].active = false;
            edges[second].active = false;
            changed = true;
            break;
        }
        if !changed {
            return Ok(());
        }
    }
}

fn orient_shell(model: &mut Model, faces: &[FaceId], dec: &mut Decider) -> Result<bool, Refusal> {
    let face = faces[0];
    let f = model.face(face).clone();
    crate::points::find(model, face, dec, |model, point, p, dec| {
        for &other in faces {
            if other != face
                && coincident(
                    model.geom().side(model.face(other).surface, point, dec),
                    "boolean.classify-ambiguous",
                )?
                && crate::region::locate(model, other, point, dec)?
                    != crate::region::Location::Outside
            {
                return Ok(None);
            }
        }
        let gradient = oriented(
            &model.geom().surface(f.surface).gradient(&p.to_real()),
            f.sense,
        );
        let normal = RVector3::new(
            gradient
                .x
                .as_rational()
                .expect("polynomial at a rational point"),
            gradient
                .y
                .as_rational()
                .expect("polynomial at a rational point"),
            gradient
                .z
                .as_rational()
                .expect("polynomial at a rational point"),
        );
        let tangent = (0..3)
            .map(|i| normal.cross(&unit(i)))
            .find(|v| !v.is_zero())
            .expect("regular supported surface");
        let second = normal.cross(&tangent);
        for k in 0..crate::geometry::DIRECTION_ATTEMPTS {
            let (a, b) = crate::geometry::probe_parameters(k);
            let direction = normal.plus(&tangent.scaled(&a)).plus(&second.scaled(&b));
            if let Some(inside) = crate::containment::ray_parity(model, faces, p, &direction, dec)?
            {
                return Ok(Some(!inside));
            }
        }
        Ok(None)
    })
}

fn shell_inside(
    model: &mut Model,
    probe_faces: &[FaceId],
    container: &[FaceId],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    crate::points::find(model, probe_faces[0], dec, |model, _, p, dec| {
        match crate::containment::faces(model, container, p, dec) {
            Outcome::Proven(Containment::Inside) => Ok(Some(true)),
            Outcome::Proven(Containment::Outside) => Ok(Some(false)),
            Outcome::ProvenCoincident(_) => Ok(None),
            Outcome::NotExact(n) => Err(refusal("boolean.classify-ambiguous", n.detail)),
            Outcome::Refused(e) => Err(e),
        }
    })
}

pub(crate) fn make_maximal(model: &mut Model, body: BodyId, dec: &mut Decider) -> Answered<BodyId> {
    if model.body(body).kind == BodyKind::Empty {
        return Ok(Answer::exact(body));
    }
    let selection = Selection {
        faces: model
            .faces_of_body(body)
            .into_iter()
            .map(|f| (f, Sense::Same))
            .collect(),
        coincident: BTreeMap::new(),
        marks: MarkLog::new(),
    };
    let marks = model.body(body).marks.clone();
    rebuild(model, &BTreeSet::new(), &selection, marks, Some(body), dec)
}

pub(crate) fn sew(
    model: &mut Model,
    imprinted: &Imprinted,
    selection: &Selection,
    dec: &mut Decider,
) -> Answered<BodyId> {
    let mut marks = selection.marks.clone();
    marks.extend(&model.body(imprinted.a).marks);
    marks.extend(&model.body(imprinted.b).marks);
    rebuild(model, &imprinted.singular, selection, marks, None, dec)
}

fn rebuild(
    model: &mut Model,
    nodes: &BTreeSet<PointId>,
    selection: &Selection,
    marks: MarkLog,
    source: Option<BodyId>,
    dec: &mut Decider,
) -> Answered<BodyId> {
    if selection.faces.is_empty() {
        return Ok(Answer::exact(model.add_body(Body {
            kind: BodyKind::Empty,
            lumps: Vec::new(),
            marks,
        })));
    }
    let (original, uses, by_edge) = read(model, selection);
    let twins = pair(model, &original, &uses, &by_edge, dec)?;
    let (mut faces, removed, next, mut components) =
        merge_faces(model, &original, &uses, &twins, dec)?;
    for list in by_edge.values() {
        if list.len() > 2 {
            let mut seen = BTreeSet::new();
            for &h in list {
                if h < twins[h] && !seen.insert(components.root(uses[h].face)) {
                    return Err(refusal(
                        "boolean.self-touching-shell",
                        "one edge-connected boundary touches itself along an edge",
                    ));
                }
            }
        }
    }
    let mut fans = Groups::new(uses.len());
    for h in 0..uses.len() {
        if !removed[h] && model.edge(uses[h].edge).start.is_some() {
            fans.join(next[h], twins[h]);
        }
    }
    let mut fan_points = BTreeMap::new();
    let mut fan_components = BTreeMap::new();
    let mut by_point: BTreeMap<PointId, BTreeSet<usize>> = BTreeMap::new();
    for (h, u) in uses.iter().enumerate() {
        if removed[h] {
            continue;
        }
        let edge = model.edge(u.edge);
        let vertex = if u.sense == Sense::Same {
            edge.start
        } else {
            edge.end
        };
        if let Some(v) = vertex {
            let root = fans.root(h);
            let point = model.vertex(v).point;
            fan_points.insert(root, point);
            fan_components.insert(root, components.root(u.face));
            by_point.entry(point).or_default().insert(root);
        }
    }
    for list in by_point.values() {
        let mut seen = BTreeSet::new();
        for fan in list {
            if !seen.insert(fan_components[fan]) {
                return Err(refusal(
                    "boolean.self-touching-shell",
                    "kept faces at a vertex form multiple cycles in one connected boundary",
                ));
            }
        }
    }
    let mut face_for_use = BTreeMap::new();
    for (i, f) in faces.iter().enumerate() {
        for l in &f.loops {
            for &(h, _) in l {
                face_for_use.insert(h, i);
            }
        }
    }
    let mut edges = Vec::new();
    let mut edge_for_use = BTreeMap::new();
    for h in 0..uses.len() {
        if removed[h] || h > twins[h] {
            continue;
        }
        let u = &uses[h];
        let edge = model.edge(u.edge);
        let partner = twins[h];
        let (from, to) = if edge.start.is_none() {
            (None, None)
        } else if u.sense == Sense::Same {
            (Some(fans.root(h)), Some(fans.root(next[h])))
        } else {
            (Some(fans.root(next[h])), Some(fans.root(h)))
        };
        let mut pair = [face_for_use[&h], face_for_use[&partner]];
        pair.sort();
        let index = edges.len();
        edges.push(EdgeData {
            curve: edge.curve,
            sense: edge.sense,
            from,
            to,
            faces: pair,
            active: true,
        });
        edge_for_use.insert(h, index);
        edge_for_use.insert(partner, index);
    }
    for face in &mut faces {
        for l in &mut face.loops {
            for item in l {
                let h = item.0;
                item.0 = edge_for_use[&h];
                item.1 = uses[h].sense;
            }
        }
    }
    maximal_edges(model, &mut faces, &mut edges, &fan_points, nodes, dec)?;
    let mut vertices = BTreeMap::new();
    for edge in &edges {
        if edge.active {
            for fan in [edge.from, edge.to].into_iter().flatten() {
                vertices.entry(fan).or_insert_with(|| {
                    model.add_vertex(Vertex {
                        point: fan_points[&fan],
                    })
                });
            }
        }
    }
    let mut stored_edges = BTreeMap::new();
    for (i, edge) in edges.iter().enumerate() {
        if !edge.active {
            continue;
        }
        let start = edge.from.map(|v| vertices[&v]);
        let end = edge.to.map(|v| vertices[&v]);
        for fan in [edge.from, edge.to].into_iter().flatten() {
            record_incident(
                model,
                GeomKey::Point(fan_points[&fan]),
                GeomKey::Curve(edge.curve),
            );
        }
        stored_edges.insert(
            i,
            model.add_edge(Edge {
                curve: edge.curve,
                sense: edge.sense,
                start,
                end,
            }),
        );
    }
    let mut stored_faces = Vec::new();
    for face in &faces {
        let mut loops = Vec::new();
        for l in &face.loops {
            let mut cs = Vec::new();
            for &(edge, sense) in l {
                record_incident(
                    model,
                    GeomKey::Curve(edges[edge].curve),
                    GeomKey::Surface(face.surface),
                );
                cs.push(model.add_coedge(Coedge {
                    edge: stored_edges[&edge],
                    sense,
                }));
            }
            if !cs.is_empty() {
                loops.push(model.add_loop(Loop::Coedges(cs)));
            }
        }
        for &vertex in &face.vertices {
            let point = model.vertex(vertex).point;
            record_incident(model, GeomKey::Point(point), GeomKey::Surface(face.surface));
            loops.push(model.add_loop(Loop::Vertex(vertex)));
        }
        stored_faces.push(model.add_named_face(
            Face {
                surface: face.surface,
                sense: face.sense,
                loops,
            },
            face.origins.clone(),
        ));
    }
    if let Some(source) = source {
        // Removing redundant boundaries does not change shell orientation or
        // containment. Retain the constructor's proven lump/void partition.
        let mut owner = BTreeMap::new();
        for &lump in &model.body(source).lumps {
            for &shell in &model.lump(lump).shells {
                for &face in &model.shell(shell).faces {
                    owner.insert(face, shell);
                }
            }
        }
        let mut grouped: BTreeMap<_, Vec<FaceId>> = BTreeMap::new();
        for (face, &stored) in faces.iter().zip(&stored_faces) {
            let original = selection.faces[face.original].0;
            grouped.entry(owner[&original]).or_default().push(stored);
        }
        let mut lumps = Vec::with_capacity(model.body(source).lumps.len());
        for i in 0..model.body(source).lumps.len() {
            let original_lump = model.body(source).lumps[i];
            let mut shells = Vec::with_capacity(model.lump(original_lump).shells.len());
            for j in 0..model.lump(original_lump).shells.len() {
                let original_shell = model.lump(original_lump).shells[j];
                shells.push(
                    model.add_shell(Shell {
                        faces: grouped
                            .remove(&original_shell)
                            .expect("a valid shell retains a boundary"),
                    }),
                );
            }
            lumps.push(model.add_lump(Lump { shells }));
        }
        let body = model.add_body(Body {
            kind: model.body(source).kind,
            lumps,
            marks,
        });
        return checked(model, body);
    }
    let mut shell_faces: BTreeMap<usize, Vec<FaceId>> = BTreeMap::new();
    for (i, face) in faces.iter().enumerate() {
        shell_faces
            .entry(components.root(face.original))
            .or_default()
            .push(stored_faces[i]);
    }
    let shell_faces: Vec<_> = shell_faces.into_iter().collect();
    let mut outer = Vec::new();
    let mut voids = Vec::new();
    for (i, (_, fs)) in shell_faces.iter().enumerate() {
        if orient_shell(model, fs, dec)? {
            outer.push(i);
        } else {
            voids.push(i);
        }
    }
    let mut lump_shells: BTreeMap<usize, Vec<usize>> =
        outer.iter().map(|&i| (i, vec![i])).collect();
    let mut shell_lump: BTreeMap<usize, usize> = outer.iter().map(|&i| (i, i)).collect();
    for void in voids {
        let mut enclosing: Option<usize> = None;
        for &candidate in &outer {
            if shell_inside(model, &shell_faces[void].1, &shell_faces[candidate].1, dec)? {
                if let Some(current) = enclosing {
                    if shell_inside(
                        model,
                        &shell_faces[candidate].1,
                        &shell_faces[current].1,
                        dec,
                    )? {
                        enclosing = Some(candidate);
                    }
                } else {
                    enclosing = Some(candidate);
                }
            }
        }
        let outer = enclosing.ok_or_else(|| {
            refusal(
                "boolean.classify-ambiguous",
                "void shell has no proven enclosing outer shell",
            )
        })?;
        lump_shells
            .get_mut(&outer)
            .expect("outer shell starts a lump")
            .push(void);
        shell_lump.insert(void, outer);
    }
    let component_lump: BTreeMap<_, _> = shell_faces
        .iter()
        .enumerate()
        .map(|(i, (component, _))| (*component, shell_lump[&i]))
        .collect();
    for list in by_point.values() {
        let mut seen = BTreeSet::new();
        for fan in list {
            if !seen.insert(component_lump[&fan_components[fan]]) {
                return Err(refusal(
                    "boolean.self-touching-shell",
                    "vertex contact joins two boundary cycles of one lump",
                ));
            }
        }
    }
    for list in by_edge.values() {
        if list.len() > 2 {
            let mut seen = BTreeSet::new();
            for &h in list {
                if h < twins[h] && !seen.insert(component_lump[&components.root(uses[h].face)]) {
                    return Err(refusal(
                        "boolean.self-touching-shell",
                        "curve contact belongs to one connected material lump",
                    ));
                }
            }
        }
    }
    let mut lumps = Vec::new();
    for shells in lump_shells.into_values() {
        let shells = shells
            .into_iter()
            .map(|i| {
                model.add_shell(Shell {
                    faces: shell_faces[i].1.clone(),
                })
            })
            .collect();
        lumps.push(model.add_lump(Lump { shells }));
    }
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps,
        marks,
    });
    checked(model, body)
}

fn checked(model: &Model, body: BodyId) -> Answered<BodyId> {
    let faults = model.check(body);
    assert!(
        faults.is_empty(),
        "boolean.sew produced invalid topology: {faults:?}"
    );
    Ok(Answer::exact(body))
}
