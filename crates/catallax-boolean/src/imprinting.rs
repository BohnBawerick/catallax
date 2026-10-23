//! A face-local halfedge arrangement. Boundary walks are grouped by winding
//! signatures, so holes, cylinder bands and noncontractible cuts use one rule.

use crate::geometry::{
    canonical_point, coincident, natural_ends, on_arc, oriented, refusal, reverse, sign,
};
use crate::graph::sort_points;
use crate::points::{curve_sample, left_probe};
use crate::region::{Location, locate, signatures};
use crate::stages::{Imprinted, IntersectionGraph};
use catallax_intersect::{curve_is_closed, same_curve, tangent_at};
use catallax_numbers::{
    Answer, Answered, Decider, Fact, FactKind, FactOrigin, Refusal, Sign, ZeroProof,
};
use catallax_shapes::{CurveId, GeomKey, PointId, Surface};
use catallax_topology::{
    Coedge, Edge, EdgeId, Face, FaceId, Loop, LoopId, Model, Sense, Vertex, VertexId,
};
use std::collections::{BTreeMap, BTreeSet};

struct Builder {
    points: Vec<PointId>,
    vertices: BTreeMap<PointId, VertexId>,
    edges: Vec<EdgeId>,
}
impl Builder {
    fn point(
        &mut self,
        model: &mut Model,
        p: PointId,
        dec: &mut Decider,
    ) -> Result<PointId, Refusal> {
        canonical_point(model, &mut self.points, p, dec, "boolean.imprint-ambiguous")
    }
    fn vertex(&mut self, model: &mut Model, p: PointId) -> VertexId {
        *self
            .vertices
            .entry(p)
            .or_insert_with(|| model.add_vertex(Vertex { point: p }))
    }
    fn edge(
        &mut self,
        model: &mut Model,
        curve: CurveId,
        a: Option<PointId>,
        b: Option<PointId>,
        dec: &mut Decider,
    ) -> Result<(EdgeId, Sense), Refusal> {
        for &edge in &self.edges {
            let other = model.edge(edge).clone();
            if !coincident(
                same_curve(model.geom(), curve, other.curve, dec),
                "boolean.imprint-ambiguous",
            )? {
                continue;
            }
            let (x, y) = natural_ends(model, edge);
            let sample = curve_sample(model, curve, a, b, dec)?;
            let ours = tangent_at(model.geom(), curve, sample, dec)?[1].clone();
            let theirs = tangent_at(model.geom(), other.curve, sample, dec)?[1].clone();
            let same = match sign(&ours.dot(&theirs), dec, "boolean.imprint-ambiguous")? {
                Sign::Positive => true,
                Sign::Negative => false,
                Sign::Zero => {
                    return Err(refusal(
                        "boolean.imprint-ambiguous",
                        "coincident curve directions vanished",
                    ));
                }
            };
            if (same && a == x && b == y) || (!same && a == y && b == x) {
                model.geom_mut().record_fact(Fact {
                    kind: FactKind::Coincident,
                    a: GeomKey::Curve(curve),
                    b: GeomKey::Curve(other.curve),
                    at: None,
                    origin: FactOrigin::Proof(ZeroProof::Algebraic),
                });
                return Ok((edge, if same { Sense::Same } else { Sense::Reversed }));
            }
        }
        let start = a.map(|p| self.vertex(model, p));
        let end = b.map(|p| self.vertex(model, p));
        let edge = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start,
            end,
        });
        self.edges.push(edge);
        Ok((edge, Sense::Same))
    }
    fn split(
        &mut self,
        model: &mut Model,
        edge: EdgeId,
        dec: &mut Decider,
    ) -> Result<Vec<(EdgeId, Sense)>, Refusal> {
        let original = model.edge(edge).clone();
        let (from, to) = natural_ends(model, edge);
        let mut result = self.split_arc(model, original.curve, from, to, dec)?;
        if original.sense == Sense::Reversed {
            result.reverse();
            for (_, s) in &mut result {
                *s = reverse(*s);
            }
        }
        Ok(result)
    }
    fn split_arc(
        &mut self,
        model: &mut Model,
        curve: CurveId,
        from: Option<PointId>,
        to: Option<PointId>,
        dec: &mut Decider,
    ) -> Result<Vec<(EdgeId, Sense)>, Refusal> {
        let from = from.map(|p| self.point(model, p, dec)).transpose()?;
        let to = to.map(|p| self.point(model, p, dec)).transpose()?;
        let mut cuts = Vec::new();
        if let Some(p) = from {
            cuts.push(p);
        }
        for &p in &self.points {
            if Some(p) != from
                && on_arc(model, curve, from, to, p, dec, "boolean.imprint-ambiguous")?
            {
                cuts.push(p);
            }
        }
        sort_points(model, curve, &mut cuts, dec)?;
        let whole = curve_is_closed(model.geom(), curve, dec)? && (from.is_none() || from == to);
        let arcs = if whole && cuts.is_empty() {
            vec![(None, None)]
        } else if whole {
            (0..cuts.len())
                .map(|i| (Some(cuts[i]), Some(cuts[(i + 1) % cuts.len()])))
                .collect()
        } else {
            cuts.windows(2).map(|w| (Some(w[0]), Some(w[1]))).collect()
        };
        let mut result = Vec::new();
        for (a, b) in arcs {
            result.push(self.edge(model, curve, a, b, dec)?);
        }
        Ok(result)
    }
}

fn ends(model: &Model, edge: EdgeId, sense: Sense) -> (Option<VertexId>, Option<VertexId>) {
    let e = model.edge(edge);
    if sense == Sense::Same {
        (e.start, e.end)
    } else {
        (e.end, e.start)
    }
}

fn walks(
    model: &mut Model,
    face: FaceId,
    halfedges: &[(EdgeId, Sense)],
    dec: &mut Decider,
) -> Result<Vec<LoopId>, Refusal> {
    let f = model.face(face).clone();
    let mut outgoing: BTreeMap<VertexId, Vec<usize>> = BTreeMap::new();
    for (i, &(edge, sense)) in halfedges.iter().enumerate() {
        if let Some(start) = ends(model, edge, sense).0 {
            outgoing.entry(start).or_default().push(i);
        }
    }
    let mut successors = Vec::with_capacity(halfedges.len());
    for (index, &(edge, sense)) in halfedges.iter().enumerate() {
        let Some(end) = ends(model, edge, sense).1 else {
            successors.push(index);
            continue;
        };
        let point = model.vertex(end).point;
        let cone_axis = if let Surface::Cone(cone) = model.geom().surface(f.surface) {
            let displacement = model.geom().position(point).minus(&cone.apex);
            if sign(
                &displacement.norm_squared(),
                dec,
                "boolean.imprint-ambiguous",
            )? == Sign::Zero
            {
                Some(cone.axis.clone())
            } else {
                None
            }
        } else {
            None
        };
        let normal = if let Some(axis) = &cone_axis {
            // The apex link is the radial disk, oriented opposite the axis.
            oriented(&axis.scaled(&catallax_numbers::Real::from_i64(-1)), f.sense)
        } else {
            oriented(&model.geom().normal(f.surface, point, dec)?, f.sense)
        };
        let candidates = outgoing.get(&end).ok_or_else(|| {
            refusal(
                "boolean.imprint-ambiguous",
                "imprint has an open face boundary",
            )
        })?;
        let mut germs = Vec::new();
        for &candidate in candidates {
            let (e, s) = halfedges[candidate];
            germs.push(crate::germ::outgoing(
                model,
                e,
                s,
                point,
                candidate,
                dec,
                "boolean.imprint-ambiguous",
            )?);
        }
        let reference = candidates
            .iter()
            .copied()
            .find(|&i| halfedges[i] == (edge, reverse(sense)));
        let reference = if let Some(i) = reference {
            i
        } else {
            let key = halfedges.len();
            germs.push(crate::germ::outgoing(
                model,
                edge,
                reverse(sense),
                point,
                key,
                dec,
                "boolean.imprint-ambiguous",
            )?);
            key
        };
        if let Some(axis) = &cone_axis {
            let inverse = catallax_numbers::Real::one() / axis.norm_squared();
            for germ in &mut germs {
                germ.tangent = germ
                    .tangent
                    .minus(&axis.scaled(&(&germ.tangent.dot(axis) * &inverse)));
                germ.curvature = germ
                    .curvature
                    .minus(&axis.scaled(&(&germ.curvature.dot(axis) * &inverse)));
            }
        }
        crate::germ::sort(model, &mut germs, &normal, dec, "boolean.imprint-ambiguous")?;
        let at = germs
            .iter()
            .position(|g| g.key == reference)
            .expect("incoming germ is present");
        let next = germs[(at + germs.len() - 1) % germs.len()].key;
        if next >= halfedges.len() {
            return Err(refusal(
                "boolean.imprint-ambiguous",
                "a face boundary has no outgoing halfedge",
            ));
        }
        successors.push(next);
    }
    let mut visited = BTreeSet::new();
    let mut loops = Vec::new();
    for seed in 0..halfedges.len() {
        if visited.contains(&seed) {
            continue;
        }
        let mut cs = Vec::new();
        let mut current = seed;
        loop {
            if !visited.insert(current) {
                if current == seed {
                    break;
                }
                return Err(refusal(
                    "boolean.imprint-ambiguous",
                    "radial order did not form disjoint boundary walks",
                ));
            }
            let (edge, sense) = halfedges[current];
            cs.push(model.add_coedge(Coedge { edge, sense }));
            current = successors[current];
        }
        loops.push(model.add_loop(Loop::Coedges(cs)));
    }
    Ok(loops)
}

fn pieces(
    model: &mut Model,
    face: FaceId,
    boundary: &BTreeMap<EdgeId, Vec<(EdgeId, Sense)>>,
    sections: &[(EdgeId, Sense)],
    dec: &mut Decider,
) -> Result<Vec<FaceId>, Refusal> {
    let f = model.face(face).clone();
    let mut halfedges = Vec::new();
    let mut vertex_loops = Vec::new();
    for &l in &f.loops {
        match model.loop_(l) {
            Loop::Coedges(cs) => {
                for &c in cs {
                    let coedge = model.coedge(c);
                    let mut split = boundary[&coedge.edge].clone();
                    if coedge.sense == Sense::Reversed {
                        split.reverse();
                        for (_, s) in &mut split {
                            *s = reverse(*s);
                        }
                    }
                    halfedges.extend(split);
                }
            }
            Loop::Vertex(_) => vertex_loops.push(l),
        }
    }
    let boundary_edges: BTreeSet<_> = halfedges.iter().map(|x| x.0).collect();
    for &(edge, _) in sections {
        if !boundary_edges.contains(&edge) {
            halfedges.push((edge, Sense::Same));
            halfedges.push((edge, Sense::Reversed));
        }
    }
    if halfedges.is_empty() {
        return Ok(vec![face]);
    }
    let mut distinct = BTreeSet::new();
    halfedges.retain(|x| distinct.insert(*x));
    let loops = walks(model, face, &halfedges, dec)?;
    let groups: Vec<_> = loops.iter().map(|&l| vec![l]).collect();
    let barriers: Vec<_> = halfedges
        .iter()
        .map(|x| x.0)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut cells: BTreeMap<Vec<i64>, Vec<LoopId>> = BTreeMap::new();
    for &l in &loops {
        let Loop::Coedges(original_coedges) = model.loop_(l) else {
            unreachable!("arrangement walks use coedges")
        };
        let uses: Vec<_> = original_coedges
            .iter()
            .map(|&c| model.coedge(c).clone())
            .collect();
        let temporary_coedges = uses.into_iter().map(|c| model.add_coedge(c)).collect();
        let temporary_loop = model.add_loop(Loop::Coedges(temporary_coedges));
        let temporary = model.add_face(Face {
            surface: f.surface,
            sense: f.sense,
            loops: vec![temporary_loop],
        });
        let surface = model.geom().surface(f.surface).clone();
        let (inside, signature) =
            left_probe(model, temporary, &barriers, dec, |model, point, _, dec| {
                let inside = locate(model, face, point, dec)? == Location::Inside;
                let signature = signatures(model, &surface, &groups, point, dec)?;
                Ok(Some((inside, signature)))
            })?;
        if inside {
            cells.entry(signature).or_default().push(l);
        }
    }
    let origins = model.face_origins(face).to_vec();
    let mut result = Vec::with_capacity(cells.len());
    for mut loops in cells.into_values() {
        if !vertex_loops.is_empty() {
            let provisional = model.add_face(Face {
                surface: f.surface,
                sense: f.sense,
                loops: loops.clone(),
            });
            for &l in &vertex_loops {
                let Loop::Vertex(vertex) = model.loop_(l) else {
                    unreachable!("point-only boundary")
                };
                let point = model.vertex(*vertex).point;
                if locate(model, provisional, point, dec)? == Location::Inside {
                    loops.push(l);
                }
                // A boundary hit is already represented by an incident edge.
            }
        }
        result.push(model.add_named_face(
            Face {
                surface: f.surface,
                sense: f.sense,
                loops,
            },
            origins.clone(),
        ));
    }
    if result.is_empty() {
        return Err(refusal(
            "boolean.imprint-ambiguous",
            "face arrangement lost every interior cell",
        ));
    }
    Ok(result)
}

pub(crate) fn imprint(
    model: &mut Model,
    graph: &IntersectionGraph,
    dec: &mut Decider,
) -> Answered<Imprinted> {
    let mut builder = Builder {
        points: Vec::new(),
        vertices: BTreeMap::new(),
        edges: Vec::new(),
    };
    for body in [graph.a, graph.b] {
        for vertex in model.vertices_of_body(body) {
            builder.point(model, model.vertex(vertex).point, dec)?;
        }
    }
    for &point in graph.points.iter().chain(graph.singular.iter()) {
        builder.point(model, point, dec)?;
    }
    for edge in &graph.edges {
        for point in [edge.start, edge.end].into_iter().flatten() {
            builder.point(model, point, dec)?;
        }
    }
    let mut boundary = BTreeMap::new();
    for body in [graph.a, graph.b] {
        for edge in model.edges_of_body(body) {
            if let std::collections::btree_map::Entry::Vacant(entry) = boundary.entry(edge) {
                let parts = builder.split(model, edge, dec)?;
                entry.insert(parts);
            }
        }
    }
    let mut sections: BTreeMap<FaceId, Vec<(EdgeId, Sense)>> = BTreeMap::new();
    let mut common = Vec::new();
    for edge in &graph.edges {
        for stored in builder.split_arc(model, edge.curve, edge.start, edge.end, dec)? {
            if !common.contains(&stored.0) {
                common.push(stored.0);
            }
            for &face in edge.faces_of_first.iter().chain(&edge.faces_of_second) {
                sections.entry(face).or_default().push(stored);
            }
        }
    }
    let mut all = Vec::new();
    for body in [graph.a, graph.b] {
        let mut body_pieces = Vec::new();
        for face in model.faces_of_body(body) {
            body_pieces.extend(pieces(
                model,
                face,
                &boundary,
                sections.get(&face).map(Vec::as_slice).unwrap_or(&[]),
                dec,
            )?);
        }
        all.push(body_pieces);
    }
    let second = all.pop().expect("two operands");
    let first = all.pop().expect("two operands");
    Ok(Answer::exact(Imprinted {
        a: graph.a,
        b: graph.b,
        first,
        second,
        common,
        singular: graph.singular.clone(),
        marks: graph.marks.clone(),
    }))
}
