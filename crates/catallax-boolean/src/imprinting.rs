//! A face-local halfedge arrangement. Boundary walks are grouped by winding
//! signatures, so holes, cylinder bands and noncontractible cuts use one rule.

use crate::geometry::{
    canonical_point, coincident, natural_ends, on_arc, oriented, refusal, reverse, sign,
};
use crate::graph::sort_points;
use crate::points::{curve_sample, left_probe};
use crate::region::{Location, RegionCache, locate};
use crate::stages::{Imprinted, IntersectionGraph};
use catallax_intersect::{curve_is_closed, same_curve, tangent_at};
use catallax_numbers::{
    Answer, Answered, Decider, Fact, FactKind, FactOrigin, Outcome, Refusal, Sign, ZeroProof,
};
use catallax_shapes::{Curve, CurveId, GeomKey, PointId, Surface};
use catallax_topology::{
    Coedge, Edge, EdgeId, Face, FaceId, Loop, LoopId, Model, Sense, Vertex, VertexId,
};
use std::collections::{BTreeMap, BTreeSet};

mod bezier;

/// The sense in which stored `edge` is the arc of `curve` between `a` and
/// `b`, when it is that arc.
fn reusable(
    model: &mut Model,
    fitted: bool,
    edge: EdgeId,
    curve: CurveId,
    (a, b): (Option<PointId>, Option<PointId>),
    dec: &mut Decider,
) -> Result<Option<Sense>, Refusal> {
    let other = model.edge(edge).clone();
    let verdict = same_curve(model.geom(), curve, other.curve, dec);
    let (x, y) = natural_ends(model, edge);
    // Only an edge with these two ends is reused below. Without a
    // spline equivalence proof, a different-ended edge stays distinct.
    if fitted
        && let Outcome::Refused(refused) = &verdict
        && refused.code == "intersect.parametric-coincidence"
        && !((a == x && b == y) || (a == y && b == x))
    {
        return Ok(None);
    }
    if !coincident(verdict, "boolean.imprint-ambiguous")? {
        return Ok(None);
    }
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
        return Ok(Some(if same { Sense::Same } else { Sense::Reversed }));
    }
    Ok(None)
}

#[derive(Clone)]
struct Builder {
    // Whether an operand has a fitted band, whose spline edges may lack a
    // coincidence proof against an edge with other ends.
    fitted: bool,
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
        // The first reusable edge ends the search; workers compare against
        // later edges ahead. Merging an edge into the shared list is an
        // ordered append.
        #[cfg(feature = "work-span")]
        let _append = catallax_numbers::work_span::append();
        #[cfg(feature = "work-span")]
        let span_tasks = catallax_numbers::work_span::search("edge merge");
        let mut ahead = catallax_numbers::run_ahead!(dec, self.edges.len());
        for (index, &edge) in self.edges.iter().enumerate() {
            #[cfg(feature = "work-span")]
            let _span_task = span_tasks.task();
            ahead.at(index, || {
                let (model, edges, fitted) = (model.clone(), self.edges.clone(), self.fitted);
                move |index: usize, dec: &mut Decider| {
                    let _ = reusable(&mut model.clone(), fitted, edges[index], curve, (a, b), dec);
                }
            });
            if let Some(sense) = reusable(model, self.fitted, edge, curve, (a, b), dec)? {
                return Ok((edge, sense));
            }
        }
        drop(ahead);
        #[cfg(feature = "work-span")]
        drop(span_tasks);
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
        // Whether a point lies on the arc is a value of what every thread
        // shares, so other threads test later points whole.
        let mut ahead = catallax_numbers::values_ahead!(dec, self.points.len());
        for (index, &p) in self.points.iter().enumerate() {
            if Some(p) == from {
                continue;
            }
            let adopted = ahead.take(index, dec, || {
                let (model, points) = (model.clone(), self.points.clone());
                move |index: usize, dec: &mut Decider| {
                    let p = points[index];
                    // The loop skips the arc's own start; its value is not read.
                    if Some(p) == from {
                        return Ok(false);
                    }
                    on_arc(&model, curve, from, to, p, dec, "boolean.imprint-ambiguous")
                }
            });
            let on = match adopted {
                Some(on) => on,
                None => on_arc(model, curve, from, to, p, dec, "boolean.imprint-ambiguous"),
            }?;
            if on {
                cuts.push(p);
            }
        }
        drop(ahead);
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

fn named_germ(model: &Model, edge: EdgeId, mut refusal: Refusal) -> Refusal {
    if refusal.code == "boolean.parametric-query"
        && refusal.detail.starts_with("curve equation unavailable")
    {
        let kind = match model.geom().curve(model.edge(edge).curve) {
            Curve::BSpline(_) => "spline",
            Curve::PatchSection(_) => "patch-section",
            Curve::PatchPair(_) => "patch-pair",
            Curve::Intersection(_) => "traced",
            _ => "parametric",
        };
        refusal.detail = format!("outgoing germ of a {kind} edge: {}", refusal.detail);
    }
    refusal
}

fn walks(
    model: &mut Model,
    face: FaceId,
    halfedges: &[(EdgeId, Sense)],
    dec: &mut Decider,
    context: &mut crate::germ::Context,
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
            if sign(&displacement.norm_squared(), dec, context.code)? == Sign::Zero {
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
        let candidates = outgoing
            .get(&end)
            .ok_or_else(|| refusal(context.code, "imprint has an open face boundary"))?;
        let mut germs = Vec::new();
        for &candidate in candidates {
            let (e, s) = halfedges[candidate];
            germs.push(
                crate::germ::outgoing(model, e, s, point, candidate, dec, context)
                    .map_err(|refusal| named_germ(model, e, refusal))?,
            );
        }
        let reference = candidates
            .iter()
            .copied()
            .find(|&i| halfedges[i] == (edge, reverse(sense)));
        let reference = if let Some(i) = reference {
            i
        } else {
            let key = halfedges.len();
            germs.push(
                crate::germ::outgoing(model, edge, reverse(sense), point, key, dec, context)
                    .map_err(|refusal| named_germ(model, edge, refusal))?,
            );
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

/// Whether the cell left of walk `l` lies in the face, and its winding
/// signature against every walk.
fn cell(
    model: &mut Model,
    (face, f): (FaceId, &Face),
    l: LoopId,
    barriers: &[EdgeId],
    groups: &[Vec<LoopId>],
    regions: &mut RegionCache,
    dec: &mut Decider,
) -> Result<(bool, Vec<i64>), Refusal> {
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
    left_probe(model, temporary, barriers, dec, |model, point, _, dec| {
        let inside = locate(model, face, point, dec)? == Location::Inside;
        let signature = regions.signatures_prepared(model, f.surface, groups, point, dec)?;
        Ok(Some((inside, signature)))
    })
}

fn pieces(
    model: &mut Model,
    face: FaceId,
    boundary: &BTreeMap<EdgeId, Vec<(EdgeId, Sense)>>,
    sections: &[(EdgeId, Sense)],
    dec: &mut Decider,
    context: &mut crate::germ::Context,
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
    let loops = walks(model, face, &halfedges, dec, context)?;
    if let Some(pieces) = bezier::pieces(
        model,
        face,
        &loops,
        &boundary_edges,
        sections,
        &vertex_loops,
        dec,
    )? {
        return Ok(pieces);
    }
    let groups: Vec<_> = loops.iter().map(|&l| vec![l]).collect();
    let barriers: Vec<_> = halfedges
        .iter()
        .map(|x| x.0)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut cells: BTreeMap<Vec<i64>, Vec<LoopId>> = BTreeMap::new();
    let mut regions = RegionCache::default();
    // Each walk's cell is found on its own. Other threads find later cells
    // whole, each on a copy of the model and with the region cache as they
    // were when the loop was handed over. This thread takes such a result
    // when the cell would have been found the same on the model as it is
    // now, and the other thread's copy of the cache held all the cell asked
    // for. This thread's cache holds that too, so the cell costs the same.
    #[cfg(feature = "work-span")]
    let span_tasks = catallax_numbers::work_span::parallel("face cells");
    let mut ahead = catallax_numbers::values_ahead!(dec, loops.len());
    for (index, &l) in loops.iter().enumerate() {
        #[cfg(feature = "work-span")]
        let _span_task = span_tasks.task();
        let planned = ahead.take_if(
            index,
            dec,
            || {
                let (base, f, loops, barriers, groups, regions) = (
                    std::sync::Arc::new(model.clone()),
                    f.clone(),
                    loops.clone(),
                    barriers.clone(),
                    groups.clone(),
                    regions.clone(),
                );
                move |index: usize, dec: &mut Decider| {
                    regions.ahead(&base, |model, regions| {
                        cell(
                            model,
                            (face, &f),
                            loops[index],
                            &barriers,
                            &groups,
                            regions,
                            dec,
                        )
                    })
                }
            },
            |ahead| ahead.holds(model),
        );
        let (inside, signature) = match planned {
            Some(ahead) => ahead.commit(model),
            None => crate::geometry::on_a_copy(model, |model| {
                cell(model, (face, &f), l, &barriers, &groups, &mut regions, dec)
            }),
        }?;
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
    imprint_graph(
        model,
        graph,
        dec,
        &mut crate::germ::Context::new("boolean.imprint-ambiguous"),
    )
}

pub(crate) fn imprint_for_boolean(
    model: &mut Model,
    graph: &IntersectionGraph,
    dec: &mut Decider,
    proof: &mut Option<crate::trim_bounds::BodyApartProof>,
) -> Answered<Imprinted> {
    let mut context = crate::germ::Context::new("boolean.imprint-ambiguous");
    context.fitted_operands = Some((graph.a, graph.b));
    let answer = imprint_graph(model, graph, dec, &mut context);
    *proof = context.proven_apart;
    answer
}

fn imprint_graph(
    model: &mut Model,
    graph: &IntersectionGraph,
    dec: &mut Decider,
    context: &mut crate::germ::Context,
) -> Answered<Imprinted> {
    let first = model.faces_of_body(graph.a);
    let second = model.faces_of_body(graph.b);
    imprint_groups(
        model,
        graph,
        &[(graph.a, &first), (graph.b, &second)],
        dec,
        context,
        |_, _| {},
    )
}

pub(crate) fn arrangement(
    model: &mut Model,
    graph: &IntersectionGraph,
    faces: &[FaceId],
    dec: &mut Decider,
) -> Answered<(Imprinted, BTreeMap<FaceId, FaceId>)> {
    let mut owners = BTreeMap::new();
    let mut context = crate::germ::Context::new("boolean.imprint-ambiguous");
    let answer = imprint_groups(
        model,
        graph,
        &[(graph.a, faces)],
        dec,
        &mut context,
        |face, pieces| {
            owners.extend(pieces.iter().map(|&piece| (piece, face)));
        },
    )?;
    Ok(Answer::with_marks((answer.value, owners), answer.marks))
}

fn imprint_groups(
    model: &mut Model,
    graph: &IntersectionGraph,
    groups: &[(catallax_topology::BodyId, &[FaceId])],
    dec: &mut Decider,
    context: &mut crate::germ::Context,
    mut visit: impl FnMut(FaceId, &[FaceId]),
) -> Answered<Imprinted> {
    let mut builder = Builder {
        fitted: groups.iter().any(|(_, faces)| {
            faces
                .iter()
                .any(|face| crate::trim_bounds::fitted_band(model, *face))
        }),
        points: Vec::new(),
        vertices: BTreeMap::new(),
        edges: Vec::new(),
    };
    let mut additional_edges = None;
    if !builder.fitted || graph.a == graph.b {
        context.fitted_operands = None;
    }
    if graph.a == graph.b {
        let mut vertices = BTreeSet::new();
        let mut edges = BTreeSet::new();
        for &(_, faces) in groups {
            for &face in faces {
                for lp in &model.face(face).loops {
                    match model.loop_(*lp) {
                        Loop::Vertex(vertex) => {
                            vertices.insert(*vertex);
                        }
                        Loop::Coedges(coedges) => {
                            for coedge in coedges {
                                let edge = model.coedge(*coedge).edge;
                                edges.insert(edge);
                                let edge = model.edge(edge);
                                vertices.extend([edge.start, edge.end].into_iter().flatten());
                            }
                        }
                    }
                }
            }
        }
        for vertex in vertices {
            builder.point(model, model.vertex(vertex).point, dec)?;
        }
        additional_edges = Some(edges.into_iter().collect::<Vec<_>>());
    } else {
        for body in [graph.a, graph.b] {
            for vertex in model.vertices_of_body(body) {
                builder.point(model, model.vertex(vertex).point, dec)?;
            }
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
    // Each boundary edge once, under the first body that has it. A body's own
    // edges do not change while its edges are split.
    let mut boundary_edges = Vec::new();
    let mut listed = BTreeSet::new();
    for &(body, _) in groups {
        for edge in additional_edges
            .take()
            .unwrap_or_else(|| model.edges_of_body(body))
        {
            if listed.insert(edge) {
                boundary_edges.push((body, edge));
            }
        }
    }
    let mut boundary = BTreeMap::new();
    // A split merges its cuts into the builder's lists, so workers split later
    // edges ahead against the lists as they stood.
    #[cfg(feature = "work-span")]
    let span_tasks = catallax_numbers::work_span::parallel("split boundary edges");
    let mut ahead = catallax_numbers::run_ahead!(dec, boundary_edges.len());
    for (index, &(body, edge)) in boundary_edges.iter().enumerate() {
        #[cfg(feature = "work-span")]
        let _span_task = span_tasks.task();
        ahead.at(index, || {
            let (model, builder, edges) = (model.clone(), builder.clone(), boundary_edges.clone());
            move |index: usize, dec: &mut Decider| {
                let _ = builder
                    .clone()
                    .split(&mut model.clone(), edges[index].1, dec);
            }
        });
        let parts = builder.split(model, edge, dec)?;
        let mut names = model.edge_constructions(edge).to_vec();
        if names.is_empty() {
            names.extend(model.edge_name(body, edge));
        }
        for &(part, _) in &parts {
            model.record_edge_constructions(part, names.clone());
        }
        boundary.insert(edge, parts);
    }
    drop(ahead);
    #[cfg(feature = "work-span")]
    drop(span_tasks);
    let mut sections: BTreeMap<FaceId, Vec<(EdgeId, Sense)>> = BTreeMap::new();
    let mut common = Vec::new();
    #[cfg(feature = "work-span")]
    let span_tasks = catallax_numbers::work_span::parallel("split section edges");
    let mut ahead = catallax_numbers::run_ahead!(dec, graph.edges.len());
    for (index, edge) in graph.edges.iter().enumerate() {
        #[cfg(feature = "work-span")]
        let _span_task = span_tasks.task();
        ahead.at(index, || {
            let (model, builder, edges) = (model.clone(), builder.clone(), graph.edges.clone());
            move |index: usize, dec: &mut Decider| {
                let edge = &edges[index];
                let _ = builder.clone().split_arc(
                    &mut model.clone(),
                    edge.curve,
                    edge.start,
                    edge.end,
                    dec,
                );
            }
        });
        for stored in builder.split_arc(model, edge.curve, edge.start, edge.end, dec)? {
            if !common.contains(&stored.0) {
                common.push(stored.0);
            }
            for &face in edge.faces_of_first.iter().chain(&edge.faces_of_second) {
                sections.entry(face).or_default().push(stored);
            }
        }
    }
    drop(ahead);
    #[cfg(feature = "work-span")]
    drop(span_tasks);
    let mut all = Vec::with_capacity(groups.len());
    // Each face is cut on its own. Other threads cut later faces whole, each
    // on a copy of the model as it was when the loop was handed over. This
    // thread takes such a result when the face would have been cut the same
    // on the model as it is now, and cuts the face itself otherwise.
    let every_face: Vec<FaceId> = groups
        .iter()
        .flat_map(|(_, faces)| faces.iter().copied())
        .collect();
    #[cfg(feature = "work-span")]
    let span_tasks = catallax_numbers::work_span::parallel("face pieces");
    let mut ahead = catallax_numbers::values_ahead!(dec, every_face.len());
    let mut index = 0;
    for &(_, faces) in groups {
        let mut group_pieces = Vec::new();
        for &face in faces {
            #[cfg(feature = "work-span")]
            let _span_task = span_tasks.task();
            let planned = ahead.take_if(
                index,
                dec,
                || {
                    // Each worker owns its context and keeps any refusal proof there.
                    let (base, faces, boundary, sections, context) = (
                        std::sync::Arc::new(model.clone()),
                        every_face.clone(),
                        boundary.clone(),
                        sections.clone(),
                        context.clone(),
                    );
                    move |index: usize, dec: &mut Decider| {
                        let face = faces[index];
                        let mut context = context.clone();
                        crate::geometry::plan(&base, |model| {
                            pieces(
                                model,
                                face,
                                &boundary,
                                sections.get(&face).map(Vec::as_slice).unwrap_or(&[]),
                                dec,
                                &mut context,
                            )
                        })
                    }
                },
                // A refusal is the writer's to give, from its own run.
                |planned| planned.value().is_ok() && planned.holds(model),
            );
            index += 1;
            let split = match planned {
                Some(planned) => planned.commit(model),
                None => pieces(
                    model,
                    face,
                    &boundary,
                    sections.get(&face).map(Vec::as_slice).unwrap_or(&[]),
                    dec,
                    context,
                ),
            }?;
            visit(face, &split);
            group_pieces.extend(split);
        }
        all.push(group_pieces);
    }
    let mut all = all.into_iter();
    Ok(Answer::exact(Imprinted {
        a: graph.a,
        b: graph.b,
        first: all.next().unwrap_or_default(),
        second: all.next().unwrap_or_default(),
        common,
        singular: graph.singular.clone(),
        marks: graph.marks.clone(),
    }))
}
