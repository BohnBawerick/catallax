//! The one checker of a body, whether its entities are exact, tolerant or
//! both (section 4.5 of the Checkpoint 2 decisions). It reads what builders
//! recorded and measures nothing.
//!
//! Structure comes first: kinds, sharing, loop closure, edge uses,
//! connectivity and the Euler count. None of it depends on a tolerance.
//!
//! Then each pair that must meet is asked for its record. The pairs are a
//! vertex on its edge's curve, an edge's curve on its face's surface and a
//! vertex loop's point on its face's surface. Two exact entities need their
//! `Incident` fact in the fact book, which is never recomputed here. A
//! tolerant entity has its tolerance as the record. Its builder measured the
//! distance and stored the result on the entity (rule T4), as an exact
//! builder records its fact.
//!
//! Last, the order. A vertex's tolerance is at least that of every edge that
//! ends at it, and an edge's at least that of every face it bounds, with
//! exact counting as zero. Parasolid keeps this order (Jackson 1995,
//! "Boundary representation modelling with local tolerances") and so does
//! OpenCASCADE (`BRepLib::UpdateTolerances`). A tolerant edge therefore
//! forces its vertices up and never its faces, and an exact face beside it
//! stays exact.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::FactKind;
use catallax_shapes::GeomKey;

use crate::{
    BodyId, BodyKind, CoedgeId, EdgeId, FaceId, Loop, LoopId, LumpId, Model, Sense, ShellId,
    Tolerance, TopologyFault, VertexId,
};

pub(super) fn check(model: &Model, body: BodyId) -> Vec<TopologyFault> {
    check_with_incidence(model, body, true)
}

/// The structure part of [`check`] alone. Numerical admission runs it on a
/// body whose tolerances it has not measured yet, so that a broken shell
/// refuses before any distance is taken; the admitted body then goes through
/// [`check`] like every other.
pub(super) fn structure(model: &Model, body: BodyId) -> Vec<TopologyFault> {
    check_with_incidence(model, body, false)
}

fn check_with_incidence(model: &Model, body: BodyId, exact_incidence: bool) -> Vec<TopologyFault> {
    let mut check = Checker {
        exact_incidence,
        surface: model.body(body).kind == BodyKind::Surface,
        model,
        faults: Vec::new(),
        lumps: BTreeSet::new(),
        shells: BTreeSet::new(),
        faces: BTreeSet::new(),
        loops: BTreeSet::new(),
        coedges: BTreeSet::new(),
        edge_shells: BTreeMap::new(),
        vertex_shells: BTreeMap::new(),
    };
    let b = model.body(body);
    if (matches!(b.kind, BodyKind::Solid | BodyKind::Surface)
        && (b.lumps.is_empty() || !b.wires.is_empty()))
        || (b.kind == BodyKind::Empty && (!b.lumps.is_empty() || !b.wires.is_empty()))
        || (b.kind == BodyKind::Wire && (!b.lumps.is_empty() || b.wires.is_empty()))
    {
        check.fault("topology.body-kind", || {
            format!("body {body:?} has incompatible lumps")
        });
    }
    let mut wire_edges = BTreeSet::new();
    for wire in &b.wires {
        if wire.coedges.is_empty() {
            check.fault("topology.empty-part", || "wire has no coedge".into());
        }
        for (index, coedge) in wire.coedges.iter().enumerate() {
            if !check.coedges.insert(*coedge) {
                check.shared(coedge);
            }
            let edge = model.coedge(*coedge).edge;
            if !wire_edges.insert(edge) {
                check.shared(&edge);
            }
            check.edge_geometry(edge);
            let (start, end) = model.coedge_ends(*coedge);
            if (start.is_none() && end.is_none() && wire.coedges.len() != 1)
                || (index > 0 && model.coedge_ends(wire.coedges[index - 1]).1 != start)
            {
                check.fault("topology.wire-not-connected", || {
                    format!("wire is broken at coedge {coedge:?}")
                });
            }
        }
    }
    for l in &b.lumps {
        if !check.lumps.insert(*l) {
            check.shared(l);
        }
        let lump = model.lump(*l);
        if lump.shells.is_empty() {
            check.fault("topology.empty-part", || format!("lump {l:?} has no shell"));
        }
        for s in &lump.shells {
            if !check.shells.insert(*s) {
                check.shared(s);
            }
            check.shell(*s);
        }
    }
    check.faults
}

struct Checker<'a> {
    exact_incidence: bool,
    surface: bool,
    model: &'a Model,
    faults: Vec<TopologyFault>,
    lumps: BTreeSet<LumpId>,
    shells: BTreeSet<ShellId>,
    faces: BTreeSet<FaceId>,
    loops: BTreeSet<LoopId>,
    coedges: BTreeSet<CoedgeId>,
    edge_shells: BTreeMap<EdgeId, ShellId>,
    vertex_shells: BTreeMap<VertexId, ShellId>,
}

#[derive(Default)]
struct EdgeUse {
    same: usize,
    reversed: usize,
    first_face: usize,
}

struct ShellScan {
    edges: BTreeMap<EdgeId, EdgeUse>,
    vertices: BTreeSet<VertexId>,
    // Union by size and path compression keep connectivity work near-linear,
    // even for an invalid edge used by many faces.
    parents: Vec<usize>,
    sizes: Vec<usize>,
    face_chi: i128,
    rings: usize,
}

impl ShellScan {
    fn root(&mut self, mut f: usize) -> usize {
        while self.parents[f] != f {
            self.parents[f] = self.parents[self.parents[f]];
            f = self.parents[f];
        }
        f
    }

    fn join(&mut self, a: usize, b: usize) {
        let (mut a, mut b) = (self.root(a), self.root(b));
        if a != b {
            if self.sizes[a] < self.sizes[b] {
                std::mem::swap(&mut a, &mut b);
            }
            self.parents[b] = a;
            self.sizes[a] += self.sizes[b];
        }
    }
}

impl Checker<'_> {
    fn fault(&mut self, code: &'static str, detail: impl FnOnce() -> String) {
        self.faults.push(TopologyFault {
            code,
            detail: detail(),
        });
    }

    fn shared(&mut self, id: &impl std::fmt::Debug) {
        self.fault("topology.shared-entity", || {
            format!("entity {id:?} is shared within the body")
        });
    }

    /// The incidence of `a` on `b`, where `a` belongs to the entity with
    /// tolerance `inner` and `b` to the one with tolerance `outer`: a vertex
    /// on an edge, an edge on a face, a vertex on a face.
    ///
    /// Two exact entities need the recorded fact. A tolerant one has its
    /// tolerance as the record. The inner entity's tolerance may not be less
    /// than the outer's.
    fn incidence(
        &mut self,
        inner: (&dyn std::fmt::Debug, Option<&Tolerance>, GeomKey),
        outer: (&dyn std::fmt::Debug, Option<&Tolerance>, GeomKey),
    ) {
        let (inner_id, inner_tolerance, a) = inner;
        let (outer_id, outer_tolerance, b) = outer;
        if inner_tolerance < outer_tolerance {
            self.fault("topology.tolerance-order", || {
                format!(
                    "{inner_id:?} has {} and bounds {outer_id:?}, which has {}",
                    stated(inner_tolerance),
                    stated(outer_tolerance)
                )
            });
        }
        if self.exact_incidence
            && inner_tolerance.is_none()
            && outer_tolerance.is_none()
            && !self.model.geom().facts().holds(FactKind::Incident, a, b)
        {
            self.fault("topology.missing-incidence", || {
                format!("missing Incident({a:?}, {b:?})")
            });
        }
    }

    fn vertex(&mut self, v: VertexId, shell: ShellId, scan: &mut ShellScan) {
        if scan.vertices.insert(v) {
            if let Some(previous) = self.vertex_shells.get(&v) {
                if *previous != shell {
                    self.shared(&v);
                }
            } else {
                self.vertex_shells.insert(v, shell);
            }
        }
    }

    fn shell(&mut self, s: ShellId) {
        let m = self.model;
        let shell = m.shell(s);
        if shell.faces.is_empty() {
            self.fault("topology.empty-part", || format!("shell {s:?} has no face"));
        }
        let mut scan = ShellScan {
            edges: BTreeMap::new(),
            vertices: BTreeSet::new(),
            parents: (0..shell.faces.len()).collect(),
            sizes: vec![1; shell.faces.len()],
            face_chi: 0,
            rings: 0,
        };
        for (index, f) in shell.faces.iter().copied().enumerate() {
            if !self.faces.insert(f) {
                self.shared(&f);
            }
            self.face(f, index, s, &mut scan);
        }
        for (e, uses) in &scan.edges {
            let boundary = self.surface && uses.same + uses.reversed == 1;
            if !boundary && (uses.same != 1 || uses.reversed != 1) {
                self.fault("topology.edge-use", || {
                    format!(
                        "edge {e:?} in shell {s:?} has {} same and {} reversed uses",
                        uses.same, uses.reversed
                    )
                });
            }
        }
        if !shell.faces.is_empty() {
            let root = scan.root(0);
            if (1..shell.faces.len()).any(|f| scan.root(f) != root) {
                self.fault("topology.disconnected-shell", || {
                    format!("shell {s:?} is not edge-connected")
                });
            }
        }
        let chi = scan.vertices.len() as i128 + scan.rings as i128 - scan.edges.len() as i128
            + scan.face_chi;
        let boundaries = if self.surface {
            self.surface_boundaries(s, &scan)
        } else {
            0
        };
        if chi + boundaries as i128 > 2 || (chi + boundaries as i128) % 2 != 0 {
            self.fault("topology.euler", || {
                format!("shell {s:?} has Euler characteristic {chi}")
            });
        }
    }

    fn surface_boundaries(&mut self, shell: ShellId, scan: &ShellScan) -> usize {
        let mut boundary: BTreeMap<VertexId, Vec<VertexId>> = BTreeMap::new();
        let mut rings = 0;
        for (edge, uses) in &scan.edges {
            let edge = self.model.edge(*edge);
            if uses.same + uses.reversed != 1 {
                continue;
            }
            match (edge.start, edge.end) {
                (Some(a), Some(b)) => {
                    boundary.entry(a).or_default().push(b);
                    boundary.entry(b).or_default().push(a);
                }
                (None, None) => rings += 1,
                _ => {}
            }
        }
        for (vertex, neighbors) in &boundary {
            if neighbors.len() != 2 {
                self.fault("topology.boundary-vertex", || {
                    format!(
                        "vertex {vertex:?} in shell {shell:?} has {} boundary uses",
                        neighbors.len()
                    )
                });
            }
        }
        let mut seen = BTreeSet::new();
        for vertex in boundary.keys() {
            if seen.contains(vertex) {
                continue;
            }
            rings += 1;
            let mut pending = vec![*vertex];
            while let Some(vertex) = pending.pop() {
                if seen.insert(vertex) {
                    pending.extend(boundary[&vertex].iter().copied());
                }
            }
        }
        self.surface_fans(shell);
        rings
    }

    fn surface_fans(&mut self, shell: ShellId) {
        let m = self.model;
        let mut fans: BTreeMap<VertexId, BTreeMap<CoedgeId, BTreeSet<CoedgeId>>> = BTreeMap::new();
        let mut endpoints: BTreeMap<(EdgeId, bool, VertexId), Vec<CoedgeId>> = BTreeMap::new();
        for face in &m.shell(shell).faces {
            for boundary in &m.face(*face).loops {
                let Loop::Coedges(coedges) = m.loop_(*boundary) else {
                    continue;
                };
                for (index, corner) in coedges.iter().copied().enumerate() {
                    let outgoing = m.coedge(corner);
                    let edge = m.edge(outgoing.edge);
                    let start = outgoing.sense == Sense::Same;
                    let Some(vertex) = (if start { edge.start } else { edge.end }) else {
                        continue;
                    };
                    fans.entry(vertex).or_default().entry(corner).or_default();
                    endpoints
                        .entry((outgoing.edge, start, vertex))
                        .or_default()
                        .push(corner);
                    let incoming = m.coedge(coedges[(index + coedges.len() - 1) % coedges.len()]);
                    endpoints
                        .entry((incoming.edge, incoming.sense == Sense::Reversed, vertex))
                        .or_default()
                        .push(corner);
                }
            }
        }
        // A node is a corner occurrence, so two wedges of one face stay distinct.
        for ((_, _, vertex), corners) in endpoints {
            if let [a, b] = corners.as_slice() {
                let fan = fans.get_mut(&vertex).expect("endpoint has a corner");
                fan.get_mut(a).expect("corner was registered").insert(*b);
                fan.get_mut(b).expect("corner was registered").insert(*a);
            }
        }
        for (vertex, fan) in fans {
            let mut seen = BTreeSet::new();
            let mut pending = vec![*fan.keys().next().expect("vertex has a corner")];
            while let Some(corner) = pending.pop() {
                if seen.insert(corner) {
                    pending.extend(fan[&corner].iter().copied());
                }
            }
            if seen.len() != fan.len() {
                self.fault("topology.vertex-fan", || {
                    format!("vertex {vertex:?} in shell {shell:?} has disconnected corner fans")
                });
            }
        }
    }

    fn face(&mut self, f: FaceId, index: usize, s: ShellId, scan: &mut ShellScan) {
        let m = self.model;
        let face = m.face(f);
        let surface = m.geom().surface(face.surface);
        if face.loops.is_empty()
            && matches!(
                surface,
                catallax_shapes::Surface::Plane(_)
                    | catallax_shapes::Surface::Cylinder(_)
                    | catallax_shapes::Surface::Cone(_)
            )
        {
            self.fault("topology.unbounded-face", || {
                format!("face {f:?} has no boundary")
            });
        }
        scan.face_chi +=
            if face.loops.is_empty() && matches!(surface, catallax_shapes::Surface::Torus(_)) {
                0
            } else {
                2 - face.loops.len() as i128
            };
        let has_coedge_loop = face
            .loops
            .iter()
            .any(|l| matches!(m.loop_(*l), Loop::Coedges(_)));
        let mut vertex_loops = BTreeSet::new();
        let mut first_loop = BTreeMap::new();
        let mut multiple_loops = BTreeSet::new();
        for (li, l) in face.loops.iter().enumerate() {
            if !self.loops.insert(*l) {
                self.shared(l);
            }
            let mut visit = |v| {
                self.vertex(v, s, scan);
                if let Some(previous) = first_loop.get(&v) {
                    if *previous != li {
                        multiple_loops.insert(v);
                    }
                } else {
                    first_loop.insert(v, li);
                }
            };
            match m.loop_(*l) {
                Loop::Vertex(v) => {
                    visit(*v);
                    vertex_loops.insert(*v);
                    let vertex = m.vertex(*v);
                    self.incidence(
                        (v, vertex.tolerance.as_ref(), GeomKey::Point(vertex.point)),
                        (&f, face.tolerance.as_ref(), GeomKey::Surface(face.surface)),
                    );
                }
                Loop::Coedges(cs) => {
                    for c in cs {
                        let edge = m.edge(m.coedge(*c).edge);
                        for v in [edge.start, edge.end].into_iter().flatten() {
                            visit(v);
                        }
                    }
                    if cs.is_empty() {
                        self.fault("topology.empty-part", || {
                            format!("loop {l:?} has no coedge")
                        });
                    }
                    for (ci, c) in cs.iter().enumerate() {
                        if !self.coedges.insert(*c) {
                            self.shared(c);
                        }
                        let coedge = m.coedge(*c);
                        let edge = m.edge(coedge.edge);
                        let ring = edge.start.is_none() && edge.end.is_none();
                        let next = m.coedge(cs[(ci + 1) % cs.len()]);
                        let next_edge = m.edge(next.edge);
                        let end = if coedge.sense == Sense::Same {
                            edge.end
                        } else {
                            edge.start
                        };
                        let next_start = if next.sense == Sense::Same {
                            next_edge.start
                        } else {
                            next_edge.end
                        };
                        if (ring && cs.len() != 1) || end != next_start {
                            self.fault("topology.loop-not-closed", || {
                                format!("loop {l:?} is broken at coedge {c:?}")
                            });
                        }
                        self.incidence(
                            (
                                &coedge.edge,
                                edge.tolerance.as_ref(),
                                GeomKey::Curve(edge.curve),
                            ),
                            (&f, face.tolerance.as_ref(), GeomKey::Surface(face.surface)),
                        );
                        let first_face = if let Some(uses) = scan.edges.get_mut(&coedge.edge) {
                            match coedge.sense {
                                Sense::Same => uses.same += 1,
                                Sense::Reversed => uses.reversed += 1,
                            }
                            uses.first_face
                        } else {
                            scan.edges.insert(
                                coedge.edge,
                                EdgeUse {
                                    same: usize::from(coedge.sense == Sense::Same),
                                    reversed: usize::from(coedge.sense == Sense::Reversed),
                                    first_face: index,
                                },
                            );
                            scan.rings += usize::from(ring);
                            self.edge(coedge.edge, s);
                            index
                        };
                        scan.join(first_face, index);
                    }
                }
            }
        }
        for v in vertex_loops {
            if !has_coedge_loop || multiple_loops.contains(&v) {
                self.fault("topology.vertex-loop", || {
                    format!("vertex loop {v:?} is not an isolated boundary of face {f:?}")
                });
            }
        }
    }

    fn edge(&mut self, e: EdgeId, shell: ShellId) {
        if let Some(previous) = self.edge_shells.get(&e) {
            if *previous != shell {
                self.shared(&e);
            }
        } else {
            self.edge_shells.insert(e, shell);
        }
        self.edge_geometry(e);
    }

    fn edge_geometry(&mut self, e: EdgeId) {
        let m = self.model;
        let edge = m.edge(e);
        if edge.start.is_none() != edge.end.is_none()
            || (edge.start == edge.end && m.geom().curve(edge.curve).is_closed() == Some(false))
        {
            self.fault("topology.edge-ends", || {
                format!("edge {e:?} has invalid endpoints for its curve")
            });
        }
        for v in [edge.start, edge.end].into_iter().flatten() {
            let vertex = m.vertex(v);
            self.incidence(
                (&v, vertex.tolerance.as_ref(), GeomKey::Point(vertex.point)),
                (&e, edge.tolerance.as_ref(), GeomKey::Curve(edge.curve)),
            );
        }
    }
}

fn stated(tolerance: Option<&Tolerance>) -> String {
    match tolerance {
        Some(tolerance) => format!("tolerance {} mm", tolerance.mm()),
        None => "no tolerance".into(),
    }
}
