//! The model: write-once arenas of topology over a geometry store.

use im::{OrdMap, Vector};
use std::collections::{BTreeMap, BTreeSet};

use catallax_recipe::{EdgeName, Origin, VertexName};
use catallax_shapes::{Arena, CurveKind, GeomStore, SurfaceKind};

use crate::{
    Body, BodyId, Coedge, CoedgeId, Edge, EdgeId, Face, FaceId, Loop, LoopId, Lump, LumpId, Sense,
    Shell, ShellId, Vertex, VertexId,
};

mod journal;
pub(crate) use journal::TopoRead;
pub use journal::{Journal, Reads, Relocate, Remap};
use journal::{TopoJournal, TopoOp};

/// Entity counts of one body.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct EntityCounts {
    /// Lumps.
    pub lumps: usize,
    /// Shells.
    pub shells: usize,
    /// Faces.
    pub faces: usize,
    /// Loops.
    pub loops: usize,
    /// Edges, ring edges included.
    pub edges: usize,
    /// Ring edges.
    pub ring_edges: usize,
    /// Vertices.
    pub vertices: usize,
    /// Faces by surface kind.
    pub faces_by_kind: BTreeMap<SurfaceKind, usize>,
    /// Edges by curve kind.
    pub edges_by_kind: BTreeMap<CurveKind, usize>,
}

/// One topological fault found by `Model::check`.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct TopologyFault {
    /// Stable fault code.
    pub code: &'static str,
    /// Free text for people.
    pub detail: String,
}

/// Topology over a geometry store. Write-once.
#[derive(Clone, Debug, Default)]
pub struct Model {
    geom: GeomStore,
    vertices: Arena<Vertex>,
    edges: Arena<Edge>,
    coedges: Arena<Coedge>,
    loops: Arena<Loop>,
    faces: Arena<Face>,
    shells: Arena<Shell>,
    lumps: Arena<Lump>,
    bodies: Arena<Body>,
    edge_coedges: Vector<Vec<CoedgeId>>,
    coedge_owners: Vector<Option<(LoopId, FaceId)>>,
    face_origins: OrdMap<FaceId, Vec<Origin>>,
    edge_constructions: OrdMap<EdgeId, Vec<EdgeName>>,
    blend_bases: OrdMap<BodyId, BodyId>,
    pub(crate) chart_expansions: OrdMap<BodyId, std::sync::Arc<crate::ChartExpansion>>,
    pub(crate) face_expansions: OrdMap<FaceId, std::sync::Arc<[FaceId]>>,
    // What has been added since `journal_begin`, while a journal is kept.
    journal: Option<TopoJournal>,
    // What this model and its copies have been asked since `watch`.
    watch: catallax_numbers::Watch<TopoRead>,
    // Work-span observation only: which copy of the model last changed each
    // side-table entry, tagged by table in the high bits.
    #[cfg(feature = "work-span")]
    span: catallax_numbers::work_span::Entries,
}

impl Model {
    /// An empty model with an empty geometry store.
    pub fn new() -> Self {
        Self::default()
    }

    /// Exact numerical chart topology for a coalesced logical body. This cover
    /// is not the body's selectable face, edge or vertex inventory.
    pub fn chart_expansion(&self, body: BodyId) -> Option<&crate::ChartExpansion> {
        self.body(body);
        self.read(|| TopoRead::ChartExpansion(body));
        self.chart_expansions.get(&body).map(std::sync::Arc::as_ref)
    }

    /// Exact chart faces covering one published logical face. These are
    /// integration/exchange cells, not additional selectable model faces.
    pub fn face_chart_cover(&self, face: FaceId) -> Option<&[FaceId]> {
        self.face(face);
        self.read(|| TopoRead::FaceCover(face));
        self.face_expansions.get(&face).map(std::sync::Arc::as_ref)
    }

    /// The geometry store.
    pub fn geom(&self) -> &GeomStore {
        &self.geom
    }

    /// The geometry store, for adding geometry and facts.
    pub fn geom_mut(&mut self) -> &mut GeomStore {
        &mut self.geom
    }

    // Writes the change into the journal, while one is kept.
    fn note(&mut self, op: impl FnOnce() -> TopoOp) {
        if let Some(journal) = &mut self.journal {
            journal.ops.push(op());
        }
    }

    // Writes into the list a read of something later additions can change,
    // while a list is kept.
    pub(crate) fn read(&self, read: impl FnOnce() -> TopoRead) {
        self.watch.note(read);
    }

    /// A change the journal does not describe. None of the work that runs on
    /// a journalled copy makes one: these belong to building and moving whole
    /// bodies.
    pub(crate) fn not_journaled(&self, change: &str) {
        assert!(
            self.journal.is_none(),
            "{change} cannot be replayed from a journal"
        );
    }

    /// Adds a vertex.
    pub fn add_vertex(&mut self, v: Vertex) -> VertexId {
        self.note(|| TopoOp::Vertex(v.clone()));
        self.geom.point(v.point);
        self.vertices.push(v)
    }

    /// Adds an edge.
    pub fn add_edge(&mut self, e: Edge) -> EdgeId {
        self.note(|| TopoOp::Edge(e.clone()));
        self.geom.curve(e.curve);
        for v in [e.start, e.end].into_iter().flatten() {
            self.vertex(v);
        }
        let id = self.edges.push(e);
        self.edge_coedges.push_back(Vec::new());
        id
    }

    /// Adds a coedge.
    pub fn add_coedge(&mut self, c: Coedge) -> CoedgeId {
        self.note(|| TopoOp::Coedge(c.clone()));
        self.edge(c.edge);
        let edge = c.edge;
        let id = self.coedges.push(c);
        self.edge_coedges[edge.index() as usize].push(id);
        #[cfg(feature = "work-span")]
        self.span.write("edge coedges", u64::from(edge.index()));
        self.coedge_owners.push_back(None);
        id
    }

    /// Adds a loop.
    pub fn add_loop(&mut self, l: Loop) -> LoopId {
        self.note(|| TopoOp::Loop(l.clone()));
        match &l {
            Loop::Coedges(cs) => {
                for c in cs {
                    self.coedge(*c);
                }
            }
            Loop::Vertex(v) => {
                self.vertex(*v);
            }
        }
        self.loops.push(l)
    }

    /// Adds a face and records which loops and coedges it owns.
    pub fn add_face(&mut self, f: Face) -> FaceId {
        self.note(|| TopoOp::Face(f.clone()));
        self.push_face(f)
    }

    fn push_face(&mut self, f: Face) -> FaceId {
        self.geom.surface(f.surface);
        for l in &f.loops {
            self.loop_(*l);
        }
        let id = self.faces.push(f);
        for l in &self.faces.get(id).loops {
            if let Loop::Coedges(cs) = self.loops.get(*l) {
                for c in cs {
                    self.coedge_owners[c.index() as usize].get_or_insert((*l, id));
                    #[cfg(feature = "work-span")]
                    self.span
                        .write("coedge owner", (1 << 40) | u64::from(c.index()));
                }
            }
        }
        id
    }

    /// Adds a face with its canonical, duplicate-free set of recipe origins.
    /// Unnamed exchange geometry continues to use [`Model::add_face`].
    pub fn add_named_face(&mut self, face: Face, mut origins: Vec<Origin>) -> FaceId {
        origins.sort();
        origins.dedup();
        self.note(|| TopoOp::NamedFace(face.clone(), origins.clone()));
        let id = self.push_face(face);
        if !origins.is_empty() {
            self.face_origins.insert(id, origins);
            #[cfg(feature = "work-span")]
            self.span
                .write("face origins", (2 << 40) | u64::from(id.index()));
        }
        id
    }

    /// Origins of a face in compact canonical JSON order.
    pub fn face_origins(&self, face: FaceId) -> &[Origin] {
        self.face(face);
        #[cfg(feature = "work-span")]
        self.span
            .read("face origins", (2 << 40) | u64::from(face.index()));
        self.read(|| TopoRead::FaceOrigins(face));
        self.face_origins.get(&face).map_or(&[], Vec::as_slice)
    }

    /// Rewrites a body's face and edge provenance without changing its geometry or incidence.
    /// Call before publishing a newly built body. Shared faces retain this rewrite.
    pub fn map_face_origins(&mut self, body: BodyId, mut map: impl FnMut(&Origin) -> Origin) {
        self.not_journaled("a rewrite of face origins");
        for face in self.faces_of_body(body) {
            if let Some(origins) = self.face_origins.get_mut(&face) {
                #[cfg(feature = "work-span")]
                self.span
                    .write("face origins", (2 << 40) | u64::from(face.index()));
                for origin in origins.iter_mut() {
                    *origin = map(origin);
                }
                origins.sort();
                origins.dedup();
            }
        }
        if self.edge_constructions.is_empty() {
            return;
        }
        for edge in self.edges_of_body(body) {
            if let Some(names) = self.edge_constructions.remove(&edge) {
                #[cfg(feature = "work-span")]
                self.span
                    .write("edge constructions", (3 << 40) | u64::from(edge.index()));
                let names = names
                    .into_iter()
                    .filter_map(|name| {
                        let [a, b] = name.between();
                        EdgeName::new(map(a), map(b)).ok()
                    })
                    .collect();
                self.record_edge_constructions(edge, names);
            }
        }
    }

    /// Resolves an origin on the result, including every split piece and a
    /// merged face whose origin set contains it.
    pub fn faces_with_origin(&self, body: BodyId, origin: &Origin) -> Vec<FaceId> {
        self.unique_faces(body)
            .filter(|f| self.face_origins(*f).contains(origin))
            .collect()
    }

    /// Resolves a saved edge name against the complete incident origin sets.
    /// A merge can change a face's first origin without losing this reference.
    pub fn edges_with_name(&self, body: BodyId, name: &EdgeName) -> Vec<EdgeId> {
        let [a, b] = name.between();
        self.edges_of_body(body).into_iter().filter(|edge| {
            let mut faces = self.unique_faces(body).filter(|f| {
                self.face(*f).loops.iter().any(|lp| {
                    matches!(self.loop_(*lp), Loop::Coedges(cs) if cs.iter().any(|c| self.coedge(*c).edge == *edge))
                })
            });
            let (Some(first), Some(second)) = (faces.next(), faces.next()) else { return false; };
            if faces.next().is_some() { return false; }
            let first = self.face_origins(first);
            let second = self.face_origins(second);
            (first.contains(a) && second.contains(b)) || (first.contains(b) && second.contains(a))
        }).collect()
    }

    /// Construction names retained through edge splitting and maximal merging.
    /// Unlike a merged face's origins, these pairs were incident before merging.
    pub fn edge_constructions(&self, edge: EdgeId) -> &[EdgeName] {
        self.edge(edge);
        #[cfg(feature = "work-span")]
        self.span
            .read("edge constructions", (3 << 40) | u64::from(edge.index()));
        self.read(|| TopoRead::EdgeConstructions(edge));
        self.edge_constructions
            .get(&edge)
            .map_or(&[], Vec::as_slice)
    }

    /// Records the canonical union of construction names for a newly built edge.
    pub fn record_edge_constructions(&mut self, edge: EdgeId, mut names: Vec<EdgeName>) {
        self.note(|| TopoOp::EdgeConstructions(edge, names.clone()));
        self.edge(edge);
        if names.is_empty() {
            return;
        }
        #[cfg(feature = "work-span")]
        self.span
            .write("edge constructions", (3 << 40) | u64::from(edge.index()));
        let stored = self.edge_constructions.entry(edge).or_default();
        stored.append(&mut names);
        stored.sort_by(|a, b| a.between().cmp(b.between()));
        stored.dedup();
    }

    /// Names an edge by a retained construction pair still incident on this body.
    /// Without construction history, uses the first origins of its two faces.
    pub fn edge_name(&self, body: BodyId, edge: EdgeId) -> Option<EdgeName> {
        let mut faces = self.unique_faces(body).filter(|face| {
            self.face(*face).loops.iter().any(|lp| {
                matches!(self.loop_(*lp), Loop::Coedges(cs) if cs.iter().any(|c| self.coedge(*c).edge == edge))
            })
        });
        let (Some(a), Some(b)) = (faces.next(), faces.next()) else {
            return None;
        };
        if faces.next().is_some() {
            return None;
        }
        let a = self.face_origins(a);
        let b = self.face_origins(b);
        if let Some(name) = self.edge_constructions(edge).iter().find(|name| {
            let [x, y] = name.between();
            (a.contains(x) && b.contains(y)) || (a.contains(y) && b.contains(x))
        }) {
            return Some(name.clone());
        }
        EdgeName::new(a.first()?.clone(), b.first()?.clone()).ok()
    }

    /// Names a vertex by the first origins of all incident faces.
    pub fn vertex_name(&self, body: BodyId, vertex: VertexId) -> Option<VertexName> {
        let origins = self
            .unique_faces(body)
            .filter_map(|face| {
                self.face(face)
                    .loops
                    .iter()
                    .any(|lp| match self.loop_(*lp) {
                        Loop::Vertex(v) => *v == vertex,
                        Loop::Coedges(cs) => cs.iter().any(|c| {
                            let edge = self.edge(self.coedge(*c).edge);
                            edge.start == Some(vertex) || edge.end == Some(vertex)
                        }),
                    })
                    .then(|| self.face_origins(face).first().cloned())
                    .flatten()
            })
            .collect();
        VertexName::new(origins).ok()
    }

    /// Records the unfilleted body a fitted blend changed. Measurement reads a
    /// fitted junction as this body plus the blend's exact ring.
    pub fn record_blend_base(&mut self, blend: BodyId, base: BodyId) {
        self.note(|| TopoOp::BlendBase(blend, base));
        self.body(blend);
        self.body(base);
        self.blend_bases.insert(blend, base);
        #[cfg(feature = "work-span")]
        self.span
            .write("blend base", (4 << 40) | u64::from(blend.index()));
    }

    /// The unfilleted body behind a fitted blend.
    pub fn blend_base(&self, blend: BodyId) -> Option<BodyId> {
        #[cfg(feature = "work-span")]
        self.span
            .read("blend base", (4 << 40) | u64::from(blend.index()));
        self.read(|| TopoRead::BlendBase(blend));
        self.blend_bases.get(&blend).copied()
    }

    /// Adds a shell.
    pub fn add_shell(&mut self, s: Shell) -> ShellId {
        self.note(|| TopoOp::Shell(s.clone()));
        for f in &s.faces {
            self.face(*f);
        }
        self.shells.push(s)
    }

    /// Adds a lump.
    pub fn add_lump(&mut self, l: Lump) -> LumpId {
        self.note(|| TopoOp::Lump(l.clone()));
        for s in &l.shells {
            self.shell(*s);
        }
        self.lumps.push(l)
    }

    /// Adds a body.
    pub fn add_body(&mut self, b: Body) -> BodyId {
        self.note(|| TopoOp::Body(b.clone()));
        for l in &b.lumps {
            self.lump(*l);
        }
        self.bodies.push(b)
    }

    /// The vertex behind a handle.
    pub fn vertex(&self, id: VertexId) -> &Vertex {
        self.vertices.get(id)
    }

    /// The edge behind a handle.
    pub fn edge(&self, id: EdgeId) -> &Edge {
        self.edges.get(id)
    }

    /// The coedge behind a handle.
    pub fn coedge(&self, id: CoedgeId) -> &Coedge {
        self.coedges.get(id)
    }

    /// The loop behind a handle.
    pub fn loop_(&self, id: LoopId) -> &Loop {
        self.loops.get(id)
    }

    /// The face behind a handle.
    pub fn face(&self, id: FaceId) -> &Face {
        self.faces.get(id)
    }

    /// The shell behind a handle.
    pub fn shell(&self, id: ShellId) -> &Shell {
        self.shells.get(id)
    }

    /// The lump behind a handle.
    pub fn lump(&self, id: LumpId) -> &Lump {
        self.lumps.get(id)
    }

    /// The body behind a handle.
    pub fn body(&self, id: BodyId) -> &Body {
        self.bodies.get(id)
    }

    /// Every coedge that uses `e`, in creation order.
    pub fn coedges_of_edge(&self, e: EdgeId) -> Vec<CoedgeId> {
        self.edge(e);
        #[cfg(feature = "work-span")]
        self.span.read("edge coedges", u64::from(e.index()));
        self.read(|| TopoRead::EdgeCoedges(e));
        self.edge_coedges[e.index() as usize].clone()
    }

    /// The loop and face that own a coedge, once its face has been added.
    pub fn owner_of_coedge(&self, c: CoedgeId) -> Option<(LoopId, FaceId)> {
        self.coedge(c);
        #[cfg(feature = "work-span")]
        self.span
            .read("coedge owner", (1 << 40) | u64::from(c.index()));
        self.read(|| TopoRead::CoedgeOwner(c));
        self.coedge_owners[c.index() as usize]
    }

    fn body_faces(&self, b: BodyId) -> impl Iterator<Item = FaceId> + '_ {
        self.body(b).lumps.iter().flat_map(|l| {
            self.lump(*l)
                .shells
                .iter()
                .flat_map(|s| self.shell(*s).faces.iter().copied())
        })
    }

    fn unique_faces(&self, b: BodyId) -> impl Iterator<Item = FaceId> + '_ {
        let mut seen = BTreeSet::new();
        self.body_faces(b).filter(move |f| seen.insert(*f))
    }

    /// Every body with this face, in creation order. No face records its body.
    pub fn bodies_with_face(&self, f: FaceId) -> Vec<BodyId> {
        self.read(|| TopoRead::Bodies);
        self.bodies
            .handles()
            .into_iter()
            .filter(|&b| self.body_faces(b).any(|face| face == f))
            .collect()
    }

    /// Faces of a body in the fixed order of the piece description.
    pub fn faces_of_body(&self, b: BodyId) -> Vec<FaceId> {
        self.unique_faces(b).collect()
    }

    /// Edges of a body in the fixed order of the piece description.
    pub fn edges_of_body(&self, b: BodyId) -> Vec<EdgeId> {
        let mut seen = BTreeSet::new();
        let mut result = Vec::new();
        for f in self.unique_faces(b) {
            for l in &self.face(f).loops {
                if let Loop::Coedges(cs) = self.loop_(*l) {
                    for c in cs {
                        let e = self.coedge(*c).edge;
                        if seen.insert(e) {
                            result.push(e);
                        }
                    }
                }
            }
        }
        result
    }

    /// Vertices of a body in the fixed order of the piece description.
    pub fn vertices_of_body(&self, b: BodyId) -> Vec<VertexId> {
        let mut seen = BTreeSet::new();
        let mut result = Vec::new();
        let mut visit = |v| {
            if seen.insert(v) {
                result.push(v);
            }
        };
        for f in self.unique_faces(b) {
            for l in &self.face(f).loops {
                match self.loop_(*l) {
                    Loop::Vertex(v) => visit(*v),
                    Loop::Coedges(cs) => {
                        for c in cs {
                            let (start, end) = self.coedge_ends(*c);
                            for v in [start, end].into_iter().flatten() {
                                visit(v);
                            }
                        }
                    }
                }
            }
        }
        result
    }

    fn coedge_ends(&self, c: CoedgeId) -> (Option<VertexId>, Option<VertexId>) {
        let c = self.coedge(c);
        let e = self.edge(c.edge);
        match c.sense {
            Sense::Same => (e.start, e.end),
            Sense::Reversed => (e.end, e.start),
        }
    }

    /// Entity counts of a body.
    pub fn counts(&self, b: BodyId) -> EntityCounts {
        let body = self.body(b);
        let mut edges = BTreeSet::new();
        let mut vertices = BTreeSet::new();
        let mut counts = EntityCounts {
            lumps: body.lumps.len(),
            shells: body.lumps.iter().map(|l| self.lump(*l).shells.len()).sum(),
            ..EntityCounts::default()
        };
        for f in self.unique_faces(b) {
            let face = self.face(f);
            counts.faces += 1;
            counts.loops += face.loops.len();
            *counts
                .faces_by_kind
                .entry(self.geom.surface(face.surface).kind())
                .or_default() += 1;
            for l in &face.loops {
                match self.loop_(*l) {
                    Loop::Vertex(v) => {
                        vertices.insert(*v);
                    }
                    Loop::Coedges(cs) => {
                        for c in cs {
                            let e = self.coedge(*c).edge;
                            if edges.insert(e) {
                                let edge = self.edge(e);
                                counts.ring_edges +=
                                    usize::from(edge.start.is_none() && edge.end.is_none());
                                *counts
                                    .edges_by_kind
                                    .entry(self.geom.curve(edge.curve).kind())
                                    .or_default() += 1;
                                vertices.extend([edge.start, edge.end].into_iter().flatten());
                            }
                        }
                    }
                }
            }
        }
        counts.edges = edges.len();
        counts.vertices = vertices.len();
        counts
    }

    /// Topological faults of a body; empty when sound. Reads entities and
    /// recorded facts only and makes no numeric decision.
    pub fn check(&self, b: BodyId) -> Vec<TopologyFault> {
        crate::check::check(self, b)
    }
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
