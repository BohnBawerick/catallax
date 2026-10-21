//! The model: write-once arenas of topology over a geometry store.

use std::collections::{BTreeMap, BTreeSet};

use catallax_shapes::{Arena, CurveKind, GeomStore, SurfaceKind};

use crate::{
    Body, BodyId, Coedge, CoedgeId, Edge, EdgeId, Face, FaceId, Loop, LoopId, Lump, LumpId, Sense,
    Shell, ShellId, Vertex, VertexId,
};

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
    edge_coedges: Vec<Vec<CoedgeId>>,
    coedge_owners: Vec<Option<(LoopId, FaceId)>>,
}

impl Model {
    /// An empty model with an empty geometry store.
    pub fn new() -> Self {
        Self::default()
    }

    /// The geometry store.
    pub fn geom(&self) -> &GeomStore {
        &self.geom
    }

    /// The geometry store, for adding geometry and facts.
    pub fn geom_mut(&mut self) -> &mut GeomStore {
        &mut self.geom
    }

    /// Adds a vertex.
    pub fn add_vertex(&mut self, v: Vertex) -> VertexId {
        self.geom.point(v.point);
        self.vertices.push(v)
    }

    /// Adds an edge.
    pub fn add_edge(&mut self, e: Edge) -> EdgeId {
        self.geom.curve(e.curve);
        for v in [e.start, e.end].into_iter().flatten() {
            self.vertex(v);
        }
        let id = self.edges.push(e);
        self.edge_coedges.push(Vec::new());
        id
    }

    /// Adds a coedge.
    pub fn add_coedge(&mut self, c: Coedge) -> CoedgeId {
        self.edge(c.edge);
        let edge = c.edge;
        let id = self.coedges.push(c);
        self.edge_coedges[edge.index() as usize].push(id);
        self.coedge_owners.push(None);
        id
    }

    /// Adds a loop.
    pub fn add_loop(&mut self, l: Loop) -> LoopId {
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
        self.geom.surface(f.surface);
        for l in &f.loops {
            self.loop_(*l);
        }
        let id = self.faces.push(f);
        for l in &self.faces.get(id).loops {
            if let Loop::Coedges(cs) = self.loops.get(*l) {
                for c in cs {
                    self.coedge_owners[c.index() as usize].get_or_insert((*l, id));
                }
            }
        }
        id
    }

    /// Adds a shell.
    pub fn add_shell(&mut self, s: Shell) -> ShellId {
        for f in &s.faces {
            self.face(*f);
        }
        self.shells.push(s)
    }

    /// Adds a lump.
    pub fn add_lump(&mut self, l: Lump) -> LumpId {
        for s in &l.shells {
            self.shell(*s);
        }
        self.lumps.push(l)
    }

    /// Adds a body.
    pub fn add_body(&mut self, b: Body) -> BodyId {
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
        self.edge_coedges[e.index() as usize].clone()
    }

    /// The loop and face that own a coedge, once its face has been added.
    pub fn owner_of_coedge(&self, c: CoedgeId) -> Option<(LoopId, FaceId)> {
        self.coedge(c);
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
