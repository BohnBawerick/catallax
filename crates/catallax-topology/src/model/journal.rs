//! A journal of what a piece of work adds to its own copy of a model, and its
//! replay onto the model the copy was taken from.
//!
//! A model only grows: entities are appended, facts recorded, names attached.
//! So work that ran on a copy can be done again on the original in the same
//! order, even after the original has grown, provided the handles the work
//! made are moved past what the original gained meanwhile. The replay goes
//! through the model's own adding methods, so the tables those methods keep,
//! which coedges an edge has and which face owns a coedge, come out as if the
//! work had run on the original.

use catallax_recipe::{EdgeName, Origin};
use catallax_shapes::{CurveId, GeomJournal, GeomRemap, Handle, PointId, SurfaceId};

use super::Model;
use crate::{
    Body, BodyId, Coedge, CoedgeId, Edge, EdgeId, Face, FaceId, Loop, LoopId, Lump, LumpId, Shell,
    ShellId, Vertex, VertexId,
};

/// One change to a model's topology, as its caller asked for it.
#[derive(Clone, Debug)]
pub(super) enum TopoOp {
    Vertex(Vertex),
    Edge(Edge),
    Coedge(Coedge),
    Loop(Loop),
    Face(Face),
    NamedFace(Face, Vec<Origin>),
    EdgeConstructions(EdgeId, Vec<EdgeName>),
    BlendBase(BodyId, BodyId),
    Shell(Shell),
    Lump(Lump),
    Body(Body),
}

/// The topology half of a journal while it is being kept.
#[derive(Clone, Debug)]
pub(super) struct TopoJournal {
    // Vertices, edges, coedges, loops, faces, shells, lumps and bodies the
    // copy held when the journal began.
    pub(super) base: [u32; 8],
    pub(super) ops: Vec<TopoOp>,
}

/// What was added to a copy of a model between [`Model::journal_begin`] and
/// [`Model::journal_end`], in the order it was added.
#[derive(Clone, Debug)]
pub struct Journal {
    geom: GeomJournal,
    topo: TopoJournal,
}

impl Journal {
    /// Whether nothing was added.
    pub fn is_empty(&self) -> bool {
        self.geom.is_empty() && self.topo.ops.is_empty()
    }
}

/// Where the handles a journal made on its copy landed in the model it was
/// replayed onto. Handles older than the journal are unchanged.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Remap {
    geom: GeomRemap,
    base: [u32; 8],
    shift: [u32; 8],
}

impl Remap {
    fn moved<T>(&self, arena: usize, id: Handle<T>) -> Handle<T> {
        if id.index() >= self.base[arena] {
            id.moved(self.shift[arena])
        } else {
            id
        }
    }

    /// The replayed handles of geometry.
    pub fn geom(&self) -> &GeomRemap {
        &self.geom
    }

    /// The replayed handle of a vertex.
    pub fn vertex(&self, id: VertexId) -> VertexId {
        self.moved(0, id)
    }

    /// The replayed handle of an edge.
    pub fn edge(&self, id: EdgeId) -> EdgeId {
        self.moved(1, id)
    }

    /// The replayed handle of a coedge.
    pub fn coedge(&self, id: CoedgeId) -> CoedgeId {
        self.moved(2, id)
    }

    /// The replayed handle of a loop.
    pub fn loop_(&self, id: LoopId) -> LoopId {
        self.moved(3, id)
    }

    /// The replayed handle of a face.
    pub fn face(&self, id: FaceId) -> FaceId {
        self.moved(4, id)
    }

    /// The replayed handle of a shell.
    pub fn shell(&self, id: ShellId) -> ShellId {
        self.moved(5, id)
    }

    /// The replayed handle of a lump.
    pub fn lump(&self, id: LumpId) -> LumpId {
        self.moved(6, id)
    }

    /// The replayed handle of a body.
    pub fn body(&self, id: BodyId) -> BodyId {
        self.moved(7, id)
    }

    fn face_value(&self, face: &Face) -> Face {
        Face {
            surface: self.geom.surface(face.surface),
            loops: face.loops.iter().map(|l| self.loop_(*l)).collect(),
            ..face.clone()
        }
    }
}

/// A value that may hold handles a journal made, moved to where the replay
/// put them. Work that ran on a copy of a model returns its result through
/// this.
pub trait Relocate {
    /// The value with every handle the journal made replaced by its replayed
    /// handle.
    fn relocate(self, remap: &Remap) -> Self;
}

macro_rules! plain {
    ($($kind:ty),*) => {$(
        impl Relocate for $kind {
            fn relocate(self, _: &Remap) -> Self {
                self
            }
        }
    )*};
}
plain!((), bool, u32, u64, usize, i32, i64);

macro_rules! handle {
    ($($kind:ty => $via:ident),*) => {$(
        impl Relocate for $kind {
            fn relocate(self, remap: &Remap) -> Self {
                remap.$via(self)
            }
        }
    )*};
}
handle!(
    VertexId => vertex,
    EdgeId => edge,
    CoedgeId => coedge,
    LoopId => loop_,
    FaceId => face,
    ShellId => shell,
    LumpId => lump,
    BodyId => body
);

impl Relocate for SurfaceId {
    fn relocate(self, remap: &Remap) -> Self {
        remap.geom.surface(self)
    }
}

impl Relocate for CurveId {
    fn relocate(self, remap: &Remap) -> Self {
        remap.geom.curve(self)
    }
}

impl Relocate for PointId {
    fn relocate(self, remap: &Remap) -> Self {
        remap.geom.point(self)
    }
}

impl<T: Relocate> Relocate for Option<T> {
    fn relocate(self, remap: &Remap) -> Self {
        self.map(|value| value.relocate(remap))
    }
}

impl<T: Relocate> Relocate for Vec<T> {
    fn relocate(self, remap: &Remap) -> Self {
        self.into_iter()
            .map(|value| value.relocate(remap))
            .collect()
    }
}

impl<A: Relocate, B: Relocate> Relocate for (A, B) {
    fn relocate(self, remap: &Remap) -> Self {
        (self.0.relocate(remap), self.1.relocate(remap))
    }
}

// A refusal names no handle.
impl<T: Relocate, E> Relocate for Result<T, E> {
    fn relocate(self, remap: &Remap) -> Self {
        self.map(|value| value.relocate(remap))
    }
}

impl Model {
    fn extent(&self) -> [u32; 8] {
        [
            self.vertices.len(),
            self.edges.len(),
            self.coedges.len(),
            self.loops.len(),
            self.faces.len(),
            self.shells.len(),
            self.lumps.len(),
            self.bodies.len(),
        ]
        .map(|len| u32::try_from(len).expect("arena handle space exhausted"))
    }

    /// Starts a journal of everything added to this model and its geometry
    /// store from here on. Take a copy of the model first: the journal is for
    /// replaying onto the original.
    pub fn journal_begin(&mut self) {
        self.geom.journal_begin();
        self.journal = Some(TopoJournal {
            base: self.extent(),
            ops: Vec::new(),
        });
    }

    /// Ends the journal and returns it. Panics when none was begun.
    pub fn journal_end(&mut self) -> Journal {
        Journal {
            geom: self.geom.journal_end(),
            topo: self.journal.take().expect("a journal was begun"),
        }
    }

    /// Adds to this model, in order, what a journal's copy had added to it,
    /// and returns where the copy's new handles landed here. This model must
    /// be the one the copy was taken from, or have grown from it.
    ///
    /// Geometry goes first. Topology refers to geometry and never the other
    /// way, so the two orders are independent.
    pub fn replay(&mut self, journal: &Journal) -> Remap {
        let geom = self.geom.replay(&journal.geom);
        let extent = self.extent();
        let remap = Remap {
            geom,
            base: journal.topo.base,
            shift: std::array::from_fn(|arena| {
                extent[arena]
                    .checked_sub(journal.topo.base[arena])
                    .expect("a journal is replayed onto a model grown from its own")
            }),
        };
        for op in &journal.topo.ops {
            match op {
                TopoOp::Vertex(vertex) => {
                    self.add_vertex(Vertex {
                        point: geom.point(vertex.point),
                    });
                }
                TopoOp::Edge(edge) => {
                    self.add_edge(Edge {
                        curve: geom.curve(edge.curve),
                        start: edge.start.map(|v| remap.vertex(v)),
                        end: edge.end.map(|v| remap.vertex(v)),
                        ..edge.clone()
                    });
                }
                TopoOp::Coedge(coedge) => {
                    self.add_coedge(Coedge {
                        edge: remap.edge(coedge.edge),
                        ..coedge.clone()
                    });
                }
                TopoOp::Loop(l) => {
                    self.add_loop(match l {
                        Loop::Coedges(coedges) => {
                            Loop::Coedges(coedges.iter().map(|c| remap.coedge(*c)).collect())
                        }
                        Loop::Vertex(vertex) => Loop::Vertex(remap.vertex(*vertex)),
                    });
                }
                TopoOp::Face(face) => {
                    self.add_face(remap.face_value(face));
                }
                TopoOp::NamedFace(face, origins) => {
                    self.add_named_face(remap.face_value(face), origins.clone());
                }
                TopoOp::EdgeConstructions(edge, names) => {
                    self.record_edge_constructions(remap.edge(*edge), names.clone());
                }
                TopoOp::BlendBase(blend, base) => {
                    self.record_blend_base(remap.body(*blend), remap.body(*base));
                }
                TopoOp::Shell(shell) => {
                    self.add_shell(Shell {
                        faces: shell.faces.iter().map(|f| remap.face(*f)).collect(),
                    });
                }
                TopoOp::Lump(lump) => {
                    self.add_lump(Lump {
                        shells: lump.shells.iter().map(|s| remap.shell(*s)).collect(),
                    });
                }
                TopoOp::Body(body) => {
                    self.add_body(Body {
                        lumps: body.lumps.iter().map(|l| remap.lump(*l)).collect(),
                        ..body.clone()
                    });
                }
            }
        }
        remap
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Sense;
    use catallax_numbers::{Decider, Real};
    use catallax_shapes::{Curve, IntersectionCurve, Plane, Point3, PointDef, Surface, Vector3};

    fn axis(x: i64, y: i64, z: i64) -> Vector3 {
        Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
    }

    fn point(model: &mut Model, x: i64) -> PointId {
        model.geom_mut().add_point(PointDef::Explicit(Point3::new(
            Real::from_i64(x),
            Real::zero(),
            Real::zero(),
        )))
    }

    fn plane(model: &mut Model, flat: bool) -> SurfaceId {
        let mut dec = Decider::default();
        let (normal, along) = if flat {
            (axis(0, 0, 1), axis(1, 0, 0))
        } else {
            (axis(1, 0, 0), axis(0, 1, 0))
        };
        model.geom_mut().add_surface(Surface::Plane(
            Plane::new(Point3::origin(), normal, along, &mut dec).unwrap(),
        ))
    }

    // A model with one edge between two vertices on the meeting of two planes.
    fn base() -> (Model, SurfaceId, CurveId, VertexId, EdgeId) {
        let mut model = Model::new();
        let flat = plane(&mut model, true);
        let upright = plane(&mut model, false);
        let witness = point(&mut model, 0);
        let curve = model
            .geom_mut()
            .add_curve(Curve::Intersection(IntersectionCurve {
                surfaces: [flat, upright],
                witness,
            }));
        let start = model.add_vertex(Vertex { point: witness });
        let far = point(&mut model, 1);
        let end = model.add_vertex(Vertex { point: far });
        let edge = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(start),
            end: Some(end),
        });
        (model, flat, curve, start, edge)
    }

    type Made = ((VertexId, EdgeId), (LoopId, FaceId));

    // A piece of work that ties new topology to old: a vertex on a new
    // point, an edge from an old vertex to it, a coedge on the new edge and
    // one on an old edge, a loop of both, and two faces on that loop.
    fn work(
        model: &mut Model,
        surface: SurfaceId,
        curve: CurveId,
        vertex: VertexId,
        edge: EdgeId,
    ) -> Made {
        let place = point(model, 5);
        let corner = model.add_vertex(Vertex { point: place });
        let side = model.add_edge(Edge {
            curve,
            sense: Sense::Reversed,
            start: Some(vertex),
            end: Some(corner),
        });
        let out = model.add_coedge(Coedge {
            edge: side,
            sense: Sense::Same,
        });
        let back = model.add_coedge(Coedge {
            edge,
            sense: Sense::Reversed,
        });
        let ring = model.add_loop(Loop::Coedges(vec![out, back]));
        // The first face to hold a coedge owns it; the second does not.
        let face = model.add_face(Face {
            surface,
            sense: Sense::Same,
            loops: vec![ring],
        });
        model.add_named_face(
            Face {
                surface,
                sense: Sense::Reversed,
                loops: vec![ring],
            },
            Vec::new(),
        );
        ((corner, side), (ring, face))
    }

    #[test]
    fn a_journal_replays_onto_a_model_that_has_grown() {
        let (mut original, surface, curve, vertex, edge) = base();

        let mut copy = original.clone();
        copy.journal_begin();
        let made = work(&mut copy, surface, curve, vertex, edge);
        let journal = copy.journal_end();
        assert!(!journal.is_empty());

        // The original gains a point, a vertex, an edge and a coedge on the
        // old edge before the replay.
        let extra = point(&mut original, 9);
        let extra = original.add_vertex(Vertex { point: extra });
        original.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(vertex),
            end: Some(extra),
        });
        original.add_coedge(Coedge {
            edge,
            sense: Sense::Same,
        });

        let mut direct = original.clone();
        let expected = work(&mut direct, surface, curve, vertex, edge);

        let remap = original.replay(&journal);
        assert_ne!(made, expected, "the handles had to move");
        assert_eq!(made.relocate(&remap), expected);
        assert_eq!(remap.edge(edge), edge, "old handles stay");
        assert_eq!(
            original.coedges_of_edge(edge),
            direct.coedges_of_edge(edge),
            "the old edge has its coedges in the order they were added"
        );
        assert_eq!(format!("{original:?}"), format!("{direct:?}"));
    }

    #[test]
    fn a_journal_inside_a_journal_is_part_of_the_outer_one() {
        let (mut original, surface, curve, vertex, edge) = base();
        let mut direct = original.clone();
        work(&mut direct, surface, curve, vertex, edge);
        work(&mut direct, surface, curve, vertex, edge);

        let mut outer = original.clone();
        outer.journal_begin();
        work(&mut outer, surface, curve, vertex, edge);
        // Work on a copy of the copy, replayed onto the copy.
        let mut inner = outer.clone();
        inner.journal_begin();
        work(&mut inner, surface, curve, vertex, edge);
        let journal = inner.journal_end();
        outer.replay(&journal);
        let journal = outer.journal_end();
        original.replay(&journal);
        assert_eq!(format!("{original:?}"), format!("{direct:?}"));
    }
}
