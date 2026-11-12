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

use catallax_numbers::Refusal;
use catallax_recipe::{EdgeName, Origin};
use catallax_shapes::{CurveId, GeomJournal, GeomReads, GeomRemap, Handle, PointId, SurfaceId};

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

// One thing a model was asked that later additions can change the answer to.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub(crate) enum TopoRead {
    FaceOrigins(FaceId),
    EdgeConstructions(EdgeId),
    BlendBase(BodyId),
    EdgeCoedges(EdgeId),
    CoedgeOwner(CoedgeId),
    ChartExpansion(BodyId),
    FaceCover(FaceId),
    // Which bodies there are.
    Bodies,
}

/// What work on a copy of a model asked it, between [`Model::watch`] and
/// [`Model::watched`], that later additions can change the answer to.
///
/// What a handle names never changes, so reading an entity by its handle is
/// not listed. Listed are the questions whose answers grow: the coedges of an
/// edge, the owner of a coedge, the names kept beside faces, edges and
/// bodies, which bodies there are, and everything the geometry store lists.
/// [`Model::unchanged_where_read`] checks them against the model the copy
/// came from, and a journal is replayed there only when they hold: the work
/// would then have done the same on that model.
#[derive(Clone, Debug, Default)]
pub struct Reads {
    geom: GeomReads,
    topo: Vec<TopoRead>,
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

impl<T: Relocate, E: Relocate> Relocate for Result<T, E> {
    fn relocate(self, remap: &Remap) -> Self {
        match self {
            Ok(value) => Ok(value.relocate(remap)),
            Err(error) => Err(error.relocate(remap)),
        }
    }
}

// A refusal is a code, a text for people and the decision it stopped on. It
// names no handle.
impl Relocate for Refusal {
    fn relocate(self, _: &Remap) -> Self {
        self
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

    /// Starts a list of what this model and its geometry store are asked from
    /// here on that later additions can change the answer to. Take a copy of
    /// the model first: the list is for checking the original.
    pub fn watch(&mut self) {
        self.geom.watch();
        self.watch.begin();
    }

    /// Ends the list and returns it. Panics when none was begun.
    pub fn watched(&mut self) -> Reads {
        Reads {
            geom: self.geom.watched(),
            topo: self.watch.end(),
        }
    }

    /// Whether this model answers everything in `reads` as `base` does. This
    /// model must be `base` or have grown from it, and `reads` must come from
    /// work on a copy of `base`. A read about an entity the work made itself
    /// always passes.
    ///
    /// When this holds, the work would have found the same answers on this
    /// model, so it would have added the same to it: its journal can be
    /// replayed here.
    pub fn unchanged_where_read(&self, base: &Model, reads: &Reads) -> bool {
        let before = base.extent();
        let old = |arena: usize, index: u32| index < before[arena];
        reads.topo.iter().all(|read| match *read {
            TopoRead::FaceOrigins(face) => {
                !old(4, face.index())
                    || base.face_origins.get(&face) == self.face_origins.get(&face)
            }
            TopoRead::EdgeConstructions(edge) => {
                !old(1, edge.index())
                    || base.edge_constructions.get(&edge) == self.edge_constructions.get(&edge)
            }
            TopoRead::BlendBase(body) => {
                !old(7, body.index()) || base.blend_bases.get(&body) == self.blend_bases.get(&body)
            }
            TopoRead::EdgeCoedges(edge) => {
                !old(1, edge.index())
                    || base.edge_coedges[edge.index() as usize]
                        == self.edge_coedges[edge.index() as usize]
            }
            TopoRead::CoedgeOwner(coedge) => {
                !old(2, coedge.index())
                    || base.coedge_owners[coedge.index() as usize]
                        == self.coedge_owners[coedge.index() as usize]
            }
            TopoRead::ChartExpansion(body) => {
                !old(7, body.index())
                    || match (
                        base.chart_expansions.get(&body),
                        self.chart_expansions.get(&body),
                    ) {
                        (None, None) => true,
                        (Some(a), Some(b)) => std::sync::Arc::ptr_eq(a, b),
                        _ => false,
                    }
            }
            TopoRead::FaceCover(face) => {
                !old(4, face.index())
                    || match (
                        base.face_expansions.get(&face),
                        self.face_expansions.get(&face),
                    ) {
                        (None, None) => true,
                        (Some(a), Some(b)) => std::sync::Arc::ptr_eq(a, b),
                        _ => false,
                    }
            }
            TopoRead::Bodies => self.bodies.len() == base.bodies.len(),
        }) && self.geom.unchanged_where_read(&base.geom, &reads.geom)
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

    #[test]
    fn a_handle_in_an_error_moves_with_the_replay() {
        let mut original = Model::new();

        // The work makes a point on its copy and fails with that point.
        let mut copy = original.clone();
        copy.journal_begin();
        let made = point(&mut copy, 5);
        let failed: Result<(), PointId> = Err(made);
        let also: Result<VertexId, (PointId, Option<PointId>)> = Err((made, Some(made)));
        let journal = copy.journal_end();

        // The original gains another point before the replay, so the work's
        // point lands one place later than it did on the copy.
        let other = point(&mut original, 9);
        assert_eq!(
            other, made,
            "the two points share a handle until the replay"
        );
        let remap = original.replay(&journal);
        let moved = remap.geom().point(made);
        assert_ne!(moved, made);

        let Err(named) = failed.relocate(&remap) else {
            unreachable!("an error stays an error")
        };
        assert_eq!(named, moved);
        assert_eq!(
            original.geom().approx_point(named),
            [5.0, 0.0, 0.0],
            "the error names the work's point, not the one that took its place"
        );
        assert_eq!(original.geom().approx_point(other), [9.0, 0.0, 0.0]);
        assert_eq!(also.relocate(&remap), Err((moved, Some(moved))));
    }

    // What `work` asked a watched copy of `base`.
    fn asked(base: &Model, work: impl FnOnce(&mut Model)) -> Reads {
        let mut copy = base.clone();
        copy.watch();
        work(&mut copy);
        copy.watched()
    }

    #[test]
    fn work_whose_reads_still_hold_replays_as_if_it_had_run_here() {
        let (mut original, surface, curve, vertex, edge) = base();

        // The work adds to the model and asks about what it added: the
        // coedges of its own edge, the owner of its own coedge.
        let mut copy = original.clone();
        copy.journal_begin();
        copy.watch();
        let made = work(&mut copy, surface, curve, vertex, edge);
        let ((_, side), (ring, face)) = made;
        let out = copy.coedges_of_edge(side)[0];
        assert_eq!(copy.owner_of_coedge(out), Some((ring, face)));
        let reads = copy.watched();
        let journal = copy.journal_end();

        // Meanwhile the original gains an edge of its own with a coedge.
        let extra = point(&mut original, 9);
        let extra = original.add_vertex(Vertex { point: extra });
        let theirs = original.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(vertex),
            end: Some(extra),
        });
        original.add_coedge(Coedge {
            edge: theirs,
            sense: Sense::Same,
        });
        let base = base().0;
        assert!(original.unchanged_where_read(&base, &reads));

        let mut direct = original.clone();
        let expected = work(&mut direct, surface, curve, vertex, edge);
        let remap = original.replay(&journal);
        assert_eq!(made.relocate(&remap), expected);
        assert_eq!(format!("{original:?}"), format!("{direct:?}"));
    }

    #[test]
    fn a_table_entry_the_work_read_is_checked() {
        let (mut base, surface, curve, vertex, edge) = base();
        let lone = base.add_coedge(Coedge {
            edge,
            sense: Sense::Same,
        });
        let ring = base.add_loop(Loop::Coedges(vec![lone]));

        // The coedges of an old edge.
        let reads = asked(&base, |model| {
            assert_eq!(model.coedges_of_edge(edge), vec![lone]);
        });
        assert!(base.unchanged_where_read(&base, &reads));
        let mut grown = base.clone();
        grown.add_coedge(Coedge {
            edge,
            sense: Sense::Reversed,
        });
        assert!(!grown.unchanged_where_read(&base, &reads));

        // The owner of an old coedge.
        let reads = asked(&base, |model| {
            assert_eq!(model.owner_of_coedge(lone), None);
        });
        let mut grown = base.clone();
        let face = grown.add_face(Face {
            surface,
            sense: Sense::Same,
            loops: vec![ring],
        });
        assert_eq!(grown.owner_of_coedge(lone), Some((ring, face)));
        assert!(!grown.unchanged_where_read(&base, &reads));

        // Which bodies there are.
        let reads = asked(&grown, |model| {
            assert!(model.bodies_with_face(face).is_empty());
        });
        let mut later = grown.clone();
        later.add_body(Body {
            kind: crate::BodyKind::Empty,
            lumps: Vec::new(),
            marks: catallax_numbers::MarkLog::new(),
        });
        assert!(!later.unchanged_where_read(&grown, &reads));

        // A fact of the geometry store, asked through a copy of the copy.
        let start = base.vertex(vertex).point;
        let on_curve = |model: &Model| {
            model.geom().facts().holds(
                catallax_numbers::FactKind::Incident,
                catallax_shapes::GeomKey::Point(start),
                catallax_shapes::GeomKey::Curve(curve),
            )
        };
        let reads = asked(&base, |model| {
            let scratch = model.clone();
            assert!(!on_curve(&scratch));
        });
        let mut grown = base.clone();
        grown.geom_mut().record_fact(catallax_numbers::Fact {
            kind: catallax_numbers::FactKind::Incident,
            a: catallax_shapes::GeomKey::Point(start),
            b: catallax_shapes::GeomKey::Curve(curve),
            at: None,
            origin: catallax_numbers::FactOrigin::Construction("test"),
        });
        assert!(on_curve(&grown));
        assert!(!grown.unchanged_where_read(&base, &reads));

        // What the work never asked about may change.
        let reads = asked(&base, |model| {
            let _ = model.edge(edge);
            let _ = model.geom().curve(curve);
        });
        assert!(grown.unchanged_where_read(&base, &reads));
    }
}
