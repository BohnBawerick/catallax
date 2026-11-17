//! Topological entities.

use catallax_numbers::{
    MarkLog, Refusal,
    transfer::{Codec, Reader, Registry, Writer},
};
use catallax_shapes::{CurveId, Handle, PointId, SurfaceId};

/// Handle of a vertex.
pub type VertexId = Handle<Vertex>;
/// Handle of an edge.
pub type EdgeId = Handle<Edge>;
/// Handle of a coedge.
pub type CoedgeId = Handle<Coedge>;
/// Handle of a loop.
pub type LoopId = Handle<Loop>;
/// Handle of a face.
pub type FaceId = Handle<Face>;
/// Handle of a shell.
pub type ShellId = Handle<Shell>;
/// Handle of a lump.
pub type LumpId = Handle<Lump>;
/// Handle of a body.
pub type BodyId = Handle<Body>;

/// Direction relative to an underlying direction.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, catallax_numbers::Transfer)]
pub enum Sense {
    /// The same direction.
    Same,
    /// The reverse direction.
    Reversed,
}

/// A vertex: a point stored by definition.
#[derive(Clone, PartialEq, Eq, Hash, Debug, catallax_numbers::Transfer)]
pub struct Vertex {
    /// The point.
    pub point: PointId,
}

/// An edge on a curve. Both ends `None` is a ring edge over a whole closed curve.
#[derive(Clone, PartialEq, Eq, Hash, Debug, catallax_numbers::Transfer)]
pub struct Edge {
    /// The curve.
    pub curve: CurveId,
    /// Edge direction relative to the curve's natural direction.
    pub sense: Sense,
    /// Start vertex along the edge direction.
    pub start: Option<VertexId>,
    /// End vertex along the edge direction.
    pub end: Option<VertexId>,
}

/// One directed use of an edge in a face boundary or a free wire.
#[derive(Clone, PartialEq, Eq, Hash, Debug, catallax_numbers::Transfer)]
pub struct Coedge {
    /// The edge.
    pub edge: EdgeId,
    /// Direction of travel relative to the edge direction.
    pub sense: Sense,
}

/// A boundary of a face.
#[derive(Clone, PartialEq, Eq, Hash, Debug, catallax_numbers::Transfer)]
pub enum Loop {
    /// A closed cycle of coedges.
    Coedges(Vec<CoedgeId>),
    /// A single point boundary, such as a cone apex.
    Vertex(VertexId),
}

/// A face on a surface. Zero loops means the whole closed surface.
#[derive(Clone, PartialEq, Eq, Hash, Debug, catallax_numbers::Transfer)]
pub struct Face {
    /// The surface.
    pub surface: SurfaceId,
    /// Face normal relative to the surface normal; the face normal points out of the material.
    pub sense: Sense,
    /// Boundaries, none of them outer.
    pub loops: Vec<LoopId>,
}

/// A connected set of faces.
#[derive(Clone, PartialEq, Eq, Hash, Debug, catallax_numbers::Transfer)]
pub struct Shell {
    /// The faces.
    pub faces: Vec<FaceId>,
}

/// One connected piece of material: an outer shell and its voids.
#[derive(Clone, PartialEq, Eq, Hash, Debug, catallax_numbers::Transfer)]
pub struct Lump {
    /// For solids, the outer shell first, then voids. Surface shells enclose no material.
    pub shells: Vec<ShellId>,
}

/// An ordered, connected path of coedges, open or closed, without a face.
#[derive(Clone, PartialEq, Eq, Hash, Debug, catallax_numbers::Transfer)]
pub struct Wire {
    /// Directed edge uses in traversal order.
    pub coedges: Vec<CoedgeId>,
}

/// What a body is.
#[non_exhaustive]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, catallax_numbers::Transfer)]
pub enum BodyKind {
    /// A solid with zero or more lumps.
    Solid,
    /// A collection of bounded faces, with no assertion of enclosed material.
    Surface,
    /// No material at all.
    Empty,
    /// One or more connected edge paths, without faces or material.
    Wire,
}

/// A body: material/surface lumps or free wires, plus construction marks.
#[derive(Clone, PartialEq, Debug)]
pub struct Body {
    /// What the body is.
    pub kind: BodyKind,
    /// Its lumps.
    pub lumps: Vec<LumpId>,
    /// Free wires. Nonempty only for a wire body, which has no lumps.
    pub wires: Vec<Wire>,
    /// Marks of every not-exact decision behind the body.
    pub marks: MarkLog,
}

impl Codec for Body {
    fn write(&self, writer: &mut Writer<'_>) -> Result<(), Refusal> {
        self.kind.write(writer)?;
        self.lumps.write(writer)?;
        self.wires.write(writer)?;
        self.marks.write(writer)
    }

    fn read(reader: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok(Self {
            kind: Codec::read(reader)?,
            lumps: Codec::read(reader)?,
            wires: if reader.version() == 1 {
                Vec::new()
            } else {
                Codec::read(reader)?
            },
            marks: Codec::read(reader)?,
        })
    }

    fn register(registry: &mut Registry) {
        registry.add::<BodyKind>();
        registry.add::<Vec<LumpId>>();
        registry.add::<Vec<Wire>>();
        registry.add::<MarkLog>();
    }
}
