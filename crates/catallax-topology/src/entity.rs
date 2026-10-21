//! Topological entities.

use catallax_numbers::MarkLog;
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
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Sense {
    /// The same direction.
    Same,
    /// The reverse direction.
    Reversed,
}

/// A vertex: a point stored by definition.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Vertex {
    /// The point.
    pub point: PointId,
}

/// An edge on a curve. Both ends `None` is a ring edge over a whole closed curve.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
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

/// One use of an edge by a loop.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Coedge {
    /// The edge.
    pub edge: EdgeId,
    /// Direction of travel relative to the edge direction.
    pub sense: Sense,
}

/// A boundary of a face.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum Loop {
    /// A closed cycle of coedges.
    Coedges(Vec<CoedgeId>),
    /// A single point boundary, such as a cone apex.
    Vertex(VertexId),
}

/// A face on a surface. Zero loops means the whole closed surface.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Face {
    /// The surface.
    pub surface: SurfaceId,
    /// Face normal relative to the surface normal; the face normal points out of the material.
    pub sense: Sense,
    /// Boundaries, none of them outer.
    pub loops: Vec<LoopId>,
}

/// A connected set of faces.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Shell {
    /// The faces.
    pub faces: Vec<FaceId>,
}

/// One connected piece of material: an outer shell and its voids.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Lump {
    /// Outer shell first, then voids.
    pub shells: Vec<ShellId>,
}

/// What a body is.
#[non_exhaustive]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum BodyKind {
    /// A solid with zero or more lumps.
    Solid,
    /// No material at all.
    Empty,
}

/// A body: lumps plus the marks of the operations that produced it.
#[derive(Clone, PartialEq, Debug)]
pub struct Body {
    /// What the body is.
    pub kind: BodyKind,
    /// Its lumps.
    pub lumps: Vec<LumpId>,
    /// Marks of every not-exact decision behind the body.
    pub marks: MarkLog,
}
