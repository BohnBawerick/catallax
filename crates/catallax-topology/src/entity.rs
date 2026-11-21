//! Topological entities.

use catallax_numbers::{
    MarkLog, Rational, Refusal, Sign,
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

/// The tolerance of a vertex, an edge or a face, in millimetres. Always
/// greater than zero.
///
/// One model holds exact and tolerant entities (decision Q12, section 4.5 of
/// the Checkpoint 2 decisions). The structure is Parasolid's, from Jackson
/// 1995, "Boundary representation modelling with local tolerances". A
/// tolerance belongs to topology, and points, curves and surfaces are stored
/// as they are.
///
/// An entity whose tolerance is `None` is exact. It stores no tolerance and
/// has no size limit (rule T0). `Some` says how far the entity's stored
/// geometry may lie from the geometry it meets.
///
/// - A vertex's is the largest distance from its point to the end of each of
///   its edges' curves.
/// - An edge's is the largest distance between its curve and its copy on
///   each of its faces.
/// - A face's is the largest distance between its surface and the surface it
///   stands for.
///
/// The builder that makes the entity measures the value (rules T1 to T4).
/// [`Model::check`] reads it and measures nothing.
///
/// [`Model::check`]: crate::Model::check
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Tolerance(Rational);

impl Tolerance {
    /// A tolerance of `mm` millimetres. Zero and negative values refuse,
    /// because exact is `None` and never a tolerance of zero.
    pub fn new(mm: Rational) -> Result<Self, Refusal> {
        if mm.sign() == Sign::Positive {
            Ok(Self(mm))
        } else {
            Err(Refusal::new(
                "topology.tolerance-not-positive",
                format!("a tolerance must be greater than zero; {mm} mm was given"),
            ))
        }
    }

    /// The value, in millimetres.
    pub fn mm(&self) -> &Rational {
        &self.0
    }
}

impl Codec for Tolerance {
    fn write(&self, writer: &mut Writer<'_>) -> Result<(), Refusal> {
        self.0.write(writer)
    }

    fn read(reader: &mut Reader<'_>) -> Result<Self, Refusal> {
        Self::new(Codec::read(reader)?)
            .map_err(|_| catallax_numbers::transfer::invalid("a stored tolerance is not positive"))
    }

    fn register(registry: &mut Registry) {
        registry.add::<Rational>();
    }
}

/// Reads the tolerance of a vertex, an edge or a face. A record written
/// before format version 3 has none. That is right for an exact body. The
/// reader of a numerical body gives its model's entities theirs again
/// (`NumericalBody::restore_model_tolerances`).
fn read_tolerance(reader: &mut Reader<'_>) -> Result<Option<Tolerance>, Refusal> {
    if reader.version() < 3 {
        Ok(None)
    } else {
        Codec::read(reader)
    }
}

/// A vertex: a point stored by definition.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Vertex {
    /// The point.
    pub point: PointId,
    /// How far the point may lie from the ends of the vertex's edges, or
    /// `None` for an exact vertex. Never less than the tolerance of an edge
    /// that ends here.
    pub tolerance: Option<Tolerance>,
}

impl Codec for Vertex {
    fn write(&self, writer: &mut Writer<'_>) -> Result<(), Refusal> {
        self.point.write(writer)?;
        self.tolerance.write(writer)
    }

    fn read(reader: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok(Self {
            point: Codec::read(reader)?,
            tolerance: read_tolerance(reader)?,
        })
    }

    fn register(registry: &mut Registry) {
        registry.add::<PointId>();
        registry.add::<Option<Tolerance>>();
    }
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
    /// How far the curve may lie from its copy on each face of the edge, or
    /// `None` for an exact edge. Never less than the tolerance of a face the
    /// edge bounds.
    pub tolerance: Option<Tolerance>,
}

impl Codec for Edge {
    fn write(&self, writer: &mut Writer<'_>) -> Result<(), Refusal> {
        self.curve.write(writer)?;
        self.sense.write(writer)?;
        self.start.write(writer)?;
        self.end.write(writer)?;
        self.tolerance.write(writer)
    }

    fn read(reader: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok(Self {
            curve: Codec::read(reader)?,
            sense: Codec::read(reader)?,
            start: Codec::read(reader)?,
            end: Codec::read(reader)?,
            tolerance: read_tolerance(reader)?,
        })
    }

    fn register(registry: &mut Registry) {
        registry.add::<CurveId>();
        registry.add::<Sense>();
        registry.add::<Option<VertexId>>();
        registry.add::<Option<Tolerance>>();
    }
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
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Face {
    /// The surface.
    pub surface: SurfaceId,
    /// Face normal relative to the surface normal; the face normal points out of the material.
    pub sense: Sense,
    /// Boundaries, none of them outer.
    pub loops: Vec<LoopId>,
    /// How far the surface may lie from the surface it stands for, or `None`
    /// for an exact face. A face on an exact surface never takes one: a
    /// tolerant edge or vertex on it does (section 4.5 of the Checkpoint 2
    /// decisions).
    pub tolerance: Option<Tolerance>,
}

impl Codec for Face {
    fn write(&self, writer: &mut Writer<'_>) -> Result<(), Refusal> {
        self.surface.write(writer)?;
        self.sense.write(writer)?;
        self.loops.write(writer)?;
        self.tolerance.write(writer)
    }

    fn read(reader: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok(Self {
            surface: Codec::read(reader)?,
            sense: Codec::read(reader)?,
            loops: Codec::read(reader)?,
            tolerance: read_tolerance(reader)?,
        })
    }

    fn register(registry: &mut Registry) {
        registry.add::<SurfaceId>();
        registry.add::<Sense>();
        registry.add::<Vec<LoopId>>();
        registry.add::<Option<Tolerance>>();
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Decider, transfer};

    #[test]
    fn a_tolerance_is_greater_than_zero() {
        for refused in ["0", "-1/1000"] {
            assert_eq!(
                Tolerance::new(Rational::parse(refused).unwrap())
                    .unwrap_err()
                    .code,
                "topology.tolerance-not-positive"
            );
        }
        let small = Tolerance::new(Rational::parse("1/1000000000").unwrap()).unwrap();
        assert_eq!(small.mm(), &Rational::parse("1/1000000000").unwrap());
        // Exact is `None`, and it is less than every tolerance.
        assert!(None < Some(small));
    }

    #[test]
    fn a_stored_tolerance_of_zero_is_refused_on_reading() {
        let mut dec = Decider::default();
        let written = transfer::encode(&Tolerance(Rational::zero()), &mut dec).unwrap();
        assert_eq!(
            transfer::decode::<Tolerance>(&written, &mut dec)
                .unwrap_err()
                .code,
            "handover.invalid-data"
        );
        let kept = Tolerance::new(Rational::parse("1/4500").unwrap()).unwrap();
        let written = transfer::encode(&kept, &mut dec).unwrap();
        assert_eq!(
            transfer::decode::<Tolerance>(&written, &mut dec).unwrap(),
            kept
        );
    }
}
