//! Deterministic IGES 5.3 manifold-solid B-rep export.
//!
//! See the crate README for supported definitions and the whole-domain similarity
//! certificate. A geometric refusal never calls the supplied writer.

use catallax_numbers::{Decider, Rational, Refusal};
use catallax_shapes::{CurveKind, SurfaceKind};
use catallax_topology::{BodyId, Model};
use std::{fmt, io::Write};

mod exact;
mod export;
mod format;

/// Output length units. Kernel coordinates and all report bounds are millimetres.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Units {
    /// Millimetres, IGES unit flag 2.
    #[default]
    Millimetres,
}
impl Units {
    pub(crate) fn millimetres(self) -> Rational {
        match self {
            Self::Millimetres => Rational::one(),
        }
    }
}

/// Strict export settings, independent of modelling tolerance and STEP options.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IgesOptions {
    /// Inclusive bound against the retained exact definitions, in millimetres.
    pub tolerance: Rational,
    /// Significant decimal digits used to choose the common scale, 1..=100.
    /// Coordinate fields can have more digits to preserve exact incidence.
    pub scale_digits: u32,
    /// Units declared in the IGES global section.
    pub units: Units,
}
impl Default for IgesOptions {
    fn default() -> Self {
        Self {
            tolerance: Rational::from_ratio(1, 10_000_000).unwrap(),
            scale_digits: 17,
            units: Units::Millimetres,
        }
    }
}

/// A whole-domain certificate for the emitted boundary and all its trims.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IgesReport {
    /// Proven upper bound on two-sided file-to-definition boundary displacement,
    /// in millimetres, including every emitted vertex, curve and face.
    pub max_deviation: Rational,
    /// Maximum unscaled pcurve endpoint rounding bound, in millimetres.
    /// Zero for planar and cylindrical bands; conical slant distances can be irrational.
    pub trim_deviation: Rational,
    /// Exact positive similarity scale of all three-dimensional definitions.
    pub scale: Rational,
    /// Bound on `|x| + |y| + |z|` over the entire finite model boundary, in mm.
    pub coordinate_bound: Rational,
    /// Written units; bounds above remain in millimetres.
    pub units: Units,
    /// Number of manifold solids (one per lump).
    pub solids: usize,
    /// Faces, preserving the input face count.
    pub faces: usize,
    /// Shared topological edges, including ring edges.
    pub edges: usize,
    /// Export-only seam edges joining periodic ring boundaries.
    pub added_edges: usize,
    /// Written vertices, including one added vertex per vertexless ring edge.
    pub vertices: usize,
    /// Ring vertices required by the IGES edge-list representation.
    pub added_vertices: usize,
}

/// Typed failures. Unsupported definitions and unavailable certificates are
/// distinct from an accuracy request that the proven bound cannot meet.
#[derive(Debug)]
pub enum IgesError {
    /// Invalid options.
    InvalidOptions,
    /// The input does not satisfy the public closed-solid topology contract.
    InvalidTopology(String),
    /// A marked model lacks the exact-definition correspondence this writer needs.
    UncertifiedModel,
    /// An unsupported surface kind.
    UnsupportedSurface(SurfaceKind),
    /// An unsupported curve kind.
    UnsupportedCurve(CurveKind),
    /// A numeric definition cannot yet be encoded with a proved exact similarity.
    UnrepresentableDefinition(&'static str),
    /// The proven bound exceeds the requested inclusive tolerance.
    ToleranceExceeded {
        /// Proven upper bound in millimetres.
        proven: Box<Rational>,
        /// Requested maximum in millimetres.
        requested: Box<Rational>,
    },
    /// A shared exact-number decision refused.
    Decision(Refusal),
    /// Output failed after preparation; the sink may have accepted a prefix.
    Io(std::io::Error),
}
impl fmt::Display for IgesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}
impl std::error::Error for IgesError {}
impl From<Refusal> for IgesError {
    fn from(value: Refusal) -> Self {
        Self::Decision(value)
    }
}

/// Prepare and certify an entire IGES file, then write it to `out`.
///
/// The caller supplies valid closed, finite exact bodies. Empty bodies are
/// accepted. Unsupported definitions, work exhaustion and tighter requests leave
/// `out` untouched. I/O errors follow `Write::write_all` semantics.
/// No importer, healing, STEP export or kernel construction runs here.
pub fn write_iges<W: Write>(
    model: &Model,
    bodies: &[BodyId],
    options: &IgesOptions,
    dec: &mut Decider,
    out: &mut W,
) -> Result<IgesReport, IgesError> {
    if options.tolerance < Rational::zero() || !(1..=100).contains(&options.scale_digits) {
        return Err(IgesError::InvalidOptions);
    }
    if dec.policy().max_work == 0 && !bodies.is_empty() {
        return Err(IgesError::Decision(Refusal::new(
            "iges.work-limit",
            "no allowance for export proof",
        )));
    }
    let (file, report) = export::prepare(model, bodies, options, dec)?;
    out.write_all(file.as_bytes()).map_err(IgesError::Io)?;
    Ok(report)
}
