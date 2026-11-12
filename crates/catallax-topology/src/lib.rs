//! Piece 3 of catallax: faces and edges.
//!
//! The boundary representation over exact curved geometry: bodies, lumps,
//! shells, faces, loops, coedges, edges and vertices; primitive solids; and
//! topology checks.
//!
//! The contract is `docs/pieces/03-topology.md`.

mod bezier_extrusion;
mod build;
mod check;
mod entity;
mod linear_fillet;
mod model;
mod oblique_draft;
mod rational_atlas;
mod rational_patch;
mod rolling_junction;
mod rolling_mitre;
mod spline_prism;
mod transform;

pub use bezier_extrusion::BezierExtrusionBand;
pub use build::{block, cone, cylinder, empty, sphere, torus};
pub use entity::{
    Body, BodyId, BodyKind, Coedge, CoedgeId, Edge, EdgeId, Face, FaceId, Loop, LoopId, Lump,
    LumpId, Sense, Shell, ShellId, Vertex, VertexId,
};
pub use linear_fillet::LinearFilletBand;
pub use model::{EntityCounts, Journal, Model, Reads, Relocate, Remap, TopologyFault};
pub use oblique_draft::DraftedCylinder;
pub use rational_atlas::{AtlasFaceGroup, ChartExpansion, coalesce_rational_charts};
pub use rational_patch::RationalPatchFace;
pub use rolling_junction::{
    band_cap, envelope_retains_chart, envelope_rows, pipe_outward, pipe_rows, saddle_added,
    saddle_rows,
};
pub use rolling_mitre::RollingMitreBand;
pub use spline_prism::SplinePrism;
pub use transform::transformed_body;
