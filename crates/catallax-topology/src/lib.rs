//! Piece 3 of catallax: faces and edges.
//!
//! The boundary representation over exact curved geometry: bodies, lumps,
//! shells, faces, loops, coedges, edges and vertices; primitive solids; and
//! topology checks.
//!
//! The contract is `docs/pieces/03-topology.md`.

mod build;
mod check;
mod entity;
mod model;
mod rolling_junction;
mod spline_prism;
mod transform;

pub use build::{block, cone, cylinder, empty, sphere, torus};
pub use entity::{
    Body, BodyId, BodyKind, Coedge, CoedgeId, Edge, EdgeId, Face, FaceId, Loop, LoopId, Lump,
    LumpId, Sense, Shell, ShellId, Vertex, VertexId,
};
pub use model::{EntityCounts, Model, TopologyFault};
pub use rolling_junction::{envelope_retains_chart, envelope_rows};
pub use spline_prism::SplinePrism;
pub use transform::transformed_body;
