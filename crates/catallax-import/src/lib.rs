//! STEP and IGES import with bounded, recorded gap repair.
//!
//! The contract is `docs/pieces/08-step-import-repair.md`. This module holds the
//! wire types that cross the application boundary: the proposals an import
//! returns, the decisions an application passes back, the persisted change
//! record and the typed refusals. Numbers are decimal strings in millimetres.

mod wire;

pub use wire::*;
mod build;
mod cone;
mod consume;
mod curved;
mod deviation;
mod exact;
mod gap;
mod iges;
mod incidence;
mod incidence_arithmetic;
mod incidence_features;
mod incidence_repair;
mod incidence_straight;
mod intersect;
/// Closed triangle mesh import with exact validity and winding queries.
pub mod mesh;
mod part21;
mod recognize;
mod repair;
mod spline_build;

pub use build::{build_file, build_solid};
pub use consume::{import_iges_bodies, import_step_bodies, sha256_hex};
pub use gap::{ImportOptions, propose, propose_source};
pub use iges::{SourceFormat, read_iges, read_iges_with_decider};
pub use incidence::measure_incidence;
pub use intersect::check_self_intersection;
pub use recognize::*;
pub use repair::{Imported, apply, apply_source, canonical_bytes, rebuild, rebuild_source};
