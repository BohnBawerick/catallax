//! STEP import with bounded, recorded gap repair.
//!
//! The contract is `docs/pieces/08-step-import-repair.md`. This module holds the
//! wire types that cross the application boundary: the proposals an import
//! returns, the decisions an application passes back, the persisted change
//! record and the typed refusals. Numbers are decimal strings in millimetres.

mod wire;

pub use wire::*;
mod build;
mod consume;
mod curved;
mod exact;
mod gap;
mod intersect;
mod part21;
mod recognize;
mod repair;

pub use build::{build_file, build_solid};
pub use consume::{import_step_bodies, sha256_hex};
pub use gap::{ImportOptions, propose};
pub use intersect::check_self_intersection;
pub use recognize::*;
pub use repair::{Imported, apply, canonical_bytes, rebuild};
