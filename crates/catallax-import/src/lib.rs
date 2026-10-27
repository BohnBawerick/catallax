//! STEP import with bounded, recorded gap repair.
//!
//! The contract is `docs/pieces/08-step-import-repair.md`. This module holds the
//! wire types that cross the application boundary: the proposals an import
//! returns, the decisions an application passes back, the persisted change
//! record and the typed refusals. Numbers are decimal strings in millimetres.

mod wire;

pub use wire::*;
