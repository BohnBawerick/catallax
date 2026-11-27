//! Unproved B-rep previews in buffers that the renderer can draw directly.
//!
//! This crate never supplies a checked mesh, a proof, an export, or a
//! geometric query. Face IDs identify draw ranges; they do not resolve picks.

mod gpu;
pub mod mesher;
