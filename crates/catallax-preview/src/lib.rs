//! Unproved graphics-card preview routes for the field.
//!
//! Nothing in this crate is proved. Its values, meshes and pixels are display
//! estimates labelled unproved, and never feed an export, a proof or a proven
//! measurement. Exports and proven queries go through `catallax-field`'s
//! retained definitions and their proof routes.
pub mod batch;
pub mod first_paint;
pub mod gpu;
pub mod picking;
