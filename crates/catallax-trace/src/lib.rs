//! Spline-intersection tracing adapted from NVIDIA USD BRep SMLib.
//!
//! Everything this crate returns is an unproved floating-point candidate.
//! Nothing here decides a Catallax answer: a traced point or curve becomes
//! part of a result only after Catallax's own exact proofs accept it, and a
//! candidate those proofs reject is a typed refusal, never a near-enough
//! answer. SMLib's tolerances steer the search here; they are never a
//! Catallax acceptance rule.
//!
//! The port follows the unit order of the M25 scope study. Each module names
//! its SMLib source and, where one exists, the rust-brep unit it adapts. Each
//! unit is gated against the pinned C++ SMLib through saved oracle fixtures
//! under `tests/fixtures`; Cargo never builds or links C++.
//!
//! Arithmetic is IEEE binary64 without fused multiply-add, and the inverse
//! trigonometric functions come from `libm`, so a candidate is the same on
//! every machine and at every thread count.

// Ported comparisons reproduce SMLib's exact floating-point tests, including
// equality against sentinels and against 0.0 and 1.0.
#![allow(clippy::float_cmp)]

pub mod extent;
pub mod geometry;
pub mod global_solver;
pub mod local_solve;
pub mod matrix;
pub mod scalar;
mod search_math;
pub mod solution;
pub mod vector;

/// Why a ported routine produced no value.
///
/// The two kinds mirror SMLib's `SM_ERR` and `SM_ERR_INVALID_INPUT`, which
/// the differential gate compares exactly.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Fault {
    /// A degenerate or unsupported configuration (`SM_ERR`).
    Degenerate,
    /// Input outside the routine's domain (`SM_ERR_INVALID_INPUT`).
    InvalidInput,
}
