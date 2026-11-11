// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: constants and clamped scalar helpers only; error
// reporting through SMLib's global error state is dropped, and the inverse
// sine and cosine come from libm so candidates match on every machine.
// Adapted from SMLib SmMath.h, SmMath.cpp (sm_sqrt, sm_acos, sm_asin) and
// SmTol.cpp (GetScaledZero), with rust-brep crates/kernel/src/vector2d.rs and
// sm_geom_utility.rs, commit 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab.
//! SMLib's numeric constants and clamped scalar functions.
//!
//! These values steer the candidate search only. No Catallax proof reads them.

/// `SM_EFF_ZERO`, SMLib's unscaled numerical zero.
pub const EFF_ZERO: f64 = 1.0e-12;
/// `SM_EFF_ZERO_SQ`.
pub const EFF_ZERO_SQ: f64 = 1.0e-24;
/// `SM_EFF_ZERO_SQRT`.
pub const EFF_ZERO_SQRT: f64 = 1.0e-6;
/// `SM_ZONE_TOL_3D`, SMLib's default neighbourhood size in model units.
pub const ZONE_TOL_3D: f64 = 1.0e-5;
/// `SM_ANG_TOL_DEG`, also SMLib's `SM_CONTINUITY_ANGLE`, in degrees.
pub const ANGLE_TOL_DEG: f64 = 0.5;
/// `SM_BIG_DOUBLE`, SMLib's "no finite answer" value.
pub const BIG_DOUBLE: f64 = 1.0e20;
/// `SM_UNDEF_DOUBLE`, the coordinate of an uninitialised SMLib vector.
pub const UNDEFINED: f64 = 1.23456e-61;
/// `SM_PI`, which rounds to the binary64 value of pi.
pub const PI: f64 = std::f64::consts::PI;

/// `smos_Max`: `a > b ? a : b`, which differs from `f64::max` on NaN.
pub fn max(a: f64, b: f64) -> f64 {
    if a > b { a } else { b }
}

/// `smos_Min`: `a < b ? a : b`.
pub fn min(a: f64, b: f64) -> f64 {
    if a < b { a } else { b }
}

/// `sm_sqrt`: the square root, with zero for zero and negative arguments.
pub fn sqrt(x: f64) -> f64 {
    if x < -EFF_ZERO_SQ || x <= 0.0 {
        0.0
    } else {
        x.sqrt()
    }
}

/// `sm_acos`: the arc cosine, clamped to `[0, pi]` outside `[-1, 1]`.
pub fn acos(x: f64) -> f64 {
    if x >= 1.0 {
        0.0
    } else if x <= -1.0 {
        PI
    } else {
        libm::acos(x)
    }
}

/// `sm_asin`: the arc sine, clamped to `[-pi/2, pi/2]` outside `[-1, 1]`.
pub fn asin(x: f64) -> f64 {
    if x >= 1.0 {
        PI / 2.0
    } else if x <= -1.0 {
        -PI / 2.0
    } else {
        libm::asin(x)
    }
}

/// `smos_Cosine`: the cosine, from `libm` for repeatable last bits.
pub fn cos(x: f64) -> f64 {
    libm::cos(x)
}

/// `SmTol::GetScaledZero(double)`: `SM_EFF_ZERO * (1 + |value|)`.
pub fn scaled_zero(value: f64) -> f64 {
    EFF_ZERO * (1.0 + value.abs())
}

/// `SM_ARE_SAME(a, b)`: `|a - b|` within the larger scaled zero of the two.
pub fn are_same(a: f64, b: f64) -> bool {
    (a - b).abs() <= max(scaled_zero(a), scaled_zero(b))
}

/// `SM_IS_ZERO(a)`: `|a| <= SM_EFF_ZERO`.
pub fn is_zero(a: f64) -> bool {
    a.abs() <= EFF_ZERO
}

/// `SM_RAD2DEG`.
pub fn degrees(radians: f64) -> f64 {
    radians * (180.0 / PI)
}
