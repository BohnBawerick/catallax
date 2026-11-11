// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: only the vector operations the intersection chain
// enters; statuses become `Fault`; cross products are named methods rather
// than overloaded `*`; inverse trigonometry through libm.
// Adapted from rust-brep crates/kernel/src/vector2d.rs and vector3d.rs, commit
// 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab, translations of SMLib
// SmVector2d.h, SmVector3d.h and SmVector3d.cpp.
//! Two- and three-dimensional binary64 vectors with SMLib's semantics.
//!
//! Division by a scalar keeps SMLib's guards: the 2D quotient and the 3D
//! in-place quotient return the vector unchanged when the divisor is
//! negligible against every coordinate, while the 3D quotient operator
//! divides unconditionally.

use crate::Fault;
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ, PI, UNDEFINED, acos, asin, sqrt};
use std::ops::{Add, Div, DivAssign, Mul, MulAssign, Neg, Sub};

/// A 2D vector or parameter point (`SmVector2d`, `SmPoint2d`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vector2 {
    /// First coordinate, or the u parameter.
    pub x: f64,
    /// Second coordinate, or the v parameter.
    pub y: f64,
}

impl Default for Vector2 {
    /// SMLib's uninitialised vector.
    fn default() -> Self {
        Self::new(UNDEFINED, UNDEFINED)
    }
}

impl Vector2 {
    /// The vector `(x, y)`.
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// `LengthSquared`.
    pub fn length_squared(self) -> f64 {
        (self.x * self.x) + (self.y * self.y)
    }

    /// `Length`.
    pub fn length(self) -> f64 {
        sqrt(self.length_squared())
    }

    /// `Dot`.
    pub fn dot(self, other: Self) -> f64 {
        (self.x * other.x) + (self.y * other.y)
    }

    /// `DistanceBetween`.
    pub fn distance(self, other: Self) -> f64 {
        (self - other).length()
    }

    /// `DistanceBetweenSquared`.
    pub fn distance_squared(self, other: Self) -> f64 {
        (self - other).length_squared()
    }

    /// `GetMaxDimension`: the largest coordinate magnitude.
    pub fn max_dimension(self) -> f64 {
        if self.x.abs() > self.y.abs() {
            self.x.abs()
        } else {
            self.y.abs()
        }
    }

    /// `operator!=`: the coordinate gap reaches the scaled zero of `self`.
    pub fn differs(self, other: Self) -> bool {
        let scale = 1.0 + self.max_dimension();
        (self.x - other.x).abs() + (self.y - other.y).abs() >= EFF_ZERO * scale
    }
}

impl Neg for Vector2 {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y)
    }
}

impl Add for Vector2 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y)
    }
}

impl Sub for Vector2 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self::new(self.x - other.x, self.y - other.y)
    }
}

impl Mul<f64> for Vector2 {
    type Output = Self;
    fn mul(self, scale: f64) -> Self {
        Self::new(self.x * scale, self.y * scale)
    }
}

impl Mul<Vector2> for f64 {
    type Output = Vector2;
    fn mul(self, vector: Vector2) -> Vector2 {
        vector * self
    }
}

impl Div<f64> for Vector2 {
    type Output = Self;
    /// Returns `self` unchanged when the divisor is negligible against a coordinate.
    fn div(self, divisor: f64) -> Self {
        if divisor.abs() < self.x.abs() * EFF_ZERO || divisor.abs() < self.y.abs() * EFF_ZERO {
            self
        } else {
            Self::new(self.x / divisor, self.y / divisor)
        }
    }
}

/// A 3D vector or point (`SmVector3d`, `SmPoint3d`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vector3 {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
    /// Z coordinate.
    pub z: f64,
}

impl Default for Vector3 {
    /// SMLib's uninitialised vector.
    fn default() -> Self {
        Self::new(UNDEFINED, UNDEFINED, UNDEFINED)
    }
}

impl From<Vector2> for Vector3 {
    /// `SmVector3d(const SmVector2d&)`: z is zero.
    fn from(v: Vector2) -> Self {
        Self::new(v.x, v.y, 0.0)
    }
}

impl Vector3 {
    /// The vector `(x, y, z)`.
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// `IsInitialized`: no coordinate is SMLib's undefined value.
    pub fn is_initialized(self) -> bool {
        self.x != UNDEFINED && self.y != UNDEFINED && self.z != UNDEFINED
    }

    /// `IsUndef`: a coordinate is undefined or outside `[-1e18, 1e18]`.
    pub fn is_undefined(self) -> bool {
        let valid = |v: f64| (-1e18..=1e18).contains(&v);
        [self.x, self.y, self.z]
            .into_iter()
            .any(|v| v == UNDEFINED || !valid(v))
    }

    /// `LengthSquared`.
    pub fn length_squared(self) -> f64 {
        (self.x * self.x) + (self.y * self.y) + (self.z * self.z)
    }

    /// `Length`.
    pub fn length(self) -> f64 {
        sqrt(self.length_squared())
    }

    /// `Dot`.
    pub fn dot(self, b: Self) -> f64 {
        (self.x * b.x) + (self.y * b.y) + (self.z * b.z)
    }

    /// The cross product, SMLib's `operator*` between two vectors.
    pub fn cross(self, b: Self) -> Self {
        Self::new(
            (self.y * b.z) - (self.z * b.y),
            (self.z * b.x) - (self.x * b.z),
            (self.x * b.y) - (self.y * b.x),
        )
    }

    /// `Unitize`. A vector shorter than `SM_EFF_ZERO` becomes `(0, 0, 1)`
    /// and reports invalid input; a unit vector is left untouched.
    pub fn unitize(&mut self) -> Result<(), Fault> {
        let squared = self.length_squared();
        if squared < EFF_ZERO_SQ {
            *self = Self::new(0.0, 0.0, 1.0);
            return Err(Fault::InvalidInput);
        }
        if (squared - 1.0).abs() > EFF_ZERO_SQ {
            let length = sqrt(squared);
            self.x /= length;
            self.y /= length;
            self.z /= length;
        }
        Ok(())
    }

    /// `GetMaxDimension`: the largest coordinate magnitude.
    pub fn max_dimension(self) -> f64 {
        let mut result = if self.x.abs() > self.y.abs() {
            self.x.abs()
        } else {
            self.y.abs()
        };
        if self.z.abs() > result {
            result = self.z.abs();
        }
        result
    }

    /// `DistanceBetween`.
    pub fn distance(self, other: Self) -> f64 {
        sqrt(self.distance_squared(other))
    }

    /// `DistanceBetweenSquared`.
    pub fn distance_squared(self, other: Self) -> f64 {
        (self - other).length_squared()
    }

    /// `CloserThan`: every coordinate gap and the distance are within `maximum`.
    pub fn closer_than(self, maximum: f64, other: Self) -> bool {
        let (dx, dy, dz) = (self.x - other.x, self.y - other.y, self.z - other.z);
        !(dx.abs() > maximum
            || dy.abs() > maximum
            || dz.abs() > maximum
            || dx * dx + dy * dy + dz * dz > maximum * maximum)
    }

    /// `AngleBetween`: the angle in `[0, pi]`, from the sine near 0 and pi
    /// to avoid the arc cosine's loss of precision there.
    pub fn angle_between(self, other: Self) -> Result<f64, Fault> {
        let a = self.length_squared();
        let b = other.length_squared();
        if a.abs() <= EFF_ZERO_SQ || b.abs() <= EFF_ZERO_SQ {
            return Err(Fault::InvalidInput);
        }
        let dot = self.dot(other);
        let sign = if dot < 0.0 { -1.0 } else { 1.0 };
        let squared = (dot * dot) / a / b;
        let cosine = sign * sqrt(squared);
        if squared < 0.81 {
            return Ok(acos(cosine));
        }
        let sine = sqrt(self.cross(other).length_squared() / a / b);
        Ok(if cosine < 0.0 {
            PI - asin(sine)
        } else {
            asin(sine)
        })
    }

    /// `IsParallelTo`: within `degrees` of parallel or antiparallel. Returns
    /// whether parallel and, when both vectors are nonzero, whether antiparallel.
    pub fn is_parallel_to(self, other: Self, degrees: f64) -> (bool, Option<bool>) {
        if self.length_squared() > EFF_ZERO_SQ && other.length_squared() > EFF_ZERO_SQ {
            let mut angle = self.angle_between(other).unwrap_or(0.0);
            let anti = angle > PI / 2.0;
            if anti {
                angle = PI - angle;
            }
            (angle < degrees * (PI / 180.0), Some(anti))
        } else {
            (false, None)
        }
    }

    /// `IsPerpendicularTo`: within `degrees` of a right angle; zero vectors are not.
    pub fn is_perpendicular_to(self, other: Self, degrees: f64) -> bool {
        if self.length_squared() > EFF_ZERO_SQ && other.length_squared() > EFF_ZERO_SQ {
            let angle = self.angle_between(other).unwrap_or(0.0);
            (angle - PI / 2.0).abs() < degrees * (PI / 180.0)
        } else {
            false
        }
    }

    /// `MakeUnitOrthoVectors`: a right-handed unit frame with x along `self`
    /// and y in the plane of `self` and `reference` when that plane exists.
    pub fn make_unit_ortho_vectors(self, reference: Option<Self>) -> Result<[Self; 3], Fault> {
        let mut x = self;
        x.unitize()?;
        if let Some(reference) = reference {
            let mut z = x.cross(reference);
            if z.length_squared() > EFF_ZERO_SQ {
                let _ = z.unitize();
                let mut y = z.cross(x);
                let _ = y.unitize();
                return Ok([x, y, z]);
            }
        }
        let mut z = x.cross(Self::new(0.0, 1.0, 0.0));
        if z.length_squared() < EFF_ZERO_SQ {
            z = x.cross(Self::new(0.0, 0.0, 1.0));
            if z.length_squared() < EFF_ZERO_SQ {
                return Err(Fault::Degenerate);
            }
        }
        let _ = z.unitize();
        let mut y = z.cross(x);
        let _ = y.unitize();
        Ok([x, y, z])
    }

    /// `UnitizedDerivative`: the derivative of `self / |self|` given the
    /// derivative of `self`; a negligible `self` returns the derivative.
    pub fn unitized_derivative(self, derivative: Self) -> Self {
        let magnitude = self.length();
        if magnitude < EFF_ZERO {
            return derivative;
        }
        let unit = self / magnitude;
        let magnitude_t = unit.dot(derivative);
        let mut output = derivative - magnitude_t * unit;
        output /= magnitude;
        output
    }
}

impl Neg for Vector3 {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl Add for Vector3 {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self::new(self.x + b.x, self.y + b.y, self.z + b.z)
    }
}

impl Sub for Vector3 {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self::new(self.x - b.x, self.y - b.y, self.z - b.z)
    }
}

impl Mul<f64> for Vector3 {
    type Output = Self;
    fn mul(self, b: f64) -> Self {
        Self::new(self.x * b, self.y * b, self.z * b)
    }
}

impl Mul<Vector3> for f64 {
    type Output = Vector3;
    fn mul(self, b: Vector3) -> Vector3 {
        b * self
    }
}

impl Div<f64> for Vector3 {
    type Output = Self;
    /// Divides unconditionally, as SMLib's quotient operator does.
    fn div(self, b: f64) -> Self {
        Self::new(self.x / b, self.y / b, self.z / b)
    }
}

impl MulAssign<f64> for Vector3 {
    fn mul_assign(&mut self, b: f64) {
        self.x *= b;
        self.y *= b;
        self.z *= b;
    }
}

impl DivAssign<f64> for Vector3 {
    /// Leaves `self` unchanged when the divisor is negligible against every coordinate.
    fn div_assign(&mut self, b: f64) {
        let check = b.abs();
        if check <= EFF_ZERO * self.x.abs()
            && check <= EFF_ZERO * self.y.abs()
            && check <= EFF_ZERO * self.z.abs()
        {
            return;
        }
        self.x /= b;
        self.y /= b;
        self.z /= b;
    }
}
