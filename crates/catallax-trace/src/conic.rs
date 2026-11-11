// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: a plain value type in the plane; statuses become
// returned values, as every method used here always succeeds in SMLib.
// Translated from SMLib SmConic.cpp (the SmConic constructor,
// GetNaturalInterval, EvaluatePoint, ComputeImplicitEquation,
// IntersectImplicitConic, IntersectConic, InvertPointOnConic, Transform2D);
// rust-brep has no port.
//! The planar conics of the Dupin indicatrix (`SmConic`): rotated and moved
//! ellipses, hyperbola branches, parabolas and lines with SMLib's rational
//! parameterisations, and their analytic intersection.

use crate::extent::Extent1;
use crate::geometry::solve_quadratic;
use crate::polynomial::solve_quartic;
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQRT, cos, sin};
use crate::vector::Vector2;

/// `SmConicType`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConicKind {
    /// `x^2/a^2 + y^2/b^2 = 1`, parameter `[-1, 3]` (two halves).
    Ellipse,
    /// The right branch of `x^2/a^2 - y^2/b^2 = 1`, parameter inside
    /// `(-1, 1)`.
    Hyperbola,
    /// `y^2 = 4 a x`.
    Parabola,
    /// `y = a`.
    Line,
}

/// `SmConic`: a conic in natural position, then rotated by `theta`
/// (clockwise, as SMLib stores it) and moved to `origin`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Conic {
    /// The kind.
    pub kind: ConicKind,
    /// `m_dA`.
    pub a: f64,
    /// `m_dB`.
    pub b: f64,
    /// `m_vOrigin`.
    pub origin: Vector2,
    /// `m_dThetaRad`.
    pub theta: f64,
    sin_theta: f64,
    cos_theta: f64,
}

impl Conic {
    /// `SmConic(type, a, b)`: in natural position.
    pub fn new(kind: ConicKind, a: f64, b: f64) -> Self {
        Self {
            kind,
            a,
            b,
            origin: Vector2::new(0.0, 0.0),
            theta: 0.0,
            sin_theta: 0.0,
            cos_theta: 1.0,
        }
    }

    /// `GetNaturalInterval`.
    pub fn natural_interval(&self) -> Extent1 {
        match self.kind {
            ConicKind::Parabola | ConicKind::Line => Extent1::new(-999.0, 999.0),
            ConicKind::Ellipse => Extent1::new(-1.0, 3.0),
            ConicKind::Hyperbola => Extent1::new(-1.0 + EFF_ZERO_SQRT, 1.0 - EFF_ZERO_SQRT),
        }
    }

    /// `EvaluatePoint`.
    pub fn evaluate_point(&self, t: f64) -> Vector2 {
        let mut t = t;
        let mut t_sq = t * t;
        let mut point = match self.kind {
            ConicKind::Ellipse => {
                let mut sign = 1.0;
                if t > 1.0 {
                    t -= 2.0;
                    t_sq = t * t;
                    sign = -1.0;
                }
                let denominator = 1.0 + t_sq;
                Vector2::new(
                    sign * self.a * (1.0 - t_sq) / denominator,
                    sign * self.b * (2.0 * t) / denominator,
                )
            }
            ConicKind::Hyperbola => {
                let denominator = 1.0 - t_sq;
                Vector2::new(
                    self.a * (1.0 + t_sq) / denominator,
                    self.b * (2.0 * t) / denominator,
                )
            }
            ConicKind::Parabola => Vector2::new(self.a * t_sq, 2.0 * self.a * t),
            ConicKind::Line => Vector2::new(t, self.a),
        };
        if self.theta != 0.0 {
            point = Vector2::new(
                point.x * self.cos_theta + point.y * self.sin_theta,
                point.y * self.cos_theta - point.x * self.sin_theta,
            );
        }
        if self.origin.x != 0.0 || self.origin.y != 0.0 {
            point = point + self.origin;
        }
        point
    }

    /// `Transform2D`: rotates counter-clockwise by `angle` about the
    /// origin, then translates.
    pub fn transform_2d(&mut self, angle: f64, translation: Vector2) {
        self.theta -= angle;
        let (c, s) = (cos(-angle), sin(-angle));
        let rotated = Vector2::new(
            self.origin.x * c + self.origin.y * s,
            self.origin.y * c - self.origin.x * s,
        );
        self.origin = rotated + translation;
        self.cos_theta = cos(self.theta);
        self.sin_theta = sin(self.theta);
    }

    /// `ComputeImplicitEquation`: `A x^2 + B xy + C y^2 + D x + E y + F`.
    pub fn implicit_equation(&self) -> [f64; 6] {
        let (a, b, c, d, e, f) = match self.kind {
            ConicKind::Ellipse => (
                1.0 / (self.a * self.a),
                0.0,
                1.0 / (self.b * self.b),
                0.0,
                0.0,
                -1.0,
            ),
            ConicKind::Hyperbola => (
                1.0 / (self.a * self.a),
                0.0,
                -1.0 / (self.b * self.b),
                0.0,
                0.0,
                -1.0,
            ),
            ConicKind::Parabola => (0.0, 0.0, 1.0, -4.0 * self.a, 0.0, 0.0),
            ConicKind::Line => (0.0, 0.0, 0.0, 0.0, -1.0 / self.a, 1.0),
        };
        let s = sin(-self.theta);
        let co = cos(-self.theta);
        let (sin_sq, cos_sq, sin_cos) = (s * s, co * co, s * co);
        let a2 = a * cos_sq - b * sin_cos + c * sin_sq;
        let b2 = 2.0 * a * sin_cos + b * cos_sq - b * sin_sq - 2.0 * c * sin_cos;
        let c2 = a * sin_sq + b * sin_cos + c * cos_sq;
        let d2 = d * co - e * s;
        let e2 = d * s + e * co;
        let f2 = f;
        let (c1, c_2) = (-self.origin.x, -self.origin.y);
        [
            a2,
            b2,
            c2,
            2.0 * a2 * c1 + b2 * c_2 + d2,
            b2 * c1 + 2.0 * c2 * c_2 + e2,
            a2 * c1 * c1 + b2 * c1 * c_2 + c2 * c_2 * c_2 + d2 * c1 + e2 * c_2 + f2,
        ]
    }

    /// `IntersectImplicitConic`: the parameters of this conic (taken in
    /// natural position) on the implicit conic `coefficients`, from the
    /// quartic in the rational parameter, each half of an ellipse in turn.
    pub fn intersect_implicit(&self, coefficients: [f64; 6]) -> Vec<f64> {
        let [a, b, c, d, e, f] = coefficients;
        let mut found = Vec::new();
        for side in 0..2 {
            if side == 1 && self.kind != ConicKind::Ellipse {
                break;
            }
            let sign = if side == 1 { -1.0 } else { 1.0 };
            let ka = self.a * sign;
            let a_sq = ka * ka;
            let kb = self.b * sign;
            let b_sq = kb * kb;
            let q = match self.kind {
                ConicKind::Ellipse => [
                    a * a_sq + d * ka + f,
                    2.0 * (b * ka * kb + e * kb),
                    2.0 * (-a * a_sq + 2.0 * c * b_sq + f),
                    2.0 * (-b * ka * kb + e * kb),
                    a * a_sq - d * ka + f,
                ],
                ConicKind::Hyperbola => [
                    a * a_sq + d * ka + f,
                    2.0 * (b * ka * kb + e * kb),
                    2.0 * (a * a_sq + 2.0 * c * b_sq - f),
                    2.0 * (b * ka * kb - e * kb),
                    a * a_sq - d * ka + f,
                ],
                ConicKind::Parabola => [
                    a * a_sq,
                    2.0 * (b * a_sq),
                    4.0 * c * a_sq + d * ka,
                    2.0 * e * ka,
                    f,
                ],
                ConicKind::Line => [c * a_sq + e * ka + f, b * ka + d, a, 0.0, 0.0],
            };
            let (count, roots) = solve_quartic(q, EFF_ZERO);
            for &root in &roots[..count] {
                let mut root = root;
                match self.kind {
                    ConicKind::Ellipse => {
                        if root.abs() <= 1.0 + EFF_ZERO {
                            if sign < 0.0 {
                                if root.abs() < 1.0 {
                                    found.push(2.0 + root);
                                }
                            } else {
                                found.push(root);
                            }
                        }
                    }
                    ConicKind::Hyperbola => {
                        let interval = self.natural_interval();
                        if root.abs() < 1.0 + EFF_ZERO_SQRT {
                            root = interval.clamp_value(root);
                        }
                        if interval.contains_value(root, 0.0) {
                            found.push(root);
                        }
                    }
                    ConicKind::Parabola | ConicKind::Line => {
                        if self.natural_interval().contains_value(root, 0.0) {
                            found.push(root);
                        }
                    }
                }
            }
        }
        found
    }

    /// `IntersectConic`: the parameter pairs `(t, other t, distance)` where
    /// the conics meet within `tolerance`.
    pub fn intersect_conic(&self, other: &Self, tolerance: f64) -> Vec<(f64, f64, f64)> {
        let mut moved = Self::new(other.kind, other.a, other.b);
        moved.transform_2d(-other.theta, other.origin);
        moved.transform_2d(0.0, -self.origin);
        moved.transform_2d(self.theta, Vector2::new(0.0, 0.0));
        let mut found = Vec::new();
        for t in self.intersect_implicit(moved.implicit_equation()) {
            let point = self.evaluate_point(t);
            let Some(t2) = other.invert_point(point) else {
                continue;
            };
            let distance = other.evaluate_point(t2).distance(point);
            if distance > tolerance {
                continue;
            }
            found.push((t, t2, distance));
        }
        found
    }

    /// `InvertPointOnConic`: the parameter of a point on the conic, if it
    /// lies within the scaled square-root tolerance.
    pub fn invert_point(&self, point: Vector2) -> Option<f64> {
        let offset = point - self.origin;
        let (c, s) = (cos(-self.theta), sin(-self.theta));
        let local = Vector2::new(offset.x * c + offset.y * s, offset.y * c - offset.x * s);
        let (mut parameter, second) = match self.kind {
            ConicKind::Hyperbola => {
                let (count, roots) = solve_quadratic([-local.y, 2.0 * self.b, local.y], EFF_ZERO);
                if count == 0 {
                    return None;
                }
                let interval = self.natural_interval();
                if roots[0].abs() <= 1.0 + 0.0001 {
                    let second = if count > 1 && roots[1].abs() <= 1.0 + 0.0001 {
                        interval.clamp_value(roots[1])
                    } else {
                        0.0
                    };
                    (interval.clamp_value(roots[0]), second)
                } else if count > 1 && roots[1].abs() <= 1.0 + 0.0001 {
                    (interval.clamp_value(roots[1]), 0.0)
                } else {
                    return None;
                }
            }
            ConicKind::Ellipse => {
                let (count, roots) = solve_quadratic([local.y, -2.0 * self.b, local.y], EFF_ZERO);
                if count == 0 {
                    return None;
                }
                let mut parameter = if roots[0].abs() <= 1.0 + EFF_ZERO {
                    roots[0]
                } else if count > 1 && roots[1].abs() <= 1.0 + EFF_ZERO {
                    roots[1]
                } else {
                    return None;
                };
                if local.x < 0.0 {
                    parameter = 2.0 - parameter;
                }
                (parameter, 0.0)
            }
            ConicKind::Parabola => (local.y / (2.0 * self.a), 0.0),
            ConicKind::Line => (local.x, 0.0),
        };
        let test = self.evaluate_point(parameter);
        let scale = 1.0 + test.max_dimension();
        if test.distance(point) > EFF_ZERO_SQRT * scale {
            if self.kind != ConicKind::Hyperbola {
                return None;
            }
            parameter = second;
            if self.evaluate_point(parameter).distance(point) > EFF_ZERO_SQRT * scale {
                return None;
            }
        }
        Some(parameter)
    }
}
