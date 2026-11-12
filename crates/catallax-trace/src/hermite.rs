// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: statuses become `Fault` and output references
// become returned values; the curve's dimension, form and the continuity,
// offset and database services are not carried.
// Translated from SMLib SmHermiteCurve.h (the two constructors,
// GetNaturalInterval) and SmHermiteCurve.cpp (CalculateBoundingBox, Evaluate,
// EvaluatePoint) at c9c5979f04cea8c38213004d4f04df0dcfb40d05; rust-brep has
// no port of these.
//! SMLib's `SmHermiteCurve`: the cubic through two points with given
//! end derivatives, over the parameter range `[0, 1]`.

use crate::Fault;
use crate::extent::{Extent1, Extent3};
use crate::polar_box::PolarBox;
use crate::pseudo_box::PseudoBox;
use crate::scalar::{EFF_ZERO, is_zero};
use crate::vector::{Vector2, Vector3};

/// The cubic Hermite curve through `p1` and `p2` with end derivatives `d1`
/// and `d2` (`SmHermiteCurve`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hermite {
    p1: Vector3,
    d1: Vector3,
    p2: Vector3,
    d2: Vector3,
    c: Vector3,
    d: Vector3,
}

impl Hermite {
    /// The three-dimensional constructor. It stores the power-basis
    /// coefficients that the evaluators use.
    pub fn new(p1: Vector3, d1: Vector3, p2: Vector3, d2: Vector3) -> Self {
        Self {
            p1,
            d1,
            p2,
            d2,
            c: -3.0 * p1 - 2.0 * d1 + 3.0 * p2 - d2,
            d: 2.0 * p1 + d1 - 2.0 * p2 + d2,
        }
    }

    /// The two-dimensional constructor: the curve in the plane `z = 0`.
    pub fn planar(p1: Vector2, d1: Vector2, p2: Vector2, d2: Vector2) -> Self {
        let lift = |v: Vector2| Vector3::new(v.x, v.y, 0.0);
        Self::new(lift(p1), lift(d1), lift(p2), lift(d2))
    }

    /// `GetNaturalInterval`: always `[0, 1]`.
    pub fn natural_interval(&self) -> Extent1 {
        Extent1::new(0.0, 1.0)
    }

    /// `Evaluate`: the point at `u` and `derivatives` derivatives, up to the
    /// third; the fourth and beyond are [`Fault::InvalidInput`], where SMLib
    /// leaves them unset. With `nonzero_tangent`, a first derivative shorter
    /// than a scaled `SM_EFF_ZERO` takes the direction of the second.
    pub fn evaluate(
        &self,
        u: f64,
        derivatives: usize,
        nonzero_tangent: bool,
    ) -> Result<Vec<Vector3>, Fault> {
        if derivatives > 3 {
            return Err(Fault::InvalidInput);
        }
        let (a, b, c, d) = (self.p1, self.d1, self.c, self.d);
        let mut values = vec![Vector3::new(0.0, 0.0, 0.0); derivatives + 1];
        values[0] = a + u * (b + u * (c + u * d));
        if derivatives >= 1 {
            values[1] = b + u * (2.0 * c + 3.0 * u * d);
        }
        if derivatives >= 2 {
            values[2] = 2.0 * c + 6.0 * u * d;
        }
        if derivatives >= 3 {
            values[3] = 6.0 * d;
        }
        if nonzero_tangent {
            let new_length = 1.1 * (EFF_ZERO * (1.0 + values[0].max_dimension()));
            if derivatives >= 1 && values[1].length_squared() < new_length * new_length {
                if derivatives < 2 {
                    // The recursive call keeps the default of `TRUE`.
                    values[1] = self.evaluate(u, 2, true)?[1];
                } else {
                    let length = values[2].length();
                    if length > new_length {
                        values[1] = (new_length / length) * values[2];
                        if u > self.natural_interval().mid() {
                            values[1] *= -1.0;
                        }
                    }
                }
            }
        }
        Ok(values)
    }

    /// `EvaluatePoint`: the point at `u`.
    pub fn evaluate_point(&self, u: f64) -> Vector3 {
        self.p1 + u * (self.d1 + u * (self.c + u * self.d))
    }

    /// `CalculateBoundingBox`: the axis-aligned box and the pseudo box of the
    /// four Bezier control points, whatever the interval. The pseudo box
    /// points along the chord, with its second axis toward the start
    /// derivative.
    pub fn bounding_boxes(&self) -> (Extent3, PseudoBox) {
        let points = [
            self.p1,
            self.p1 + self.d1 / 3.0,
            self.p2 - self.d2 / 3.0,
            self.p2,
        ];
        let mut pseudo = PseudoBox::new();
        let chord = points[3] - points[0];
        if !is_zero(chord.length_squared())
            && let Ok([x, y, z]) = chord.make_unit_ortho_vectors(Some(self.d1))
        {
            let _ = pseudo.set_basis(x, y, z);
        }
        let mut bbox = Extent3::singleton(points[0]);
        for &point in &points[1..] {
            bbox.add_point(point);
        }
        for &point in &points {
            pseudo.add_point(point);
        }
        (bbox, pseudo)
    }

    /// The polar part of `CalculateBoundingBox`: the directions of the three
    /// Bezier control-polygon edges, over the whole curve regardless of the
    /// requested interval. This is an unproved candidate bound.
    pub fn polar_box(&self) -> PolarBox {
        let points = [
            self.p1,
            self.p1 + self.d1 / 3.0,
            self.p2 - self.d2 / 3.0,
            self.p2,
        ];
        let mut polar = PolarBox::default();
        for edge in points.windows(2) {
            polar.add_vector(edge[1] - edge[0], None);
        }
        polar
    }
}
