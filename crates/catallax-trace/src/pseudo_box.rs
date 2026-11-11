// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: only the operations the intersection chain enters;
// `AreDisjoint` without its projected form, which the chain never asks for;
// statuses become `Fault` and output references become returned values.
// Adapted from rust-brep crates/kernel/src/sm_pseudo_box.rs, commit
// 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab, a translation of SMLib
// SmPseudoBox.cpp.
//! SMLib's `SmPseudoBox`: a parallelepiped given by three basis vectors and
//! an interval along each. It bounds a Bezier span or patch more tightly
//! than an axis-aligned box.

use crate::Fault;
use crate::extent::Extent1;
use crate::geometry::{
    PlanePair, determinant_3_vectors, intersect_three_planes, intersect_two_planes,
    line_closest_point,
};
use crate::scalar::{EFF_ZERO, UNDEFINED, is_zero, sqrt};
use crate::vector::Vector3;

/// `SmPseudoBox`: the points whose coordinates in `basis` lie inside the
/// three `intervals`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PseudoBox {
    /// Whether the basis vectors are mutually perpendicular.
    pub orthogonal_basis: bool,
    /// The three unit basis vectors.
    pub basis: [Vector3; 3],
    /// The coordinate interval along each basis vector.
    pub intervals: [Extent1; 3],
}

impl Default for PseudoBox {
    /// The empty box on the XYZ basis.
    fn default() -> Self {
        Self::new()
    }
}

impl PseudoBox {
    /// `SmPseudoBox()`: the XYZ basis and three empty intervals.
    pub fn new() -> Self {
        Self {
            orthogonal_basis: true,
            basis: [
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.0, 1.0, 0.0),
                Vector3::new(0.0, 0.0, 1.0),
            ],
            intervals: [Extent1::default(); 3],
        }
    }

    /// `Init`: empties the intervals and keeps the basis.
    pub fn init(&mut self) {
        self.intervals = [Extent1::default(); 3];
    }

    /// `InitBasis`: the XYZ basis.
    pub fn init_basis(&mut self) {
        self.orthogonal_basis = true;
        self.basis = Self::new().basis;
    }

    /// `AssertDefined`: the intervals are not destroyed and the basis is set.
    fn is_defined(&self) -> bool {
        self.intervals
            .iter()
            .all(|interval| interval.min != UNDEFINED && interval.max != -UNDEFINED)
            && self.basis.iter().all(|b| b.is_initialized())
    }

    /// `InvertAbsolute`: the coordinates of `point` in the basis, by Cramer's rule.
    pub fn invert_absolute(&self, point: Vector3) -> Result<Vector3, Fault> {
        if !self.is_defined() {
            return Err(Fault::Degenerate);
        }
        let [b0, b1, b2] = self.basis;
        let det = determinant_3_vectors(b0, b1, b2);
        if is_zero(det) {
            return Err(Fault::Degenerate);
        }
        Ok(Vector3::new(
            determinant_3_vectors(point, b1, b2) / det,
            determinant_3_vectors(b0, point, b2) / det,
            determinant_3_vectors(b0, b1, point) / det,
        ))
    }

    /// `AddPoint3d`: grows the intervals to include `point`. A singular basis
    /// adds SMLib's uninitialised coordinates, as the C++ does.
    pub fn add_point(&mut self, point: Vector3) {
        let a = self.invert_absolute(point).unwrap_or_default();
        self.intervals[0].add_value(a.x);
        self.intervals[1].add_value(a.y);
        self.intervals[2].add_value(a.z);
    }

    /// `ExpandAbsolute`.
    pub fn expand_absolute(&mut self, expansion: f64) {
        for interval in &mut self.intervals {
            interval.expand_absolute(expansion);
        }
    }

    /// `SetBasis`: stores the three directions normalised; refuses a
    /// singular set.
    pub fn set_basis(
        &mut self,
        basis1: Vector3,
        basis2: Vector3,
        basis3: Vector3,
    ) -> Result<(), Fault> {
        if determinant_3_vectors(basis1, basis2, basis3).abs() < EFF_ZERO {
            return Err(Fault::InvalidInput);
        }
        self.basis = [basis1, basis2, basis3];
        for b in &mut self.basis {
            let length_squared = b.length_squared();
            if !is_zero(length_squared - 1.0) {
                *b /= sqrt(length_squared);
            }
        }
        let [b0, b1, b2] = self.basis;
        self.orthogonal_basis = b0.dot(b1).abs() < EFF_ZERO
            && b1.dot(b2).abs() < EFF_ZERO
            && b2.dot(b0).abs() < EFF_ZERO;
        Ok(())
    }

    /// `GetVolume`: interval lengths times the basis triple product.
    pub fn volume(&self) -> f64 {
        let lengths =
            self.intervals[0].length() * self.intervals[1].length() * self.intervals[2].length();
        let [b0, b1, b2] = self.basis;
        let triple = b0.x * (b1.y * b2.z - b1.z * b2.y)
            + b0.y * (b1.z * b2.x - b1.x * b2.z)
            + b0.z * (b1.x * b2.y - b1.y * b2.x);
        lengths * triple
    }

    /// `CalcPlanes`: a point and outward unit normal for each face, in the
    /// order u-min, u-max, v-min, v-max, w-min, w-max.
    pub fn planes(&self) -> ([Vector3; 6], [Vector3; 6]) {
        let [b0, b1, b2] = self.basis;
        let no_flip = determinant_3_vectors(b0, b1, b2) > 0.0;
        let origin =
            b0 * self.intervals[0].min + b1 * self.intervals[1].min + b2 * self.intervals[2].min;
        let edges = [
            b0 * self.intervals[0].length(),
            b1 * self.intervals[1].length(),
            b2 * self.intervals[2].length(),
        ];
        let mut points = [Vector3::default(); 6];
        let mut normals = [Vector3::default(); 6];
        for axis in 0..3 {
            let (next, last) = ((axis + 1) % 3, (axis + 2) % 3);
            let mut normal = self.basis[next].cross(self.basis[last]);
            let _ = normal.unitize();
            points[2 * axis] = origin + 0.5 * (edges[next] + edges[last]);
            normals[2 * axis] = if no_flip { -normal } else { normal };
            points[2 * axis + 1] = points[2 * axis] + edges[axis];
            normals[2 * axis + 1] = -normals[2 * axis];
        }
        (points, normals)
    }

    /// `CalcCorners`: the eight corners, bottom face counter-clockwise first.
    pub fn corners(&self) -> [Vector3; 8] {
        let [b0, b1, b2] = self.basis;
        let c0 =
            b0 * self.intervals[0].min + b1 * self.intervals[1].min + b2 * self.intervals[2].min;
        let del0 = b0 * self.intervals[0].length();
        let del1 = b1 * self.intervals[1].length();
        let del2 = b2 * self.intervals[2].length();
        let c1 = c0 + del0;
        let c2 = c1 + del1;
        let c3 = c0 + del1;
        [c0, c1, c2, c3, c0 + del2, c1 + del2, c2 + del2, c3 + del2]
    }

    /// `DistanceToPoint`: 0 inside; the face, edge or corner distance when
    /// outside one, two or three face planes; SMLib's `SM_ERR` value, 1001,
    /// beyond that.
    pub fn distance_to_point(&self, point: Vector3) -> f64 {
        let (points, normals) = self.planes();
        let mut outside = Vec::with_capacity(6);
        for i in 0..6 {
            let d = -points[i].dot(normals[i]);
            let height = normals[i].dot(point) + d;
            if height > EFF_ZERO {
                outside.push((i, height));
            }
        }
        match *outside.as_slice() {
            [] => 0.0,
            [(_, height)] => height,
            [(a, _), (b, _)] => {
                let (line_point, line_vector) =
                    match intersect_two_planes(points[a], normals[a], points[b], normals[b]) {
                        Ok(PlanePair::Line { point, direction }) => (point, direction),
                        _ => (Vector3::default(), Vector3::default()),
                    };
                let parameter = line_closest_point(line_point, line_vector, point);
                (line_point + parameter * line_vector).distance(point)
            }
            [(a, _), (b, _), (c, _)] => {
                let vertex = intersect_three_planes(
                    points[a], normals[a], points[b], normals[b], points[c], normals[c],
                )
                .unwrap_or_default();
                vertex.distance(point)
            }
            _ => 1001.0,
        }
    }

    /// Whether no corner lies less than `SM_EFF_ZERO` outside one face plane.
    // Negated as in SMLib, so a NaN height does not count as inside.
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    fn outside_one_plane(
        points: &[Vector3; 6],
        normals: &[Vector3; 6],
        corners: &[Vector3; 8],
    ) -> bool {
        (0..6).any(|i| {
            let d = -points[i].dot(normals[i]);
            corners
                .iter()
                .all(|corner| !(normals[i].dot(*corner) + d < EFF_ZERO))
        })
    }

    /// `AreDisjoint` without a projection: some face plane of one box
    /// separates it from all corners of the other.
    pub fn are_disjoint(&self, other: &Self) -> bool {
        let (points, normals) = self.planes();
        if Self::outside_one_plane(&points, &normals, &other.corners()) {
            return true;
        }
        let (points, normals) = other.planes();
        Self::outside_one_plane(&points, &normals, &self.corners())
    }
}
