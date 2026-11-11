// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: translated from SMLib SmExtent3d.cpp
// (MinimumDistanceSquared, MaximumDistanceSquared, MinimumDistance) and
// SmVector3d.h (RotateProjectPointToPlane,
// PerspectiveProjectToPlane, PerspectiveProjectPointToPlane), the helpers the
// global search solver's box tests call. Assertion macros are dropped.
//! Box-distance and projection helpers for the global search solver.

use crate::extent::Extent3;
use crate::scalar::sqrt;
use crate::vector::Vector3;

/// `SmExtent3d::MaximumDistanceSquared`: the largest squared distance between
/// two points of the union of the boxes.
pub(crate) fn maximum_distance_squared(a: Extent3, b: Extent3) -> f64 {
    let union = a.union(b);
    (union.max - union.min).length_squared()
}

/// `SmExtent3d::MinimumDistanceSquared`: the smallest squared distance
/// between a point of each box, zero when they overlap.
pub(crate) fn minimum_distance_squared(a: Extent3, b: Extent3) -> f64 {
    let mut gap = Vector3::new(0.0, 0.0, 0.0);
    if a.min.x > b.max.x {
        gap.x += a.min.x - b.max.x;
    }
    if a.max.x < b.min.x {
        gap.x += b.min.x - a.max.x;
    }
    if a.min.y > b.max.y {
        gap.y += a.min.y - b.max.y;
    }
    if a.max.y < b.min.y {
        gap.y += b.min.y - a.max.y;
    }
    if a.min.z > b.max.z {
        gap.z += a.min.z - b.max.z;
    }
    if a.max.z < b.min.z {
        gap.z += b.min.z - a.max.z;
    }
    gap.length_squared()
}

/// `SmExtent3d::MinimumDistance`.
pub(crate) fn minimum_distance(a: Extent3, b: Extent3) -> f64 {
    sqrt(minimum_distance_squared(a, b))
}

/// `SmVector3d::RotateProjectPointToPlane`: the point rotated about the unit
/// axis through `rotate_point` into the plane spanned by the axis and the
/// unit `x_axis`.
pub(crate) fn rotate_project_point_to_plane(
    point: Vector3,
    rotate_point: Vector3,
    rotate_axis: Vector3,
    x_axis: Vector3,
) -> Vector3 {
    let diff = point - rotate_point;
    rotate_point + rotate_axis * diff.dot(rotate_axis) + x_axis * diff.dot(x_axis)
}

/// `SmVector3d::PerspectiveProjectToPlane`: the vector `vector` based at
/// `base`, projected towards `eye` onto the plane through `plane_point` with
/// unit normal `normal`.
pub(crate) fn perspective_project_vector_to_plane(
    vector: Vector3,
    base: Vector3,
    plane_point: Vector3,
    normal: Vector3,
    eye: Vector3,
) -> Vector3 {
    let eye_plane = plane_point - eye;
    let eye_point = base - eye;
    let eye_vector = base + vector - eye;
    eye_vector * (eye_plane.dot(normal) / eye_vector.dot(normal))
        - eye_point * (eye_plane.dot(normal) / eye_point.dot(normal))
}

/// `SmVector3d::PerspectiveProjectPointToPlane`: the point projected towards
/// `eye` onto the plane through `plane_point` with unit normal `normal`.
pub(crate) fn perspective_project_point_to_plane(
    point: Vector3,
    plane_point: Vector3,
    normal: Vector3,
    eye: Vector3,
) -> Vector3 {
    let eye_plane = plane_point - eye;
    let eye_point = point - eye;
    eye + eye_point * (eye_plane.dot(normal) / eye_point.dot(normal))
}
