// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: results are returned values and `Fault`s instead
// of status codes and output references; calls that always succeed in SMLib
// return plain values. `evaluate_curve_continuity` is translated directly
// from smgu_EvaluateCurveContinuity.
// Adapted from rust-brep crates/kernel/src/sm_geom_utility.rs and
// sm_geom_intersect.rs, commit 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab,
// translations of SMLib SmGeomUtility.cpp.
//! The `smgu_*` point, line, plane, triangle and circle helpers that the
//! intersection chain enters (step 1 of the M25 port order).

use crate::Fault;
use crate::extent::Extent1;
use crate::scalar::{
    EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, PI, ZONE_TOL_3D, degrees, is_zero, max, sqrt,
};
use crate::vector::{Vector2, Vector3};

/// `smgu_TriangleAreaSquared`.
pub fn triangle_area_squared(p0: Vector3, p1: Vector3, p2: Vector3) -> f64 {
    let ay_bz = (p1.y - p0.y) * (p2.z - p0.z);
    let az_bx = (p1.z - p0.z) * (p2.x - p0.x);
    let ax_by = (p1.x - p0.x) * (p2.y - p0.y);
    let az_by = (p1.z - p0.z) * (p2.y - p0.y);
    let ax_bz = (p1.x - p0.x) * (p2.z - p0.z);
    let ay_bx = (p1.y - p0.y) * (p2.x - p0.x);
    ((ay_bz - az_by) * (ay_bz - az_by)
        + (az_bx - ax_bz) * (az_bx - ax_bz)
        + (ax_by - ay_bx) * (ax_by - ay_bx))
        / 4.0
}

/// The nearest point of a triangle (`smgu_TrianglePointDistance`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TrianglePoint {
    /// Distance from the test point to `closest`.
    pub distance: f64,
    /// The nearest point found.
    pub closest: Vector3,
    /// Parameter along `p1 - p0`.
    pub u: f64,
    /// Parameter along `p2 - p0`.
    pub v: f64,
}

/// `smgu_TrianglePointDistance`. A degenerate triangle falls back to the
/// nearer of its two sides from `p0`.
pub fn triangle_point_distance(
    p0: Vector3,
    p1: Vector3,
    p2: Vector3,
    test: Vector3,
) -> TrianglePoint {
    let x_side = p1 - p0;
    let y_side = p2 - p0;
    let Ok(uv) = vec_lin_comb_two_vectors(x_side, y_side, test - p0) else {
        let (distance1, u) = segment_point_distance(p0, p1, test);
        let (distance2, v) = segment_point_distance(p0, p2, test);
        return if distance1 <= distance2 {
            TrianglePoint {
                distance: distance1,
                closest: p0 + u * x_side,
                u,
                v,
            }
        } else {
            TrianglePoint {
                distance: distance2,
                closest: p0 + v * y_side,
                u,
                v,
            }
        };
    };
    let (mut u, mut v) = (uv.x, uv.y);
    if u <= 0.0 || v <= 0.0 {
        u = u.clamp(0.0, 1.0);
        v = v.clamp(0.0, 1.0);
    } else if u + v > 1.0 {
        v = segment_point_distance(p1, p2, test).1;
        u = 1.0 - v;
    }
    let closest = p0 + u * x_side + v * y_side;
    TrianglePoint {
        distance: test.distance(closest),
        closest,
        u,
        v,
    }
}

/// `smgu_LinePointDistance`: the distance from `test` to the line and the
/// line parameter of the foot. A negligible `vector` gives parameter 0.
pub fn line_point_distance(point: Vector3, vector: Vector3, test: Vector3) -> (f64, f64) {
    let parameter = line_closest_point(point, vector, test);
    let foot = point + parameter * vector;
    (foot.distance(test), parameter)
}

/// `smgu_LineClosestPoint`: the line parameter nearest `test`.
pub fn line_closest_point(point: Vector3, vector: Vector3, test: Vector3) -> f64 {
    let length_squared = vector.length_squared();
    if length_squared < EFF_ZERO_SQ {
        0.0
    } else {
        (test - point).dot(vector) / length_squared
    }
}

/// `smgu_SegmentPointDistance`: the distance from `test` to the segment and
/// the clamped segment parameter in `[0, 1]`.
pub fn segment_point_distance(start: Vector3, end: Vector3, test: Vector3) -> (f64, f64) {
    let vector = end - start;
    let parameter = line_closest_point(start, vector, test).clamp(0.0, 1.0);
    let found = start + parameter * vector;
    (found.distance(test), parameter)
}

/// `smgu_SegmentPointAngleDegToPerp`: how far, in degrees, the gap from the
/// nearest segment point to `test` is from perpendicular to the segment,
/// with the clamped segment parameter. Zero when the foot is interior or
/// `test` is on the segment.
pub fn segment_point_angle_deg_to_perp(start: Vector3, end: Vector3, test: Vector3) -> (f64, f64) {
    let vector = end - start;
    let parameter = line_closest_point(start, vector, test);
    let clip = Extent1::new(0.0, 1.0);
    if clip.contains_value(parameter, 0.0) {
        return (0.0, parameter);
    }
    let parameter = clip.clamp_value(parameter);
    let gap = test - (start + vector * parameter);
    if gap.length_squared() < EFF_ZERO_SQ {
        return (0.0, parameter);
    }
    match gap.angle_between(vector) {
        Ok(angle) => ((angle * 180.0 / PI - 90.0).abs(), parameter),
        Err(_) => (0.0, parameter),
    }
}

/// Up to two line parameters (`smgu_LineCoPlanarCircleIntersect`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LineHits {
    /// How many of `parameters` are meaningful: 0, 1 or 2.
    pub count: usize,
    /// Line parameters in increasing order.
    pub parameters: [f64; 2],
}

/// `smgu_LineCoPlanarCircleIntersect`: where a line in the circle's plane
/// meets the circle, in the line's own parameter. Within `tolerance` of
/// tangency the fuzzy band resolves to one touch or two crossings by which
/// is closer to the radius. Fails when the line is perpendicular to the plane.
pub fn line_coplanar_circle_intersect(
    line_point: Vector3,
    line_vector: Vector3,
    center: Vector3,
    normal: Vector3,
    radius: f64,
    tolerance: f64,
) -> Result<LineHits, Fault> {
    let mut hits = LineHits {
        count: 0,
        parameters: [0.0; 2],
    };
    let point = point_project_to_plane(line_point, center, normal);
    let mut vector = vector_project_to_plane(line_vector, normal);
    let scale = vector.length();
    if scale < EFF_ZERO * radius {
        return Err(Fault::Degenerate);
    }
    vector /= scale;
    let t = line_closest_point(point, vector, center);
    let to_foot = point + t * vector - center;
    let distance = sqrt(to_foot.dot(to_foot));
    if distance > radius + tolerance {
        return Ok(hits);
    }
    if distance >= radius {
        hits.count = 1;
        hits.parameters[0] = t / scale;
        return Ok(hits);
    }
    let half_chord = sqrt((radius + distance) * (radius - distance));
    let two = LineHits {
        count: 2,
        parameters: [(t - half_chord) / scale, (t + half_chord) / scale],
    };
    let one = LineHits {
        count: 1,
        parameters: [t / scale, 0.0],
    };
    if distance <= radius - tolerance {
        return Ok(two);
    }
    if half_chord * 2.0 <= tolerance {
        return Ok(one);
    }
    // Fuzzy band: keep whichever answer lies closer to the circle.
    let crossing = point + (t - half_chord) * vector;
    if (radius - distance).abs() < (radius - crossing.distance(center)).abs() {
        Ok(one)
    } else {
        Ok(two)
    }
}

/// The meeting of two planes (`smgu_IntersectTwoPlanes`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PlanePair {
    /// The planes cross in the line `point + s * direction`, with a unit direction.
    Line {
        /// A point on both planes.
        point: Vector3,
        /// Unit direction of the line.
        direction: Vector3,
    },
    /// The planes are parallel; `coincident` when each base point is within
    /// the scaled zero of the other plane's positive side.
    Parallel {
        /// SMLib's coincidence test.
        coincident: bool,
    },
}

/// `smgu_IntersectTwoPlanes`. Normals must be unit length.
pub fn intersect_two_planes(
    point1: Vector3,
    normal1: Vector3,
    point2: Vector3,
    normal2: Vector3,
) -> Result<PlanePair, Fault> {
    if (normal1.dot(normal1) - 1.0).abs() > EFF_ZERO
        || (normal2.dot(normal2) - 1.0).abs() > EFF_ZERO
    {
        return Err(Fault::InvalidInput);
    }
    let middle = (point1 + point2) / 2.0;
    let mut direction = normal1.cross(normal2);
    if direction.length_squared() < EFF_ZERO_SQ {
        let scaled_zero = EFF_ZERO * (1.0 + max(point1.max_dimension(), point2.max_dimension()));
        let distance1 = (point2 - point1).dot(normal1);
        let distance2 = (point1 - point2).dot(normal2);
        return Ok(PlanePair::Parallel {
            coincident: distance1 < scaled_zero && distance2 < scaled_zero,
        });
    }
    direction.unitize()?;
    let point = intersect_three_planes(point1, normal1, point2, normal2, middle, direction)?;
    Ok(PlanePair::Line { point, direction })
}

/// `smgu_IntersectThreePlanes`. Normals must be unit length and independent.
pub fn intersect_three_planes(
    point1: Vector3,
    normal1: Vector3,
    point2: Vector3,
    normal2: Vector3,
    point3: Vector3,
    normal3: Vector3,
) -> Result<Vector3, Fault> {
    if [normal1, normal2, normal3]
        .into_iter()
        .any(|n| (n.dot(n) - 1.0).abs() > EFF_ZERO)
    {
        return Err(Fault::InvalidInput);
    }
    let determinant = determinant_3_vectors(normal1, normal2, normal3);
    if determinant.abs() < EFF_ZERO {
        return Err(Fault::InvalidInput);
    }
    let work = (point1.dot(normal1) * normal2.cross(normal3))
        + (point2.dot(normal2) * normal3.cross(normal1))
        + (point3.dot(normal3) * normal1.cross(normal2));
    Ok(work / determinant)
}

/// `smgu_PlanePointDistance`: unsigned distance; refuses a normal longer than unit.
pub fn plane_point_distance(
    plane_point: Vector3,
    normal: Vector3,
    test: Vector3,
) -> Result<f64, Fault> {
    if normal.length_squared() - 1.0 > EFF_ZERO {
        return Err(Fault::Degenerate);
    }
    Ok((test - plane_point).dot(normal).abs())
}

/// `smgu_PointProjectToPlane`, for a normal of any nonzero length.
pub fn point_project_to_plane(point: Vector3, plane_point: Vector3, normal: Vector3) -> Vector3 {
    let k = ((point.x - plane_point.x) * normal.x
        + (point.y - plane_point.y) * normal.y
        + (point.z - plane_point.z) * normal.z)
        / normal.dot(normal);
    point - k * normal
}

/// `smgu_VectorProjectToPlane`, for a unit normal.
pub fn vector_project_to_plane(vector: Vector3, normal: Vector3) -> Vector3 {
    vector - vector.dot(normal) * normal
}

/// `smgu_VecLinCombTwoVectors`: `(u, v)` with `u * basis1 + v * basis2`
/// nearest `vector`. Fails when the basis is singular at binary64 precision.
pub fn vec_lin_comb_two_vectors(
    basis1: Vector3,
    basis2: Vector3,
    vector: Vector3,
) -> Result<Vector2, Fault> {
    let a00 = basis1.dot(basis1);
    let a01 = basis1.dot(basis2);
    let a10 = a01;
    let a11 = basis2.dot(basis2);
    let det = a00 * a11 - a10 * a01;
    let b0 = basis1.dot(vector);
    let b1 = basis2.dot(vector);
    let numerator0 = b0 * a11 - b1 * a01;
    let numerator1 = a00 * b1 - a10 * b0;
    if det == 0.0
        || det.abs() < numerator0.abs() * EFF_ZERO
        || det.abs() < numerator1.abs() * EFF_ZERO
    {
        return Err(Fault::Degenerate);
    }
    Ok(Vector2::new(numerator0 / det, numerator1 / det))
}

/// `smgu_Determinant3Vectors`: the determinant with rows `v1`, `v2`, `v3`.
pub fn determinant_3_vectors(v1: Vector3, v2: Vector3, v3: Vector3) -> f64 {
    v1.x * v2.y * v3.z - v1.x * v2.z * v3.y + v1.y * v2.z * v3.x - v1.y * v2.x * v3.z
        + v1.z * v2.x * v3.y
        - v1.z * v2.y * v3.x
}

/// `smgu_FrenetFromDerivatives`. From position and the first `count`
/// derivatives (`count` at most 3), writes position, first derivative,
/// curvature vector and torsion vector into `frame`. On a fault the frame
/// keeps whatever was written before it, as in SMLib.
pub fn frenet_from_derivatives(
    derivatives: &[Vector3; 4],
    count: usize,
    frame: &mut [Vector3; 4],
) -> Result<(), Fault> {
    frame[0] = derivatives[0];
    if count == 0 {
        return Ok(());
    }
    let tangent = derivatives[1];
    frame[1] = tangent;
    if count == 1 {
        return Ok(());
    }
    frame[2] = Vector3::new(0.0, 0.0, 0.0);
    if count == 3 {
        frame[3] = Vector3::new(0.0, 0.0, 0.0);
    }
    let tangent_squared = tangent.length_squared();
    if tangent_squared < EFF_ZERO_SQ {
        return Err(Fault::Degenerate);
    }
    let second = derivatives[2];
    frame[2] = tangent.cross(second).cross(tangent) / (tangent_squared * tangent_squared);
    if count == 2 {
        return Ok(());
    }
    if is_zero(second.length()) {
        return Err(Fault::Degenerate);
    }
    let cross_length = tangent.cross(second).length();
    let denominator = cross_length * cross_length;
    if is_zero(denominator) {
        return Err(Fault::Degenerate);
    }
    let torsion = determinant_3_vectors(tangent, second, derivatives[3]) / denominator;
    let mut binormal = frame[1].cross(frame[2]);
    let length = binormal.length();
    if length < EFF_ZERO {
        return Err(Fault::Degenerate);
    }
    binormal *= torsion / length;
    frame[3] = binormal;
    Ok(())
}

/// SMLib's continuity classes (`SmContinuityType`), with its numbering.
///
/// The variants are declared in SMLib's numeric order, so comparisons
/// between classes are SMLib's integer comparisons.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum Continuity {
    /// Not yet classified (`SM_CT_UNDEFINED`).
    Undefined = 0,
    /// Positions differ.
    Discontinuous = 1,
    /// Continuous in position.
    C0 = 2,
    /// Tangent directions agree.
    G1 = 3,
    /// G1 and curvature vectors agree.
    G1G2 = 5,
    /// G1, G2 and third-order agreement.
    G1G2G3 = 6,
    /// First derivatives agree.
    C1 = 7,
    /// C1 and curvature vectors agree.
    C1G2 = 8,
    /// C1, G2 and third-order agreement.
    C1G2G3 = 9,
    /// C1 and second derivatives agree.
    C1C2 = 10,
    /// C1, C2 and third-order agreement.
    C1C2G3 = 11,
    /// Third derivatives agree as well.
    C1C2C3 = 12,
    /// Unlimited continuity (`SM_CT_CINFINITY`).
    CInfinity = 13,
}

/// One end of a continuity comparison: position and derivatives.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurveJet {
    /// Position.
    pub position: Vector3,
    /// First derivative.
    pub first: Vector3,
    /// Second derivative.
    pub second: Vector3,
    /// Third derivative, when known.
    pub third: Option<Vector3>,
}

/// `smgu_EvaluateCurveContinuity`: the continuity class where two curve
/// evaluations meet. Third-order classes need both third derivatives.
///
/// Two SMLib behaviours are kept: the curvature magnitude is the length of
/// the raw second derivative, and the G2 normal-angle test compares the
/// first normal with itself, so it never rejects.
pub fn evaluate_curve_continuity(a: &CurveJet, b: &CurveJet, angle_degrees: f64) -> Continuity {
    let parameter_zero = ZONE_TOL_3D / 10.0;
    let first_max = max(a.first.max_dimension(), b.first.max_dimension());
    let second_max = max(a.second.max_dimension(), b.second.max_dimension());
    let third_max = max(
        a.third.map_or(1.0, Vector3::max_dimension),
        b.third.map_or(1.0, Vector3::max_dimension),
    );

    let scaled_zero =
        EFF_ZERO * (1.0 + max(a.position.max_dimension(), b.position.max_dimension()));
    if a.position.distance(b.position) > scaled_zero {
        return Continuity::Discontinuous;
    }
    let length_a = a.first.length();
    if length_a < EFF_ZERO * (1.0 + a.first.max_dimension()) {
        return Continuity::C0;
    }
    let length_b = b.first.length();
    if length_b < EFF_ZERO * (1.0 + b.first.max_dimension()) {
        return Continuity::C0;
    }
    let mut unit_a = a.first;
    unit_a /= length_a;
    let mut unit_b = b.first;
    unit_b /= length_b;
    if degrees(unit_a.angle_between(unit_b).unwrap_or(0.0)) > angle_degrees {
        return Continuity::C0;
    }
    let mut class = Continuity::G1;
    if a.first.distance(b.first) < parameter_zero * (1.0 + first_max) {
        class = Continuity::C1;
    }

    // SMLib computes both Frenet frames here and never reads them.
    let curvature_a = a.second.length();
    let curvature_b = b.second.length();
    let mut angle = 0.0;
    if !is_zero(curvature_a) && !is_zero(curvature_b) {
        let mut normal = a.second;
        normal /= curvature_a;
        angle = normal.angle_between(normal).unwrap_or(0.0);
    }
    // Negated as in SMLib, so a NaN curvature gap stops at the first-order class.
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    let curvature_differs = !((curvature_a - curvature_b).abs() <= EFF_ZERO_SQRT / 10.0);
    if curvature_differs || degrees(angle) > angle_degrees {
        return class;
    }
    class = if class == Continuity::G1 {
        Continuity::G1G2
    } else {
        Continuity::C1G2
    };
    if class == Continuity::C1G2
        && a.second.distance(b.second) < parameter_zero * (1.0 + second_max)
    {
        class = Continuity::C1C2;
    }

    let (Some(third_a), Some(third_b)) = (a.third, b.third) else {
        return class;
    };
    let origin = Vector3::new(0.0, 0.0, 0.0);
    let mut tolerance = ZONE_TOL_3D / 10.0 * (1.0 + first_max);
    let (gap1, a1) = line_point_distance(origin, b.first, a.first);
    if tolerance < gap1 {
        tolerance = gap1;
    }
    let (gap2, a2) = line_point_distance(origin, b.first, a.second - a1 * a1 * b.second);
    if tolerance < gap2 {
        tolerance = gap2;
    }
    let largest = if first_max > second_max {
        max(first_max, third_max)
    } else {
        max(second_max, third_max)
    };
    tolerance = max(tolerance, ZONE_TOL_3D / 10.0 * (1.0 + largest));
    let target = third_a - a1 * a1 * a1 * third_b - 3.0 * a1 * a2 * b.second;
    let (gap3, _) = line_point_distance(origin, b.first, target);
    if gap3 > tolerance {
        return class;
    }
    class = match class {
        Continuity::G1G2 => Continuity::G1G2G3,
        Continuity::C1G2 => Continuity::C1G2G3,
        Continuity::C1C2 => Continuity::C1C2G3,
        _ => Continuity::G1G2G3,
    };
    if class == Continuity::C1C2G3 && third_a.distance(third_b) < parameter_zero * (1.0 + third_max)
    {
        class = Continuity::C1C2C3;
    }
    class
}

/// `smgu_TriangleDecomposition`: for a triangle with `base` and the other
/// two sides `left` and `right`, the distance along the base from its left
/// end to the foot of the apex, and the apex height.
pub fn triangle_decomposition(base: f64, left: f64, right: f64) -> Result<(f64, f64), Fault> {
    if is_zero(base) {
        return Err(Fault::InvalidInput);
    }
    if left + right < base + (EFF_ZERO * (1.0 + base)) {
        return Err(Fault::InvalidInput);
    }
    let q = (left * left - right * right - base * base) / (-2.0 * base);
    if q.abs() > right {
        return Err(Fault::InvalidInput);
    }
    Ok((base - q, sqrt(right * right - q * q)))
}

/// `smgu_CircleCoPlanarCircleIntersect`: the meetings of two circles in one
/// plane: none, one touch, two crossings, or (count 3) coincident circles.
pub fn circle_coplanar_circle_intersect(
    normal: Vector3,
    center1: Vector3,
    radius1: f64,
    center2: Vector3,
    radius2: f64,
    tolerance: f64,
) -> Result<(usize, [Vector3; 2]), Fault> {
    let mut points = [Vector3::default(); 2];
    let projected = point_project_to_plane(center2, center1, normal);
    let distance = center1.distance(projected);
    if distance < tolerance {
        let count = if (radius1 - radius2).abs() < tolerance {
            3
        } else {
            0
        };
        return Ok((count, points));
    }
    if distance > radius1 + radius2 + tolerance
        || distance + radius1 < radius2 - tolerance
        || distance + radius2 < radius1 - tolerance
    {
        return Ok((0, points));
    }
    let zero = EFF_ZERO * (1.0 + radius1 + radius2);
    if distance > radius1 + radius2 - zero
        || radius1 > distance + radius2 - zero
        || radius2 > distance + radius1 - zero
    {
        let mut vector = projected - center1;
        vector /= distance;
        let scaled_zero = EFF_ZERO * (1.0 + center1.max_dimension());
        let (point1, point2) = if distance < radius1 + radius2 - scaled_zero {
            if radius2 < radius1 {
                (center1 + radius1 * vector, center2 + radius2 * vector)
            } else {
                (center1 - radius1 * vector, center2 - radius2 * vector)
            }
        } else {
            (center1 + radius1 * vector, center2 - radius2 * vector)
        };
        points[0] = 0.5 * (point1 + point2);
        return Ok((1, points));
    }
    let (along, height) = triangle_decomposition(distance, radius1, radius2)?;
    let mut center_vector = projected - center1;
    center_vector.unitize()?;
    let mut perpendicular = center_vector.cross(normal);
    perpendicular.unitize()?;
    points[0] = center1 + along * center_vector + height * perpendicular;
    points[1] = center1 + along * center_vector - height * perpendicular;
    Ok((2, points))
}

/// `SmPolynomial::SolveQuadraticEqn`: the real roots of
/// `c[0] + c[1] x + c[2] x^2` with SMLib's degenerate cases.
pub fn solve_quadratic(c: [f64; 3], zero: f64) -> (usize, [f64; 2]) {
    let mut roots = [0.0; 2];
    if c[2].abs() < zero {
        if c[1].abs() < zero {
            return (0, roots);
        }
        roots[0] = -c[0] / c[1];
        return (1, roots);
    }
    if c[0].abs() < zero {
        roots[0] = 0.0;
        if c[1].abs() < zero {
            return (1, roots);
        }
        roots[1] = -c[1] / c[2];
        return (2, roots);
    }
    let p = c[1] / (2.0 * c[2]);
    let q = c[0] / c[2];
    let discriminant = p * p - q;
    if discriminant.abs() < zero || (discriminant < 0.0 && discriminant.abs() < EFF_ZERO_SQRT) {
        roots[0] = -p;
        return (1, roots);
    }
    if discriminant < 0.0 {
        return (0, roots);
    }
    let root = sqrt(discriminant);
    roots[0] = root - p;
    roots[1] = -root - p;
    (2, roots)
}
