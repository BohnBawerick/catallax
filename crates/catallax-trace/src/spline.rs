// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0
//! Owned spline candidates adapted from SMLib's `SmBSplineCurve.cpp`,
//! `SmBSplineSurface.cpp`, `SmNurbsCrv.cpp`, `SmNurbsSrf.cpp`,
//! `NL_Math.cpp` and `NL_SrfTool.cpp` at
//! c9c5979f04cea8c38213004d4f04df0dcfb40d05.
//!
//! The storage replaces SMLib's contexts, inheritance and raw allocations.
//! Control-point and span arithmetic also adapts rust-brep's `curve/sm_nurbs.rs`,
//! `surface/sm_bspline_net.rs`, `surface/partial_net_decompose.rs` and
//! `basis/basic.rs`. These values carry no proof.

// Indexed recurrences and non-assignment division preserve SMLib arithmetic.
#![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]

use crate::{
    Fault,
    extent::{Extent1, Extent2, Extent3},
    geometry::triangle_area_squared,
    pseudo_box::PseudoBox,
    scalar::EFF_ZERO_SQ,
    vector::{Vector2, Vector3},
};

/// SMLib's marker for a non-rational control point.
pub const NO_WEIGHT: f64 = -f64::MAX;
/// SMLib's marker for a two-dimensional control point.
pub const NO_Z: f64 = f64::MAX;

/// Coordinates requested by SMLib's control-point accessors.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PointForm {
    /// Euclidean coordinates with weight one.
    NonRational,
    /// Weighted coordinates and their weight.
    Homogeneous,
    /// Euclidean coordinates and their original weight.
    Euclidean,
}

/// Homogeneous control point, including SMLib's missing-coordinate markers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ControlPoint {
    /// Stored coordinates, before perspective division.
    pub coordinates: [f64; 4],
}

impl ControlPoint {
    /// Converts stored coordinates by SMLib's `TO_EUCLID` operation order.
    pub fn euclidean(self) -> Vector3 {
        let [x, y, z, w] = self.coordinates;
        if w != NO_WEIGHT {
            Vector3::new(x / w, y / w, if z != NO_Z { z / w } else { 0.0 })
        } else {
            Vector3::new(x, y, if z != NO_Z { z } else { 0.0 })
        }
    }

    /// Reads a control point in one of SMLib's three output forms.
    pub fn get(self, form: PointForm) -> (Vector3, f64) {
        let [x, y, z, w] = self.coordinates;
        let weight = if w == NO_WEIGHT { 1.0 } else { w };
        let point = if form == PointForm::Homogeneous {
            Vector3::new(x, y, if z == NO_Z { 0.0 } else { z })
        } else {
            self.euclidean()
        };
        (
            point,
            if form == PointForm::NonRational {
                1.0
            } else {
                weight
            },
        )
    }
}

fn validate(degree: usize, knots: &[f64], count: usize) -> Result<(), Fault> {
    if !(1..32).contains(&degree)
        || count <= degree
        || count.checked_add(degree + 1) != Some(knots.len())
        || knots.iter().any(|v| !v.is_finite())
        || knots.windows(2).any(|w| w[0] > w[1])
        || knots[degree] >= knots[count]
    {
        return Err(Fault::InvalidInput);
    }
    Ok(())
}

fn validate_points(points: &[ControlPoint]) -> Result<(), Fault> {
    let rational = points[0].coordinates[3] != NO_WEIGHT;
    let dim2 = points[0].coordinates[2] == NO_Z;
    if points.iter().any(|p| {
        p.coordinates.iter().any(|x| !x.is_finite())
            || (p.coordinates[3] != NO_WEIGHT) != rational
            || (p.coordinates[2] == NO_Z) != dim2
            || (rational && p.coordinates[3] == 0.0)
    }) {
        return Err(Fault::InvalidInput);
    }
    Ok(())
}

/// SMLib `N_BasisFindSpan`, with the active interval and repeated-end handling.
pub(crate) fn find_span(knots: &[f64], degree: usize, u: f64, from_left: bool) -> usize {
    let a = &knots[degree..knots.len() - degree];
    let m = a.len() - 1;
    if u >= a[m] {
        let mut i = m - 1;
        while i > 0 && a[i + 1] < a[i] + 1e-12 {
            i -= 1;
        }
        return i + degree;
    }
    if u <= a[0] {
        let mut i = 0;
        while i < m - 1 && a[i + 1] < a[i] + 1e-12 {
            i += 1;
        }
        return i + degree;
    }
    let (mut low, mut high) = (0, m);
    let mut mid = (low + high) / 2;
    loop {
        let outside = if from_left {
            u <= a[mid] || u > a[mid + 1]
        } else {
            u < a[mid] || u >= a[mid + 1]
        };
        if !outside {
            break;
        }
        if if from_left { u > a[mid] } else { u >= a[mid] } {
            low = mid;
        } else {
            high = mid;
        }
        let next = (low + high) / 2;
        if next == mid {
            break;
        }
        mid = next;
    }
    mid + degree
}

fn greville(knots: &[f64], degree: usize, index: usize) -> f64 {
    if index == 0 {
        return knots[0];
    }
    if index >= knots.len() - degree - 2 {
        return knots[knots.len() - 1];
    }
    let mut sum = 0.0;
    for j in 1..=degree {
        sum += knots[index + j];
    }
    sum / degree as f64
}

/// STEP knots and multiplicities, using SMLib's absolute knot equality rule.
pub fn unique_knots(knots: &[f64], interval: Option<Extent1>) -> Vec<(f64, usize)> {
    let mut result: Vec<(f64, usize)> = Vec::new();
    for &value in knots {
        if interval.is_some_and(|iv| !iv.contains_value(value, 0.0)) {
            continue;
        }
        if let Some(last) = result
            .last_mut()
            .filter(|last| (last.0 - value).abs() <= 1e-12)
        {
            last.1 += 1;
        } else {
            result.push((value, 1));
        }
    }
    result
}

fn excess_multiplicity(knots: &[f64], degree: usize) -> bool {
    let unique = unique_knots(knots, None);
    unique
        .iter()
        .enumerate()
        .any(|(i, &(_, n))| n > degree + usize::from(i == 0 || i + 1 == unique.len()))
}

/// An owned, unproved B-spline curve. Cloning replaces SMLib's copy context.
#[derive(Clone, Debug, PartialEq)]
pub struct BSplineCurve {
    degree: usize,
    knots: Vec<f64>,
    points: Vec<ControlPoint>,
}

impl BSplineCurve {
    /// `sm_CrvSplit` with fresh owned outputs. Splits at an interior parameter;
    /// a split at an end or through an excessive knot multiplicity is refused.
    pub fn split(&self, parameter: f64) -> Result<[Self; 2], Fault> {
        let iv = self.natural_interval();
        if !parameter.is_finite() || parameter <= iv.min || parameter >= iv.max {
            return Err(Fault::InvalidInput);
        }
        let p = self.degree;
        let span = find_span(&self.knots, p, parameter, false);
        let mult = self.knots[..=span]
            .iter()
            .rev()
            .take_while(|&&u| u == parameter)
            .count();
        if mult > p {
            return Err(Fault::Degenerate);
        }
        let mut left = vec![self.points[0]; span - mult + 1];
        let mut right = vec![self.points[0]; self.points.len() + p - span];
        let mut sw = self.points[span - p..=span - mult].to_vec();
        left[..=span - p].copy_from_slice(&self.points[..=span - p]);
        right[p - mult..].copy_from_slice(&self.points[span - mult..]);
        for i in 1..=p - mult {
            let k = span - p + i;
            for j in 0..=p - i - mult {
                let alpha = (parameter - self.knots[k + j])
                    / (self.knots[span + j + 1] - self.knots[k + j]);
                sw[j] = combine(alpha, sw[j + 1], 1.0 - alpha, sw[j]);
            }
            left[k] = sw[0];
            right[p - i - mult] = sw[p - i - mult];
        }
        let mut lu = self.knots[..=span - mult].to_vec();
        lu.extend(std::iter::repeat_n(parameter, p + 1));
        let mut ru = vec![parameter; p + 1];
        ru.extend_from_slice(&self.knots[span + 1..]);
        Ok([Self::new(p, lu, left)?, Self::new(p, ru, right)?])
    }
    /// Whether non-rational Bezier pieces meet with full interior multiplicity.
    pub fn is_piecewise_c0(&self) -> bool {
        let unique = unique_knots(&self.knots, None);
        !self.is_rational()
            && unique
                .iter()
                .enumerate()
                .all(|(i, &(_, n))| n == self.degree + usize::from(i == 0 || i + 1 == unique.len()))
    }

    /// Repeated end controls of a non-degenerate curve, using SMLib's default
    /// neighbourhood when `tolerance` is zero or negative.
    pub fn has_repeated_end_control_points(&self, tolerance: f64) -> bool {
        if self.points.len() < 3 || self.is_degenerate(1e-12) {
            return false;
        }
        let tol = if tolerance > 0.0 { tolerance } else { 1e-5 };
        let n = self.points.len() - 1;
        self.points[0]
            .euclidean()
            .distance(self.points[1].euclidean())
            <= tol.abs()
            || self.points[n - 1]
                .euclidean()
                .distance(self.points[n].euclidean())
                <= tol.abs()
    }

    /// SMLib's constant-speed line test. It samples a candidate and proves nothing.
    pub fn is_line(&self, samples: usize, tolerance: f64) -> Result<bool, Fault> {
        let [start, end] = self.ends().map(|v| v.0);
        let base = end - start;
        if base.length_squared() < tolerance * tolerance {
            return Ok(false);
        }
        if !self.is_rational() {
            if self.points.len() == 2 {
                return Ok(true);
            }
            if self.degree == self.points.len() - 1 {
                let standard = base / (self.points.len() - 1) as f64;
                if self.points.windows(2).all(|w| {
                    ((w[1].euclidean() - w[0].euclidean()) - standard).length_squared()
                        <= tolerance * tolerance
                }) {
                    return Ok(true);
                }
            }
        }
        let samples = samples.max(3);
        let spans = (samples - 1) as f64;
        let standard = base / spans;
        let mut previous = start;
        for i in 1..samples {
            let parameter = if i == samples - 1 {
                self.natural_interval().max
            } else {
                self.natural_interval().evaluate(i as f64 / spans)
            };
            let values = self.evaluate(parameter, 1, false, false)?;
            if ((values[0] - previous) - standard).length_squared() > tolerance * tolerance
                || !values[1].is_parallel_to(base, 1.0).0
            {
                return Ok(false);
            }
            previous = values[0];
        }
        Ok(true)
    }

    /// Axis-aligned `CalculateBoundingBox`, including its one-step refinement
    /// for quadratic and cubic Bezier curves.
    pub fn bounding_box(&self, interval: Extent1) -> Result<Extent3, Fault> {
        self.bounding_boxes(interval).map(|(bbox, _)| bbox)
    }

    /// `CalculateBoundingBox` with both its axis-aligned box and its pseudo
    /// box. The pseudo-box basis comes from the widest triangle among up to
    /// seven control points sampled from the start of the whole polygon, as
    /// SMLib samples them, whatever the interval.
    pub fn bounding_boxes(&self, interval: Extent1) -> Result<(Extent3, PseudoBox), Fault> {
        if !interval.min.is_finite() || !interval.max.is_finite() {
            return Err(Fault::InvalidInput);
        }
        let natural = self.natural_interval();
        let interval = Extent1::new(
            natural.clamp_value(interval.min),
            natural.clamp_value(interval.max),
        );
        let (start, count) = if (natural.min - interval.min).abs() <= 1e-12
            && (natural.max - interval.max).abs() <= 1e-12
        {
            (0, self.points.len())
        } else {
            let start = self.find_span(interval.min, true)? - self.degree;
            let end = self.find_span(interval.max, false)?;
            (start, end - start + 1)
        };
        let pseudo = self.pseudo_box_basis(count);
        let mut pseudo = pseudo;
        let mut points = self.points[start..start + count].to_vec();
        if self.degree == self.points.len() - 1 && self.degree == count - 1 {
            let tolerance = 1e-12 * (1.0 + points[0].euclidean().max_dimension());
            if self.is_line(4, tolerance)? {
                let start = self.evaluate_point(interval.min)?;
                let end = self.evaluate_point(interval.max)?;
                let mut bbox = Extent3::singleton(start);
                bbox.add_point(end);
                pseudo.add_point(start);
                pseudo.add_point(end);
                return Ok((bbox, pseudo));
            }
            if self.degree == 2 || self.degree == 3 {
                let mut row = points.clone();
                let mut left = vec![row[0]];
                let mut right = vec![row[row.len() - 1]];
                while row.len() > 1 {
                    row = row.windows(2).map(|p| average(p[0], p[1])).collect();
                    left.push(row[0]);
                    right.push(row[row.len() - 1]);
                }
                // The joining point is an average of its neighbours.
                left.pop();
                right.pop();
                right.reverse();
                left.extend(right);
                points = left;
            }
        }
        let mut bbox = Extent3::singleton(points[0].euclidean());
        for point in points {
            bbox.add_point(point.euclidean());
            pseudo.add_point(point.euclidean());
        }
        Ok((bbox, pseudo))
    }

    /// The pseudo-box basis part of `CalculateBoundingBox` for a sub-polygon
    /// of `count` control points: empty intervals on an orthonormal basis
    /// spanned by the widest sampled triangle.
    fn pseudo_box_basis(&self, count: usize) -> PseudoBox {
        let mut pseudo = PseudoBox::new();
        let samples = count.min(7);
        let step = count / samples;
        let (mut i_max, mut j_max, mut k_max) = match count {
            0 | 1 => (0, 0, 0),
            2 => (0, 0, 1),
            3 => (0, 1, 2),
            _ => (1, 2, 3),
        };
        if count > 3 {
            let mut max_cross = 0.0;
            for i in 0..samples - 2 {
                let p0 = self.points[i * step].euclidean();
                for j in i + 1..samples - 1 {
                    let p1 = self.points[j * step].euclidean();
                    for k in j + 1..samples {
                        let p2 = self.points[k * step].euclidean();
                        let cross = triangle_area_squared(p0, p1, p2);
                        if cross > max_cross {
                            max_cross = cross;
                            (i_max, j_max, k_max) = (i, j, k);
                        }
                    }
                }
            }
        }
        let p0 = self.points[i_max * step].euclidean();
        let mut v0 = self.points[j_max * step].euclidean() - p0;
        let mut v1 = self.points[k_max * step].euclidean() - p0;
        if v1.length_squared() > v0.length_squared() {
            std::mem::swap(&mut v0, &mut v1);
        }
        if v0.length_squared() >= EFF_ZERO_SQ
            && let Ok([x, y, z]) = v0.make_unit_ortho_vectors(Some(v1))
        {
            let _ = pseudo.set_basis(x, y, z);
        }
        pseudo
    }
    /// Copies a raw NURBS curve into checked owned storage.
    pub fn new(degree: usize, knots: Vec<f64>, points: Vec<ControlPoint>) -> Result<Self, Fault> {
        validate(degree, &knots, points.len())?;
        validate_points(&points)?;
        Ok(Self {
            degree,
            knots,
            points,
        })
    }
    /// Polynomial degree.
    pub fn degree(&self) -> usize {
        self.degree
    }
    /// Natural knot vector, including repeated knots.
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }
    /// Stored control polygon.
    pub fn control_points(&self) -> &[ControlPoint] {
        &self.points
    }
    /// Natural interval, matching the first and last stored knots in SMLib.
    pub fn natural_interval(&self) -> Extent1 {
        Extent1::new(self.knots[0], self.knots[self.knots.len() - 1])
    }
    /// Whether the polygon carries weights, including a rational polygon with all weights one.
    pub fn is_rational(&self) -> bool {
        self.points[0].coordinates[3] != NO_WEIGHT
    }
    /// Two or three Euclidean dimensions.
    pub fn dimension(&self) -> usize {
        if self.points[0].coordinates[2] == NO_Z {
            2
        } else {
            3
        }
    }
    /// End positions and weights, read directly from the polygon.
    pub fn ends(&self) -> [(Vector3, f64); 2] {
        [
            self.points[0].get(PointForm::Euclidean),
            self.points[self.points.len() - 1].get(PointForm::Euclidean),
        ]
    }
    /// Greville parameter, clamped to the extreme knot for an extreme index.
    pub fn greville(&self, index: usize) -> f64 {
        greville(&self.knots, self.degree, index)
    }
    /// Span selected from either side of a repeated knot.
    pub fn find_span(&self, u: f64, from_left: bool) -> Result<usize, Fault> {
        if !u.is_finite() {
            return Err(Fault::InvalidInput);
        }
        Ok(find_span(&self.knots, self.degree, u, !from_left))
    }
    /// Control point in the requested form.
    pub fn control_point(&self, index: usize, form: PointForm) -> Result<(Vector3, f64), Fault> {
        self.points
            .get(index)
            .map(|p| p.get(form))
            .ok_or(Fault::InvalidInput)
    }
    /// Control-polygon degeneracy test for the natural interval.
    pub fn is_degenerate(&self, tolerance: f64) -> bool {
        let mut bbox = Extent3::default();
        for p in &self.points {
            bbox.add_point(p.euclidean());
        }
        bbox.is_point_sized(tolerance.max(1e-12 * (1.0 + bbox.max_dimension())))
    }
    /// Excessive interior or end knot multiplicity. Degenerate curves are exempt in SMLib.
    pub fn has_excess_knot_multiplicity(&self) -> bool {
        !self.is_degenerate(1e-12)
            && self.points.len() >= 3
            && excess_multiplicity(&self.knots, self.degree)
    }
    /// `EvaluatePoint`, including the direct endpoint shortcut and parameter clamp.
    pub fn evaluate_point(&self, u: f64) -> Result<Vector3, Fault> {
        if !u.is_finite() {
            return Err(Fault::InvalidInput);
        }
        let u = self.natural_interval().clamp_value(u);
        if u == self.knots[0] {
            return finite(self.points[0].euclidean());
        }
        if u == self.knots[self.knots.len() - 1] {
            return finite(self.points[self.points.len() - 1].euclidean());
        }
        // sm_BasisEval stores the active endpoint basis but returns NL_YES.
        // EvaluatePoint's natural-end shortcut handles the clamped case first.
        if u == self.knots[self.degree] || u == self.knots[self.knots.len() - self.degree - 1] {
            return Err(Fault::Degenerate);
        }
        let (span, basis) = basis_values(&self.knots, self.degree, u)?;
        let mut point = [0.0; 4];
        for (i, &b) in basis.iter().enumerate() {
            let cp = self.points[span - self.degree + i].coordinates;
            for k in 0..4 {
                if cp[k] == NO_Z || cp[k] == NO_WEIGHT {
                    point[k] = cp[k];
                } else if b == 1.0 {
                    point[k] += cp[k];
                } else if b != 0.0 {
                    point[k] += b * cp[k];
                }
            }
        }
        finite(ControlPoint { coordinates: point }.euclidean())
    }
    /// `Evaluate`, through fourth derivatives. `nonzero_tangent` reproduces
    /// SMLib's second-derivative replacement of a vanishing first derivative.
    pub fn evaluate(
        &self,
        u: f64,
        derivatives: usize,
        from_left: bool,
        nonzero_tangent: bool,
    ) -> Result<Vec<Vector3>, Fault> {
        if derivatives == 0 {
            return Ok(vec![self.evaluate_point(u)?]);
        }
        if derivatives > 4 || !u.is_finite() {
            return Err(Fault::InvalidInput);
        }
        let u = self.natural_interval().clamp_value(u);
        let (span, mut basis) =
            basis_derivatives(&self.knots, self.degree, u, !from_left, derivatives)?;
        if self.is_rational() {
            let mut denominator = vec![0.0; derivatives + 1];
            for k in 0..=derivatives {
                for j in 0..=self.degree {
                    denominator[k] +=
                        self.points[span - self.degree + j].coordinates[3] * basis[k][j];
                }
            }
            if denominator[0] == 0.0 {
                return Err(Fault::Degenerate);
            }
            let ordinary = basis.clone();
            for i in 0..=self.degree {
                let w = self.points[span - self.degree + i].coordinates[3];
                for k in 0..=derivatives {
                    let mut value = w * ordinary[k][i];
                    for j in 1..=k {
                        value -= binomial(k, j) * denominator[j] * basis[k - j][i];
                    }
                    basis[k][i] = value / denominator[0];
                }
            }
        }
        let mut result = vec![Vector3::new(0.0, 0.0, 0.0); derivatives + 1];
        for k in 0..=derivatives {
            for i in 0..=self.degree {
                let cp = self.points[span - self.degree + i];
                // sm_CrvGetEPolygon calls N_CPtToPtEuclid, which multiplies
                // by a reciprocal rather than using TO_EUCLID's divisions.
                let [x, y, z, w] = cp.coordinates;
                let p = if w == NO_WEIGHT {
                    cp.euclidean()
                } else {
                    let reciprocal = 1.0 / w;
                    Vector3::new(
                        x * reciprocal,
                        y * reciprocal,
                        if z == NO_Z { 0.0 } else { z * reciprocal },
                    )
                };
                result[k] = result[k] + basis[k][i] * p;
            }
            finite(result[k])?;
        }
        if nonzero_tangent {
            let new_length = 1.1 * (1e-12 * (1.0 + result[0].max_dimension()));
            if result[1].length_squared() < new_length * new_length {
                if derivatives < 2 {
                    // The recursive C++ call omits its fifth argument. The
                    // declaration defaults to FALSE, despite the body's note.
                    result[1] = self.evaluate(u, 2, from_left, false)?[1];
                } else {
                    let second = result[2];
                    let length = second.length();
                    if length > new_length {
                        result[1] = (new_length / length) * second;
                        if u > self.natural_interval().mid() {
                            result[1] *= -1.0;
                        }
                    }
                }
            }
        }
        Ok(result)
    }
}

fn average(a: ControlPoint, b: ControlPoint) -> ControlPoint {
    let mut coordinates = [0.0; 4];
    for (i, value) in coordinates.iter_mut().enumerate() {
        *value = if i == 2 && a.coordinates[i] == NO_Z {
            NO_Z
        } else if i == 3 && a.coordinates[i] == NO_WEIGHT {
            NO_WEIGHT
        } else {
            (a.coordinates[i] + b.coordinates[i]) / 2.0
        };
    }
    ControlPoint { coordinates }
}

fn finite(v: Vector3) -> Result<Vector3, Fault> {
    if v.x.is_finite() && v.y.is_finite() && v.z.is_finite() {
        Ok(v)
    } else {
        Err(Fault::Degenerate)
    }
}

fn binomial(n: usize, k: usize) -> f64 {
    let mut value = 1usize;
    for i in 0..k {
        value = value * (n - i) / (i + 1);
    }
    value as f64
}

/// Owned NURBS surface, with V varying fastest in the control net.
///
/// An `SmPlane` is this type with a [`PlaneFrame`]: SMLib's plane derives
/// from its B-spline surface, keeps the bilinear net of `MakeNurb` for its
/// caches and overrides evaluation and a few queries with the analytic
/// frame.
#[derive(Clone, Debug, PartialEq)]
pub struct BSplineSurface {
    degrees: [usize; 2],
    knots: [Vec<f64>; 2],
    counts: [usize; 2],
    points: Vec<ControlPoint>,
    plane: Option<PlaneFrame>,
}

/// `SM_INFINITE_PARAMETER`: the half-width of an unbounded plane's domain.
pub const INFINITE_PARAMETER: f64 = 1_234_567.0;

/// The analytic data of an `SmPlane`: the point at `(u, v)` is
/// `origin + u * u_scale * x_axis + v * v_scale * y_axis` over `domain`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlaneFrame {
    /// `m_vPosition` origin.
    pub origin: Vector3,
    /// `m_vPosition` unit X axis.
    pub x_axis: Vector3,
    /// `m_vPosition` unit Y axis, perpendicular to the X axis.
    pub y_axis: Vector3,
    /// `m_vUVScale`.
    pub uv_scale: Vector2,
    /// `m_vAnalUVDomain`.
    pub domain: Extent2,
}

impl PlaneFrame {
    /// `SmAxis2Placement::SetCanonical` for the plane's position: unit axes,
    /// made perpendicular from the given ones when they are not.
    pub(crate) fn canonical(
        origin: Vector3,
        x_axis: Vector3,
        y_axis: Vector3,
    ) -> Result<[Vector3; 3], Fault> {
        let mut x = x_axis;
        let mut y = y_axis;
        if (x.dot(x) - 1.0).abs() > crate::scalar::EFF_ZERO / 100.0 {
            x.unitize()?;
        }
        if (y.dot(y) - 1.0).abs() > crate::scalar::EFF_ZERO / 100.0 {
            y.unitize()?;
        }
        if x.dot(y).abs() > crate::scalar::EFF_ZERO {
            let [ux, uy, _] = x_axis.make_unit_ortho_vectors(Some(y_axis))?;
            x = ux;
            y = uy;
        }
        Ok([origin, x, y])
    }

    /// `SmPlane::EvaluatePointFast`: the point without clamping.
    pub fn point_fast(&self, uv: Vector2) -> Vector3 {
        self.origin + uv.x * self.x_axis * self.uv_scale.x + uv.y * self.y_axis * self.uv_scale.y
    }

    /// `SmAxis2Placement::GetZAxis`: the unit normal `x × y` (an
    /// unnormalisable product is returned as computed).
    pub fn z_axis(&self) -> Vector3 {
        let mut z = self.x_axis.cross(self.y_axis);
        let _ = z.unitize();
        z
    }

    /// `SmPlane::ProjectPointToUVDomain`.
    pub fn project(&self, point: Vector3) -> Vector2 {
        let offset = point - self.origin;
        Vector2::new(
            self.x_axis.dot(offset) / self.uv_scale.x,
            self.y_axis.dot(offset) / self.uv_scale.y,
        )
    }
}

impl BSplineSurface {
    /// `sm_SplitSrf` with owned outputs. Each net row uses the same insertion
    /// recurrence as `N_SrfSplit`, without the C++ output-allocation protocol.
    pub fn split(&self, axis: usize, parameter: f64) -> Result<[Self; 2], Fault> {
        if axis > 1 {
            return Err(Fault::InvalidInput);
        }
        let other = 1 - axis;
        let mut rows = Vec::new();
        for j in 0..self.counts[other] {
            let points = (0..self.counts[axis])
                .map(|i| {
                    if axis == 0 {
                        self.points[i * self.counts[1] + j]
                    } else {
                        self.points[j * self.counts[1] + i]
                    }
                })
                .collect();
            rows.push(
                BSplineCurve::new(self.degrees[axis], self.knots[axis].clone(), points)?
                    .split(parameter)?,
            );
        }
        let mut outputs = Vec::new();
        for side in 0..2 {
            let size = rows[0][side].points.len();
            let mut counts = self.counts;
            counts[axis] = size;
            let mut points = vec![self.points[0]; counts[0] * counts[1]];
            for j in 0..self.counts[other] {
                for i in 0..size {
                    let index = if axis == 0 {
                        i * counts[1] + j
                    } else {
                        j * counts[1] + i
                    };
                    points[index] = rows[j][side].points[i];
                }
            }
            let mut knots = self.knots.clone();
            knots[axis] = rows[0][side].knots.clone();
            outputs.push(Self::new(self.degrees, knots, counts, points)?);
        }
        let [left, right] = outputs.try_into().unwrap();
        Ok([left, right])
    }
    /// `CreateIsoParametricCurve` without an optional trim. `axis` is the
    /// constant parameter, zero for U and one for V. Extraction uses the
    /// original weighted net, including SMLib's direct boundary copy.
    pub fn iso_curve(&self, axis: usize, parameter: f64) -> Result<BSplineCurve, Fault> {
        if axis > 1 || !parameter.is_finite() {
            return Err(Fault::InvalidInput);
        }
        let knots = &self.knots[axis];
        let iv = Extent1::new(knots[0], knots[knots.len() - 1]);
        if !iv.contains_value(parameter, 0.0) && !iv.is_value_on_boundary(parameter, 1e-6).0 {
            return Err(Fault::Degenerate);
        }
        let parameter = iv.clamp_value(parameter);
        let degree = self.degrees[axis];
        let other = 1 - axis;
        let at = |i, j| {
            if axis == 0 {
                self.points[i * self.counts[1] + j]
            } else {
                self.points[j * self.counts[1] + i]
            }
        };
        let points = if parameter == knots[degree] || parameter == knots[knots.len() - degree - 1] {
            let index = if parameter == knots[degree] {
                0
            } else {
                self.counts[axis] - 1
            };
            (0..self.counts[other]).map(|j| at(index, j)).collect()
        } else {
            let span = find_span(knots, degree, parameter, false);
            let multiplicity = knots[..=span]
                .iter()
                .rev()
                .take_while(|&&u| u == parameter)
                .count();
            if multiplicity > degree {
                return Err(Fault::Degenerate);
            }
            let mut points = Vec::new();
            for j in 0..self.counts[other] {
                let mut row: Vec<_> = (0..=degree - multiplicity)
                    .map(|i| at(span - degree + i, j))
                    .collect();
                for i in 1..=degree - multiplicity {
                    let a = span - degree + i;
                    for k in 0..=degree - i - multiplicity {
                        let alpha =
                            (parameter - knots[a + k]) / (knots[span + k + 1] - knots[a + k]);
                        row[k] = combine(alpha, row[k + 1], 1.0 - alpha, row[k]);
                    }
                }
                points.push(row[0]);
            }
            points
        };
        BSplineCurve::new(self.degrees[other], self.knots[other].clone(), points)
    }
    /// Copies a control net and two full knot vectors into checked storage.
    pub fn new(
        degrees: [usize; 2],
        knots: [Vec<f64>; 2],
        counts: [usize; 2],
        points: Vec<ControlPoint>,
    ) -> Result<Self, Fault> {
        if counts[0].checked_mul(counts[1]) != Some(points.len()) {
            return Err(Fault::InvalidInput);
        }
        for axis in 0..2 {
            validate(degrees[axis], &knots[axis], counts[axis])?;
        }
        validate_points(&points)?;
        Ok(Self {
            degrees,
            knots,
            counts,
            points,
            plane: None,
        })
    }
    /// `SmPlane::SmPlane` from an origin, axes, UV scale and domain, with
    /// the bilinear net of `SmPlane::MakeNurb` over the domain's corners.
    pub fn plane(
        origin: Vector3,
        x_axis: Vector3,
        y_axis: Vector3,
        uv_scale: Vector2,
        domain: Extent2,
    ) -> Result<Self, Fault> {
        // A zero or non-finite scale collapses the plane and makes its
        // projection divide by zero; SMLib's constructor does not check.
        let finite = |v: Vector3| v.x.is_finite() && v.y.is_finite() && v.z.is_finite();
        if !(finite(origin) && finite(x_axis) && finite(y_axis))
            || !(uv_scale.x.is_finite() && uv_scale.y.is_finite())
            || uv_scale.x == 0.0
            || uv_scale.y == 0.0
            || ![domain.min.x, domain.min.y, domain.max.x, domain.max.y]
                .iter()
                .all(|c| c.is_finite())
        {
            return Err(Fault::InvalidInput);
        }
        let [origin, x_axis, y_axis] = PlaneFrame::canonical(origin, x_axis, y_axis)?;
        let frame = PlaneFrame {
            origin,
            x_axis,
            y_axis,
            uv_scale,
            domain,
        };
        let corner = |a: f64, b: f64| {
            let p = frame.point_fast(domain.evaluate(a, b));
            ControlPoint {
                coordinates: [p.x, p.y, p.z, NO_WEIGHT],
            }
        };
        let points = vec![
            corner(0.0, 0.0),
            corner(0.0, 1.0),
            corner(1.0, 0.0),
            corner(1.0, 1.0),
        ];
        let u = vec![domain.min.x, domain.min.x, domain.max.x, domain.max.x];
        let v = vec![domain.min.y, domain.min.y, domain.max.y, domain.max.y];
        let mut surface = Self::new([1, 1], [u, v], [2, 2], points)?;
        surface.plane = Some(frame);
        Ok(surface)
    }
    /// `SmPlane::CreateCanonical`: a plane through `origin` spanned by the
    /// axes, with unit scale over `±SM_INFINITE_PARAMETER`.
    pub fn canonical_plane(
        origin: Vector3,
        x_axis: Vector3,
        y_axis: Vector3,
    ) -> Result<Self, Fault> {
        Self::plane(
            origin,
            x_axis,
            y_axis,
            Vector2::new(1.0, 1.0),
            Extent2::new(
                Vector2::new(-INFINITE_PARAMETER, -INFINITE_PARAMETER),
                Vector2::new(INFINITE_PARAMETER, INFINITE_PARAMETER),
            ),
        )
    }
    /// The analytic frame of an `SmPlane`, `None` for other surfaces.
    pub fn plane_frame(&self) -> Option<&PlaneFrame> {
        self.plane.as_ref()
    }
    /// U and V degrees.
    pub fn degrees(&self) -> [usize; 2] {
        self.degrees
    }
    /// Numbers of control points in U and V.
    pub fn counts(&self) -> [usize; 2] {
        self.counts
    }
    /// Full U and V knot vectors.
    pub fn knots(&self) -> [&[f64]; 2] {
        [&self.knots[0], &self.knots[1]]
    }
    /// Stored control net, with V varying fastest.
    pub fn control_points(&self) -> &[ControlPoint] {
        &self.points
    }
    /// Parameter rectangle from the extreme knots.
    pub fn natural_domain(&self) -> Extent2 {
        Extent2::new(
            Vector2::new(self.knots[0][0], self.knots[1][0]),
            Vector2::new(
                self.knots[0][self.knots[0].len() - 1],
                self.knots[1][self.knots[1].len() - 1],
            ),
        )
    }
    /// Whether the net stores rational weights.
    pub fn is_rational(&self) -> bool {
        self.points[0].coordinates[3] != NO_WEIGHT
    }
    /// Both control-net dimensions are two, as in SMLib's `IsBilinear`.
    pub fn is_bilinear(&self) -> bool {
        self.counts == [2, 2]
    }
    /// Greville parameter in U or V.
    pub fn greville(&self, axis: usize, index: usize) -> Result<f64, Fault> {
        if axis > 1 {
            return Err(Fault::InvalidInput);
        }
        Ok(greville(&self.knots[axis], self.degrees[axis], index))
    }
    /// Knot spans from either side of a repeated knot in each direction.
    pub fn find_spans(&self, uv: Vector2, from_left: [bool; 2]) -> Result<[usize; 2], Fault> {
        if !uv.x.is_finite() || !uv.y.is_finite() {
            return Err(Fault::InvalidInput);
        }
        Ok([
            find_span(&self.knots[0], self.degrees[0], uv.x, from_left[0]),
            find_span(&self.knots[1], self.degrees[1], uv.y, from_left[1]),
        ])
    }
    /// Control point in the requested form.
    pub fn control_point(
        &self,
        index: [usize; 2],
        form: PointForm,
    ) -> Result<(Vector3, f64), Fault> {
        if index[0] >= self.counts[0] || index[1] >= self.counts[1] {
            return Err(Fault::InvalidInput);
        }
        Ok(self.points[index[0] * self.counts[1] + index[1]].get(form))
    }
    /// Copies a rectangular Euclidean control net (`sm_SurfaceGetENet`).
    pub fn euclidean_net(&self, low: [usize; 2], high: [usize; 2]) -> Result<Vec<Vector3>, Fault> {
        if (0..2).any(|a| low[a] > high[a] || high[a] >= self.counts[a]) {
            return Err(Fault::InvalidInput);
        }
        let mut result = Vec::new();
        for u in low[0]..=high[0] {
            for v in low[1]..=high[1] {
                result.push(self.control_point([u, v], PointForm::Euclidean)?.0);
            }
        }
        Ok(result)
    }
    /// Excessive knot multiplicity in either direction.
    pub fn has_excess_knot_multiplicity(&self) -> bool {
        (0..2).any(|a| excess_multiplicity(&self.knots[a], self.degrees[a]))
    }
    /// Control-net box and the scale-dependent `sm_IsNSrfDegenerate` test.
    pub fn is_degenerate(&self, tolerance: f64) -> bool {
        let mut bbox = Extent3::default();
        for p in &self.points {
            bbox.add_point(p.euclidean());
        }
        bbox.is_point_sized(tolerance.max(1000.0 * 1e-12 * (1.0 + bbox.max_dimension())))
    }
    /// `sm_FindNetExtrema`, preserving stable traversal order when distances tie.
    pub fn net_extrema(&self, test: Vector3, count: usize) -> Result<NetExtrema, Fault> {
        if count > 10 {
            return Err(Fault::InvalidInput);
        }
        finite(test)?;
        let mut result = NetExtrema {
            nearest: Vec::new(),
            farthest: Vec::new(),
        };
        for u in 0..self.counts[0] {
            for v in 0..self.counts[1] {
                let d = (test - self.points[u * self.counts[1] + v].euclidean()).length_squared();
                for (list, nearest) in [(&mut result.nearest, true), (&mut result.farthest, false)]
                {
                    let at = list
                        .iter()
                        .position(|&(old, _)| if nearest { d < old } else { d > old })
                        .unwrap_or(list.len());
                    if at < count {
                        list.insert(at, (d, [u, v]));
                        list.truncate(count);
                    }
                }
            }
        }
        Ok(result)
    }
}

pub(crate) fn combine(alpha: f64, a: ControlPoint, beta: f64, b: ControlPoint) -> ControlPoint {
    let mut coordinates = [0.0; 4];
    for (i, value) in coordinates.iter_mut().enumerate() {
        *value = if i == 2 && a.coordinates[i] == NO_Z {
            NO_Z
        } else if i == 3 && a.coordinates[i] == NO_WEIGHT {
            NO_WEIGHT
        } else {
            alpha * a.coordinates[i] + beta * b.coordinates[i]
        };
    }
    ControlPoint { coordinates }
}

/// Squared distances and control-net indexes from `sm_FindNetExtrema`.
#[derive(Clone, Debug, PartialEq)]
pub struct NetExtrema {
    /// Nearest points in ascending distance order.
    pub nearest: Vec<(f64, [usize; 2])>,
    /// Farthest points in descending distance order.
    pub farthest: Vec<(f64, [usize; 2])>,
}

pub(crate) fn basis_derivatives(
    knots: &[f64],
    degree: usize,
    u: f64,
    from_left: bool,
    derivatives: usize,
) -> Result<(usize, Vec<Vec<f64>>), Fault> {
    // Same-knot tolerance: see discussion by its first use, below.
    let d_knot_tol = 1e-24;
    let uk = knots;
    let p = degree as i64;
    let der = derivatives as i64;
    let mut nd = vec![vec![0.0; degree + 1]; derivatives + 1];

    let mut left = [0.0f64; 32];
    let mut right = [0.0f64; 32];
    let mut ndu = [[0.0f64; 32]; 32];
    let mut a = [[0.0f64; 32]; 32];

    // Find the knot span u is in
    let i = find_span(knots, p as usize, u, from_left) as i64;
    // In case of early exit:
    // Get maximum derivative index and set zero derivatives
    let mder = if p < der { p } else { der };
    for k in (p + 1)..=der {
        for j in 0..=p {
            nd[k as usize][j as usize] = 0.0;
        }
    }

    // Compute the basis functions
    ndu[0][0] = 1.0;
    for j in 1..=p as usize {
        left[j] = u - uk[i as usize + 1 - j];
        right[j] = uk[i as usize + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];

            // Divide by zero? Only a zero sum can be a floating-point problem; keep the divisions
            // inside the parentheses.
            if ndu[j][r].abs() <= d_knot_tol {
                return Err(Fault::Degenerate);
            }

            ndu[r][j] = saved + ndu[r][j - 1] * (right[r + 1] / ndu[j][r]);
            saved = ndu[r][j - 1] * (left[j - r] / ndu[j][r]);
        }
        ndu[j][j] = saved;
    }

    // Load the basis functions
    for j in 0..=p as usize {
        nd[0][j] = ndu[j][p as usize];
    }

    // Compute derivatives
    for r in 0..=p {
        let (mut s1, mut s2) = (0usize, 1usize);
        a[0][0] = 1.0;
        for k in 1..=mder {
            let mut d = 0.0;
            let rk = r - k;
            let pk = p - k;
            let row = (pk + 1) as usize;
            if r >= k {
                if ndu[row][rk as usize].abs() <= d_knot_tol {
                    return Err(Fault::Degenerate);
                }
                a[s2][0] = a[s1][0] / ndu[row][rk as usize];
                d = a[s2][0] * ndu[rk as usize][pk as usize];
            }
            let j1 = if rk >= -1 { 1 } else { -rk };
            let j2 = if (r - 1) <= pk { k - 1 } else { p - r };
            for j in j1..=j2 {
                let col = (rk + j) as usize;
                if ndu[row][col].abs() <= d_knot_tol {
                    return Err(Fault::Degenerate);
                }
                a[s2][j as usize] = (a[s1][j as usize] - a[s1][j as usize - 1]) / ndu[row][col];
                d += a[s2][j as usize] * ndu[col][pk as usize];
            }
            if r <= pk {
                if ndu[row][r as usize].abs() <= d_knot_tol {
                    return Err(Fault::Degenerate);
                }
                a[s2][k as usize] = -a[s1][k as usize - 1] / ndu[row][r as usize];
                d += a[s2][k as usize] * ndu[r as usize][pk as usize];
            }
            nd[k as usize][r as usize] = d;
            std::mem::swap(&mut s1, &mut s2);
        }
    }

    // Multiply through by the correct factors
    let mut r = p;
    for k in 1..=mder {
        for j in 0..=p as usize {
            nd[k as usize][j] *= r as f64;
        }
        r = r.wrapping_mul(p - k);
    }

    // Made it through.
    Ok((i as usize, nd))
}

/// Constants of a rectangular subset of the control net, before any proof.
#[derive(Clone, Debug, PartialEq)]
pub struct PartialNet {
    /// Corners indexed first by U, then by V.
    pub corners: [[Vector3; 2]; 2],
    /// Axis-aligned box accumulated by SMLib's row traversal.
    pub bounding_box: Extent3,
    /// Maximum chord heights in U and V.
    pub chord_heights: [f64; 2],
    /// Accumulated polygon turning angles in U and V, in degrees.
    pub angles: [f64; 2],
}

impl BSplineSurface {
    /// `CalculatePartialNetValues`, including SMLib's side choices at domain ends.
    /// This computes the normal box, chord heights and angles. Pseudo and polar
    /// boxes belong to the later cache units.
    pub fn partial_net(&self, domain: Extent2) -> Result<PartialNet, Fault> {
        let low = self.find_spans(domain.min, [false, false])?;
        let high = self.find_spans(domain.max, [true, true])?;
        self.partial_net_spans(low, high)
    }

    /// `sm_ComputePartialNetConstants` for inclusive knot-span indexes.
    pub fn partial_net_spans(
        &self,
        low: [usize; 2],
        high: [usize; 2],
    ) -> Result<PartialNet, Fault> {
        if (0..2).any(|a| low[a] > high[a] || low[a] < self.degrees[a] || high[a] >= self.counts[a])
        {
            return Err(Fault::InvalidInput);
        }
        let start = [low[0] - self.degrees[0], low[1] - self.degrees[1]];
        let point = |i: usize, j: usize| self.points[i * self.counts[1] + j].euclidean();
        let mut result = PartialNet {
            corners: [
                [point(start[0], start[1]), point(start[0], high[1])],
                [point(high[0], start[1]), point(high[0], high[1])],
            ],
            bounding_box: Extent3::default(),
            chord_heights: [0.0; 2],
            angles: [0.0; 2],
        };
        // SMLib visits V rows before U columns. A collapsed row is skipped,
        // including its box writes, and an angle break also ends its chord scan.
        for axis in [1, 0] {
            let other = 1 - axis;
            let at = |fixed, varying| {
                if axis == 1 {
                    point(fixed, varying)
                } else {
                    point(varying, fixed)
                }
            };
            for fixed in start[other]..=high[other] {
                let line_point = at(fixed, start[axis]);
                let mut line_vec = at(fixed, high[axis]) - line_point;
                let mut try_index = high[axis] - 1;
                while line_vec.length_squared() < 1e-24 {
                    line_vec = at(fixed, try_index) - line_point;
                    if try_index == 0 {
                        break;
                    }
                    try_index -= 1;
                }
                if line_vec.length_squared() < 1e-24 {
                    continue;
                }
                let mut total_angle = 0.0;
                let mut last_vec = Vector3::new(0.0, 0.0, 0.0);
                let mut current = line_point;
                if axis == 1 {
                    result.bounding_box.add_point(current);
                }
                for i in start[axis] + 1..=high[axis] {
                    let next = at(fixed, i);
                    if axis == 1 {
                        result.bounding_box.add_point(next);
                    }
                    let param = crate::geometry::line_closest_point(line_point, line_vec, next);
                    let on_line = line_point + line_vec * param;
                    let chord = (next - on_line).length_squared();
                    if chord > result.chord_heights[axis] {
                        result.chord_heights[axis] = chord;
                    }
                    let mut current_vec = next - current;
                    let distance = current_vec.length();
                    let tolerance = 1e-6 * 100.0 * (1.0 + next.max_dimension());
                    if distance < tolerance {
                        break;
                    }
                    current_vec = current_vec / distance;
                    if last_vec.length_squared() > 1e-24 {
                        let cos = current_vec.dot(last_vec).clamp(-1.0, 1.0);
                        total_angle += if cos >= 1.0 {
                            0.0
                        } else if cos <= -1.0 {
                            std::f64::consts::PI
                        } else {
                            libm::acos(cos)
                        };
                    }
                    last_vec = current_vec;
                    current = next;
                }
                if total_angle > result.angles[axis] {
                    result.angles[axis] = total_angle;
                }
            }
        }
        for axis in 0..2 {
            result.chord_heights[axis] = result.chord_heights[axis].sqrt();
            result.angles[axis] = result.angles[axis] * 180.0 / std::f64::consts::PI;
        }
        Ok(result)
    }
}

impl BSplineSurface {
    /// Repeated end controls, excluding collapsed boundary net curves.
    pub fn has_repeated_end_control_points(&self, tol: f64) -> bool {
        let [n, m] = self.counts.map(|n| n as i64 - 1);
        let pw = &self.points;
        let at = |i: i64, j: i64| pw[(i * (m + 1) + j) as usize];
        let tolerance = if tol > 0.0 { tol } else { 1e-5 };
        // for every edge of the surface polyhedron
        for ii in 0..4 {
            let (j_max, row_length, id, jd, mut i0, mut j0, mut i1, mut j1) = match ii {
                0 => (n, m, 1, 0, 0, 0, 0, 1),
                1 => (n, m, 1, 0, 0, m, 0, m - 1),
                2 => (m, n, 0, 1, 0, 0, 1, 0),
                _ => (m, n, 0, 1, n, 0, n - 1, 0),
            };
            // skip net curves with less than 3 control points
            if row_length < 2 {
                continue;
            }
            let mut jj = 0;
            while jj <= j_max {
                // ignore repeated control points on singular boundaries
                let mut degenerate = false;
                if jj == 0 || jj == j_max {
                    let mut bbox = Extent3::default();
                    let (mut ki, kid, mut kj, kjd) = if ii == 0 || ii == 1 {
                        (jj, 0, 0, 1)
                    } else {
                        (0, 1, jj, 0)
                    };
                    let mut kk = 0;
                    while kk <= row_length {
                        bbox.add_point(at(ki, kj).euclidean());
                        if kk == 1 {
                            degenerate = bbox.is_point_sized(tolerance);
                            if !degenerate {
                                break;
                            }
                        }
                        kk += 1;
                        ki += kid;
                        kj += kjd;
                    }
                    if degenerate {
                        degenerate = bbox.is_point_sized(tolerance);
                    }
                }
                if !degenerate {
                    let dist = at(i0, j0).euclidean().distance(at(i1, j1).euclidean());
                    if dist.abs() <= tolerance.abs() {
                        return true;
                    }
                }
                jj += 1;
                i0 += id;
                i1 += id;
                j0 += jd;
                j1 += jd;
            }
        }
        false
    }
}

// sm_BasisEval differs from sm_BasisDerivs at the exact knot-spacing threshold.
fn basis_values(knots: &[f64], degree: usize, u: f64) -> Result<(usize, Vec<f64>), Fault> {
    let span = find_span(knots, degree, u, false);
    let mut basis = vec![0.0; degree + 1];
    basis[0] = 1.0;
    let mut left = vec![0.0; degree + 1];
    let mut right = left.clone();
    for k in 1..=degree {
        left[k] = u - knots[span + 1 - k];
        right[k] = knots[span + k] - u;
        let mut saved = 0.0;
        for j in 0..k {
            let delta = right[j + 1] + left[k - j];
            if delta < 1e-24 {
                return Err(Fault::Degenerate);
            }
            let temp = basis[j];
            basis[j] = saved + temp * (right[j + 1] / delta);
            saved = temp * (left[k - j] / delta);
        }
        basis[k] = saved;
    }
    Ok((span, basis))
}
