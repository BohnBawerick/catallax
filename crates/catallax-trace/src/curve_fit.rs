// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: NLib stacks and matrix descriptors become owned
// vectors, statuses become `Fault`, and the vectors SMLib rewrites in place
// are returned. Translated directly from SMLib SmBSplineCurveNLib.cpp
// (`CreateInterpolatingCurve`), NL_CrvFit.cpp (`N_FitCrvFirstDeriv`,
// `N_FitCalcCrvParamValues`, `N_FitCalcKnotsDerivs`), NL_Math.cpp
// (`N_RealMatrixLuDecompose`, `N_RealMatrixForBack`) and NL_Geom.cpp
// (`N_DistPtPt`, `N_DistPolygon`). The matching rust-brep unit
// (`curve/crv_fit_tangents_first_derivs.rs`) was read as a cross-check.
//! Interpolating a curve through points with a first derivative at each.
//!
//! The curve is degree 2 or 3 with a knot at every data parameter and halfway
//! (degree 2) or a third of the way (degree 3) between them. It solves a
//! banded system whose storage SMLib addresses by raw offsets: an entry that
//! falls outside its row but inside the block lands in the neighbouring row,
//! and the port reproduces that.
//!
//! The other branches of `CreateInterpolatingCurve` (end derivatives, higher
//! derivatives, none) need `N_FitCrvHighDerivs` and are not carried; they
//! refuse with [`Fault::Degenerate`] and the gate records SMLib's result.
//!
//! Every curve is an unproved candidate.

// Indexed loops keep SMLib's index arithmetic.
#![allow(clippy::needless_range_loop)]

use crate::{
    Fault,
    curve_build::constructed,
    curve_compat::division_is_bad,
    geometry::solve_quadratic,
    scalar::{EFF_ZERO, EFF_ZERO_SQ},
    spline::{BSplineCurve, ControlPoint, NO_WEIGHT, NO_Z, basis_derivatives},
    vector::Vector3,
};

/// `NL_LUDT`, the pivot tolerance of the LU decomposition.
const LU_PIVOT_TOLERANCE: f64 = 1.0e-9;
/// `NL_PTOL`, the parameter-space tolerance.
const PARAMETER_TOLERANCE: f64 = 1.0e-12;

/// `SmCurveParameterizationType`: how the data parameters are spaced.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Parameterization {
    /// `SM_CP_UNIFORM`.
    Uniform,
    /// `SM_CP_CHORDLENGTH`.
    ChordLength,
    /// `SM_CP_CENTRIPETAL`.
    Centripetal,
}

/// An interpolating curve and the vectors as SMLib leaves them.
#[derive(Clone, Debug, PartialEq)]
pub struct Interpolation {
    /// The curve.
    pub curve: BSplineCurve,
    /// The input vectors after the call. With `ignore_vector_magnitude`
    /// SMLib normalises the caller's array and scales it by the length of
    /// the data polygon; otherwise they are unchanged.
    pub vectors: Vec<Vector3>,
}

fn distance(p: Vector3, q: Vector3) -> f64 {
    ((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) + (p.z - q.z) * (p.z - q.z)).sqrt()
}

/// `N_FitCalcCrvParamValues` for Euclidean points.
fn parameter_values(
    points: &[Vector3],
    parameterization: Parameterization,
) -> Result<Vec<f64>, Fault> {
    let n = points.len() - 1;
    let mut u = vec![0.0; n + 1];
    u[n] = 1.0;
    if parameterization == Parameterization::Uniform {
        let factor = 1.0 / n as f64;
        for i in 1..n {
            u[i] = i as f64 * factor;
        }
        return Ok(u);
    }
    let mut dist = vec![0.0; n + 1];
    let mut sum = 0.0;
    for i in 1..=n {
        dist[i] = distance(points[i - 1], points[i]);
        if parameterization == Parameterization::Centripetal {
            dist[i] = dist[i].sqrt();
        }
        sum += dist[i];
    }
    // When the data have no length SMLib writes no interior parameter. Its
    // allocator returns zeroed memory in the pinned build, so they are zero.
    if sum > PARAMETER_TOLERANCE {
        for i in 1..n {
            u[i] = u[i - 1] + dist[i] / sum;
        }
    }
    Ok(u)
}

/// `N_FitCalcKnotsDerivs`: the `2k + p + 3` knots.
fn derivative_knots(u: &[f64], k: usize, p: usize) -> Vec<f64> {
    let n = 2 * k + 1;
    let mut knots = vec![0.0; n + p + 2];
    for i in 0..=p {
        knots[i] = u[0];
        knots[n + i + 1] = u[k];
    }
    let mut j = p + 1;
    match p {
        2 => {
            knots[p + 1] = 0.5 * (u[1] + u[0]);
            for i in 1..k {
                j += 1;
                knots[j] = u[i];
                j += 1;
                knots[j] = 0.5 * (u[i] + u[i + 1]);
            }
        }
        3 => {
            if k > 1 {
                knots[p + 1] = 0.5 * (u[1] + u[0]);
                knots[n] = 0.5 * (u[k] + u[k - 1]);
            }
            let (alpha, beta) = (1.0 / 3.0, 2.0 / 3.0);
            for i in 1..k.saturating_sub(1) {
                j += 1;
                knots[j] = beta * u[i] + alpha * u[i + 1];
                j += 1;
                knots[j] = alpha * u[i] + beta * u[i + 1];
            }
        }
        _ => {}
    }
    knots
}

/// A banded `NL_RMATRIX` in SMLib's block layout: `rows` rows of `width`
/// entries, the entry for column `j` of row `i` at offset `j - i + width / 2`.
struct Band {
    values: Vec<f64>,
    rows: usize,
    width: usize,
    half: usize,
}

impl Band {
    fn new(rows: usize, width: usize) -> Self {
        Self {
            values: vec![0.0; rows * width],
            rows,
            width,
            half: width / 2,
        }
    }

    /// The flat position of `row`, offset `offset` (which may spill into a
    /// neighbouring row, as SMLib's pointer arithmetic does).
    fn slot(&self, row: usize, offset: i64) -> Result<usize, Fault> {
        let flat = (row * self.width) as i64 + offset;
        if flat < 0 || flat as usize >= self.values.len() {
            return Err(Fault::Degenerate);
        }
        Ok(flat as usize)
    }

    fn at(&self, row: usize, offset: i64) -> f64 {
        self.values[self
            .slot(row, offset)
            .expect("banded entry inside the block")]
    }

    fn set(&mut self, row: usize, offset: i64, value: f64) -> Result<(), Fault> {
        let slot = self.slot(row, offset)?;
        self.values[slot] = value;
        Ok(())
    }

    /// `N_RealMatrixLuDecompose` for a banded matrix.
    fn decompose(&mut self) -> Result<(), Fault> {
        let n = self.rows as i64 - 1;
        let half = self.half as i64;
        for i in 0..=n {
            let high = n.min(i + half);
            for j in i..=high {
                let low = 0.max(j - half);
                let mut sum = 0.0;
                for k in low..i {
                    sum += self.at(i as usize, k - i + half) * self.at(k as usize, j - k + half);
                }
                let slot = self.slot(i as usize, j - i + half)?;
                self.values[slot] -= sum;
            }
            for j in i + 1..=high {
                let low = 0.max(j - half);
                let mut sum = 0.0;
                for k in low..i {
                    sum += self.at(j as usize, k - j + half) * self.at(k as usize, i - k + half);
                }
                if self.at(i as usize, half).abs() < LU_PIVOT_TOLERANCE {
                    return Err(Fault::Degenerate);
                }
                let value = (self.at(j as usize, i - j + half) - sum) / self.at(i as usize, half);
                let slot = self.slot(j as usize, i - j + half)?;
                self.values[slot] = value;
            }
        }
        Ok(())
    }

    /// `N_RealMatrixForBack` for Euclidean points.
    fn solve(&self, rhs: &[Vector3]) -> Result<Vec<Vector3>, Fault> {
        let n = self.rows as i64 - 1;
        let half = self.half as i64;
        let mut z = vec![Vector3::new(0.0, 0.0, 0.0); self.rows];
        for i in 0..=n {
            let mut s = Vector3::new(0.0, 0.0, 0.0);
            for j in 0.max(i - half)..i {
                let b = self.at(i as usize, j - i + half);
                let zj = z[j as usize];
                s = Vector3::new(s.x + b * zj.x, s.y + b * zj.y, s.z + b * zj.z);
            }
            let p = rhs[i as usize];
            z[i as usize] = Vector3::new(p.x - s.x, p.y - s.y, p.z - s.z);
        }
        let mut q = vec![Vector3::new(0.0, 0.0, 0.0); self.rows];
        for i in (0..=n).rev() {
            let pivot = self.at(i as usize, half);
            if division_is_bad(1.0, pivot) {
                return Err(Fault::Degenerate);
            }
            let factor = 1.0 / pivot;
            let mut s = Vector3::new(0.0, 0.0, 0.0);
            for j in i + 1..=n.min(i + half) {
                let b = self.at(i as usize, j - i + half);
                let qj = q[j as usize];
                s = Vector3::new(s.x + b * qj.x, s.y + b * qj.y, s.z + b * qj.z);
            }
            let zi = z[i as usize];
            q[i as usize] = Vector3::new(
                factor * zi.x + -factor * s.x,
                factor * zi.y + -factor * s.y,
                factor * zi.z + -factor * s.z,
            );
        }
        Ok(q)
    }
}

/// `N_FitCrvFirstDeriv`: knots and control points of the interpolating
/// curve. `vectors` is rewritten as SMLib rewrites it for `tangents`.
fn fit_first_derivatives(
    points: &[Vector3],
    vectors: &mut [Vector3],
    p: usize,
    tangents: bool,
    parameterization: Parameterization,
) -> Result<(Vec<f64>, Vec<Vector3>), Fault> {
    if !(2..=3).contains(&p) {
        return Err(Fault::Degenerate);
    }
    let kk = points.len() - 1;
    let n = 2 * kk + 1;
    let u = parameter_values(points, parameterization)?;
    let knots = derivative_knots(&u, kk, p);
    let width = 2 * p - 1;
    let half = (p - 1) as i64;
    let mut a = Band::new(n + 1, width);
    let sbw = half;
    a.set(0, sbw, 1.0)?;
    a.set(1, sbw - 1, -1.0)?;
    a.set(1, sbw, 1.0)?;
    a.set(n - 1, sbw, -1.0)?;
    a.set(n - 1, sbw + 1, 1.0)?;
    a.set(n, sbw, 1.0)?;
    let high = if p == 3 { p } else { p - 1 };
    for i in 1..kk {
        let (span, derivatives) = basis_derivatives(&knots, p, u[i], false, 1)?;
        let row = 2 * i;
        let l = span as i64 - row as i64 - 1;
        for k in 0..=high {
            a.set(row, l + k as i64, derivatives[0][k])?;
            a.set(row + 1, l + k as i64 - 1, derivatives[1][k])?;
        }
    }
    a.decompose()?;
    if tangents {
        let length: f64 = (1..=kk)
            .map(|i| distance(points[i - 1], points[i]))
            .fold(0.0, |s, d| s + d);
        for vector in vectors.iter_mut() {
            // `N_VectorNormalizeRef`, then `N_VectorScale`.
            let magnitude =
                (vector.x * vector.x + vector.y * vector.y + vector.z * vector.z).sqrt();
            if division_is_bad(1.0, magnitude) {
                return Err(Fault::Degenerate);
            }
            let inverse = 1.0 / magnitude;
            let unit = Vector3::new(inverse * vector.x, inverse * vector.y, inverse * vector.z);
            *vector = Vector3::new(length * unit.x, length * unit.y, length * unit.z);
        }
    }
    let scale = |factor: f64, v: Vector3| Vector3::new(factor * v.x, factor * v.y, factor * v.z);
    let mut rhs = vec![Vector3::new(0.0, 0.0, 0.0); n + 1];
    rhs[0] = points[0];
    rhs[1] = scale(knots[p + 1] / p as f64, vectors[0]);
    for i in 1..kk {
        rhs[2 * i] = points[i];
        rhs[2 * i + 1] = vectors[i];
    }
    rhs[n - 1] = scale((1.0 - knots[n]) / p as f64, vectors[kk]);
    rhs[n] = points[kk];
    Ok((knots, a.solve(&rhs)?))
}

impl BSplineCurve {
    /// `SmBSplineCurve::CreateInterpolatingCurve` for points with a vector
    /// at each (more than two points), and for a degree-one curve through two
    /// points. The vectors are first derivatives, or with
    /// `ignore_vector_magnitude` directions that SMLib scales by the length of
    /// the data polygon. A planar result has the planar marker; SMLib keeps
    /// the points' z and labels the curve two-dimensional.
    ///
    /// Every other combination of vectors needs `N_FitCrvHighDerivs`, which is
    /// not ported, and refuses with [`Fault::Degenerate`].
    pub fn interpolating_curve(
        parameterization: Parameterization,
        dimension: usize,
        degree: usize,
        points: &[Vector3],
        vectors: &[Vector3],
        ignore_vector_magnitude: bool,
    ) -> Result<Interpolation, Fault> {
        if degree == 1 && points.len() == 2 {
            let curve = BSplineCurve::line_segment(dimension, points[0], points[1])?;
            return Ok(Interpolation {
                curve,
                vectors: vectors.to_vec(),
            });
        }
        if dimension != 2 && dimension != 3 {
            return Err(Fault::Degenerate);
        }
        if vectors.len() != points.len() || points.len() <= 2 {
            // The end-vector, start-vector and no-vector branches.
            return Err(Fault::Degenerate);
        }
        if degree != 2 && degree != 3 {
            return Err(Fault::Degenerate);
        }
        let mut vectors = vectors.to_vec();
        let (knots, controls) = fit_first_derivatives(
            points,
            &mut vectors,
            degree,
            ignore_vector_magnitude,
            parameterization,
        )?;
        let points = controls
            .into_iter()
            .map(|c| ControlPoint {
                coordinates: [c.x, c.y, if dimension == 2 { NO_Z } else { c.z }, NO_WEIGHT],
            })
            .collect();
        let curve = constructed(degree, knots, points)?;
        Ok(Interpolation { curve, vectors })
    }
}

/// The tangents and higher derivatives after [`fair_blend_curve_derivatives`].
#[derive(Clone, Debug, PartialEq)]
pub struct FairedDerivatives {
    /// The two end tangents, scaled.
    pub tangents: [Vector3; 2],
    /// The higher-order derivatives, scaled when there are two or four.
    pub higher: Vec<Vector3>,
}

/// `smgu_LinePlaneIntersect`: the parameter at which `start + t * direction`
/// meets the plane through `plane_point` with normal `normal`.
fn line_plane_intersect(
    start: Vector3,
    direction: Vector3,
    plane_point: Vector3,
    normal: Vector3,
) -> Result<f64, Fault> {
    // `ZERO_VEC_ER`: a NaN length refuses too.
    for vector in [direction, normal] {
        let squared = vector.length_squared();
        if squared.is_nan() || squared < EFF_ZERO_SQ {
            return Err(Fault::InvalidInput);
        }
    }
    let plane_d = -normal.dot(plane_point);
    let tolerance = EFF_ZERO * direction.max_dimension() + EFF_ZERO_SQ;
    let denominator = direction.dot(normal);
    if denominator.abs() < tolerance {
        return Err(Fault::Degenerate);
    }
    Ok(-(plane_d + start.dot(normal)) / denominator)
}

/// `SmBSplineCurve::FairBlendCurveDerivatives`: scales the end tangents of a
/// blend between two points so the curve's control polygon is evenly spaced,
/// optionally keeping it free of inflections, and scales `higher` to match.
pub fn fair_blend_curve_derivatives(
    points: [Vector3; 2],
    tangents: [Vector3; 2],
    higher: &[Vector3],
    prevent_inflection: bool,
) -> Result<FairedDerivatives, Fault> {
    let (mut scale1, mut scale2) = (1.0, 1.0);
    let a = points[1] - points[0];
    let mut v1 = -tangents[1];
    v1.unitize()?;
    let mut v0 = tangents[0];
    v0.unitize()?;
    let b = v1 - v0;
    let coefficients = [a.dot(a), 2.0 * a.dot(b), b.dot(b) - 1.0];
    let (count, roots) = solve_quadratic(coefficients, EFF_ZERO);
    let mut t = 0.0;
    if count == 1 {
        t = 3.0 * roots[0];
    }
    if count == 2 {
        t = if roots[0] < 0.0 {
            3.0 * roots[1]
        } else if roots[1] < 0.0 {
            3.0 * roots[0]
        } else {
            3.0 * if roots[0] < roots[1] {
                roots[0]
            } else {
                roots[1]
            }
        };
    }
    if t > 0.0 {
        scale1 = t / tangents[0].length();
        scale2 = t / tangents[1].length();
    }
    if prevent_inflection {
        let z = tangents[0].cross(tangents[1]);
        let mut normal1 = z.cross(tangents[1]);
        let mut offset = points[0] - points[1];
        normal1.unitize()?;
        offset.unitize()?;
        let dot1 = offset.dot(normal1);
        let control1 = points[0] + tangents[0] * scale1 / 3.0;
        let line1 = control1 - points[0];
        let to_second = control1 - points[1];
        let dot2 = to_second.dot(normal1);
        let (mut t1, mut t2) = (1.0, 1.0);
        if dot1 * dot2 < 0.0 {
            t1 = line_plane_intersect(points[0], line1, points[1], normal1)?;
            if t1 < 0.03 {
                t1 = 1.0;
            }
        }
        let mut normal2 = z.cross(tangents[0]);
        normal2.unitize()?;
        offset = -offset;
        let dot3 = offset.dot(normal2);
        let control2 = points[1] - tangents[1] * scale2 / 3.0;
        let to_first = control2 - points[0];
        let dot4 = to_first.dot(normal2);
        let line2 = control2 - points[1];
        if dot3 * dot4 < 0.0 {
            t2 = line_plane_intersect(points[1], line2, points[0], normal2)?;
            if t2 < 0.03 {
                t2 = 1.0;
            }
        }
        t = if t1 < t2 { t1 } else { t2 };
        scale1 *= t;
        scale2 *= t;
    }
    let scaled = [tangents[0] * scale1, tangents[1] * scale2];
    let mut higher = higher.to_vec();
    if higher.len() == 2 {
        higher[0] = higher[0] * scale1 * scale1;
        higher[1] = higher[1] * scale2 * scale2;
    }
    if higher.len() == 4 {
        higher[0] = higher[0] * scale1 * scale1;
        higher[1] = higher[1] * scale2 * scale2;
        higher[2] = higher[2] * scale1 * scale1 * scale1;
        higher[3] = higher[3] * scale2 * scale2 * scale2;
    }
    Ok(FairedDerivatives {
        tangents: scaled,
        higher,
    })
}
