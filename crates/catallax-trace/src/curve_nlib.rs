// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: NLib stacks and in-place aliasing become owned
// values and statuses become `Fault`. Translated directly from SMLib
// NL_CrvNurb.cpp (`N_CrvDerivs`), NL_CrvGeom.cpp (`N_CrvClosestPt`),
// NL_CrvTool.cpp (`N_CrvExtractCrvSeg`), NL_CrvBasic.cpp (`N_CrvIsClosed`)
// and NL_Math.cpp (`N_FloatOpIsBad`). rust-brep has no port as a unit.
//! NLib's own curve evaluation, closest-point search and segment extraction,
//! which SMLib's arc builders call directly (not through `SmBSplineCurve`).
//!
//! Every result is an unproved candidate.

#![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]

use crate::{
    Fault,
    curve_build::span_and_multiplicity,
    curve_compat::division_is_bad,
    spline::{BSplineCurve, ControlPoint, NO_WEIGHT, NO_Z, combine, find_span},
    vector::Vector3,
};

/// `NL_MAXDER`, the highest derivative NLib evaluates.
const MAX_DERIVATIVE: usize = 4;
/// `NL_ITLIM`, the iteration limit of the closest-point Newton search.
const ITERATION_LIMIT: usize = 50;
/// `NL_MTOL`, the model-space tolerance.
pub(crate) const MODEL_TOLERANCE: f64 = 1.0e-7;

/// `N_CPtToPtEuclid`: the Euclidean point, by the weight's reciprocal.
fn euclid(point: ControlPoint) -> Vector3 {
    let [x, y, z, w] = point.coordinates;
    let z_or_zero = |value: f64| if z == NO_Z { 0.0 } else { value };
    if w == NO_WEIGHT {
        Vector3::new(x, y, z_or_zero(z))
    } else {
        let reciprocal = 1.0 / w;
        Vector3::new(x * reciprocal, y * reciprocal, z_or_zero(z * reciprocal))
    }
}

/// `N_CrvDerivs`: the point and derivatives of the curve at `u`, spans taken
/// from the left or the right. Refuses a parameter outside the knots and
/// more than seventeen derivatives.
pub fn curve_derivatives(
    curve: &BSplineCurve,
    u: f64,
    from_left: bool,
    order: usize,
) -> Result<Vec<Vector3>, Fault> {
    let p = curve.degree();
    let knots = curve.knots();
    let points = curve.control_points();
    let m = knots.len() - 1;
    if u < knots[0] || u > knots[m] || order > MAX_DERIVATIVE {
        return Err(Fault::Degenerate);
    }
    // `NL_LEFT` is the span search without the right-continuous flag.
    let span = find_span(knots, p, u, !from_left) as i64;
    let degree = p as i64;
    let mut ndu = vec![vec![0.0; p + 1]; p + 1];
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = u - knots[(span + 1 - j as i64) as usize];
        right[j] = knots[span as usize + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            if division_is_bad(ndu[r][j - 1], ndu[j][r]) {
                return Err(Fault::Degenerate);
            }
            let temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    let mut basis = vec![vec![0.0; p + 1]; order.max(p) + 1];
    for j in 0..=p {
        basis[0][j] = ndu[j][p];
    }
    let mder = p.min(order);
    let mut a = [vec![0.0; p + 2], vec![0.0; p + 2]];
    for r in 0..=p {
        let (mut s1, mut s2) = (0usize, 1usize);
        a[0][0] = 1.0;
        for k in 1..=mder {
            let mut d = 0.0;
            let rk = r as i64 - k as i64;
            let pk = degree - k as i64;
            if r >= k {
                a[s2][0] = a[s1][0] / ndu[(pk + 1) as usize][rk as usize];
                d = a[s2][0] * ndu[rk as usize][pk as usize];
            }
            let j1 = if rk >= -1 { 1 } else { -rk };
            let j2 = if r as i64 - 1 <= pk {
                k as i64 - 1
            } else {
                degree - r as i64
            };
            for j in j1..=j2 {
                let ju = j as usize;
                a[s2][ju] = (a[s1][ju] - a[s1][ju - 1]) / ndu[(pk + 1) as usize][(rk + j) as usize];
                d += a[s2][ju] * ndu[(rk + j) as usize][pk as usize];
            }
            if r as i64 <= pk {
                a[s2][k] = -a[s1][k - 1] / ndu[(pk + 1) as usize][r];
                d += a[s2][k] * ndu[r][pk as usize];
            }
            basis[k][r] = d;
            std::mem::swap(&mut s1, &mut s2);
        }
    }
    let mut factor = p as i64;
    for k in 1..=mder {
        for j in 0..=p {
            basis[k][j] *= factor as f64;
        }
        factor *= degree - k as i64;
    }
    let first = (span - degree) as usize;
    if points[0].coordinates[3] != NO_WEIGHT {
        let mut dd = vec![0.0; order + 1];
        for k in 0..=order {
            for j in first..=span as usize {
                dd[k] += points[j].coordinates[3] * basis[k][j - first];
            }
        }
        // Pascal's triangle by NLib's recurrence.
        let mut tri = vec![vec![0i64; order + 2]; order + 2];
        tri[0][0] = 1;
        if order >= 1 {
            tri[1][0] = 1;
            tri[1][1] = 1;
        }
        if order > 1 {
            tri[2][0] = 1;
            tri[2][1] = 2;
            tri[2][2] = 1;
            for k in 3..=order {
                tri[k][0] = 1;
                let half = k / 2;
                let mut previous = 1;
                for j in 1..=half {
                    let held = tri[k - 1][j];
                    tri[k][j] = tri[k - 1][j] + previous;
                    tri[k][k - j] = tri[k][j];
                    previous = held;
                }
                tri[k][k] = 1;
            }
        }
        for i in 0..=p {
            let w = points[first + i].coordinates[3];
            for k in 0..=order {
                let mut d = w * basis[k][i];
                for j in 1..=k {
                    d -= tri[k][j] as f64 * dd[j] * basis[k - j][i];
                }
                basis[k][i] = d / dd[0];
            }
        }
    }
    let euclidean: Vec<Vector3> = (first..=span as usize).map(|i| euclid(points[i])).collect();
    let mut result = Vec::with_capacity(order + 1);
    for k in 0..=order {
        let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
        for i in 0..=p {
            let d = basis[k][i];
            x += d * euclidean[i].x;
            y += d * euclidean[i].y;
            z += d * euclidean[i].z;
        }
        result.push(Vector3::new(x, y, z));
    }
    Ok(result)
}

/// `N_CrvIsClosed`: the end points lie within `NL_MTOL`.
pub fn curve_is_closed(curve: &BSplineCurve) -> bool {
    let points = curve.control_points();
    let (first, last) = (euclid(points[0]), euclid(points[points.len() - 1]));
    let distance = ((first.x - last.x) * (first.x - last.x)
        + (first.y - last.y) * (first.y - last.y)
        + (first.z - last.z) * (first.z - last.z))
        .sqrt();
    distance < MODEL_TOLERANCE
}

/// `N_CrvClosestPt`: Newton's search from `start` for the parameter of the
/// point on the curve nearest `point`, stopping when the distance is within
/// `distance_tolerance` or the tangent is within `cosine_tolerance` of
/// perpendicular to the offset. Returns the parameter and the curve point.
pub fn closest_point_on_curve(
    curve: &BSplineCurve,
    point: Vector3,
    start: f64,
    distance_tolerance: f64,
    cosine_tolerance: f64,
) -> Result<(f64, Vector3), Fault> {
    let knots = curve.knots();
    let m = knots.len() - 1;
    let closed = curve_is_closed(curve);
    let mut unew = start;
    let mut best_distance = f64::MAX;
    let mut crossings = 0;
    let (mut found, mut found_point) = (0.0, Vector3::new(0.0, 0.0, 0.0));
    let mut iteration = 0;
    while iteration < ITERATION_LIMIT {
        let d = curve_derivatives(curve, unew, true, 2)?;
        let v = d[0] - point;
        let numerator = d[1].dot(v);
        let derivative = d[1].length_squared().sqrt();
        let distance = v.length_squared().sqrt();
        if distance < best_distance {
            found = unew;
            found_point = d[0];
            best_distance = distance;
        }
        if distance <= distance_tolerance {
            break;
        }
        if division_is_bad(numerator, derivative * distance) {
            return Err(Fault::Degenerate);
        }
        let cosine = numerator.abs() / (derivative * distance);
        if cosine <= cosine_tolerance && distance == best_distance {
            break;
        }
        let dot = d[2].dot(v);
        let denominator = dot + derivative * derivative;
        if division_is_bad(numerator, denominator) {
            return Err(Fault::Degenerate);
        }
        let uold = unew;
        unew = uold - numerator / denominator;
        if closed {
            let range = knots[m] - knots[0];
            if unew < knots[0] {
                if knots[0] - unew < 100.0 * range {
                    while unew < knots[0] {
                        unew = knots[m] - knots[0] + unew;
                    }
                    crossings = 0;
                } else {
                    if crossings >= 2 {
                        return Err(Fault::Degenerate);
                    }
                    unew = uold;
                    crossings += 1;
                }
            } else if unew > knots[m] {
                if unew - knots[m] < 100.0 * range {
                    while unew > knots[m] {
                        unew = knots[0] - knots[m] + unew;
                    }
                    crossings = 0;
                } else {
                    if crossings >= 2 {
                        return Err(Fault::Degenerate);
                    }
                    unew = uold;
                    crossings += 1;
                }
            } else {
                crossings = 0;
            }
        } else {
            if unew < knots[0] {
                unew = knots[0];
            }
            if unew > knots[m] {
                unew = knots[m];
            }
        }
        if ((unew - uold) * derivative).abs() <= distance_tolerance && distance == best_distance {
            break;
        }
        iteration += 1;
    }
    if iteration >= ITERATION_LIMIT {
        return Err(Fault::Degenerate);
    }
    Ok((found, found_point))
}

/// `N_CrvExtractCrvSeg`: the part of the curve between `left` and `right`,
/// on a clamped knot vector, by knot insertion.
pub fn extract_curve_segment(
    curve: &BSplineCurve,
    left: f64,
    right: f64,
) -> Result<BSplineCurve, Fault> {
    // A NaN bound makes SMLib's comparisons all fail; so it does here.
    if right <= left {
        return Err(Fault::Degenerate);
    }
    let (p, knots, points) = (curve.degree(), curve.knots(), curve.control_points());
    let degree = p as i64;
    let m = knots.len() as i64 - 1;
    let (spl, mll) = span_and_multiplicity(knots, p, left);
    let (mut spr, mut mlr) = span_and_multiplicity(knots, p, right);
    let first = spl - degree;
    if right == knots[(m - degree) as usize] {
        spr = m;
        mlr = degree + 1;
    }
    let last = spr - mlr;
    let n = last - first;
    let count = spr - spl - mlr + 2 * degree + 1;
    if n < 0 {
        return Err(Fault::Degenerate);
    }
    let mut q: Vec<ControlPoint> = (first..=last).map(|i| points[i as usize]).collect();
    let ll = spl - degree;
    for i in 1..=degree - mll {
        for j in 0..=degree - i - mll {
            let leftmost = knots[(ll + i + j) as usize];
            let alpha = (left - leftmost) / (knots[(spl + j + 1) as usize] - leftmost);
            q[j as usize] = combine(alpha, q[j as usize + 1], 1.0 - alpha, q[j as usize]);
        }
    }
    let lr = spr - degree;
    let lk = n - degree + mlr;
    for i in 1..=degree - mlr {
        for j in (0..=degree - i - mlr).rev() {
            let k = lk + i + j;
            let mut leftmost = knots[(lr + i + j) as usize];
            if leftmost < left {
                leftmost = left;
            }
            let alpha = (right - leftmost) / (knots[(spr + j + 1) as usize] - leftmost);
            q[k as usize] = combine(alpha, q[k as usize], 1.0 - alpha, q[k as usize - 1]);
        }
    }
    let mut uq = Vec::with_capacity(count as usize + 1);
    uq.extend(std::iter::repeat_n(left, p + 1));
    uq.extend((spl + 1..=spr - mlr).map(|i| knots[i as usize]));
    uq.extend(std::iter::repeat_n(right, p + 1));
    BSplineCurve::new(p, uq, q)
}
