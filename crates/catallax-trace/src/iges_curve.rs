// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: the curve storage is an owned `RawCurve`, error
// codes saved in SMLib's global error state become `Fault`, and the C++
// memory stacks are dropped. Reads that the C++ would make outside its
// arrays are refused instead of performed.
// Adapted from rust-brep crates/kernel/src/iges_crv_plane_sweep_rotation.rs,
// commit 0c40280847447020b9371410cd03393b8f07fb5f, a translation of SMLib
// NL_Iges.cpp (N_Iges126Crv, N_IgesClampCrv, N_IgesCorrectCrvKnots,
// N_IgesValidate) at c9c5979f04cea8c38213004d4f04df0dcfb40d05.
//! IGES entity 126 curve construction: validate, clamp to the bounds, and
//! correct the knot multiplicities.
//!
//! The result of [`iges_126_curve`] is a clamped curve with no interior knot
//! of multiplicity above its degree and no end knot above `degree + 1`. It is
//! an unproved candidate like everything in this crate.

use crate::Fault;
use crate::spline::{ControlPoint, NO_WEIGHT, NO_Z, combine};
use crate::vector::Vector3;

/// `NL_DMAX` as `N_InitNurbs` sets it: the largest degree accepted.
pub const MAX_DEGREE: usize = 17;
/// `NL_WMIN`: the smallest weight accepted.
pub const MIN_WEIGHT: f64 = 1.0e-3;
/// `NL_WMAX`: the largest weight accepted.
pub const MAX_WEIGHT: f64 = 100.0;

/// SMLib's `NL_CURVE`: a degree, a knot vector and a control polygon, not
/// yet checked against the relations a valid B-spline curve satisfies.
#[derive(Clone, Debug, PartialEq)]
pub struct RawCurve {
    /// Polynomial degree.
    pub degree: usize,
    /// Every knot, repeats included.
    pub knots: Vec<f64>,
    /// Stored control points, weighted when the curve is rational.
    pub points: Vec<ControlPoint>,
}

/// A knot or control point by signed index, refusing an index outside the
/// array where the C++ would read beyond it.
fn at<T: Copy>(values: &[T], index: isize) -> Result<T, Fault> {
    usize::try_from(index)
        .ok()
        .and_then(|i| values.get(i).copied())
        .ok_or(Fault::InvalidInput)
}

/// Writes `value` at a signed index, refusing an index outside the array.
fn put<T>(values: &mut [T], index: isize, value: T) -> Result<(), Fault> {
    let slot = usize::try_from(index)
        .ok()
        .and_then(|i| values.get_mut(i))
        .ok_or(Fault::InvalidInput)?;
    *slot = value;
    Ok(())
}

/// `N_Sum2CPts`: the sum of two control points, with SMLib's handling of the
/// missing-coordinate markers.
fn sum_points(a: ControlPoint, b: ControlPoint) -> ControlPoint {
    let (a, b) = (a.coordinates, b.coordinates);
    let z = match (a[2] != NO_Z, b[2] != NO_Z) {
        (true, true) => a[2] + b[2],
        (true, false) => a[2],
        (false, true) => b[2],
        (false, false) => NO_Z,
    };
    let w = if a[3] != NO_WEIGHT && b[3] != NO_WEIGHT {
        a[3] + b[3]
    } else {
        NO_WEIGHT
    };
    ControlPoint {
        coordinates: [a[0] + b[0], a[1] + b[1], z, w],
    }
}

/// `N_ScaleCPt`: a control point scaled by `alpha`, markers kept.
fn scale_point(alpha: f64, a: ControlPoint) -> ControlPoint {
    let [x, y, z, w] = a.coordinates;
    ControlPoint {
        coordinates: [
            alpha * x,
            alpha * y,
            if z != NO_Z { alpha * z } else { NO_Z },
            if w != NO_WEIGHT { alpha * w } else { NO_WEIGHT },
        ],
    }
}

/// `N_IgesValidate`: the bounds lie inside the active knot range, in
/// increasing order, and the knots do not decrease.
fn validate(knots: &[f64], degree: usize, bounds: [f64; 2]) -> Result<(), Fault> {
    let last = knots.len() - 1;
    if knots[degree] > bounds[0] || bounds[0] >= bounds[1] || knots[last - degree] < bounds[1] {
        return Err(Fault::Degenerate);
    }
    if knots.windows(2).any(|pair| pair[0] > pair[1]) {
        return Err(Fault::Degenerate);
    }
    Ok(())
}

/// `N_Iges126Crv`: the clamped curve of IGES entity 126 data.
///
/// `points` are Euclidean and may carry [`NO_Z`] in every `z` for a planar
/// curve; `weights` are one per point. Equal weights make a non-rational
/// curve whose coordinates are left unscaled, as SMLib does. A degree above
/// [`MAX_DEGREE`], a weight outside `[MIN_WEIGHT, MAX_WEIGHT]` or bounds
/// outside the knot range return [`Fault::Degenerate`], SMLib's `SM_ERR`.
/// Counts that do not satisfy `knots = points + degree + 1`, and non-finite
/// bounds, are [`Fault::InvalidInput`]; SMLib reads past its arrays there.
pub fn iges_126_curve(
    degree: usize,
    knots: &[f64],
    weights: &[f64],
    points: &[Vector3],
    bounds: [f64; 2],
) -> Result<RawCurve, Fault> {
    if points.is_empty()
        || weights.len() != points.len()
        || knots.len() != points.len() + degree + 1
        || !bounds.iter().all(|v| v.is_finite())
    {
        return Err(Fault::InvalidInput);
    }
    if !(1..=MAX_DEGREE).contains(&degree) {
        return Err(Fault::Degenerate);
    }
    let mut rational = false;
    for &w in weights {
        if w < MIN_WEIGHT || w > MAX_WEIGHT {
            return Err(Fault::Degenerate);
        }
        if w != weights[0] {
            rational = true;
        }
    }
    validate(knots, degree, bounds)?;
    let control = points
        .iter()
        .zip(weights)
        .map(|(point, &w)| {
            let (x, y, z) = (point.x, point.y, point.z);
            let coordinates = if rational {
                [w * x, w * y, if z == NO_Z { z } else { w * z }, w]
            } else {
                [x, y, z, NO_WEIGHT]
            };
            ControlPoint { coordinates }
        })
        .collect();
    let mut curve = RawCurve {
        degree,
        knots: knots.to_vec(),
        points: control,
    };
    clamp_curve(&mut curve, bounds)?;
    correct_curve_knots(&mut curve)?;
    Ok(curve)
}

/// `N_IgesClampCrv`: clamps `curve` in place to `bounds`, which must satisfy
/// `U[p] <= bounds[0] < bounds[1] <= U[m-p]` for a nondecreasing knot vector.
pub fn clamp_curve(curve: &mut RawCurve, bounds: [f64; 2]) -> Result<(), Fault> {
    if curve.knots.is_empty() || curve.points.is_empty() {
        return Err(Fault::InvalidInput);
    }
    let p = curve.degree as isize;
    let m = curve.knots.len() as isize - 1;
    let mut up = curve.knots.clone();
    let pw = curve.points.clone();
    let [v0, v1] = bounds;

    // Spans and multiplicities of the bounds, as knots.
    let (mut span1, mut mult1) = (0isize, 0isize);
    for i in 0..=m {
        let u = up[i as usize];
        if u == v0 {
            mult1 += 1;
        }
        if u > v0 {
            span1 = i - 1;
            break;
        }
    }
    let (mut span2, mut mult2) = (0isize, 0isize);
    for i in (0..=m).rev() {
        let u = up[i as usize];
        if u == v1 {
            mult2 += 1;
        }
        if u < v1 {
            span2 = i + mult2;
            break;
        }
    }

    // Junk in the first and last knots. The writes land in the snapshot the
    // insertion loops read, and in the caller's curve.
    if span1 >= p && mult1 == span1 {
        up[0] = v0;
        mult1 += 1;
    }
    if span2 == m - 1 && mult2 >= p {
        up[m as usize] = v1;
        mult2 += 1;
        span2 += 1;
    }
    curve.knots.clone_from(&up);

    // Already clamped at both ends.
    if span1 >= p && mult1 == span1 + 1 && span2 == m && mult2 > p {
        return Ok(());
    }

    let start = span1 - p;
    let end = span2 - mult2;
    let n = end - start;
    let new_m = span2 - span1 - mult2 + 2 * p + 1;
    if n < 0 {
        return Err(Fault::Degenerate);
    }

    let mut qw: Vec<ControlPoint> = (start..=end)
        .map(|i| at(&pw, i))
        .collect::<Result<_, _>>()?;

    // Insert the left knot.
    let ll = span1 - p;
    for i in 1..=p - mult1 {
        for j in 0..=p - i - mult1 {
            let left = at(&up, ll + i + j)?;
            let alpha = (v0 - left) / (at(&up, span1 + j + 1)? - left);
            let (next, here) = (at(&qw, j + 1)?, at(&qw, j)?);
            put(&mut qw, j, combine(alpha, next, 1.0 - alpha, here))?;
        }
    }

    // Insert the right knot.
    let lr = span2 - p;
    let lk = n - p + mult2;
    for i in 1..=p - mult2 {
        for j in (0..=p - i - mult2).rev() {
            let k = lk + i + j;
            let mut left = at(&up, lr + i + j)?;
            if left < v0 {
                left = v0;
            }
            let alpha = (v1 - left) / (at(&up, span2 + j + 1)? - left);
            let (here, below) = (at(&qw, k)?, at(&qw, k - 1)?);
            put(&mut qw, k, combine(alpha, here, 1.0 - alpha, below))?;
        }
    }

    // The knot vector: bounds with full multiplicity around the interior.
    let mut uq = vec![v0; p as usize + 1];
    for i in span1 + 1..=span2 - mult2 {
        uq.push(at(&up, i)?);
    }
    uq.extend(std::iter::repeat_n(v1, p as usize + 1));
    if uq.len() as isize != new_m + 1 {
        return Err(Fault::Degenerate);
    }
    curve.points = qw;
    curve.knots = uq;
    Ok(())
}

/// `N_IgesCorrectCrvKnots`: reduces an interior knot of multiplicity above
/// `p` to `p`, and an end knot above `p + 1` to `p + 1`, merging the
/// control points that the removed knots carried.
pub fn correct_curve_knots(curve: &mut RawCurve) -> Result<(), Fault> {
    let p = curve.degree as isize;
    let n0 = curve.points.len() as isize - 1;
    let mut nn = n0;
    let mut mm = curve.knots.len() as isize - 1;
    let mut pw = std::mem::take(&mut curve.points);
    let mut up = std::mem::take(&mut curve.knots);
    let result = correct(&mut pw, &mut up, p, &mut nn, &mut mm);
    if nn < n0 {
        pw.truncate(usize::try_from(nn + 1).unwrap_or(0));
        up.truncate(usize::try_from(mm + 1).unwrap_or(0));
    }
    curve.points = pw;
    curve.knots = up;
    result
}

fn correct(
    pw: &mut [ControlPoint],
    up: &mut [f64],
    p: isize,
    nn: &mut isize,
    mm: &mut isize,
) -> Result<(), Fault> {
    let us = at(up, 0)?;
    let ue = at(up, *mm)?;

    // Interior multiplicities above p.
    let mut ii = p + 1;
    while at(up, ii)? == at(up, ii - 1)? {
        ii += 1;
    }
    while at(up, ii)? < ue {
        let mut mult = 1;
        while at(up, ii)? == at(up, ii + mult)? {
            mult += 1;
        }
        if mult > p {
            let mut k1 = ii - 1;
            let k2 = k1 + mult - p;
            let k3 = k2 - k1;
            for jj in k1 + 1..=k2 {
                put(pw, k1, sum_points(at(pw, k1)?, at(pw, jj)?))?;
            }
            let d = 1.0 / (k3 as f64 + 1.0);
            put(pw, k1, scale_point(d, at(pw, k1)?))?;
            k1 += 1;
            for jj in k2 + 1..=*nn {
                put(pw, k1, at(pw, jj)?)?;
                k1 += 1;
            }
            *nn -= k3;
            k1 = ii + p;
            for jj in ii + mult..=*mm {
                put(up, k1, at(up, jj)?)?;
                k1 += 1;
            }
            *mm -= k3;
            mult = p;
        }
        ii += mult;
    }

    // The first knot's multiplicity is at most p + 1.
    let mut mult = 0;
    while us == at(up, p + mult + 1)? {
        mult += 1;
    }
    if mult > 0 {
        let mut k1 = p + 1;
        for jj in p + mult + 1..=*mm {
            put(up, k1, at(up, jj)?)?;
            k1 += 1;
        }
        *mm -= mult;
        k1 = 1;
        for jj in mult + 1..=*nn {
            put(pw, k1, at(pw, jj)?)?;
            k1 += 1;
        }
        *nn -= mult;
    }

    // The last knot's multiplicity is at most p + 1.
    let mut mult = 0;
    while ue == at(up, *mm - p - mult - 1)? {
        mult += 1;
    }
    if mult > 0 {
        *mm -= mult;
        put(pw, *nn - mult, at(pw, *nn)?)?;
        *nn -= mult;
    }
    Ok(())
}
