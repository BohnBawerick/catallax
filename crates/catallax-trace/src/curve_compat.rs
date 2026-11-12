// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: NLib stacks, in-place aliasing and the legacy error
// flag become owned values and `Fault`, and curves are checked by
// `BSplineCurve::new`. Translated directly from SMLib NL_CrvTool.cpp
// (`N_CrvElevateDegree`, `N_CrvRemoveKnot`, `N_CrvsMakeCompatibleAdjKnots`),
// NL_Bezier.cpp (`N_BezGetDegreeElevationMatrix`, `N_BezElevateDegree`),
// NL_BasisBasic.cpp (`N_BasisIDerivs`, `N_BasisFindGlobalMax`,
// `N_KnotVectorIsValid`, `N_BasisGetSpanCount`, `N_GetCompatibleKnotVectorToTol`),
// NL_CrvBasic.cpp (`N_CrvsMakeRatCompatible`, `N_CrvReparamToInterval`),
// NL_CrvNurb.cpp (`N_CrvGetMinMaxWeightsAndPts`), NL_Geom.cpp
// (`N_DistCptCptHomo`), NL_Math.cpp (`N_VectorBlendCPt`, `N_FloatOpIsBad`)
// and SmBSplineCurveNLib.cpp (`SyncronizeKnotsOfCurves`). The matching
// rust-brep units (`curve/compatible_adjusted_knots.rs`, `curve/knot_removal.rs`,
// `basis/basic.rs`, `curve/sm_nurbs.rs`) were read as a cross-check.
//! Making curves compatible: degree elevation, knot removal and knot
//! synchronisation, the steps that skinned and ruled surfaces need before
//! their curves share one knot vector.
//!
//! Every result is an unproved candidate. Where SMLib returns a curve the
//! checked storage cannot hold, the port refuses with [`Fault::InvalidInput`];
//! the gate records each such case against SMLib's result. A legacy error in
//! SMLib's working copies becomes a [`Fault::Degenerate`] and leaves the
//! caller's curves unchanged, as `SyncronizeKnotsOfCurves` does.

#![allow(clippy::needless_range_loop, clippy::explicit_counter_loop)]

use crate::{
    Fault,
    curve_build::refine_parts,
    spline::{BSplineCurve, ControlPoint, NO_WEIGHT, NO_Z, combine},
};

/// `NL_DMAX`: the highest degree of an NLib curve (`NL_MAXDEG`).
const MAX_DEGREE: i64 = 17;
/// `NL_WMIN` and `NL_WMAX`: the weights a knot removal may leave.
pub(crate) const MIN_WEIGHT: f64 = 1.0e-3;
pub(crate) const MAX_WEIGHT: f64 = 100.0;
/// The knot-removal basis tolerance factor `cto` of NL_CrvTool.cpp.
pub(crate) const BASIS_TOLERANCE: f64 = 1.0e-5;

/// A curve's raw storage, without the checks of [`BSplineCurve::new`].
#[derive(Clone, Debug)]
struct Parts {
    degree: usize,
    knots: Vec<f64>,
    points: Vec<ControlPoint>,
}

impl Parts {
    fn of(curve: &BSplineCurve) -> Self {
        Self {
            degree: curve.degree(),
            knots: curve.knots().to_vec(),
            points: curve.control_points().to_vec(),
        }
    }

    fn into_curve(self) -> Result<BSplineCurve, Fault> {
        BSplineCurve::new(self.degree, self.knots, self.points)
    }

    /// `N_IsCrvRat`.
    fn is_rational(&self) -> bool {
        self.points[0].coordinates[3] != NO_WEIGHT
    }
}

const ZERO: ControlPoint = ControlPoint {
    coordinates: [0.0; 4],
};

/// `N_VectorBlendCPt`: `target += alpha * source`.
fn blend(alpha: f64, source: ControlPoint, target: &mut ControlPoint) {
    let (a, b) = (source.coordinates, &mut target.coordinates);
    b[0] += alpha * a[0];
    b[1] += alpha * a[1];
    b[2] = if a[2] != NO_Z {
        b[2] + alpha * a[2]
    } else {
        NO_Z
    };
    b[3] = if a[3] != NO_WEIGHT {
        b[3] + alpha * a[3]
    } else {
        NO_WEIGHT
    };
}

/// `N_DistCptCptHomo`.
pub(crate) fn homogeneous_distance(p: ControlPoint, q: ControlPoint) -> f64 {
    let (p, q) = (p.coordinates, q.coordinates);
    let mut discriminant = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]);
    if p[2] != NO_Z && q[2] != NO_Z {
        discriminant += (p[2] - q[2]) * (p[2] - q[2]);
    }
    if p[3] != NO_WEIGHT && q[3] != NO_WEIGHT {
        discriminant += (p[3] - q[3]) * (p[3] - q[3]);
    }
    discriminant.sqrt()
}

/// `N_FloatOpIsBad` for a division `x / y`.
pub(crate) fn division_is_bad(mut x: f64, mut y: f64) -> bool {
    if y == 0.0 {
        return true;
    }
    if x == 0.0 {
        return false;
    }
    x = x.abs();
    y = y.abs();
    if y <= 1.0 && x < f64::MAX * y {
        return false;
    }
    if y >= 1.0 && x > f64::MIN_POSITIVE * y {
        return false;
    }
    true
}

/// `NL_PascalTri[n][k]`, which SMLib tabulates up to degree 17.
fn pascal(n: i64, k: i64) -> f64 {
    let mut value = 1u64;
    for i in 0..k as u64 {
        value = value * (n as u64 - i) / (i + 1);
    }
    value as f64
}

/// `N_KnotVectorIsValid`: end multiplicities are `p + 1`, knots do not
/// decrease and no interior knot is repeated more than `p` times.
fn knot_vector_is_valid(u: &[f64], p: i64) -> bool {
    let m = u.len() as i64 - 1;
    let at = |i: i64| u[i as usize];
    let mut i = 0;
    while i < m && at(i) == at(i + 1) {
        i += 1;
    }
    if i != p {
        return false;
    }
    i = m;
    while i > 0 && at(i) == at(i - 1) {
        i -= 1;
    }
    if m - i != p {
        return false;
    }
    let (mut start, mut end) = (p, p + 1);
    while end <= m - p {
        if at(start) > at(end) {
            return false;
        }
        i = end;
        while end < m - p && at(end) == at(end + 1) {
            end += 1;
        }
        if end - i + 1 > p {
            return false;
        }
        start = end;
        end += 1;
    }
    true
}

/// `N_BasisGetSpanCount`.
fn span_count(u: &[f64], p: i64) -> i64 {
    let m = u.len() as i64 - 1;
    (p..m - p)
        .filter(|&i| u[i as usize] != u[i as usize + 1])
        .count() as i64
}

/// `N_BasisIDerivs` from the left: the value and derivatives, up to `order`,
/// of the `i`th basis function at `x`. `None` is SMLib's error flag.
fn basis_function_derivatives(u: &[f64], i: i64, p: i64, x: f64, order: i64) -> Option<[f64; 3]> {
    let m = u.len() as i64 - 1;
    let n = m - p - 1;
    let at = |j: i64| u[j as usize];
    if x < at(0) || x > at(m) || i < 0 || i > n {
        return None;
    }
    let mut result = [0.0; 3];
    let mut table = vec![vec![0.0; p as usize + 1]; p as usize + 1];
    let right_closed = x == at(m - p) && i >= n - p;
    if !right_closed && (x < at(i) || x >= at(i + p + 1)) {
        return Some(result);
    }
    for j in 0..=p {
        let inside = if right_closed {
            x > at(i + j) && x <= at(i + j + 1)
        } else {
            x >= at(i + j) && x < at(i + j + 1)
        };
        table[j as usize][0] = if inside { 1.0 } else { 0.0 };
    }
    for k in 1..=p {
        let (ku, mut saved) = (k as usize, 0.0);
        if table[0][ku - 1] != 0.0 {
            saved = ((x - at(i)) * table[0][ku - 1]) / (at(i + k) - at(i));
        }
        for j in 0..p - k + 1 {
            let ju = j as usize;
            let (right, left) = (at(i + j + k + 1), at(i + j + 1));
            if table[ju + 1][ku - 1] == 0.0 {
                table[ju][ku] = saved;
                saved = 0.0;
            } else {
                let temp = table[ju + 1][ku - 1] / (right - left);
                table[ju][ku] = saved + (right - x) * temp;
                saved = (x - left) * temp;
            }
        }
    }
    result[0] = table[0][p as usize];
    let mut derivative = [0.0; 20];
    for k in 1..=p.min(order) {
        for j in 0..=k {
            derivative[j as usize] = table[j as usize][(p - k) as usize];
        }
        for l in 1..=k {
            let mut saved = 0.0;
            if derivative[0] != 0.0 {
                saved = derivative[0] / (at(i + p - k + l) - at(i));
            }
            for j in 0..k - l + 1 {
                let ju = j as usize;
                let (right, left) = (at(p - k + l + i + j + 1), at(i + j + 1));
                let factor = (p - k + l) as f64;
                if derivative[ju + 1] == 0.0 {
                    derivative[ju] = factor * saved;
                    saved = 0.0;
                } else {
                    let temp = derivative[ju + 1] / (right - left);
                    derivative[ju] = factor * (saved - temp);
                    saved = temp;
                }
            }
        }
        result[k as usize] = derivative[0];
    }
    Some(result)
}

/// `N_BasisFindGlobalMax`: the largest value of the `i`th basis function and
/// where it is taken. `None` is SMLib's error flag, which `N_CrvRemoveKnot`
/// answers with 1.
pub(crate) fn basis_global_max(u: &[f64], i: i64, p: i64, tolerance: f64) -> Option<(f64, f64)> {
    const ITERATIONS: i64 = 20;
    const SAMPLES: i64 = 10;
    let m = u.len() as i64 - 1;
    let at = |j: i64| u[j as usize];
    if i < 0 || i > m - p - 1 {
        return None;
    }
    if p == 1 {
        return Some((1.0, at(i + 1)));
    }
    let mut repeated = 1;
    for k in i + 1..i + p {
        if at(k) == at(k + 1) {
            repeated += 1;
        }
    }
    if repeated == p {
        return Some((1.0, at(i + 1)));
    }
    let steps = p * SAMPLES;
    let mut du = (at(i + p + 1) - at(i)) / steps as f64;
    let (mut right, mut left) = (at(i), at(i));
    while right < at(i + p + 1) {
        left = right;
        right += du;
        if right > at(i + p + 1) {
            right = at(i + p + 1);
        }
        // A step too small to move the parameter would never end the scan;
        // SMLib would loop, the port reports its error flag.
        if right <= left {
            return None;
        }
        let derivatives = basis_function_derivatives(u, i, p, right, 1)?;
        if derivatives[1] < 0.0 {
            break;
        }
    }
    let mut iteration = 0;
    let mut guess = 0.0;
    let mut converged = false;
    let mut derivatives = [0.0; 3];
    while iteration < ITERATIONS {
        guess = 0.5 * (left + right);
        let mut k = 0;
        while k < ITERATIONS {
            derivatives = basis_function_derivatives(u, i, p, guess, 2)?;
            if derivatives[1].abs() < tolerance && derivatives[0] > tolerance {
                converged = true;
                break;
            }
            if division_is_bad(derivatives[1], derivatives[2]) {
                return None;
            }
            guess -= derivatives[1] / derivatives[2];
            if guess <= at(i) || guess >= at(i + p + 1) {
                break;
            }
            k += 1;
        }
        if converged {
            break;
        }
        du = (right - left) / SAMPLES as f64;
        right = left;
        for _ in 1..=SAMPLES {
            left = right;
            right += du;
            let derivatives = basis_function_derivatives(u, i, p, right, 1)?;
            if derivatives[1] < 0.0 {
                break;
            }
        }
        iteration += 1;
    }
    if iteration >= ITERATIONS {
        return None;
    }
    Some((derivatives[0], guess))
}

/// `N_BezGetDegreeElevationMatrix`: the `(p + t + 1) x (p + 1)` matrix.
fn elevation_matrix(p: i64, t: i64) -> Vec<Vec<f64>> {
    let q = p + t;
    let mut matrix = vec![vec![0.0; p as usize + 1]; q as usize + 1];
    let half = q / 2;
    matrix[0][0] = 1.0;
    matrix[q as usize][p as usize] = 1.0;
    for i in 1..=half {
        let inverse = 1.0 / pascal(q, i);
        let low = if i - t > 0 { i - t } else { 0 };
        let high = if i > p { p } else { i };
        for j in low..=high {
            matrix[i as usize][j as usize] = inverse * pascal(p, j) * pascal(t, i - j);
        }
    }
    for i in half + 1..q {
        let low = if i - t > 0 { i - t } else { 0 };
        let high = if i > p { p } else { i };
        for j in low..=high {
            matrix[i as usize][j as usize] = matrix[(q - i) as usize][(p - j) as usize];
        }
    }
    matrix
}

/// `N_CrvElevateDegree`: the curve of degree `p + t` that the same curve
/// has on the knots with each interior multiplicity raised by `t`.
fn elevate(curve: &Parts, t: i64) -> Result<Parts, Fault> {
    let p = curve.degree as i64;
    let (np, mp) = (curve.points.len() as i64 - 1, curve.knots.len() as i64 - 1);
    if t < 0 || t > MAX_DEGREE - p || !knot_vector_is_valid(&curve.knots, p) {
        return Err(Fault::Degenerate);
    }
    let up = &curve.knots;
    let pw = &curve.points;
    let s = span_count(up, p);
    let nq = np + t * s;
    let q = p + t;
    let mq = mp + t * (s + 1);
    let mut qw = vec![ZERO; nq as usize + 1];
    let mut uq = vec![0.0; mq as usize + 1];
    if t == 0 {
        qw.copy_from_slice(pw);
        uq.copy_from_slice(up);
        return Ok(Parts {
            degree: q as usize,
            knots: uq,
            points: qw,
        });
    }
    let mut bw = vec![ZERO; p as usize + 1];
    let mut nw = vec![ZERO; p as usize + 1];
    let mut dw = vec![ZERO; q as usize + 1];
    let mut alphas = vec![0.0; p as usize + 1];
    let mut complements = vec![0.0; p as usize + 1];
    let (mut a, mut b, mut r, mut cind, mut kind) = (p, p + 1, -1i64, 1i64, q + 1);
    for i in 0..=q {
        uq[i as usize] = up[a as usize];
    }
    qw[0] = pw[0];
    for i in 0..=p {
        bw[i as usize] = pw[i as usize];
    }
    let matrix = elevation_matrix(p, t);
    let upi = |i: i64| up[i as usize];
    while b < mp {
        let mut i = b;
        while b < mp && upi(b) == upi(b + 1) {
            b += 1;
        }
        let multiplicity = b - i + 1;
        let oldr = r;
        r = p - multiplicity;
        let lbz = if oldr > 0 { (oldr + 2) / 2 } else { 1 };
        let rbz = if r > 0 { q - (r + 1) / 2 } else { q };
        if r > 0 {
            let numerator = upi(b) - upi(a);
            for k in (multiplicity + 1..=p).rev() {
                let slot = (k - multiplicity - 1) as usize;
                alphas[slot] = numerator / (upi(a + k) - upi(a));
                complements[slot] = 1.0 - alphas[slot];
            }
            for j in 1..=r {
                let save = r - j;
                let s = multiplicity + j;
                for k in (s..=p).rev() {
                    bw[k as usize] = combine(
                        alphas[(k - s) as usize],
                        bw[k as usize],
                        complements[(k - s) as usize],
                        bw[(k - 1) as usize],
                    );
                }
                nw[save as usize] = bw[p as usize];
            }
        }
        // `N_BezElevateDegree` for the entries lbz..=q.
        for i in lbz..=q {
            let low = if i - t > 0 { i - t } else { 0 };
            let high = if i > p { p } else { i };
            dw[i as usize] = ZERO;
            for j in low..=high {
                blend(
                    matrix[i as usize][j as usize],
                    bw[j as usize],
                    &mut dw[i as usize],
                );
            }
        }
        if oldr > 1 {
            let mut first = kind - 2;
            let mut last = kind;
            let den = upi(b) - upi(a);
            let beta = (upi(b) - uq[(kind - 1) as usize]) / den;
            let complement = 1.0 - beta;
            for k in 1..oldr {
                i = first;
                let mut j = last;
                let mut l = j - kind + 1;
                while j - i > k {
                    if i < cind {
                        let alpha = (upi(b) - uq[i as usize]) / (upi(a) - uq[i as usize]);
                        qw[i as usize] =
                            combine(alpha, qw[i as usize], 1.0 - alpha, qw[(i - 1) as usize]);
                    }
                    if j >= lbz {
                        if j - k <= kind - q + oldr {
                            let gamma = (upi(b) - uq[(j - k) as usize]) / den;
                            dw[l as usize] =
                                combine(gamma, dw[l as usize], 1.0 - gamma, dw[(l + 1) as usize]);
                        } else {
                            dw[l as usize] =
                                combine(beta, dw[l as usize], complement, dw[(l + 1) as usize]);
                        }
                    }
                    i += 1;
                    j -= 1;
                    l -= 1;
                }
                first -= 1;
                last += 1;
            }
        }
        if a != p {
            for _ in 0..q - oldr {
                uq[kind as usize] = upi(a);
                kind += 1;
            }
        }
        for i in lbz..=rbz {
            qw[cind as usize] = dw[i as usize];
            cind += 1;
        }
        if b < mp {
            for i in 0..r {
                bw[i as usize] = nw[i as usize];
            }
            for i in r..=p {
                bw[i as usize] = pw[(b - p + i) as usize];
            }
            a = b;
            b += 1;
        } else {
            for i in 0..=q {
                uq[(kind + i) as usize] = upi(b);
            }
        }
    }
    Ok(Parts {
        degree: q as usize,
        knots: uq,
        points: qw,
    })
}

/// `N_CrvRemoveKnot`: removes up to `times` copies of the knot `u` while
/// the curve moves by no more than `tolerance`. Returns the curve and the
/// number of copies removed.
fn remove_knot(
    curve: &Parts,
    u: f64,
    times: i64,
    mut tolerance: f64,
) -> Result<(Parts, i64), Fault> {
    let p = curve.degree as i64;
    let (n, m) = (curve.points.len() as i64 - 1, curve.knots.len() as i64 - 1);
    let rational = curve.is_rational();
    if rational {
        let wmin = curve
            .points
            .iter()
            .map(|point| point.coordinates[3])
            .reduce(|a, b| if b < a { b } else { a })
            .unwrap();
        // `N_CrvGetMinMaxWeightsAndPts` always reports a largest position of 0.
        tolerance = (tolerance * wmin) / (1.0 + 0.0);
    }
    let basis_tolerance = BASIS_TOLERANCE * (curve.knots[m as usize] - curve.knots[0]).abs();
    // `N_KnotVectorIsEndParam`.
    if u <= curve.knots[0] || u >= curve.knots[m as usize] {
        return Err(Fault::Degenerate);
    }
    let span = crate::spline::find_span(&curve.knots, curve.degree, u, false) as i64;
    let multiplicity = if u == curve.knots[(m - p) as usize] {
        p + 1
    } else {
        let mut count = 0;
        let mut k = span;
        while k >= 0 {
            let equal = u == curve.knots[k as usize];
            k -= 1;
            if !equal {
                break;
            }
            count += 1;
        }
        count
    };
    let times = times.min(multiplicity);
    let mut out = curve.clone();
    let mut work = vec![ZERO; 2 * p as usize + 1];
    let fout = (2 * span - multiplicity - p) / 2;
    let mut first = span - p;
    let mut last = span - multiplicity;
    let mut removed = 0;
    // `N_CrvRemoveKnot` returns the error flag of its last call to
    // `N_BasisFindGlobalMax`, though it carries on with a maximum of 1.
    let mut flagged = false;
    let mut global_max = |knots: &[f64], i: i64| {
        let found = basis_global_max(knots, i, p, basis_tolerance);
        flagged = found.is_none();
        found.map_or(1.0, |(max, _)| max)
    };
    while removed < times {
        let k = removed;
        let offset = first - 1;
        if offset < 0 || last + 1 > n || last + 1 - offset > 2 * p {
            return Err(Fault::Degenerate);
        }
        let (mut i, mut j, mut ii, mut jj) = (first, last, 1i64, last - offset);
        let uq = |out: &Parts, index: i64| out.knots[index as usize];
        work[0] = out.points[offset as usize];
        work[(last + 1 - offset) as usize] = out.points[(last + 1) as usize];
        let (mut alpha, mut beta) = (0.5, 0.5);
        while j - i > k {
            alpha = (uq(&out, i + p + 1) - uq(&out, i)) / (u - uq(&out, i));
            let beta_numerator = uq(&out, j + p - k + 1) - uq(&out, j - k);
            beta = beta_numerator / (uq(&out, j + p - k + 1) - u);
            work[ii as usize] = combine(
                alpha,
                out.points[i as usize],
                1.0 - alpha,
                work[(ii - 1) as usize],
            );
            work[jj as usize] = combine(
                beta,
                out.points[j as usize],
                1.0 - beta,
                work[(jj + 1) as usize],
            );
            i += 1;
            j -= 1;
            ii += 1;
            jj -= 1;
        }
        let delta = (u - uq(&out, i)) / (uq(&out, i + p + 1) - uq(&out, i));
        let delta_complement = 1.0 - delta;
        if j - i < k {
            let lambda = alpha / (alpha + beta);
            let lambda_complement = 1.0 - lambda;
            let distance = homogeneous_distance(work[(ii - 1) as usize], work[(jj + 1) as usize]);
            let left = global_max(&out.knots, i - 1);
            let right = global_max(&out.knots, i);
            let (a, b) = (
                lambda * alpha * left,
                lambda_complement * delta_complement * right,
            );
            let max = if a > b { a } else { b };
            if max * distance > tolerance {
                break;
            }
            work[(jj + 1) as usize] = combine(
                lambda,
                work[(jj + 1) as usize],
                lambda_complement,
                work[(ii - 1) as usize],
            );
        } else {
            let average = combine(
                delta,
                work[(jj + 1) as usize],
                delta_complement,
                work[(ii - 1) as usize],
            );
            let distance = homogeneous_distance(out.points[i as usize], average);
            let max = global_max(&out.knots, i);
            if max * distance > tolerance {
                break;
            }
        }
        if rational {
            let (mut low, mut high) = (first, last);
            let (mut wmin, mut wmax) = (f64::MAX, f64::MIN_POSITIVE);
            while high - low > k {
                let wi = work[(low - offset) as usize].coordinates[3];
                let wj = work[(high - offset) as usize].coordinates[3];
                if wi < wmin {
                    wmin = wi;
                }
                if wj < wmin {
                    wmin = wj;
                }
                if wi > wmax {
                    wmax = wi;
                }
                if wj > wmax {
                    wmax = wj;
                }
                low += 1;
                high -= 1;
            }
            if wmin < MIN_WEIGHT || wmax > MAX_WEIGHT {
                break;
            }
        }
        let (mut low, mut high) = (first, last);
        while high - low > k {
            out.points[low as usize] = work[(low - offset) as usize];
            out.points[high as usize] = work[(high - offset) as usize];
            low += 1;
            high -= 1;
        }
        for l in span - k..=m - k - 1 {
            out.knots[l as usize] = out.knots[(l + 1) as usize];
        }
        first -= 1;
        last += 1;
        removed += 1;
    }
    if removed == 0 {
        return if flagged {
            Err(Fault::Degenerate)
        } else {
            Ok((out, 0))
        };
    }
    let mut j = fout;
    let mut i = j;
    for l in 1..removed {
        if l % 2 != 0 {
            i += 1;
        } else {
            j -= 1;
        }
    }
    for l in i + 1..=n {
        out.points[j as usize] = out.points[l as usize];
        j += 1;
    }
    out.points.truncate((n - removed + 1) as usize);
    out.knots.truncate((m - removed + 1) as usize);
    if flagged {
        return Err(Fault::Degenerate);
    }
    Ok((out, removed))
}

/// What `N_GetCompatibleKnotVectorToTol` leaves: each vector with its knots
/// moved to the shared averages, and the knots each one lacks.
#[derive(Clone, Debug, PartialEq)]
pub struct CompatibleKnots {
    /// The input vector, with knots within the tolerance averaged.
    pub adjusted: Vec<f64>,
    /// The knots to insert so the vector matches all the others.
    pub additions: Vec<f64>,
}

/// `N_GetCompatibleKnotVectorToTol`: merges the knot vectors in `vectors`,
/// treating knots within `tolerance` as one.
pub fn compatible_knot_vectors(
    vectors: &[Vec<f64>],
    tolerance: f64,
) -> Result<Vec<CompatibleKnots>, Fault> {
    // SMLib's merge makes no progress for a tolerance that is not positive
    // and finite or for knots that are not finite; refuse both.
    if !tolerance.is_finite()
        || tolerance <= 0.0
        || vectors.is_empty()
        || vectors.iter().flatten().any(|knot| !knot.is_finite())
    {
        return Err(Fault::InvalidInput);
    }
    let count = vectors.len();
    let mut knots: Vec<Vec<f64>> = vectors.to_vec();
    let top: Vec<i64> = knots.iter().map(|u| u.len() as i64 - 1).collect();
    let mut index = vec![0i64; count];
    let mut multiplicity = vec![0i64; count];
    let mut additions: Vec<Vec<f64>> = vec![Vec::new(); count];
    // `ST_notdone`.
    let not_done = |index: &[i64]| (0..count).any(|r| index[r] < top[r]);
    while not_done(&index) {
        if index[0] > top[0] {
            return Err(Fault::Degenerate);
        }
        let gap_after = |u: &[f64], at: i64, top: i64| {
            let mut next = at + 1;
            while next <= top && u[next as usize] <= u[at as usize] {
                next += 1;
            }
            if next <= top {
                u[next as usize] - u[at as usize]
            } else {
                -1.0
            }
        };
        let mut smallest = knots[0][index[0] as usize];
        let mut gap = gap_after(&knots[0], index[0], top[0]);
        for r in 1..count {
            if index[r] > top[r] {
                return Err(Fault::Degenerate);
            }
            if knots[r][index[r] as usize] < smallest {
                smallest = knots[r][index[r] as usize];
                gap = gap_after(&knots[r], index[r], top[r]);
            }
        }
        let limit = if gap < tolerance && gap != -1.0 {
            0.99999 * gap
        } else {
            tolerance
        };
        for r in 0..count {
            if (knots[r][index[r] as usize] - smallest).abs() < limit {
                let start = index[r];
                while index[r] < top[r]
                    && knots[r][index[r] as usize + 1] <= knots[r][index[r] as usize]
                {
                    index[r] += 1;
                }
                multiplicity[r] = index[r] - start + 1;
            } else {
                multiplicity[r] = 0;
                index[r] -= 1;
            }
        }
        if multiplicity.iter().all(|&m| m == 0) {
            // No vector holds the smallest knot: the scan cannot advance.
            return Err(Fault::Degenerate);
        }
        let mut contributors = 0;
        let mut sum = 0.0;
        for r in 0..count {
            if multiplicity[r] > 0 && (knots[r][index[r] as usize] - smallest).abs() < limit {
                sum += knots[r][index[r] as usize];
                contributors += 1;
            }
        }
        let average = sum / f64::from(contributors);
        for r in 0..count {
            if index[r] >= 0 && (knots[r][index[r] as usize] - average).abs() < limit {
                for s in 1..=multiplicity[r] {
                    knots[r][(index[r] - multiplicity[r] + s) as usize] = average;
                }
            }
        }
        let highest = multiplicity.iter().copied().max().unwrap_or(0);
        for r in 0..count {
            for _ in multiplicity[r] + 1..=highest {
                additions[r].push(average);
            }
            index[r] += 1;
        }
    }
    Ok(knots
        .into_iter()
        .zip(additions)
        .map(|(adjusted, additions)| CompatibleKnots {
            adjusted,
            additions,
        })
        .collect())
}

/// `N_CrvsMakeRatCompatible`: every curve rational if any is, and every
/// curve three-dimensional if any is.
fn make_rational_compatible(parts: &mut [Parts]) {
    if parts.iter().any(Parts::is_rational) {
        for curve in parts.iter_mut().filter(|curve| !curve.is_rational()) {
            for point in &mut curve.points {
                point.coordinates[3] = 1.0;
            }
        }
    }
    if parts
        .iter()
        .any(|curve| curve.points[0].coordinates[2] != NO_Z)
    {
        for curve in parts.iter_mut() {
            if curve.points[0].coordinates[2] == NO_Z {
                for point in &mut curve.points {
                    point.coordinates[2] = 0.0;
                }
            }
        }
    }
}

/// `N_CrvsMakeCompatibleAdjKnots` on copies of `curves`: the same curves
/// with one rationality and dimension, one degree, one parameter range and
/// one knot vector, knots within `tolerance` counting as equal.
fn compatible_parts(curves: &[BSplineCurve], tolerance: f64) -> Result<Vec<Parts>, Fault> {
    let mut parts: Vec<Parts> = curves.iter().map(Parts::of).collect();
    make_rational_compatible(&mut parts);
    let ends = |curve: &Parts| (curve.knots[0], curve.knots[curve.knots.len() - 1]);
    let (lower, upper) = ends(&parts[0]);
    let mut last_knot = upper;
    let mut scale = false;
    for curve in &parts[1..] {
        let (start, end) = ends(curve);
        if (start - lower).abs() > tolerance || (end - upper).abs() > tolerance {
            scale = true;
            break;
        }
        if end < last_knot {
            last_knot = end;
        }
    }
    if scale {
        for curve in &mut parts {
            rescale(curve, 0.0, 1.0);
        }
        last_knot = 1.0;
    }
    let highest = parts.iter().map(|curve| curve.degree).max().unwrap();
    for curve in &mut parts {
        let raise = (highest - curve.degree) as i64;
        if raise > 0 {
            *curve = elevate(curve, raise)?;
        }
    }
    for curve in &mut parts {
        let m = curve.knots.len() - 1;
        let knot = curve.knots[m - highest - 1];
        if last_knot - knot <= tolerance {
            if m == 2 * highest + 1 {
                return Err(Fault::Degenerate);
            }
            *curve = remove_knot(curve, knot, highest as i64, f64::MAX)?.0;
        }
    }
    let vectors: Vec<Vec<f64>> = parts.iter().map(|curve| curve.knots.clone()).collect();
    let merged = compatible_knot_vectors(&vectors, tolerance)?;
    for (curve, merge) in parts.iter_mut().zip(merged) {
        curve.knots = merge.adjusted;
        if !merge.additions.is_empty() {
            let (knots, points) =
                refine_parts(curve.degree, &curve.knots, &curve.points, &merge.additions)?;
            curve.knots = knots;
            curve.points = points;
        }
    }
    Ok(parts)
}

/// `N_CrvsMakeCompatibleKnotTol` on copies of `curves`: one rationality,
/// dimension, degree, parameter range (unit when they differ) and knot
/// vector. A `knot_tolerance` of 0 is chosen from the range as SMLib does.
fn compatible_knot_tol_parts(
    curves: &[BSplineCurve],
    mut knot_tolerance: f64,
) -> Result<Vec<Parts>, Fault> {
    let mut parts: Vec<Parts> = curves.iter().map(Parts::of).collect();
    make_rational_compatible(&mut parts);
    let ends = |curve: &Parts| (curve.knots[0], curve.knots[curve.knots.len() - 1]);
    let (start, end) = ends(&parts[0]);
    let differs = parts[1..].iter().any(|curve| {
        let (u1, u2) = ends(curve);
        (u1 - start).abs() > 1.0e-12 || (u2 - end).abs() > 1.0e-12
    });
    if differs {
        for curve in &mut parts {
            rescale(curve, 0.0, 1.0);
        }
        // All curves now run over the unit interval.
        if knot_tolerance == 0.0 {
            knot_tolerance = 0.001;
        }
    }
    if knot_tolerance == 0.0 {
        knot_tolerance = (end - start) / 1000.0;
    }
    let highest = parts.iter().map(|curve| curve.degree).max().unwrap();
    for curve in &mut parts {
        let raise = (highest - curve.degree) as i64;
        if raise > 0 {
            *curve = elevate(curve, raise)?;
        }
    }
    let vectors: Vec<Vec<f64>> = parts.iter().map(|curve| curve.knots.clone()).collect();
    // `N_GetCompatibleKnotArray` is the same merge with a floor of 1e-8.
    let tolerance = if knot_tolerance <= 1.0e-8 {
        1.0e-8
    } else {
        knot_tolerance
    };
    let merged = compatible_knot_vectors(&vectors, tolerance)?;
    for (curve, merge) in parts.iter_mut().zip(merged) {
        curve.knots = merge.adjusted;
        if !merge.additions.is_empty() {
            let (knots, points) =
                refine_parts(curve.degree, &curve.knots, &curve.points, &merge.additions)?;
            curve.knots = knots;
            curve.points = points;
        }
    }
    Ok(parts)
}

/// `N_CrvsMakeCompatible` (a knot tolerance of 0) and
/// `N_CrvsMakeCompatibleKnotTol`: copies of `curves` on one knot vector.
pub fn make_compatible(
    curves: &[BSplineCurve],
    knot_tolerance: f64,
) -> Result<Vec<BSplineCurve>, Fault> {
    if curves.is_empty() {
        return Err(Fault::InvalidInput);
    }
    compatible_knot_tol_parts(curves, knot_tolerance)?
        .into_iter()
        .map(Parts::into_curve)
        .collect()
}

/// `N_CrvReparamToInterval`.
fn rescale(curve: &mut Parts, a: f64, b: f64) {
    let p = curve.degree;
    let m = curve.knots.len() - 1;
    let u = &mut curve.knots;
    if a != u[0] || b != u[m] {
        let u0 = u[0];
        let factor = (b - a) / (u[m] - u[0]);
        for i in 0..=p {
            u[i] = a;
        }
        for i in p + 1..=m - p - 1 {
            u[i] = factor * (u[i] - u0) + a;
        }
        for i in m - p..=m {
            u[i] = b;
        }
    }
}

impl BSplineCurve {
    /// `N_CrvElevateDegree`: the same curve at degree `degree() + times`.
    /// Refuses a degree above 16 and a knot vector whose end knots are not
    /// of multiplicity `degree + 1`.
    pub fn elevate_degree(&self, times: usize) -> Result<Self, Fault> {
        // A count beyond the signed range is beyond degree 17 as well.
        let times = i64::try_from(times).map_err(|_| Fault::Degenerate)?;
        elevate(&Parts::of(self), times)?.into_curve()
    }

    /// `N_CrvRemoveKnot`: removes up to `times` copies of the interior knot
    /// `parameter` while the curve moves by no more than `tolerance`, and
    /// reports how many were removed.
    pub fn remove_knot(
        &self,
        parameter: f64,
        times: usize,
        tolerance: f64,
    ) -> Result<(Self, usize), Fault> {
        // Any count above the knot's multiplicity removes all of them.
        let times = i64::try_from(times).unwrap_or(i64::MAX);
        let (parts, removed) = remove_knot(&Parts::of(self), parameter, times, tolerance)?;
        Ok((parts.into_curve()?, removed as usize))
    }

    /// `N_CrvReparamToInterval`: the same curve with its knots rescaled to
    /// run from `start` to `end`.
    pub fn rescaled(&self, start: f64, end: f64) -> Result<Self, Fault> {
        let mut parts = Parts::of(self);
        rescale(&mut parts, start, end);
        parts.into_curve()
    }
}

/// `SmBSplineCurve::SyncronizeKnotsOfCurves`, through
/// `N_CrvsMakeCompatibleAdjKnots`: copies of `curves` on one knot vector,
/// one degree, one rationality and one dimension. On an error the input
/// curves stand unchanged, as SMLib works on copies.
pub fn synchronize_knots(
    curves: &[BSplineCurve],
    tolerance: f64,
) -> Result<Vec<BSplineCurve>, Fault> {
    if curves.is_empty() {
        return Err(Fault::InvalidInput);
    }
    compatible_parts(curves, tolerance)?
        .into_iter()
        .map(Parts::into_curve)
        .collect()
}
