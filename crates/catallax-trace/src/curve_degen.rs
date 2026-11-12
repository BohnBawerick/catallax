// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline curves of `spline.rs`;
// NLib stacks and in-place aliasing become owned vectors, statuses become
// `Fault`, and NLib's error record is dropped (these routines only clear it
// or fail). Adapted from rust-brep crates/kernel/src/curve/arc_length.rs
// (`N_CrvArcLength`, `N_CrvGetDegenSegs`, `N_BezCalcLength`,
// `ST_BezCalcArcLength`, `N_BezSplit`, `N_DistCPolygon`),
// curve/geometry.rs (`N_CrvEvalCurvature`, `N_CrvEvalFrenetFrame`) and
// surface/srf_type.rs (`N_CrvGetType`), commit
// 50404914a96c3855cb21cd314d1a3b870fcb0310, translations of SMLib
// NL_CrvGeom.cpp, NL_CrvNurb.cpp, NL_Bezier.cpp and NL_Geom.cpp;
// `N_CrvRemoveDegenSegs` (NL_CrvTool.cpp), `N_CrvEval`,
// `N_CrvGetBBoxMaxDiagDist` (NL_CrvNurb.cpp) and `N_ProjectPtPlane`,
// `N_IsectLinePlane`, `N_VectorsAngle` (NL_Geom.cpp) are translated
// directly at c9c5979f04cea8c38213004d4f04df0dcfb40d05.
// Indexed loops follow NLib's index arithmetic.
#![allow(clippy::needless_range_loop, clippy::explicit_counter_loop)]
//! The short spans a trim can leave: NLib's search for curve segments no
//! longer than a tolerance, the arc length it measures them with, and their
//! removal, which `SmBSplineCurve::Trim` applies to a cut curve.
//!
//! Every result is an unproved candidate.

use crate::Fault;
use crate::curve_build::span_and_multiplicity;
use crate::curve_compat::division_is_bad;
use crate::curve_edit::TrimRefusal;
use crate::curve_nlib::{curve_derivatives, extract_curve_segment};
use crate::scalar::PI;
use crate::spline::{BSplineCurve, ControlPoint, NO_WEIGHT, NO_Z, combine, find_span};
use crate::surface_eval::{MAX_DEGREE, basis_eval};
use crate::vector::Vector3;

/// `NL_PTOL`, the parameter-space tolerance.
const PARAMETER_TOLERANCE: f64 = 1.0e-12;
/// `NL_LTOL`, the parallel-line tolerance.
const LINE_TOLERANCE: f64 = 1.0e-9;
/// `NL_RECURSELIM`, the recursion limit of `ST_BezCalcArcLength`.
const RECURSION_LIMIT: usize = 20;
/// `ctl` of `N_CrvArcLength`: the line and circle test tolerance relative
/// to the segment's box diagonal.
const CURVE_TYPE_TOLERANCE: f64 = 1.0e-9;

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

/// `N_DistPtPt`.
fn distance(p: Vector3, q: Vector3) -> f64 {
    (p - q).length()
}

/// `N_VectorNormalize`: the unit vector and the length.
fn normalize(a: Vector3) -> Result<(Vector3, f64), Fault> {
    let d = a.length();
    if division_is_bad(1.0, d) {
        return Err(Fault::Degenerate);
    }
    let inverse = 1.0 / d;
    Ok((a * inverse, d))
}

/// `N_CrvEval` from `NL_LEFT`: the curve point at `u`.
fn curve_point(curve: &BSplineCurve, u: f64) -> Result<Vector3, Fault> {
    let (p, knots, points) = (curve.degree(), curve.knots(), curve.control_points());
    // `basis_eval` holds `NL_MAXDEG` degrees, while `BSplineCurve::new`
    // accepts more; a higher degree would index past its work arrays.
    if p > MAX_DEGREE {
        return Err(Fault::InvalidInput);
    }
    // `N_KnotVectorIsParamOutOfBounds`.
    if u < knots[0] || u > knots[knots.len() - 1] {
        return Err(Fault::InvalidInput);
    }
    let (span, basis) = basis_eval(knots, p, u, true)?;
    let mut cw = [0.0; 4];
    let j = span - p;
    for i in j..=span {
        let alpha = basis[i - j];
        let [x, y, z, w] = points[i].coordinates;
        cw[0] += alpha * x;
        cw[1] += alpha * y;
        cw[2] = if z != NO_Z { cw[2] + alpha * z } else { NO_Z };
        cw[3] = if w != NO_WEIGHT {
            cw[3] + alpha * w
        } else {
            NO_WEIGHT
        };
    }
    Ok(euclid(ControlPoint { coordinates: cw }))
}

/// `N_CrvGetBBoxMaxDiagDist`: the diagonal of the box of the Euclidean
/// control points.
fn box_diagonal(curve: &BSplineCurve) -> f64 {
    let points = curve.control_points();
    let first = euclid(points[0]);
    let (mut low, mut high) = ([first.x, first.y, first.z], [first.x, first.y, first.z]);
    for point in &points[1..] {
        let e = euclid(*point);
        for (k, value) in [e.x, e.y, e.z].into_iter().enumerate() {
            if value < low[k] {
                low[k] = value;
            }
            if value > high[k] {
                high[k] = value;
            }
        }
    }
    distance(
        Vector3::new(low[0], low[1], low[2]),
        Vector3::new(high[0], high[1], high[2]),
    )
}

/// `N_BasisGetSpanCount`: the knot spans of nonzero length.
fn span_count(knots: &[f64], p: usize) -> usize {
    (p..knots.len() - 1 - p)
        .filter(|&i| knots[i] != knots[i + 1])
        .count()
}

/// `N_DistCPolygon`: the length of the control polygon.
fn polygon_length(points: &[ControlPoint]) -> f64 {
    let mut length = 0.0;
    for i in 1..points.len() {
        length += distance(euclid(points[i - 1]), euclid(points[i]));
    }
    length
}

/// `N_BezSplit` at one half, the only split the arc length takes.
fn split_in_half(points: &[ControlPoint]) -> (Vec<ControlPoint>, Vec<ControlPoint>) {
    let u = 0.5;
    let p = points.len() - 1;
    let mut work = points.to_vec();
    let mut left = points.to_vec();
    let mut right = points.to_vec();
    left[0] = points[0];
    right[p] = points[p];
    for k in 1..=p {
        for i in 0..=p - k {
            work[i] = combine(1.0 - u, work[i], u, work[i + 1]);
        }
        left[k] = work[0];
        right[p - k] = work[p - k];
    }
    (left, right)
}

/// `ST_BezCalcArcLength`: `eps` is the tolerance to beat on entry and the
/// tolerance reached on exit; `length` accumulates the sums of the polygon
/// and chord lengths. Past the recursion limit it fails.
fn bezier_arc_length(
    points: &[ControlPoint],
    eps: &mut f64,
    length: &mut f64,
    depth: usize,
) -> Result<(), Fault> {
    if depth > RECURSION_LIMIT {
        return Err(Fault::Degenerate);
    }
    let p = points.len() - 1;
    let upper = polygon_length(points);
    let lower = distance(euclid(points[0]), euclid(points[p]));
    let del = (upper - lower).abs();
    if del <= *eps {
        *length += upper + lower;
        *eps = del;
    } else {
        let (left, right) = split_in_half(points);
        let mut eps_left = *eps / 2.0;
        bezier_arc_length(&left, &mut eps_left, length, depth + 1)?;
        let mut eps_right = *eps - eps_left;
        bezier_arc_length(&right, &mut eps_right, length, depth + 1)?;
        *eps = eps_left + eps_right;
    }
    Ok(())
}

/// `N_BezCalcLength`: the arc length of a Bezier curve to `tol`.
fn bezier_length(points: &[ControlPoint], tol: f64) -> Result<f64, Fault> {
    let mut eps = 2.0 * tol;
    let mut total = 0.0;
    bezier_arc_length(points, &mut eps, &mut total, 0)?;
    Ok(0.5 * total)
}

/// The shape `N_CrvGetType` reports.
#[derive(Clone, Copy, Debug, PartialEq)]
enum CurveType {
    /// A line from `start` along `direction` (the chord, not unit).
    Line { start: Vector3, direction: Vector3 },
    /// A circular arc of `radius`, sweeping `angle` degrees.
    Circle { radius: f64, angle: f64 },
    /// Neither.
    FreeForm,
}

/// `N_CrvEvalCurvature` with its osculating circle, from `NL_LEFT`: the
/// curvature, the point, the circle's centre and the binormal (SMLib's
/// circle-plane normal).
fn curvature_circle(
    curve: &BSplineCurve,
    t: f64,
) -> Result<(f64, Vector3, Vector3, Vector3), Fault> {
    let d = curve_derivatives(curve, t, true, 2)?;
    let num = d[1].cross(d[2]).length();
    let mut den = d[1].length();
    den = den * den * den;
    if division_is_bad(num, den) {
        return Err(Fault::Degenerate);
    }
    let kappa = num / den;
    if division_is_bad(1.0, kappa) {
        return Err(Fault::Degenerate);
    }
    let radius = 1.0 / kappa;
    // `N_CrvEvalFrenetFrame`; its normal and binormal outputs are passed
    // swapped, so the binormal is returned as the circle's normal.
    if curve.degree() <= 1 {
        return Err(Fault::Degenerate);
    }
    let cd = curve_derivatives(curve, t, true, 2)?;
    let (tangent, _) = normalize(cd[1])?;
    let (binormal, _) = normalize(cd[1].cross(cd[2]))?;
    let normal = binormal.cross(tangent);
    let center = d[0] + normal * radius;
    Ok((kappa, cd[0], center, binormal))
}

/// `N_ProjectPtPlane` through `N_IsectLinePlane`: `point` moved along the
/// plane normal onto the plane.
fn project_point_plane(origin: Vector3, normal: Vector3, point: Vector3) -> Result<Vector3, Fault> {
    let (v, length_v) = normalize(normal)?;
    let (n, _) = normalize(normal)?;
    let pl = origin - point;
    let den = v.dot(n);
    if den.abs() < LINE_TOLERANCE {
        // SMLib leaves the output unwritten. The line runs along the
        // normal, so `den` is one to rounding and this is not reached.
        return Err(Fault::Degenerate);
    }
    let num = pl.dot(n);
    if division_is_bad(num, den) {
        return Err(Fault::Degenerate);
    }
    let alpha = num / den;
    if division_is_bad(alpha, length_v) {
        return Err(Fault::Degenerate);
    }
    Ok(point + v * alpha)
}

/// `N_VectorsAngle`: the angle between two vectors in degrees.
fn vectors_angle(a: Vector3, b: Vector3) -> Result<f64, Fault> {
    let dot = a.dot(b);
    let (ma, mb) = (a.length(), b.length());
    if division_is_bad(dot, ma * mb) {
        return Err(Fault::Degenerate);
    }
    // Like SMLib's two comparisons, `clamp` keeps a NaN quotient.
    let arg = (dot / (ma * mb)).clamp(-1.0, 1.0);
    let rad = crate::scalar::acos(arg);
    Ok(180.0 * (rad / PI))
}

/// `N_CrvGetType`: whether the curve is a line or a circular arc within
/// `tol`. A failed evaluation leaves it free-form, as SMLib clears the error.
fn curve_type(curve: &BSplineCurve, tol: f64) -> Result<CurveType, Fault> {
    let (p, knots, points) = (curve.degree(), curve.knots(), curve.control_points());
    let n = points.len() - 1;
    let m = knots.len() - 1;
    let a = euclid(points[0]);
    let b = euclid(points[n]);
    let chord = b - a;
    let dab = chord.length();
    let mut line = true;
    // SMLib tests `NL_MTOL` and `NL_PTOL`; the smaller one decides.
    if dab <= PARAMETER_TOLERANCE {
        // Closed curves are not lines.
        line = false;
    } else if n > 1 {
        let dab2 = dab * dab;
        let mut kh = 2 * p;
        let mut tl = knots[p];
        let mut i = p + 1;
        'spans: while i < m {
            while i < m && knots[i] == knots[i + 1] {
                i += 1;
            }
            let tr = knots[i];
            let ti = (tr - tl) / kh as f64;
            if i >= m {
                kh -= 1;
            }
            for k in 1..=kh {
                let t = if k == kh && i < m {
                    tr
                } else {
                    tl + k as f64 * ti
                };
                let Ok(q) = curve_point(curve, t) else {
                    return Ok(CurveType::FreeForm);
                };
                let mut d = (q - a).dot(chord);
                d /= dab2;
                if !(0.0..=1.0).contains(&d) {
                    line = false;
                    break 'spans;
                }
                let r = a + chord * d;
                if distance(q, r) > tol || distance(q, a) > dab {
                    line = false;
                    break 'spans;
                }
            }
            tl = tr;
            i += 1;
        }
    }
    if line {
        return Ok(CurveType::Line {
            start: a,
            direction: chord,
        });
    }

    let closed = distance(a, b) < tol;
    let kh = 2 * p;
    let mut j: i64 = -1;
    let mut tl = knots[p];
    let mut i = p + 1;
    let tol2 = 2.0 * tol;
    let mut max_radius = 0.0;
    let mut min_radius = f64::MAX;
    let mut center_sum = Vector3::new(0.0, 0.0, 0.0);
    let mut normal_sum = Vector3::new(0.0, 0.0, 0.0);
    let mut samples = Vec::new();
    while i < m {
        let kl = if i == p + 1 { 0 } else { 1 };
        while i < m && knots[i] == knots[i + 1] {
            i += 1;
        }
        let kr = if closed && i == m { kh - 1 } else { kh };
        let tr = knots[i];
        let ti = (tr - tl) / kh as f64;
        for k in kl..=kr {
            let t = if k == kh { tr } else { tl + k as f64 * ti };
            j += 1;
            let Ok((kappa, point, center, normal)) = curvature_circle(curve, t) else {
                return Ok(CurveType::FreeForm);
            };
            if division_is_bad(1.0, kappa) {
                return Ok(CurveType::FreeForm);
            }
            let rd = 1.0 / kappa;
            if min_radius > rd {
                min_radius = rd;
            }
            if max_radius < rd {
                max_radius = rd;
            }
            if max_radius - min_radius > tol2 {
                return Ok(CurveType::FreeForm);
            }
            center_sum = center_sum + center;
            normal_sum = normal_sum + normal;
            samples.push(point);
        }
        tl = tr;
        i += 1;
    }
    let fac = 1.0 / (j as f64 + 1.0);
    let center = center_sum * fac;
    let normal = normal_sum * fac;
    let mut radius_sum = 0.0;
    let mut max_radius = 0.0;
    let mut min_radius = f64::MAX;
    for &sample in &samples {
        let q = project_point_plane(center, normal, sample)?;
        let d = distance(sample, center);
        radius_sum += d;
        if min_radius > d {
            min_radius = d;
        }
        if max_radius < d {
            max_radius = d;
        }
        if distance(sample, q) > tol {
            return Ok(CurveType::FreeForm);
        }
    }
    let radius = fac * radius_sum;
    if (max_radius - radius).abs() > tol || (radius - min_radius).abs() > tol {
        return Ok(CurveType::FreeForm);
    }
    let (x_axis, _) = normalize(a - center)?;
    let y_axis = normal.cross(x_axis);
    let angle = if closed {
        360.0
    } else {
        let to_end = b - center;
        let fac = vectors_angle(x_axis, to_end)?;
        if y_axis.dot(to_end) < 0.0 {
            360.0 - fac
        } else {
            fac
        }
    };
    Ok(CurveType::Circle { radius, angle })
}

/// `N_CrvArcLength` with an absolute tolerance: the length of the curve
/// over `[ul, ur]` to `tol`.
fn arc_length(curve: &BSplineCurve, ul: f64, ur: f64, tol: f64) -> Result<f64, Fault> {
    if (ur - ul).abs() < PARAMETER_TOLERANCE {
        return Ok(0.0);
    }
    let segment = extract_curve_segment(curve, ul, ur)?;
    let lctl = box_diagonal(&segment) * CURVE_TYPE_TOLERANCE;
    match curve_type(&segment, lctl)? {
        CurveType::Line { start, direction } => return Ok(distance(start, start + direction)),
        CurveType::Circle { radius, angle } => return Ok((angle / 180.0) * PI * radius),
        CurveType::FreeForm => {}
    }
    let (p, knots, points) = (segment.degree(), segment.knots(), segment.control_points());
    let m = knots.len() - 1;
    let eps = tol / span_count(knots, p) as f64;
    let mut length = 0.0;
    let mut bezier = points[..=p].to_vec();
    let mut next = bezier.clone();
    let mut alphas = vec![0.0; p + 1];
    let mut complements = vec![0.0; p + 1];
    let (mut a, mut b) = (p, p + 1);
    let mut r: i64;
    while b < m {
        let i = b;
        while b < m && knots[b] == knots[b + 1] {
            b += 1;
        }
        let multiplicity = b - i + 1;
        r = p as i64 - multiplicity as i64;
        if r > 0 {
            let num = knots[b] - knots[a];
            for k in (multiplicity + 1..=p).rev() {
                alphas[k - multiplicity - 1] = num / (knots[a + k] - knots[a]);
                complements[k - multiplicity - 1] = 1.0 - alphas[k - multiplicity - 1];
            }
            for j in 1..=r as usize {
                let save = r as usize - j;
                let s = multiplicity + j;
                for k in (s..=p).rev() {
                    bezier[k] =
                        combine(alphas[k - s], bezier[k], complements[k - s], bezier[k - 1]);
                }
                next[save] = bezier[p];
            }
        }
        length += bezier_length(&bezier, eps)?;
        if b < m {
            if r < 0 {
                return Err(Fault::Degenerate);
            }
            let r = r as usize;
            bezier[..r].copy_from_slice(&next[..r]);
            for i in r..=p {
                bezier[i] = points[b - p + i];
            }
            a = b;
            b += 1;
        }
    }
    Ok(length)
}

/// `N_CrvGetDegenSegs`: the parameter intervals of the curve's segments no
/// longer than `tol`, neighbours merged while their summed length stays
/// within `tol`.
fn degenerate_segments(curve: &BSplineCurve, tol: f64) -> Result<Vec<[f64; 2]>, Fault> {
    let (p, knots) = (curve.degree(), curve.knots());
    let m = knots.len() - 1;
    let mut segments: Vec<[f64; 2]> = Vec::new();
    let mut length = 0.0;
    let mut ii = p;
    let mut start = curve_point(curve, knots[ii])?;
    loop {
        let end = curve_point(curve, knots[ii + 1])?;
        let mut d1 = distance(start, end);
        if d1 <= tol {
            if m <= 2 * p + 1 {
                let middle = curve_point(curve, 0.5 * (knots[ii] + knots[ii + 1]))?;
                let d2 = distance(start, middle);
                d1 = distance(middle, end);
                if d2 + d1 > tol {
                    break;
                }
            }
            d1 = arc_length(curve, knots[ii], knots[ii + 1], tol)?;
            if d1 > tol {
                length = 0.0;
            } else if length > 0.0 && length + d1 <= tol {
                if let Some(last) = segments.last_mut() {
                    last[1] = knots[ii + 1];
                }
                length += d1;
            } else {
                length = d1;
                segments.push([knots[ii], knots[ii + 1]]);
            }
        } else {
            length = 0.0;
        }
        if ii as i64 >= m as i64 - p as i64 - 1 {
            break;
        }
        start = end;
        ii += 1;
        while knots[ii] == knots[ii + 1] {
            ii += 1;
            // SMLib reads past the knots of a curve whose last span has no
            // length; refuse it instead.
            if ii + 1 > m {
                return Err(Fault::InvalidInput);
            }
        }
    }
    Ok(segments)
}

/// `N_CrvInsertKnot` on the working knots and control points.
fn insert_knot(
    p: usize,
    knots: &mut Vec<f64>,
    points: &mut Vec<ControlPoint>,
    u: f64,
    times: usize,
) -> Result<(), Fault> {
    let curve = BSplineCurve::new(p, knots.clone(), points.clone())?.insert_knot(u, times)?;
    *knots = curve.knots().to_vec();
    *points = curve.control_points().to_vec();
    Ok(())
}

/// The last index of the run of degenerate segments that starts at `ii`.
fn run_end(segments: &[[f64; 2]], ii: usize) -> usize {
    let ns = segments.len() - 1;
    let mut jj = ii;
    while jj < ns {
        if segments[jj][1] != segments[jj + 1][0] {
            break;
        }
        jj += 1;
    }
    jj
}

/// `N_BasisFindSpan` with `NL_LEFT` (the span that starts at `u`).
fn span_left(knots: &[f64], p: usize, u: f64) -> usize {
    find_span(knots, p, u, false)
}

/// `N_BasisFindSpan` with `NL_RIGHT` (the span that ends at `u`).
fn span_right(knots: &[f64], p: usize, u: f64) -> usize {
    find_span(knots, p, u, true)
}

/// `N_CrvRemoveDegenSegs`: the curve without its segments no longer than
/// `tol` (arc length), or `None` when it has none to remove. Each run of
/// such segments is first given a parameter length of at most `0.00005` of
/// the curve's, then cut out. Runs covering the whole curve need
/// `N_CrvRemoveKnots`, which is not ported: [`TrimRefusal::ShortSpan`].
pub(crate) fn remove_degenerate_segments(
    curve: &BSplineCurve,
    tol: f64,
) -> Result<Option<BSplineCurve>, TrimRefusal> {
    let (p, up, pw) = (curve.degree(), curve.knots(), curve.control_points());
    let (np, mp) = (pw.len() - 1, up.len() - 1);
    if np <= p {
        return Ok(None);
    }
    let mut segs = degenerate_segments(curve, tol)?;
    if segs.is_empty() {
        return Ok(None);
    }
    let ns = segs.len() - 1;
    if ns == 0 && segs[0][0] == up[0] && segs[0][1] == up[mp] {
        return Ok(None);
    }
    let dlen = 0.00005 * (up[mp] - up[0]);
    let mut uq = up.to_vec();
    let mut qw = pw.to_vec();

    // Make each run a Bezier piece and shrink its parameter length.
    let mut ii = 0;
    while ii <= ns {
        let jj = run_end(&segs, ii);
        if segs[ii][0] == up[0] && segs[jj][1] == up[mp] {
            // `N_CrvRemoveKnots` on the whole curve.
            return Err(TrimRefusal::ShortSpan);
        }
        let (_, mult1) = span_and_multiplicity(&uq, p, segs[ii][0]);
        if mult1 < p as i64 {
            insert_knot(p, &mut uq, &mut qw, segs[ii][0], p - mult1 as usize)?;
        }
        let (_, mult2) = span_and_multiplicity(&uq, p, segs[jj][1]);
        if mult2 < p as i64 {
            insert_knot(p, &mut uq, &mut qw, segs[jj][1], p - mult2 as usize)?;
        }
        let mq = uq.len() - 1;
        let du = segs[jj][1] - segs[ii][0];
        if du > dlen {
            if segs[ii][0] == uq[0] {
                // A run at the start.
                let span2 = span_right(&uq, p, segs[jj][1]) + 1;
                let u2 = uq[0] + dlen;
                let d2 = dlen / du;
                for kk in ii..jj {
                    segs[kk][1] = d2 * (segs[kk][1] - uq[0]) + uq[0];
                    segs[kk + 1][0] = d2 * (segs[kk + 1][0] - uq[0]) + uq[0];
                }
                for kk in p + 1..span2 {
                    uq[kk] = d2 * (uq[kk] - uq[0]) + uq[0];
                }
                shrink_right(&mut segs, &mut uq, p, jj, span2, u2);
            } else if segs[jj][1] == uq[mq] {
                // A run at the end.
                let span1 = span_left(&uq, p, segs[ii][0]);
                let u2 = uq[mq] - dlen;
                let d2 = dlen / du;
                for kk in ii..jj {
                    segs[kk][1] = d2 * (segs[kk][1] - uq[span1]) + u2;
                    segs[kk + 1][0] = d2 * (segs[kk + 1][0] - uq[span1]) + u2;
                }
                for kk in span1 + 1..mq - p {
                    uq[kk] = d2 * (uq[kk] - uq[span1]) + u2;
                }
                shrink_left(&mut segs, &mut uq, p, ii, span1, u2);
            } else {
                // An interior run.
                let span1 = span_left(&uq, p, segs[ii][0]);
                let span2 = span_right(&uq, p, segs[jj][1]) + 1;
                let d1 = 0.5 * (segs[ii][0] + segs[jj][1]);
                let u1 = d1 - 0.499 * dlen;
                let u2 = d1 + 0.499 * dlen;
                let d2 = (u2 - u1) / du;
                for kk in ii..jj {
                    segs[kk][1] = d2 * (segs[kk][1] - uq[span1]) + u1;
                    segs[kk + 1][0] = d2 * (segs[kk + 1][0] - uq[span1]) + u1;
                }
                for kk in span1 + 1..span2 {
                    uq[kk] = d2 * (uq[kk] - uq[span1]) + u1;
                }
                shrink_left(&mut segs, &mut uq, p, ii, span1, u1);
                shrink_right(&mut segs, &mut uq, p, jj, span2, u2);
            }
        }
        ii = jj + 1;
    }

    // Cut the runs out: end runs by dropping their knots and control
    // points, interior runs by collapsing them to one knot of multiplicity
    // `p` at their middle.
    let mut ii = 0;
    while ii <= ns {
        let jj = run_end(&segs, ii);
        let (nq, mq) = (qw.len() - 1, uq.len() - 1);
        if segs[ii][0] == up[0] {
            let span2 = span_right(&uq, p, segs[jj][1]);
            let mut k1 = p + 1;
            for kk in span2 + p + 1..=mq {
                uq[k1] = uq[kk];
                k1 += 1;
            }
            uq.truncate(k1);
            let mut k1 = 1;
            for kk in span2 + 1..=nq {
                qw[k1] = qw[kk];
                k1 += 1;
            }
            qw.truncate(k1);
        } else if segs[jj][1] == up[mp] {
            let span1 = span_left(&uq, p, segs[ii][0]) + 1;
            let mut k1 = span1 - p;
            for _ in 0..=p {
                uq[k1] = uq[mq];
                k1 += 1;
            }
            uq.truncate(span1 + 1);
            let last = nq + span1 - mq;
            qw[last] = qw[nq];
            qw.truncate(last + 1);
        } else {
            let u2 = 0.5 * (segs[ii][0] + segs[jj][1]);
            let (_, mult2) = span_and_multiplicity(&uq, p, u2);
            if mult2 < p as i64 {
                insert_knot(p, &mut uq, &mut qw, u2, p - mult2 as usize)?;
            }
            let (nq, mq) = (qw.len() - 1, uq.len() - 1);
            let span1 = span_left(&uq, p, segs[ii][0]);
            let span2 = span_left(&uq, p, u2);
            let span3 = span_left(&uq, p, segs[jj][1]);
            let mut k1 = span1 - p + 1;
            for _ in 0..p {
                uq[k1] = u2;
                k1 += 1;
            }
            for kk in span3 + 1..=mq {
                uq[k1] = uq[kk];
                k1 += 1;
            }
            uq.truncate(k1);
            let mut k1 = span1 - p;
            qw[k1] = qw[span2 - p];
            k1 += 1;
            for kk in span3 - p + 1..=nq {
                qw[k1] = qw[kk];
                k1 += 1;
            }
            qw.truncate(k1);
        }
        ii = jj + 1;
    }
    Ok(Some(BSplineCurve::new(p, uq, qw)?))
}

/// The knots and later runs right of a shrunk run ending at `segs[jj][1]`
/// (the knot run from `span2`), rescaled from `u2` to the curve's end.
fn shrink_right(segs: &mut [[f64; 2]], uq: &mut [f64], p: usize, jj: usize, span2: usize, u2: f64) {
    let ns = segs.len() - 1;
    let mq = uq.len() - 1;
    let d3 = (uq[mq] - u2) / (uq[mq] - segs[jj][1]);
    for kk in jj + 1..=ns {
        segs[kk][0] = d3 * (segs[kk][0] - segs[jj][1]) + u2;
        if segs[kk][1] != uq[mq] {
            segs[kk][1] = d3 * (segs[kk][1] - segs[jj][1]) + u2;
        }
    }
    for kk in span2 + p..mq - p {
        uq[kk] = d3 * (uq[kk] - uq[span2]) + u2;
    }
    segs[jj][1] = u2;
    for kk in 0..p {
        uq[span2 + kk] = u2;
    }
}

/// The knots and earlier runs left of a shrunk run starting at
/// `segs[ii][0]` (the knot run ending at `span1`), rescaled from the curve's
/// start to `u1`.
fn shrink_left(segs: &mut [[f64; 2]], uq: &mut [f64], p: usize, ii: usize, span1: usize, u1: f64) {
    let d3 = (u1 - uq[0]) / (uq[span1] - uq[0]);
    for kk in 0..ii {
        segs[kk][1] = d3 * (segs[kk][1] - uq[0]) + uq[0];
        if segs[kk][0] != uq[0] {
            segs[kk][0] = d3 * (segs[kk][0] - uq[0]) + uq[0];
        }
    }
    for kk in p + 1..=span1.saturating_sub(p) {
        uq[kk] = d3 * (uq[kk] - uq[0]) + uq[0];
    }
    segs[ii][0] = u1;
    for kk in 0..p {
        uq[span1 - kk] = u1;
    }
}
