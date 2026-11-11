// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline surface candidates of
// `spline.rs` in their non-polynomial form, with SMLib's operation order
// kept; statuses become `Fault`. The power-basis form is not carried, since
// no Catallax candidate uses it; the out-of-bounds switch is a separate entry
// point.
// Adapted from rust-brep crates/kernel/src/surface/evaluation.rs,
// surface/sm_surface.rs and basis/basic.rs, commit
// 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab, translations of SMLib
// NL_SrfNurb.cpp (N_SrfEvalPt, N_SrfBasisDerivs, N_SrfDerivs),
// NL_BasisBasic.cpp (N_BasisEval) and SmBSplineSurfaceNLib.cpp
// (SmBSplineSurface::Evaluate, EvaluatePoint); EvaluateSimple is translated
// directly.
//! Point and derivative evaluation of B-spline surfaces as SMLib does it,
//! including its repair of vanishing tangents.

use crate::Fault;
use crate::extent::Extent2;
use crate::scalar::{EFF_ZERO, ZONE_TOL_3D};
use crate::spline::{BSplineSurface, NO_WEIGHT, NO_Z, find_span};
use crate::vector::{Vector2, Vector3};

/// `NL_MAXDER`: the highest derivative order evaluated.
const MAX_DERIVATIVE: usize = 4;

/// `NL_MAXDEG`: the highest degree NLib's basis work arrays hold.
pub const MAX_DEGREE: usize = 17;

/// The surface forms NLib evaluates: degrees up to [`MAX_DEGREE`], and knot
/// vectors clamped at both ends (degree + 1 equal end knots), which NLib's
/// end-value shortcut in `N_BasisEval` assumes. SMLib's surfaces have these
/// forms; others the constructor accepts are refused with
/// `Fault::InvalidInput` rather than indexed past the work arrays or
/// evaluated with the shortcut where it does not hold.
pub fn check_supported(surface: &BSplineSurface) -> Result<(), Fault> {
    for (degree, knots) in surface.degrees().into_iter().zip(surface.knots()) {
        let n = knots.len();
        if degree > MAX_DEGREE
            || n < 2 * (degree + 1)
            || knots[..=degree].iter().any(|&k| k != knots[0])
            || knots[n - 1 - degree..].iter().any(|&k| k != knots[n - 1])
        {
            return Err(Fault::InvalidInput);
        }
    }
    Ok(())
}

/// The control point `[i][j]` (u index, v index) as stored, weighted.
fn net(surface: &BSplineSurface, i: usize, j: usize) -> [f64; 4] {
    surface.control_points()[i * surface.counts()[1] + j].coordinates
}

/// `N_BasisFindSpan`; `nl_left` is SMLib's `NL_LEFT` side.
fn span(knots: &[f64], degree: usize, u: f64, nl_left: bool) -> usize {
    find_span(knots, degree, u, !nl_left)
}

/// `N_BasisEval`: the nonzero basis functions at `u` and their span, exact
/// at the ends of the active knot range.
fn basis_eval(knots: &[f64], p: usize, u: f64, nl_left: bool) -> Result<(usize, [f64; 18]), Fault> {
    let mut output = [0.0; 18];
    let last = knots.len() - 1;
    if u == knots[p] {
        output[0] = 1.0;
        return Ok((p, output));
    }
    if u == knots[last - p] {
        output[p] = 1.0;
        return Ok((last - p - 1, output));
    }
    let i = span(knots, p, u, nl_left);
    let mut left = [0.0; 19];
    let mut right = [0.0; 19];
    output[0] = 1.0;
    for k in 1..=p {
        left[k] = u - knots[i + 1 - k];
        right[k] = knots[i + k] - u;
        let mut saved = 0.0;
        for l in 0..k {
            let delta = right[l + 1] + left[k - l];
            if delta < 1e-24 {
                return Err(Fault::Degenerate);
            }
            let temp = output[l];
            output[l] = saved + temp * (right[l + 1] / delta);
            saved = temp * (left[k - l] / delta);
        }
        output[k] = saved;
    }
    Ok((i, output))
}

fn blend(no_z: bool, no_w: bool, alpha: f64, pt: [f64; 4], output: &mut [f64; 4]) {
    output[0] += alpha * pt[0];
    output[1] += alpha * pt[1];
    if no_z {
        output[2] = NO_Z;
    } else {
        output[2] += alpha * pt[2];
    }
    if no_w {
        output[3] = NO_WEIGHT;
    } else {
        output[3] += alpha * pt[3];
    }
}

fn euclid(no_z: bool, no_w: bool, pt: [f64; 4]) -> [f64; 3] {
    if no_w {
        return [pt[0], pt[1], if no_z { 0.0 } else { pt[2] }];
    }
    let inv = 1.0 / pt[3];
    [
        pt[0] * inv,
        pt[1] * inv,
        if no_z { 0.0 } else { pt[2] * inv },
    ]
}

fn markers(surface: &BSplineSurface) -> (bool, bool) {
    let first = surface.control_points()[0].coordinates;
    (first[2] == NO_Z, !surface.is_rational())
}

/// `N_SrfEvalPt`: the surface point, rows in v first and then in u.
pub fn eval_point(
    surface: &BSplineSurface,
    uv: Vector2,
    sides: [bool; 2],
) -> Result<Vector3, Fault> {
    check_supported(surface)?;
    let [p, q] = surface.degrees();
    let [ku, kv] = surface.knots();
    let (usp, nu) = basis_eval(ku, p, uv.x, sides[0])?;
    let (vsp, nv) = basis_eval(kv, q, uv.y, sides[1])?;
    let (no_z, no_w) = markers(surface);
    let mut sum = [0.0; 4];
    for (i, &basis_u) in nu[..=p].iter().enumerate() {
        let mut row = [0.0; 4];
        for (j, &basis_v) in nv[..=q].iter().enumerate() {
            blend(
                no_z,
                no_w,
                basis_v,
                net(surface, usp - p + i, vsp - q + j),
                &mut row,
            );
        }
        blend(no_z, no_w, basis_u, row, &mut sum);
    }
    let [x, y, z] = euclid(no_z, no_w, sum);
    Ok(Vector3::new(x, y, z))
}

/// The univariate basis derivatives of `N_SrfBasisDerivs`, scaled by the
/// falling factorials of the degree.
fn univariate(knots: &[f64], p: usize, u: f64, span: usize, order: usize) -> [[f64; 18]; 5] {
    let mut ndu = [[0.0; 18]; 18];
    let mut left = [0.0; 18];
    let mut right = [0.0; 18];
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    let mut du = [[0.0; 18]; 5];
    for j in 0..=p {
        du[0][j] = ndu[j][p];
    }
    let mut work = [[0.0; 18]; 2];
    for r in 0..=p {
        let (mut a, mut b) = (0, 1);
        work[a][0] = 1.0;
        for k in 1..=p.min(order) {
            let mut value = 0.0;
            let rk = r as isize - k as isize;
            let pk = p - k;
            if r >= k {
                work[b][0] = work[a][0] / ndu[pk + 1][rk as usize];
                value = work[b][0] * ndu[rk as usize][pk];
            }
            let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
            let j2 = if r as isize - 1 <= pk as isize {
                k - 1
            } else {
                p - r
            };
            for j in j1..=j2 {
                let row = (rk + j as isize) as usize;
                work[b][j] = (work[a][j] - work[a][j - 1]) / ndu[pk + 1][row];
                value += work[b][j] * ndu[row][pk];
            }
            if r <= pk {
                work[b][k] = -work[a][k - 1] / ndu[pk + 1][r];
                value += work[b][k] * ndu[r][pk];
            }
            du[k][r] = value;
            std::mem::swap(&mut a, &mut b);
        }
    }
    let mut factor = p;
    for (k, derivative) in du.iter_mut().enumerate().take(p.min(order) + 1).skip(1) {
        for value in &mut derivative[..=p] {
            *value *= factor as f64;
        }
        factor *= p - k;
    }
    du
}

const PASCAL: [[f64; 5]; 5] = [
    [1.0, 0.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 2.0, 1.0, 0.0, 0.0],
    [1.0, 3.0, 3.0, 1.0, 0.0],
    [1.0, 4.0, 6.0, 4.0, 1.0],
];

/// `NL_IsRecipBad`: whether `1/value` would overflow or underflow.
fn reciprocal_is_bad(value: f64) -> bool {
    let y = value.abs();
    !((y <= 1.0 && 1.0 < f64::MAX * y) || (y >= 1.0 && 1.0 > f64::MIN_POSITIVE * y))
}

/// `N_SrfDerivs`: derivatives `S^(k,l)` for `k <= hu`, `l <= hv` (or
/// `k + l <= hu` with `upper`), in cells `k * (hv + 1) + l`, through the
/// rational basis derivatives as SMLib forms them. Cells not computed stay
/// zero.
pub fn derivatives(
    surface: &BSplineSurface,
    uv: Vector2,
    sides: [bool; 2],
    upper: bool,
    hu: usize,
    hv: usize,
) -> Result<Vec<Vector3>, Fault> {
    if hu > MAX_DERIVATIVE || hv > MAX_DERIVATIVE {
        return Err(Fault::Degenerate);
    }
    check_supported(surface)?;
    let [p, q] = surface.degrees();
    let [knots_u, knots_v] = surface.knots();
    let columns = |k: usize| {
        if upper && hu == hv {
            hu - k + 1
        } else {
            hv + 1
        }
    };
    let ku = span(knots_u, p, uv.x, sides[0]);
    let du = univariate(knots_u, p, uv.x, ku, hu);
    let kv = span(knots_v, q, uv.y, sides[1]);
    let dv = univariate(knots_v, q, uv.y, kv, hv);
    let index =
        |k: usize, l: usize, i: usize, j: usize| ((k * (hv + 1) + l) * (p + 1) + i) * (q + 1) + j;
    let mut basis = vec![0.0; (hu + 1) * (hv + 1) * (p + 1) * (q + 1)];
    if surface.is_rational() {
        let mut denominator = [[0.0; 5]; 5];
        let mut temporary = [0.0; 18];
        for l in 0..=hv {
            for (i, item) in temporary[..=p].iter_mut().enumerate() {
                *item = 0.0;
                for (j, &basis_v) in dv[l][..=q].iter().enumerate() {
                    *item += net(surface, ku - p + i, kv - q + j)[3] * basis_v;
                }
            }
            for k in 0..=hu {
                denominator[k][l] = 0.0;
                for i in 0..=p {
                    denominator[k][l] += temporary[i] * du[k][i];
                }
            }
        }
        if reciprocal_is_bad(denominator[0][0]) {
            return Err(Fault::Degenerate);
        }
        let inv = 1.0 / denominator[0][0];
        for n in 0..=p {
            for m in 0..=q {
                let w = net(surface, ku - p + n, kv - q + m)[3];
                for k in 0..=hu {
                    for l in 0..columns(k) {
                        let mut value = w * du[k][n] * dv[l][m];
                        for j in 1..=l {
                            value -=
                                PASCAL[l][j] * denominator[0][j] * basis[index(k, l - j, n, m)];
                        }
                        for i in 1..=k {
                            value -=
                                PASCAL[k][i] * denominator[i][0] * basis[index(k - i, l, n, m)];
                            let mut value2 = 0.0;
                            for j in 1..=l {
                                value2 += PASCAL[l][j]
                                    * denominator[i][j]
                                    * basis[index(k - i, l - j, n, m)];
                            }
                            value -= PASCAL[k][i] * value2;
                        }
                        basis[index(k, l, n, m)] = value * inv;
                    }
                }
            }
        }
    } else {
        for k in 0..=hu {
            for l in 0..columns(k) {
                for i in 0..=p {
                    for j in 0..=q {
                        basis[index(k, l, i, j)] = du[k][i] * dv[l][j];
                    }
                }
            }
        }
    }
    let (no_z, no_w) = markers(surface);
    let mut points = vec![[0.0; 3]; (p + 1) * (q + 1)];
    for i in 0..=p {
        for j in 0..=q {
            points[i * (q + 1) + j] = euclid(no_z, no_w, net(surface, ku - p + i, kv - q + j));
        }
    }
    let mut out = vec![Vector3::new(0.0, 0.0, 0.0); (hu + 1) * (hv + 1)];
    for k in 0..=hu {
        for l in 0..columns(k) {
            let mut sum = [0.0; 3];
            for i in 0..=p {
                for j in 0..=q {
                    let pt = points[i * (q + 1) + j];
                    let alpha = basis[index(k, l, i, j)];
                    sum[0] += alpha * pt[0];
                    sum[1] += alpha * pt[1];
                    if !no_z {
                        sum[2] += alpha * pt[2];
                    }
                }
            }
            out[k * (hv + 1) + l] = Vector3::new(sum[0], sum[1], sum[2]);
        }
    }
    Ok(out)
}

/// `SmSurface::GetDerivativeLimit`: below this a derivative counts as zero.
pub fn derivative_limit(position: Vector3) -> f64 {
    10.0 * EFF_ZERO * (1.0 + position.max_dimension())
}

/// `CheckZeroSurfaceDerivs`: replaces a vanishing first derivative by a
/// short vector along a second derivative, oriented into the domain.
#[allow(clippy::too_many_arguments)]
fn repair_zero_derivatives(
    uv: Vector2,
    domain: Extent2,
    su: &mut Vector3,
    sv: &mut Vector3,
    suu: Vector3,
    suv: Vector3,
    svv: Vector3,
    scaled_zero: f64,
) -> bool {
    let mut result = true;
    let new_length = 1.1 * scaled_zero;
    let mut flip = false;
    if su.length() < scaled_zero {
        result = false;
        let length = suv.length();
        if length > scaled_zero {
            *su = suv * (new_length / length);
            flip = uv.y > domain.v_interval().mid();
            result = true;
        } else {
            let length = suu.length();
            if length > scaled_zero {
                *su = suu * (new_length / length);
                flip = uv.x > domain.u_interval().mid();
                result = true;
            }
        }
        if result && flip {
            *su *= -1.0;
        }
    }
    if sv.length() < scaled_zero {
        result = false;
        let length = suv.length();
        if length > scaled_zero {
            *sv = suv * (new_length / length);
            flip = uv.x > domain.u_interval().mid();
            result = true;
        } else {
            let length = svv.length();
            if length > scaled_zero {
                *sv = svv * (new_length / length);
                flip = uv.y > domain.v_interval().mid();
                result = true;
            }
        }
        if result && flip {
            *sv *= -1.0;
        }
    }
    result
}

/// `SmBSplineSurface::Evaluate`: derivatives at `uv` clamped into the
/// natural domain, with `non_zero_tangents` replacing vanishing first
/// derivatives as SMLib does. `sides` are SMLib's `bUFromLeft` and
/// `bVFromLeft`.
#[allow(clippy::too_many_arguments)]
pub fn evaluate(
    surface: &BSplineSurface,
    uv_in: Vector2,
    hu: usize,
    hv: usize,
    sides: [bool; 2],
    upper: bool,
    non_zero_tangents: bool,
    zero_sampling: bool,
) -> Result<Vec<Vector3>, Fault> {
    evaluate_bounded(
        surface,
        uv_in,
        [hu, hv],
        sides,
        upper,
        non_zero_tangents,
        zero_sampling,
        true,
    )
}

/// [`evaluate`] with SMLib's out-of-bounds switch on
/// (`SetOutOfBoundsEnabled(TRUE)`): parameters outside the natural domain
/// extrapolate the end spans instead of being clamped.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_out_of_bounds(
    surface: &BSplineSurface,
    uv: Vector2,
    hu: usize,
    hv: usize,
    sides: [bool; 2],
    upper: bool,
    non_zero_tangents: bool,
    zero_sampling: bool,
) -> Result<Vec<Vector3>, Fault> {
    evaluate_bounded(
        surface,
        uv,
        [hu, hv],
        sides,
        upper,
        non_zero_tangents,
        zero_sampling,
        false,
    )
}

#[allow(clippy::too_many_arguments)]
fn evaluate_bounded(
    surface: &BSplineSurface,
    uv_in: Vector2,
    [hu, hv]: [usize; 2],
    sides: [bool; 2],
    upper: bool,
    non_zero_tangents: bool,
    zero_sampling: bool,
    clamp: bool,
) -> Result<Vec<Vector3>, Fault> {
    if let Some(frame) = surface.plane_frame() {
        return plane_evaluate(frame, uv_in, [hu, hv], upper, clamp);
    }
    if (upper && hu != hv) || hu > MAX_DERIVATIVE || hv > MAX_DERIVATIVE {
        return Err(Fault::InvalidInput);
    }
    let domain = surface.natural_domain();
    let uv = if !clamp || domain.contains_point(uv_in, 0.0) {
        uv_in
    } else {
        domain.clamp_point(uv_in)
    };
    let mut d = derivatives(surface, uv, sides, upper, hu, hv)?;
    if !non_zero_tangents || (hu < 1 && hv < 1) {
        return Ok(d);
    }
    let original_du = if hu >= 1 { d[hv + 1] } else { d[1] };
    let original_dv = if hv >= 1 { d[1] } else { d[hv + 1] };
    let (du_length, dv_length) = (original_du.length(), original_dv.length());
    let limit = derivative_limit(d[0]);
    let (zero_du, zero_dv) = (du_length < limit, dv_length < limit);
    if !zero_du && !zero_dv {
        return Ok(d);
    }
    if hu < 2 || hv < 2 {
        let again = evaluate_bounded(surface, uv_in, [2, 2], sides, true, true, true, clamp)?;
        if hu >= 1 {
            d[hv + 1] = again[3];
        }
        if hv >= 1 {
            d[1] = again[1];
        }
        return Ok(d);
    }
    let (mut su, mut sv) = (d[hv + 1], d[1]);
    let repaired = repair_zero_derivatives(
        uv,
        domain,
        &mut su,
        &mut sv,
        d[2 * hv + 2],
        d[hv + 2],
        d[2],
        limit,
    );
    d[hv + 1] = su;
    d[1] = sv;
    if repaired {
        return Ok(d);
    }
    let step_u = if uv.x > domain.u_interval().mid() {
        -1.0
    } else {
        1.0
    };
    let step_v = if uv.y > domain.v_interval().mid() {
        -1.0
    } else {
        1.0
    };
    let mut step = ZONE_TOL_3D / 10.0;
    while zero_sampling && step < 0.1 {
        step *= 10.0;
        let first = if du_length < dv_length {
            Vector2::new(uv.x, uv.y + step_v * step)
        } else {
            Vector2::new(uv.x + step_u * step, uv.y)
        };
        let samples = [
            first,
            Vector2::new(uv.x + step, uv.y + step),
            Vector2::new(uv.x - step, uv.y - step),
            Vector2::new(uv.x - step, uv.y + step),
            Vector2::new(uv.x + step, uv.y - step),
        ];
        for sample in samples {
            let sample = domain.clamp_point(sample);
            let Ok(mut again) = evaluate(surface, sample, 1, 1, sides, true, false, false) else {
                continue;
            };
            if again[2].length() >= limit && again[1].length() >= limit {
                if zero_du {
                    let _ = again[2].unitize();
                    again[2] *= 1.1 * limit;
                    d[hv + 1] = again[2];
                }
                if zero_dv {
                    let _ = again[1].unitize();
                    again[1] *= 1.1 * limit;
                    d[1] = again[1];
                }
                return Ok(d);
            }
        }
    }
    Err(Fault::Degenerate)
}

/// `SmBSplineSurface::EvaluatePoint`: the point at `uv`, unclamped.
pub fn evaluate_point(surface: &BSplineSurface, uv: Vector2) -> Result<Vector3, Fault> {
    if let Some(frame) = surface.plane_frame() {
        // `SmPlane::EvaluatePoint` scales the parameter before the axis.
        let uv = if frame.domain.contains_point(uv, 0.0) {
            uv
        } else {
            frame.domain.clamp_point(uv)
        };
        return Ok(frame.origin
            + uv.x * frame.uv_scale.x * frame.x_axis
            + uv.y * frame.uv_scale.y * frame.y_axis);
    }
    eval_point(surface, uv, [true, true])
}

/// `SmPlane::Evaluate`: the point and its two constant first derivatives in
/// SMLib's cell layout, every other cell zero.
fn plane_evaluate(
    frame: &crate::spline::PlaneFrame,
    uv: Vector2,
    [hu, hv]: [usize; 2],
    upper: bool,
    clamp: bool,
) -> Result<Vec<Vector3>, Fault> {
    if upper && hu > hv {
        return Err(Fault::InvalidInput);
    }
    let uv = if clamp && !frame.domain.contains_point(uv, 0.0) {
        frame.domain.clamp_point(uv)
    } else {
        uv
    };
    // SMLib zeroes the cells it writes.
    let mut cells = vec![Vector3::new(0.0, 0.0, 0.0); (hu + 1) * (hv + 1)];
    cells[0] = frame.point_fast(uv);
    if hv > 0 {
        cells[1] = frame.y_axis * frame.uv_scale.y;
    }
    if hu > 0 {
        cells[hv + 1] = frame.x_axis * frame.uv_scale.x;
    }
    Ok(cells)
}

/// `SmBSplineSurface::EvaluateSimple`: derivatives at `uv` clamped into the
/// natural domain, without tangent repair.
pub fn evaluate_simple(
    surface: &BSplineSurface,
    uv: Vector2,
    hu: usize,
    hv: usize,
    sides: [bool; 2],
    upper: bool,
) -> Result<Vec<Vector3>, Fault> {
    if let Some(frame) = surface.plane_frame() {
        // `SmPlane::EvaluateSimple` is its `Evaluate`.
        return plane_evaluate(frame, uv, [hu, hv], upper, true);
    }
    if hu > 17 || hv > 17 {
        return Err(Fault::InvalidInput);
    }
    let domain = surface.natural_domain();
    let uv = if domain.contains_point(uv, 0.0) {
        uv
    } else {
        domain.clamp_point(uv)
    };
    derivatives(surface, uv, sides, upper, hu, hv)
}

/// `SmSurface::GetNormalLimit`: below this a cross product counts as zero.
const NORMAL_LIMIT: f64 = 1.0e-8;

/// `SmSurface::Evaluate1stDerivatives`: position, `Su` and `Sv`, with
/// tangent repair.
pub fn first_derivatives(
    surface: &BSplineSurface,
    uv: Vector2,
    sides: [bool; 2],
) -> Result<[Vector3; 3], Fault> {
    let d = evaluate(surface, uv, 1, 1, sides, true, true, true)?;
    Ok([d[0], d[2], d[1]])
}

/// Unitizes a derivative only when it is longer than the limit.
fn shrink(vector: &mut Vector3, limit: f64) {
    let length = vector.length();
    if length > limit {
        *vector /= length;
    }
}

/// `SmSurface::EvaluateNormal`: the unit cross product of the first
/// derivatives, with stepping, higher-derivative and sampling fallbacks. It
/// always answers; an undefined normal is the zero vector.
pub fn evaluate_normal(surface: &BSplineSurface, uv: Vector2, sides: [bool; 2]) -> Vector3 {
    if let Some(frame) = surface.plane_frame() {
        return frame.z_axis();
    }
    let mut point = Vector3::default();
    let mut su = Vector3::default();
    let mut sv = Vector3::default();
    let (mut len_su, mut len_sv);
    let first = first_derivatives(surface, uv, sides);
    if let Ok(values) = first {
        [point, su, sv] = values;
    }
    let limit = derivative_limit(point);
    if first.is_ok() {
        let mut normal = su.cross(sv);
        let length = normal.length();
        if length > NORMAL_LIMIT {
            normal /= length;
            return normal;
        }
        len_su = su.length();
        len_sv = sv.length();
        if len_su > limit {
            su /= len_su;
            len_su = 1.0;
        }
        if len_sv > limit {
            sv /= len_sv;
            len_sv = 1.0;
        }
        if len_su > limit && len_sv > limit {
            let mut normal = su.cross(sv);
            let length = normal.length();
            if length > NORMAL_LIMIT {
                normal /= length;
                return normal;
            }
        }
    }
    let domain = surface.natural_domain();
    let mut step = 1.0e-8;
    len_su = su.length();
    len_sv = sv.length();
    if len_su >= limit && len_sv >= limit {
        let mut step_uv = Vector2::new(0.0, 0.0);
        let size = domain.size();
        if uv.x < domain.min.x + size.x / 20.0 {
            step_uv.x = size.x;
        } else if uv.x > domain.max.x - size.x / 20.0 {
            step_uv.x = -size.x;
        }
        if uv.y < domain.min.y + size.y / 20.0 {
            step_uv.y = size.y;
        } else if uv.y > domain.max.y - size.y / 20.0 {
            step_uv.y = -size.y;
        }
        if step_uv.length_squared() < EFF_ZERO {
            let mid = Vector2::new(domain.u_interval().mid(), domain.v_interval().mid());
            step_uv.x = if uv.x < mid.x { size.x } else { -size.x };
            step_uv.y = if uv.y < mid.y { size.y } else { -size.y };
        }
        while step < 0.1 {
            let Ok([_, a, b]) = first_derivatives(surface, uv + step * step_uv, sides) else {
                step *= 10.0;
                continue;
            };
            (su, sv) = (a, b);
            len_su = su.length();
            len_sv = sv.length();
            shrink(&mut su, limit);
            shrink(&mut sv, limit);
            let mut normal = su.cross(sv);
            let length = normal.length();
            if length > NORMAL_LIMIT {
                normal /= length;
                let Ok([_, mut a, mut b]) =
                    first_derivatives(surface, uv + (2.0 * step) * step_uv, sides)
                else {
                    return normal;
                };
                let first = normal;
                shrink(&mut a, limit);
                shrink(&mut b, limit);
                let mut second = a.cross(b);
                let length = second.length();
                if length > NORMAL_LIMIT {
                    second /= length;
                    normal = 2.0 * first - second;
                    let _ = normal.unitize();
                }
                return normal;
            }
            step *= 10.0;
        }
    }
    if let Ok(matrix) = evaluate(surface, uv, 2, 2, sides, true, false, true) {
        if len_su < 0.9 {
            su = matrix[3];
        }
        if len_sv < 0.9 {
            sv = matrix[1];
        }
        let (mut suu, mut suv, mut svv) = (matrix[6], matrix[4], matrix[2]);
        let flip_u = if uv.x > domain.u_interval().mid() {
            -1.0
        } else {
            1.0
        };
        let flip_v = if uv.y > domain.v_interval().mid() {
            -1.0
        } else {
            1.0
        };
        let tolerance = NORMAL_LIMIT * 1000.0;
        shrink(&mut suv, limit);
        for (mut candidate, flip) in [(suv.cross(sv), flip_v), (su.cross(suv), flip_u)] {
            let length = candidate.length();
            if length > tolerance {
                candidate *= flip / length;
                return candidate;
            }
        }
        shrink(&mut suu, limit);
        let mut candidate = suu.cross(sv);
        let length = candidate.length();
        if length > tolerance {
            candidate *= flip_u / length;
            return candidate;
        }
        shrink(&mut svv, limit);
        for (mut candidate, flip) in [
            (su.cross(svv), flip_v),
            (suu.cross(suv), flip_u * flip_u),
            (suv.cross(svv), flip_v * flip_v),
            (suu.cross(svv), flip_u * flip_v),
        ] {
            let length = candidate.length();
            if length > tolerance {
                candidate *= flip / length;
                return candidate;
            }
        }
    }
    step = 1.0e-8;
    let (u_size, v_size) = (domain.u_interval().length(), domain.v_interval().length());
    while step < 0.1 {
        let (u_step, v_step) = (step * u_size, step * v_size);
        let first = Vector2::new(
            if sides[0] {
                uv.x + u_step
            } else {
                uv.x - u_step
            },
            if sides[1] {
                uv.y + v_step
            } else {
                uv.y - v_step
            },
        );
        for sample in [
            first,
            Vector2::new(uv.x + u_step, uv.y + v_step),
            Vector2::new(uv.x - u_step, uv.y - v_step),
            Vector2::new(uv.x - u_step, uv.y + v_step),
            Vector2::new(uv.x + u_step, uv.y - v_step),
        ] {
            let Ok([_, mut a, mut b]) = first_derivatives(surface, sample, sides) else {
                continue;
            };
            shrink(&mut a, limit);
            shrink(&mut b, limit);
            let mut normal = a.cross(b);
            if normal.length() > NORMAL_LIMIT {
                let _ = normal.unitize();
                return normal;
            }
        }
        step *= 10.0;
    }
    Vector3::new(0.0, 0.0, 0.0)
}
