// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: NLib's stack-allocated surfaces become owned
// storage, its error flags become `Fault`, and where a failing NLib routine
// leaves an in-place surface half edited the port returns that surface
// beside the flag. Translated directly from SMLib NL_SrfTool.cpp
// (`N_SrfRemoveKnotConditional`), NL_SrfAdv.cpp (`N_SrfExtendByParamDist`
// for `NL_G1R`, `N_SrfEvalPtDerivsUnbounded`) and NL_SrfNurb.cpp
// (`N_SrfIsClosedSmooth`); rust-brep has no port.
//! NLib's surface knot removal and surface extension.
//!
//! Extension is the one that `SmOffsetSurface` asks for, `NL_G1R`: a ruled
//! strip with tangent-plane continuity, one span of the surface's degree, with
//! an attempt to remove one knot at the old end. The other continuities are
//! not carried.
//!
//! The routines keep NLib's quirks, which the gates pin: a removal that fails
//! on a later column or row leaves the earlier ones already edited, and a
//! failed basis maximum is reported although the removal went on with a
//! maximum of 1.
//!
//! Everything returned is an unproved candidate.

use crate::Fault;
use crate::curve_compat::division_is_bad;
use crate::curve_compat::{
    BASIS_TOLERANCE, MAX_WEIGHT, MIN_WEIGHT, basis_global_max, homogeneous_distance,
};
use crate::curve_nlib::MODEL_TOLERANCE;
use crate::scalar::{PI, acos};
use crate::spline::{BSplineSurface, ControlPoint, NO_WEIGHT, NO_Z, combine};
use crate::surface_eval::{derivatives, eval_point};
use crate::vector::{Vector2, Vector3};

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

/// A surface's raw storage, without the checks of [`BSplineSurface::new`].
#[derive(Clone, Debug)]
pub(crate) struct SurfaceParts {
    pub(crate) degrees: [usize; 2],
    pub(crate) knots: [Vec<f64>; 2],
    pub(crate) counts: [usize; 2],
    pub(crate) points: Vec<ControlPoint>,
}

impl SurfaceParts {
    pub(crate) fn of(surface: &BSplineSurface) -> Self {
        Self {
            degrees: surface.degrees(),
            knots: [0, 1].map(|axis| surface.knots()[axis].to_vec()),
            counts: surface.counts(),
            points: surface.control_points().to_vec(),
        }
    }

    pub(crate) fn into_surface(self) -> Result<BSplineSurface, Fault> {
        BSplineSurface::new(self.degrees, self.knots, self.counts, self.points)
    }

    /// `N_IsSrfRat`.
    fn is_rational(&self) -> bool {
        self.points[0].coordinates[3] != NO_WEIGHT
    }

    fn at(&self, u: i64, v: i64) -> ControlPoint {
        self.points[u as usize * self.counts[1] + v as usize]
    }

    fn set(&mut self, u: i64, v: i64, point: ControlPoint) {
        let columns = self.counts[1];
        self.points[u as usize * columns + v as usize] = point;
    }
}

/// What `N_SrfRemoveKnotConditional` leaves.
#[derive(Clone, Debug)]
pub struct SurfaceKnotRemoval {
    /// The surface after the removal. When NLib fails part way, the net is
    /// as far edited as NLib got.
    pub surface: BSplineSurface,
    /// How many copies of the knot were removed.
    pub removed: usize,
    /// NLib's error flag, set when its last call for a basis maximum failed
    /// and the removal went on with a maximum of 1.
    pub flagged: bool,
}

/// `N_SrfRemoveKnotConditional`: removes up to `times` copies of the interior
/// knot `parameter` in direction `axis` (0 for U, 1 for V) while the surface
/// moves by no more than `tolerance`. The net is edited in place, as NLib's
/// is when its two surfaces are the same: a copy whose removal fails on a
/// later row or column keeps the earlier rows and columns edited.
pub fn remove_knot_conditional(
    surface: &BSplineSurface,
    axis: usize,
    parameter: f64,
    times: usize,
    tolerance: f64,
) -> Result<SurfaceKnotRemoval, Fault> {
    let mut parts = SurfaceParts::of(surface);
    let (removed, flagged) = remove_knot_in_place(&mut parts, axis, parameter, times, tolerance)?;
    Ok(SurfaceKnotRemoval {
        surface: parts.into_surface()?,
        removed,
        flagged,
    })
}

/// The removal itself on raw storage; `(removed, flagged)`.
pub(crate) fn remove_knot_in_place(
    s: &mut SurfaceParts,
    axis: usize,
    t: f64,
    times: usize,
    mut tol: f64,
) -> Result<(usize, bool), Fault> {
    if axis > 1 {
        return Err(Fault::InvalidInput);
    }
    let (p, q) = (s.degrees[0] as i64, s.degrees[1] as i64);
    let (n, m) = (s.counts[0] as i64 - 1, s.counts[1] as i64 - 1);
    let (r, sn) = (s.knots[0].len() as i64 - 1, s.knots[1].len() as i64 - 1);
    let rational = s.is_rational();
    if rational {
        let wmin = s
            .points
            .iter()
            .map(|point| point.coordinates[3])
            .reduce(|a, b| if b < a { b } else { a })
            .unwrap_or(1.0);
        let pmax = s
            .points
            .iter()
            .map(|point| euclid(*point).length())
            .reduce(|a, b| if b > a { b } else { a })
            .unwrap_or(0.0);
        tol = (tol * wmin) / (1.0 + pmax);
    }
    let (degree, last_knot) = if axis == 0 { (p, r) } else { (q, sn) };
    let knots = &s.knots[axis];
    // `N_KnotVectorIsEndParam`.
    if t <= knots[0] || t >= knots[last_knot as usize] {
        return Err(Fault::Degenerate);
    }
    // `N_BasisFindSpanAndMult` with `NL_LEFT`.
    let span = crate::spline::find_span(knots, degree as usize, t, false) as i64;
    let multiplicity = if t == knots[(last_knot - degree) as usize] {
        degree + 1
    } else {
        let mut count = 0;
        let mut k = span;
        while k >= 0 {
            let equal = t == knots[k as usize];
            k -= 1;
            if !equal {
                break;
            }
            count += 1;
        }
        count
    };
    let times = i64::try_from(times).unwrap_or(i64::MAX).min(multiplicity);
    let a = p.max(q);
    let size = 2 * a;
    let zero = ControlPoint {
        coordinates: [0.0; 4],
    };
    let mut alf = vec![0.0; size as usize + 1];
    let mut oma = vec![0.0; size as usize + 1];
    let mut bet = vec![0.0; size as usize + 1];
    let mut omb = vec![0.0; size as usize + 1];
    let mut work = vec![zero; size as usize + 1];
    let along = if axis == 0 { n } else { m };
    let across = if axis == 0 { m } else { n };
    let fout = (2 * span - multiplicity - degree) / 2;
    let mut first = span - degree;
    let mut last = span - multiplicity;
    let mut lambda = 0.0;
    let mut lambda_complement = 0.0;
    let mut removable = true;
    let mut flagged = false;
    let basis_tolerance =
        BASIS_TOLERANCE * (s.knots[axis][last_knot as usize] - s.knots[axis][0]).abs();
    let mut k = 0;
    while k < times {
        let off = first - 1;
        // NLib's scratch arrays run to the index `2 * max(p, q)`.
        if off < 0 || last + 1 > along || last + 1 - off > size || last - first > size {
            return Err(Fault::Degenerate);
        }
        let knot = |s: &SurfaceParts, index: i64| s.knots[axis][index as usize];
        let (mut i, mut j) = (first, last);
        while j - i > k {
            let (fi, fj) = ((i - first) as usize, (j - first) as usize);
            alf[fi] = (knot(s, i + degree + 1) - knot(s, i)) / (t - knot(s, i));
            oma[fi] = 1.0 - alf[fi];
            bet[fj] =
                (knot(s, j + degree - k + 1) - knot(s, j - k)) / (knot(s, j + degree - k + 1) - t);
            omb[fj] = 1.0 - bet[fj];
            i += 1;
            j -= 1;
        }
        let delta = (t - knot(s, i)) / (knot(s, i + degree + 1) - knot(s, i));
        let delta_complement = 1.0 - delta;
        let max;
        if j - i < k {
            let (left, right) = (i - first - 1, j - first + 1);
            if left >= 0 && right >= 0 {
                if left as usize >= alf.len() || right as usize >= bet.len() {
                    return Err(Fault::Degenerate);
                }
                lambda = alf[left as usize] / (alf[left as usize] + bet[right as usize]);
                lambda_complement = 1.0 - lambda;
                let mut global_max = |index: i64| {
                    let found = basis_global_max(&s.knots[axis], index, degree, basis_tolerance);
                    flagged = found.is_none();
                    found.map_or(1.0, |(max, _)| max)
                };
                let max_left = global_max(i - 1);
                let max_right = global_max(i);
                let (x, y) = (
                    lambda * alf[left as usize] * max_left,
                    lambda_complement * delta_complement * max_right,
                );
                max = if x > y { x } else { y };
            } else {
                lambda = 0.5;
                lambda_complement = 0.5;
                max = 1.0;
            }
        } else {
            let found = basis_global_max(&s.knots[axis], i, degree, basis_tolerance);
            flagged = found.is_none();
            max = found.map_or(1.0, |(max, _)| max);
        }
        // Try to remove one knot along each row or column.
        for line in 0..=across {
            let get = |s: &SurfaceParts, index: i64| {
                if axis == 0 {
                    s.at(index, line)
                } else {
                    s.at(line, index)
                }
            };
            let (mut i, mut j) = (first, last);
            let (mut ii, mut jj) = (1i64, last - off);
            work[0] = get(s, off);
            work[(last + 1 - off) as usize] = get(s, last + 1);
            while j - i > k {
                let (fi, fj) = ((i - first) as usize, (j - first) as usize);
                work[ii as usize] = combine(alf[fi], get(s, i), oma[fi], work[(ii - 1) as usize]);
                work[jj as usize] = combine(bet[fj], get(s, j), omb[fj], work[(jj + 1) as usize]);
                i += 1;
                j -= 1;
                ii += 1;
                jj -= 1;
            }
            if j - i < k {
                let distance =
                    homogeneous_distance(work[(ii - 1) as usize], work[(jj + 1) as usize]);
                if max * distance > tol {
                    removable = false;
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
                let distance = homogeneous_distance(get(s, i), average);
                if max * distance > tol {
                    removable = false;
                    break;
                }
            }
            if rational {
                let (mut low, mut high) = (first, last);
                let (mut wmin, mut wmax) = (f64::MAX, f64::MIN_POSITIVE);
                while high - low > k {
                    let wi = work[(low - off) as usize].coordinates[3];
                    let wj = work[(high - off) as usize].coordinates[3];
                    wmin = if wi < wmin { wi } else { wmin };
                    wmin = if wj < wmin { wj } else { wmin };
                    wmax = if wi > wmax { wi } else { wmax };
                    wmax = if wj > wmax { wj } else { wmax };
                    low += 1;
                    high -= 1;
                }
                if wmin < MIN_WEIGHT || wmax > MAX_WEIGHT {
                    removable = false;
                    break;
                }
            }
            let (mut low, mut high) = (first, last);
            while high - low > k {
                let (a, b) = (work[(low - off) as usize], work[(high - off) as usize]);
                if axis == 0 {
                    s.set(low, line, a);
                    s.set(high, line, b);
                } else {
                    s.set(line, low, a);
                    s.set(line, high, b);
                }
                low += 1;
                high -= 1;
            }
        }
        if !removable {
            break;
        }
        // Successful removal: shift the knots down.
        for l in span - k..=last_knot - k - 1 {
            let next = s.knots[axis][(l + 1) as usize];
            s.knots[axis][l as usize] = next;
        }
        first -= 1;
        last += 1;
        k += 1;
    }
    if k == 0 {
        return Ok((0, flagged));
    }
    // Shift the control points down.
    let mut j = fout;
    let mut i = j;
    for l in 1..k {
        if l % 2 != 0 {
            i += 1;
        } else {
            j -= 1;
        }
    }
    for line in 0..=across {
        for (target, l) in (j..).zip(i + 1..=along) {
            let point = if axis == 0 {
                s.at(l, line)
            } else {
                s.at(line, l)
            };
            if axis == 0 {
                s.set(target, line, point);
            } else {
                s.set(line, target, point);
            }
        }
    }
    // `N_SrfSetSizeIndices`: the net and the knot vector shrink.
    let new_count = (along + 1 - k) as usize;
    let mut points = Vec::with_capacity(s.points.len());
    let counts = if axis == 0 {
        [new_count, s.counts[1]]
    } else {
        [s.counts[0], new_count]
    };
    for u in 0..counts[0] as i64 {
        for v in 0..counts[1] as i64 {
            points.push(s.at(u, v));
        }
    }
    s.points = points;
    s.counts = counts;
    s.knots[axis].truncate((last_knot - k + 1) as usize);
    Ok((k as usize, flagged))
}

/// NLib's closure flag of a surface: not closed, or smoothly closed in one or
/// both directions (`NL_NO`, `NL_UDIR`, `NL_VDIR`, `NL_UVDIR`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Closed {
    /// Closed in neither direction.
    No,
    /// Closed in U.
    U,
    /// Closed in V.
    V,
    /// Closed in both.
    UV,
}

impl Closed {
    fn in_u(self) -> bool {
        matches!(self, Self::U | Self::UV)
    }

    fn in_v(self) -> bool {
        matches!(self, Self::V | Self::UV)
    }
}

/// `N_VectorsAngle`: the angle between two vectors in degrees; `None` is
/// NLib's error flag, for a quotient that would overflow or vanish.
fn vectors_angle(a: Vector3, b: Vector3) -> Option<f64> {
    let (dot, ma, mb) = (a.dot(b), a.length(), b.length());
    if division_is_bad(dot, ma * mb) {
        return None;
    }
    let arg = (dot / (ma * mb)).clamp(-1.0, 1.0);
    Some(180.0 * (acos(arg) / PI))
}

/// `N_SrfIsClosedSmooth` with `n = 1` over both directions: whether the
/// opposite edges of the surface meet with tangent-plane continuity, within
/// `tolerances` (distance, then angle in degrees). The returned error flag
/// is NLib's, which also stands when a tangent angle could not be formed and
/// no later evaluation reset it.
pub(crate) fn is_closed_smooth(
    surface: &BSplineSurface,
    tolerances: [f64; 2],
) -> (Result<(), Fault>, Closed) {
    let [p, q] = surface.degrees();
    let [knots_u, knots_v] = surface.knots();
    let (r, s) = (knots_u.len() - 1, knots_v.len() - 1);
    let mut error: Result<(), Fault> = Ok(());
    let mut closed_u = false;
    let mut closed_v = false;
    let result = |uv: Vector2, sides: [bool; 2]| derivatives(surface, uv, sides, true, 1, 1);
    // U closure: the edges at the first and last U knots, along V.
    let mut ii = q;
    'u: loop {
        let dv = (knots_v[ii + 1] - knots_v[ii]) / q as f64;
        let mut jj = 0;
        while jj <= q {
            let (first, last) = if jj == q {
                let v = knots_v[ii + 1];
                (
                    result(Vector2::new(knots_u[0], v), [true, false]),
                    result(Vector2::new(knots_u[r], v), [false, false]),
                )
            } else {
                let v = knots_v[ii] + jj as f64 * dv;
                (
                    result(Vector2::new(knots_u[0], v), [true, true]),
                    result(Vector2::new(knots_u[r], v), [false, true]),
                )
            };
            let (d1, d2) = match (first, last) {
                (Ok(d1), Ok(d2)) => (d1, d2),
                (Err(fault), _) | (_, Err(fault)) => return (Err(fault), Closed::No),
            };
            error = Ok(());
            if d1[0].distance(d2[0]) > tolerances[0] {
                break;
            }
            match vectors_angle(d1[2], d2[2]) {
                None => {
                    error = Err(Fault::Degenerate);
                    break;
                }
                Some(angle) if angle > tolerances[1] => break,
                Some(_) => {}
            }
            jj += 1;
        }
        if jj <= q {
            break 'u;
        }
        ii += 1;
        if ii + q >= s {
            closed_u = true;
            break;
        }
        while ii < s && knots_v[ii] == knots_v[ii + 1] {
            ii += 1;
        }
        if ii >= s {
            return (Err(Fault::Degenerate), Closed::No);
        }
    }
    // V closure: the edges at the first and last V knots, along U.
    let mut ii = p;
    'v: loop {
        let du = (knots_u[ii + 1] - knots_u[ii]) / p as f64;
        let mut jj = 0;
        while jj <= p {
            let (first, last) = if jj == p {
                let u = knots_u[ii + 1];
                (
                    result(Vector2::new(u, knots_v[0]), [false, true]),
                    result(Vector2::new(u, knots_v[s]), [false, false]),
                )
            } else {
                let u = knots_u[ii] + jj as f64 * du;
                (
                    result(Vector2::new(u, knots_v[0]), [true, true]),
                    result(Vector2::new(u, knots_v[s]), [true, false]),
                )
            };
            let (d1, d2) = match (first, last) {
                (Ok(d1), Ok(d2)) => (d1, d2),
                (Err(fault), _) | (_, Err(fault)) => return (Err(fault), Closed::No),
            };
            error = Ok(());
            if d1[0].distance(d2[0]) > tolerances[0] {
                break;
            }
            match vectors_angle(d1[1], d2[1]) {
                None => {
                    error = Err(Fault::Degenerate);
                    break;
                }
                Some(angle) if angle > tolerances[1] => break,
                Some(_) => {}
            }
            jj += 1;
        }
        if jj <= p {
            break 'v;
        }
        ii += 1;
        if ii + p >= r {
            closed_v = true;
            break;
        }
        while ii < r && knots_u[ii] == knots_u[ii + 1] {
            ii += 1;
        }
        if ii >= r {
            return (Err(Fault::Degenerate), Closed::No);
        }
    }
    let closed = match (closed_u, closed_v) {
        (true, true) => Closed::UV,
        (true, false) => Closed::U,
        (false, true) => Closed::V,
        (false, false) => Closed::No,
    };
    (error, closed)
}

/// `N_SrfExtendByParamDist` for `NL_G1R`: the surface extended in direction
/// `axis` to the parameter `par` at its start or its end, by a ruled strip of
/// one span whose end control points continue the first or last two rows of
/// the net, then one knot removed from the old end if it can be. `s` is
/// changed in place, as NLib's surface is when its two arguments are the same,
/// and stays as far extended as the routine got when it returns an error.
pub(crate) fn extend_by_parameter_distance(
    s: &mut SurfaceParts,
    par: f64,
    axis: usize,
    at_start: bool,
) -> Result<(), Fault> {
    let degree = s.degrees[axis] as i64;
    let (n, m) = (s.counts[0] as i64 - 1, s.counts[1] as i64 - 1);
    let old = s.clone();
    let along = if axis == 0 { n } else { m };
    let new_along = along + degree;
    let mut extended = SurfaceParts {
        degrees: old.degrees,
        knots: old.knots.clone(),
        counts: if axis == 0 {
            [(new_along + 1) as usize, old.counts[1]]
        } else {
            [old.counts[0], (new_along + 1) as usize]
        },
        points: vec![
            ControlPoint {
                coordinates: [0.0; 4]
            };
            (if axis == 0 { new_along + 1 } else { n + 1 }) as usize
                * (if axis == 0 { m + 1 } else { new_along + 1 }) as usize
        ],
    };
    let last_knot = old.knots[axis].len() as i64 - 1;
    let shift = if at_start { degree } else { 0 };
    // The new knot vector: `degree + 1` copies of `par` at the new end, the
    // old knots beyond the old end.
    let mut knots = Vec::with_capacity((last_knot + degree + 1) as usize);
    if at_start {
        knots.extend(std::iter::repeat_n(par, degree as usize + 1));
        knots.extend_from_slice(&old.knots[axis][1..]);
    } else {
        knots.extend_from_slice(&old.knots[axis][..last_knot as usize]);
        knots.extend(std::iter::repeat_n(par, degree as usize + 1));
    }
    extended.knots[axis] = knots;
    // The control points of the old surface move along by `shift`.
    let across = if axis == 0 { m } else { n };
    let net_at = |parts: &SurfaceParts, line: i64, index: i64| {
        if axis == 0 {
            parts.at(index, line)
        } else {
            parts.at(line, index)
        }
    };
    let put = |parts: &mut SurfaceParts, line: i64, index: i64, point: ControlPoint| {
        if axis == 0 {
            parts.set(index, line, point);
        } else {
            parts.set(line, index, point);
        }
    };
    let knot = |index: i64| old.knots[axis][index as usize];
    let elevation = degree as f64;
    for line in 0..=across {
        for index in 0..=along {
            put(
                &mut extended,
                line,
                index + shift,
                net_at(&old, line, index),
            );
        }
        let (edge, inner, boundary, step) = if at_start {
            (
                0,
                1,
                0,
                -elevation * ((knot(0) - par) / (knot(degree + 1) - knot(0))),
            )
        } else {
            (
                along,
                along - 1,
                new_along,
                elevation
                    * ((par - knot(last_knot)) / (knot(last_knot) - knot(last_knot - degree - 1))),
            )
        };
        // The boundary row continues the line through the first (or last) two rows.
        let (p1, p2) = if at_start {
            (euclid(net_at(&old, line, 0)), euclid(net_at(&old, line, 1)))
        } else {
            (
                euclid(net_at(&old, line, along - 1)),
                euclid(net_at(&old, line, along)),
            )
        };
        let direction = p2 - p1;
        let origin = if at_start { p1 } else { p2 };
        let moved = origin + step * direction;
        let weight = net_at(&old, line, if at_start { 1 } else { inner }).coordinates[3];
        let point = if weight == NO_WEIGHT {
            ControlPoint {
                coordinates: [moved.x, moved.y, moved.z, NO_WEIGHT],
            }
        } else {
            ControlPoint {
                coordinates: [weight * moved.x, weight * moved.y, weight * moved.z, weight],
            }
        };
        put(&mut extended, line, boundary, point);
        // Degree elevation between the boundary row and the old edge.
        for k in 1..degree {
            let fraction = k as f64 / elevation;
            let (a, b, target) = if at_start {
                (boundary, edge + shift, k)
            } else {
                (edge, boundary, along + k)
            };
            let combined = combine(
                1.0 - fraction,
                net_at(&extended, line, a),
                fraction,
                if at_start {
                    net_at(&old, line, edge)
                } else {
                    net_at(&extended, line, b)
                },
            );
            put(&mut extended, line, target, combined);
        }
    }
    *s = extended;
    // Try to remove one knot from the old end.
    let old_end = if at_start { knot(0) } else { knot(last_knot) };
    remove_knot_in_place(s, axis, old_end, 1, MODEL_TOLERANCE).and_then(|(_, flagged)| {
        if flagged {
            Err(Fault::Degenerate)
        } else {
            Ok(())
        }
    })
}

/// The extension of a surface that `N_SrfEvalPtDerivsUnbounded` keeps in its
/// `esur`, and the closure flag it computed.
#[derive(Clone, Debug)]
pub struct ExtendedSurface {
    pub(crate) parts: SurfaceParts,
    /// `gflg`: where the surface closes smoothly, which it wraps rather than
    /// extends.
    pub closed: Closed,
}

impl ExtendedSurface {
    /// The extended surface.
    pub fn surface(&self) -> Result<BSplineSurface, Fault> {
        self.parts.clone().into_surface()
    }
}

/// How far, in a fraction of the knot range, one extension step goes.
const EXTENSION_STEP: f64 = 0.9;
/// `N_SrfEvalPtDerivsUnbounded`'s limit on a parameter.
const PARAMETER_LIMIT: f64 = 1.0e12;

/// `N_SrfEvalPtDerivsUnbounded` for a point with `NL_G1R` extension: the
/// surface point at `uv`, where a parameter outside the domain is wrapped
/// when the surface closes smoothly in that direction and otherwise reached by
/// extending it, in steps of at most nine tenths of its range. `extension`
/// is the routine's `esur` and `gflg`: empty before the first extension, it
/// is then the surface extended so far, and every later call extends that
/// one. On an error it holds what the routine reached, which a later call
/// goes on from.
pub fn evaluate_unbounded(
    surface: &BSplineSurface,
    extension: &mut Option<ExtendedSurface>,
    uv: Vector2,
) -> Result<Vector3, Fault> {
    let (mut u, mut v) = (uv.x, uv.y);
    let bounds = |parts: &SurfaceParts| {
        let [ku, kv] = &parts.knots;
        (ku[0], ku[ku.len() - 1], kv[0], kv[kv.len() - 1])
    };
    let current = match extension {
        Some(extended) => bounds(&extended.parts),
        None => {
            let [ku, kv] = surface.knots();
            (ku[0], ku[ku.len() - 1], kv[0], kv[kv.len() - 1])
        }
    };
    let (mut us, mut ue, mut vs, mut ve) = current;
    if u.abs() > PARAMETER_LIMIT || v.abs() > PARAMETER_LIMIT {
        return Err(Fault::Degenerate);
    }
    if u < us || u > ue || v < vs || v > ve {
        if extension.is_none() {
            // `N_SrfCopy`, then the closure test. A failing test leaves `gflg`
            // unset; the port reads it as not closed.
            let parts = SurfaceParts::of(surface);
            let (tested, closed) = is_closed_smooth(surface, [MODEL_TOLERANCE, 0.1]);
            *extension = Some(ExtendedSurface { parts, closed });
            tested?;
        }
        let extended = extension.as_mut().ok_or(Fault::Degenerate)?;
        let closed = extended.closed;
        let wrap = |value: &mut f64, low: f64, high: f64| -> Result<(), Fault> {
            // SMLib loops; a surface of no extent would loop for ever.
            let ordered = high > low;
            if !ordered || !value.is_finite() {
                return Err(Fault::Degenerate);
            }
            let mut turns = 0u32;
            while *value < low {
                *value += high - low;
                turns += 1;
                if turns > 1_000_000 {
                    return Err(Fault::Degenerate);
                }
            }
            while *value > high {
                *value += low - high;
                turns += 1;
                if turns > 1_000_000 {
                    return Err(Fault::Degenerate);
                }
            }
            Ok(())
        };
        if u < us || u > ue {
            if closed.in_u() {
                wrap(&mut u, us, ue)?;
            } else if u < us {
                loop {
                    let target = if us - u <= EXTENSION_STEP * (ue - us) {
                        u
                    } else {
                        us - EXTENSION_STEP * (ue - us)
                    };
                    extend_by_parameter_distance(&mut extended.parts, target, 0, true)?;
                    us = target;
                    if u >= us {
                        break;
                    }
                }
            } else {
                loop {
                    let target = if u - ue <= EXTENSION_STEP * (ue - us) {
                        u
                    } else {
                        ue + EXTENSION_STEP * (ue - us)
                    };
                    extend_by_parameter_distance(&mut extended.parts, target, 0, false)?;
                    ue = target;
                    if u <= ue {
                        break;
                    }
                }
            }
        }
        if v < vs || v > ve {
            if closed.in_v() {
                wrap(&mut v, vs, ve)?;
            } else if v < vs {
                loop {
                    let target = if vs - v <= EXTENSION_STEP * (ve - vs) {
                        v
                    } else {
                        vs - EXTENSION_STEP * (ve - vs)
                    };
                    extend_by_parameter_distance(&mut extended.parts, target, 1, true)?;
                    vs = target;
                    if v >= vs {
                        break;
                    }
                }
            } else {
                loop {
                    let target = if v - ve <= EXTENSION_STEP * (ve - vs) {
                        v
                    } else {
                        ve + EXTENSION_STEP * (ve - vs)
                    };
                    extend_by_parameter_distance(&mut extended.parts, target, 1, false)?;
                    ve = target;
                    if v <= ve {
                        break;
                    }
                }
            }
        }
    }
    // Now evaluate.
    match extension {
        Some(extended) => eval_point(&extended.surface()?, Vector2::new(u, v), [true, true]),
        None => eval_point(surface, Vector2::new(u, v), [true, true]),
    }
}
