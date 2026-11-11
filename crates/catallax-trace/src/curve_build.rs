// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: SMLib's contexts, NLib stacks and in-place swaps
// become owned values, statuses become `Fault`, and every curve is checked by
// `BSplineCurve::new`. Translated directly from SMLib SmBSplineCurve.cpp
// (`CreateLineSegment`, `CreateDegenerateCurve`, `CreatePointCurve`,
// `EditParameterization`, `EditEndPoint`), SmBSplineCurveNLib.cpp (`InsertOneKnot`,
// `InsertKnots`), SmNurbsCrv.cpp (`sm_FixupKnotVector`), NL_CrvNurb.cpp
// (`N_CrvLineFromPtAndVector`), NL_CrvTool.cpp (`N_CrvInsertKnot`,
// `N_CrvRefine`), NL_BasisBasic.cpp (`N_KnotVectorIsEndParam`,
// `N_BasisFindSpanAndMult`, `N_BasisReparam`) and NL_Math.cpp
// (`N_PtToCPt`, `N_Combine2CPts`). rust-brep has no port of these routines.
//! Curve construction and knot editing for the blend and fillet surface
//! generators.
//!
//! Every curve is an unproved candidate. SMLib builds a curve through
//! `SmBSplineCurve` constructors that run `sm_FixupKnotVector`; the
//! constructor-built results here do the same. A refinement keeps SMLib's
//! raw copy and does not.
//!
//! Where SMLib returns a curve the checked storage cannot hold, the port
//! refuses with [`Fault::InvalidInput`]: a knot vector whose active interval
//! is empty (a reparameterisation to a single value) and a result with a
//! non-finite control point (a refinement whose knots are not ordered, or
//! repeat a knot past the degree). The gate records each such case against
//! SMLib's result.

#![allow(clippy::needless_range_loop)]

use crate::{
    Fault,
    extent::Extent1,
    scalar::{EFF_ZERO, are_same, max},
    spline::{BSplineCurve, ControlPoint, NO_WEIGHT, NO_Z, combine, find_span},
    vector::Vector3,
};

/// `sm_FixupKnotVector` on `knots` with parameter-space tolerance
/// `tolerance`: knots closer than the scaled zero become equal, and distinct
/// knots closer than `tolerance` are spread by `1.01 * tolerance`.
pub(crate) fn fixup_knot_vector(knots: &mut [f64], tolerance: f64) {
    let count = knots.len();
    let (u0, u1) = (knots[0], knots[count - 1]);
    // `SmTol::GetScaledZero(u0, u1)`.
    let scaled_zero = max(EFF_ZERO * (1.0 + u0.abs()), EFF_ZERO * (1.0 + u1.abs()));
    let increment = tolerance * 1.01;
    let mut distinct = Vec::new();
    let mut multiplicities = Vec::new();
    let mut multiplicity = 1usize;
    for i in 0..count - 1 {
        if knots[i + 1] - knots[i] < scaled_zero {
            multiplicity += 1;
            knots[i + 1] = knots[i];
        } else {
            multiplicities.push(multiplicity);
            distinct.push(knots[i]);
            multiplicity = 1;
        }
    }
    multiplicities.push(multiplicity);
    distinct.push(knots[count - 1]);
    let last = distinct.len();
    for i in 0..last - 1 {
        if distinct[i + 1] - distinct[i] < tolerance {
            if i < 1 {
                distinct[i + 1] = distinct[i] + increment;
            } else if i + 2 >= last {
                distinct[i] = distinct[i + 1] - increment;
            } else {
                let mid = (distinct[i] + distinct[i + 1]) / 2.0;
                distinct[i] = mid - increment / 2.0;
                distinct[i + 1] = mid + increment / 2.0;
            }
        }
    }
    let mut at = 0;
    for (value, multiplicity) in distinct.into_iter().zip(multiplicities) {
        for _ in 0..multiplicity {
            knots[at] = value;
            at += 1;
        }
    }
}

/// `N_PtToCPt`, with `N_Crv3dTo2d` applied for a planar result.
fn control_point(point: Vector3, dimension: usize) -> ControlPoint {
    ControlPoint {
        coordinates: [
            point.x,
            point.y,
            if dimension == 2 { NO_Z } else { point.z },
            NO_WEIGHT,
        ],
    }
}

/// A constructor-built curve: SMLib's `SmBSplineCurve` constructors run
/// `sm_FixupKnotVector` with `SM_EFF_ZERO` on the new knots.
pub(crate) fn constructed(
    degree: usize,
    mut knots: Vec<f64>,
    points: Vec<ControlPoint>,
) -> Result<BSplineCurve, Fault> {
    fixup_knot_vector(&mut knots, EFF_ZERO);
    BSplineCurve::new(degree, knots, points)
}

impl BSplineCurve {
    /// `SmBSplineCurve::CreateLineSegment`: a degree-one curve from `start`
    /// to `end` in two or three dimensions. The second control point is
    /// `start + (end - start)`, as `N_CrvLineFromPtAndVector` forms it.
    /// Refuses a dimension outside 2 to 3 and ends closer than `SM_EFF_ZERO`.
    pub fn line_segment(dimension: usize, start: Vector3, end: Vector3) -> Result<Self, Fault> {
        if !(2..=3).contains(&dimension) {
            return Err(Fault::InvalidInput);
        }
        // `SAME_PNT_ER`: a NaN distance refuses too.
        let distance = start.distance(end);
        if distance.is_nan() || distance < EFF_ZERO {
            return Err(Fault::InvalidInput);
        }
        let vector = end - start;
        constructed(
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                control_point(start, dimension),
                control_point(start + vector, dimension),
            ],
        )
    }

    /// `SmBSplineCurve::CreateDegenerateCurve`, and with dimension three
    /// `CreatePointCurve`: a degree-one curve of two equal control points on
    /// the knots `0, 0, 1, 1`.
    pub fn degenerate(dimension: usize, point: Vector3) -> Result<Self, Fault> {
        if !(2..=3).contains(&dimension) {
            return Err(Fault::InvalidInput);
        }
        constructed(
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                control_point(point, dimension),
                control_point(point, dimension),
            ],
        )
    }

    /// `SmBSplineCurve::InsertOneKnot` through `N_CrvInsertKnot`: inserts
    /// `parameter` `times` more times. Refuses a parameter on or outside the
    /// end knots, and an insertion that would raise the knot's multiplicity
    /// above the degree.
    pub fn insert_knot(&self, parameter: f64, times: usize) -> Result<Self, Fault> {
        let (p, knots, points) = (self.degree(), self.knots(), self.control_points());
        let (n, m) = (points.len() as i64 - 1, knots.len() as i64 - 1);
        let (degree, r) = (p as i64, times as i64);
        // `N_KnotVectorIsEndParam`.
        if parameter <= knots[0] || parameter >= knots[m as usize] {
            return Err(Fault::Degenerate);
        }
        // `N_BasisFindSpanAndMult` from the left.
        let span = find_span(knots, p, parameter, false) as i64;
        let multiplicity = if parameter == knots[(m - degree) as usize] {
            degree + 1
        } else {
            let mut count = 0;
            let mut k = span;
            while k >= 0 {
                let equal = parameter == knots[k as usize];
                k -= 1;
                if !equal {
                    break;
                }
                count += 1;
            }
            count
        };
        if multiplicity + r > degree {
            return Err(Fault::Degenerate);
        }
        let mut q = vec![
            ControlPoint {
                coordinates: [0.0; 4]
            };
            (n + r + 1) as usize
        ];
        let mut uq = vec![0.0; (m + r + 1) as usize];
        let mut k = span - degree;
        let mut auxiliary: Vec<ControlPoint> = (0..=degree - multiplicity)
            .map(|i| points[(k + i) as usize])
            .collect();
        auxiliary.resize(
            (degree + 1) as usize,
            ControlPoint {
                coordinates: [0.0; 4],
            },
        );
        for i in 0..=k {
            q[i as usize] = points[i as usize];
        }
        for i in span - multiplicity..=n {
            q[(i + r) as usize] = points[i as usize];
        }
        for i in 1..=r {
            k = span - degree + i;
            for j in 0..=degree - i - multiplicity {
                let alpha = (parameter - knots[(k + j) as usize])
                    / (knots[(span + j + 1) as usize] - knots[(k + j) as usize]);
                let beta = 1.0 - alpha;
                auxiliary[j as usize] = combine(
                    alpha,
                    auxiliary[(j + 1) as usize],
                    beta,
                    auxiliary[j as usize],
                );
            }
            q[k as usize] = auxiliary[0];
            q[(span + r - i - multiplicity) as usize] =
                auxiliary[(degree - i - multiplicity) as usize];
        }
        for i in k + 1..span - multiplicity {
            q[i as usize] = auxiliary[(i - k) as usize];
        }
        for i in 0..=span {
            uq[i as usize] = knots[i as usize];
        }
        for i in 1..=r {
            uq[(i + span) as usize] = parameter;
        }
        for i in span + 1..=m {
            uq[(i + r) as usize] = knots[i as usize];
        }
        constructed(p, uq, q)
    }

    /// `SmBSplineCurve::InsertKnots` through `N_CrvRefine`: inserts every
    /// value of `new_knots`, which must be in ascending order and lie inside
    /// the end knots. An empty list returns the curve unchanged.
    pub fn refine(&self, new_knots: &[f64]) -> Result<Self, Fault> {
        if new_knots.is_empty() {
            return Ok(self.clone());
        }
        let (knots, points) = refine_parts(
            self.degree(),
            self.knots(),
            self.control_points(),
            new_knots,
        )?;
        BSplineCurve::new(self.degree(), knots, points)
    }

    /// `SmBSplineCurve::EditEndPoint`: moves the first or last control point
    /// to `point`. With a `weight` the point is homogeneous, `x`, `y` and `z`
    /// already multiplied by it; without one a rational end point is scaled by
    /// its own weight. A weight for a polynomial curve is refused. The port
    /// returns a new curve and leaves this one as it was; SMLib has already
    /// moved the coordinates when it refuses.
    pub fn edit_end_point(
        &self,
        point: Vector3,
        start: bool,
        weight: Option<f64>,
    ) -> Result<Self, Fault> {
        let mut points = self.control_points().to_vec();
        let end = if start { 0 } else { points.len() - 1 };
        let c = &mut points[end].coordinates;
        c[0] = point.x;
        c[1] = point.y;
        if c[2] != NO_Z {
            c[2] = point.z;
        }
        match weight {
            Some(weight) => {
                if c[3] == NO_WEIGHT {
                    return Err(Fault::Degenerate);
                }
                c[3] = weight;
            }
            None if c[3] != NO_WEIGHT => {
                c[0] *= c[3];
                c[1] *= c[3];
                if c[2] != NO_Z {
                    c[2] *= c[3];
                }
            }
            None => {}
        }
        BSplineCurve::new(self.degree(), self.knots().to_vec(), points)
    }

    /// `SmCurve::EditParameterization` through `N_BasisReparam`: the same
    /// curve on the natural interval `interval`. A request within SMLib's
    /// scaled zero of the current interval returns the curve unchanged.
    pub fn reparameterized(&self, interval: Extent1) -> Result<Self, Fault> {
        let old = self.natural_interval();
        if are_same(old.min, interval.min) && are_same(old.max, interval.max) {
            return Ok(self.clone());
        }
        let (p, mut knots) = (self.degree(), self.knots().to_vec());
        let m = knots.len() - 1;
        let (left, right) = (interval.min, interval.max);
        if knots[0] != left || knots[m] != right {
            let u0 = knots[0];
            let factor = (right - left) / (knots[m] - knots[0]);
            for i in 0..=p {
                knots[i] = left;
            }
            for i in p + 1..=m - p - 1 {
                knots[i] = left + factor * (knots[i] - u0);
            }
            for i in m - p..=m {
                knots[i] = right;
            }
        }
        BSplineCurve::new(p, knots, self.control_points().to_vec())
    }
}

/// `N_CrvRefine` on raw parts: the knots and control points after inserting
/// `x`, a non-empty ascending list inside the end knots.
pub(crate) fn refine_parts(
    p: usize,
    knots: &[f64],
    points: &[ControlPoint],
    x: &[f64],
) -> Result<(Vec<f64>, Vec<ControlPoint>), Fault> {
    let (n, m) = (points.len() as i64 - 1, knots.len() as i64 - 1);
    let degree = p as i64;
    let r = x.len() as i64 - 1;

    for end in [x[0], x[r as usize]] {
        if end <= knots[0] || end >= knots[m as usize] {
            return Err(Fault::Degenerate);
        }
    }
    let mut q = vec![
        ControlPoint {
            coordinates: [0.0; 4]
        };
        (n + r + 2) as usize
    ];
    let mut uq = vec![0.0; (m + r + 2) as usize];
    let a = find_span(knots, p, x[0], false) as i64;
    let b = find_span(knots, p, x[r as usize], false) as i64 + 1;
    for j in 0..=a {
        uq[j as usize] = knots[j as usize];
    }
    for j in b + degree..=m {
        uq[(j + r + 1) as usize] = knots[j as usize];
    }
    for j in 0..=a - degree {
        q[j as usize] = points[j as usize];
    }
    for j in b - 1..=n {
        q[(j + r + 1) as usize] = points[j as usize];
    }
    let mut i = b + degree - 1;
    let mut k = b + degree + r;
    for j in (0..=r).rev() {
        while x[j as usize] <= knots[i as usize] && i > a {
            q[(k - degree - 1) as usize] = points[(i - degree - 1) as usize];
            uq[k as usize] = knots[i as usize];
            k -= 1;
            i -= 1;
        }
        q[(k - degree - 1) as usize] = q[(k - degree) as usize];
        for l in 1..=degree {
            let t = k - degree + l;
            let mut alpha = uq[(k + l) as usize] - x[j as usize];
            if alpha.abs() <= 0.0 {
                q[(t - 1) as usize] = q[t as usize];
            } else {
                alpha /= uq[(k + l) as usize] - knots[(i - degree + l) as usize];
                let beta = 1.0 - alpha;
                q[(t - 1) as usize] = combine(alpha, q[(t - 1) as usize], beta, q[t as usize]);
            }
        }
        uq[k as usize] = x[j as usize];
        k -= 1;
    }
    Ok((uq, q))
}
