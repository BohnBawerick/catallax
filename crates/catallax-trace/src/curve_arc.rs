// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: statuses become `Fault` and the output curve is an
// owned, checked `BSplineCurve`. Translated directly from SMLib
// SmBSplineCurve.cpp (`CreateCircleSegment`), NL_CrvCommon.cpp
// (`N_CreateCircArc`, `N_CreateCircArcTol`, `N_CreateQuadraticArcTol`,
// `N_CreateQuinticArcTol`), NL_CrvGeom.cpp (`N_CrvClosestPt`), NL_CrvTool.cpp
// (`N_CrvExtractCrvSeg`), NL_CrvBasic.cpp (`N_CrvGetConicData`), NL_Geom.cpp
// (`N_IsectLineLine`, `N_VectorNormalize`, `N_VectorNormalizeRef`) and
// SmNurbsCrv.cpp (`sm_ExtractFrame`). rust-brep has no port of these.
//! Exact circular arcs as rational B-splines.
//!
//! The quadratic arc is a chain of up to four quarter-or-smaller conic
//! spans with double interior knots. The quintic circle is a single span; a
//! quintic arc that is not a whole circle is cut from it by
//! `N_CrvClosestPt` and `N_CrvExtractCrvSeg`.
//!
//! Sines and cosines come from `libm` so a candidate is the same on every
//! machine; SMLib's come from the C library and can differ in the last bits.

use crate::{
    Fault,
    curve_build::constructed,
    curve_compat::division_is_bad,
    curve_nlib::{MODEL_TOLERANCE, closest_point_on_curve, extract_curve_segment},
    scalar::{EFF_ZERO, PI},
    spline::{BSplineCurve, ControlPoint, NO_Z, PlaneFrame},
    vector::Vector3,
};

/// `NL_ZCTL`, the zero-coefficient tolerance of `N_IsectLineLine`.
const ZERO_COEFFICIENT: f64 = 1.0e-12;
/// The angle tolerance `N_CreateCircArc` passes, in degrees.
const ANGLE_TOLERANCE: f64 = 0.01;

/// `SmNurbCircleParam`: the parameterisation of a circle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CircleParameterization {
    /// `SM_CO_QUADRATIC`: degree 2 with double interior knots.
    Quadratic,
    /// `SM_CO_QUINTIC`: degree 5 with no interior knots.
    Quintic,
}

fn dot(a: Vector3, b: Vector3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

/// `N_VectorNormalizeRef`: scales by the reciprocal of the length.
fn normalize_ref(v: Vector3) -> Result<Vector3, Fault> {
    let length = dot(v, v).sqrt();
    if division_is_bad(1.0, length) {
        return Err(Fault::Degenerate);
    }
    let inverse = 1.0 / length;
    Ok(Vector3::new(inverse * v.x, inverse * v.y, inverse * v.z))
}

/// `N_Combine2Pts`.
fn combine(alpha: f64, a: Vector3, beta: f64, b: Vector3) -> Vector3 {
    Vector3::new(
        alpha * a.x + beta * b.x,
        alpha * a.y + beta * b.y,
        alpha * a.z + beta * b.z,
    )
}

/// `N_TranslateSum2Pts`.
fn translate_sum(c: Vector3, alpha: f64, a: Vector3, beta: f64, b: Vector3) -> Vector3 {
    Vector3::new(
        c.x + alpha * a.x + beta * b.x,
        c.y + alpha * a.y + beta * b.y,
        c.z + alpha * a.z + beta * b.z,
    )
}

/// `N_Weight` with a weight other than `NL_NOW`.
fn weighted(point: Vector3, weight: f64) -> ControlPoint {
    ControlPoint {
        coordinates: [weight * point.x, weight * point.y, weight * point.z, weight],
    }
}

/// `N_CrvGetConicData`: the highest control-point index `n` and the end
/// angle brought within one turn of `start`.
fn conic_data(kind: CircleParameterization, start: f64, mut end: f64) -> (usize, f64) {
    const TOLERANCE: f64 = 1.0e-6;
    while end < start + TOLERANCE {
        end += 360.0;
    }
    while end - start > 360.0 - TOLERANCE {
        end -= 360.0;
    }
    if end - start < TOLERANCE {
        end = start + 360.0;
    }
    let theta = end - start;
    let n = match kind {
        CircleParameterization::Quadratic => {
            if theta < 91.0 {
                2
            } else if theta < 181.0 {
                4
            } else if theta < 271.0 {
                6
            } else {
                8
            }
        }
        CircleParameterization::Quintic => 5,
    };
    (n, end)
}

/// `N_IsectLineLine` for two unbounded lines: the point midway between the
/// closest points, `None` for parallel lines, `Err` for a numeric failure.
fn intersect_lines(
    p1: Vector3,
    v1: Vector3,
    p2: Vector3,
    v2: Vector3,
) -> Result<Option<Vector3>, Fault> {
    let length = |v: Vector3| dot(v, v).sqrt();
    let (m1, m2) = (length(v1), length(v2));
    if division_is_bad(1.0, m1) || division_is_bad(1.0, m2) {
        return Err(Fault::Degenerate);
    }
    let (inverse1, inverse2) = (1.0 / m1, 1.0 / m2);
    let v1 = Vector3::new(inverse1 * v1.x, inverse1 * v1.y, inverse1 * v1.z);
    let v2 = Vector3::new(inverse2 * v2.x, inverse2 * v2.y, inverse2 * v2.z);
    let p21 = p2 - p1;
    let (a00, a01, a11) = (dot(v1, v1), dot(v1, v2), dot(v2, v2));
    let (b0, b1) = (dot(v1, p21), dot(v2, p21));
    let denominator = a00 * a11 - a01 * a01;
    if denominator.abs() < ZERO_COEFFICIENT {
        return Ok(None);
    }
    let numerator_a = b0 * a11 - b1 * a01;
    let numerator_b = b0 * a01 - b1 * a00;
    if division_is_bad(numerator_a, denominator) || division_is_bad(numerator_b, denominator) {
        return Err(Fault::Degenerate);
    }
    let (alpha, beta) = (numerator_a / denominator, numerator_b / denominator);
    if division_is_bad(alpha, m1) || division_is_bad(beta, m2) {
        return Err(Fault::Degenerate);
    }
    let q1 = Vector3::new(
        p1.x + alpha * v1.x,
        p1.y + alpha * v1.y,
        p1.z + alpha * v1.z,
    );
    let q2 = Vector3::new(p2.x + beta * v2.x, p2.y + beta * v2.y, p2.z + beta * v2.z);
    Ok(Some(combine(0.5, q1, 0.5, q2)))
}

/// `N_CreateQuadraticArcTol`: knots and weighted control points.
fn quadratic_arc(
    c: Vector3,
    x: Vector3,
    y: Vector3,
    r: f64,
    start: f64,
    end: f64,
) -> Result<(Vec<f64>, Vec<ControlPoint>), Fault> {
    let (n, end) = conic_data(CircleParameterization::Quadratic, start, end);
    let arcs = n / 2;
    let theta = end - start;
    let delta = theta / arcs as f64;
    let (x, y) = (normalize_ref(x)?, normalize_ref(y)?);
    let mut angle = (PI * delta) / 360.0;
    let w1 = libm::cos(angle);
    angle = (PI * start) / 180.0;
    let (mut cosa, mut sina) = (libm::cos(angle), libm::sin(angle));
    let mut p0 = translate_sum(c, r * cosa, x, r * sina, y);
    let mut t0 = combine(-sina, x, cosa, y);
    let mut points = vec![weighted(p0, 1.0)];
    let step = (PI * delta) / 180.0;
    for i in 1..=arcs {
        angle += step;
        cosa = libm::cos(angle);
        sina = libm::sin(angle);
        let p2 = translate_sum(c, r * cosa, x, r * sina, y);
        let t2 = combine(-sina, x, cosa, y);
        // Lines that do not meet make a very flat span: use its midpoint.
        let p1 = match intersect_lines(p0, t0, p2, t2)? {
            Some(point) => point,
            None => combine(0.5, p0, 0.5, p2),
        };
        points.push(weighted(p1, w1));
        points.push(weighted(p2, 1.0));
        if i < arcs {
            p0 = p2;
            t0 = t2;
        }
    }
    if theta > 360.0 - ANGLE_TOLERANCE {
        points[n] = points[0];
    }
    let mut knots = vec![0.0; n + 4];
    match arcs {
        2 => knots[3..5].fill(0.5),
        3 => {
            knots[3..5].fill(1.0 / 3.0);
            knots[5..7].fill(2.0 / 3.0);
        }
        4 => {
            knots[3..5].fill(0.25);
            knots[5..7].fill(0.50);
            knots[7..9].fill(0.75);
        }
        _ => {}
    }
    knots[n + 1..].fill(1.0);
    Ok((knots, points))
}

/// `N_CreateQuinticArcTol`: the quintic circle, cut to the arc between the
/// angles when they span less than a whole turn.
fn quintic_arc(
    c: Vector3,
    x: Vector3,
    y: Vector3,
    r: f64,
    start: f64,
    end: f64,
) -> Result<(Vec<f64>, Vec<ControlPoint>), Fault> {
    let (_, end) = conic_data(CircleParameterization::Quintic, start, end);
    let middle = 0.5 * (start + end) + 180.0;
    let (x, y) = (normalize_ref(x)?, normalize_ref(y)?);
    let on_circle = |degrees: f64| {
        let angle = (PI * degrees) / 180.0;
        translate_sum(c, r * libm::cos(angle), x, r * libm::sin(angle), y)
    };
    let (start_point, end_point) = (on_circle(start), on_circle(end));
    let angle = (PI * middle) / 180.0;
    let (cosa, sina) = (libm::cos(angle), libm::sin(angle));
    let pm = translate_sum(c, r * cosa, x, r * sina, y);
    let t = normalize_ref(combine(-sina, x, cosa, y))?;
    let s = normalize_ref(combine(-cosa, x, -sina, y))?;
    let mut points = vec![weighted(pm, 5.0); 6];
    points[1] = weighted(combine(1.0, pm, 4.0 * r, t), 1.0);
    points[4] = weighted(combine(1.0, pm, -4.0 * r, t), 1.0);
    points[2] = weighted(translate_sum(c, 3.0 * r, s, 2.0 * r, t), 1.0);
    points[3] = weighted(translate_sum(c, 3.0 * r, s, -2.0 * r, t), 1.0);
    let mut knots = vec![0.0; 12];
    knots[6..].fill(1.0);
    if end - start > 360.0 - ANGLE_TOLERANCE {
        points[5] = points[0];
        return Ok((knots, points));
    }
    let circle = BSplineCurve::new(5, knots, points)?;
    let guess = |degrees: f64| degrees / 360.0;
    let first = 360.0 + (start - middle);
    let (lower, _) = closest_point_on_curve(
        &circle,
        start_point,
        guess(first),
        MODEL_TOLERANCE,
        MODEL_TOLERANCE,
    )?;
    let mut second = end - middle;
    if second < 0.0 {
        second += 360.0;
    }
    let (upper, _) = closest_point_on_curve(
        &circle,
        end_point,
        guess(second),
        MODEL_TOLERANCE,
        MODEL_TOLERANCE,
    )?;
    let arc = extract_curve_segment(&circle, lower, upper)?;
    Ok((arc.knots().to_vec(), arc.control_points().to_vec()))
}

impl BSplineCurve {
    /// `SmBSplineCurve::CreateCircleSegment`: a rational curve for the arc of
    /// the circle of `radius` in the XY plane of the frame, from `start` to
    /// `end` degrees counter-clockwise from the X axis. The frame is taken
    /// as `SmAxis2Placement::SetCanonical` leaves it: unit, perpendicular
    /// axes made from the ones given.
    ///
    /// Refuses a radius below `SM_EFF_ZERO`, a dimension outside 2 to 3 and an
    /// angle outside `[-360, 360]`. A quintic arc short of a whole circle also
    /// refuses; see the module notes.
    #[allow(clippy::too_many_arguments)]
    pub fn circle_segment(
        dimension: usize,
        origin: Vector3,
        x_axis: Vector3,
        y_axis: Vector3,
        radius: f64,
        start_degrees: f64,
        end_degrees: f64,
        parameterization: CircleParameterization,
    ) -> Result<Self, Fault> {
        // `LE_ZERO_ER`, then `RANGE_ER`: a NaN fails each.
        if radius.is_nan() || radius < EFF_ZERO {
            return Err(Fault::InvalidInput);
        }
        let in_range = |value: f64, low: f64, high: f64| value >= low && value <= high;
        if !in_range(dimension as f64, 2.0, 3.0)
            || !in_range(start_degrees, -360.0, 360.0)
            || !in_range(end_degrees, -360.0, 360.0)
        {
            return Err(Fault::InvalidInput);
        }
        let [c, x, y] = PlaneFrame::canonical(origin, x_axis, y_axis)?;
        let (knots, mut points) = match parameterization {
            CircleParameterization::Quadratic => {
                quadratic_arc(c, x, y, radius, start_degrees, end_degrees)?
            }
            CircleParameterization::Quintic => {
                quintic_arc(c, x, y, radius, start_degrees, end_degrees)?
            }
        };
        if dimension == 2 {
            for point in &mut points {
                point.coordinates[2] = NO_Z;
            }
        }
        let degree = match parameterization {
            CircleParameterization::Quadratic => 2,
            CircleParameterization::Quintic => 5,
        };
        constructed(degree, knots, points)
    }
}
