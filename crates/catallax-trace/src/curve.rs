// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: the `SmCurve` services the intersection chain
// enters, written for the B-spline candidates of `spline.rs`; statuses
// become `Fault` and output references become returned values. A failed
// evaluation keeps SMLib's fallback values (the previous sample, or an
// unclassified continuity).
// Translated from SMLib SmCurve.cpp (ApproximateLength, EvaluateContinuity,
// EvaluateGeometric, sm_UnitizeTangent, SmFindCurveDistanceExtremaEFO) and
// SmBSplineCurve.cpp (CalculateContinuities, GetPolygon,
// GetNumberNaturalKnots); rust-brep has no port of these.
//! Curve services shared by the curve cache and the curve solvers.

use crate::Fault;
use crate::extent::Extent1;
use crate::geometry::{Continuity, CurveJet, evaluate_curve_continuity, frenet_from_derivatives};
use crate::local_solve::EvalFunction;
use crate::scalar::{ANGLE_TOL_DEG, EFF_ZERO};
use crate::spline::{BSplineCurve, unique_knots};
use crate::vector::Vector3;

/// `GetNumberNaturalKnots`: every knot, repeats included.
pub fn natural_knot_count(curve: &BSplineCurve) -> usize {
    curve.knots().len()
}

/// `ApproximateLength` for 3D use: the polyline length through `samples`
/// equally spaced parameters after the start. A sample that cannot be
/// evaluated repeats the previous point, as in SMLib.
pub fn approximate_length(curve: &BSplineCurve, interval: Extent1, samples: usize) -> f64 {
    let mut length = 0.0;
    let mut last = curve.evaluate_point(interval.min).unwrap_or_default();
    let mut current = Vector3::default();
    for ii in 1..=samples {
        let parameter = if ii < samples {
            interval.evaluate(ii as f64 / samples as f64)
        } else {
            interval.max
        };
        if let Ok(point) = curve.evaluate_point(parameter) {
            current = point;
        }
        length += current.distance(last);
        last = current;
    }
    length
}

/// `GetPolygon`: the Euclidean control points that influence `interval`.
pub fn polygon(curve: &BSplineCurve, interval: Extent1) -> Result<Vec<Vector3>, Fault> {
    let natural = curve.natural_interval();
    let points = curve.control_points();
    if points.len() > 1 && !interval.is_contained_by(natural, EFF_ZERO) {
        return Err(Fault::InvalidInput);
    }
    let (start, count) =
        if (natural.min == interval.min && natural.max == interval.max) || points.len() == 1 {
            (0, points.len())
        } else {
            let start = curve.find_span(interval.min, true)? - curve.degree();
            let end = curve.find_span(interval.max, false)?;
            (start, end - start + 1)
        };
    Ok(points[start..start + count]
        .iter()
        .map(|p| p.euclidean())
        .collect())
}

/// `EvaluateGeometric`: position and Frenet vectors from exact derivatives,
/// for `count` up to 3 (tangent, normal and binormal); a larger count is
/// refused.
pub fn evaluate_geometric(
    curve: &dyn crate::curve_solve::PointSolveCurve,
    parameter: f64,
    count: usize,
    from_left: bool,
) -> Result<[Vector3; 4], Fault> {
    if count > 3 {
        return Err(Fault::InvalidInput);
    }
    let values = curve.evaluate(parameter, count, from_left, false)?;
    if values.len() != count + 1 {
        return Err(Fault::InvalidInput);
    }
    let mut derivatives = [Vector3::default(); 4];
    derivatives[..values.len()].copy_from_slice(&values);
    let mut frame = [Vector3::default(); 4];
    frenet_from_derivatives(&derivatives, count, &mut frame)?;
    Ok(frame)
}

/// One side of [`evaluate_continuity`]: a curve, a parameter and which side
/// of a knot to evaluate from.
#[derive(Clone, Copy, Debug)]
pub struct CurveSide<'a> {
    /// The curve.
    pub curve: &'a BSplineCurve,
    /// The parameter.
    pub parameter: f64,
    /// Evaluate from the left of a knot.
    pub from_left: bool,
}

/// `EvaluateContinuity`: the continuity class where two curve evaluations
/// meet, with both position, first and second derivative triples.
pub fn evaluate_continuity(
    this: CurveSide<'_>,
    other: CurveSide<'_>,
    angle_degrees: f64,
) -> Result<(Continuity, [Vector3; 3], [Vector3; 3]), Fault> {
    let a = this
        .curve
        .evaluate(this.parameter, 3, this.from_left, false)?;
    let b = other
        .curve
        .evaluate(other.parameter, 3, other.from_left, false)?;
    let jet = |v: &[Vector3]| CurveJet {
        position: v[0],
        first: v[1],
        second: v[2],
        third: Some(v[3]),
    };
    let class = evaluate_curve_continuity(&jet(&a), &jet(&b), angle_degrees);
    Ok((class, [a[0], a[1], a[2]], [b[0], b[1], b[2]]))
}

/// `CalculateContinuities` with SMLib's default angle: the class at every
/// unique knot (the ends are discontinuous) and the smallest interior class.
/// A knot that cannot be evaluated stays unclassified, as SMLib leaves it.
pub fn calculate_continuities(curve: &BSplineCurve) -> (Continuity, Vec<Continuity>) {
    let knots = unique_knots(curve.knots(), None);
    let mut classes = vec![Continuity::Undefined; knots.len()];
    classes[0] = Continuity::Discontinuous;
    let last = knots.len() - 1;
    classes[last] = Continuity::Discontinuous;
    let mut minimum = Continuity::CInfinity;
    for i in 1..last {
        let knot = knots[i].0;
        let side = |from_left| CurveSide {
            curve,
            parameter: knot,
            from_left,
        };
        if let Ok((class, _, _)) = evaluate_continuity(side(true), side(false), ANGLE_TOL_DEG) {
            classes[i] = class;
        }
        if classes[i] < minimum {
            minimum = classes[i];
        }
    }
    (minimum, classes)
}

/// `sm_UnitizeTangent`: the unit tangent from position, first and second
/// derivatives. A vanishing tangent falls back to the second derivative when
/// it points along `gap`, then to `gap` itself, else to zero. Also returns
/// the first-derivative length.
pub fn unitize_tangent(pvv: &[Vector3; 3], gap: Vector3) -> (Vector3, f64) {
    let length = pvv[1].length();
    if length > EFF_ZERO {
        return (pvv[1] / length, length);
    }
    if pvv[2].dot(gap) > EFF_ZERO {
        let second = pvv[2].length();
        if second > EFF_ZERO {
            return (pvv[2] / second, length);
        }
        let gap_length = gap.length();
        if gap_length > EFF_ZERO {
            return (gap / gap_length, length);
        }
    }
    (Vector3::new(0.0, 0.0, 0.0), length)
}

/// `SmFindCurveDistanceExtremaEFO`: for two curves sharing a parameter,
/// `F(t) = (C2 - C1) · (C2' - C1')`, which vanishes where their distance is
/// extreme.
pub struct CurveDistanceExtrema<'a> {
    /// First curve.
    pub first: &'a BSplineCurve,
    /// Second curve.
    pub second: &'a BSplineCurve,
    /// `|F|` at or below this is reported as an answer.
    pub tolerance: f64,
}

impl EvalFunction for CurveDistanceExtrema<'_> {
    fn evaluate(
        &mut self,
        t: f64,
        value: &mut f64,
        derivative: &mut f64,
        found: &mut bool,
        _signal_errors: bool,
    ) -> Result<(), Fault> {
        *value = 0.0;
        *derivative = 0.0;
        *found = false;
        let a = self.first.evaluate(t, 2, true, true)?;
        let b = self.second.evaluate(t, 2, true, true)?;
        *value = (b[0] - a[0]).dot(b[1] - a[1]);
        *derivative = (b[1] - a[1]).dot(b[1] - a[1]) + (b[0] - a[0]).dot(b[2] - a[2]);
        *found = value.abs() <= self.tolerance;
        Ok(())
    }
}
