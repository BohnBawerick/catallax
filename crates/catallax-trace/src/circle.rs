// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: the circle owns its knot-doubled spline instead of
// inheriting from `SmBSplineCurve`; contexts, owners and notifications are
// dropped; statuses become `Fault`; the polar converter is built once, as
// the constructor's notification builds it; only the pieces
// `CreateAnalyticCurve` and `DropPointFast` run for minimising and meeting
// are carried. A circle whose doubled spline cannot be converted, which SMLib
// rebuilds with `MakeNurb`, is reported rather than rebuilt.
// Translated from SMLib SmCurve.cpp (IsArc; LocalPropertyAnalysis and
// SmCurvePropertyEFO::Evaluate for point inversion), SmBSplineCurve.cpp
// (CreateAnalyticCurve's arc), SmCircle.cpp and SmEllipse.cpp (the
// constructor from a NURBS curve, SmPolarConversion::SetUpPolarConversion,
// ConvertToNURBSParameter and AngleParamToNurbParam, STEPInversion,
// EvaluateSTEP, TransformPointToStartPlane, IsCircle, DropPointFast),
// SmGeomUtility.cpp (smgu_CircleCenterFrom3Points,
// smgu_TransformPointToStartPlane, smgu_PlaneDropPoint,
// smgu_PlaneEvaluatePoint) and SmExtent1d.cpp (ContainsPeriodicValue,
// ClampPeriodicValue, IsValueOnPeriodicBoundary) at
// c9c5979f04cea8c38213004d4f04df0dcfb40d05; rust-brep has no port of these.
//! SMLib's analytic circle for a rational quadratic arc: the arc test, the
//! `SmCircle` that `SmCurveCache::GetAnalyticCurve` builds for it, and the
//! circle's fast point drop.
//!
//! The circle keeps the arc's spline with every interior knot doubled, so
//! each span is one conic piece, and SMLib's polar converter between its
//! angle in degrees and that spline's parameter. Drops answer spline
//! parameters; they propose candidates and prove nothing.
//!
//! Sines, cosines and arc tangents come from `libm`, so a candidate is the
//! same on every machine; SMLib's come from the C library and can differ in
//! the last bits.

use crate::Fault;
use crate::extent::Extent1;
use crate::geometry::{line_point_distance, plane_point_distance};
use crate::local_solve::{EvalFunction, LocalSolve1d};
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, PI, ZONE_TOL_3D, are_same, cos, sin};
use crate::spline::{BSplineCurve, PlaneFrame, PointForm, unique_knots};
use crate::surface_curves::ccw_angle_between;
use crate::vector::{Vector2, Vector3};

/// `SM_EFF_ZERO_DEG`.
const EFF_ZERO_DEG: f64 = 1.0e-4;
/// `ANALYTIC_TOL_SCALE`.
const ANALYTIC_TOL_SCALE: f64 = 100.0;
/// The sample count `CreateAnalyticCurve` gives the arc test.
const ARC_SAMPLES: usize = 20;

/// `SM_RAD2DEG`.
fn rad2deg(radians: f64) -> f64 {
    radians * (180.0 / PI)
}

/// `SM_ARE_SAME_TO_TOL`.
fn are_same_to_tol(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

/// `SmExtent1d::ContainsPeriodicValue`'s outputs.
#[derive(Clone, Copy, Debug)]
struct Periodic {
    inside: bool,
    value: f64,
    to_min: f64,
    to_max: f64,
}

/// `SmExtent1d::ContainsPeriodicValue`: whether `value`, moved by whole
/// periods, lies in `extent` within `tol` (`SM_EFF_ZERO_DEG` when not
/// positive), with that value and its distances to the ends.
fn contains_periodic(extent: Extent1, value: f64, period: f64, tol: f64) -> Periodic {
    let zero = if tol > 0.0 { tol } else { EFF_ZERO_DEG };
    let inside = extent.contains_value(value, zero);
    if inside || period == 0.0 {
        return Periodic {
            inside,
            value,
            to_min: (value - extent.min).abs(),
            to_max: (value - extent.max).abs(),
        };
    }
    let base = Extent1::new(extent.min, extent.min + period);
    let mut mapped = value;
    if !base.contains_value(value, zero) {
        mapped = crate::local_solve::periodic_wrap(base, value);
    }
    let mut inside = mapped > extent.min - zero && mapped < extent.max + zero;
    if !inside
        && are_same_to_tol(mapped, base.max, zero)
        && are_same_to_tol(mapped - period, extent.min, zero)
    {
        mapped = extent.min;
        inside = true;
    }
    let periodic =
        if inside && ((value - extent.min).abs() < zero || (value - extent.max).abs() < zero) {
            value
        } else if inside && mapped > extent.max + zero {
            mapped - period
        } else {
            mapped
        };
    let mut to_min = mapped - extent.min;
    if to_min > period / 2.0 {
        to_min = period - to_min;
    }
    let mut to_max = (extent.max - mapped).abs();
    if to_max > period / 2.0 {
        to_max = period - to_max;
    }
    Periodic {
        inside,
        value: periodic,
        to_min,
        to_max,
    }
}

/// `SmExtent1d::ClampPeriodicValue`: the periodic value inside `extent`, or
/// the nearer end.
fn clamp_periodic(extent: Extent1, value: f64, period: f64, tol: f64) -> f64 {
    let found = contains_periodic(extent, value, period, tol);
    if found.inside {
        found.value
    } else if found.to_min < found.to_max {
        extent.min
    } else {
        extent.max
    }
}

/// `SmExtent1d::IsValueOnPeriodicBoundary`.
fn is_on_periodic_boundary(extent: Extent1, value: f64, period: f64, tol: f64) -> bool {
    let found = contains_periodic(extent, value, period, tol);
    found.to_min < tol || found.to_max < tol
}

/// `smgu_CircleCenterFrom3Points`: the centre of the circle through three
/// plane points, `None` for coincident or nearly collinear points.
fn circle_center_from_3_points(a: Vector2, b: Vector2, c: Vector2) -> Option<Vector2> {
    let tol = EFF_ZERO_SQRT * (1.0 + a.max_dimension() + c.max_dimension());
    let tol2 = tol * tol;
    if a.distance_squared(b) < tol2 || b.distance_squared(c) < tol2 || a.distance_squared(c) < tol2
    {
        return None;
    }
    // `smgu_LinePointDistance` in two dimensions.
    let along = c - a;
    let length_squared = along.length_squared();
    let parameter = if length_squared < EFF_ZERO_SQ {
        0.0
    } else {
        (b - a).dot(along) / length_squared
    };
    let to_line = (a + parameter * along).distance(b);
    if to_line < 2.0 * EFF_ZERO * a.distance_squared(c) {
        return None;
    }
    let a00 = 2.0 * (a.x - b.x);
    let a01 = 2.0 * (a.y - b.y);
    let a10 = 2.0 * (b.x - c.x);
    let a11 = 2.0 * (b.y - c.y);
    let b0 = a.dot(a) - b.dot(b);
    let b1 = b.dot(b) - c.dot(c);
    let determinant = a00 * a11 - a10 * a01;
    let numerator0 = b0 * a11 - b1 * a01;
    let numerator1 = a00 * b1 - a10 * b0;
    if determinant.abs() < numerator0.abs() * EFF_ZERO
        || determinant.abs() < numerator1.abs() * EFF_ZERO
    {
        return None;
    }
    Some(Vector2::new(
        numerator0 / determinant,
        numerator1 / determinant,
    ))
}

/// `SmAxis2Placement::GetZAxis`: `x × y`, unitised when it can be.
fn z_axis(x: Vector3, y: Vector3) -> Vector3 {
    let mut z = x.cross(y);
    let _ = z.unitize();
    z
}

/// The arc `SmCurve::IsArc` finds: its frame (centre, unit X axis towards
/// the start, unit Y axis), mean radius and end angle in degrees from a
/// start angle of zero.
#[derive(Clone, Copy, Debug)]
struct ArcFrame {
    origin: Vector3,
    x_axis: Vector3,
    y_axis: Vector3,
    radius: f64,
    end_degrees: f64,
}

/// `SmCurve::IsArc`: whether `samples` points of `curve` lie within
/// `tolerance` of the circle through three of them and of its plane.
fn is_arc(curve: &BSplineCurve, samples: usize, tolerance: f64) -> Result<Option<ArcFrame>, Fault> {
    if curve.is_degenerate(tolerance) {
        return Ok(None);
    }
    let interval = curve.natural_interval();
    let p1 = curve.evaluate_point(interval.evaluate(0.00))?;
    let p2 = curve.evaluate_point(interval.evaluate(0.45))?;
    let p3 = curve.evaluate_point(interval.evaluate(0.90))?;
    let v1 = p2 - p1;
    let v2 = p3 - p1;
    let approximate_length = p1.distance(p2) + p2.distance(p3);
    if v1.length_squared() < EFF_ZERO_SQ {
        return Ok(None);
    }
    let Ok([x, y, z]) = v1.make_unit_ortho_vectors(Some(v2)) else {
        return Ok(None);
    };
    let [origin, x, y] = PlaneFrame::canonical(p1, x, y)?;
    let drop = |p: Vector3| {
        let offset = p - origin;
        Vector2::new(x.dot(offset), y.dot(offset))
    };
    let Some(uv) = circle_center_from_3_points(drop(p1), drop(p2), drop(p3)) else {
        return Ok(None);
    };
    let center = origin + uv.x * x + uv.y * y;
    let radius = center.distance(p1);
    if radius > approximate_length / EFF_ZERO_SQRT {
        return Ok(None);
    }
    let mut total = radius;
    for i in 1..samples {
        let parameter = if i == samples - 1 {
            interval.max
        } else {
            interval.evaluate(i as f64 / (samples as f64 - 1.0))
        };
        let Ok(point) = curve.evaluate_point(parameter) else {
            return Ok(None);
        };
        let test = point.distance(center);
        if (test - radius).abs() > tolerance {
            return Ok(None);
        }
        total += test;
        let Ok(to_plane) = plane_point_distance(center, z, point) else {
            return Ok(None);
        };
        if to_plane > tolerance {
            return Ok(None);
        }
    }
    let radius = total / samples as f64;
    let Ok([x, y, z]) = (p1 - center).make_unit_ortho_vectors(Some(p2 - center)) else {
        return Ok(None);
    };
    let [origin, x_axis, y_axis] = PlaneFrame::canonical(center, x, y)?;
    let end = curve.evaluate_point(interval.max)?;
    let mut end_degrees = ccw_angle_between(z, x, end - center)? * 180.0 / PI;
    if end_degrees.abs() < EFF_ZERO_SQRT {
        end_degrees = 360.0;
    } else if end_degrees < 0.0 {
        end_degrees += 360.0;
    }
    Ok(Some(ArcFrame {
        origin,
        x_axis,
        y_axis,
        radius,
        end_degrees,
    }))
}

/// `SmPolarConversion`'s tables: for every distinct knot of the circle's
/// spline its angle in degrees, and for every span whether its weights make
/// the conic conversion exact.
#[derive(Clone, Debug, PartialEq)]
struct PolarConversion {
    knots: Vec<f64>,
    angles: Vec<f64>,
    exact: Vec<bool>,
}

/// What `SmBSplineCurve::CreateAnalyticCurve` makes of a rational quadratic.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum AnalyticArc {
    /// Not a circular arc, or one whose circle does not reproduce the
    /// curve's ends: SMLib keeps the spline.
    Spline,
    /// The analytic circle.
    Circle(Box<Circle>),
    /// A circle whose knot-doubled spline the polar converter cannot take,
    /// which SMLib rebuilds from the circle with `MakeNurb`.
    Rebuilt,
}

/// `SmCircle` built from a rational quadratic arc.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Circle {
    origin: Vector3,
    x_axis: Vector3,
    y_axis: Vector3,
    /// `m_vAnalDomain`, in degrees.
    domain: Extent1,
    radius: f64,
    /// The arc's spline with every interior knot doubled.
    nurb: BSplineCurve,
    polar: PolarConversion,
}

/// One answer of [`Circle::drop_point_fast`]: a spline parameter and its
/// distance from the test point.
pub(crate) type CircleAnswer = (f64, f64);

/// The outcome of [`Circle::drop_point_fast`]. `complete` is false where
/// SMLib's drop fails part-way (`SM_ERR`), keeping the answers it added.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct CircleDrop {
    pub answers: Vec<CircleAnswer>,
    pub complete: bool,
}

/// `SmBSplineCurve::CreateAnalyticCurve` for a rational quadratic arc over
/// the cache `interval`.
pub(crate) fn analytic_arc(curve: &BSplineCurve, interval: Extent1) -> Result<AnalyticArc, Fault> {
    if curve.degree() != 2 || !curve.is_rational() {
        return Ok(AnalyticArc::Spline);
    }
    let middle = curve.evaluate_point(interval.evaluate(0.456_789))?;
    let scaled_zero = EFF_ZERO * (1.0 + middle.max_dimension());
    let Some(arc) = is_arc(curve, ARC_SAMPLES, scaled_zero)? else {
        return Ok(AnalyticArc::Spline);
    };
    let (mut start, mut end) = (0.0, arc.end_degrees);
    if !curve.natural_interval().is_contained_by(interval, 0.0) {
        // Unreached by `GlobalPointSolve`, whose cache spans the natural
        // interval. SMLib's degrees for a positive angle are a full turn up.
        let low = curve.evaluate_point(interval.min)? - arc.origin;
        let high = curve.evaluate_point(interval.max)? - arc.origin;
        let z = z_axis(arc.x_axis, arc.y_axis);
        let low = ccw_angle_between(z, arc.x_axis, low)?;
        let high = ccw_angle_between(z, arc.x_axis, high)?;
        let degrees = |radians: f64| {
            if radians < 0.0 {
                radians * 180.0 / PI
            } else {
                radians * 180.0 / PI + 360.0
            }
        };
        start = degrees(low);
        end = degrees(high);
        if end < start {
            end += 360.0;
        }
    }
    let mut nurb = curve.clone();
    let knots = unique_knots(curve.knots(), None);
    for &(knot, multiplicity) in knots.iter().take(knots.len().saturating_sub(1)).skip(1) {
        if multiplicity != 2 {
            // SMLib ignores a failed insertion and keeps the spline.
            if let Ok(doubled) = nurb.insert_knot(knot, 1) {
                nurb = doubled;
            }
        }
    }
    let [origin, x_axis, y_axis] = PlaneFrame::canonical(arc.origin, arc.x_axis, arc.y_axis)?;
    let mut circle = Circle {
        origin,
        x_axis,
        y_axis,
        domain: Extent1::new(start, end),
        radius: arc.radius,
        nurb,
        polar: PolarConversion {
            knots: Vec::new(),
            angles: Vec::new(),
            exact: Vec::new(),
        },
    };
    let Some(polar) = circle.set_up_polar_conversion() else {
        return Ok(AnalyticArc::Rebuilt);
    };
    circle.polar = polar;
    for t in [interval.min, interval.max] {
        if curve
            .evaluate_point(t)?
            .distance(circle.nurb.evaluate_point(t)?)
            > scaled_zero
        {
            return Ok(AnalyticArc::Spline);
        }
    }
    Ok(AnalyticArc::Circle(Box::new(circle)))
}

/// `SmCurvePropertyEFO` for `SM_CP_POINT_INVERSION`: `|C - P|^2` and
/// `C' · (C - P)`, evaluated from the right of a knot with nonzero tangents.
struct PointInversion<'a> {
    curve: &'a BSplineCurve,
    target: Vector3,
}

impl EvalFunction for PointInversion<'_> {
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
        let c = self.curve.evaluate(t, 1, true, true)?;
        if c[1].length_squared() < EFF_ZERO {
            return Err(Fault::Degenerate);
        }
        let gap = c[0] - self.target;
        *value = gap.dot(gap);
        *derivative = gap.dot(c[1]);
        Ok(())
    }
}

impl Circle {
    /// The center, unit normal and radius used by analytic intersections.
    pub(crate) fn geometry(&self) -> (Vector3, Vector3, f64) {
        (self.origin, z_axis(self.x_axis, self.y_axis), self.radius)
    }

    /// The inherited B-spline representation of the analytic circle.
    pub(crate) fn spline(&self) -> &BSplineCurve {
        &self.nurb
    }

    /// `STEPInversion` against `domain`: the angle in degrees of `point`,
    /// moved into the domain by whole turns, and the domain's start for the
    /// end of a closed circle.
    fn step_inversion(&self, domain: Extent1, point: Vector3) -> f64 {
        let offset = point - self.origin;
        let x = offset.dot(self.x_axis);
        let y = offset.dot(self.y_axis);
        // SMLib scales `y` by the ratio of the radii, exactly one here.
        let mut angle = libm::atan2(y, x);
        if angle < 0.0 {
            angle += 2.0 * PI;
        }
        let mut degrees = contains_periodic(domain, rad2deg(angle), 360.0, -1.0).value;
        let approximate = ZONE_TOL_3D / 100.0;
        let end = domain.max * PI / 180.0;
        let end = self.origin
            + self.radius * cos(end) * self.x_axis
            + self.radius * sin(end) * self.y_axis;
        if end.distance(point) <= approximate
            && are_same_to_tol(domain.max, degrees, EFF_ZERO_DEG)
            && are_same_to_tol(domain.length(), 360.0, EFF_ZERO_DEG)
        {
            degrees = domain.min;
        }
        degrees
    }

    /// `EvaluateSTEP`'s point: the circle at `degrees`, moved into the
    /// domain by whole turns or clamped to its nearer end.
    fn evaluate_step(&self, degrees: f64) -> Vector3 {
        let found = contains_periodic(self.domain, degrees, 360.0, -1.0);
        let mut degrees = found.value;
        if !found.inside {
            degrees = clamp_periodic(self.domain, degrees, 360.0, -1.0);
        }
        let radians = degrees * PI / 180.0;
        self.origin
            + self.radius * cos(radians) * self.x_axis
            + self.radius * sin(radians) * self.y_axis
    }

    /// `SmPolarConversion::SetUpPolarConversion` for the circle's spline;
    /// may move the domain to the angles of the spline's ends. `None` where
    /// SMLib leaves the conversion impossible.
    fn set_up_polar_conversion(&mut self) -> Option<PolarConversion> {
        let unique = unique_knots(self.nurb.knots(), None);
        let count = unique.len();
        if count < 2 || unique[0].1 != 3 || unique[count - 1].1 != 3 {
            return None;
        }
        if unique[1..count - 1].iter().any(|&(_, m)| m != 2) {
            return None;
        }
        let knots: Vec<f64> = unique.iter().map(|&(k, _)| k).collect();
        let mut domain = self.domain;
        if domain.length() < 360.0 - EFF_ZERO {
            let full = Extent1::new(0.0, 360.0);
            self.domain = full;
            let low = self.nurb.evaluate_point(knots[0]).ok()?;
            let mut lowest = self.step_inversion(full, low);
            let high = self.nurb.evaluate_point(knots[count - 1]).ok()?;
            let mut highest = self.step_inversion(full, high);
            if highest <= lowest + EFF_ZERO {
                highest += 360.0;
            }
            let input_mid = domain.mid();
            let mut this_mid = (lowest + highest) / 2.0;
            while this_mid > input_mid + 180.0 + EFF_ZERO {
                lowest -= 360.0;
                highest -= 360.0;
                this_mid -= 360.0;
            }
            while this_mid < input_mid - 180.0 - EFF_ZERO {
                lowest += 360.0;
                highest += 360.0;
                this_mid += 360.0;
            }
            let low_change = lowest - domain.min;
            if low_change < EFF_ZERO * 1000.0 && low_change > -EFF_ZERO * 1000.0 {
                lowest = domain.min;
            }
            let high_change = highest - domain.max;
            if high_change < EFF_ZERO * 1000.0 && high_change > -EFF_ZERO * 1000.0 {
                highest = domain.max;
            }
            // `SetMinMax` keeps the interval when the pair is reversed.
            let _ = domain.set_min_max(lowest, highest);
            self.domain = domain;
        }
        let weights: Vec<f64> = self
            .nurb
            .control_points()
            .iter()
            .map(|p| p.get(PointForm::Euclidean).1)
            .collect();
        let tol = EFF_ZERO_DEG;
        let mut angles = Vec::with_capacity(count);
        let mut exact = Vec::with_capacity(count - 1);
        let mut span_end = 0.0;
        for (ii, &knot) in knots.iter().enumerate() {
            let point = self.nurb.evaluate_point(knot).ok()?;
            let mut angle = self.step_inversion(self.domain, point);
            let mut end_point = false;
            if ii == 0 && ((angle - domain.min).abs() < tol || (angle - domain.max).abs() < tol) {
                angle = domain.min;
                end_point = true;
            }
            if ii == count - 1
                && ((angle - domain.max).abs() < tol || (angle - domain.min).abs() < tol)
            {
                angle = domain.max;
                end_point = true;
            }
            if !end_point {
                if (angle + 360.0).abs() < tol {
                    angle = -360.0;
                } else if angle.abs() < tol {
                    angle = 0.0;
                } else if (angle - 360.0).abs() < tol {
                    angle = 360.0;
                }
            }
            let span_start = span_end;
            span_end = angle;
            angles.push(angle);
            if ii == 0 {
                continue;
            }
            let wa = *weights.get(ii * 2 - 2)?;
            let wb = *weights.get(ii * 2 - 1)?;
            let wc = *weights.get(ii * 2)?;
            let cos_half = cos((span_end - span_start) / 2.0 * PI / 180.0);
            exact.push(are_same(wa, wc) && are_same(wb / wc, cos_half));
        }
        self.domain = Extent1::new(angles[0], angles[count - 1]);
        Some(PolarConversion {
            knots,
            angles,
            exact,
        })
    }

    /// `SmPolarConversion::AngleParamToNurbParam`: the spline parameter of
    /// `degrees` within one conic span.
    fn angle_to_nurb(
        span: Extent1,
        start_degrees: f64,
        end_degrees: f64,
        degrees: f64,
    ) -> Result<f64, Fault> {
        let theta = (end_degrees - start_degrees) / 2.0 * PI / 180.0;
        let tol = EFF_ZERO_DEG;
        let range = Extent1::new(
            start_degrees.min(end_degrees),
            start_degrees.max(end_degrees),
        );
        let inside = range.contains_value(degrees, 0.0);
        let degrees = range.clamp_value(degrees);
        if !inside {
            return Err(Fault::Degenerate);
        }
        if (degrees - start_degrees).abs() < tol {
            return Ok(span.min);
        }
        if (degrees - end_degrees).abs() < tol {
            return Ok(span.max);
        }
        let angle = (degrees - start_degrees) * PI / 180.0;
        let cos_angle = cos(angle);
        let cos_theta = cos(theta);
        let a = 2.0 * cos_angle * (1.0 - cos_theta) - 1.0 + 2.0 * cos_theta - cos(2.0 * theta);
        let b = 2.0 * cos_angle * (-1.0 + cos_theta) + 2.0 - 2.0 * cos_theta;
        let c = cos_angle - 1.0;
        let discriminant = b * b - 4.0 * a * c;
        if discriminant < 0.0 {
            return Err(Fault::Degenerate);
        }
        let root = discriminant.sqrt();
        let mut t = (-b + root) / (2.0 * a);
        if !(0.0..=1.0).contains(&t) {
            t = (-b - root) / (2.0 * a);
            if !(0.0..=1.0).contains(&t) {
                return Err(Fault::Degenerate);
            }
        }
        Ok((span.max - span.min) * t + span.min)
    }

    /// `SmPolarConversion::ConvertToNURBSParameter`: the spline parameter of
    /// `degrees`, from the conic span that holds it, kept when its point is
    /// on the circle's point at that angle and otherwise refined by point
    /// inversion on the spline.
    fn convert_to_nurbs(&self, degrees: f64) -> Result<f64, Fault> {
        let polar = &self.polar;
        let (first, last) = (polar.knots[0], polar.knots[polar.knots.len() - 1]);
        let tol = EFF_ZERO_DEG;
        let mut degrees = degrees;
        let mut guess = first
            + (degrees - self.domain.min) / (self.domain.max - self.domain.min) * (last - first);
        let range = Extent1::new(polar.angles[0], polar.angles[polar.angles.len() - 1]);
        if !range.contains_value(degrees, tol) {
            degrees = range.clamp_value(degrees);
        }
        if are_same_to_tol(degrees, polar.angles[0], tol) {
            guess = first;
        } else {
            for ii in 1..polar.angles.len() {
                if are_same_to_tol(degrees, polar.angles[ii], tol) {
                    guess = polar.knots[ii];
                    break;
                }
                if degrees <= polar.angles[ii] {
                    guess = Self::angle_to_nurb(
                        Extent1::new(polar.knots[ii - 1], polar.knots[ii]),
                        polar.angles[ii - 1],
                        polar.angles[ii],
                        degrees,
                    )?;
                    break;
                }
            }
        }
        let target = self.evaluate_step(degrees);
        let point = self.nurb.evaluate_point(guess)?;
        if (point - target).length() < EFF_ZERO * (1.0 + point.max_dimension()) {
            return Ok(guess);
        }
        let natural = self.nurb.natural_interval();
        let mut solver = LocalSolve1d::with_interval(
            PointInversion {
                curve: &self.nurb,
                target,
            },
            natural,
            false,
        );
        solver.set_desired_accuracy(EFF_ZERO_SQ);
        let root = solver.solve_it(guess, EFF_ZERO_SQRT)?;
        if !root.found {
            return Err(Fault::Degenerate);
        }
        Ok(root.t)
    }

    /// `smgu_TransformPointToStartPlane` through
    /// `SmEllipse::TransformPointToStartPlane`, angles only: none for a
    /// point on the axis, two for a point on the seam of a whole circle
    /// (one exact, the far one snapped), else one.
    fn angles_of(&self, test: Vector3, tolerance: f64) -> Result<Vec<f64>, Fault> {
        let z = z_axis(self.x_axis, self.y_axis);
        let (to_axis, _) = line_point_distance(self.origin, z, test);
        let scaled_zero = EFF_ZERO * (1.0 + test.max_dimension());
        if to_axis < scaled_zero {
            return Ok(Vec::new());
        }
        let sweep = self.domain;
        let angle_tol = ANALYTIC_TOL_SCALE * EFF_ZERO * (1.0 + 360.0);
        let mut radians = ccw_angle_between(z, self.x_axis, test - self.origin)?;
        if radians < EFF_ZERO {
            if radians > -EFF_ZERO {
                radians = 0.0;
            } else {
                radians += 2.0 * PI;
            }
        }
        let found = contains_periodic(sweep, rad2deg(radians), 360.0, -1.0);
        let mut degrees = found.value;
        if !found.inside && degrees < 0.0 {
            degrees += 360.0;
        }
        let boundary_tol = rad2deg(tolerance / to_axis);
        if are_same_to_tol(sweep.length(), 360.0, angle_tol)
            && is_on_periodic_boundary(sweep, degrees, 360.0, boundary_tol)
        {
            if (sweep.min - degrees).abs() < (sweep.max - degrees).abs() {
                return Ok(vec![degrees, sweep.max]);
            }
            return Ok(vec![sweep.min, degrees]);
        }
        Ok(vec![degrees])
    }

    /// `SmEllipse::DropPointFast` for minimising (`meeting` false) and
    /// meeting within `tolerance`, with no target distance and no pointing
    /// vector, as `GlobalPointSolve` calls it: the spline parameters, within
    /// the search `interval`, of the circle's points at the test point's
    /// angle, or at both domain ends for a point on the axis. `single`
    /// keeps the nearer of two.
    pub(crate) fn drop_point_fast(
        &self,
        interval: Extent1,
        meeting: bool,
        test: Vector3,
        tolerance: f64,
        single: bool,
    ) -> CircleDrop {
        let mut answers = Vec::new();
        let Ok(mut angles) = self.angles_of(test, tolerance) else {
            return CircleDrop {
                answers,
                complete: false,
            };
        };
        if angles.is_empty() {
            angles = vec![self.domain.min, self.domain.max];
        }
        let lone = angles.len() == 1;
        let mut distances = [0.0; 2];
        for (ii, &angle) in angles.iter().enumerate() {
            let clamped = if lone {
                clamp_periodic(self.domain, angle, 360.0, EFF_ZERO_DEG)
            } else {
                self.domain.clamp_value(angle)
            };
            // The conversion reports every success as exact, so SMLib's
            // inexact branch never runs.
            let Ok(parameter) = self.convert_to_nurbs(clamped) else {
                return CircleDrop {
                    answers,
                    complete: false,
                };
            };
            let parameter = interval.clamp_value(parameter);
            let Ok(point) = self.nurb.evaluate_point(parameter) else {
                return CircleDrop {
                    answers,
                    complete: false,
                };
            };
            distances[ii] = point.distance(test);
            if meeting && distances[ii] > tolerance {
                continue;
            }
            answers.push((parameter, distances[ii]));
        }
        if single && answers.len() > 1 {
            if distances[1] < distances[0] {
                answers.remove(0);
            } else {
                answers.truncate(1);
            }
        }
        CircleDrop {
            answers,
            complete: true,
        }
    }
}
