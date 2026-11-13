// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: a line owns its two-point curve directly instead of
// inheriting from `SmBSplineCurve`, contexts, owner surfaces and
// notifications are dropped, statuses become `Fault`, and output references
// become returned values. Only the bounded two-point forms are built.
// Translated from SMLib SmLine.cpp (the bounded constructor,
// CreateLineSegment, ConvertTFromNURBSToSTEP, ConvertTFromSTEPToNURBS,
// DropPointFast for minimising and meeting, Evaluate, EvaluateSTEP,
// EvaluatePoint, IsBounded, IsDegenerate, GetNaturalInterval,
// CalculateBoundingBox, MakeNurb) and SmBSplineCurve.cpp (EditParameterization
// through N_BasisReparam) at c9c5979f04cea8c38213004d4f04df0dcfb40d05;
// rust-brep has no port of these.
//! SMLib's `SmLine`: a straight curve with an analytic (STEP) parameter and a
//! NURBS parameter that may differ.
//!
//! The line is `point + s * scale * direction` for the STEP parameter `s` in
//! the analytic domain. Its NURBS curve has degree one, two control points
//! and knots at the ends of the analytic domain until the parameter range
//! is edited. A line of no length keeps a collapsed domain, which
//! [`BSplineCurve`] cannot hold; [`Line::make_nurb`] refuses it.

use crate::Fault;
use crate::extent::{Extent1, Extent3};
use crate::geometry::line_point_distance;
use crate::pseudo_box::PseudoBox;
use crate::scalar::{EFF_ZERO, are_same, is_zero, max};
use crate::spline::{BSplineCurve, ControlPoint, NO_WEIGHT, NO_Z};
use crate::vector::Vector3;

/// A straight curve with analytic and NURBS parameters (`SmLine`).
#[derive(Clone, Debug, PartialEq)]
pub struct Line {
    point: Vector3,
    direction: Vector3,
    domain: Extent1,
    scale: f64,
    inside_out: bool,
    planar: bool,
    knots: [f64; 4],
    ends: [Vector3; 2],
}

/// The most derivatives the evaluators return, as for [`BSplineCurve`].
pub const MAX_DERIVATIVES: usize = 4;

impl Line {
    /// The analytic origin and scaled STEP direction used by intersections.
    pub(crate) fn geometry(&self) -> (Vector3, Vector3) {
        (self.point, self.direction * self.scale)
    }

    /// The bounded `SmLine` constructor: `point + s * scale * direction` for
    /// `s` in `domain`. `direction` is unitised; a direction shorter than
    /// `SM_EFF_ZERO` becomes the x axis and collapses the domain to its
    /// start. `planar` makes the NURBS control points two-dimensional.
    /// `inside_out` runs the NURBS parameter against the STEP parameter, as
    /// a cone's generator does.
    pub fn new(
        point: Vector3,
        direction: Vector3,
        domain: Extent1,
        scale: f64,
        planar: bool,
        inside_out: bool,
    ) -> Self {
        let length = direction.length();
        let (direction, domain) = if length > EFF_ZERO {
            (direction.divided(length), domain)
        } else {
            (
                Vector3::new(1.0, 0.0, 0.0),
                Extent1::new(domain.min, domain.min),
            )
        };
        let mut line = Self {
            point,
            direction,
            domain,
            scale,
            inside_out,
            planar,
            knots: [0.0; 4],
            ends: [Vector3::default(); 2],
        };
        line.build_nurb();
        line
    }

    /// `CreateLineSegment`: the segment from `start` to `end` over the
    /// analytic domain `[0, 1]`, with its NURBS parameter set to `range`
    /// when given. Coincident ends are [`Fault::InvalidInput`].
    pub fn segment(
        start: Vector3,
        end: Vector3,
        planar: bool,
        range: Option<Extent1>,
    ) -> Result<Self, Fault> {
        let mut direction = end - start;
        let scale = direction.length();
        direction.unitize()?;
        let mut line = Self::new(
            start,
            direction,
            Extent1::new(0.0, 1.0),
            scale,
            planar,
            false,
        );
        if let Some(range) = range {
            line.edit_parameterization(range);
        }
        Ok(line)
    }

    /// `MakeNurb`: the two NURBS end points and the knots at the analytic
    /// domain's ends. An inside-out line starts at its far end. A planar
    /// line is then projected to the xy plane, by `ConvertTo2D`: its point
    /// and direction lose their z, and the direction is not unitised again.
    fn build_nurb(&mut self) {
        let (first, last) = if self.inside_out {
            (self.domain.max, self.domain.min)
        } else {
            (self.domain.min, self.domain.max)
        };
        self.ends = [
            self.evaluate_step_point(first),
            self.evaluate_step_point(last),
        ];
        self.knots = [
            self.domain.min,
            self.domain.min,
            self.domain.max,
            self.domain.max,
        ];
        if self.planar {
            self.point.z = 0.0;
            self.direction.z = 0.0;
        }
    }

    /// `EvaluateSTEP` at no derivatives: `point + s * direction * scale`.
    fn evaluate_step_point(&self, s: f64) -> Vector3 {
        self.point + s * self.direction * self.scale
    }

    /// The knots of the NURBS curve.
    pub fn nurb_knots(&self) -> [f64; 4] {
        self.knots
    }

    /// The control points of the NURBS curve. A planar line has [`NO_Z`]
    /// for every `z`, and no weights.
    pub fn nurb_points(&self) -> [ControlPoint; 2] {
        self.ends.map(|end| ControlPoint {
            coordinates: [
                end.x,
                end.y,
                if self.planar { NO_Z } else { end.z },
                NO_WEIGHT,
            ],
        })
    }

    /// `MakeNurb`: the NURBS curve as a [`BSplineCurve`]. A line of no
    /// length has no valid knot vector and is [`Fault::InvalidInput`].
    pub fn make_nurb(&self) -> Result<BSplineCurve, Fault> {
        BSplineCurve::new(1, self.knots.to_vec(), self.nurb_points().to_vec())
    }

    /// `GetNaturalInterval`: the range of the NURBS parameter.
    pub fn natural_interval(&self) -> Extent1 {
        Extent1::new(self.knots[0], self.knots[3])
    }

    /// `IsBounded`: neither end of the analytic domain is unbounded.
    pub fn is_bounded(&self) -> bool {
        self.domain.is_bounded()
    }

    /// `EditParameterization`: the same line over a new NURBS parameter
    /// range. A range within the scaled zero of the current one changes
    /// nothing. The analytic domain is kept.
    pub fn edit_parameterization(&mut self, range: Extent1) {
        let old = self.natural_interval();
        if are_same(old.min, range.min) && are_same(old.max, range.max) {
            return;
        }
        // `N_BasisReparam` on `[a, a, b, b]` sets both ends and no interior knot.
        self.knots = [range.min, range.min, range.max, range.max];
    }

    /// `ConvertTFromNURBSToSTEP`: the analytic parameter of a NURBS parameter.
    pub fn nurbs_to_step(&self, t: f64) -> f64 {
        let natural = self.natural_interval();
        if self.is_bounded() {
            if is_zero(natural.length()) {
                t
            } else {
                let fraction = if self.inside_out {
                    (t - natural.max) / -natural.length() * self.domain.length()
                } else {
                    (t - natural.min) / natural.length() * self.domain.length()
                };
                self.domain.min + fraction
            }
        } else if self.inside_out {
            -t
        } else {
            t
        }
    }

    /// `ConvertTFromSTEPToNURBS`: the NURBS parameter of an analytic one.
    pub fn step_to_nurbs(&self, s: f64) -> f64 {
        let natural = self.natural_interval();
        if self.is_bounded() {
            // A collapsed analytic domain, or a collapsed NURBS parameter
            // range, has one parameter for every analytic one. Dividing by
            // its zero length would turn the candidate into NaN.
            if is_zero(self.domain.length()) || is_zero(natural.length()) {
                natural.min
            } else if self.inside_out {
                natural.max + (s - self.domain.min) / self.domain.length() * -natural.length()
            } else {
                natural.min + (s - self.domain.min) / self.domain.length() * natural.length()
            }
        } else if self.inside_out {
            -s
        } else {
            s
        }
    }

    /// `DropPointFast` for minimising and meeting (others are `None`, SMLib's
    /// `SM_ERR`): the foot of `test` on the line, clamped to the NURBS
    /// `interval`, as the NURBS parameter and its distance. A meeting farther
    /// than `tolerance` has no answer (`Some(None)`).
    /// A collapsed parameter map measures from the retained curve point,
    /// not from the infinite support line.
    pub fn drop_point_fast(
        &self,
        interval: Extent1,
        meeting: bool,
        test: Vector3,
        tolerance: f64,
    ) -> Option<(f64, f64)> {
        let (mut distance, step) = line_point_distance(self.point, self.direction, test);
        let raw = self.step_to_nurbs(step / self.scale);
        let parameter = interval.clamp_value(raw);
        if parameter != raw
            || (self.is_bounded()
                && (is_zero(self.domain.length()) || is_zero(self.natural_interval().length())))
        {
            distance = self.evaluate_point(parameter).distance(test);
        }
        (!meeting || distance <= tolerance).then_some((parameter, distance))
    }

    /// `Evaluate`: the point at NURBS parameter `t` and `derivatives`
    /// derivatives. Every derivative above the first is zero. More than
    /// [`MAX_DERIVATIVES`] is [`Fault::InvalidInput`].
    pub fn evaluate(&self, t: f64, derivatives: usize) -> Result<Vec<Vector3>, Fault> {
        if derivatives > MAX_DERIVATIVES {
            return Err(Fault::InvalidInput);
        }
        let step = self.nurbs_to_step(t);
        let mut values = vec![Vector3::new(0.0, 0.0, 0.0); derivatives + 1];
        values[0] = self.evaluate_step_point(step);
        if derivatives >= 1 {
            let jacobian = self.domain.length() / self.natural_interval().length();
            values[1] = if self.inside_out {
                -self.direction * self.scale * jacobian
            } else {
                self.direction * self.scale * jacobian
            };
        }
        Ok(values)
    }

    /// `EvaluateSTEP`: the point at analytic parameter `s` and `derivatives`
    /// derivatives, which are the constant velocity and zeros. More than
    /// [`MAX_DERIVATIVES`] is [`Fault::InvalidInput`].
    pub fn evaluate_step(&self, s: f64, derivatives: usize) -> Result<Vec<Vector3>, Fault> {
        if derivatives > MAX_DERIVATIVES {
            return Err(Fault::InvalidInput);
        }
        let mut values = vec![Vector3::new(0.0, 0.0, 0.0); derivatives + 1];
        values[0] = self.evaluate_step_point(s);
        if derivatives >= 1 {
            values[1] = self.direction * self.scale;
        }
        Ok(values)
    }

    /// `EvaluatePoint`: the point at NURBS parameter `t`.
    pub fn evaluate_point(&self, t: f64) -> Vector3 {
        let step = self.nurbs_to_step(t);
        self.point + step * self.scale * self.direction
    }

    /// `IsDegenerate`: whether the line, or the part over a NURBS
    /// `interval`, is shorter than `tolerance` or `SM_EFF_ZERO`.
    pub fn is_degenerate(&self, tolerance: f64, interval: Option<Extent1>) -> bool {
        let tolerance = max(tolerance, EFF_ZERO);
        let distance = match interval {
            Some(interval) => {
                let low = self.nurbs_to_step(interval.min);
                let high = self.nurbs_to_step(interval.max);
                self.scale * (high - low).abs()
            }
            None => self.scale * self.domain.length(),
        };
        distance < tolerance
    }

    /// `CalculateBoundingBox`: the axis-aligned box and the pseudo box of
    /// the part over a NURBS `interval`, both through its end points. The
    /// polar box, which SMLib also offers, is not ported.
    pub fn bounding_boxes(&self, interval: Extent1) -> (Extent3, PseudoBox) {
        let start = self.evaluate_point(interval.min);
        let end = self.evaluate_point(interval.max);
        let mut pseudo = PseudoBox::new();
        let across = end - start;
        if !is_zero(across.length_squared())
            && let Ok([x, y, z]) = across.make_unit_ortho_vectors(None)
        {
            let _ = pseudo.set_basis(x, y, z);
        }
        let mut bbox = Extent3::singleton(start);
        bbox.add_point(end);
        pseudo.add_point(start);
        pseudo.add_point(end);
        (bbox, pseudo)
    }
}

impl crate::curve_solve::PointSolveCurve for Line {
    fn evaluate(
        &self,
        t: f64,
        derivatives: usize,
        _from_left: bool,
        _nonzero_tangent: bool,
    ) -> Result<Vec<Vector3>, Fault> {
        Line::evaluate(self, t, derivatives)
    }
    fn evaluate_point(&self, t: f64) -> Result<Vector3, Fault> {
        Ok(Line::evaluate_point(self, t))
    }
    fn is_degenerate(&self, tolerance: f64) -> bool {
        Line::is_degenerate(self, tolerance, None)
    }
}

impl crate::curve_solve::SolverCurve for Line {
    fn natural_interval(&self) -> Extent1 {
        Line::natural_interval(self)
    }
    fn degree(&self) -> usize {
        1
    }
    fn natural_knot_count(&self) -> usize {
        4
    }
    fn knot_values(&self) -> Vec<f64> {
        vec![self.knots[0], self.knots[3]]
    }
}
