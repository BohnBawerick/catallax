// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: the intersector borrows a `RollingBall` instead of
// SMLib's fillet solver and offset-surface objects, statuses become `Fault`
// or `StepRefusal`, output references become returned values, and
// `SetSurfacesToZeroOffset` becomes evaluation of zero-offset copies of the
// base surfaces. Translated directly from SMLib SmFilletIntersector.cpp
// (the constructor, `ComputePointValues`, `ComputePointValuesDelta`,
// `TestSpanAccuracy`, `ComputeNextPoint`, `EvaluateLawPoint`) and
// SmCrvOnSurf.cpp (`Evaluate` to one derivative); rust-brep has no port.
//! SMLib's `SmFilletIntersector`: the intersector of two offset surfaces
//! that the rolling-ball fillet traces its spine with. Its stepping methods
//! are here, and `DoPointIntersection` (the trace through given points)
//! running the `surface_trace.rs` tracer over them, which takes any stepper
//! that implements `TraceSteps` where SMLib calls virtual methods.
//!
//! `FlushCurve`, which builds the fillet surface and its rails from the
//! traced points, is not carried here; [`FilletTrace`] holds what
//! `FlushCurve` would be given.
//!
//! One refusal stands where SMLib goes on: a step that a domain clips is
//! re-solved by SMLib with `LocalCurveIntersect` on an offset surface
//! ([`StepRefusal::ClippedStep`]), and the port refuses it. An evaluation
//! outside a base surface's domain extends the surface, as SMLib's does.
//!
//! Everything returned is an unproved candidate.

use crate::{
    Fault,
    curve_solve::local_point_solve,
    curve_surface::drop_vector,
    extent::Extent2,
    geometry::line_closest_point,
    global_solver::{SolverOperation, SolverRefusal},
    hermite::Hermite,
    rolling_ball::{RollingBall, StepRefusal},
    scalar::{EFF_ZERO, EFF_ZERO_SQ, PI, max, sqrt},
    surface_curves::IntersectionCurve,
    surface_offset::OffsetSurface,
    surface_trace::{PointKind, StepOutcome, TraceState, TraceSteps, Tracer, TsectPoint},
    vector::{Vector2, Vector3},
};
use std::cell::Cell;

/// Where `TestSpanAccuracy` and `ComputeNextPoint` test a span.
const TEST_PARAMETERS: [f64; 5] = [0.5, 0.17, 0.83, 0.33, 0.67];

/// SMLib's `SmFilletIntersector` over a rolling ball.
#[derive(Clone, Copy, Debug)]
pub struct FilletIntersector<'a> {
    ball: &'a RollingBall<'a>,
    /// `m_vUVDomain`: the domains as given, which SMLib lets exceed the
    /// natural ones; the port's offset surfaces clamp to the natural ones.
    pub domains: [Extent2; 2],
    /// `m_dThisAngTolRad`.
    pub angle_tolerance: f64,
    /// `m_dThisApproxTol3d`.
    pub approx_tolerance: f64,
    /// `m_bClippedByThroughPoint`.
    pub clipped_by_through_point: bool,
    /// `m_dCurveTraceDirection`.
    pub trace_direction: f64,
    /// `m_dExtensionDistance`: how far past the last through point a trace
    /// may extend, in curve parameter, when it is allowed to extend.
    pub extension_distance: f64,
}

/// What `DoPointIntersection` leaves for `FlushCurve`: the traced points,
/// and the through points as the trace left them.
#[derive(Clone, Debug, PartialEq)]
pub struct FilletTrace {
    /// Whether `TraceCurve` succeeded; SMLib flushes nothing, and still
    /// succeeds, when it does not.
    pub traced: bool,
    /// `m_vCurvePoints` at the call of `FlushCurve`.
    pub points: Vec<TsectPoint>,
    /// The caller's through points: `m_ePointType` and `m_dCurveParameter`
    /// updated, the parameter still `-1` for a point the trace did not reach.
    pub through_points: Vec<TsectPoint>,
    /// `m_bCurveIsClosed`.
    pub curve_is_closed: bool,
}

/// What `DoPointIntersection` is asked: SMLib's optional arguments.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TraceRequest {
    /// `bExtendBeforeStart`.
    pub extend_before_start: bool,
    /// `bExtendAfterEnd`.
    pub extend_after_end: bool,
    /// `pOptStartDirection`.
    pub start_direction: Option<Vector3>,
    /// `pOptEndDirection`.
    pub end_direction: Option<Vector3>,
    /// `pOptApproxTol3d`.
    pub approx_tolerance: Option<f64>,
    /// `pdOptAngTolRad`.
    pub angle_tolerance: Option<f64>,
}

/// What `TestSpanAccuracy` found: SMLib's `rbFoundGoodPoint`,
/// `rdAngleFoundRad` and whether it wrote them. The angle is left at its
/// input value when the next point has no tangent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpanTest {
    /// The span met the angle and distance tolerances.
    pub good: bool,
    /// The turn of the curve tangent over the span.
    pub angle: f64,
}

/// `SmCrvOnSurf::Evaluate` of a parameter-space curve on a surface to one
/// derivative, from the left: the point, and the derivative lifted into
/// space. SMLib falls back to no derivatives when the surface cannot give
/// one and then returns the tangent unset, which a caller reads
/// uninitialised; the port refuses that.
fn lift(surface: &OffsetSurface, uv_curve: &Hermite, t: f64) -> Result<[Vector3; 2], StepRefusal> {
    let uv = uv_curve.evaluate(t, 1, true)?;
    let at = Vector2::new(uv[0].x, uv[0].y);
    let cells = surface.evaluate(at, 1, 1, [true, true], true, true)?;
    Ok([cells[0], uv[1].x * cells[2] + uv[1].y * cells[1]])
}

impl<'a> FilletIntersector<'a> {
    /// `SmFilletIntersector::SmFilletIntersector`: the domains are used as
    /// given, the angle tolerance is twenty degrees and the 3D tolerance a
    /// thousandth of the first offset surface's size over its domain.
    pub fn new(ball: &'a RollingBall<'a>, domain_1: Extent2, domain_2: Extent2) -> Self {
        let surface = ball.surface(0);
        let p1 = surface.evaluate_point(domain_1.min).unwrap_or_default();
        let middle = surface
            .evaluate_point(domain_1.evaluate(0.5, 0.5))
            .unwrap_or_default();
        let p2 = surface.evaluate_point(domain_1.max).unwrap_or_default();
        let size = p1.distance(middle) + middle.distance(p2);
        Self {
            ball,
            domains: [domain_1, domain_2],
            angle_tolerance: 20.0 * PI / 180.0,
            approx_tolerance: size / 1000.0,
            clipped_by_through_point: false,
            trace_direction: 1.0,
            extension_distance: 0.0,
        }
    }

    /// `ComputePointValues`: the solver's, with this intersector's domains
    /// and trace direction. Any failure is `SM_ERR`, whatever the solver's
    /// status was.
    pub fn compute_point_values(
        &self,
        uv: [Vector2; 2],
        point: &mut TsectPoint,
        previous: Option<&TsectPoint>,
    ) -> Result<(), StepRefusal> {
        self.ball
            .compute_point_values(uv, self.domains, self.trace_direction, point, previous)
            .map_err(|refusal| match refusal {
                StepRefusal::Failed(_) => StepRefusal::Failed(Fault::Degenerate),
                other => other,
            })
    }

    /// `ComputePointValuesDelta`: the surface values at `uv` into `point` and
    /// at `uv_delta` into `delta`, the unit chord between the two curve
    /// points as `point`'s tangent (negated when the delta has the smaller
    /// curve parameter), `point` a crossing, and its parameter tangents.
    pub fn compute_point_values_delta(
        &self,
        uv: [Vector2; 2],
        uv_delta: [Vector2; 2],
        point: &mut TsectPoint,
        delta: &mut TsectPoint,
    ) -> Result<(), Fault> {
        self.ball.compute_surface_values(uv, point)?;
        self.ball.compute_surface_values(uv_delta, delta)?;
        let mut chord = delta.curve_point - point.curve_point;
        if delta.curve_parameter < point.curve_parameter {
            chord = -chord;
        }
        if chord.length() < EFF_ZERO {
            return Err(Fault::Degenerate);
        }
        chord.unitize()?;
        point.curve_tangent = chord;
        point.kind = PointKind::Crossing;
        self.ball.compute_uv_vectors(point)
    }

    /// `EvaluateLawPoint` for a constant radius: the surface parameters on
    /// the plane through `origin` with `normal`, from the guesses. SMLib
    /// takes the curve parameter from the first user double of the solver's
    /// point, which a constant-radius solve (four unknowns) never sets: it is
    /// the point's zero. `None` when no solution is found, where SMLib returns
    /// the unset parameters.
    pub fn evaluate_law_point(
        &self,
        origin: Vector3,
        normal: Vector3,
        guesses: [Vector2; 2],
    ) -> Result<Option<([Vector2; 2], f64)>, StepRefusal> {
        Ok(self
            .ball
            .point_on_plane_solve(origin, normal, self.domains, guesses)?
            .map(|uvs| (uvs, 0.0)))
    }

    /// The Hermite curves of one step: the parameter-space curves through
    /// both ends' parameters and tangents, and the 3D rail curves through
    /// their lifts onto the base surfaces. `adjusted` takes the 3D tangents
    /// from points a thousandth in from each end.
    fn rail_curves(
        &self,
        start: &TsectPoint,
        next: &TsectPoint,
        step: f64,
        adjusted: bool,
    ) -> Result<([Hermite; 2], [Hermite; 2]), StepRefusal> {
        let uv_curves = [0, 1].map(|s| {
            Hermite::planar(
                start.uv[s],
                start.uv_deriv[s] * step,
                next.uv[s],
                next.uv_deriv[s] * step,
            )
        });
        let zero = self.ball.zero_offset_surfaces();
        let mut rails = Vec::with_capacity(2);
        for s in 0..2 {
            let [p0, mut d0] = lift(&zero[s], &uv_curves[s], 0.0)?;
            let [p1, mut d1] = lift(&zero[s], &uv_curves[s], 1.0)?;
            if adjusted {
                d0 = lift(&zero[s], &uv_curves[s], 0.001)?[1];
                d1 = lift(&zero[s], &uv_curves[s], 0.999)?[1];
            }
            rails.push(Hermite::new(p0, d0, p1, d1));
        }
        Ok((uv_curves, [rails[0], rails[1]]))
    }

    /// The five-point test the span tests share: the largest squared
    /// distance between a base-surface point on the parameter curves and the
    /// matching rail point, or `None` when one exceeds `tolerance_sq` on a
    /// step long enough to reject.
    fn test_points(
        &self,
        uv_curves: &[Hermite; 2],
        rails: &[Hermite; 2],
        step: f64,
        tolerance_sq: f64,
    ) -> Result<Option<f64>, StepRefusal> {
        let zero = self.ball.zero_offset_surfaces();
        let mut max_sq: f64 = 0.0;
        for t in TEST_PARAMETERS {
            let mut distances_sq = [0.0; 2];
            for s in 0..2 {
                let q = uv_curves[s].evaluate_point(t);
                let uv = self.domains[s].clamp_point(Vector2::new(q.x, q.y));
                let on_surface = zero[s].evaluate_point(uv)?;
                distances_sq[s] = on_surface.distance_squared(rails[s].evaluate_point(t));
            }
            let distance_sq = max(distances_sq[0], distances_sq[1]);
            if distance_sq > tolerance_sq && step > EFF_ZERO * 100.0 {
                return Ok(None);
            }
            if distance_sq > max_sq {
                max_sq = distance_sq;
            }
        }
        Ok(Some(max_sq))
    }

    /// `TestSpanAccuracy`: the span from `start` to `next` is good when the
    /// tangent turns less than the angle tolerance and, at five points of the
    /// Hermite curves through both ends (and the optional `middle` point's
    /// nearness to the rails), the base surfaces and the rails stay within
    /// the tolerance. A failing span demotes a singular `next` to crossing.
    pub fn test_span_accuracy(
        &self,
        start: &TsectPoint,
        next: &mut TsectPoint,
        middle: Option<&TsectPoint>,
    ) -> Result<SpanTest, StepRefusal> {
        let step = next.curve_parameter - start.curve_parameter;
        if next.curve_tangent.length_squared() < EFF_ZERO_SQ {
            return Err(Fault::Degenerate.into());
        }
        let angle = start.curve_tangent.angle_between(next.curve_tangent)?;
        let demote = |next: &mut TsectPoint| {
            if next.kind == PointKind::Singularity {
                next.kind = PointKind::Crossing;
            }
        };
        let fail = |next: &mut TsectPoint| {
            demote(next);
            Ok(SpanTest { good: false, angle })
        };
        if angle > self.angle_tolerance && step > EFF_ZERO * 100.0 {
            return fail(next);
        }
        let (uv_curves, rails) = self.rail_curves(start, next, step, false)?;
        // SMLib's adjusted tolerance never leaves the tolerance itself.
        let tolerance = self.approx_tolerance;
        if let Some(middle) = middle {
            let zero = self.ball.zero_offset_surfaces();
            let big_sq = (tolerance * 100.0) * (tolerance * 100.0);
            for s in 0..2 {
                let test = zero[s].evaluate_point(middle.uv[s])?;
                let found = local_point_solve(
                    &rails[s],
                    rails[s].natural_interval(),
                    SolverOperation::Minimize,
                    test,
                    0.5,
                )?;
                let at = found.map_or(0.5, |hit| hit.parameter);
                if rails[s].evaluate_point(at).distance_squared(test) > big_sq {
                    return Ok(SpanTest { good: false, angle });
                }
            }
        }
        let tolerance_sq = tolerance * tolerance;
        match self.test_points(&uv_curves, &rails, step, tolerance_sq)? {
            Some(_) => Ok(SpanTest { good: true, angle }),
            None => fail(next),
        }
    }

    /// `ComputeNextPoint`: steps `step` from `current` with the solver,
    /// computes the next point's values and tests the span. A through-point
    /// trace may change `current`'s tangents, see [`Self::adjust_derivatives`].
    /// A clipped step is refused, see [`StepRefusal::ClippedStep`].
    pub fn compute_next_point(
        &self,
        current: &mut TsectPoint,
        step: f64,
        next: &mut TsectPoint,
    ) -> Result<StepOutcome, StepRefusal> {
        let mut outcome = StepOutcome::default();
        if current.kind != PointKind::Crossing {
            return Ok(outcome);
        }
        let stepped = match self
            .ball
            .step_solve(current, step, self.domains[0], self.domains[1])
        {
            Ok(stepped) => stepped,
            // SMLib reads any other failure as no next point.
            Err(StepRefusal::Failed(_)) => return Ok(outcome),
            Err(refusal) => return Err(refusal),
        };
        outcome.clipped = stepped.clipped;
        outcome.boundary_hit = stepped.boundary_hit;
        if !stepped.found {
            return Ok(outcome);
        }
        next.uv = stepped.next.uv;
        let uvs = next.uv;
        match self.compute_point_values(uvs, next, Some(current)) {
            Ok(()) => {}
            Err(StepRefusal::Failed(_)) => {
                // The derivatives fail leaving a surface boundary: the trace
                // is done.
                outcome.boundary_hit = true;
                return Ok(outcome);
            }
            Err(refusal) => return Err(refusal),
        }
        if next.kind == PointKind::TangentPoint || next.kind == PointKind::Coincidence {
            return Ok(outcome);
        }
        if next.curve_tangent.length_squared() < EFF_ZERO_SQ {
            return Ok(outcome);
        }
        let angle = current.curve_tangent.angle_between(next.curve_tangent)?;
        outcome.angle = angle;
        if angle > self.angle_tolerance && step > EFF_ZERO * 1000.0 {
            // Split the span unless the two points coincide.
            let tolerance = EFF_ZERO * 1000.0 * (1.0 + next.curve_point.max_dimension());
            if current.curve_point.distance(next.curve_point) > tolerance {
                if next.kind == PointKind::Singularity {
                    next.kind = PointKind::Crossing;
                }
                return Ok(outcome);
            }
        }
        // The next point must lie ahead of the current one.
        let ahead = line_closest_point(
            Vector3::new(current.uv[0].x, current.uv[0].y, 0.0),
            Vector3::new(current.uv_deriv[0].x, current.uv_deriv[0].y, 0.0),
            Vector3::new(next.uv[0].x, next.uv[0].y, 0.0),
        );
        if ahead < 0.0 {
            if next.kind == PointKind::Singularity {
                next.kind = PointKind::Crossing;
            }
            return Ok(outcome);
        }
        next.curve_parameter = current.curve_parameter + step;
        let adjusted = self.adjust_derivatives(current, next)?;
        let (uv_curves, rails) = self.rail_curves(current, next, step, adjusted)?;
        let mut tolerance = self.approx_tolerance;
        if adjusted {
            tolerance *= 10.0;
        }
        match self.test_points(&uv_curves, &rails, step, tolerance * tolerance)? {
            Some(max_sq) => {
                outcome.deviation = sqrt(max_sq);
                outcome.deviation_written = true;
                outcome.good = true;
            }
            None => {
                if next.kind == PointKind::Singularity {
                    next.kind = PointKind::Crossing;
                }
            }
        }
        Ok(outcome)
    }

    /// The `m_bClippedByThroughPoint` branch of `ComputeNextPoint`: when a
    /// parameter tangent's length has changed threefold over the span, take
    /// the tangents from a step of `1e-5`. Whether it did.
    ///
    /// The step's point values are filled in at the *full* step's UVs, as
    /// SMLib does: its `ComputeNextPoint` passes `rNextTSP.UVPos` to
    /// `ComputePointValues` after `StepSolve`. The resulting chord is
    /// therefore the full step's direction, divided by `1e-5` and unitized,
    /// and the surface derivatives are the full endpoint's; the short step
    /// supplies no direction of its own. This is kept for donor parity and
    /// is pinned by the curved through-point cases of the gate.
    fn adjust_derivatives(
        &self,
        current: &mut TsectPoint,
        next: &TsectPoint,
    ) -> Result<bool, StepRefusal> {
        if !self.clipped_by_through_point {
            return Ok(false);
        }
        for surface in 0..2 {
            let start_length = current.uv_deriv[surface].length();
            let end_length = next.uv_deriv[surface].length();
            if start_length > 3.0 * end_length || end_length > 3.0 * start_length {
                let small = 1.0e-5;
                let stepped =
                    self.ball
                        .step_solve(current, small, self.domains[0], self.domains[1])?;
                if !stepped.found {
                    return Ok(false);
                }
                let mut small_step = stepped.next;
                if self
                    .compute_point_values(next.uv, &mut small_step, None)
                    .is_err()
                {
                    return Ok(false);
                }
                let mut tangent = (small_step.curve_point - current.curve_point) / small;
                tangent.unitize()?;
                current.curve_tangent = tangent;
                for s in 0..2 {
                    let new_deriv = drop_vector(small_step.du(s), small_step.dv(s), tangent)?;
                    if new_deriv.length() < current.uv_deriv[s].length() {
                        current.uv_deriv[s] = new_deriv;
                    }
                }
                return Ok(true);
            }
        }
        Ok(false)
    }
}

/// The fillet intersector as [`Tracer::trace_curve`] steps it, remembering
/// why a step was refused where the trace itself only sees a refusal.
struct FilletSteps<'a, 'b> {
    intersector: &'b FilletIntersector<'a>,
    refusal: Cell<Option<StepRefusal>>,
}

impl FilletSteps<'_, '_> {
    /// The intersector as the trace's state has it.
    fn with(&self, state: TraceState) -> FilletIntersector<'_> {
        FilletIntersector {
            trace_direction: state.trace_direction,
            clipped_by_through_point: state.clipped_by_through_point,
            ..*self.intersector
        }
    }

    fn solver_refusal(&self, refusal: StepRefusal) -> SolverRefusal {
        match refusal {
            StepRefusal::Failed(fault) => SolverRefusal::Failed(fault),
            other => {
                self.refusal.set(Some(other));
                SolverRefusal::NotPorted(SolverOperation::Intersect)
            }
        }
    }
}

impl TraceSteps for FilletSteps<'_, '_> {
    fn domains(&self) -> [Extent2; 2] {
        self.intersector.domains
    }

    fn approx_tolerance(&self) -> f64 {
        self.intersector.approx_tolerance
    }

    fn angle_tolerance(&self) -> f64 {
        self.intersector.angle_tolerance
    }

    fn stepping_parallel_tolerance(&self) -> f64 {
        EFF_ZERO
    }

    /// `SmSurfaceIntersector::ComputeStepSize`, which the fillet intersector
    /// keeps: the old step, or a hundred tolerances to start.
    fn compute_step_size(&self, old_step: f64) -> f64 {
        if crate::scalar::is_zero(old_step) {
            100.0 * self.intersector.approx_tolerance
        } else {
            old_step
        }
    }

    fn compute_next_point(
        &self,
        last: &mut TsectPoint,
        step: f64,
        next: &mut TsectPoint,
        _curve_end: Option<Vector3>,
        state: TraceState,
    ) -> Result<StepOutcome, SolverRefusal> {
        self.with(state)
            .compute_next_point(last, step, next)
            .map_err(|refusal| self.solver_refusal(refusal))
    }

    fn compute_point_values(
        &self,
        uvs: [Vector2; 2],
        point: &mut TsectPoint,
        previous: Option<&TsectPoint>,
        _curve_end: Option<Vector3>,
        _step_size: Option<&mut f64>,
        state: TraceState,
    ) -> Result<(), SolverRefusal> {
        self.with(state)
            .compute_point_values(uvs, point, previous)
            .map_err(|refusal| self.solver_refusal(refusal))
    }

    /// The fillet trace has no earlier curves to be on.
    fn is_point_on_curve(
        &self,
        _point: Vector3,
        curves: &[IntersectionCurve],
        _skip: Option<usize>,
    ) -> Result<bool, SolverRefusal> {
        if curves.is_empty() {
            Ok(false)
        } else {
            Err(SolverRefusal::NotPorted(SolverOperation::Intersect))
        }
    }
}

impl FilletIntersector<'_> {
    /// `SmSurfaceIntersector::DoPointIntersection` for the fillet
    /// intersector: the through points' values, the trace through them from
    /// the first (and, when asked, back before it), up to the call of
    /// `FlushCurve`, which would build the fillet surface and its rails and
    /// is not carried. SMLib's `FlushCurve` of the fillet intersector adds no
    /// curve to the outputs, so its singularity retry never runs and is
    /// not ported either.
    ///
    /// The through points are updated as SMLib updates the caller's: kind
    /// and curve parameter (`-1` for a point not reached).
    pub fn do_point_intersection(
        &self,
        through_points: &mut [TsectPoint],
        request: TraceRequest,
    ) -> Result<FilletTrace, StepRefusal> {
        if through_points.is_empty() {
            return Err(Fault::Degenerate.into());
        }
        let mut fi = *self;
        fi.trace_direction = 1.0;
        if let Some(tolerance) = request.approx_tolerance {
            fi.approx_tolerance = tolerance;
        }
        if let Some(tolerance) = request.angle_tolerance {
            fi.angle_tolerance = tolerance;
        }
        let steps = FilletSteps {
            intersector: &fi,
            refusal: Cell::new(None),
        };
        let from_direction = |direction: Option<Vector3>| {
            direction.map(|direction| TsectPoint {
                kind: PointKind::Unknown,
                curve_tangent: direction,
                ..TsectPoint::default()
            })
        };
        let mut previous = from_direction(request.start_direction);
        let end_previous = from_direction(request.end_direction);
        let mut tracer = Tracer {
            trace_direction: 1.0,
            extension_distance: fi.extension_distance,
            ..Tracer::default()
        };
        let mut start = None;
        // Where each point of the tracer's through list came from.
        let mut origins: Vec<usize> = Vec::new();
        let last_index = through_points.len() - 1;
        for (index, point) in through_points.iter_mut().enumerate() {
            point.kind = PointKind::Unknown;
            point.curve_parameter = -1.0;
            let uvs = point.uv;
            if index > 0 {
                previous = None;
            }
            if index == last_index {
                previous = end_previous;
            }
            fi.compute_point_values(uvs, point, previous.as_ref())
                .map_err(|refusal| match refusal {
                    StepRefusal::Failed(_) => StepRefusal::Failed(Fault::Degenerate),
                    other => other,
                })?;
            if index == 0
                && !point
                    .position(0)
                    .closer_than(fi.approx_tolerance, point.position(1))
            {
                return Err(Fault::Degenerate.into());
            }
            if index == 0 {
                start = Some(*point);
                point.curve_parameter = 0.0;
                if !(request.start_direction.is_some() && last_index > 0) {
                    tracer.has_start_point = true;
                }
            } else if index == last_index
                && request.start_direction.is_some()
                && !request.extend_after_end
            {
                tracer.end_point = Some(*point);
            } else {
                let duplicate = tracer.through_points.iter().any(|other| {
                    (point.uv[0] - other.uv[0]).length_squared() < EFF_ZERO_SQ
                        && (point.uv[1] - other.uv[1]).length_squared() < EFF_ZERO_SQ
                });
                if !duplicate {
                    tracer.through_points.push(*point);
                    origins.push(index);
                }
            }
        }
        let mut start = start.ok_or(Fault::Degenerate)?;
        // `AddStartPoint` on an empty list adds the point; the trace takes it
        // back out.
        tracer.start_points.clear();
        tracer.trace_direction = 1.0;
        if fi.extension_distance > 0.0 && request.extend_after_end && tracer.end_point.is_none() {
            tracer.extension_clipping = true;
        }
        if let Some(direction) = request.start_direction
            && direction.dot(start.curve_tangent) < -EFF_ZERO
        {
            start.curve_tangent = -start.curve_tangent;
            start.uv_deriv = [-start.uv_deriv[0], -start.uv_deriv[1]];
            tracer.trace_direction = -1.0;
        }
        let input_direction = tracer.trace_direction;
        let refused = |steps: &FilletSteps, refusal: SolverRefusal| match refusal {
            SolverRefusal::Failed(fault) => StepRefusal::Failed(fault),
            SolverRefusal::NotPorted(_) => steps
                .refusal
                .take()
                .unwrap_or(StepRefusal::NotPorted(SolverOperation::Intersect)),
        };
        let traced = match tracer.trace_curve(&steps, &start, &[]) {
            Ok(()) => true,
            Err(SolverRefusal::Failed(_)) => false,
            Err(refusal) => return Err(refused(&steps, refusal)),
        };
        if traced
            && request.extend_before_start
            && !tracer.curve_is_closed
            && !tracer.curve_points.is_empty()
            && input_direction == tracer.trace_direction
        {
            tracer.reverse_curve_direction()?;
            tracer.trace_direction *= -1.0;
            tracer
                .trace_curve(&steps, &start, &[])
                .map_err(|refusal| refused(&steps, refusal))?;
            tracer.reverse_curve_direction()?;
            tracer.trace_direction *= -1.0;
        }
        // The trace wrote the parameters it reached into its own copies.
        for (origin, parameter) in origins.iter().zip(&tracer.through_parameters) {
            through_points[*origin].curve_parameter = *parameter;
        }
        if let Some(end) = tracer.end_point {
            through_points[last_index].curve_parameter = end.curve_parameter;
        }
        Ok(FilletTrace {
            traced,
            points: if traced {
                std::mem::take(&mut tracer.curve_points)
            } else {
                Vec::new()
            },
            through_points: through_points.to_vec(),
            curve_is_closed: tracer.curve_is_closed,
        })
    }
}
