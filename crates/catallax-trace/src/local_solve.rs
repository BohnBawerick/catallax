// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: statuses become `Fault`; `SmExtentNd` becomes a
// vector of intervals; the bisection step and oscillation queue are folded
// in; a periodic wrap needing more than a million period steps finishes
// with one floating remainder instead of looping, so work stays bounded.
// Adapted from rust-brep crates/kernel/src/local_solvers.rs, bisection.rs
// and iteration_queue.rs, commit 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab,
// translations of SMLib SmLocalSolve1d.cpp, SmLocalSolveNd.cpp and
// SmExtentNd.cpp.
//! SMLib's local root finders: damped one-variable Newton iteration with
//! bracketing, bisection and Brent minimisation (`SmLocalSolve1d`), and
//! multi-variable Newton-Raphson with clamping to parameter boxes
//! (`SmLocalSolveNd`).
//!
//! The solvers propose candidates for the global search. Their accuracies
//! steer the search, never what Catallax accepts.

use crate::Fault;
use crate::extent::Extent1;
use crate::matrix::Matrix;
use crate::scalar::{BIG_DOUBLE, EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, UNDEFINED, is_zero};

/// `SM_UNDEF_ULONG`.
const UNDEFINED_COUNT: u64 = 123_456_789;
/// The Brent minimiser's absolute tolerance floor (`ZEPS`).
const ZEPS: f64 = 1.0e-12;
/// Periodic wraps longer than this finish with a floating remainder.
const WRAP_STEPS: u32 = 1_000_000;

/// `SMSIGN(a, b)`: `|a|` with the sign of `b`.
fn sign(a: f64, b: f64) -> f64 {
    if b >= 0.0 { a.abs() } else { -a.abs() }
}

fn snap_ends(interval: Extent1, mut value: f64) -> f64 {
    if (value - interval.min).abs() < EFF_ZERO {
        value = interval.min;
    }
    if (value - interval.max).abs() < EFF_ZERO {
        value = interval.max;
    }
    value
}

/// `SmExtent1d::PeriodicWrap`: moves `value` into the interval by whole
/// periods, snapping to an end within `SM_EFF_ZERO`.
pub fn periodic_wrap(interval: Extent1, value: f64) -> f64 {
    let period = interval.max - interval.min;
    if period <= EFF_ZERO {
        return value;
    }
    let mut value = snap_ends(interval, value);
    if value + period == value {
        return value;
    }
    let mut steps = 0;
    while value < interval.min {
        if steps == WRAP_STEPS {
            return snap_ends(
                interval,
                interval.min + (value - interval.min).rem_euclid(period),
            );
        }
        steps += 1;
        value = snap_ends(interval, value + period);
    }
    while value > interval.max {
        if steps == WRAP_STEPS {
            return snap_ends(
                interval,
                interval.min + (value - interval.min).rem_euclid(period),
            );
        }
        steps += 1;
        value = snap_ends(interval, value - period);
    }
    value
}

/// `SmSolverAlgorithmType`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolverAlgorithm {
    /// Damped Newton iteration (`SM_SA_NEWTON`).
    Newton,
    /// Bisection (`SM_SA_BISECTION`).
    Bisection,
    /// Newton on the derivative with a Brent maximiser fallback.
    BrentMaximize,
    /// Newton on the derivative with a Brent minimiser fallback.
    BrentMinimize,
}

/// `SmTerminationReasonType`, with SMLib's numbering.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum TerminationReason {
    /// Not yet run.
    Undefined = 0,
    /// Converged to the desired accuracy.
    FoundAnswerConverged = 1,
    /// Within the acceptable accuracy only.
    FoundAnswerClose = 2,
    /// Stopped against the parameter bounds.
    OutOfBounds = 3,
    /// The Newton system was singular.
    BadJacobian = 4,
    /// No convergence.
    UnableToConverge = 5,
    /// Still iterating.
    StillIterating = 6,
}

/// `SmBoundaryHandlerType`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundaryHandler {
    /// Stop after more than the allowed number of significant bound hits.
    TotalBoundaryHits,
    /// Stop when every variable is clamped three times.
    HitsAllBoundaries,
    /// Never stop for bound hits.
    IgnoreBoundaryHits,
}

/// `SmEvalFunctionObject`: a scalar function of one variable.
///
/// Outputs are written through references, as in SMLib, so an evaluator
/// that leaves an output alone keeps the caller's value.
pub trait EvalFunction {
    /// Evaluates the function and its derivative at `t`; `found` reports an
    /// exact answer. The default is SMLib's base class, which fails.
    fn evaluate(
        &mut self,
        _t: f64,
        _value: &mut f64,
        _derivative: &mut f64,
        _found: &mut bool,
        _signal_errors: bool,
    ) -> Result<(), Fault> {
        Err(Fault::Degenerate)
    }
}

/// The outcome of [`LocalSolve1d::solve_it`]: SMLib's `bFoundSolution`
/// and `dFoundT`. The parameter is the last iterate even when no root was
/// found, as SMLib reports it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Root1d {
    /// Whether the function was within the acceptable accuracy there.
    pub found: bool,
    /// The parameter reached.
    pub t: f64,
}

/// `SmIterationValue`.
#[derive(Clone, Copy, Debug)]
struct IterationValue {
    t: f64,
    f_of_t: f64,
    f_prime_of_t: f64,
    prev_step: f64,
}

impl Default for IterationValue {
    fn default() -> Self {
        Self {
            t: UNDEFINED,
            f_of_t: UNDEFINED,
            f_prime_of_t: UNDEFINED,
            prev_step: UNDEFINED,
        }
    }
}

/// `SmIterQueue1d`: the last eight Newton steps, to detect oscillation.
struct IterationQueue {
    steps: [f64; 8],
    count: usize,
    next: usize,
    full: bool,
}

impl IterationQueue {
    fn new() -> Self {
        Self {
            steps: [0.0; 8],
            count: 0,
            next: 0,
            full: false,
        }
    }

    fn add(&mut self, previous_step: f64) {
        self.steps[self.next] = previous_step;
        self.next += 1;
        if self.next >= self.steps.len() {
            self.next = 0;
            self.full = true;
        }
        if self.count < self.steps.len() {
            self.count += 1;
        }
    }

    fn is_oscillating(&self, step: u64) -> bool {
        if step != 2 || self.count < 6 {
            return false;
        }
        let top = if self.full {
            self.steps.len()
        } else {
            self.next
        };
        let (mut sum1, mut sum2) = (0.0, 0.0);
        for (index, value) in self.steps[..top].iter().enumerate() {
            if index % 2 == 0 {
                sum1 += value;
            } else {
                sum2 += value;
            }
        }
        (sum1 + sum2).abs() < (sum1.abs() + sum2.abs()) / 100.0
    }
}

/// A sign-change bracket: the function is negative at `lower_t` and
/// positive at `upper_t`, in either parameter order.
#[derive(Clone, Copy, Debug)]
struct Bracket {
    lower_t: f64,
    upper_t: f64,
    lower_value: f64,
    upper_value: f64,
}

/// `SmLocalSolve1d`.
pub struct LocalSolve1d<E: EvalFunction> {
    use_interval: bool,
    interval: Extent1,
    algorithm: SolverAlgorithm,
    /// The function being solved.
    pub evaluator: E,
    desired_accuracy: f64,
    acceptable_accuracy: f64,
    max_iter: u64,
    max_boundary_hits: u64,
    is_periodic: bool,
    have_lower_bracket: bool,
    lower_bracket_f: f64,
    lower_bracket_t: f64,
    have_upper_bracket: bool,
    upper_bracket_f: f64,
    upper_bracket_t: f64,
    prev2: IterationValue,
    prev1: IterationValue,
    curr: IterationValue,
    iter: u64,
    found_accuracy: f64,
    lower_bound_hits: u64,
    upper_bound_hits: u64,
    termination_reason: TerminationReason,
}

impl<E: EvalFunction> LocalSolve1d<E> {
    /// A solver over the whole real line.
    pub fn new(evaluator: E) -> Self {
        Self {
            use_interval: false,
            interval: Extent1::default(),
            algorithm: SolverAlgorithm::Newton,
            evaluator,
            desired_accuracy: EFF_ZERO,
            acceptable_accuracy: UNDEFINED,
            max_iter: 100,
            max_boundary_hits: UNDEFINED_COUNT,
            is_periodic: false,
            have_lower_bracket: false,
            lower_bracket_f: UNDEFINED,
            lower_bracket_t: UNDEFINED,
            have_upper_bracket: false,
            upper_bracket_f: UNDEFINED,
            upper_bracket_t: UNDEFINED,
            prev2: IterationValue::default(),
            prev1: IterationValue::default(),
            curr: IterationValue::default(),
            iter: 0,
            found_accuracy: UNDEFINED,
            lower_bound_hits: UNDEFINED_COUNT,
            upper_bound_hits: UNDEFINED_COUNT,
            termination_reason: TerminationReason::Undefined,
        }
    }

    /// A solver limited to `interval`, wrapped when `is_periodic`.
    pub fn with_interval(evaluator: E, interval: Extent1, is_periodic: bool) -> Self {
        let mut solver = Self::new(evaluator);
        solver.use_interval = true;
        solver.interval = interval;
        solver.is_periodic = is_periodic;
        solver
    }

    /// `GetFoundAccuracy`: the function value at the answer.
    pub fn found_accuracy(&self) -> f64 {
        self.found_accuracy
    }

    /// `GetTerminationReason`.
    pub fn termination_reason(&self) -> TerminationReason {
        self.termination_reason
    }

    /// `SetSolutionInterval`: replaces the interval without enabling it.
    pub fn set_solution_interval(&mut self, interval: Extent1) {
        self.interval = interval;
    }

    /// `SetSolverAlgorithm`.
    pub fn set_solver_algorithm(&mut self, algorithm: SolverAlgorithm) {
        self.algorithm = algorithm;
    }

    /// `SetDesiredAccuracy`.
    pub fn set_desired_accuracy(&mut self, accuracy: f64) {
        self.desired_accuracy = accuracy;
    }

    /// `SolveByBisection`: false position clamped to `[0.15, 0.85]`.
    fn bisect_bracket(
        &mut self,
        mut levels: u64,
        mut bracket: Bracket,
    ) -> Result<(f64, f64), Fault> {
        loop {
            let Bracket {
                lower_t,
                upper_t,
                lower_value,
                upper_value,
            } = bracket;
            if levels == 0 {
                return Ok(if lower_value + upper_value >= 0.0 {
                    (lower_t, lower_value)
                } else {
                    (upper_t, upper_value)
                });
            }
            let fraction = (-lower_value / (upper_value - lower_value)).clamp(0.15, 0.85);
            let mid_t = lower_t + fraction * (upper_t - lower_t);
            let scale = mid_t.abs();
            if (mid_t - upper_t).abs() < EFF_ZERO_SQ * scale
                || (mid_t - lower_t).abs() < EFF_ZERO_SQ * scale
            {
                return Ok(if lower_value.abs() < upper_value.abs() {
                    (lower_t, lower_value)
                } else if upper_value.abs() < lower_value.abs() {
                    (upper_t, upper_value)
                } else {
                    ((lower_t + upper_t) / 2.0, (lower_value + upper_value) / 2.0)
                });
            }
            let (mut value, mut derivative, mut found) = (0.0, 0.0, false);
            self.evaluator
                .evaluate(mid_t, &mut value, &mut derivative, &mut found, true)
                .map_err(|_| Fault::Degenerate)?;
            if found || value == 0.0 {
                return Ok((mid_t, value));
            }
            if value > 0.0 {
                bracket.upper_t = mid_t;
                bracket.upper_value = value;
            } else {
                bracket.lower_t = mid_t;
                bracket.lower_value = value;
            }
            levels -= 1;
        }
    }

    fn current_bracket(&self) -> Bracket {
        Bracket {
            lower_t: self.lower_bracket_t,
            upper_t: self.upper_bracket_t,
            lower_value: self.lower_bracket_f,
            upper_value: self.upper_bracket_f,
        }
    }

    fn evaluate_signed(&mut self, t: f64, maximize: bool) -> Result<(f64, f64), Fault> {
        let (mut value, mut derivative, mut found) = (0.0, 0.0, false);
        self.evaluator
            .evaluate(t, &mut value, &mut derivative, &mut found, true)
            .map_err(|_| Fault::Degenerate)?;
        Ok(if maximize {
            (-value, -derivative)
        } else {
            (value, derivative)
        })
    }

    /// `SolveByBrent`: Brent minimisation (or maximisation) with derivatives
    /// inside `[ax, cx]` from `bx`. Returns the parameter and value found.
    pub fn solve_by_brent(
        &mut self,
        ax: f64,
        bx: f64,
        cx: f64,
        tol: f64,
        maximize: bool,
    ) -> Result<(f64, f64), Fault> {
        let mut bracket_min = if ax < cx { ax } else { cx };
        let mut bracket_max = if ax > cx { ax } else { cx };
        let mut u;
        let (mut fu, mut du);
        let mut x = bx;
        let (mut fx, mut dx) = self.evaluate_signed(x, maximize)?;
        let (mut v, mut fv, mut dv) = (bx, 0.0, 0.0);
        let (mut w, mut fw, mut dw) = (bx, 0.0, 0.0);
        let mut delta_x: f64 = 0.0;
        let mut e: f64 = 0.0;
        let mut iter: u64 = 1;
        while iter <= self.max_iter {
            let tol1 = tol * x.abs() + ZEPS;
            let tol2 = 2.0 * tol1;
            let mut fall_back = true;
            if e.abs() > tol1 {
                let mut d1 = 2.0 * (bracket_max - bracket_min);
                let mut d2 = d1;
                if dw != dx && w != x {
                    d1 = (w - x) * dx / (dx - dw);
                }
                if dv != dx && v != x {
                    d2 = (v - x) * dx / (dx - dv);
                }
                let u1 = x + d1;
                let u2 = x + d2;
                let ok1 = (bracket_min - u1) * (u1 - bracket_max) > 0.0 && dx * d1 <= 0.0;
                let ok2 = (bracket_min - u2) * (u2 - bracket_max) > 0.0 && dx * d2 <= 0.0;
                let olde = e;
                e = delta_x;
                if ok1 || ok2 {
                    delta_x = if ok1 && ok2 {
                        if d1.abs() < d2.abs() { d1 } else { d2 }
                    } else if ok1 {
                        d1
                    } else {
                        d2
                    };
                    if delta_x.abs() <= (0.5 * olde).abs() {
                        fall_back = false;
                        u = x + delta_x;
                        if u - bracket_min < tol2 || bracket_max - u < tol2 {
                            delta_x = sign(tol1, delta_x);
                        }
                    }
                }
            }
            if fall_back {
                e = if dx >= 0.0 {
                    bracket_min - x
                } else {
                    bracket_max - x
                };
                delta_x = 0.5 * e;
            }
            if delta_x.abs() >= tol1 {
                u = x + delta_x;
                if u == x {
                    return Ok((x, fx));
                }
                (fu, du) = self.evaluate_signed(u, maximize)?;
            } else {
                u = x + sign(tol1, delta_x);
                if u <= bracket_min {
                    u = (x + bracket_min) / 2.0;
                }
                if u >= bracket_max {
                    u = (x + bracket_max) / 2.0;
                }
                (fu, du) = self.evaluate_signed(u, maximize)?;
                if u == x || fu >= fx {
                    return Ok((x, fx));
                }
            }
            if fu <= fx {
                if u >= x {
                    bracket_min = x;
                } else {
                    bracket_max = x;
                }
                (v, fv, dv) = (w, fw, dw);
                (w, fw, dw) = (x, fx, dx);
                (x, fx, dx) = (u, fu, du);
            } else {
                if u < x {
                    bracket_min = u;
                } else {
                    bracket_max = u;
                }
                if fu <= fw || w == x {
                    (v, fv, dv) = (w, fw, dw);
                    (w, fw, dw) = (u, fu, du);
                } else if fu < fv || v == x || v == w {
                    (v, fv, dv) = (u, fu, du);
                }
            }
            iter += 1;
        }
        Err(Fault::Degenerate)
    }

    /// `SolveOscillating`: bracket and bisect, run Brent for the minimising
    /// and maximising algorithms, or settle on the best interior value.
    fn solve_oscillating(&mut self, guess_t: f64) -> Result<Root1d, Fault> {
        const GUESSES: [f64; 8] = [0.25, 0.75, 0.13, 0.38, 0.82, 0.87, 0.01, 0.99];
        let mut found = false;
        let mut temp = IterationValue {
            t: 0.5 * (self.prev1.t + self.curr.t),
            ..IterationValue::default()
        };
        self.evaluator.evaluate(
            temp.t,
            &mut temp.f_of_t,
            &mut temp.f_prime_of_t,
            &mut found,
            true,
        )?;
        let mut best_interior_t = temp.t;
        let mut best_interior_f = temp.f_of_t;

        if !self.have_lower_bracket || !self.have_upper_bracket {
            let (curr_t, prev_t, curr_f) = (self.curr.t, self.prev1.t, self.curr.f_of_t);
            for guess in GUESSES {
                let mid_t = prev_t + guess * (curr_t - prev_t);
                temp.t = mid_t;
                if self
                    .evaluator
                    .evaluate(
                        temp.t,
                        &mut temp.f_of_t,
                        &mut temp.f_prime_of_t,
                        &mut found,
                        true,
                    )
                    .is_err()
                {
                    continue;
                }
                if temp.f_of_t.abs() < best_interior_f.abs() {
                    best_interior_f = temp.f_of_t;
                    best_interior_t = mid_t;
                }
                if temp.f_of_t * curr_f <= 0.0 {
                    break;
                }
            }
            if temp.f_of_t * curr_f <= 0.0 {
                let previous = (prev_t, self.prev1.f_of_t);
                let current = (curr_t, curr_f);
                // The bracket's other end: the iterate on the guess's side.
                let other = if guess_t < temp.t {
                    if curr_t < prev_t { current } else { previous }
                } else if curr_t > prev_t {
                    current
                } else {
                    previous
                };
                let ((mut lo_t, mut lo_f), (mut hi_t, mut hi_f)) = if temp.f_of_t > 0.0 {
                    (other, (temp.t, temp.f_of_t))
                } else {
                    ((temp.t, temp.f_of_t), other)
                };
                if self.have_lower_bracket
                    && (self.lower_bracket_t - hi_t).abs() < (lo_t - hi_t).abs()
                {
                    lo_t = self.lower_bracket_t;
                    lo_f = self.lower_bracket_f;
                }
                if self.have_upper_bracket
                    && (self.upper_bracket_t - lo_t).abs() < (hi_t - lo_t).abs()
                {
                    hi_t = self.upper_bracket_t;
                    hi_f = self.upper_bracket_f;
                }
                self.have_lower_bracket = true;
                self.lower_bracket_t = lo_t;
                self.lower_bracket_f = lo_f;
                self.have_upper_bracket = true;
                self.upper_bracket_t = hi_t;
                self.upper_bracket_f = hi_f;
            }
        }

        if self.have_lower_bracket
            && self.have_upper_bracket
            && let Ok((t, f)) = self.bisect_bracket(40, self.current_bracket())
        {
            self.curr.t = t;
            self.curr.f_of_t = f;
            self.found_accuracy = f;
            self.termination_reason = if f.abs() < self.desired_accuracy {
                TerminationReason::FoundAnswerConverged
            } else {
                TerminationReason::FoundAnswerClose
            };
            return Ok(Root1d { found: true, t });
        }

        let brent = matches!(
            self.algorithm,
            SolverAlgorithm::BrentMinimize | SolverAlgorithm::BrentMaximize
        );
        if brent {
            let maximize = self.algorithm == SolverAlgorithm::BrentMaximize;
            let (t, _) = self.solve_by_brent(
                self.prev1.t,
                best_interior_t,
                self.curr.t,
                self.desired_accuracy,
                maximize,
            )?;
            self.curr.t = t;
        } else if best_interior_f.abs() < self.curr.f_of_t.abs() {
            self.curr.t = best_interior_t;
        } else {
            self.curr.t = 0.5 * (self.prev1.t + self.curr.t);
        }
        self.evaluator.evaluate(
            self.curr.t,
            &mut self.curr.f_of_t,
            &mut self.curr.f_prime_of_t,
            &mut found,
            true,
        )?;
        Ok(Root1d {
            found: brent || found || self.curr.f_of_t.abs() <= self.acceptable_accuracy,
            t: self.curr.t,
        })
    }

    fn classify(&self, magnitude: f64) -> TerminationReason {
        if magnitude < self.desired_accuracy {
            TerminationReason::FoundAnswerConverged
        } else if magnitude < self.acceptable_accuracy {
            TerminationReason::FoundAnswerClose
        } else {
            TerminationReason::UnableToConverge
        }
    }

    fn converged(&mut self) -> Result<Root1d, Fault> {
        self.found_accuracy = self.curr.f_of_t;
        self.termination_reason = TerminationReason::FoundAnswerConverged;
        Ok(Root1d {
            found: true,
            t: self.curr.t,
        })
    }

    /// `SolveIt`: damped Newton iteration with bracketing, oscillation
    /// handling, Brent hand-off and a bisection fallback. Fails only when the
    /// function cannot be evaluated.
    pub fn solve_it(&mut self, guess_t: f64, acceptable_accuracy: f64) -> Result<Root1d, Fault> {
        let mut guess_t = guess_t;
        self.have_lower_bracket = false;
        self.have_upper_bracket = false;
        self.lower_bracket_f = -BIG_DOUBLE;
        self.upper_bracket_f = BIG_DOUBLE;
        self.iter = 0;
        self.lower_bound_hits = 0;
        self.upper_bound_hits = 0;
        self.max_boundary_hits = 2;
        self.termination_reason = TerminationReason::UnableToConverge;
        self.acceptable_accuracy = acceptable_accuracy;
        self.curr.prev_step = 0.0;
        if self.use_interval && !self.interval.contains_value(guess_t, 0.0) {
            guess_t = self.interval.clamp_value(guess_t);
        }
        self.curr.t = guess_t;
        let mut found = false;
        let mut status = self.evaluator.evaluate(
            self.curr.t,
            &mut self.curr.f_of_t,
            &mut self.curr.f_prime_of_t,
            &mut found,
            false,
        );
        if status.is_err() && self.use_interval {
            let size = self.interval.length();
            self.curr.t = guess_t + size / 1000.0;
            if !self.interval.contains_value(self.curr.t, 0.0) {
                self.curr.t = guess_t - size / 1000.0;
            }
            status = self.evaluator.evaluate(
                self.curr.t,
                &mut self.curr.f_of_t,
                &mut self.curr.f_prime_of_t,
                &mut found,
                true,
            );
        }
        status.map_err(|_| Fault::Degenerate)?;
        if found {
            return self.converged();
        }

        let mut cut_factor = 1.0;
        let mut minimum = self.curr.f_of_t.abs();
        let mut min_count: u64 = 0;
        let mut queue = IterationQueue::new();
        self.iter = 0;
        while self.iter < self.max_iter {
            if self.curr.f_of_t < 0.0 {
                self.have_lower_bracket = true;
                if self.curr.f_of_t > self.lower_bracket_f {
                    self.lower_bracket_t = self.curr.t;
                    self.lower_bracket_f = self.curr.f_of_t;
                }
            }
            if self.curr.f_of_t > 0.0 {
                self.have_upper_bracket = true;
                if self.curr.f_of_t < self.upper_bracket_f {
                    self.upper_bracket_t = self.curr.t;
                    self.upper_bracket_f = self.curr.f_of_t;
                }
            }
            if self.curr.f_of_t.abs() < self.desired_accuracy {
                return self.converged();
            }
            if queue.is_oscillating(2) {
                return self.solve_oscillating(guess_t);
            }
            if is_zero(self.curr.f_prime_of_t) {
                self.curr.f_prime_of_t = if self.curr.f_prime_of_t < 0.0 {
                    -1.0
                } else {
                    1.0
                };
            }
            let mut delta_t = -self.curr.f_of_t / self.curr.f_prime_of_t;
            let mut step_limit = self.interval.length() * 0.08 + EFF_ZERO;
            if delta_t.abs() > step_limit {
                step_limit *= 0.99;
                if delta_t < 0.0 {
                    step_limit = -step_limit;
                }
                delta_t = step_limit;
            }
            if self.iter > 20 && self.iter >= self.max_iter.wrapping_sub(1) {
                cut_factor *= 0.5;
                if cut_factor > 0.0001 {
                    self.iter -= 20;
                }
            }
            delta_t *= cut_factor;
            let mut new_t = self.curr.t + delta_t;
            if self.use_interval && !self.interval.contains_value(new_t, 0.0) {
                if self.is_periodic {
                    new_t = periodic_wrap(self.interval, new_t);
                } else {
                    if new_t < self.interval.min {
                        new_t = self.interval.min;
                        self.lower_bound_hits += 1;
                        self.upper_bound_hits = 0;
                    } else if new_t > self.interval.max {
                        new_t = self.interval.max;
                        self.upper_bound_hits += 1;
                        self.lower_bound_hits = 0;
                    }
                    if self.lower_bound_hits > self.max_boundary_hits
                        || self.upper_bound_hits > self.max_boundary_hits
                    {
                        self.termination_reason = TerminationReason::OutOfBounds;
                        break;
                    }
                }
            } else {
                self.lower_bound_hits = 0;
                self.upper_bound_hits = 0;
            }
            self.curr.prev_step = if self.iter > 0 {
                self.curr.t - self.prev1.t
            } else {
                0.0
            };
            self.prev2 = if self.iter > 1 { self.prev1 } else { self.curr };
            self.prev1 = self.curr;
            queue.add(self.prev1.prev_step);
            self.curr.t = new_t;
            self.evaluator
                .evaluate(
                    self.curr.t,
                    &mut self.curr.f_of_t,
                    &mut self.curr.f_prime_of_t,
                    &mut found,
                    true,
                )
                .map_err(|_| Fault::Degenerate)?;
            if found || self.curr.f_of_t.abs() < self.desired_accuracy {
                return self.converged();
            }
            if self.curr.f_of_t.abs() >= minimum {
                min_count += 1;
                if min_count >= 4 && self.curr.f_of_t.abs() < self.acceptable_accuracy {
                    break;
                }
            } else {
                min_count = 0;
                minimum = self.curr.f_of_t.abs();
            }
            for maximize in [true, false] {
                let wanted = if maximize {
                    SolverAlgorithm::BrentMaximize
                } else {
                    SolverAlgorithm::BrentMinimize
                };
                let beyond = if maximize {
                    self.curr.f_of_t > self.prev1.f_of_t && self.curr.f_of_t > self.prev2.f_of_t
                } else {
                    self.curr.f_of_t < self.prev1.f_of_t && self.curr.f_of_t < self.prev2.f_of_t
                };
                if self.algorithm == wanted && self.iter > 1 && beyond {
                    let mut span = Extent1::new(self.prev1.t, self.prev1.t);
                    span.add_value(self.prev2.t);
                    if span.contains_value(self.curr.t, 0.0) {
                        let (t, _) = self.solve_by_brent(
                            self.prev1.t,
                            self.curr.t,
                            self.prev2.t,
                            self.desired_accuracy,
                            maximize,
                        )?;
                        return Ok(Root1d { found: true, t });
                    }
                }
            }
            if self.curr.f_of_t.abs() >= self.prev1.f_of_t.abs()
                && self.curr.f_of_t * self.prev1.f_of_t > 0.0
            {
                if self.curr.f_of_t != self.prev1.f_of_t
                    && self.curr.t != self.prev1.t
                    && self.have_lower_bracket
                    && self.have_upper_bracket
                {
                    break;
                }
                if delta_t.abs() < EFF_ZERO {
                    let magnitude = self.curr.f_of_t.abs();
                    self.found_accuracy = self.curr.f_of_t;
                    self.termination_reason = self.classify(magnitude);
                    return Ok(Root1d {
                        found: magnitude < self.acceptable_accuracy,
                        t: self.curr.t,
                    });
                }
            }
            self.iter += 1;
        }

        if self.have_lower_bracket
            && self.have_upper_bracket
            && let Ok((t, f)) = self.bisect_bracket(20, self.current_bracket())
        {
            self.curr.t = t;
            self.curr.f_of_t = f;
        }
        self.found_accuracy = self.curr.f_of_t;
        let found = self.curr.f_of_t.abs() < self.acceptable_accuracy;
        if found {
            self.termination_reason = self.classify(self.curr.f_of_t.abs());
        }
        Ok(Root1d {
            found,
            t: self.curr.t,
        })
    }
}

/// `SmExtentNd`: one parameter interval per variable.
#[derive(Clone, Debug, PartialEq)]
pub struct ExtentN {
    /// The intervals, in variable order.
    pub intervals: Vec<Extent1>,
}

impl ExtentN {
    /// `ClampVector`: clamps (or, for periodic variables, wraps) `vector`
    /// into the bounded intervals; unbounded variables leave `clamped`
    /// unchanged. With `out_of_bounds`, only values outside their interval
    /// are clamped, and a variable moved by more than `SM_EFF_ZERO` adds one
    /// to its count. Returns whether any variable was moved that far.
    pub fn clamp_vector(
        &self,
        vector: &[f64],
        periodic: Option<&[bool]>,
        clamped: &mut Vec<f64>,
        mut out_of_bounds: Option<&mut Vec<u64>>,
    ) -> bool {
        let dimension = self.intervals.len();
        let mut moved = false;
        clamped.resize(dimension, 0.0);
        if let Some(counts) = out_of_bounds.as_deref_mut() {
            counts.resize(dimension, 0);
        }
        for (i, extent) in self.intervals.iter().enumerate() {
            let value = vector[i];
            if !extent.is_bounded() {
                continue;
            }
            if periodic.is_some_and(|flags| flags[i]) {
                clamped[i] = periodic_wrap(*extent, value);
            } else if let Some(counts) = out_of_bounds.as_deref_mut() {
                if extent.contains_value(value, 0.0) {
                    clamped[i] = value;
                } else {
                    let limited = extent.clamp_value(value);
                    clamped[i] = limited;
                    if (limited - value).abs() > EFF_ZERO {
                        counts[i] += 1;
                        moved = true;
                    }
                }
            } else {
                let limited = extent.clamp_value(value);
                if (limited - value).abs() > EFF_ZERO {
                    moved = true;
                }
                clamped[i] = limited;
            }
        }
        moved
    }
}

/// `SmEvalNFunctionsObject`: `n` functions of `n` variables.
pub trait EvalNFunctions {
    /// `m_dSolutionSize` after the last evaluation.
    fn solution_size(&self) -> f64 {
        0.0
    }

    /// `m_dScaledZero`.
    fn scaled_zero(&self) -> f64 {
        1.0
    }

    /// `IsConverged`: the solution size is below the scaled zero.
    fn is_converged(&self) -> bool {
        self.solution_size() < self.scaled_zero()
    }

    /// `HasIntervals`: whether the solver should pass its intervals in.
    fn has_intervals(&mut self) -> bool {
        false
    }

    /// `SetIntervals`; SMLib's base class fails.
    fn set_intervals(
        &mut self,
        _intervals: Option<&ExtentN>,
        _periodic: Option<&[bool]>,
    ) -> Result<(), Fault> {
        Err(Fault::Degenerate)
    }

    /// `Evaluate`: the function values at `x`, the Jacobian when asked, and
    /// whether `x` is an exact answer. SMLib's base class fails.
    fn evaluate(
        &mut self,
        _x: &[f64],
        _f: &mut [f64],
        _jacobian: Option<&mut Matrix>,
        _found: &mut bool,
    ) -> Result<(), Fault> {
        Err(Fault::Degenerate)
    }
}

/// `SmEvalNFunctionsObject::ComputeJacobian`: central differences, clamped
/// to optional intervals. A failure leaves the current variable perturbed,
/// as in SMLib.
pub fn compute_jacobian<E: EvalNFunctions + ?Sized>(
    evaluator: &mut E,
    x: &mut [f64],
    f: &[f64],
    jacobian: &mut Matrix,
    intervals: Option<&[Extent1]>,
) -> Result<(), Fault> {
    let count = f.len();
    let mut low = vec![0.0; count];
    let mut high = vec![0.0; count];
    let mut found = false;
    for jj in 0..count {
        let base = x[jj];
        let delta = EFF_ZERO_SQRT * (1.0 + base.abs());
        let mut low_t = base - delta;
        let mut high_t = base + delta;
        if let Some(intervals) = intervals {
            low_t = intervals[jj].clamp_value(low_t);
            high_t = intervals[jj].clamp_value(high_t);
        }
        x[jj] = low_t;
        evaluator.evaluate(x, &mut low, None, &mut found)?;
        x[jj] = high_t;
        evaluator.evaluate(x, &mut high, None, &mut found)?;
        let delta = high_t - low_t;
        if delta.abs() < EFF_ZERO {
            return Err(Fault::Degenerate);
        }
        for ii in 0..count {
            jacobian[ii][jj] = (high[ii] - low[ii]) / delta;
        }
        x[jj] = base;
    }
    Ok(())
}

/// The outcome of [`LocalSolveNd::solve_it`]: SMLib's `bFoundSolution` and
/// solution vector, which holds the best point reached even without a root.
#[derive(Clone, Debug, PartialEq)]
pub struct RootNd {
    /// Whether the residual sum was within the acceptable accuracy.
    pub found: bool,
    /// The variables reached.
    pub x: Vec<f64>,
}

/// `SmLocalSolveNd`.
pub struct LocalSolveNd<E: EvalNFunctions> {
    /// The equations being solved.
    pub evaluator: E,
    desired_accuracy: f64,
    acceptable_accuracy: f64,
    maximum_iterations: u64,
    intervals: Option<ExtentN>,
    periodic: Option<Vec<bool>>,
    boundary_handler: BoundaryHandler,
    max_boundary_hits: u64,
    step_cutting_enabled: bool,
    current_iteration: u64,
    termination_reason: TerminationReason,
}

impl<E: EvalNFunctions> LocalSolveNd<E> {
    /// A solver with SMLib's defaults: desired accuracy `SM_EFF_ZERO/10`,
    /// 100 iterations, at most three hits on all bounds.
    pub fn new(evaluator: E, intervals: Option<ExtentN>, periodic: Option<Vec<bool>>) -> Self {
        Self {
            evaluator,
            desired_accuracy: EFF_ZERO / 10.0,
            acceptable_accuracy: EFF_ZERO,
            maximum_iterations: 100,
            intervals,
            periodic,
            boundary_handler: BoundaryHandler::HitsAllBoundaries,
            max_boundary_hits: 3,
            step_cutting_enabled: false,
            current_iteration: 0,
            termination_reason: TerminationReason::Undefined,
        }
    }

    /// `GetTerminationReason`.
    pub fn termination_reason(&self) -> TerminationReason {
        self.termination_reason
    }

    /// `SetStepCuttingEnabled`.
    pub fn set_step_cutting_enabled(&mut self, enabled: bool) {
        self.step_cutting_enabled = enabled;
    }

    /// `SetDesiredAccuracy`.
    pub fn set_desired_accuracy(&mut self, accuracy: f64) {
        self.desired_accuracy = accuracy;
    }

    /// `SetMaximumIterations`.
    pub fn set_maximum_iterations(&mut self, iterations: u64) {
        self.maximum_iterations = iterations;
    }

    /// `SetBoundaryHandler`.
    pub fn set_boundary_handler(&mut self, handler: BoundaryHandler, max_hits: u64) {
        self.boundary_handler = handler;
        self.max_boundary_hits = max_hits;
    }

    fn done(
        &mut self,
        found: bool,
        reason: TerminationReason,
        x: Vec<f64>,
    ) -> Result<RootNd, Fault> {
        self.termination_reason = reason;
        Ok(RootNd { found, x })
    }

    /// `SolveIt`: Newton-Raphson with clamping to the intervals, the better
    /// of the clamped and the scaled clamped step, optional step cutting,
    /// boundary-hit limits and a best-point fallback. Fails only when the
    /// equations reject the intervals.
    pub fn solve_it(&mut self, guess: &[f64], acceptable_accuracy: f64) -> Result<RootNd, Fault> {
        use TerminationReason::{
            BadJacobian, FoundAnswerClose, FoundAnswerConverged, OutOfBounds, UnableToConverge,
        };
        self.termination_reason = UnableToConverge;
        let count = guess.len();
        let mut found = false;
        let mut jacobian = Matrix::new(count, count);
        let mut equations = vec![BIG_DOUBLE; count];
        let mut next_guess = vec![0.0; count];
        let mut step_unclamped = vec![0.0; count];
        if self.evaluator.has_intervals() {
            self.evaluator
                .set_intervals(self.intervals.as_ref(), self.periodic.as_deref())?;
        }
        let mut found2 = false;
        let mut jacobian2 = jacobian.clone();
        let mut next_clamped = vec![0.0; count];
        let mut next_clamped_scaled = vec![0.0; count];
        let mut equations2 = vec![0.0; count];
        self.acceptable_accuracy = acceptable_accuracy;
        let mut solution = guess.to_vec();

        if self
            .evaluator
            .evaluate(&solution, &mut equations, Some(&mut jacobian), &mut found)
            .is_err()
        {
            return self.done(true, UnableToConverge, solution);
        }
        if found {
            return self.done(true, FoundAnswerConverged, solution);
        }
        let mut sum: f64 = equations.iter().map(|v| v.abs()).sum();
        if sum < self.desired_accuracy && self.evaluator.is_converged() {
            return self.done(true, FoundAnswerConverged, solution);
        }
        let mut best_sum = sum;
        let mut best_guess = guess.to_vec();
        let mut boundary_hits: u64 = 0;
        let mut hit_all_bounds_count: u64 = 0;
        let mut clamped = false;
        let mut out_of_bounds: Vec<u64> = Vec::new();
        let mut next_sum = 0.0;

        self.current_iteration = 0;
        while self.current_iteration < self.maximum_iterations {
            for value in equations.iter_mut() {
                *value = -*value;
            }
            if jacobian
                .solve_linear_system(&equations, &mut step_unclamped)
                .is_err_and(|fault| fault == Fault::Degenerate)
            {
                self.termination_reason = BadJacobian;
                break;
            }
            let mut have_good_step = false;
            let mut step_scale = 1.0;
            'line_search: while !have_good_step {
                let mut step_size = 0.0;
                have_good_step = true;
                for ii in 0..count {
                    let this_step = step_scale * step_unclamped[ii];
                    step_size += this_step.abs();
                    next_guess[ii] = solution[ii] + this_step;
                }
                if step_size < self.desired_accuracy * self.desired_accuracy {
                    if sum < self.acceptable_accuracy {
                        return self.done(true, FoundAnswerClose, solution);
                    }
                    let reason = if clamped {
                        OutOfBounds
                    } else {
                        UnableToConverge
                    };
                    return self.done(false, reason, solution);
                }
                next_clamped.clone_from(&next_guess);
                if let Some(intervals) = &self.intervals {
                    clamped = intervals.clamp_vector(
                        &next_guess,
                        self.periodic.as_deref(),
                        &mut next_clamped,
                        Some(&mut out_of_bounds),
                    );
                }
                let mut first_failed = false;
                if self
                    .evaluator
                    .evaluate(
                        &next_clamped,
                        &mut equations,
                        Some(&mut jacobian),
                        &mut found,
                    )
                    .is_err()
                {
                    self.termination_reason = UnableToConverge;
                    first_failed = true;
                    if !clamped {
                        break 'line_search;
                    }
                }
                next_sum = equations.iter().map(|v| v.abs()).sum();
                if clamped {
                    let mut scale_factor = 1.0;
                    for ii in 0..count {
                        let base_step = next_guess[ii] - solution[ii];
                        if base_step.abs() > EFF_ZERO {
                            let this_scale = (next_clamped[ii] - solution[ii]) / base_step;
                            if this_scale < scale_factor {
                                scale_factor = this_scale;
                            }
                        }
                    }
                    if scale_factor < 1.0 - EFF_ZERO && scale_factor > 0.0 {
                        for ii in 0..count {
                            let this_step = step_scale * step_unclamped[ii];
                            next_clamped_scaled[ii] = solution[ii] + scale_factor * this_step;
                        }
                        if let Some(intervals) = &self.intervals {
                            let input = next_clamped_scaled.clone();
                            intervals.clamp_vector(
                                &input,
                                self.periodic.as_deref(),
                                &mut next_clamped_scaled,
                                Some(&mut out_of_bounds),
                            );
                        }
                        if self
                            .evaluator
                            .evaluate(
                                &next_clamped_scaled,
                                &mut equations2,
                                Some(&mut jacobian2),
                                &mut found2,
                            )
                            .is_err()
                            && first_failed
                        {
                            self.termination_reason = UnableToConverge;
                            break 'line_search;
                        }
                        let next_sum2: f64 = equations2.iter().map(|v| v.abs()).sum();
                        if next_sum2 < next_sum || first_failed {
                            next_clamped.clone_from(&next_clamped_scaled);
                            equations.clone_from(&equations2);
                            next_sum = next_sum2;
                        }
                    }
                }
                next_guess.clone_from(&next_clamped);
                found |= next_sum < self.desired_accuracy;
                if found {
                    return self.done(true, FoundAnswerConverged, next_guess);
                }
                if next_sum >= sum
                    && self.step_cutting_enabled
                    && sum > self.acceptable_accuracy
                    && step_scale.abs() > 0.0001
                {
                    step_scale /= 10.0;
                    if step_scale > 0.0 && step_scale < 0.001 {
                        step_scale = -1.0;
                    }
                    have_good_step = false;
                }
            }
            if self.current_iteration > 2
                && next_sum >= sum
                && sum < self.acceptable_accuracy / 100.0
            {
                self.termination_reason = FoundAnswerClose;
                break;
            }
            let change_factor = 0.000001;
            if self.current_iteration > 2
                && (next_sum - sum).abs() < sum * change_factor
                && sum > self.acceptable_accuracy / 100.0
            {
                self.termination_reason = UnableToConverge;
                break;
            }
            if next_sum < best_sum {
                best_sum = next_sum;
                best_guess.clone_from(&next_guess);
            }
            if next_sum > sum && sum < self.acceptable_accuracy && sum < best_sum {
                best_sum = sum;
                best_guess.clone_from(&solution);
            }
            sum = next_sum;
            let mut step_size = 0.0;
            for ii in 0..count {
                step_size += (solution[ii] - next_guess[ii]).abs();
                solution[ii] = next_guess[ii];
            }
            if step_size < EFF_ZERO_SQ {
                self.termination_reason = UnableToConverge;
                break;
            }
            if matches!(
                self.boundary_handler,
                BoundaryHandler::TotalBoundaryHits | BoundaryHandler::HitsAllBoundaries
            ) {
                let mut hits_all_boundaries = true;
                for ii in 0..out_of_bounds.len() {
                    if out_of_bounds[ii] != 0 {
                        out_of_bounds[ii] = 0;
                        let extent = self
                            .intervals
                            .as_ref()
                            .map_or_else(Extent1::default, |intervals| intervals.intervals[ii]);
                        if step_unclamped[ii].abs() > 1.0e-3 * extent.length() {
                            boundary_hits += 1;
                        }
                    } else {
                        hits_all_boundaries = false;
                        if self.boundary_handler == BoundaryHandler::HitsAllBoundaries {
                            break;
                        }
                    }
                }
                if hits_all_boundaries && !out_of_bounds.is_empty() {
                    hit_all_bounds_count += 1;
                }
                if (self.boundary_handler == BoundaryHandler::TotalBoundaryHits
                    && boundary_hits > self.max_boundary_hits)
                    || (self.boundary_handler == BoundaryHandler::HitsAllBoundaries
                        && hit_all_bounds_count > 2)
                {
                    self.termination_reason = OutOfBounds;
                    solution.copy_from_slice(&next_guess);
                    break;
                }
            }
            self.current_iteration += 1;
        }
        if sum > best_sum && best_sum.abs() < BIG_DOUBLE {
            sum = best_sum;
            solution.clone_from(&best_guess);
            let _ =
                self.evaluator
                    .evaluate(&solution, &mut equations, Some(&mut jacobian), &mut found);
        }
        if sum < self.desired_accuracy {
            return self.done(true, FoundAnswerConverged, solution);
        }
        if sum < self.acceptable_accuracy {
            return self.done(true, FoundAnswerClose, solution);
        }
        Ok(RootNd {
            found: false,
            x: solution,
        })
    }
}
