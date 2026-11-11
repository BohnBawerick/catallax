//! Parameter-coupled degree-seven UV fits, with outward Taylor bounds.

use catallax_numbers::{Ball, Decider, Rational, Refusal};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap};
use std::sync::{Arc, Condvar, Mutex, MutexGuard, PoisonError};

use crate::{fit, number};
use number::Number;

/// Taylor coefficients, not derivatives. Constant and polynomial jets retain
/// their degree so affine coordinate changes do not multiply known zero tails.
#[derive(Clone)]
pub(crate) struct Jet<const N: usize> {
    pub coefficients: [Ball; N],
    degree: usize,
}

fn scalar(value: i64, precision: u32) -> Ball {
    Ball::from_rational(&Rational::from_i64(value), precision)
}

impl<const N: usize> Jet<N> {
    pub fn constant(value: Ball) -> Self {
        let precision = value.precision();
        Self {
            coefficients: std::array::from_fn(|i| {
                if i == 0 {
                    value.clone()
                } else {
                    scalar(0, precision)
                }
            }),
            degree: 0,
        }
    }

    pub fn variable(value: Ball) -> Self {
        let precision = value.precision();
        let mut result = Self::constant(value);
        if N > 1 {
            result.coefficients[1] = scalar(1, precision);
            result.degree = 1;
        }
        result
    }

    pub fn from_coefficients(coefficients: [Ball; N]) -> Self {
        Self {
            coefficients,
            degree: N - 1,
        }
    }

    pub fn add(&self, other: &Self) -> Self {
        let mut result = Self::constant(scalar(0, self.coefficients[0].precision()));
        result.degree = self.degree.max(other.degree);
        for i in 0..=result.degree {
            result.coefficients[i] = &self.coefficients[i] + &other.coefficients[i];
        }
        result
    }

    pub fn sub(&self, other: &Self) -> Self {
        let mut result = Self::constant(scalar(0, self.coefficients[0].precision()));
        result.degree = self.degree.max(other.degree);
        for i in 0..=result.degree {
            result.coefficients[i] = &self.coefficients[i] - &other.coefficients[i];
        }
        result
    }

    pub fn scale(&self, factor: &Ball) -> Self {
        let mut result = Self::constant(scalar(0, self.coefficients[0].precision()));
        if factor.is_exact_zero() {
            return result;
        }
        result.degree = self.degree;
        for i in 0..=self.degree {
            result.coefficients[i] = &self.coefficients[i] * factor;
        }
        result
    }

    pub fn mul(&self, other: &Self) -> Self {
        let mut result = Self::constant(scalar(0, self.coefficients[0].precision()));
        result.degree = (self.degree + other.degree).min(N - 1);
        for n in 0..=result.degree {
            for k in n.saturating_sub(other.degree)..=n.min(self.degree) {
                result.coefficients[n] =
                    &result.coefficients[n] + &(&self.coefficients[k] * &other.coefficients[n - k]);
            }
        }
        result
    }

    /// The quotient series, for a divisor whose constant term is separated
    /// from zero.
    pub fn div(&self, other: &Self) -> Result<Self, Refusal> {
        let inverse = other.coefficients[0]
            .recip()
            .ok_or_else(|| number::budget("UV jet divisor is not separated from zero"))?;
        let mut result = Self::constant(&self.coefficients[0] * &inverse);
        if self.degree == 0 && other.degree == 0 {
            return Ok(result);
        }
        result.degree = N - 1;
        for n in 1..N {
            let mut residual = self.coefficients[n].clone();
            for k in 0..n {
                if n - k <= other.degree {
                    residual = residual - &result.coefficients[k] * &other.coefficients[n - k];
                }
            }
            result.coefficients[n] = &residual * &inverse;
        }
        Ok(result)
    }

    pub fn square(&self) -> Self {
        let precision = self.coefficients[0].precision();
        let mut result = Self::constant(self.coefficients[0].square());
        result.degree = (self.degree * 2).min(N - 1);
        for n in 1..=result.degree {
            for k in n.saturating_sub(self.degree)..=self.degree.min(n / 2) {
                let other = n - k;
                let term = if k == other {
                    self.coefficients[k].square()
                } else {
                    &(&self.coefficients[k] * &self.coefficients[other]) * &scalar(2, precision)
                };
                result.coefficients[n] = &result.coefficients[n] + &term;
            }
        }
        result
    }

    pub fn sqrt(&self) -> Result<Self, Refusal> {
        if self.coefficients[0].lo() <= Rational::zero() {
            return Err(number::budget("UV radial jet is not separated from a pole"));
        }
        let root = self.coefficients[0]
            .sqrt()
            .expect("positive UV radial square");
        let inverse = (&root * &scalar(2, root.precision()))
            .recip()
            .ok_or_else(|| number::budget("UV radial derivative is not separated"))?;
        let mut result = Self::constant(root);
        if self.degree == 0 {
            return Ok(result);
        }
        result.degree = N - 1;
        for n in 1..N {
            let mut residual = self.coefficients[n].clone();
            for k in 1..n {
                residual = residual - &result.coefficients[k] * &result.coefficients[n - k];
            }
            result.coefficients[n] = &residual * &inverse;
        }
        Ok(result)
    }

    /// Integrate d atan2(y,x) = (x dy - y dx)/(x²+y²). The caller selects
    /// the continuous angular lift for coefficient zero; derivatives do not
    /// depend on that lift and need no transcendental series for remainder work.
    pub fn angle(y: &Self, x: &Self, constant: Ball) -> Result<Self, Refusal> {
        let mut result = Self::constant(constant);
        if N == 1 || (x.degree == 0 && y.degree == 0) {
            return Ok(result);
        }
        let precision = result.coefficients[0].precision();
        let mut denominator: [Ball; N] = std::array::from_fn(|_| scalar(0, precision));
        let mut quotient: [Ball; N] = std::array::from_fn(|_| scalar(0, precision));
        let inverse = (x.coefficients[0].square() + y.coefficients[0].square())
            .recip()
            .ok_or_else(|| number::budget("UV angular jet is not separated from its axis"))?;
        for n in 0..N - 1 {
            let mut numerator = scalar(0, precision);
            for k in 0..=n {
                if n != 0 && k <= x.degree && n - k <= x.degree {
                    denominator[n] =
                        &denominator[n] + &(&x.coefficients[k] * &x.coefficients[n - k]);
                }
                if n != 0 && k <= y.degree && n - k <= y.degree {
                    denominator[n] =
                        &denominator[n] + &(&y.coefficients[k] * &y.coefficients[n - k]);
                }
                let order = n - k + 1;
                let factor = scalar(order as i64, precision);
                if k <= x.degree && order <= y.degree {
                    numerator =
                        numerator + &(&x.coefficients[k] * &y.coefficients[order]) * &factor;
                }
                if k <= y.degree && order <= x.degree {
                    numerator =
                        numerator - &(&y.coefficients[k] * &x.coefficients[order]) * &factor;
                }
            }
            for k in 1..=n {
                numerator = numerator - &denominator[k] * &quotient[n - k];
            }
            quotient[n] = &numerator * &inverse;
            result.coefficients[n + 1] = &quotient[n]
                * &Ball::from_rational(
                    &Rational::from_ratio(1, (n + 1) as i64).unwrap(),
                    precision,
                );
        }
        result.degree = N - 1;
        Ok(result)
    }
}

/// A local curve has an actual STEP parameter interval. Unlike the existing
/// space-curve fit, UV refinement must never reparameterize its Bezier spans.
pub(crate) trait Parametric {
    fn range(&self) -> [Rational; 2];
    fn reuse_endpoints(&self) -> bool {
        false
    }
    /// `false` requests interval subdivision before evaluating derivatives;
    /// arithmetic refusals still propagate without being swallowed or retried.
    fn regular(&self, _parameter: &Ball, _dec: &Decider) -> Result<bool, Refusal> {
        Ok(true)
    }
    fn coefficients<const N: usize>(
        &self,
        parameter: &Ball,
        constants: bool,
        dec: &Decider,
    ) -> Result<[[Ball; 2]; N], Refusal>;
    fn value(
        &self,
        parameter: &Rational,
        precision: u32,
        dec: &Decider,
    ) -> Result<[Ball; 2], Refusal> {
        Ok(
            self.coefficients::<1>(&Ball::from_rational(parameter, precision), true, dec)?[0]
                .clone(),
        )
    }
}

pub(crate) struct Fit {
    pub points: Vec<[Number; 2]>,
    pub knots: Vec<Rational>,
    pub multiplicities: Vec<usize>,
    pub error: Rational,
}

#[cfg(test)]
impl Fit {
    /// Exact degree-seven controls, joined spans sharing their end points,
    /// with unit knot spacing, for chart boundary tests.
    pub fn exact(points: Vec<[Rational; 2]>) -> Self {
        assert_eq!((points.len() - 1) % fit::DEGREE, 0, "whole spans");
        let spans = (points.len() - 1) / fit::DEGREE;
        Self {
            points: points
                .into_iter()
                .map(|point| {
                    point.map(|value| Number {
                        value,
                        text: String::new(),
                        error: Rational::zero(),
                    })
                })
                .collect(),
            knots: (0..=spans as i64).map(Rational::from_i64).collect(),
            multiplicities: (0..=spans)
                .map(|k| {
                    if k == 0 || k == spans {
                        fit::DEGREE + 1
                    } else {
                        fit::DEGREE
                    }
                })
                .collect(),
            error: Rational::zero(),
        }
    }

    /// An exact polyline of integer vertices over denominator `scale`, one
    /// straight span per segment.
    pub fn polyline(vertices: &[[i64; 2]], scale: i64) -> Self {
        let degree = fit::DEGREE as i64;
        let mut points = Vec::new();
        for (index, pair) in vertices.windows(2).enumerate() {
            for i in usize::from(index != 0)..=fit::DEGREE {
                let i = i as i64;
                points.push(std::array::from_fn(|axis| {
                    Rational::from_ratio(
                        pair[0][axis] * (degree - i) + pair[1][axis] * i,
                        degree * scale,
                    )
                    .unwrap()
                }));
            }
        }
        Self::exact(points)
    }
}

struct Span {
    range: [Rational; 2],
    points: Vec<[Number; 2]>,
    ideals: Vec<[[Rational; 2]; 2]>,
    remainder: Rational,
}

impl Span {
    fn error(&self) -> Rational {
        self.remainder.clone()
            + self
                .points
                .iter()
                .zip(&self.ideals)
                .map(|(p, bounds)| {
                    (0..2).fold(Rational::zero(), |sum, i| {
                        sum + (&p[i].value - &bounds[i][0])
                            .abs()
                            .max((&p[i].value - &bounds[i][1]).abs())
                    })
                })
                .max()
                .expect("nonempty UV controls")
    }
}

/// A span's Taylor controls and remainder, before the curve's own values at
/// the span's ends replace the end controls.
#[derive(Clone)]
struct Core {
    points: Vec<[Number; 2]>,
    ideals: Vec<[[Rational; 2]; 2]>,
    remainder: Rational,
}

/// The controls of `piece` over `range`, or `None` when the chart is not
/// regular there and the range must be cut before derivatives are taken.
/// Depends on nothing but its arguments.
fn core<P: Parametric>(
    piece: &P,
    range: &[Rational; 2],
    digits: u32,
    precision: u32,
    dec: &Decider,
) -> Result<Option<Core>, Refusal> {
    let [lo, hi] = range;
    let domain = Ball::from_rational(lo, precision).hull(&Ball::from_rational(hi, precision));
    if !piece.regular(&domain, dec)? {
        return Ok(None);
    }
    let width = hi - lo;
    let center = Ball::from_rational(&(lo + hi).mul_pow2(-1), precision);
    let coefficients = piece.coefficients::<{ fit::DEGREE + 1 }>(&center, true, dec)?;
    let remainder_coefficient = piece.coefficients::<{ fit::DEGREE + 2 }>(&domain, false, dec)?;
    let remainder = remainder_coefficient[fit::DEGREE + 1]
        .iter()
        .fold(Rational::zero(), |sum, value| {
            sum + number::magnitude(value)
        })
        * width.mul_pow2(-1).pow((fit::DEGREE + 1) as i32).unwrap();
    let mut points = Vec::with_capacity(fit::DEGREE + 1);
    let mut ideals = Vec::with_capacity(fit::DEGREE + 1);
    fit::taylor_controls(&coefficients, &width, digits, |point, ideal| {
        points.push(point);
        ideals.push(ideal);
    });
    Ok(Some(Core {
        points,
        ideals,
        remainder,
    }))
}

/// Whether a span must be cut in two: it misses the tolerance, and its
/// remainder is a large enough share of the miss for a cut to help.
fn needs_cut(span: &Span, tolerance: &Rational) -> bool {
    let error = span.error();
    let useful = &span.remainder * &Rational::from_i64(4) > error.clone().max(tolerance.clone());
    error > *tolerance && useful
}

fn cut(range: &[Rational; 2], digits: u32) -> Option<Rational> {
    let middle = number::rational(&(&range[0] + &range[1]).mul_pow2(-1), digits).value;
    (middle > range[0] && middle < range[1]).then_some(middle)
}

pub(crate) struct Budget {
    pub precision: u32,
    pieces: usize,
    cost: u64,
    remaining: u64,
}

impl Budget {
    pub fn new(pieces: usize, digits: u32, dec: &Decider) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        assert!(pieces != 0, "empty UV trace");
        // Two coordinate jets, two endpoints and Taylor-to-Bezier controls.
        // The request-wide arithmetic meter remains active as well.
        const OPERATIONS: u64 = 32_768;
        let count = u64::try_from(pieces).unwrap();
        let initial = count
            .checked_mul(OPERATIONS)
            .ok_or_else(|| number::budget("too many initial UV spans"))?;
        let precision = number::work_precision(
            dec.policy(),
            digits.max(1).saturating_mul(4).saturating_add(32),
            initial,
        )
        .ok_or_else(|| number::budget("no work available for an initial UV fit"))?;
        let cost = u64::from(precision.div_ceil(64)).pow(2) * OPERATIONS;
        Ok(Self {
            precision,
            pieces,
            cost,
            remaining: dec.policy().max_work - cost * count,
        })
    }
}

/// A fit runs alone this long before it asks for helper threads. Most fits
/// finish sooner, and starting threads costs more than they would save.
const ALONE: std::time::Duration = std::time::Duration::from_millis(2);

pub(crate) fn fit<P: Parametric + Sync>(
    pieces: &[P],
    digits: u32,
    tolerance: &Rational,
    budget: Budget,
    dec: &Decider,
) -> Result<Fit, Refusal> {
    dec.check_deadline()?;
    assert_eq!(pieces.len(), budget.pieces, "changed UV reservation");
    let Budget {
        precision,
        cost,
        remaining,
        ..
    } = budget;
    let mut walk = Walk {
        pieces,
        digits,
        tolerance,
        precision,
        cost,
        remaining,
        piece: 0,
        pending: vec![pieces[0].range()],
        endpoints: BTreeMap::new(),
        spans: Vec::new(),
    };
    // The walk below is the whole fit, in order, on this thread. Helper
    // threads only compute spans the walk is likely to reach; the walk adopts
    // a finished one with the work it charged, or computes the span itself.
    let began = std::time::Instant::now();
    let mut helpers = None;
    let mut asked = false;
    // A test that compares complete replays asks for helpers at once and lets
    // them finish first, so that the walk adopts every span.
    let wait = catallax_numbers::run_ahead_waits();
    while !walk.done() {
        if !asked && (wait || began.elapsed() >= ALONE) {
            asked = true;
            helpers = catallax_numbers::Helpers::new(dec);
            if helpers.is_some() {
                break;
            }
        }
        walk.step(None, dec)?;
    }
    if let Some(helpers) = helpers {
        let ahead = Ahead::new(&walk, &helpers);
        std::thread::scope(|scope| {
            let mut started = 0;
            for _ in 0..helpers.count() {
                let spawned = std::thread::Builder::new()
                    .name("catallax-uv-fit".into())
                    .stack_size(helpers.stack_size())
                    .spawn_scoped(scope, || helpers.run(|dec| ahead.explore(dec)));
                if spawned.is_err() {
                    break;
                }
                started += 1;
            }
            if wait && started != 0 {
                ahead.settle();
            }
            let walked = (|| {
                while !walk.done() {
                    walk.step(Some(&ahead), dec)?;
                }
                Ok::<_, Refusal>(())
            })();
            ahead.stop();
            walked
        })?;
    }
    let spans = walk.spans;
    let mut points: Vec<[Number; 2]> = Vec::with_capacity(spans.len() * fit::DEGREE + 1);
    let mut knots = Vec::with_capacity(spans.len() + 1);
    let mut error = Rational::zero();
    for mut span in spans {
        if let Some(previous) = points.last() {
            span.points[0] = previous.clone();
            assert_eq!(
                knots.last(),
                Some(&span.range[0]),
                "disjoint UV parameter pieces"
            );
        } else {
            knots.push(span.range[0].clone());
        }
        error = error.max(span.error());
        let skip = usize::from(!points.is_empty());
        points.extend(span.points.into_iter().skip(skip));
        let [_, hi] = span.range;
        knots.push(hi);
    }
    let multiplicities = (0..knots.len())
        .map(|i| {
            if i == 0 || i + 1 == knots.len() {
                fit::DEGREE + 1
            } else {
                fit::DEGREE
            }
        })
        .collect();
    Ok(Fit {
        points,
        knots,
        multiplicities,
        error,
    })
}

/// The fit's subdivision, depth first and left first, one piece after another.
struct Walk<'a, P> {
    pieces: &'a [P],
    digits: u32,
    tolerance: &'a Rational,
    precision: u32,
    cost: u64,
    remaining: u64,
    /// The piece being cut, its ranges still to do and its end values so far.
    piece: usize,
    pending: Vec<[Rational; 2]>,
    endpoints: BTreeMap<Rational, [Ball; 2]>,
    spans: Vec<Span>,
}

impl<P: Parametric + Sync> Walk<'_, P> {
    fn done(&self) -> bool {
        self.pending.is_empty() && self.piece + 1 >= self.pieces.len()
    }

    /// Does the next range: keeps its span, or cuts it in two.
    fn step(&mut self, ahead: Option<&Ahead<'_, P>>, dec: &Decider) -> Result<(), Refusal> {
        let range = loop {
            if let Some(range) = self.pending.pop() {
                break range;
            }
            self.piece += 1;
            self.endpoints = BTreeMap::new();
            self.pending.push(self.pieces[self.piece].range());
        };
        let (index, piece) = (self.piece, &self.pieces[self.piece]);
        dec.check_deadline()?;
        let core = match ahead {
            Some(ahead) => ahead.core(index, &range, dec)?,
            None => core(piece, &range, self.digits, self.precision, dec)?,
        };
        let mut current = None;
        let refine = if let Some(Core {
            mut points,
            ideals,
            remainder,
        }) = core
        {
            for (slot, parameter) in [(0, &range[0]), (fit::DEGREE, &range[1])] {
                let value = |dec: &Decider| match ahead {
                    Some(ahead) => ahead.value(index, parameter, dec),
                    None => piece.value(parameter, self.precision, dec),
                };
                let value = if piece.reuse_endpoints() {
                    dec.check_deadline()?;
                    match self.endpoints.entry(parameter.clone()) {
                        std::collections::btree_map::Entry::Occupied(entry) => entry.get().clone(),
                        std::collections::btree_map::Entry::Vacant(entry) => {
                            entry.insert(value(dec)?).clone()
                        }
                    }
                } else {
                    value(dec)?
                };
                points[slot] = value
                    .each_ref()
                    .map(|value| number::ball_number(value, self.digits));
            }
            let span = Span {
                range: range.clone(),
                points,
                ideals,
                remainder,
            };
            let refine = needs_cut(&span, self.tolerance);
            current = Some(span);
            refine
        } else {
            true
        };
        if refine
            && self.remaining >= self.cost.saturating_mul(2)
            && let Some(middle) = cut(&range, self.digits)
        {
            self.remaining -= self.cost * 2;
            self.pending.push([middle.clone(), range[1].clone()]);
            self.pending.push([range[0].clone(), middle]);
            return Ok(());
        }
        self.spans.push(
            current.ok_or_else(|| {
                number::budget("UV chart domain cannot be separated from its axis")
            })?,
        );
        Ok(())
    }
}

/// What a helper or the writer knows of one computation.
enum Slot<T> {
    /// A helper is computing it.
    Pending,
    /// The writer is computing it itself.
    Taken,
    /// A helper's result, with the work it charged.
    Recorded(Arc<(T, catallax_numbers::Charges)>),
    /// The writer's own result. Its work is already charged.
    Known(T),
    /// A helper's computation that left no usable record.
    Unusable,
}

type CoreResult = Result<Option<Core>, Refusal>;
type ValueResult = Result<[Ball; 2], Refusal>;

#[derive(Default)]
struct Explored {
    cores: BTreeMap<(usize, Rational, Rational), Slot<CoreResult>>,
    values: BTreeMap<(usize, Rational), Slot<ValueResult>>,
    /// Ranges to explore; the walk reaches the least first.
    queue: BinaryHeap<Reverse<(usize, Rational, Rational)>>,
    /// Helpers exploring a range now.
    busy: usize,
}

/// Spans computed ahead of the walk by helper threads.
struct Ahead<'a, P> {
    pieces: &'a [P],
    digits: u32,
    tolerance: &'a Rational,
    precision: u32,
    helpers: &'a catallax_numbers::Helpers,
    explored: Mutex<Explored>,
    changed: Condvar,
}

fn locked<T>(mutex: &Mutex<T>) -> MutexGuard<'_, T> {
    mutex.lock().unwrap_or_else(PoisonError::into_inner)
}

impl<'a, P: Parametric + Sync> Ahead<'a, P> {
    fn new(walk: &Walk<'a, P>, helpers: &'a catallax_numbers::Helpers) -> Self {
        let mut explored = Explored::default();
        for range in &walk.pending {
            explored
                .queue
                .push(Reverse((walk.piece, range[0].clone(), range[1].clone())));
        }
        for (index, piece) in walk.pieces.iter().enumerate().skip(walk.piece + 1) {
            let [lo, hi] = piece.range();
            explored.queue.push(Reverse((index, lo, hi)));
        }
        Self {
            pieces: walk.pieces,
            digits: walk.digits,
            tolerance: walk.tolerance,
            precision: walk.precision,
            helpers,
            explored: Mutex::new(explored),
            changed: Condvar::new(),
        }
    }

    /// Waits until the helpers have explored everything they can reach.
    fn settle(&self) {
        let mut explored = locked(&self.explored);
        while explored.busy != 0 || !explored.queue.is_empty() {
            explored = self.wait(explored);
        }
    }

    fn stop(&self) {
        self.helpers.stop();
        drop(locked(&self.explored));
        self.changed.notify_all();
    }

    fn wait<'g>(&self, guard: MutexGuard<'g, Explored>) -> MutexGuard<'g, Explored> {
        self.changed
            .wait(guard)
            .unwrap_or_else(PoisonError::into_inner)
    }

    /// The writer's span controls: a helper's, with its work charged here, or
    /// the writer's own.
    fn core(&self, index: usize, range: &[Rational; 2], dec: &Decider) -> CoreResult {
        let key = (index, range[0].clone(), range[1].clone());
        let cold = || core(&self.pieces[index], range, self.digits, self.precision, dec);
        let mut explored = locked(&self.explored);
        loop {
            match explored.cores.get(&key) {
                Some(Slot::Pending) => explored = self.wait(explored),
                Some(Slot::Recorded(recorded)) => {
                    let recorded = recorded.clone();
                    drop(explored);
                    return if recorded.1.adopt() {
                        recorded.0.clone()
                    } else {
                        cold()
                    };
                }
                None => {
                    explored.cores.insert(key.clone(), Slot::Taken);
                    drop(explored);
                    let result = cold();
                    locked(&self.explored)
                        .cores
                        .insert(key, Slot::Known(result.clone()));
                    self.changed.notify_all();
                    return result;
                }
                Some(Slot::Taken | Slot::Known(_) | Slot::Unusable) => {
                    drop(explored);
                    return cold();
                }
            }
        }
    }

    /// The writer's value of a piece at one parameter, likewise.
    fn value(&self, index: usize, parameter: &Rational, dec: &Decider) -> ValueResult {
        let key = (index, parameter.clone());
        let cold = || self.pieces[index].value(parameter, self.precision, dec);
        let mut explored = locked(&self.explored);
        loop {
            match explored.values.get(&key) {
                Some(Slot::Pending) => explored = self.wait(explored),
                Some(Slot::Recorded(recorded)) => {
                    let recorded = recorded.clone();
                    drop(explored);
                    return if recorded.1.adopt() {
                        recorded.0.clone()
                    } else {
                        cold()
                    };
                }
                None => {
                    explored.values.insert(key.clone(), Slot::Taken);
                    drop(explored);
                    let result = cold();
                    locked(&self.explored)
                        .values
                        .insert(key, Slot::Known(result.clone()));
                    self.changed.notify_all();
                    return result;
                }
                Some(Slot::Taken | Slot::Known(_) | Slot::Unusable) => {
                    drop(explored);
                    return cold();
                }
            }
        }
    }

    /// A helper thread: explores queued ranges until none is left or the
    /// writer stops it.
    fn explore(&self, dec: &Decider) {
        loop {
            let (index, lo, hi) = {
                let mut explored = locked(&self.explored);
                loop {
                    if self.helpers.stopped() {
                        return;
                    }
                    if let Some(Reverse(task)) = explored.queue.pop() {
                        explored.busy += 1;
                        break task;
                    }
                    if explored.busy == 0 {
                        return;
                    }
                    explored = self.wait(explored);
                }
            };
            // A panic in a span computed ahead is the writer's to meet.
            let range = [lo, hi];
            let halves = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                self.explore_range(index, &range, dec)
            }))
            .unwrap_or(None);
            let mut explored = locked(&self.explored);
            explored.busy -= 1;
            if let Some(middle) = halves {
                let [lo, hi] = range;
                explored.queue.push(Reverse((index, lo, middle.clone())));
                explored.queue.push(Reverse((index, middle, hi)));
            }
            drop(explored);
            self.changed.notify_all();
        }
    }

    /// Computes one range's span as the walk would, and returns where the
    /// walk would cut it, if it would.
    fn explore_range(
        &self,
        index: usize,
        range: &[Rational; 2],
        dec: &Decider,
    ) -> Option<Rational> {
        let piece = &self.pieces[index];
        let key = (index, range[0].clone(), range[1].clone());
        let core = self.helped(
            |explored| &mut explored.cores,
            key,
            || core(piece, range, self.digits, self.precision, dec),
        )?;
        let refine = match core {
            Err(_) => return None,
            Ok(None) => true,
            Ok(Some(Core {
                mut points,
                ideals,
                remainder,
            })) => {
                for (slot, parameter) in [(0, &range[0]), (fit::DEGREE, &range[1])] {
                    let value = self
                        .helped(
                            |explored| &mut explored.values,
                            (index, parameter.clone()),
                            || piece.value(parameter, self.precision, dec),
                        )?
                        .ok()?;
                    points[slot] = value
                        .each_ref()
                        .map(|value| number::ball_number(value, self.digits));
                }
                needs_cut(
                    &Span {
                        range: range.clone(),
                        points,
                        ideals,
                        remainder,
                    },
                    self.tolerance,
                )
            }
        };
        refine.then(|| cut(range, self.digits)).flatten()
    }

    /// A helper's result for one computation: its own, recorded for the
    /// writer, or one already known. `None` when no usable result exists.
    fn helped<K: Ord + Clone, T: Clone>(
        &self,
        slots: impl Fn(&mut Explored) -> &mut BTreeMap<K, Slot<T>>,
        key: K,
        compute: impl FnOnce() -> T,
    ) -> Option<T> {
        let mut explored = locked(&self.explored);
        loop {
            if self.helpers.stopped() {
                return None;
            }
            match slots(&mut explored).get(&key) {
                Some(Slot::Recorded(recorded)) => return Some(recorded.0.clone()),
                Some(Slot::Known(value)) => return Some(value.clone()),
                Some(Slot::Unusable) => return None,
                Some(Slot::Pending | Slot::Taken) => explored = self.wait(explored),
                None => break,
            }
        }
        slots(&mut explored).insert(key.clone(), Slot::Pending);
        drop(explored);
        // The claim is settled on every path out, including a panic.
        struct Claim<'c, K: Ord + Clone, T, S: Fn(&mut Explored) -> &mut BTreeMap<K, Slot<T>>> {
            explored: &'c Mutex<Explored>,
            changed: &'c Condvar,
            slots: &'c S,
            key: K,
            settled: Option<Slot<T>>,
        }
        impl<K: Ord + Clone, T, S: Fn(&mut Explored) -> &mut BTreeMap<K, Slot<T>>> Drop
            for Claim<'_, K, T, S>
        {
            fn drop(&mut self) {
                let slot = self.settled.take().unwrap_or(Slot::Unusable);
                (self.slots)(&mut locked(self.explored)).insert(self.key.clone(), slot);
                self.changed.notify_all();
            }
        }
        let mut claim = Claim {
            explored: &self.explored,
            changed: &self.changed,
            slots: &slots,
            key,
            settled: None,
        };
        let (value, charges) = catallax_numbers::record(compute);
        let charges = charges?;
        claim.settled = Some(Slot::Recorded(Arc::new((value.clone(), charges))));
        Some(value)
    }
}

pub(crate) type Atan = fn(&Ball, &Decider) -> Result<Ball, Refusal>;

/// A locally continuous atan2 branch. Periodic callers choose its integer lift.
pub(crate) fn angle(
    y: &Ball,
    x: &Ball,
    pi: &Ball,
    atan: Atan,
    dec: &Decider,
) -> Result<Ball, Refusal> {
    dec.check_deadline()?;
    if x.lo() > Rational::zero() {
        return atan(&y.checked_div(x).unwrap(), dec);
    }
    if x.hi() < Rational::zero() {
        return Ok(pi - &atan(&y.checked_div(&(-x)).unwrap(), dec)?);
    }
    if y.lo() > Rational::zero() {
        return Ok(
            pi * &scalar_ratio(1, 2, pi.precision()) - atan(&x.checked_div(y).unwrap(), dec)?
        );
    }
    if y.hi() < Rational::zero() {
        return Ok(
            -(pi * &scalar_ratio(1, 2, pi.precision())) - atan(&x.checked_div(y).unwrap(), dec)?
        );
    }
    Err(number::budget("UV angle is not separated from its axis"))
}

pub(crate) fn positive_angle(
    y: &Ball,
    x: &Ball,
    pi: &Ball,
    atan: Atan,
    dec: &Decider,
) -> Result<Ball, Refusal> {
    dec.check_deadline()?;
    if x.lo() > Rational::zero() {
        if y.lo() >= Rational::zero() {
            return atan(&y.checked_div(x).unwrap(), dec);
        }
        if y.hi() < Rational::zero() {
            return Ok(pi * &scalar(2, pi.precision()) + atan(&y.checked_div(x).unwrap(), dec)?);
        }
    }
    if x.hi() < Rational::zero() {
        return Ok(pi - &atan(&y.checked_div(&(-x)).unwrap(), dec)?);
    }
    if y.lo() > Rational::zero() {
        return Ok(
            pi * &scalar_ratio(1, 2, pi.precision()) - atan(&x.checked_div(y).unwrap(), dec)?
        );
    }
    if y.hi() < Rational::zero() {
        return Ok(
            pi * &scalar_ratio(3, 2, pi.precision()) + atan(&x.checked_div(&(-y)).unwrap(), dec)?
        );
    }
    Err(number::budget(
        "UV angular lift is not separated from its cut",
    ))
}

fn scalar_ratio(numerator: i64, denominator: i64, precision: u32) -> Ball {
    Ball::from_rational(
        &Rational::from_ratio(numerator, denominator).unwrap(),
        precision,
    )
}
