//! Budgeted interval filtering, refinement, and exact algebraic zero proofs.

use crate::{
    Algebraic, Ball, BoxRoots, Coincidence, MPoly, NotExact, Outcome, PolySystem, Rational, Real,
    RealRoot, Refusal, RootBox, Sign, StrictOrder, StrictSign, UPoly, ZeroProof,
};
use crate::{
    algebraic::{self, AlgebraicRepr, RootEnclosure},
    real::Expr,
    system,
};
use std::{
    collections::{HashMap, VecDeque},
    rc::Rc,
    sync::Arc,
    time::Instant,
};

// Only test builds replace the clock. Normal builds read Instant directly.
fn deadline_passed(limit: Instant) -> bool {
    #[cfg(any(test, feature = "test-deadline"))]
    if let Some(now) = test_deadline::now(limit) {
        return now >= limit;
    }
    Instant::now() >= limit
}

#[cfg(any(test, feature = "test-deadline"))]
pub(crate) mod test_deadline {
    use std::{
        cell::Cell,
        time::{Duration, Instant},
    };

    #[derive(Clone, Copy, PartialEq, Eq)]
    pub enum Point {
        Enclosure,
        Work,
        Definedness,
    }

    #[derive(Clone, Copy)]
    struct Clock {
        point: Point,
        remaining: usize,
    }

    thread_local! {
        static CLOCK: Cell<Option<Clock>> = Cell::new(process_clock());
    }

    fn process_clock() -> Option<Clock> {
        #[cfg(feature = "test-deadline")]
        if let Ok(value) = std::env::var("CATALLAX_TEST_DEADLINE_CHECKS") {
            let checks = value.parse::<usize>().expect("test deadline check count");
            assert!(checks > 1);
            return Some(Clock {
                point: Point::Work,
                remaining: checks,
            });
        }
        None
    }

    // Freeze time until a chosen inner loop has made progress, then cross the
    // deadline. Thread-local state and RAII keep parallel tests and panics isolated.
    #[cfg(test)]
    pub struct Deadline;
    #[cfg(test)]
    impl Deadline {
        pub fn after(point: Point, checks: usize) -> Self {
            assert!(checks > 1);
            CLOCK.with(|slot| {
                assert!(slot.get().is_none());
                slot.set(Some(Clock {
                    point,
                    remaining: checks,
                }));
            });
            Self
        }

        pub fn instant(&self) -> Instant {
            Instant::now()
        }

        pub fn assert_expired(&self) {
            CLOCK.with(|slot| assert_eq!(slot.get().unwrap().remaining, 0));
        }
    }
    #[cfg(test)]
    impl Drop for Deadline {
        fn drop(&mut self) {
            CLOCK.with(|slot| slot.set(None));
        }
    }

    pub(super) fn now(limit: Instant) -> Option<Instant> {
        CLOCK.with(|slot| {
            slot.get().map(|clock| {
                if clock.remaining == 0 {
                    limit
                } else {
                    limit - Duration::from_nanos(1)
                }
            })
        })
    }

    pub(super) fn checkpoint(point: Point) {
        CLOCK.with(|slot| {
            if let Some(mut clock) = slot.get()
                && clock.point == point
                && clock.remaining > 0
            {
                clock.remaining -= 1;
                slot.set(Some(clock));
            }
        });
    }
}

/// Per-decision deterministic work and working-precision limits.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct BudgetPolicy {
    /// Maximum working precision during refinement.
    pub max_precision_bits: u32,
    /// Maximum deterministic work units.
    pub max_work: u64,
}
impl Default for BudgetPolicy {
    fn default() -> Self {
        Self {
            max_precision_bits: 16384,
            max_work: 1_000_000_000,
        }
    }
}
/// Cumulative decision counts, work, and highest working precision.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct DecisionStats {
    /// Total decisions.
    pub decisions: u64,
    /// Proven by the cached interval filter.
    pub filter: u64,
    /// Proven after refinement.
    pub refined: u64,
    /// Proven nonzero by exact arithmetic.
    pub exact: u64,
    /// Zeros proven by rational or algebraic arithmetic.
    pub coincident_exact: u64,
    /// Zeros proven by separation.
    pub coincident_separation: u64,
    /// Decisions settled by construction facts.
    pub coincident_fact: u64,
    /// Bounded but unresolved decisions.
    pub not_exact: u64,
    /// Refusals.
    pub refused: u64,
    /// Highest working precision used.
    pub max_bits: u32,
    /// Total deterministic work.
    pub work: u64,
}
const ROOT_CACHE_ENTRIES: usize = 128;

#[derive(Debug)]
struct CachedRoots {
    roots: Vec<RealRoot>,
    work: u64,
}

/// Certified decisions with fresh per-call budgets and cumulative statistics.
#[derive(Debug)]
pub struct Decider {
    policy: BudgetPolicy,
    deadline: Option<Instant>,
    stats: DecisionStats,
    root_cache: HashMap<Arc<UPoly>, CachedRoots>,
    root_cache_order: VecDeque<Arc<UPoly>>,
}
impl Default for Decider {
    fn default() -> Self {
        Self::new(BudgetPolicy::default())
    }
}
impl Decider {
    /// Creates a decider without a wall-clock deadline.
    pub fn new(policy: BudgetPolicy) -> Self {
        Self {
            policy,
            deadline: None,
            stats: DecisionStats::default(),
            root_cache: HashMap::new(),
            root_cache_order: VecDeque::new(),
        }
    }
    /// Current policy.
    pub fn policy(&self) -> &BudgetPolicy {
        &self.policy
    }
    /// Sets or removes the operation deadline.
    pub fn set_deadline(&mut self, deadline: Option<Instant>) {
        self.deadline = deadline;
    }
    /// Certified sign, exact coincidence, bounded uncertainty, or named refusal.
    pub fn sign(&mut self, x: &Real) -> Outcome<StrictSign> {
        self.stats.decisions += 1;
        if self.expired() {
            return self.refused(deadline());
        }
        let mut work = Work::new(self.policy, self.deadline);
        let proof_counts = (self.stats.filter, self.stats.refined, self.stats.exact);
        let result = self.sign_inner(x, &mut work);
        self.account(&work);
        if self.expired() {
            (self.stats.filter, self.stats.refined, self.stats.exact) = proof_counts;
            return self.refused(deadline());
        }
        match &result {
            Outcome::Proven(_) => {}
            Outcome::ProvenCoincident(c) => {
                if c.proof == ZeroProof::SeparationBound {
                    self.stats.coincident_separation += 1;
                } else {
                    self.stats.coincident_exact += 1;
                }
            }
            Outcome::NotExact(_) => self.stats.not_exact += 1,
            Outcome::Refused(_) => self.stats.refused += 1,
        }
        result
    }
    /// Exactly the decision and accounting of `sign(a-b)`.
    pub fn compare(&mut self, a: &Real, b: &Real) -> Outcome<StrictOrder> {
        self.sign(&(a - b)).map(|s| match s {
            StrictSign::Negative => StrictOrder::Less,
            StrictSign::Positive => StrictOrder::Greater,
        })
    }
    /// Attempts an exact rational coefficient reduction. `None` is inconclusive,
    /// not a proof of irrationality or definedness. Domain and resource failures
    /// retain their numeric reasons; original expression definitions are unchanged.
    pub fn rational_value(&mut self, x: &Real) -> Result<Option<Rational>, Refusal> {
        self.stats.decisions += 1;
        if self.expired() {
            self.stats.refused += 1;
            return Err(deadline());
        }
        let mut work = Work::new(self.policy, self.deadline);
        let result = (|| {
            if let Some(value) = x.rational_ref() {
                work.exact(value, value, false)
                    .map_err(|code| Refusal::new(code, "rational coefficient reduction stopped"))?;
                return Ok(Some(value.clone()));
            }
            if !x.node.fast_defined {
                match evaluate(
                    x,
                    64u32.min(self.policy.max_precision_bits),
                    &mut work,
                    false,
                ) {
                    Ok(_) => {}
                    Err(EvalError::Pending(_)) => return Ok(None),
                    Err(EvalError::Domain(error)) => return Err(error),
                    Err(EvalError::Limit(code)) => {
                        return Err(Refusal::new(
                            code,
                            "rational coefficient domain proof stopped",
                        ));
                    }
                }
            }
            let reduction = crate::identity::rational(x, &mut work).and_then(|value| match value {
                Some(value) => Ok(Some(value)),
                None => crate::tower::rational(x, &mut work).and_then(|value| match value {
                    Some(value) => Ok(Some(value)),
                    None => crate::generator::rational(x, &mut work),
                }),
            });
            reduction.map_err(|code| Refusal::new(code, "rational coefficient reduction stopped"))
        })();
        self.account(&work);
        if self.expired() {
            self.stats.refused += 1;
            return Err(deadline());
        }
        match result {
            Ok(Some(value)) => {
                self.stats.exact += 1;
                Ok(Some(value))
            }
            Ok(None) => Ok(None),
            Err(error) => {
                self.stats.refused += 1;
                Err(error)
            }
        }
    }
    pub(crate) fn angle_degrees(&mut self, angle: &Rational) -> Result<[Real; 2], Refusal> {
        self.stats.decisions += 1;
        let mut work = Work::new(self.policy, self.deadline);
        let result = crate::angle::construct(angle, &mut work);
        self.account(&work);
        self.check_deadline()?;
        match result {
            Ok(values) => {
                self.stats.exact += 1;
                Ok(values)
            }
            Err(code) => {
                self.stats.refused += 1;
                Err(Refusal::new(
                    code,
                    "rational-degree angle construction stopped",
                ))
            }
        }
    }
    pub(crate) fn real_poly_eliminant(
        &mut self,
        polynomial: &crate::RealPoly,
    ) -> Result<UPoly, Refusal> {
        self.stats.decisions += 1;
        let mut work = Work::new(self.policy, self.deadline);
        let result = crate::real_poly_roots::eliminant(polynomial, &mut work);
        self.account(&work);
        self.check_deadline()?;
        match result {
            Ok(polynomial) => {
                self.stats.exact += 1;
                Ok(polynomial)
            }
            Err(code) => {
                self.stats.refused += 1;
                Err(Refusal::new(
                    code,
                    "algebraic coefficient elimination stopped",
                ))
            }
        }
    }
    pub(crate) fn real_poly_resultant(
        &mut self,
        left: &crate::RealPoly,
        right: &crate::RealPoly,
        var: usize,
    ) -> Result<crate::RealPoly, Refusal> {
        self.stats.decisions += 1;
        let mut work = Work::new(self.policy, self.deadline);
        let result = crate::real_poly_roots::resultant(left, right, var, &mut work);
        self.account(&work);
        self.check_deadline()?;
        match result {
            Ok(polynomial) => {
                self.stats.exact += 1;
                Ok(polynomial)
            }
            Err(code) => {
                self.stats.refused += 1;
                Err(Refusal::new(
                    code,
                    "algebraic resultant construction stopped",
                ))
            }
        }
    }
    /// Refines to the requested width, or returns the tightest attained ball.
    pub fn enclose(&mut self, x: &Real, max_width: &Rational) -> Result<Ball, Refusal> {
        assert!(max_width.sign() == Sign::Positive);
        self.stats.decisions += 1;
        if self.expired() {
            self.stats.refused += 1;
            return Err(deadline());
        }
        let mut work = Work::new(self.policy, self.deadline);
        let mut best = x.node.coarse.clone();
        let mut defined = x.node.fast_defined;
        let mut p = 64u32.min(self.policy.max_precision_bits);
        let mut error = None;
        loop {
            if p == 0 {
                break;
            }
            match evaluate(x, p, &mut work, false) {
                Ok(b) => {
                    defined = true;
                    best = Some(tighter(best, b));
                    if best.as_ref().unwrap().rad().mul_pow2(1) <= *max_width {
                        break;
                    }
                }
                Err(EvalError::Domain(r)) => {
                    defined = false;
                    error = Some(r);
                    break;
                }
                Err(EvalError::Pending(r)) => {
                    error = Some(r);
                }
                Err(EvalError::Limit(_)) => break,
            }
            if p == self.policy.max_precision_bits {
                break;
            }
            p = p.saturating_mul(2).min(self.policy.max_precision_bits);
        }
        self.account(&work);
        if self.expired() {
            self.stats.refused += 1;
            return Err(deadline());
        }
        if defined && let Some(best) = best {
            if best.rad().mul_pow2(1) <= *max_width {
                self.stats.refined += 1;
            } else {
                self.stats.not_exact += 1;
            }
            Ok(best)
        } else {
            self.stats.refused += 1;
            Err(error
                .or_else(|| unresolved(x, self.deadline))
                .unwrap_or_else(|| {
                    Refusal::new(
                        "numbers.divisor-not-separated",
                        "no finite defined enclosure within budget",
                    )
                }))
        }
    }
    /// Encloses a known-defined value at one working precision, without making
    /// a decision or changing statistics. A divisor containing zero requires a
    /// higher precision; square roots clip their arguments' nonnegative parts.
    pub fn enclose_at(&self, x: &Real, precision: u32) -> Result<Ball, Refusal> {
        self.check_deadline()?;
        let precision = precision.max(1);
        let policy = BudgetPolicy {
            max_precision_bits: precision,
            max_work: u64::MAX,
        };
        let mut work = Work::new(policy, self.deadline);
        let result = evaluate(x, precision, &mut work, true).map_err(|error| match error {
            EvalError::Domain(refusal) | EvalError::Pending(refusal) => refusal,
            EvalError::Limit(code) => Refusal::new(code, "working-precision enclosure stopped"),
        });
        self.check_deadline()?;
        result
    }
    /// Encloses a known-defined value under the caller's numeric limits.
    /// Unlike `enclose_at`, this charges proof work and rejects a requested
    /// working precision above the caller's policy. It creates no default policy.
    pub fn enclose_budgeted_at(&mut self, x: &Real, precision: u32) -> Result<Ball, Refusal> {
        self.stats.decisions += 1;
        let precision = precision.max(1);
        let policy = BudgetPolicy {
            max_precision_bits: precision.min(self.policy.max_precision_bits),
            max_work: self.policy.max_work,
        };
        let mut work = Work::new(policy, self.deadline);
        let result = (|| {
            self.check_deadline()?;
            if precision > self.policy.max_precision_bits {
                return Err(Refusal::new(
                    "numbers.precision-limit",
                    "requested working precision exceeds the caller's policy",
                ));
            }
            evaluate(x, precision, &mut work, true).map_err(|error| match error {
                EvalError::Domain(refusal) | EvalError::Pending(refusal) => refusal,
                EvalError::Limit(code) => Refusal::new(code, "working-precision enclosure stopped"),
            })
        })();
        self.account(&work);
        if self.expired() {
            self.stats.refused += 1;
            return Err(deadline());
        }
        match result {
            Ok(ball) => {
                self.stats.refined += 1;
                Ok(ball)
            }
            Err(error) => {
                self.stats.refused += 1;
                Err(error)
            }
        }
    }
    /// All real roots in increasing order, with exact rational roots and multiplicities.
    pub fn real_roots(&mut self, p: &UPoly) -> Outcome<Vec<RealRoot>> {
        self.stats.decisions += 1;
        if self.expired() {
            return self.refused(deadline());
        }
        if p.is_zero() {
            return self.refused(Refusal::new(
                "numbers.zero-polynomial",
                "every real number is a root",
            ));
        }
        let mut work = Work::new(self.policy, self.deadline);
        // Linear roots need only one rational division; memoize the nonlinear
        // isolation certificates that geometry predicates repeatedly request.
        let memoize = p.degree().is_some_and(|degree| degree > 1);
        let cached = if memoize {
            self.root_cache.get(p)
        } else {
            None
        };
        let cache_hit = cached.is_some();
        let result = match cached {
            // Certificate reuse must not change per-decision work or budgets.
            Some(cached) => work.charge(cached.work).map(|()| cached.roots.clone()),
            None => algebraic::isolate(p, &mut work),
        };
        self.account(&work);
        if self.expired() {
            return self.refused(deadline());
        }
        match result {
            Ok(roots) => {
                if memoize && !cache_hit {
                    if self.root_cache.len() == ROOT_CACHE_ENTRIES {
                        let oldest = self.root_cache_order.pop_front().unwrap();
                        self.root_cache.remove(oldest.as_ref());
                    }
                    let key = Arc::new(p.clone());
                    self.root_cache_order.push_back(key.clone());
                    self.root_cache.insert(
                        key,
                        CachedRoots {
                            roots: roots.clone(),
                            work: work.used,
                        },
                    );
                    if self.expired() {
                        return self.refused(deadline());
                    }
                }
                self.stats.exact += 1;
                Outcome::Proven(roots)
            }
            Err(_) => self.refused(Refusal::new(
                "numbers.budget-exhausted",
                "root isolation exhausted its budget",
            )),
        }
    }
    /// Certifies a unique nonsingular root or absence of roots in a closed box.
    pub fn certify_root(
        &mut self,
        system: &Arc<PolySystem>,
        candidate: &RootBox,
    ) -> Outcome<BoxRoots> {
        assert_eq!(system.dimension(), candidate.dimension());
        self.stats.decisions += 1;
        if self.expired() {
            return self.refused(deadline());
        }
        let mut work = Work::new(self.policy, self.deadline);
        let mut p = 64u32.min(self.policy.max_precision_bits);
        let mut answer = None;
        while p != 0 {
            match system::certify(system, candidate, p, &mut work) {
                Ok(Some(v)) => {
                    answer = Some(v);
                    break;
                }
                Ok(None) => {}
                Err(_) => break,
            }
            if p == self.policy.max_precision_bits {
                break;
            }
            p = p.saturating_mul(2).min(self.policy.max_precision_bits);
        }
        self.account(&work);
        if self.expired() {
            return self.refused(deadline());
        }
        match answer {
            Some(v) => {
                self.stats.exact += 1;
                Outcome::Proven(v)
            }
            None => self.refused(Refusal::new(
                "numbers.root-not-certified",
                "neither a unique nonsingular root nor an empty box was certified",
            )),
        }
    }
    /// Counts a fact-based decision made by another piece.
    pub fn note_fact_decision(&mut self) {
        self.stats.decisions += 1;
        self.stats.coincident_fact += 1;
    }
    /// Cumulative statistics.
    pub fn stats(&self) -> &DecisionStats {
        &self.stats
    }
    /// Takes and resets all statistics, retaining policy and deadline.
    pub fn take_stats(&mut self) -> DecisionStats {
        std::mem::take(&mut self.stats)
    }
    /// Refuses work, including a decision in progress, once its deadline passes.
    pub fn check_deadline(&self) -> Result<(), Refusal> {
        #[cfg(any(test, feature = "test-deadline"))]
        test_deadline::checkpoint(test_deadline::Point::Enclosure);
        if self.expired() {
            Err(deadline())
        } else {
            Ok(())
        }
    }
    fn expired(&self) -> bool {
        crate::operation::charge(1) || self.deadline.is_some_and(deadline_passed)
    }
    fn refused<T>(&mut self, r: Refusal) -> Outcome<T> {
        self.stats.refused += 1;
        Outcome::Refused(r)
    }
    fn account(&mut self, w: &Work) {
        self.stats.work = self.stats.work.saturating_add(w.used);
        self.stats.max_bits = self.stats.max_bits.max(w.bits);
    }
    fn sign_inner(&mut self, x: &Real, work: &mut Work) -> Outcome<StrictSign> {
        let mut stopped = work.charge(1).err();
        if stopped.is_none()
            && x.node.fast_defined
            && let Some(s) = x.approx().sign()
        {
            self.stats.filter += 1;
            return Outcome::Proven(s);
        }
        if stopped.is_none()
            && x.node.fast_defined
            && let Expr::Algebraic(value) = &x.node.expr
            && let AlgebraicRepr::Root(root) = &value.repr
        {
            // The root lies strictly inside its certified interval. A zero
            // endpoint still proves a strict sign without changing its enclosure.
            let sign = if root.lo.sign() != Sign::Negative {
                Some(StrictSign::Positive)
            } else if root.hi.sign() != Sign::Positive {
                Some(StrictSign::Negative)
            } else {
                None
            };
            if let Some(sign) = sign {
                self.stats.exact += 1;
                return Outcome::Proven(sign);
            }
        }
        if stopped.is_none()
            && let Some(r) = &x.node.rational
        {
            if let Some(refusal) = unresolved(x, self.deadline)
                && (refusal.code == "numbers.division-by-zero"
                    || refusal.code == "numbers.sqrt-of-negative")
            {
                return Outcome::Refused(refusal);
            }
            if work.exact(r, &Rational::zero(), false).is_ok() {
                match r.sign() {
                    Sign::Zero => return coincident(ZeroProof::ExactArithmetic),
                    Sign::Negative => {
                        self.stats.exact += 1;
                        return Outcome::Proven(StrictSign::Negative);
                    }
                    Sign::Positive => {
                        self.stats.exact += 1;
                        return Outcome::Proven(StrictSign::Positive);
                    }
                }
            } else {
                stopped = Some("numbers.budget-exhausted");
            }
        }
        // A singleton certificate is sufficient once its entire definition
        // is known valid; an undefined factor must never disappear into zero.
        if stopped.is_none()
            && x.node.fast_defined
            && x.node.coarse.as_ref().is_some_and(Ball::is_exact_zero)
        {
            return coincident(ZeroProof::Algebraic);
        }
        let mut best = x.node.coarse.clone();
        let mut defined = x.node.fast_defined;
        let mut pending = None;
        // Bounds and exact identities are independent of interval precision.
        let mut separation = None;
        let mut tried_exact = false;
        let mut p = 64u32.min(work.policy.max_precision_bits);
        while stopped.is_none() {
            if p == 0 {
                stopped = Some("numbers.precision-limit");
                break;
            }
            match evaluate(x, p, work, false) {
                Ok(b) => {
                    defined = true;
                    best = Some(tighter(best, b.clone()));
                    pending = None;
                    if best.as_ref().is_some_and(Ball::is_exact_zero) {
                        return coincident(ZeroProof::Algebraic);
                    }
                    if let Some(s) = b.sign() {
                        self.stats.refined += 1;
                        return Outcome::Proven(s);
                    }
                    if !tried_exact {
                        tried_exact = true;
                        let exact =
                            crate::identity::rational(x, work).and_then(|value| match value {
                                Some(value) => Ok(Some(value.sign())),
                                None => {
                                    crate::tower::sign(x, work).and_then(|sign| match sign {
                                        Some(s) => Ok(Some(s)),
                                        None => crate::quadratic::sign(x, work).and_then(|sign| {
                                            match sign {
                                                Some(s) => Ok(Some(s)),
                                                None => crate::generator::rational(x, work)
                                                    .map(|r| r.map(|r| r.sign())),
                                            }
                                        }),
                                    })
                                }
                            });
                        match exact {
                            Ok(Some(Sign::Zero)) => return coincident(ZeroProof::Algebraic),
                            Ok(Some(s)) => {
                                self.stats.exact += 1;
                                return Outcome::Proven(if s == Sign::Negative {
                                    StrictSign::Negative
                                } else {
                                    StrictSign::Positive
                                });
                            }
                            Ok(None) | Err(_) => {}
                        }
                    }
                    if separation.is_none() {
                        // A failed reservation spends no work. Keep the work
                        // already performed, but never retry a failed attempt.
                        separation = Some(separation_bound(x, work).unwrap_or(None));
                    }
                    if separation
                        .as_ref()
                        .and_then(Option::as_ref)
                        .is_some_and(|bound| b.magnitude() < *bound)
                    {
                        return coincident(ZeroProof::SeparationBound);
                    }
                }
                Err(EvalError::Domain(r)) => return Outcome::Refused(r),
                Err(EvalError::Pending(r)) => pending = Some(r),
                Err(EvalError::Limit(code)) => {
                    stopped = Some(code);
                    break;
                }
            }
            if p == work.policy.max_precision_bits {
                stopped = Some("numbers.precision-limit");
                break;
            }
            p = p.saturating_mul(2).min(work.policy.max_precision_bits);
        }
        if !defined {
            return Outcome::Refused(
                pending
                    .or_else(|| unresolved(x, self.deadline))
                    .unwrap_or_else(|| {
                        Refusal::new(
                            "numbers.divisor-not-separated",
                            "definedness was not established",
                        )
                    }),
            );
        }
        if let Some(b) = best {
            Outcome::NotExact(NotExact {
                bound: b.magnitude(),
                code: stopped.unwrap_or("numbers.budget-exhausted"),
                detail: "certified enclosure retained at the first exhausted limit".into(),
            })
        } else {
            Outcome::Refused(
                pending
                    .or_else(|| unresolved(x, self.deadline))
                    .unwrap_or_else(|| {
                        Refusal::new(
                            "numbers.divisor-not-separated",
                            "no finite enclosure within budget",
                        )
                    }),
            )
        }
    }
}
fn deadline() -> Refusal {
    if crate::operation::exhausted() {
        return Refusal::new(
            "numbers.operation-budget",
            "operation work allowance exhausted",
        );
    }
    Refusal::new("numbers.deadline", "operation deadline passed")
}
fn coincident<T>(proof: ZeroProof) -> Outcome<T> {
    Outcome::ProvenCoincident(Coincidence { proof })
}
fn tighter(old: Option<Ball>, b: Ball) -> Ball {
    match old {
        Some(a) if a.rad() < b.rad() => a,
        _ => b,
    }
}

// Logical arithmetic costs. A rational's size includes numerator and denominator
// limbs. Ball add/sub cost 4L, multiply/divide/sqrt cost 8L^2, L=ceil(p/64).
// Polynomial division/GCD reserve dense coefficient-operation bounds before
// invoking exact routines. These conservative reservations are deterministic;
// failed reservations spend no work. Construction-time caches spend no decision
// budget and a cached interval decision costs one interval operation.
// Refinement memoizes successful evaluations of shared expression nodes at one
// precision, and annihilator construction memoizes their exact polynomials.
// Cache hits do no arithmetic work; errors are never cached. Node addresses are
// used only for lookups, never iteration or a numeric decision.
pub(crate) struct Work {
    pub policy: BudgetPolicy,
    deadline: Option<Instant>,
    pub used: u64,
    pub bits: u32,
}
impl Work {
    pub fn new(policy: BudgetPolicy, deadline: Option<Instant>) -> Self {
        Self {
            policy,
            deadline,
            used: 0,
            bits: 0,
        }
    }
    pub fn check_deadline(&self) -> Result<(), &'static str> {
        #[cfg(any(test, feature = "test-deadline"))]
        test_deadline::checkpoint(test_deadline::Point::Work);
        if crate::operation::charge(1) {
            Err("numbers.operation-budget")
        } else if self.deadline.is_some_and(deadline_passed) {
            Err("numbers.deadline")
        } else {
            Ok(())
        }
    }
    pub fn charge(&mut self, n: u64) -> Result<(), &'static str> {
        self.check_deadline()?;
        if n > self.policy.max_work.saturating_sub(self.used) {
            Err("numbers.budget-exhausted")
        } else {
            if crate::operation::charge(n) {
                return Err("numbers.operation-budget");
            }
            self.used += n;
            Ok(())
        }
    }
    pub fn exact(&mut self, a: &Rational, b: &Rational, mul: bool) -> Result<(), &'static str> {
        let a = limbs(a);
        let b = limbs(b);
        self.charge(if mul {
            a.saturating_mul(b)
        } else {
            a.saturating_add(b)
        })
    }
    pub fn ball(&mut self, p: u32, mul: bool) -> Result<(), &'static str> {
        if p > self.policy.max_precision_bits {
            return Err("numbers.precision-limit");
        }
        let l = u64::from(p).div_ceil(64).max(1);
        self.charge(if mul {
            l.saturating_mul(l).saturating_mul(8)
        } else {
            l.saturating_mul(4)
        })?;
        self.bits = self.bits.max(p);
        Ok(())
    }
    pub fn poly(&mut self, a: &UPoly, b: &UPoly) -> Result<(), &'static str> {
        let l = a
            .coefficients()
            .iter()
            .chain(b.coefficients())
            .map(limbs)
            .max()
            .unwrap_or(1);
        let n = (a.coefficients().len() + b.coefficients().len()).max(1) as u64;
        self.charge(
            n.saturating_mul(n)
                .saturating_mul(l)
                .saturating_mul(l)
                .saturating_mul(4),
        )
    }
    pub fn eval_poly(&mut self, p: &UPoly, x: &Rational) -> Result<Rational, &'static str> {
        let mut value = Rational::zero();
        for c in p.coefficients().iter().rev() {
            self.exact(&value, x, true)?;
            value = &value * x;
            self.exact(&value, c, false)?;
            value = &value + c;
        }
        Ok(value)
    }
    pub fn mpoly(&mut self, a: &MPoly, b: &MPoly) -> Result<(), &'static str> {
        let at = a.terms();
        let bt = b.terms();
        let l = at
            .iter()
            .chain(&bt)
            .map(|(_, c)| limbs(c))
            .max()
            .unwrap_or(1);
        let d = u64::from(a.total_degree().unwrap_or(0))
            .saturating_add(u64::from(b.total_degree().unwrap_or(0)))
            .max(1);
        self.charge(
            (at.len().max(1) as u64)
                .saturating_mul(bt.len().max(1) as u64)
                .saturating_mul(d.saturating_pow(3))
                .saturating_mul(l)
                .saturating_mul(l),
        )
    }
    pub fn system(&mut self, s: &PolySystem, p: u32) -> Result<(), &'static str> {
        self.ball(p, true)?;
        let n = s.dimension() as u64;
        let terms = s
            .equations()
            .iter()
            .map(|f| f.terms().len() as u64)
            .sum::<u64>();
        let l = u64::from(p).div_ceil(64).max(1);
        self.charge(
            n.saturating_pow(3)
                .saturating_add(terms.saturating_mul(n).saturating_mul(8))
                .saturating_mul(l)
                .saturating_mul(l)
                .saturating_mul(16),
        )
    }
}
fn limbs(r: &Rational) -> u64 {
    r.numerator()
        .bit_length()
        .div_ceil(64)
        .max(1)
        .saturating_add(r.denominator().bit_length().div_ceil(64).max(1))
}

enum EvalError {
    Limit(&'static str),
    Domain(Refusal),
    Pending(Refusal),
}
impl From<&'static str> for EvalError {
    fn from(code: &'static str) -> Self {
        Self::Limit(code)
    }
}
fn evaluate(x: &Real, p: u32, w: &mut Work, known_defined: bool) -> Result<Ball, EvalError> {
    let mut cache = HashMap::new();
    let value = evaluate_cached(x, p, w, &mut cache, known_defined)?;
    // No cached alias survives this evaluation; reclaim the root without copying it.
    drop(cache);
    Ok(Rc::unwrap_or_clone(value))
}
fn evaluate_cached(
    x: &Real,
    p: u32,
    w: &mut Work,
    cache: &mut HashMap<*const crate::real::Node, Rc<Ball>>,
    known_defined: bool,
) -> Result<Rc<Ball>, EvalError> {
    let key = Arc::as_ptr(&x.node);
    if let Some(b) = cache.get(&key) {
        return Ok(b.clone());
    }
    let b = Rc::new(evaluate_uncached(x, p, w, cache, known_defined)?);
    cache.insert(key, b.clone());
    Ok(b)
}
fn evaluate_uncached(
    x: &Real,
    p: u32,
    w: &mut Work,
    cache: &mut HashMap<*const crate::real::Node, Rc<Ball>>,
    known_defined: bool,
) -> Result<Ball, EvalError> {
    if x.node.fast_defined
        && let Some(r) = &x.node.rational
    {
        // The exact construction cache includes every domain dependency.
        // Do not reevaluate an entire defined rational subtree.
        w.ball(p, true)?;
        w.exact(r, r, true)?;
        return Ok(Ball::from_rational(r, p));
    }
    match &x.node.expr {
        Expr::Rational(r) => {
            w.ball(p, true)?;
            w.exact(r, r, true)?;
            Ok(Ball::from_rational(r, p))
        }
        Expr::Algebraic(a) => enclose_algebraic(a, p, w).map_err(Into::into),
        Expr::Coordinate(root, axis) => Ok(root.enclose_budgeted(p, w)?[*axis].clone()),
        Expr::Neg(a) => {
            let b = evaluate_cached(a, p, w, cache, known_defined)?;
            w.ball(p, false)?;
            Ok(-b.as_ref())
        }
        Expr::Square(a) => {
            let b = evaluate_cached(a, p, w, cache, known_defined)?;
            w.ball(p, true)?;
            Ok(b.square())
        }
        Expr::Pow(a, n) => {
            let b = evaluate_cached(a, p, w, cache, known_defined)?;
            for _ in 0..(64 - n.leading_zeros() * 2) {
                w.ball(p, true)?;
            }
            Ok(b.pow_u32(*n))
        }
        Expr::Sqrt(a) => {
            let mut b = evaluate_cached(a, p, w, cache, known_defined)?;
            if b.hi().sign() == Sign::Negative {
                return Err(EvalError::Domain(Refusal::new(
                    "numbers.sqrt-of-negative",
                    "negative square-root argument",
                )));
            }
            // Clipping alone loses half the working bits at a zero radicand,
            // and another half at each nested root. Try only bounded field
            // arithmetic here, never an annihilator or separation bound.
            if known_defined
                && b.contains_zero()
                && !b.rad().is_zero()
                && crate::tower::sign(a, w)? == Some(Sign::Zero)
            {
                b = Rc::new(Ball::from_rational(&Rational::zero(), p));
            }
            if !known_defined && b.lo().sign() == Sign::Negative {
                if zero_proof(a, &b, w)? {
                    b = Rc::new(Ball::from_rational(&Rational::zero(), p));
                } else {
                    return Err(EvalError::Pending(Refusal::new(
                        "numbers.sqrt-not-separated",
                        "square-root argument not proved nonnegative",
                    )));
                }
            }
            w.ball(p, true)?;
            Ok(b.sqrt().unwrap())
        }
        Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
            let lhs = evaluate_cached(a, p, w, cache, known_defined)?;
            let rhs = evaluate_cached(b, p, w, cache, known_defined)?;
            if matches!(&x.node.expr, Expr::Div(..)) && rhs.contains_zero() {
                if !known_defined && zero_proof(b, &rhs, w)? {
                    return Err(EvalError::Domain(Refusal::new(
                        "numbers.division-by-zero",
                        "exactly zero divisor",
                    )));
                }
                return Err(EvalError::Pending(Refusal::new(
                    "numbers.divisor-not-separated",
                    "divisor not separated from zero",
                )));
            }
            w.ball(p, matches!(&x.node.expr, Expr::Mul(..) | Expr::Div(..)))?;
            Ok(match &x.node.expr {
                Expr::Add(..) => lhs.as_ref() + rhs.as_ref(),
                Expr::Sub(..) => lhs.as_ref() - rhs.as_ref(),
                Expr::Mul(..) => lhs.as_ref() * rhs.as_ref(),
                Expr::Div(..) => lhs.checked_div(&rhs).unwrap(),
                _ => unreachable!(),
            })
        }
    }
}
const ROOT_SIGN_FILTER_BITS: u32 = 512;

fn root_polynomial_sign(
    polynomial: &UPoly,
    x: &Rational,
    precision: u32,
    work: &mut Work,
) -> Result<Sign, &'static str> {
    let denominator = x.denominator();
    let bits = denominator.bit_length();
    if bits <= u64::from(ROOT_SIGN_FILTER_BITS)
        || bits
            - denominator
                .0
                .trailing_zeros()
                .expect("positive denominator")
            <= u64::from(ROOT_SIGN_FILTER_BITS)
    {
        return Ok(work.eval_poly(polynomial, x)?.sign());
    }
    work.ball(precision, true)?;
    let argument = Ball::from_rational(x, precision);
    let value = algebraic::eval_ball_checked(polynomial, &argument, work)?;
    match value.sign() {
        Some(StrictSign::Negative) => Ok(Sign::Negative),
        Some(StrictSign::Positive) => Ok(Sign::Positive),
        None => Ok(work.eval_poly(polynomial, x)?.sign()),
    }
}

fn enclose_algebraic(a: &Algebraic, p: u32, w: &mut Work) -> Result<Ball, &'static str> {
    w.ball(p, true)?;
    match &a.repr {
        AlgebraicRepr::Rational(r) => Ok(Ball::from_rational(r, p)),
        AlgebraicRepr::Root(root) => {
            let filter_precision = p.min(ROOT_SIGN_FILTER_BITS);
            let start = w.used;
            let (mut lo, mut hi, mut lower_sign) = {
                let cache = root.enclosures.lock().unwrap();
                let allowance = w.policy.max_work.saturating_sub(w.used);
                if let Some(certificate) = cache.get(&p)
                    && certificate.cost <= allowance
                {
                    w.charge(certificate.cost)?;
                    return Ok(certificate.ball.clone());
                }
                if let Some((prefix_precision, prefix)) = cache.range(..p).next_back()
                    && (*prefix_precision).min(ROOT_SIGN_FILTER_BITS) == filter_precision
                    && prefix.cost <= allowance
                {
                    // Only identical sign-filter precision has identical
                    // prefix work. Replay it before continuing the same
                    // bisections; never substitute a later, tighter bracket.
                    w.charge(prefix.cost)?;
                    (prefix.lo.clone(), prefix.hi.clone(), prefix.lower_sign)
                } else {
                    (root.lo.clone(), root.hi.clone(), None)
                }
            };
            let target = Rational::one().mul_pow2(-i64::from(p));
            let ball = if lower_sign == Some(Sign::Zero) {
                Ball::from_rational(&lo, p)
            } else if &hi - &lo <= target {
                Ball::bounds(lo.clone(), hi.clone(), p)
            } else {
                let s = match lower_sign {
                    Some(s) => s,
                    None => root_polynomial_sign(&root.polynomial, &lo, filter_precision, w)?,
                };
                lower_sign = Some(s);
                loop {
                    w.exact(&lo, &hi, false)?;
                    let mid = (&lo + &hi).mul_pow2(-1);
                    let t = root_polynomial_sign(&root.polynomial, &mid, filter_precision, w)?;
                    if t == Sign::Zero {
                        lo = mid.clone();
                        hi = mid;
                        lower_sign = Some(Sign::Zero);
                        break Ball::from_rational(&lo, p);
                    }
                    if t == s {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                    if &hi - &lo <= target {
                        break Ball::bounds(lo.clone(), hi.clone(), p);
                    }
                }
            };
            let certificate = RootEnclosure {
                ball: ball.clone(),
                cost: w.used - start,
                lo,
                hi,
                lower_sign,
            };
            let mut cache = root.enclosures.lock().unwrap();
            if cache.len() == 8 && !cache.contains_key(&p) {
                cache.pop_first();
            }
            cache.insert(p, certificate);
            Ok(ball)
        }
    }
}
// Definedness is checked left-first in postorder, including under multiplication
// by zero and exponent zero. This fallback reads only immutable caches.
fn unresolved(x: &Real, limit: Option<Instant>) -> Option<Refusal> {
    #[cfg(any(test, feature = "test-deadline"))]
    test_deadline::checkpoint(test_deadline::Point::Definedness);
    if limit.is_some_and(deadline_passed) {
        return Some(deadline());
    }
    if x.node.fast_defined {
        return None;
    }
    match &x.node.expr {
        Expr::Rational(_) | Expr::Algebraic(_) | Expr::Coordinate(..) => None,
        Expr::Neg(a) | Expr::Square(a) | Expr::Pow(a, _) => unresolved(a, limit),
        Expr::Sqrt(a) => unresolved(a, limit).or_else(|| {
            if let Some(r) = &a.node.rational {
                return if r.sign() == Sign::Negative {
                    Some(Refusal::new(
                        "numbers.sqrt-of-negative",
                        "negative square-root argument",
                    ))
                } else {
                    None
                };
            }
            if a.approx().lo() < 0.0 {
                Some(Refusal::new(
                    "numbers.sqrt-not-separated",
                    "square-root argument not proved nonnegative",
                ))
            } else {
                None
            }
        }),
        Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) => {
            unresolved(a, limit).or_else(|| unresolved(b, limit))
        }
        Expr::Div(a, b) => unresolved(a, limit)
            .or_else(|| unresolved(b, limit))
            .or_else(|| {
                if b.node.rational.as_ref().is_some_and(Rational::is_zero) {
                    Some(Refusal::new(
                        "numbers.division-by-zero",
                        "exactly zero divisor",
                    ))
                } else if b.approx().sign().is_none() {
                    Some(Refusal::new(
                        "numbers.divisor-not-separated",
                        "divisor not separated from zero",
                    ))
                } else {
                    None
                }
            }),
    }
}
fn zero_proof(x: &Real, b: &Ball, w: &mut Work) -> Result<bool, &'static str> {
    if let Some(r) = &x.node.rational {
        w.exact(r, r, false)?;
        return Ok(r.is_zero());
    }
    if let Some(value) = crate::identity::rational(x, w)? {
        return Ok(value.is_zero());
    }
    if let Some(sign) = crate::tower::sign(x, w)? {
        return Ok(sign == Sign::Zero);
    }
    if crate::quadratic::zero(x, w)? {
        return Ok(true);
    }
    if crate::generator::rational(x, w)?.is_some_and(|r| r.is_zero()) {
        return Ok(true);
    }
    Ok(separation_bound(x, w)?.is_some_and(|bound| b.magnitude() < bound))
}
fn separation_bound(x: &Real, w: &mut Work) -> Result<Option<Rational>, &'static str> {
    if let Some(bound) = crate::height::separation(x, w)? {
        return Ok(Some(bound));
    }
    let Some(p) = annihilator(x, w)? else {
        return Ok(None);
    };
    let Some(k) = p.coefficients().iter().position(|c| !c.is_zero()) else {
        return Ok(None);
    };
    if k == 0 {
        return Ok(None);
    }
    let a = &p.coefficients()[k];
    let mut max = Rational::zero();
    for c in &p.coefficients()[k + 1..] {
        w.exact(c, a, true)?;
        max = max.max((c / a).abs());
    }
    // Cauchy's polynomial root bound, applied to the reciprocal polynomial.
    // For 0 < |x| < 1/(1+M), M=max_i>k |a_i/a_k|, the nonconstant
    // terms of P(x)/x^k have relative sum <= M|x|/(1-|x|) < 1.
    // Thus every nonzero root has |x| >= 1/(1+M). This proof uses
    // the actual annihilating polynomial, including system-coordinate leaves.
    let bound = (&Rational::one() + &max).recip().unwrap();
    Ok(Some(bound))
}
fn annihilator(x: &Real, w: &mut Work) -> Result<Option<UPoly>, &'static str> {
    annihilator_cached(x, w, &mut HashMap::new())
}
fn annihilator_cached(
    x: &Real,
    w: &mut Work,
    cache: &mut HashMap<*const crate::real::Node, Option<UPoly>>,
) -> Result<Option<UPoly>, &'static str> {
    let key = Arc::as_ptr(&x.node);
    if let Some(p) = cache.get(&key) {
        return Ok(p.clone());
    }
    let p = annihilator_uncached(x, w, cache)?;
    cache.insert(key, p.clone());
    Ok(p)
}
fn annihilator_uncached(
    x: &Real,
    w: &mut Work,
    cache: &mut HashMap<*const crate::real::Node, Option<UPoly>>,
) -> Result<Option<UPoly>, &'static str> {
    if let Some(r) = &x.node.rational {
        w.exact(r, r, false)?;
        return Ok(Some(UPoly::new(vec![-r, Rational::one()])));
    }
    if !matches!(&x.node.expr, Expr::Algebraic(_) | Expr::Coordinate(..))
        && let Some(p) = crate::generator::annihilator(x, w)?
    {
        return Ok(Some(p));
    }
    let p = match &x.node.expr {
        Expr::Rational(r) => UPoly::new(vec![-r, Rational::one()]),
        Expr::Algebraic(a) => a.defining_polynomial(),
        Expr::Coordinate(root, axis) => {
            return system::coordinate_polynomial(root.system(), *axis, w);
        }
        Expr::Neg(a) => {
            let Some(p) = annihilator_cached(a, w, cache)? else {
                return Ok(None);
            };
            UPoly::new(
                p.coefficients()
                    .iter()
                    .enumerate()
                    .map(|(i, c)| if i % 2 == 0 { c.clone() } else { -c })
                    .collect(),
            )
        }
        Expr::Sqrt(a) => {
            let Some(p) = annihilator_cached(a, w, cache)? else {
                return Ok(None);
            };
            let mut c = vec![Rational::zero(); 2 * p.coefficients().len()];
            for (i, a) in p.coefficients().iter().enumerate() {
                c[2 * i] = a.clone();
            }
            UPoly::new(c)
        }
        Expr::Square(a) => {
            let Some(p) = annihilator_cached(a, w, cache)? else {
                return Ok(None);
            };
            square_annihilator(&p, w)?
        }
        Expr::Pow(a, n) => {
            let Some(mut p) = annihilator_cached(a, w, cache)? else {
                return Ok(None);
            };
            let mut result = UPoly::new(vec![Rational::from_i64(-1), Rational::one()]);
            let mut n = *n;
            while n != 0 {
                if n & 1 != 0 {
                    result = combine(&result, &p, 2, w)?;
                }
                n >>= 1;
                if n != 0 {
                    p = square_annihilator(&p, w)?;
                }
            }
            result
        }
        Expr::Mul(a, b) if Arc::ptr_eq(&a.node, &b.node) => {
            let Some(p) = annihilator_cached(a, w, cache)? else {
                return Ok(None);
            };
            square_annihilator(&p, w)?
        }
        Expr::Sub(a, b) if Arc::ptr_eq(&a.node, &b.node) => UPoly::x(),
        Expr::Div(a, b) if Arc::ptr_eq(&a.node, &b.node) => {
            // Definedness has already proved the common operand nonzero.
            UPoly::new(vec![-Rational::one(), Rational::one()])
        }
        Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
            let Some(a) = annihilator_cached(a, w, cache)? else {
                return Ok(None);
            };
            let Some(b) = annihilator_cached(b, w, cache)? else {
                return Ok(None);
            };
            let op = match &x.node.expr {
                Expr::Add(..) => 0,
                Expr::Sub(..) => 1,
                Expr::Mul(..) => 2,
                _ => 3,
            };
            combine(&a, &b, op, w)?
        }
    };
    if p.is_zero() {
        return Ok(None);
    }
    w.poly(&p, &p)?;
    let gcd = p.gcd_checked(&p.derivative(), &|| w.check_deadline())?;
    let p = p.div_rem_checked(&gcd, &|| w.check_deadline())?.unwrap().0;
    Ok(Some(UPoly::new(
        p.primitive_integer_checked(&|| w.check_deadline())?
            .into_iter()
            .map(Rational::from_integer)
            .collect(),
    )))
}
fn square_annihilator(p: &UPoly, w: &mut Work) -> Result<UPoly, &'static str> {
    // P(t) = E(t²) + t O(t²). Therefore P(t)=0 implies
    // E(t²)² - t² O(t²)² = 0. Unlike treating two copies of t as
    // independent conjugates in a resultant, this keeps degree <= deg(P).
    let even = UPoly::new(p.coefficients().iter().step_by(2).cloned().collect());
    let odd = UPoly::new(
        p.coefficients()
            .iter()
            .skip(1)
            .step_by(2)
            .cloned()
            .collect(),
    );
    w.poly(&even, &even)?;
    w.poly(&odd, &odd)?;
    let even_squared = &even * &even;
    let odd_squared = &odd * &odd;
    w.poly(&UPoly::x(), &odd_squared)?;
    let shifted = &UPoly::x() * &odd_squared;
    w.poly(&even_squared, &shifted)?;
    Ok(&even_squared - &shifted)
}
fn combine(a: &UPoly, b: &UPoly, op: u8, w: &mut Work) -> Result<UPoly, &'static str> {
    // A degree-one annihilator fixes that operand to one rational value.
    // Substitute the inverse affine map rather than constructing a Sylvester
    // matrix. The resulting polynomial still vanishes at the exact result.
    if b.degree() == Some(1) {
        w.exact(&b.coefficients()[0], &b.coefficients()[1], true)?;
        let d = -(&b.coefficients()[0] / &b.coefficients()[1]);
        w.exact(&d, &Rational::one(), true)?;
        let inner = match op {
            0 => Some(UPoly::new(vec![-&d, Rational::one()])),
            1 => Some(UPoly::new(vec![d.clone(), Rational::one()])),
            2 if d.is_zero() => return Ok(UPoly::x()),
            2 => Some(UPoly::new(vec![Rational::zero(), d.recip().unwrap()])),
            3 if !d.is_zero() => Some(UPoly::new(vec![Rational::zero(), d])),
            _ => None,
        };
        if let Some(inner) = inner {
            return compose_budgeted(a, &inner, w);
        }
    }
    if a.degree() == Some(1) && op != 3 {
        w.exact(&a.coefficients()[0], &a.coefficients()[1], true)?;
        let c = -(&a.coefficients()[0] / &a.coefficients()[1]);
        w.exact(&c, &Rational::one(), true)?;
        let inner = match op {
            0 => UPoly::new(vec![-&c, Rational::one()]),
            1 => UPoly::new(vec![c, -Rational::one()]),
            2 if c.is_zero() => return Ok(UPoly::x()),
            2 => UPoly::new(vec![Rational::zero(), c.recip().unwrap()]),
            _ => unreachable!(),
        };
        return compose_budgeted(b, &inner, w);
    }
    let y = MPoly::var(2, 1);
    let x = MPoly::var(2, 0);
    let lift = |p: &UPoly| {
        MPoly::from_terms(
            2,
            p.coefficients()
                .iter()
                .enumerate()
                .map(|(i, c)| (vec![0, i as u32], c.clone()))
                .collect(),
        )
    };
    let (f, g) = match op {
        0 => (lift(a), lift(b).substitute(1, &(&x - &y))),
        1 => (lift(a), lift(b).substitute(1, &(&y - &x))),
        2 => {
            let n = b.degree().unwrap();
            (
                lift(a),
                MPoly::from_terms(
                    2,
                    b.coefficients()
                        .iter()
                        .enumerate()
                        .map(|(i, c)| (vec![i as u32, (n - i) as u32], c.clone()))
                        .collect(),
                ),
            )
        }
        _ => {
            let k = b.coefficients().iter().position(|c| !c.is_zero()).unwrap();
            (
                MPoly::from_terms(
                    2,
                    a.coefficients()
                        .iter()
                        .enumerate()
                        .map(|(i, c)| (vec![i as u32, i as u32], c.clone()))
                        .collect(),
                ),
                lift(&UPoly::new(b.coefficients()[k..].to_vec())),
            )
        }
    };
    w.mpoly(&f, &g)?;
    Ok(f.resultant_checked(&g, 1, &|| w.check_deadline())?
        .to_upoly(0)
        .unwrap())
}
fn compose_budgeted(p: &UPoly, inner: &UPoly, w: &mut Work) -> Result<UPoly, &'static str> {
    let mut result = UPoly::zero();
    for c in p.coefficients().iter().rev() {
        w.poly(&result, inner)?;
        result = &result * inner;
        let constant = UPoly::constant(c.clone());
        w.poly(&result, &constant)?;
        result = &result + &constant;
    }
    Ok(result)
}
