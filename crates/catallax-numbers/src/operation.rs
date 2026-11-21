//! Request-scoped deterministic accounting, independent of wall time.

use crate::{
    Ball, BudgetPolicy, Real, UPoly,
    algebraic::{AlgebraicRoot, RootEnclosure},
    decider::CachedRoots,
};
use std::{
    cell::{Cell, RefCell},
    cmp::Ordering,
    collections::{BTreeMap, HashMap},
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, AtomicU64, Ordering as AtomicOrdering},
    },
};

static NEXT_REQUEST: AtomicU64 = AtomicU64::new(1);

#[derive(Clone)]
struct State {
    request: u64,
    aborted: Arc<AtomicBool>,
    limit: u64,
    used: u64,
    // Charges that `refund_on_refusal` took back out of `used`.
    refunded: u64,
    // Arithmetic counted at the decider's rate and never charged: what was
    // added by hand (`count_uncharged`) and what the number types counted
    // while an earlier holder of this thread had the state set aside. What
    // they have counted since it was installed is in `ARITHMETIC`. It wraps,
    // because a reservation takes its price off before the operation it pays
    // for adds it (`reserved_arithmetic`); `work` reports it floored at zero.
    uncharged: u64,
    exhausted: bool,
    reduced_angle_series: bool,
    reuse_enclosures: bool,
    // A worker running a loop item ahead of the writer (`ahead.rs`). Its
    // counter is its own and is never reported.
    worker: bool,
    ticks: u32,
}

#[derive(Clone)]
pub(crate) struct RootProofKey {
    root: Arc<AlgebraicRoot>,
    precision: u32,
}
impl RootProofKey {
    pub(crate) fn new(root: &Arc<AlgebraicRoot>, precision: u32) -> Self {
        Self {
            root: root.clone(),
            precision,
        }
    }
}
impl Ord for RootProofKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.precision
            .cmp(&other.precision)
            .then_with(|| {
                self.root
                    .polynomial
                    .coefficients()
                    .cmp(other.root.polynomial.coefficients())
            })
            .then_with(|| self.root.lo.cmp(&other.root.lo))
            .then_with(|| self.root.hi.cmp(&other.root.hi))
    }
}
impl PartialOrd for RootProofKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl PartialEq for RootProofKey {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}
impl Eq for RootProofKey {}

// Persistent maps: a loop handed to workers takes a copy for each of them
// without copying the entries.
type RootProofs = im::OrdMap<RootProofKey, Arc<RootEnclosure>>;
type RootListProofs = im::OrdMap<(u32, u64, bool), im::HashMap<Arc<UPoly>, Arc<CachedRoots>>>;
type EnclosureKey = (*const crate::real::Node, u32, u32, u64, bool);
type EnclosureProofs = HashMap<EnclosureKey, (Arc<crate::real::Node>, EnclosureProof)>;

#[derive(Clone)]
pub(crate) struct EnclosureProof {
    pub ball: Rc<Ball>,
    pub work: u64,
    pub bits: u32,
}

thread_local! {
    static ACTIVE: RefCell<Option<State>> = const { RefCell::new(None) };
    static RESERVATIONS: RefCell<Option<Vec<u64>>> = const { RefCell::new(None) };
    static ROOT_PROOFS: RefCell<RootProofs> = RefCell::new(RootProofs::new());
    static ROOT_LIST_PROOFS: RefCell<RootListProofs> = RefCell::new(RootListProofs::new());
    static ENCLOSURE_PROOFS: RefCell<EnclosureProofs> = RefCell::new(HashMap::new());
    static PI_PROOFS: RefCell<BTreeMap<u32, Ball>> = const { RefCell::new(BTreeMap::new()) };
    static ARITHMETIC: Cell<Arithmetic> = const { Cell::new(Arithmetic { closed: 0, units: 0 }) };
}

// The count the number types keep on this thread (`count_arithmetic`).
#[derive(Clone, Copy)]
struct Arithmetic {
    // Above zero while operations are not counted: inside a decision, inside
    // an operation counted whole, or inside the fill of a shared memo.
    closed: u32,
    // Counted since this thread's state was installed. It wraps, as
    // `State::uncharged` does.
    units: u64,
}

/// Accounting for one complete operation, including all its deciders.
#[derive(Clone, Copy, Debug)]
pub struct OperationWork {
    /// Accepted arithmetic reservations and cancellation checkpoints.
    /// An aborted task conservatively consumes the request's remaining cap.
    pub used: u64,
    /// Whether a reservation exceeded the allowance or a task aborted.
    pub exhausted: bool,
    /// Arithmetic of refused optional attempts ([`refund_on_refusal`]). The
    /// allowance got these units back, so they are not in `used`, but the
    /// machine did the work. Every refused attempt inside the operation is
    /// counted once, including one made by an [`OperationTask`] item, whether
    /// that item is then admitted or discarded. A discarded item's other
    /// charges were never admitted and are in neither counter.
    pub refunded: u64,
    /// Ball and rational arithmetic that ran outside a decision, priced at the
    /// decider's rate and not charged: what `used` would have risen by had it
    /// been charged. It is outside `used` and outside the allowance and never
    /// exhausts a request. `used + refunded + uncharged` is the request's
    /// arithmetic counted once.
    ///
    /// The number types count it themselves as they run, so no caller can
    /// leave an operation out: every sum, difference, product, square,
    /// reciprocal and square root of a [`Ball`], and every sum, difference,
    /// product, quotient and decimal conversion of a [`crate::Rational`]. The
    /// rules are at `count_arithmetic` in `operation.rs`.
    ///
    /// The count is that of one thread. Arithmetic a helper thread did reaches
    /// it only with a record the request adopts ([`crate::Charges::adopt`]),
    /// as that record's charges reach `used`. An [`OperationTask`] item
    /// reports its count as it reports `refunded`.
    pub uncharged: u64,
}

impl State {
    fn work(&self) -> OperationWork {
        let aborted = self.aborted.load(AtomicOrdering::Acquire);
        OperationWork {
            used: if aborted { self.limit } else { self.used },
            exhausted: self.exhausted || aborted,
            refunded: self.refunded,
            uncharged: floor(self.counted()),
        }
    }

    // The wrapping count of the installed state: what it holds and what the
    // number types have counted since it was installed.
    fn counted(&self) -> u64 {
        self.uncharged
            .wrapping_add(ARITHMETIC.with(|count| count.get().units))
    }
}

// A wrapped count below zero reads as zero: reservations whose operations
// never ran took off more than was counted.
fn floor(counted: u64) -> u64 {
    if counted > u64::MAX / 2 { 0 } else { counted }
}

// Makes `state` this thread's request state. The number types' count starts
// again from zero for it.
fn install(state: State) {
    ARITHMETIC.with(|count| {
        count.set(Arithmetic {
            units: 0,
            ..count.get()
        });
    });
    ACTIVE.with(|slot| *slot.borrow_mut() = Some(state));
}

// Sets this thread's request state aside, with what the number types counted
// for it folded in.
fn uninstall() -> Option<State> {
    let mut state = ACTIVE.with(|slot| slot.take())?;
    state.uncharged = state.counted();
    ARITHMETIC.with(|count| {
        count.set(Arithmetic {
            units: 0,
            ..count.get()
        });
    });
    Some(state)
}

/// Whether the caller already owns a request-scoped work allowance.
/// Consumers must not replace it with a fresh allowance for a nested build.
pub fn operation_budget_active() -> bool {
    ACTIVE.with(|slot| slot.borrow().is_some())
}

/// Reuse numeric enclosure proofs for the rest of this operation.
/// Each operation starts with reuse disabled, including reserved allowances.
pub fn reuse_operation_enclosures() {
    ACTIVE.with(|slot| {
        if let Some(state) = slot.borrow_mut().as_mut() {
            state.reuse_enclosures = true;
        }
    });
}

pub(crate) fn enclosure_reuse_active() -> bool {
    ACTIVE.with(|slot| {
        slot.borrow()
            .as_ref()
            .is_some_and(|state| state.reuse_enclosures)
    })
}

/// Prefer repeated half-angle reductions for the rest of this operation.
/// Enclosing balls retain their working precision and outward guarantee.
/// This preference is inactive outside an operation budget.
pub fn prefer_reduced_angle_series() {
    ACTIVE.with(|slot| {
        if let Some(state) = slot.borrow_mut().as_mut() {
            state.reduced_angle_series = true;
        }
    });
}

pub(crate) fn angle_reductions() -> u32 {
    ACTIVE.with(|slot| {
        if slot
            .borrow()
            .as_ref()
            .is_some_and(|state| state.reduced_angle_series)
        {
            4
        } else {
            1
        }
    })
}

/// Observes request accounting without charging it.
/// Outside a request budget every counter is zero/false.
pub fn current_operation_work() -> OperationWork {
    ACTIVE.with(|slot| {
        slot.borrow().as_ref().map_or_else(
            || OperationWork {
                used: 0,
                exhausted: false,
                refunded: 0,
                uncharged: 0,
            },
            State::work,
        )
    })
}

/// Runs a synchronous operation with a deterministic work allowance.
///
/// All deciders on this thread share the allowance. Other threads use
/// [`OperationTask`]'s ordered admission or `RunAhead`/`Helpers`' replay to charge
/// the coordinator. Nested scopes are rejected. Panic unwinding
/// restores the thread, so a later request cannot inherit another request's work.
/// Callers must discard the returned value when `OperationWork::exhausted` is
/// true: lower layers may have handled a refusal and continued collecting data.
pub fn with_operation_budget<T>(limit: u64, run: impl FnOnce() -> T) -> (T, OperationWork) {
    budget_scope(limit, false, run)
}

/// Runs an optional attempt whose refusal must not spend the caller's allowance.
///
/// A success keeps every charge, exactly as if the attempt had run inline, so
/// a request that the attempt completes costs what it always did. A refusal
/// returns the attempt's charges and clears an exhaustion the attempt latched:
/// the construction that follows starts with the allowance the attempt found.
/// Request-local proof caches stay warm. An exhaustion that predates the
/// attempt, and an aborted task, stay exhausted. Outside a request allowance
/// this only runs the attempt.
///
/// Use it for a cheaper construction tried ahead of the one the request is
/// funded for. The arithmetic of a refused attempt is real machine time. Its
/// returned charges move to [`OperationWork::refunded`], which only grows, so
/// the cost of refused attempts stays on the record.
pub fn refund_on_refusal<T, E>(run: impl FnOnce() -> Result<T, E>) -> Result<T, E> {
    let before = ACTIVE.with(|slot| {
        slot.borrow()
            .as_ref()
            .map(|state| (state.request, state.used, state.exhausted))
    });
    let recorded = RESERVATIONS.with(|slot| slot.borrow().as_ref().map(Vec::len));
    // A tape of one decision cannot describe charges that are later returned.
    crate::ahead::unrecordable();
    let result = run();
    if result.is_err()
        && let Some((request, used, exhausted)) = before
    {
        ACTIVE.with(|slot| {
            if let Some(state) = slot.borrow_mut().as_mut()
                && state.request == request
                && !state.aborted.load(AtomicOrdering::Acquire)
            {
                state.refunded = state
                    .refunded
                    .saturating_add(state.used.saturating_sub(used));
                state.used = used;
                state.exhausted = exhausted;
            }
        });
        RESERVATIONS.with(|slot| {
            if let (Some(reservations), Some(length)) = (slot.borrow_mut().as_mut(), recorded) {
                reservations.truncate(length);
            }
        });
    }
    result
}

// `worker` says whose scope this is: a thread running ahead of a writer stays
// one inside an allowance it runs, so it still leaves the shared caches alone.
fn budget_scope<T>(limit: u64, worker: bool, run: impl FnOnce() -> T) -> (T, OperationWork) {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            ACTIVE.with(|slot| *slot.borrow_mut() = None);
            ROOT_PROOFS.with(|proofs| proofs.borrow_mut().clear());
            ROOT_LIST_PROOFS.with(|proofs| proofs.borrow_mut().clear());
            ENCLOSURE_PROOFS.with(|proofs| proofs.borrow_mut().clear());
            PI_PROOFS.with(|proofs| proofs.borrow_mut().clear());
        }
    }
    assert!(
        ACTIVE.with(|slot| slot.borrow().is_none()),
        "operation budgets cannot nest"
    );
    install(State {
        request: NEXT_REQUEST.fetch_add(1, AtomicOrdering::Relaxed),
        aborted: Arc::new(AtomicBool::new(false)),
        limit,
        used: 0,
        refunded: 0,
        uncharged: 0,
        exhausted: false,
        reduced_angle_series: false,
        reuse_enclosures: false,
        worker,
        ticks: 0,
    });
    ROOT_PROOFS.with(|proofs| proofs.borrow_mut().clear());
    ROOT_LIST_PROOFS.with(|proofs| proofs.borrow_mut().clear());
    ENCLOSURE_PROOFS.with(|proofs| proofs.borrow_mut().clear());
    PI_PROOFS.with(|proofs| proofs.borrow_mut().clear());
    let _reset = Reset;
    let result = run();
    (result, current_operation_work())
}

/// Single-use snapshot for numeric jobs with cold per-item proof caches.
/// Capture one wave, split it between inputs, then commit all computed work
/// before capturing another wave or retry. Workers never duplicate its cap.
pub struct OperationTask {
    state: Option<State>,
}

/// A speculative value whose reservations have not yet been admitted.
/// Admitted or discarded, it reports its refused attempts and its uncharged
/// arithmetic to its request.
pub struct RecordedWork<T> {
    value: T,
    reservations: Vec<u64>,
    request: Option<u64>,
    work: Option<OperationWork>,
    // Never read: it reports when the record is dropped.
    _unreserved: Unreserved,
}

// What one recorded item did outside the allowance: its refused attempts
// (`refund_on_refusal`) and its uncharged arithmetic (`count_uncharged`).
// Neither left a reservation to replay, and the arithmetic was done whether the
// item is then admitted, refused at admission or discarded. So both counts go
// to the request when the record leaves, by any of those ways, and so exactly
// once. A record that leaves outside its request has no counter to report to.
struct Unreserved {
    request: Option<u64>,
    refunded: u64,
    // The item's wrapping count (`uncharged`), not its floored report.
    uncharged: u64,
}

impl Drop for Unreserved {
    fn drop(&mut self) {
        if self.refunded == 0 && self.uncharged == 0 {
            return;
        }
        // A record can be dropped while its thread unwinds or ends.
        let _ = ACTIVE.try_with(|slot| {
            if let Ok(mut active) = slot.try_borrow_mut()
                && let Some(state) = active.as_mut()
                && Some(state.request) == self.request
            {
                state.refunded = state.refunded.saturating_add(self.refunded);
                state.uncharged = state.uncharged.wrapping_add(self.uncharged);
            }
        });
    }
}

impl OperationTask {
    /// Captures the remaining request cap and arithmetic preferences.
    /// The coordinator must finish and account for this wave before another
    /// capture. Input shares, not the number of workers, determine allocation.
    pub fn capture() -> Self {
        Self {
            state: ACTIVE.with(|slot| {
                slot.borrow().as_ref().map(|state| {
                    let work = state.work();
                    State {
                        limit: state.limit.saturating_sub(work.used),
                        used: 0,
                        refunded: 0,
                        uncharged: 0,
                        exhausted: work.exhausted,
                        ..state.clone()
                    }
                })
            }),
        }
    }

    /// Divides the captured cap between fixed input items. The first remainder
    /// items receive one extra unit. Shares sum to exactly the captured cap.
    pub fn split(self, count: std::num::NonZeroUsize) -> Vec<Self> {
        let count = count.get();
        (0..count)
            .map(|index| Self {
                state: self.state.as_ref().map(|state| State {
                    limit: state.limit / count as u64
                        + u64::from((index as u64) < state.limit % count as u64),
                    ..state.clone()
                }),
            })
            .collect()
    }

    /// Computes one independent item. Each item starts with empty proof caches,
    /// even when it executes on the caller or a reused pool thread.
    /// A panic exhausts the originating request, even if the caller catches it.
    pub fn run<T>(self, run: impl FnOnce() -> T) -> RecordedWork<T> {
        let Some(state) = self.state else {
            return RecordedWork {
                value: run(),
                reservations: Vec::new(),
                request: None,
                work: None,
                _unreserved: Unreserved {
                    request: None,
                    refunded: 0,
                    uncharged: 0,
                },
            };
        };
        struct AbortOnPanic(Arc<AtomicBool>);
        impl Drop for AbortOnPanic {
            fn drop(&mut self) {
                if std::thread::panicking() {
                    self.0.store(true, AtomicOrdering::Release);
                }
            }
        }
        // No record can be committed after unwind. Spend the originating
        // request's remaining cap, including when this task runs on a worker.
        let _abort = AbortOnPanic(Arc::clone(&state.aborted));
        let request = state.request;
        let state = State {
            worker: state.worker || worker(),
            ticks: 0,
            ..state
        };
        struct Restore(Option<Vec<u64>>);
        impl Drop for Restore {
            fn drop(&mut self) {
                RESERVATIONS.with(|slot| *slot.borrow_mut() = self.0.take());
            }
        }
        let ((value, reservations, uncharged), work) = OperationAllowance { limit: state.limit }
            .run(|| {
                install(state);
                let _restore = Restore(RESERVATIONS.with(|slot| slot.replace(Some(Vec::new()))));
                let value = run();
                let reservations = RESERVATIONS.with(|slot| slot.take().unwrap());
                (value, reservations, uncharged())
            });
        RecordedWork {
            value,
            reservations,
            request: Some(request),
            work: Some(work),
            _unreserved: Unreserved {
                request: Some(request),
                refunded: work.refunded,
                uncharged,
            },
        }
    }
}

impl<T> RecordedWork<T> {
    fn check_request(&self) -> Result<(), crate::Refusal> {
        if ACTIVE.with(|slot| slot.borrow().as_ref().map(|state| state.request)) != self.request {
            return Err(crate::Refusal::new(
                "numbers.operation-context",
                "numeric batch result belongs to a different request",
            ));
        }
        Ok(())
    }

    /// Charges the actual accepted reservations of one input share.
    /// `None` means the item exhausted its share and its value is discarded.
    /// Commit every share before retrying incomplete items in input order with
    /// the remaining request cap. Failed attempts are paid, never refunded.
    pub fn commit_share(self) -> Result<Option<T>, crate::Refusal> {
        self.check_request()?;
        let mut accepted = self.work.map_or(0, |work| work.used);
        for units in self.reservations {
            if units > accepted {
                break;
            }
            accepted -= units;
            if charge(units) {
                return Err(crate::Refusal::new(
                    "numbers.operation-budget",
                    "numeric batch item exceeds the request allowance",
                ));
            }
        }
        if exhausted() {
            return Err(crate::Refusal::new(
                "numbers.operation-budget",
                "numeric batch item exceeds the request allowance",
            ));
        }
        Ok(if self.work.is_some_and(|work| work.exhausted) {
            None
        } else {
            Some(self.value)
        })
    }

    /// Admits this item on the original request thread, in input order.
    /// Reservations replay individually, preserving the exact exhausted prefix.
    /// Later speculative items must be discarded after the first refusal.
    /// A discarded item's charges are not admitted; its refused attempts are
    /// still reported in [`OperationWork::refunded`] and its uncharged
    /// arithmetic in [`OperationWork::uncharged`].
    pub fn commit(self) -> Result<T, crate::Refusal> {
        self.check_request()?;
        for units in self.reservations {
            if charge(units) {
                return Err(crate::Refusal::new(
                    "numbers.operation-budget",
                    "numeric batch item exceeds the request allowance",
                ));
            }
        }
        if exhausted() {
            return Err(crate::Refusal::new(
                "numbers.operation-budget",
                "numeric batch item exceeds the request allowance",
            ));
        }
        Ok(self.value)
    }
}

/// A single-use allowance reserved for work that may run after its caller returns.
/// Not cloneable: splitting or repeating a request cannot duplicate its work.
#[derive(Debug)]
pub struct OperationAllowance {
    limit: u64,
}

impl OperationAllowance {
    /// Reserves all units from the active request, or owns `limit` units when
    /// there is no enclosing request. Reservations are charged immediately and
    /// never refunded, including when the job is cancelled.
    pub fn reserve(limit: u64) -> Result<Self, crate::Refusal> {
        if charge(limit) {
            return Err(crate::Refusal::new(
                "numbers.operation-budget",
                "deferred work allocation exceeds the caller allowance",
            ));
        }
        Ok(Self { limit })
    }

    /// Moves units to a separate single-use job without creating work.
    pub fn split_off(&mut self, limit: u64) -> Option<Self> {
        self.limit = self.limit.checked_sub(limit)?;
        Some(Self { limit })
    }

    /// The reserved cap, not the eventual amount of arithmetic performed.
    pub fn limit(&self) -> u64 {
        self.limit
    }

    /// Consumes the reservation. An enclosing request has already paid for
    /// these units; the returned counters report actual work separately.
    /// Request-local proof caches are not carried across scheduling boundaries.
    pub fn run<T>(self, run: impl FnOnce() -> T) -> (T, OperationWork) {
        struct Restore {
            state: Option<State>,
            reservations: Option<Vec<u64>>,
            roots: RootProofs,
            lists: RootListProofs,
            enclosures: EnclosureProofs,
            pi: BTreeMap<u32, Ball>,
        }
        impl Drop for Restore {
            fn drop(&mut self) {
                ACTIVE.with(|slot| *slot.borrow_mut() = None);
                if let Some(state) = self.state.take() {
                    install(state);
                }
                RESERVATIONS.with(|slot| *slot.borrow_mut() = self.reservations.take());
                ROOT_PROOFS.with(|slot| *slot.borrow_mut() = std::mem::take(&mut self.roots));
                ROOT_LIST_PROOFS.with(|slot| *slot.borrow_mut() = std::mem::take(&mut self.lists));
                ENCLOSURE_PROOFS
                    .with(|slot| *slot.borrow_mut() = std::mem::take(&mut self.enclosures));
                PI_PROOFS.with(|slot| *slot.borrow_mut() = std::mem::take(&mut self.pi));
            }
        }
        // A record cannot describe a counter that starts again from zero or
        // proof caches set aside part-way, so a thread running ahead of a
        // writer records nothing across this scope. Declared first: it ends
        // after the outer counter is back.
        let worker = worker();
        let _unrecorded = crate::ahead::suspend(worker);
        let _restore = Restore {
            state: uninstall(),
            reservations: RESERVATIONS.with(|slot| slot.take()),
            roots: ROOT_PROOFS.with(|slot| std::mem::take(&mut *slot.borrow_mut())),
            lists: ROOT_LIST_PROOFS.with(|slot| std::mem::take(&mut *slot.borrow_mut())),
            enclosures: ENCLOSURE_PROOFS.with(|slot| std::mem::take(&mut *slot.borrow_mut())),
            pi: PI_PROOFS.with(|slot| std::mem::take(&mut *slot.borrow_mut())),
        };
        budget_scope(self.limit, worker, run)
    }
}

// A cached node retains its own cold arithmetic cost. Its children are visited
// separately so a decision still counts each node of its DAG exactly once.
pub(crate) fn enclosure_proof(
    value: &Real,
    precision: u32,
    policy: BudgetPolicy,
    known_defined: bool,
) -> Option<EnclosureProof> {
    if !enclosure_reuse_active() {
        return None;
    }
    crate::ahead::unrecordable();
    let found = ENCLOSURE_PROOFS.with(|proofs| {
        proofs
            .borrow()
            .get(&enclosure_key(value, precision, policy, known_defined))
            .map(|(_, proof)| proof.clone())
    });
    #[cfg(feature = "work-span")]
    if let Some(proof) = &found {
        crate::work_span::read("proof cache: enclosure", Rc::as_ptr(&proof.ball) as u64);
    }
    found
}

pub(crate) fn record_enclosure_proof(
    value: &Real,
    precision: u32,
    policy: BudgetPolicy,
    known_defined: bool,
    proof: EnclosureProof,
) {
    if !enclosure_reuse_active() {
        return;
    }
    crate::ahead::unrecordable();
    #[cfg(feature = "work-span")]
    crate::work_span::write("proof cache: enclosure", Rc::as_ptr(&proof.ball) as u64);
    ENCLOSURE_PROOFS.with(|proofs| {
        let mut proofs = proofs.borrow_mut();
        if proofs.len() >= 8192 {
            proofs.clear();
        }
        // Pin the node so its address cannot name another definition.
        proofs.insert(
            enclosure_key(value, precision, policy, known_defined),
            (value.node.clone(), proof),
        );
    });
}

fn enclosure_key(
    value: &Real,
    precision: u32,
    policy: BudgetPolicy,
    known_defined: bool,
) -> EnclosureKey {
    (
        Arc::as_ptr(&value.node),
        precision,
        policy.max_precision_bits,
        policy.max_work,
        known_defined,
    )
}

// Unlike expression enclosures, this constant has no definition or policy key:
// its completed ball depends only on precision. Reuse it in every request,
// with the same bounded cache and request/deferred-allowance lifetime as before.
pub(crate) fn pi_proof(precision: u32) -> Option<Ball> {
    if !operation_budget_active() {
        return None;
    }
    crate::ahead::unrecordable();
    let found = PI_PROOFS.with(|proofs| proofs.borrow().get(&precision).cloned());
    #[cfg(feature = "work-span")]
    if found.is_some() {
        crate::work_span::read("proof cache: pi", u64::from(precision));
    }
    found
}

pub(crate) fn record_pi_proof(precision: u32, ball: &Ball) {
    if operation_budget_active() {
        PI_PROOFS.with(|proofs| {
            let mut proofs = proofs.borrow_mut();
            if proofs.len() >= 32 {
                proofs.pop_first();
            }
            proofs.insert(precision, ball.clone());
        });
        #[cfg(feature = "work-span")]
        crate::work_span::write("proof cache: pi", u64::from(precision));
    }
}

pub(crate) fn charge(units: u64) -> bool {
    #[cfg(feature = "work-span")]
    let mut advanced = false;
    RESERVATIONS.with(|slot| {
        if let Some(reservations) = &mut *slot.borrow_mut() {
            reservations.push(units);
        }
    });
    let exhausted = ACTIVE.with(|slot| {
        let mut active = slot.borrow_mut();
        let Some(state) = active.as_mut() else {
            return false;
        };
        if state.aborted.load(AtomicOrdering::Acquire) {
            state.used = state.limit;
            state.exhausted = true;
        } else if units > state.limit.saturating_sub(state.used) {
            state.exhausted = true;
        } else if !state.exhausted {
            state.used += units;
            #[cfg(feature = "work-span")]
            {
                advanced = true;
            }
        }
        if state.worker {
            // An item the writer no longer wants stops at its next charges.
            state.ticks = state.ticks.wrapping_add(1);
            if state.ticks % 64 == 0 && crate::ahead::worker_tick() {
                state.exhausted = true;
            }
        }
        state.exhausted
    });
    #[cfg(feature = "work-span")]
    if advanced {
        // The observer reads worker state from ACTIVE, so release its borrow first.
        crate::work_span::advance(units);
    }
    exhausted
}

/// Adds `units` to [`OperationWork::uncharged`]: a count taken on another
/// thread, or for another item, that this request now owns. Unlike a charge
/// it takes nothing from the allowance, latches nothing and cannot refuse;
/// outside a request it does nothing. `units` is a wrapping count, as
/// [`uncharged`] returns.
pub(crate) fn count_uncharged(units: u64) {
    ACTIVE.with(|slot| {
        if let Some(state) = slot.borrow_mut().as_mut() {
            state.uncharged = state.uncharged.wrapping_add(units);
        }
    });
}

/// This thread's count as it stands, wrapping: the difference of two readings
/// (`wrapping_sub`) is what was counted between them, also where a
/// reservation took its price off before its operation ran.
pub(crate) fn uncharged() -> u64 {
    ACTIVE.with(|slot| slot.borrow().as_ref().map_or(0, State::counted))
}

// ---------------------------------------------------------------------------
// The count the number types keep
// ---------------------------------------------------------------------------
//
// Checkpoint 2 decisions, section 4.5, predictor step 1, asks for a complete
// work meter, and its done-when line reads wall time against work. The
// owner's ruling of 21 November 2026 is to count every ball and rational
// operation and to charge none of it. So the number types count themselves:
//
// - `Ball`: a sum or a difference counts 4L and a product, a square, a
//   reciprocal or a square root 8L^2, with L = ceil(p / 64) of the working
//   precision. A quotient is its reciprocal and its product.
// - `Rational`: a sum or a difference counts the two operands' limbs added, a
//   product or a quotient their limbs multiplied, and a decimal conversion the
//   value's limbs squared.
//
// Those are the decider's own prices for the same operations (`Work::ball`
// and `Work::exact` in `decider.rs`, under "Logical arithmetic costs"). They
// are the limb operations of schoolbook addition and multiplication (Knuth,
// The Art of Computer Programming, volume 2, section 4.3.1, algorithms A and
// M); a decimal conversion converts and divides, each quadratic in the value's
// limbs (the same section, algorithm D, and section 4.4). Every operation
// counts one unit more, because every charge of the decider also pays one
// unit for its deadline check (`Work::charge`).
//
// What is not counted, and why:
//
// - Arithmetic inside a decision, which is while a `Work` is alive. The
//   decision charges its arithmetic to `used` as it goes, and that charge is
//   its price. A decision's charges do not depend on which proofs are cached
//   or on which thread proved them (`ahead.rs`); the arithmetic it runs does.
// - Arithmetic inside an operation that is counted whole: a ball operation
//   and a decimal conversion are priced as one operation each.
// - The fill of a memo that lives on a shared value and is not charged
//   (`uncounted_memo`). Whichever thread or request reads the value first
//   fills it, once for the process, so a count of it would differ with the
//   thread count and with what ran before. The decider prices a read of such
//   a memo at nothing ("Construction-time caches spend no decision budget").
// - Comparisons, negations, shifts by a power of two and conversions between
//   the types. The decider has no price for them.
// - `Integer` arithmetic done directly. The decider prices rationals and
//   balls only.
//
// `Decider::reserve_exact_arithmetic` charges `used` for a rational operation
// its caller runs next, outside the decision. That operation counts itself,
// so the reservation takes the same amount off (`reserved_arithmetic`), and
// what was charged is not counted a second time.
//
// A `Decider::charge_work` is a caller's own tally of its loop and pays for no
// operation of the number types here. Arithmetic its caller runs is counted.

/// Counts one operation of the number types that costs `price` at the
/// decider's rate, unless the count is closed on this thread.
#[inline]
pub(crate) fn count_arithmetic(price: u64) {
    ARITHMETIC.with(|count| {
        let mut now = count.get();
        if now.closed == 0 {
            now.units = now.units.wrapping_add(price.saturating_add(1));
            count.set(now);
        }
    });
}

/// Takes `units` off the count: the price a reservation charged for an
/// operation that will count itself when it runs.
pub(crate) fn reserved_arithmetic(units: u64) {
    ARITHMETIC.with(|count| {
        let mut now = count.get();
        if now.closed == 0 {
            now.units = now.units.wrapping_sub(units);
            count.set(now);
        }
    });
}

/// Keeps the count closed on this thread while it lives.
pub(crate) struct Uncounted(());

impl Uncounted {
    #[inline]
    pub(crate) fn enter() -> Self {
        ARITHMETIC.with(|count| {
            let mut now = count.get();
            now.closed += 1;
            count.set(now);
        });
        Self(())
    }

    /// Counts one operation that costs `price` and closes the count for the
    /// arithmetic it is made of.
    #[inline]
    pub(crate) fn whole(price: u64) -> Self {
        count_arithmetic(price);
        Self::enter()
    }
}

impl Drop for Uncounted {
    #[inline]
    fn drop(&mut self) {
        ARITHMETIC.with(|count| {
            let mut now = count.get();
            now.closed -= 1;
            count.set(now);
        });
    }
}

/// Runs the fill of a memo that lives on a value shared between threads or
/// requests and whose fill is not charged. Its arithmetic is not counted in
/// [`OperationWork::uncharged`].
///
/// Such a memo is filled once, by whichever thread or request reads the value
/// first. A count of the fill would make a request's count depend on its
/// thread count and on what ran before it in the process, and the count is to
/// be the same at any thread count. Use it at every such memo whose fill runs
/// ball or rational arithmetic; a memo whose fill is charged is guarded with
/// [`crate::shared_caches`] instead, and its arithmetic is counted where it
/// is charged.
pub fn uncounted_memo<T>(fill: impl FnOnce() -> T) -> T {
    let _closed = Uncounted::enter();
    fill()
}

/// Whether `units` more fit the allowance. Unlike a charge, asking latches
/// nothing.
pub(crate) fn affordable(units: u64) -> bool {
    ACTIVE.with(|slot| {
        slot.borrow().as_ref().is_none_or(|state| {
            let work = state.work();
            !work.exhausted && units <= state.limit.saturating_sub(work.used)
        })
    })
}

/// Whether this thread is a worker running a loop item ahead of the writer.
pub(crate) fn worker() -> bool {
    ACTIVE.with(|slot| slot.borrow().as_ref().is_some_and(|state| state.worker))
}

/// The request's standing preferences, for its workers to copy.
pub(crate) fn preferences() -> (bool, bool) {
    ACTIVE.with(|slot| {
        slot.borrow().as_ref().map_or((false, false), |state| {
            (state.reduced_angle_series, state.reuse_enclosures)
        })
    })
}

fn clear_proofs() {
    ROOT_PROOFS.with(|proofs| proofs.borrow_mut().clear());
    ROOT_LIST_PROOFS.with(|proofs| proofs.borrow_mut().clear());
    ENCLOSURE_PROOFS.with(|proofs| proofs.borrow_mut().clear());
    PI_PROOFS.with(|proofs| proofs.borrow_mut().clear());
}

/// The writer's root proofs when it hands work to other threads. A worker
/// starts from them, so it does not prove again what the writer already holds.
#[derive(Clone, Default)]
pub(crate) struct Proofs {
    pub roots: RootProofs,
    lists: RootListProofs,
}

/// This thread's root proofs, for its workers to start from.
pub(crate) fn proofs() -> Proofs {
    Proofs {
        roots: ROOT_PROOFS.with(|proofs| proofs.borrow().clone()),
        lists: ROOT_LIST_PROOFS.with(|proofs| proofs.borrow().clone()),
    }
}

impl Proofs {
    /// The writer's certificate for this root at this precision, if it had one.
    pub(crate) fn root(
        &self,
        root: &Arc<AlgebraicRoot>,
        precision: u32,
    ) -> Option<Arc<RootEnclosure>> {
        self.roots.get(&RootProofKey::new(root, precision)).cloned()
    }
}

/// Gives a worker thread its own counter for one item, and proof caches of
/// its own that start as the writer's root proofs.
pub(crate) fn enter_worker(reduced_angle_series: bool, reuse_enclosures: bool, proofs: &Proofs) {
    install(State {
        request: NEXT_REQUEST.fetch_add(1, AtomicOrdering::Relaxed),
        aborted: Arc::new(AtomicBool::new(false)),
        limit: u64::MAX,
        used: 0,
        refunded: 0,
        uncharged: 0,
        exhausted: false,
        reduced_angle_series,
        reuse_enclosures,
        worker: true,
        ticks: 0,
    });
    clear_proofs();
    ROOT_PROOFS.with(|slot| *slot.borrow_mut() = proofs.roots.clone());
    ROOT_LIST_PROOFS.with(|slot| *slot.borrow_mut() = proofs.lists.clone());
}

pub(crate) fn leave_worker() {
    ACTIVE.with(|slot| *slot.borrow_mut() = None);
    clear_proofs();
}

pub(crate) fn exhausted() -> bool {
    current_operation_work().exhausted
}

pub(crate) fn used() -> u64 {
    current_operation_work().used
}

/// Only a certificate already charged in this request may bypass another
/// operation-wide reservation. Definition keys, not addresses or cache warmth,
/// identify the same certified root and the same published precision.
pub(crate) fn root_proof(root: &Arc<AlgebraicRoot>, precision: u32) -> Option<Arc<RootEnclosure>> {
    if ACTIVE.with(|slot| slot.borrow().is_none()) {
        return None;
    }
    let found = ROOT_PROOFS.with(|proofs| {
        proofs
            .borrow()
            .get(&RootProofKey {
                root: root.clone(),
                precision,
            })
            .cloned()
    });
    #[cfg(feature = "work-span")]
    if let Some(certificate) = &found {
        crate::work_span::read("proof cache: root", Arc::as_ptr(certificate) as u64);
    }
    found
}

pub(crate) fn record_root_proof(
    root: &Arc<AlgebraicRoot>,
    precision: u32,
    certificate: &Arc<RootEnclosure>,
) {
    if certificate.cost != 0 && ACTIVE.with(|slot| slot.borrow().is_some()) {
        ROOT_PROOFS.with(|proofs| {
            proofs.borrow_mut().insert(
                RootProofKey {
                    root: root.clone(),
                    precision,
                },
                certificate.clone(),
            );
        });
        #[cfg(feature = "work-span")]
        crate::work_span::write("proof cache: root", Arc::as_ptr(certificate) as u64);
    }
}

// Policy-bounded exact shortcuts can change isolation certificates and cold
// fees. Reuse only the same polynomial, policy and bracket-publication strategy.
pub(crate) fn root_list_proof(
    polynomial: &UPoly,
    policy: BudgetPolicy,
    coarse_brackets: bool,
) -> Option<Arc<CachedRoots>> {
    if ACTIVE.with(|slot| slot.borrow().is_none()) {
        return None;
    }
    crate::ahead::unrecordable();
    let found = ROOT_LIST_PROOFS.with(|proofs| {
        proofs
            .borrow()
            .get(&(policy.max_precision_bits, policy.max_work, coarse_brackets))
            .and_then(|polynomials| polynomials.get(polynomial))
            .cloned()
    });
    #[cfg(feature = "work-span")]
    if let Some(certificate) = &found {
        crate::work_span::read("proof cache: root list", Arc::as_ptr(certificate) as u64);
    }
    found
}

pub(crate) fn record_root_list_proof(
    polynomial: &Arc<UPoly>,
    policy: BudgetPolicy,
    coarse_brackets: bool,
    certificate: &Arc<CachedRoots>,
) {
    if ACTIVE.with(|slot| slot.borrow().is_some()) {
        ROOT_LIST_PROOFS.with(|proofs| {
            proofs
                .borrow_mut()
                .entry((policy.max_precision_bits, policy.max_work, coarse_brackets))
                .or_default()
                .insert(polynomial.clone(), certificate.clone());
        });
        #[cfg(feature = "work-span")]
        crate::work_span::write("proof cache: root list", Arc::as_ptr(certificate) as u64);
    }
}

#[cfg(test)]
mod allowance_tests {
    use super::*;

    fn record_nested_allowance() -> RecordedWork<()> {
        OperationTask::capture().run(|| {
            let allowance = OperationAllowance::reserve(5).unwrap();
            let (_, child) = allowance.run(|| crate::Decider::default().charge_work(3).unwrap());
            assert_eq!(child.used, 4);
            crate::Decider::default().charge_work(1).unwrap();
        })
    }

    #[test]
    fn prepaid_child_is_charged_once_when_committing_a_share() {
        let (_, work) = with_operation_budget(10, || {
            assert_eq!(record_nested_allowance().commit_share().unwrap(), Some(()));
            assert_eq!(current_operation_work().used, 7);
        });
        assert_eq!(work.used, 7);
        assert!(!work.exhausted);
    }

    #[test]
    fn prepaid_child_is_charged_once_when_committing_a_task() {
        let (_, work) = with_operation_budget(10, || {
            record_nested_allowance().commit().unwrap();
            assert_eq!(current_operation_work().used, 7);
        });
        assert_eq!(work.used, 7);
        assert!(!work.exhausted);
    }

    #[test]
    fn run_ahead_task_abort_stays_in_the_worker_budget() {
        let (_, writer) = with_operation_budget(10, || {
            let snapshot = proofs();
            std::thread::spawn(move || {
                enter_worker(false, false, &snapshot);
                assert!(
                    std::panic::catch_unwind(|| {
                        OperationTask::capture().run(|| {
                            crate::Decider::default().charge_work(3).unwrap();
                            panic!("abort a speculative worker task");
                        })
                    })
                    .is_err()
                );
                assert!(current_operation_work().exhausted);
                leave_worker();
            })
            .join()
            .unwrap();
            assert_eq!(current_operation_work().used, 0);
            assert!(!current_operation_work().exhausted);
            crate::Decider::default().charge_work(1).unwrap();
        });
        assert_eq!(writer.used, 2);
        assert!(!writer.exhausted);
    }

    #[test]
    fn a_recorded_result_cannot_escape_its_request() {
        let (recorded, _) = with_operation_budget(10, || OperationTask::capture().run(|| 7));
        assert_eq!(
            recorded.commit().unwrap_err().code,
            "numbers.operation-context"
        );
        let unbudgeted = OperationTask::capture().run(|| 8);
        let (_, work) = with_operation_budget(10, || {
            assert_eq!(
                unbudgeted.commit().unwrap_err().code,
                "numbers.operation-context"
            );
        });
        assert_eq!(work.used, 0);
    }

    #[test]
    fn shares_charge_actual_work_and_keep_the_total_cap() {
        let (_, work) = with_operation_budget(10, || {
            let tasks = OperationTask::capture().split(std::num::NonZeroUsize::new(4).unwrap());
            assert_eq!(
                tasks
                    .iter()
                    .map(|task| task.state.as_ref().unwrap().limit)
                    .sum::<u64>(),
                10
            );
            let results = std::thread::scope(|scope| {
                let handles = tasks
                    .into_iter()
                    .map(|task| {
                        scope.spawn(move || {
                            task.run(|| {
                                assert!(!charge(2));
                                42
                            })
                        })
                    })
                    .collect::<Vec<_>>();
                handles
                    .into_iter()
                    .map(|h| h.join().unwrap())
                    .collect::<Vec<_>>()
            });
            assert_eq!(current_operation_work().used, 0);
            for result in results {
                assert_eq!(result.commit_share().unwrap(), Some(42));
            }
        });
        assert_eq!(work.used, 8);
        assert!(!work.exhausted);
    }

    #[test]
    fn an_exhausted_share_is_paid_before_an_ordered_retry() {
        let (_, work) = with_operation_budget(10, || {
            let tasks = OperationTask::capture().split(std::num::NonZeroUsize::new(2).unwrap());
            for task in tasks {
                let recorded = task.run(|| {
                    assert!(!charge(3));
                    assert!(charge(3));
                    42
                });
                assert_eq!(recorded.commit_share().unwrap(), None);
            }
            assert_eq!(current_operation_work().used, 6);
            assert!(!current_operation_work().exhausted);
            assert_eq!(
                OperationTask::capture()
                    .run(|| {
                        assert!(!charge(3));
                        7
                    })
                    .commit()
                    .unwrap(),
                7
            );
            assert_eq!(current_operation_work().used, 9);
            let recorded = OperationTask::capture().run(|| {
                assert!(charge(3));
                8
            });
            assert_eq!(
                recorded.commit().unwrap_err().code,
                "numbers.operation-budget"
            );
        });
        assert_eq!(work.used, 9);
        assert!(work.exhausted);
    }

    fn caught_task_panics_exhaust_the_request(on_worker: bool) {
        let physical = Arc::new(AtomicU64::new(0));
        let (_, work) = with_operation_budget(10, || {
            for _ in 0..3 {
                let task = OperationTask::capture();
                let physical = Arc::clone(&physical);
                let run = move || {
                    task.run(|| {
                        if crate::Decider::default().charge_work(3).is_ok() {
                            physical.fetch_add(4, AtomicOrdering::Relaxed);
                        }
                        panic!("abort a numeric task");
                    })
                };
                let caught = if on_worker {
                    std::thread::spawn(run).join().is_err()
                } else {
                    std::panic::catch_unwind(run).is_err()
                };
                assert!(caught);
            }
            assert_eq!(physical.load(AtomicOrdering::Relaxed), 4);
            assert_eq!(current_operation_work().used, 10);
            assert!(current_operation_work().exhausted);
            assert!(crate::Decider::default().charge_work(1).is_err());
        });
        assert_eq!(work.used, 10);
        assert!(work.exhausted);
        let (_, fresh) = with_operation_budget(10, || {
            crate::Decider::default().charge_work(3).unwrap();
        });
        assert_eq!(fresh.used, 4);
        assert!(!fresh.exhausted);
    }

    #[test]
    fn caught_task_panic_on_the_caller_exhausts_the_request() {
        caught_task_panics_exhaust_the_request(false);
    }

    #[test]
    fn caught_task_panic_on_a_worker_exhausts_the_request() {
        caught_task_panics_exhaust_the_request(true);
    }

    #[test]
    fn split_reservations_conserve_work_and_restore_the_parent() {
        let (_, parent) = with_operation_budget(10, || {
            let mut allowance = OperationAllowance::reserve(7).unwrap();
            assert!(allowance.split_off(8).is_none());
            let second = allowance.split_off(4).unwrap();
            assert_eq!(allowance.limit(), 3);
            let (refused, first) = allowance.run(|| charge(4));
            assert!(refused);
            assert!(first.exhausted);
            assert_eq!(first.used, 0);
            assert_eq!(current_operation_work().used, 7);
            assert!(!current_operation_work().exhausted);
            let (_, actual) = second.run(|| assert!(!charge(4)));
            assert_eq!(actual.used, 4);
            assert!(!charge(3));
            assert!(OperationAllowance::reserve(1).is_err());
        });
        assert_eq!(parent.used, 10);
        assert!(parent.exhausted);
    }

    #[test]
    fn panic_in_reserved_work_restores_the_paid_parent() {
        let (_, parent) = with_operation_budget(10, || {
            let allowance = OperationAllowance::reserve(7).unwrap();
            let panic = std::panic::catch_unwind(|| {
                allowance.run(|| {
                    assert!(!charge(2));
                    panic!("unwind deferred work");
                })
            });
            assert!(panic.is_err());
            assert_eq!(current_operation_work().used, 7);
            assert!(!charge(3));
        });
        assert_eq!(parent.used, 10);
        assert!(!parent.exhausted);
        assert!(!operation_budget_active());
    }

    #[test]
    fn a_refused_attempt_returns_its_charges_and_a_success_keeps_them() {
        let (_, work) = with_operation_budget(10, || {
            assert!(!charge(2));
            // Exhausting the allowance inside a refused attempt leaves it whole.
            let refused: Result<(), ()> = refund_on_refusal(|| {
                assert!(!charge(5));
                assert!(charge(9));
                Err(())
            });
            assert!(refused.is_err());
            assert_eq!(current_operation_work().used, 2);
            assert!(!current_operation_work().exhausted);
            // The refused charge of 9 never fitted, so only the 5 was work.
            assert_eq!(current_operation_work().refunded, 5);
            let kept: Result<(), ()> = refund_on_refusal(|| {
                assert!(!charge(3));
                Ok(())
            });
            assert!(kept.is_ok());
            assert_eq!(current_operation_work().used, 5);
            assert_eq!(current_operation_work().refunded, 5);
            assert!(!charge(5));
        });
        assert_eq!(work.used, 10);
        assert_eq!(work.refunded, 5);
        assert!(!work.exhausted);
        assert_eq!(current_operation_work().refunded, 0);
    }

    #[test]
    fn nested_refused_attempts_are_each_counted_once() {
        let (_, work) = with_operation_budget(20, || {
            let outer: Result<(), ()> = refund_on_refusal(|| {
                assert!(!charge(2));
                let inner: Result<(), ()> = refund_on_refusal(|| {
                    assert!(!charge(3));
                    Err(())
                });
                assert!(inner.is_err());
                assert_eq!(current_operation_work().used, 2);
                assert_eq!(current_operation_work().refunded, 3);
                assert!(!charge(4));
                Err(())
            });
            assert!(outer.is_err());
        });
        assert_eq!(work.used, 0);
        assert_eq!(work.refunded, 9);
    }

    #[test]
    fn a_refused_attempt_cannot_revive_an_exhausted_request() {
        let (_, work) = with_operation_budget(3, || {
            assert!(charge(4));
            let refused: Result<(), ()> = refund_on_refusal(|| Err(()));
            assert!(refused.is_err());
        });
        assert!(work.exhausted);
        assert_eq!(work.refunded, 0);
    }

    #[test]
    fn a_refused_attempt_leaves_no_reservation_to_replay() {
        let (_, parent) = with_operation_budget(10, || {
            let recorded = OperationTask::capture().run(|| {
                assert!(!charge(1));
                let refused: Result<(), ()> = refund_on_refusal(|| {
                    assert!(!charge(6));
                    Err(())
                });
                assert!(refused.is_err());
                assert!(!charge(2));
            });
            recorded.commit().unwrap();
        });
        assert_eq!(parent.used, 3);
        assert_eq!(parent.refunded, 6);
        assert!(!parent.exhausted);
    }

    #[test]
    fn a_share_carries_its_refused_attempts_to_the_request() {
        let (_, parent) = with_operation_budget(10, || {
            let tasks = OperationTask::capture().split(std::num::NonZeroUsize::new(2).unwrap());
            for task in tasks {
                let recorded = task.run(|| {
                    let refused: Result<(), ()> = refund_on_refusal(|| {
                        assert!(!charge(4));
                        Err(())
                    });
                    assert!(refused.is_err());
                    assert!(!charge(1));
                });
                assert_eq!(recorded.commit_share().unwrap(), Some(()));
            }
        });
        assert_eq!(parent.used, 2);
        assert_eq!(parent.refunded, 8);
        assert!(!parent.exhausted);
    }

    // Two items are computed against a request of 10, five units each. Each
    // charges 1 and makes one refused attempt, of 2 and of 4. Before they are
    // admitted the request spends its 10, so the first admission refuses and
    // the second item is discarded, as `RecordedWork::commit` requires.
    fn discard_after_a_refused_admission(on_threads: bool) -> OperationWork {
        let (_, parent) = with_operation_budget(10, || {
            let item = |attempt: u64| {
                move || {
                    assert!(!charge(1));
                    let refused: Result<(), ()> = refund_on_refusal(|| {
                        assert!(!charge(attempt));
                        Err(())
                    });
                    assert!(refused.is_err());
                }
            };
            let tasks = OperationTask::capture().split(std::num::NonZeroUsize::new(2).unwrap());
            let mut recorded = tasks
                .into_iter()
                .zip([2, 4])
                .map(|(task, attempt)| {
                    if on_threads {
                        std::thread::scope(|scope| {
                            scope.spawn(move || task.run(item(attempt))).join().unwrap()
                        })
                    } else {
                        task.run(item(attempt))
                    }
                })
                .collect::<Vec<_>>()
                .into_iter();
            // Computing the items reported nothing to the request yet.
            assert_eq!(current_operation_work().refunded, 0);
            assert!(!charge(10));
            let first = recorded.next().unwrap();
            assert_eq!(first.commit().unwrap_err().code, "numbers.operation-budget");
            assert_eq!(current_operation_work().refunded, 2);
            drop(recorded);
        });
        parent
    }

    #[test]
    fn a_discarded_item_still_reports_its_refused_attempts() {
        for on_threads in [false, true] {
            let parent = discard_after_a_refused_admission(on_threads);
            // The refused attempts of both items: 2 from the one refused at
            // admission and 4 from the one that was discarded.
            assert_eq!(parent.refunded, 6);
            // Control: neither item's own charge of 1 was admitted, and the
            // request is exhausted exactly as it was before refunds were kept.
            assert_eq!(parent.used, 10);
            assert!(parent.exhausted);
        }
    }

    #[test]
    fn uncharged_arithmetic_is_counted_and_takes_nothing_from_the_allowance() {
        // No request, no counter.
        count_uncharged(7);
        assert_eq!(current_operation_work().uncharged, 0);
        let (_, work) = with_operation_budget(10, || {
            // Far more than the allowance holds, and nothing is exhausted.
            count_uncharged(1_000);
            assert!(!charge(4));
            count_uncharged(5);
            assert_eq!(current_operation_work().uncharged, 1_005);
            // A refused attempt gives its charges back. What it counted was
            // never charged, so it stays counted.
            let refused: Result<(), ()> = refund_on_refusal(|| {
                assert!(!charge(3));
                count_uncharged(20);
                Err(())
            });
            assert!(refused.is_err());
            // Counting goes on after the allowance is spent.
            assert!(charge(7));
            count_uncharged(100);
        });
        assert_eq!(work.used, 4);
        assert_eq!(work.refunded, 3);
        assert_eq!(work.uncharged, 1_125);
        assert!(work.exhausted);
        assert_eq!(current_operation_work().uncharged, 0);
    }

    #[test]
    fn a_reserved_allowance_reports_its_own_uncharged_arithmetic() {
        let (_, parent) = with_operation_budget(10, || {
            count_uncharged(2);
            let (_, child) = OperationAllowance::reserve(5).unwrap().run(|| {
                count_uncharged(30);
                assert_eq!(current_operation_work().uncharged, 30);
            });
            assert_eq!(child.uncharged, 30);
            assert_eq!(current_operation_work().uncharged, 2);
        });
        assert_eq!(parent.used, 5);
        assert_eq!(parent.uncharged, 2);
    }

    // Three items of one request count 2, 4 and 8. The first is admitted, the
    // second is refused at admission and the third is discarded.
    #[test]
    fn a_recorded_item_reports_its_uncharged_arithmetic_once() {
        for on_threads in [false, true] {
            let (_, parent) = with_operation_budget(10, || {
                let item = |units: u64| {
                    move || {
                        assert!(!charge(1));
                        count_uncharged(units);
                    }
                };
                let tasks = OperationTask::capture().split(std::num::NonZeroUsize::new(3).unwrap());
                let mut recorded = tasks
                    .into_iter()
                    .zip([2, 4, 8])
                    .map(|(task, units)| {
                        if on_threads {
                            std::thread::scope(|scope| {
                                scope.spawn(move || task.run(item(units))).join().unwrap()
                            })
                        } else {
                            task.run(item(units))
                        }
                    })
                    .collect::<Vec<_>>()
                    .into_iter();
                // Computing the items reported nothing to the request yet.
                assert_eq!(current_operation_work().uncharged, 0);
                recorded.next().unwrap().commit().unwrap();
                assert_eq!(current_operation_work().uncharged, 2);
                assert!(!charge(9));
                let second = recorded.next().unwrap();
                assert_eq!(
                    second.commit().unwrap_err().code,
                    "numbers.operation-budget"
                );
                assert_eq!(current_operation_work().uncharged, 6);
                drop(recorded);
            });
            assert_eq!(parent.uncharged, 14, "on threads: {on_threads}");
            assert_eq!(parent.used, 10);
            assert!(parent.exhausted);
        }
    }

    #[test]
    fn a_share_carries_its_uncharged_arithmetic_to_the_request() {
        let (_, parent) = with_operation_budget(10, || {
            let tasks = OperationTask::capture().split(std::num::NonZeroUsize::new(2).unwrap());
            for task in tasks {
                let recorded = task.run(|| {
                    count_uncharged(6);
                    assert!(!charge(1));
                });
                assert_eq!(recorded.commit_share().unwrap(), Some(()));
            }
        });
        assert_eq!(parent.used, 2);
        assert_eq!(parent.uncharged, 12);
    }

    #[test]
    fn a_record_outside_its_request_reports_nothing() {
        let attempt = || {
            let refused: Result<(), ()> = refund_on_refusal(|| {
                assert!(!charge(4));
                Err(())
            });
            assert!(refused.is_err());
        };
        let (recorded, first) = with_operation_budget(10, || OperationTask::capture().run(attempt));
        assert_eq!(first.refunded, 0);
        let (_, other) = with_operation_budget(10, || drop(recorded));
        assert_eq!(other.refunded, 0);
        let (recorded, first) =
            with_operation_budget(10, || OperationTask::capture().run(|| count_uncharged(4)));
        assert_eq!(first.uncharged, 0);
        let (_, other) = with_operation_budget(10, || drop(recorded));
        assert_eq!(other.uncharged, 0);
    }
}
