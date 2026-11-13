//! Request-scoped deterministic accounting, independent of wall time.

use crate::{
    Ball, BudgetPolicy, Real, UPoly,
    algebraic::{AlgebraicRoot, RootEnclosure},
    decider::CachedRoots,
};
use std::{
    cell::RefCell,
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
type RootListProofs = im::OrdMap<(u32, u64), im::HashMap<Arc<UPoly>, Arc<CachedRoots>>>;
type EnclosureProofs = HashMap<(*const crate::real::Node, u32), (Arc<crate::real::Node>, Rc<Ball>)>;

thread_local! {
    static ACTIVE: RefCell<Option<State>> = const { RefCell::new(None) };
    static RESERVATIONS: RefCell<Option<Vec<u64>>> = const { RefCell::new(None) };
    static ROOT_PROOFS: RefCell<RootProofs> = RefCell::new(RootProofs::new());
    static ROOT_LIST_PROOFS: RefCell<RootListProofs> = RefCell::new(RootListProofs::new());
    static ENCLOSURE_PROOFS: RefCell<EnclosureProofs> = RefCell::new(HashMap::new());
    static PI_PROOFS: RefCell<BTreeMap<u32, Ball>> = const { RefCell::new(BTreeMap::new()) };
}

/// Accounting for one complete operation, including all its deciders.
#[derive(Clone, Copy, Debug)]
pub struct OperationWork {
    /// Accepted arithmetic reservations and cancellation checkpoints.
    /// An aborted task conservatively consumes the request's remaining cap.
    pub used: u64,
    /// Whether a reservation exceeded the allowance or a task aborted.
    pub exhausted: bool,
}

impl State {
    fn work(&self) -> OperationWork {
        let aborted = self.aborted.load(AtomicOrdering::Acquire);
        OperationWork {
            used: if aborted { self.limit } else { self.used },
            exhausted: self.exhausted || aborted,
        }
    }
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
/// Outside a request budget both counters are zero/false.
pub fn current_operation_work() -> OperationWork {
    ACTIVE.with(|slot| {
        slot.borrow().as_ref().map_or_else(
            || OperationWork {
                used: 0,
                exhausted: false,
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
    ACTIVE.with(|slot| {
        let mut active = slot.borrow_mut();
        assert!(active.is_none(), "operation budgets cannot nest");
        *active = Some(State {
            request: NEXT_REQUEST.fetch_add(1, AtomicOrdering::Relaxed),
            aborted: Arc::new(AtomicBool::new(false)),
            limit,
            used: 0,
            exhausted: false,
            reduced_angle_series: false,
            reuse_enclosures: false,
            worker,
            ticks: 0,
        });
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
pub struct RecordedWork<T> {
    value: T,
    reservations: Vec<u64>,
    request: Option<u64>,
    work: Option<OperationWork>,
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
        let ((value, reservations), work) = OperationAllowance { limit: state.limit }.run(|| {
            ACTIVE.with(|slot| *slot.borrow_mut() = Some(state));
            let _restore = Restore(RESERVATIONS.with(|slot| slot.replace(Some(Vec::new()))));
            let value = run();
            let reservations = RESERVATIONS.with(|slot| slot.take().unwrap());
            (value, reservations)
        });
        RecordedWork {
            value,
            reservations,
            request: Some(request),
            work: Some(work),
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
                ACTIVE.with(|slot| *slot.borrow_mut() = self.state.take());
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
            state: ACTIVE.with(|slot| slot.take()),
            reservations: RESERVATIONS.with(|slot| slot.take()),
            roots: ROOT_PROOFS.with(|slot| std::mem::take(&mut *slot.borrow_mut())),
            lists: ROOT_LIST_PROOFS.with(|slot| std::mem::take(&mut *slot.borrow_mut())),
            enclosures: ENCLOSURE_PROOFS.with(|slot| std::mem::take(&mut *slot.borrow_mut())),
            pi: PI_PROOFS.with(|slot| std::mem::take(&mut *slot.borrow_mut())),
        };
        budget_scope(self.limit, worker, run)
    }
}

// Only known-defined, unbounded enclosures use this cache. Bounded decisions
// still pay their complete cold proof, and each request starts empty.
pub(crate) fn enclosure_proof(value: &Real, precision: u32) -> Option<Rc<Ball>> {
    if !enclosure_reuse_active() {
        return None;
    }
    crate::ahead::unrecordable();
    let found = ENCLOSURE_PROOFS.with(|proofs| {
        proofs
            .borrow()
            .get(&(Arc::as_ptr(&value.node), precision))
            .map(|(_, ball)| ball.clone())
    });
    #[cfg(feature = "work-span")]
    if let Some(ball) = &found {
        crate::work_span::read("proof cache: enclosure", Rc::as_ptr(ball) as u64);
    }
    found
}

pub(crate) fn record_enclosure_proof(value: &Real, precision: u32, ball: &Rc<Ball>) {
    if !enclosure_reuse_active() {
        return;
    }
    crate::ahead::unrecordable();
    ENCLOSURE_PROOFS.with(|proofs| {
        let mut proofs = proofs.borrow_mut();
        if proofs.len() >= 8192 {
            proofs.clear();
        }
        // Pin the node so its address cannot name another definition.
        proofs.insert(
            (Arc::as_ptr(&value.node), precision),
            (value.node.clone(), ball.clone()),
        );
    });
    #[cfg(feature = "work-span")]
    crate::work_span::write("proof cache: enclosure", Rc::as_ptr(ball) as u64);
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
    let mut observed = None;
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
                observed = Some(units);
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
    // The observer queries worker state, so release the counter's mutable borrow first.
    #[cfg(feature = "work-span")]
    if let Some(units) = observed {
        crate::work_span::advance(units);
    }
    exhausted
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
    ACTIVE.with(|slot| {
        *slot.borrow_mut() = Some(State {
            request: NEXT_REQUEST.fetch_add(1, AtomicOrdering::Relaxed),
            aborted: Arc::new(AtomicBool::new(false)),
            limit: u64::MAX,
            used: 0,
            exhausted: false,
            reduced_angle_series,
            reuse_enclosures,
            worker: true,
            ticks: 0,
        });
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
// fees. Reuse a complete root list only for the same polynomial and policy.
pub(crate) fn root_list_proof(
    polynomial: &UPoly,
    policy: BudgetPolicy,
) -> Option<Arc<CachedRoots>> {
    if ACTIVE.with(|slot| slot.borrow().is_none()) {
        return None;
    }
    crate::ahead::unrecordable();
    let found = ROOT_LIST_PROOFS.with(|proofs| {
        proofs
            .borrow()
            .get(&(policy.max_precision_bits, policy.max_work))
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
    certificate: &Arc<CachedRoots>,
) {
    if ACTIVE.with(|slot| slot.borrow().is_some()) {
        ROOT_LIST_PROOFS.with(|proofs| {
            proofs
                .borrow_mut()
                .entry((policy.max_precision_bits, policy.max_work))
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
}
