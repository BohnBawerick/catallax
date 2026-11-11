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
        atomic::{AtomicU64, Ordering as AtomicOrdering},
    },
};

static NEXT_REQUEST: AtomicU64 = AtomicU64::new(1);

#[derive(Clone, Copy)]
struct State {
    request: u64,
    limit: u64,
    used: u64,
    exhausted: bool,
    reduced_angle_series: bool,
    reuse_enclosures: bool,
}

struct RootProofKey {
    root: Arc<AlgebraicRoot>,
    precision: u32,
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

type RootListProofs = BTreeMap<(u32, u64), HashMap<Arc<UPoly>, Arc<CachedRoots>>>;
type EnclosureProofs = HashMap<(*const crate::real::Node, u32), (Arc<crate::real::Node>, Rc<Ball>)>;

thread_local! {
    static ACTIVE: Cell<Option<State>> = const { Cell::new(None) };
    static RESERVATIONS: RefCell<Option<Vec<u64>>> = const { RefCell::new(None) };
    static ROOT_PROOFS: RefCell<BTreeMap<RootProofKey, Arc<RootEnclosure>>> =
        const { RefCell::new(BTreeMap::new()) };
    static ROOT_LIST_PROOFS: RefCell<RootListProofs> =
        const { RefCell::new(BTreeMap::new()) };
    static ENCLOSURE_PROOFS: RefCell<EnclosureProofs> = RefCell::new(HashMap::new());
    static PI_PROOFS: RefCell<BTreeMap<u32, Ball>> = const { RefCell::new(BTreeMap::new()) };
}

/// Accounting for one complete operation, including all its deciders.
#[derive(Clone, Copy, Debug)]
pub struct OperationWork {
    /// Accepted arithmetic reservations and cancellation checkpoints.
    pub used: u64,
    /// Whether any reservation exceeded the operation allowance.
    pub exhausted: bool,
}

/// Whether the caller already owns a request-scoped work allowance.
/// Consumers must not replace it with a fresh allowance for a nested build.
pub fn operation_budget_active() -> bool {
    ACTIVE.with(|slot| slot.get().is_some())
}

/// Reuse numeric enclosure proofs for the rest of this operation.
/// Each operation starts with reuse disabled, including reserved allowances.
pub fn reuse_operation_enclosures() {
    ACTIVE.with(|slot| {
        if let Some(mut state) = slot.get() {
            state.reuse_enclosures = true;
            slot.set(Some(state));
        }
    });
}

pub(crate) fn enclosure_reuse_active() -> bool {
    ACTIVE.with(|slot| slot.get().is_some_and(|state| state.reuse_enclosures))
}

/// Prefer repeated half-angle reductions for the rest of this operation.
/// Enclosing balls retain their working precision and outward guarantee.
/// This preference is inactive outside an operation budget.
pub fn prefer_reduced_angle_series() {
    ACTIVE.with(|slot| {
        if let Some(mut state) = slot.get() {
            state.reduced_angle_series = true;
            slot.set(Some(state));
        }
    });
}

pub(crate) fn angle_reductions() -> u32 {
    ACTIVE.with(|slot| {
        if slot.get().is_some_and(|state| state.reduced_angle_series) {
            4
        } else {
            1
        }
    })
}

/// Observes accepted work in the current request without charging it.
/// Outside a request budget both counters are zero/false.
pub fn current_operation_work() -> OperationWork {
    ACTIVE.with(|slot| {
        slot.get().map_or(
            OperationWork {
                used: 0,
                exhausted: false,
            },
            |state| OperationWork {
                used: state.used,
                exhausted: state.exhausted,
            },
        )
    })
}

/// Runs a synchronous operation with a deterministic work allowance.
///
/// All deciders on this thread share the allowance. Numeric work on other
/// threads must use [`OperationTask`] and commit in input order. Nested scopes are rejected. Panic unwinding
/// restores the thread, so a later request cannot inherit another request's work.
/// Callers must discard the returned value when `OperationWork::exhausted` is
/// true: lower layers may have handled a refusal and continued collecting data.
pub fn with_operation_budget<T>(limit: u64, run: impl FnOnce() -> T) -> (T, OperationWork) {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            ACTIVE.with(|slot| slot.set(None));
            ROOT_PROOFS.with(|proofs| proofs.borrow_mut().clear());
            ROOT_LIST_PROOFS.with(|proofs| proofs.borrow_mut().clear());
            ENCLOSURE_PROOFS.with(|proofs| proofs.borrow_mut().clear());
            PI_PROOFS.with(|proofs| proofs.borrow_mut().clear());
        }
    }
    ACTIVE.with(|slot| {
        assert!(slot.get().is_none(), "operation budgets cannot nest");
        slot.set(Some(State {
            request: NEXT_REQUEST.fetch_add(1, AtomicOrdering::Relaxed),
            limit,
            used: 0,
            exhausted: false,
            reduced_angle_series: false,
            reuse_enclosures: false,
        }));
    });
    ROOT_PROOFS.with(|proofs| proofs.borrow_mut().clear());
    ROOT_LIST_PROOFS.with(|proofs| proofs.borrow_mut().clear());
    ENCLOSURE_PROOFS.with(|proofs| proofs.borrow_mut().clear());
    PI_PROOFS.with(|proofs| proofs.borrow_mut().clear());
    let _reset = Reset;
    let result = run();
    let state = ACTIVE.with(|slot| slot.get().unwrap());
    (
        result,
        OperationWork {
            used: state.used,
            exhausted: state.exhausted,
        },
    )
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
            state: ACTIVE.with(|slot| slot.get()).map(|state| State {
                limit: state.limit.saturating_sub(state.used),
                used: 0,
                ..state
            }),
        }
    }

    /// Divides the captured cap between fixed input items. The first remainder
    /// items receive one extra unit. Shares sum to exactly the captured cap.
    pub fn split(self, count: std::num::NonZeroUsize) -> Vec<Self> {
        let count = count.get();
        (0..count)
            .map(|index| Self {
                state: self.state.map(|state| State {
                    limit: state.limit / count as u64
                        + u64::from((index as u64) < state.limit % count as u64),
                    ..state
                }),
            })
            .collect()
    }

    /// Computes one independent item. Each item starts with empty proof caches,
    /// even when it executes on the caller or a reused pool thread.
    pub fn run<T>(self, run: impl FnOnce() -> T) -> RecordedWork<T> {
        let Some(state) = self.state else {
            return RecordedWork {
                value: run(),
                reservations: Vec::new(),
                request: None,
                work: None,
            };
        };
        struct Restore(Option<Vec<u64>>);
        impl Drop for Restore {
            fn drop(&mut self) {
                RESERVATIONS.with(|slot| *slot.borrow_mut() = self.0.take());
            }
        }
        let _restore = Restore(RESERVATIONS.with(|slot| slot.replace(Some(Vec::new()))));
        let (value, work) = OperationAllowance { limit: state.limit }.run(|| {
            ACTIVE.with(|slot| slot.set(Some(state)));
            run()
        });
        let reservations = RESERVATIONS.with(|slot| slot.take().unwrap());
        RecordedWork {
            value,
            reservations,
            request: Some(state.request),
            work: Some(work),
        }
    }
}

impl<T> RecordedWork<T> {
    fn check_request(&self) -> Result<(), crate::Refusal> {
        if ACTIVE.with(|slot| slot.get().map(|state| state.request)) != self.request {
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
            roots: BTreeMap<RootProofKey, Arc<RootEnclosure>>,
            lists: RootListProofs,
            enclosures: EnclosureProofs,
            pi: BTreeMap<u32, Ball>,
        }
        impl Drop for Restore {
            fn drop(&mut self) {
                ACTIVE.with(|slot| slot.set(self.state));
                ROOT_PROOFS.with(|slot| *slot.borrow_mut() = std::mem::take(&mut self.roots));
                ROOT_LIST_PROOFS.with(|slot| *slot.borrow_mut() = std::mem::take(&mut self.lists));
                ENCLOSURE_PROOFS
                    .with(|slot| *slot.borrow_mut() = std::mem::take(&mut self.enclosures));
                PI_PROOFS.with(|slot| *slot.borrow_mut() = std::mem::take(&mut self.pi));
            }
        }
        let _restore = Restore {
            state: ACTIVE.with(|slot| slot.take()),
            roots: ROOT_PROOFS.with(|slot| std::mem::take(&mut *slot.borrow_mut())),
            lists: ROOT_LIST_PROOFS.with(|slot| std::mem::take(&mut *slot.borrow_mut())),
            enclosures: ENCLOSURE_PROOFS.with(|slot| std::mem::take(&mut *slot.borrow_mut())),
            pi: PI_PROOFS.with(|slot| std::mem::take(&mut *slot.borrow_mut())),
        };
        with_operation_budget(self.limit, run)
    }
}

// Only known-defined, unbounded enclosures use this cache. Bounded decisions
// still pay their complete cold proof, and each request starts empty.
pub(crate) fn enclosure_proof(value: &Real, precision: u32) -> Option<Rc<Ball>> {
    if !enclosure_reuse_active() {
        return None;
    }
    ENCLOSURE_PROOFS.with(|proofs| {
        proofs
            .borrow()
            .get(&(Arc::as_ptr(&value.node), precision))
            .map(|(_, ball)| ball.clone())
    })
}

pub(crate) fn record_enclosure_proof(value: &Real, precision: u32, ball: &Rc<Ball>) {
    if !enclosure_reuse_active() {
        return;
    }
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
}

// Unlike expression enclosures, this constant has no definition or policy key:
// its completed ball depends only on precision. Reuse it in every request,
// with the same bounded cache and request/deferred-allowance lifetime as before.
pub(crate) fn pi_proof(precision: u32) -> Option<Ball> {
    if !operation_budget_active() {
        return None;
    }
    PI_PROOFS.with(|proofs| proofs.borrow().get(&precision).cloned())
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
    }
}

pub(crate) fn charge(units: u64) -> bool {
    RESERVATIONS.with(|slot| {
        if let Some(reservations) = &mut *slot.borrow_mut() {
            reservations.push(units);
        }
    });
    ACTIVE.with(|slot| {
        let Some(mut state) = slot.get() else {
            return false;
        };
        if units > state.limit.saturating_sub(state.used) {
            state.exhausted = true;
        } else if !state.exhausted {
            state.used += units;
        }
        slot.set(Some(state));
        state.exhausted
    })
}

pub(crate) fn exhausted() -> bool {
    ACTIVE.with(|slot| slot.get().is_some_and(|state| state.exhausted))
}

pub(crate) fn used() -> u64 {
    ACTIVE.with(|slot| slot.get().map_or(0, |state| state.used))
}

/// Only a certificate already charged in this request may bypass another
/// operation-wide reservation. Definition keys, not addresses or cache warmth,
/// identify the same certified root and the same published precision.
pub(crate) fn root_proof(root: &Arc<AlgebraicRoot>, precision: u32) -> Option<Arc<RootEnclosure>> {
    if ACTIVE.with(|slot| slot.get().is_none()) {
        return None;
    }
    ROOT_PROOFS.with(|proofs| {
        proofs
            .borrow()
            .get(&RootProofKey {
                root: root.clone(),
                precision,
            })
            .cloned()
    })
}

pub(crate) fn record_root_proof(
    root: &Arc<AlgebraicRoot>,
    precision: u32,
    certificate: &Arc<RootEnclosure>,
) {
    if certificate.cost != 0 && ACTIVE.with(|slot| slot.get().is_some()) {
        ROOT_PROOFS.with(|proofs| {
            proofs.borrow_mut().insert(
                RootProofKey {
                    root: root.clone(),
                    precision,
                },
                certificate.clone(),
            );
        });
    }
}

// Policy-bounded exact shortcuts can change isolation certificates and cold
// fees. Reuse a complete root list only for the same polynomial and policy.
pub(crate) fn root_list_proof(
    polynomial: &UPoly,
    policy: BudgetPolicy,
) -> Option<Arc<CachedRoots>> {
    if ACTIVE.with(|slot| slot.get().is_none()) {
        return None;
    }
    ROOT_LIST_PROOFS.with(|proofs| {
        proofs
            .borrow()
            .get(&(policy.max_precision_bits, policy.max_work))
            .and_then(|polynomials| polynomials.get(polynomial))
            .cloned()
    })
}

pub(crate) fn record_root_list_proof(
    polynomial: &Arc<UPoly>,
    policy: BudgetPolicy,
    certificate: &Arc<CachedRoots>,
) {
    if ACTIVE.with(|slot| slot.get().is_some()) {
        ROOT_LIST_PROOFS.with(|proofs| {
            proofs
                .borrow_mut()
                .entry((policy.max_precision_bits, policy.max_work))
                .or_default()
                .insert(polynomial.clone(), certificate.clone());
        });
    }
}

#[cfg(test)]
mod allowance_tests {
    use super::*;

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
                    .map(|task| task.state.unwrap().limit)
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

    #[test]
    fn recording_on_the_caller_restores_its_allowance_after_panic() {
        let (_, work) = with_operation_budget(10, || {
            assert!(!charge(1));
            let task = OperationTask::capture();
            assert!(
                std::panic::catch_unwind(|| task.run(|| {
                    assert!(!charge(3));
                    panic!("unwind a speculative item");
                }))
                .is_err()
            );
            assert_eq!(current_operation_work().used, 1);
            assert_eq!(
                OperationTask::capture()
                    .run(|| {
                        assert!(!charge(4));
                        7
                    })
                    .commit()
                    .unwrap(),
                7
            );
            assert_eq!(current_operation_work().used, 5);
        });
        assert_eq!(work.used, 5);
        assert!(!work.exhausted);
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
