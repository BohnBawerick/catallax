//! Request-scoped deterministic accounting, independent of wall time.

use crate::{
    BudgetPolicy, UPoly,
    algebraic::{AlgebraicRoot, RootEnclosure},
    decider::CachedRoots,
};
use std::{
    cell::{Cell, RefCell},
    cmp::Ordering,
    collections::{BTreeMap, HashMap},
    sync::Arc,
};

#[derive(Clone, Copy)]
struct State {
    limit: u64,
    used: u64,
    exhausted: bool,
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

thread_local! {
    static ACTIVE: Cell<Option<State>> = const { Cell::new(None) };
    static ROOT_PROOFS: RefCell<BTreeMap<RootProofKey, Arc<RootEnclosure>>> =
        const { RefCell::new(BTreeMap::new()) };
    static ROOT_LIST_PROOFS: RefCell<RootListProofs> =
        const { RefCell::new(BTreeMap::new()) };
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
/// All deciders on this thread share the allowance. The operation must not spawn
/// numeric work on other threads. Nested scopes are rejected. Panic unwinding
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
        }
    }
    ACTIVE.with(|slot| {
        assert!(slot.get().is_none(), "operation budgets cannot nest");
        slot.set(Some(State {
            limit,
            used: 0,
            exhausted: false,
        }));
    });
    ROOT_PROOFS.with(|proofs| proofs.borrow_mut().clear());
    ROOT_LIST_PROOFS.with(|proofs| proofs.borrow_mut().clear());
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
            return Err(crate::Refusal {
                code: "numbers.operation-budget",
                detail: "deferred work allocation exceeds the caller allowance".into(),
            });
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
        }
        impl Drop for Restore {
            fn drop(&mut self) {
                ACTIVE.with(|slot| slot.set(self.state));
                ROOT_PROOFS.with(|slot| *slot.borrow_mut() = std::mem::take(&mut self.roots));
                ROOT_LIST_PROOFS.with(|slot| *slot.borrow_mut() = std::mem::take(&mut self.lists));
            }
        }
        let _restore = Restore {
            state: ACTIVE.with(|slot| slot.take()),
            roots: ROOT_PROOFS.with(|slot| std::mem::take(&mut *slot.borrow_mut())),
            lists: ROOT_LIST_PROOFS.with(|slot| std::mem::take(&mut *slot.borrow_mut())),
        };
        with_operation_budget(self.limit, run)
    }
}

pub(crate) fn charge(units: u64) -> bool {
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
