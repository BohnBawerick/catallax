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
