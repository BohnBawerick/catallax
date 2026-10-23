//! Request-scoped deterministic accounting, independent of wall time.

use std::cell::Cell;

#[derive(Clone, Copy)]
struct State {
    limit: u64,
    used: u64,
    exhausted: bool,
}

thread_local! {
    static ACTIVE: Cell<Option<State>> = const { Cell::new(None) };
}

/// Accounting for one complete operation, including all its deciders.
#[derive(Clone, Copy, Debug)]
pub struct OperationWork {
    /// Accepted arithmetic reservations and cancellation checkpoints.
    pub used: u64,
    /// Whether any reservation exceeded the operation allowance.
    pub exhausted: bool,
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

/// Remaining request allowance, when an operation scope is active.
///
/// Optional certificates use this to stop before they exhaust the caller's
/// budget, leaving room for a complete fallback proof path.
pub(crate) fn remaining() -> Option<u64> {
    ACTIVE.with(|slot| {
        slot.get().map(|state| {
            if state.exhausted {
                0
            } else {
                state.limit.saturating_sub(state.used)
            }
        })
    })
}
