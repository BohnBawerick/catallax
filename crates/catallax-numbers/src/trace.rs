//! Opt-in operation observations. Recording never asks a numeric question.
use crate::{Decider, DecisionStats, MarkLog, Refusal};
use std::{cell::RefCell, time::Instant};

thread_local! {
    static LOG: RefCell<Option<Vec<OperationEntry>>> = const { RefCell::new(None) };
}

/// One completed operation, or an interrupted operation during panic unwinding.
#[derive(Debug)]
pub struct OperationEntry {
    /// Implementation that executed the operation.
    pub kernel: &'static str,
    /// Caller-defined stable operation identifier.
    pub id: String,
    /// Schema operation name.
    pub operation: String,
    /// Decision counters before execution.
    pub before: DecisionStats,
    /// Decision counters after execution; absent when unwinding.
    pub after: Option<DecisionStats>,
    /// Marks on the operation's result, including inherited marks.
    pub marks: MarkLog,
    /// Geometric refusal, if any.
    pub refusal: Option<Refusal>,
    /// Request work count at operation entry.
    pub work_before: u64,
    /// Request work count at operation exit.
    pub work_after: u64,
    /// Observed elapsed time; never used for a decision.
    pub wall_ms: f64,
}

/// Records synchronous operations on this thread. Nested recorders are rejected.
/// Unwinding restores the recorder before the next request.
pub fn with_operation_log<T>(run: impl FnOnce() -> T) -> (T, Vec<OperationEntry>) {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            LOG.with(|log| *log.borrow_mut() = None);
        }
    }
    LOG.with(|log| {
        assert!(log.borrow().is_none(), "operation log scopes cannot nest");
        *log.borrow_mut() = Some(Vec::new());
    });
    let _reset = Reset;
    let result = run();
    let entries = LOG.with(|log| log.borrow_mut().take().unwrap());
    (result, entries)
}

/// A no-allocation disabled span. Only an active recorder copies names or stats.
pub struct OperationSpan(Option<(OperationEntry, Instant)>);
impl OperationSpan {
    /// Starts observing an operation without changing its arithmetic or budget.
    pub fn begin(kernel: &'static str, id: &impl ToString, operation: &str, dec: &Decider) -> Self {
        Self(LOG.with(|log| {
            log.borrow().as_ref().map(|_| {
                (
                    OperationEntry {
                        kernel,
                        id: id.to_string(),
                        operation: operation.to_owned(),
                        before: dec.stats().clone(),
                        after: None,
                        marks: MarkLog::new(),
                        refusal: None,
                        wall_ms: 0.0,
                        work_before: crate::operation::used(),
                        work_after: crate::operation::used(),
                    },
                    Instant::now(),
                )
            })
        }))
    }

    /// Completes the observation. No geometry or decision is recomputed.
    pub fn finish(mut self, dec: &Decider, refusal: Option<&Refusal>, marks: Option<&MarkLog>) {
        if let Some((entry, _)) = &mut self.0 {
            entry.after = Some(dec.stats().clone());
            entry.refusal = refusal.cloned();
            if let Some(marks) = marks {
                entry.marks = marks.clone();
            }
        }
    }
}
impl Drop for OperationSpan {
    fn drop(&mut self) {
        if let Some((mut entry, started)) = self.0.take() {
            entry.wall_ms = started.elapsed().as_secs_f64() * 1000.0;
            entry.work_after = crate::operation::used();
            LOG.with(|log| {
                if let Some(entries) = log.borrow_mut().as_mut() {
                    entries.push(entry);
                }
            });
        }
    }
}
