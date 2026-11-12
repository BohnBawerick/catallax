//! Opt-in operation observations. Recording never asks a numeric question.
use crate::{Decider, DecisionStats, MarkLog, Refusal};
use std::{cell::RefCell, fmt, time::Instant};

thread_local! {
    static LOG: RefCell<Option<Vec<OperationEntry>>> = const { RefCell::new(None) };
    static PATH: RefCell<Option<&'static str>> = const { RefCell::new(None) };
}

/// One completed operation, or an interrupted operation during panic unwinding.
pub struct OperationEntry {
    /// Implementation that executed the operation.
    pub kernel: &'static str,
    /// Caller-defined stable operation identifier.
    pub id: String,
    /// Schema operation name.
    pub operation: String,
    /// Stable implementation path, when the operation records one.
    pub path: Option<&'static str>,
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

// Preserve the established operation-outcome dump. Path annotations remain
// available to the structured operation log through the public field.
impl fmt::Debug for OperationEntry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("OperationEntry")
            .field("kernel", &self.kernel)
            .field("id", &self.id)
            .field("operation", &self.operation)
            .field("before", &self.before)
            .field("after", &self.after)
            .field("marks", &self.marks)
            .field("refusal", &self.refusal)
            .field("work_before", &self.work_before)
            .field("work_after", &self.work_after)
            .field("wall_ms", &self.wall_ms)
            .finish()
    }
}

/// Records synchronous operations on this thread. Nested recorders are rejected.
/// Unwinding restores the recorder before the next request.
pub fn with_operation_log<T>(run: impl FnOnce() -> T) -> (T, Vec<OperationEntry>) {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            LOG.with(|log| *log.borrow_mut() = None);
            PATH.with(|path| *path.borrow_mut() = None);
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
pub struct OperationSpan(
    Option<(OperationEntry, Instant, Option<&'static str>)>,
    // Held only so the work-span region closes with the operation.
    #[cfg(feature = "work-span")]
    #[allow(dead_code)]
    Option<crate::work_span::Region>,
);
impl OperationSpan {
    /// Starts observing an operation without changing its arithmetic or budget.
    pub fn begin(kernel: &'static str, id: &impl ToString, operation: &str, dec: &Decider) -> Self {
        #[cfg(feature = "work-span")]
        let region = LOG.with(|log| {
            log.borrow().is_some().then(|| {
                crate::work_span::region(format!("{kernel} {operation} {}", id.to_string()))
            })
        });
        let entry = LOG.with(|log| {
            log.borrow().as_ref().map(|_| {
                #[cfg(feature = "work-profile")]
                crate::profile::enter_section(operation);
                (
                    OperationEntry {
                        kernel,
                        id: id.to_string(),
                        operation: operation.to_owned(),
                        path: None,
                        before: dec.stats().clone(),
                        after: None,
                        marks: MarkLog::new(),
                        refusal: None,
                        wall_ms: 0.0,
                        work_before: crate::operation::used(),
                        work_after: crate::operation::used(),
                    },
                    Instant::now(),
                    PATH.with(|path| path.replace(None)),
                )
            })
        });
        #[cfg(feature = "work-span")]
        return Self(entry, region);
        #[cfg(not(feature = "work-span"))]
        Self(entry)
    }

    /// Completes the observation. No geometry or decision is recomputed.
    pub fn finish(mut self, dec: &Decider, refusal: Option<&Refusal>, marks: Option<&MarkLog>) {
        if let Some((entry, _, _)) = &mut self.0 {
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
        if let Some((mut entry, started, previous)) = self.0.take() {
            #[cfg(feature = "work-profile")]
            crate::profile::leave_section();
            entry.path = PATH.with(|path| path.replace(previous));
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

/// Records the path used by the active operation without asking a numeric question.
/// Nested operation spans keep their paths separate. Disabled recording allocates nothing.
pub fn record_operation_path(path: &'static str) {
    LOG.with(|log| {
        if log.borrow().is_some() {
            PATH.with(|slot| *slot.borrow_mut() = Some(path));
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn path_annotations_are_nested_thread_local_and_spend_no_work() {
        let (((), entries), work) = crate::with_operation_budget(1, || {
            with_operation_log(|| {
                let dec = Decider::default();
                let outer = OperationSpan::begin("brep", &"outer", "shell", &dec);
                record_operation_path("distance.shell.sewn");
                std::thread::scope(|scope| {
                    scope.spawn(|| record_operation_path("other-thread"));
                });
                let inner = OperationSpan::begin("brep", &"inner", "boolean", &dec);
                record_operation_path("boolean.exact");
                inner.finish(&dec, None, None);
                outer.finish(&dec, None, None);
            })
        });
        assert_eq!(work.used, 0);
        assert_eq!(
            entries.iter().map(|entry| entry.path).collect::<Vec<_>>(),
            vec![Some("boolean.exact"), Some("distance.shell.sewn")]
        );
        record_operation_path("disabled");
        let ((), entries) = with_operation_log(|| {
            let dec = Decider::default();
            OperationSpan::begin("brep", &"next", "offset", &dec).finish(&dec, None, None);
        });
        assert_eq!(entries[0].path, None);
    }
}
