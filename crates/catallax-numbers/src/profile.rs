//! Development-only attribution of request work to decision kinds.
//!
//! Compiled only with the `work-profile` feature. Observation reads the
//! request counter and never charges work, asks a numeric question or changes
//! an outcome. Profile one request at a time: rows from concurrent requests in
//! one process would be merged.

use crate::{Ball, Coincidence, DecisionStats, Outcome, Refusal, ZeroProof};
use std::{
    cell::RefCell,
    collections::BTreeMap,
    panic::Location,
    sync::{Mutex, PoisonError},
};

/// Work attributed to one decision kind inside one logged operation.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct WorkProfileRow {
    /// Innermost logged operation name, or empty outside any logged operation.
    pub section: String,
    /// Decider entry point, refined by outcome for `sign`.
    pub label: &'static str,
    /// Source location that called the decider, as `file:line`.
    pub caller: String,
    /// Completed calls.
    pub calls: u64,
    /// Request work spent inside the call, excluding nested decider calls.
    pub work: u64,
    /// For exact-zero signs, work up to the first completed refinement
    /// enclosure, or all of it when none was needed. Any bounded rule still
    /// needs that enclosure.
    pub first_enclosure_work: u64,
}

type Rows = BTreeMap<(String, &'static str, &'static str, u32), WorkProfileRow>;

static ROWS: Mutex<Option<Rows>> = Mutex::new(None);

struct Frame {
    label: &'static str,
    caller: &'static Location<'static>,
    start: u64,
    nested: u64,
    first_enclosure: Option<u64>,
}

thread_local! {
    static FRAMES: RefCell<Vec<Frame>> = const { RefCell::new(Vec::new()) };
    static SECTIONS: RefCell<Vec<String>> = const { RefCell::new(Vec::new()) };
}

/// Collects decider work attribution from every thread while `run` executes.
/// Nested profiles are rejected.
pub fn with_work_profile<T>(run: impl FnOnce() -> T) -> (T, Vec<WorkProfileRow>) {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            *ROWS.lock().unwrap_or_else(PoisonError::into_inner) = None;
        }
    }
    {
        let mut rows = ROWS.lock().unwrap_or_else(PoisonError::into_inner);
        assert!(rows.is_none(), "work profiles cannot nest");
        *rows = Some(Rows::new());
    }
    let _reset = Reset;
    let result = run();
    let rows = ROWS
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .take()
        .unwrap_or_default();
    (result, rows.into_values().collect())
}

/// One decider call. Dropping it records its exclusive work.
pub(crate) struct Span;

impl Span {
    #[track_caller]
    pub(crate) fn enter(label: &'static str) -> Self {
        let caller = Location::caller();
        FRAMES.with(|frames| {
            frames.borrow_mut().push(Frame {
                label,
                caller,
                start: crate::operation::used(),
                nested: 0,
                first_enclosure: None,
            });
        });
        Self
    }

    /// Replaces a `sign` label with its outcome class.
    pub(crate) fn classify<T>(
        &self,
        before: &DecisionStats,
        after: &DecisionStats,
        outcome: &Outcome<T>,
    ) {
        let label = match outcome {
            Outcome::Proven(_) if after.filter > before.filter => "sign.filter",
            Outcome::Proven(_) if after.refined > before.refined => "sign.refined",
            Outcome::Proven(_) => "sign.exact-nonzero",
            Outcome::ProvenCoincident(Coincidence { proof }) => match proof {
                ZeroProof::ExactArithmetic => "sign.zero.rational",
                ZeroProof::Algebraic => "sign.zero.algebraic",
                ZeroProof::SeparationBound => "sign.zero.separation",
                ZeroProof::ConstructionFact => "sign.zero.fact",
            },
            Outcome::NotExact(_) => "sign.not-exact",
            Outcome::Refused(_) => "sign.refused",
        };
        relabel(label);
    }

    /// Replaces an `enclose` label with whether the ball separated zero.
    /// A ball still containing zero is work spent near a possible zero.
    pub(crate) fn classify_enclosure(&self, result: &Result<Ball, Refusal>) {
        relabel(match result {
            Ok(ball) if ball.sign().is_some() => "enclose.separated",
            Ok(_) => "enclose.contains-zero",
            Err(_) => "enclose.refused",
        });
    }
}

fn relabel(label: &'static str) {
    FRAMES.with(|frames| {
        if let Some(frame) = frames.borrow_mut().last_mut() {
            frame.label = label;
        }
    });
}

impl Drop for Span {
    fn drop(&mut self) {
        let used = crate::operation::used();
        let Some(frame) = FRAMES.with(|frames| {
            let mut frames = frames.borrow_mut();
            let frame = frames.pop()?;
            let total = used.saturating_sub(frame.start);
            if let Some(parent) = frames.last_mut() {
                parent.nested += total;
            }
            Some(frame)
        }) else {
            return;
        };
        let total = used.saturating_sub(frame.start);
        // A zero settled before any refinement enclosure had nothing to skip.
        let first = if frame.label.starts_with("sign.zero.") {
            frame
                .first_enclosure
                .map_or(total, |at| at.saturating_sub(frame.start).min(total))
        } else {
            0
        };
        let section =
            SECTIONS.with(|sections| sections.borrow().last().cloned().unwrap_or_default());
        let mut rows = ROWS.lock().unwrap_or_else(PoisonError::into_inner);
        if let Some(rows) = rows.as_mut() {
            let row = rows
                .entry((
                    section.clone(),
                    frame.label,
                    frame.caller.file(),
                    frame.caller.line(),
                ))
                .or_insert_with(|| WorkProfileRow {
                    section,
                    label: frame.label,
                    caller: format!("{}:{}", frame.caller.file(), frame.caller.line()),
                    ..WorkProfileRow::default()
                });
            row.calls += 1;
            row.work += total.saturating_sub(frame.nested);
            row.first_enclosure_work += first;
        }
    }
}

/// Notes the first completed refinement enclosure of the innermost call.
pub(crate) fn first_enclosure() {
    let used = crate::operation::used();
    FRAMES.with(|frames| {
        if let Some(frame) = frames.borrow_mut().last_mut() {
            frame.first_enclosure.get_or_insert(used);
        }
    });
}

pub(crate) fn enter_section(operation: &str) {
    SECTIONS.with(|sections| sections.borrow_mut().push(operation.to_owned()));
}

pub(crate) fn leave_section() {
    SECTIONS.with(|sections| {
        sections.borrow_mut().pop();
    });
}

#[cfg(test)]
mod tests {
    use crate::{Decider, Rational, Real, with_operation_budget};

    fn real(value: i64) -> Real {
        Real::from_rational(Rational::from(value))
    }

    #[test]
    fn profile_attributes_all_request_work_without_charging_any() {
        let two = real(2).sqrt();
        let zero = &(&two * &two) - &real(2);
        let positive = &two - &real(1);
        let run = || {
            let mut dec = Decider::default();
            (dec.sign(&zero), dec.sign(&positive), dec.stats().clone())
        };
        let (plain, plain_work) = with_operation_budget(1_000_000, run);
        let ((profiled, work), rows) =
            super::with_work_profile(|| with_operation_budget(1_000_000, run));
        assert_eq!(plain, profiled);
        assert_eq!(plain_work.used, work.used);
        // Concurrent tests may record their own callers into the shared rows.
        let rows: Vec<_> = rows
            .into_iter()
            .filter(|row| row.caller.contains("profile.rs"))
            .collect();
        assert_eq!(rows.iter().map(|row| row.work).sum::<u64>(), work.used);
        let zero_row = rows
            .iter()
            .find(|row| row.label.starts_with("sign.zero."))
            .expect("exact zero row");
        assert_eq!(zero_row.calls, 1);
        assert!(zero_row.first_enclosure_work <= zero_row.work);
    }
}
