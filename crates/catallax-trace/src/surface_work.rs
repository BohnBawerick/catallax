//! Work checkpoints enabled only for a metered surface intersection.

use std::{cell::Cell, time::Instant};

use crate::Fault;

/// Checkpoint units count visits, not floating-point instructions or time.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct WorkReport {
    /// Requested units, including a rejected reservation.
    pub attempted: u64,
    /// Units accepted by the caller's budget callback.
    pub accepted: u64,
    /// Search branches visited, including culled branches.
    pub search_branches: u64,
    /// Newton, Brent and line-search iterations visited.
    pub local_iterations: u64,
    /// Surface evaluation requests.
    pub surface_evaluations: u64,
    /// Curve evaluation requests.
    pub curve_evaluations: u64,
    /// Source control points and knots reserved before cache construction.
    pub source_items: u64,
    /// Trace steps and output-point visits.
    pub trace_visits: u64,
    /// The callback refused a reservation. No result may be retained.
    pub cancelled: bool,
}

#[derive(Clone, Copy)]
pub(crate) enum Kind {
    Search,
    Local,
    Surface,
    Curve,
    Source,
    Trace,
}

#[derive(Clone, Copy)]
struct State {
    callback: fn(u64, Option<Instant>) -> bool,
    deadline: Option<Instant>,
    report: WorkReport,
}

thread_local! {
    static ACTIVE: Cell<Option<State>> = const { Cell::new(None) };
}

pub(crate) fn checkpoint(kind: Kind, units: u64) -> Result<(), Fault> {
    ACTIVE.with(|active| {
        let Some(mut state) = active.get() else {
            return Ok(());
        };
        if state.report.cancelled {
            return Err(Fault::Degenerate);
        }
        let count = match kind {
            Kind::Search => &mut state.report.search_branches,
            Kind::Local => &mut state.report.local_iterations,
            Kind::Surface => &mut state.report.surface_evaluations,
            Kind::Curve => &mut state.report.curve_evaluations,
            Kind::Source => &mut state.report.source_items,
            Kind::Trace => &mut state.report.trace_visits,
        };
        *count = count.saturating_add(units);
        state.report.attempted = state.report.attempted.saturating_add(units);
        state.report.cancelled = !(state.callback)(units, state.deadline);
        if !state.report.cancelled {
            state.report.accepted = state.report.accepted.saturating_add(units);
        }
        active.set(Some(state));
        if state.report.cancelled {
            Err(Fault::Degenerate)
        } else {
            Ok(())
        }
    })
}

pub(crate) struct Scope;

impl Scope {
    pub(crate) fn enter(
        callback: fn(u64, Option<Instant>) -> bool,
        deadline: Option<Instant>,
    ) -> Result<Self, Fault> {
        ACTIVE.with(|active| {
            if active.get().is_some() {
                return Err(Fault::InvalidInput);
            }
            active.set(Some(State {
                callback,
                deadline,
                report: WorkReport::default(),
            }));
            Ok(Self)
        })
    }

    pub(crate) fn report(&self) -> WorkReport {
        ACTIVE.with(|active| active.get().map(|state| state.report).unwrap_or_default())
    }
}

impl Drop for Scope {
    fn drop(&mut self) {
        ACTIVE.with(|active| active.set(None));
    }
}
