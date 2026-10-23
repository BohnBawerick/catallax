//! Observations only. No counter read changes arithmetic or the request allowance.
use crate::FaceId;
use catallax_numbers::{Decider, current_operation_work};

/// Deterministic cost of a set of executed display operations.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Cost {
    /// Visited cells, or calls when this is a stage/query cost.
    pub cells: usize,
    /// Highest parameter-cell depth visited by these operations.
    pub max_depth: u32,
    /// Decisions made by the caller's decider.
    pub decisions: u64,
    /// Complete-request work, including checks and work from internal deciders.
    pub work: u64,
}
/// Cost of one source face's boundary queries. Shared support work is separate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FaceWork {
    /// Registry handle in the evaluated plan.
    pub face: FaceId,
    /// Cells visited by the shared support cover containing this source face.
    pub candidate_cells: usize,
    /// Maximum depth in that shared support cover.
    pub max_depth: u32,
    /// Executed boundary-query costs attributed to this face, without shared work.
    pub queries: Cost,
}
/// Work on one canonical support, shared by every listed source face.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SupportWork {
    /// Coactive source-face candidates. Do not duplicate `total` for each name.
    pub faces: Vec<FaceId>,
    /// All traversal work for this support, including failed cells.
    pub total: Cost,
    /// Regional trim-classification work, included in `total`.
    pub classification: Cost,
    /// Boundary-query work, included in `total` and attributed by `FaceWork`.
    pub queries: Cost,
    /// Coordinate generation and triangle publication work, included in `total`.
    pub emission: Cost,
}
/// Partial deterministic observations survive a meshing refusal, never geometry.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Diagnostics {
    /// Full mesh attempt, including setup and a refusal's executed prefix.
    pub total: Cost,
    /// Common preparation and shared indexing, not assigned to a face.
    pub common: Cost,
    /// Registry-order face costs. Unvisited faces have zero cost.
    pub faces: Vec<FaceWork>,
    /// Canonical support order, with shared costs recorded once.
    pub supports: Vec<SupportWork>,
}
#[derive(Clone, Copy)]
pub(super) struct Snapshot {
    decisions: u64,
    work: u64,
}
impl Snapshot {
    pub(super) fn new(dec: &Decider) -> Self {
        Self {
            decisions: dec.stats().decisions,
            work: current_operation_work().used,
        }
    }
}
impl Cost {
    pub(super) fn record(&mut self, before: Snapshot, dec: &Decider, depth: u32, count_cell: bool) {
        self.cells += usize::from(count_cell);
        self.max_depth = self.max_depth.max(depth);
        self.decisions += dec.stats().decisions.saturating_sub(before.decisions);
        self.work += current_operation_work().used.saturating_sub(before.work);
    }
}
