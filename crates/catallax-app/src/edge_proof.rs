use std::sync::atomic::AtomicBool;

use catallax_topology::Tolerance;
pub use catallax_topology::numerical::edge_proof::{
    EdgeProof, EdgeProofFailure, EdgeProofKey, EdgeProofOutcome, EdgeProofReport,
};
use catallax_topology::numerical::edge_proof::{prove_exact_edges, prove_numerical_edges};

use crate::EvaluatedApplicationOutcome;

/// Requests one whole-edge gap certificate per edge of a retained part.
///
/// A certificate bounds the maximum distance to all of the edge's face copies
/// and is no greater than `limit`. `work_cap` counts face uses on an exact part
/// and interval visits or projected endpoints on a numerical part. Completed
/// edges keep their certificates on cancellation; the rest report why they
/// could not be proved. Run this synchronous, read-only operation on a worker
/// thread when interactive. No request or report is retained by the kernel.
pub fn prove_edge_gaps(
    part: &EvaluatedApplicationOutcome,
    limit: &Tolerance,
    work_cap: u64,
    cancel: &AtomicBool,
) -> EdgeProofReport {
    match part {
        EvaluatedApplicationOutcome::Exact(body) => {
            let built = body.checkpoint.built();
            prove_exact_edges(&built.model, built.result, work_cap, cancel)
        }
        EvaluatedApplicationOutcome::Numerical(body) => {
            prove_numerical_edges(body.body(), limit, work_cap, cancel)
        }
    }
}
