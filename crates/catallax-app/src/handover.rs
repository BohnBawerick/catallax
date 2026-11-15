//! Lossless transfer of built bodies between independent native or wasm kernels.
//!
//! These bytes retain a producer's state and evidence; they do not prove geometry
//! supplied by an untrusted producer. Version and checksum checks detect incompatible
//! or damaged transfers. Display meshes and incremental evaluation caches are separate.

use catallax_features::Built;
use catallax_numbers::{
    Decider, Refusal, operation_budget_active, transfer, with_operation_budget,
};
use catallax_topology::numerical::NumericalBody;

fn bounded<T>(run: impl FnOnce() -> Result<T, Refusal>) -> Result<T, Refusal> {
    if operation_budget_active() {
        return run();
    }
    let (answer, work) = with_operation_budget(crate::DEFAULT_BUILD_WORK, run);
    if work.exhausted {
        Err(Refusal::new(
            "handover.resource-limit",
            "body transfer exhausted its request allowance",
        ))
    } else {
        answer
    }
}

/// Writes the full model, persistent names, recipe and retained feature/evidence records.
/// Exact values retain their integer numerator/denominator and shared expression graph.
/// Charges the caller's active request allowance, or the application default.
pub fn write_body(body: &Built, decider: &mut Decider) -> Result<Vec<u8>, Refusal> {
    bounded(|| transfer::encode(body, decider))
}

/// Restores a build into independent stores without rebuilding, meshing or STEP import.
/// Picking handles have the same indices and belong to the returned model.
/// Incompatible, truncated or corrupt bytes return a typed refusal.
pub fn read_body(bytes: &[u8], decider: &mut Decider) -> Result<Built, Refusal> {
    bounded(|| transfer::decode(bytes, decider))
}

/// Writes numerical geometry with its definitions, tolerance, checks and repair history.
pub fn write_numerical(body: &NumericalBody, decider: &mut Decider) -> Result<Vec<u8>, Refusal> {
    bounded(|| transfer::encode(body, decider))
}

/// Restores numerical admission and evidence without promoting the result to exact.
pub fn read_numerical(bytes: &[u8], decider: &mut Decider) -> Result<NumericalBody, Refusal> {
    bounded(|| transfer::decode(bytes, decider))
}
