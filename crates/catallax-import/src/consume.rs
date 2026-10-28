//! The entry a recipe evaluator calls for an `import_step` step.

use catallax_recipe::StepId;
use catallax_topology::{BodyId, Model};
use sha2::{Digest, Sha256};

use crate::{ImportRefusal, build_file, read_step};

/// Lowercase hex SHA-256 of `bytes`.
pub fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

/// Builds the bodies of an unrepaired STEP file for the recipe step `step`.
///
/// `expected_sha256` is the hash saved in the step's `file`. Bytes that do not
/// match it are a recipe error. The geometry is exactly the file's: a file whose
/// faces do not meet exactly is `not_solid` until repaired, and nothing moves here.
pub fn import_step_bodies(
    model: &mut Model,
    bytes: &[u8],
    expected_sha256: &str,
    step: &StepId,
) -> Result<Vec<BodyId>, ImportRefusal> {
    let actual = sha256_hex(bytes);
    if actual != expected_sha256 {
        return Err(ImportRefusal::HashMismatch {
            expected: expected_sha256.to_string(),
            actual,
        });
    }
    build_file(model, &read_step(bytes)?, Some(step))
}
