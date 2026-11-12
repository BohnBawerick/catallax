//! The entry a recipe evaluator calls for an `import_step` step.

use catallax_recipe::StepId;
use catallax_topology::{BodyId, Model};
use sha2::{Digest, Sha256};

use crate::{ImportRefusal, build_file, read_iges, read_step};

/// Lowercase hex SHA-256 of `bytes`.
pub fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

/// Builds the bodies of a STEP file for the recipe step `step`.
///
/// `expected_sha256` is the hash saved in the step's `file`. Bytes that do not
/// match it are a recipe error. Supported planar and capped-cylinder spline incidence repairs build
/// exact derivatives within the 10 um and neighbouring-feature limits, with
/// attained source-deviation bounds retained in the result bodies' marks.
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

/// Builds the bodies of an IGES file's manifold solid B-rep objects for the
/// recipe step `step`, exactly as [`import_step_bodies`] builds a STEP file's.
///
/// Face keys are directory-entry sequence numbers such as `"D17"`. The file is
/// read through the same exact incidence checks, and pairs of faces that pass
/// through each other are refused as for STEP. Nothing is repaired.
pub fn import_iges_bodies(
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
    let file = read_iges(bytes)?;
    crate::check_self_intersection(&file)?;
    build_file(model, &file, Some(step))
}
