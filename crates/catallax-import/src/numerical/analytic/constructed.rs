//! Constructed faces share import's sewing and finite-parameter reconstruction.

use super::*;
use crate::tolerant::{ToleranceOptions, healing::HealedFile};

/// A separately admitted numerical solid and the original sewing evidence.
/// A diagnostic definition fingerprint is not an imported STEP source identity.
#[derive(Debug)]
pub struct SewnSurfaceBody {
    body: n::NumericalBody,
    faces: std::sync::Arc<[crate::Face]>,
    healing: HealedFile,
    pair_checks: Vec<NumericalSourcePairCheck>,
    topology_repairs: Vec<NumericalTopologyRepair>,
}

impl SewnSurfaceBody {
    /// Checked closed material, always numerical even when all gaps are zero.
    pub fn body(&self) -> &n::NumericalBody {
        &self.body
    }
    /// Canonical 3D inputs used by the completed numerical reconstruction.
    /// Parameter-only seams remain described by `topology_repairs`.
    pub(crate) fn shared_faces(&self) -> std::sync::Arc<[crate::Face]> {
        self.faces.clone()
    }
    /// Immutable original surface definitions, seam matches and repair caps.
    pub fn healing(&self) -> &HealedFile {
        &self.healing
    }
    /// Actual source-pair checks, including explicitly unresolved conditions.
    pub fn source_pair_checks(&self) -> &[NumericalSourcePairCheck] {
        &self.pair_checks
    }
    /// Finite parameter reconstruction changes, separate from seam matching.
    pub fn topology_repairs(&self) -> &[NumericalTopologyRepair] {
        &self.topology_repairs
    }
}

/// Sew constructed finite analytic faces into a checked numerical solid.
/// This is the existing import sewer and carrier mapper, not a second matcher.
/// Sewing's displacement cap and reconstruction's numerical tolerance are
/// independent caller controls. Unsupported carriers retain their typed refusal.
/// Open, ambiguous or disconnected candidates never acquire solid authority;
/// the borrowed input faces remain available unchanged on every failure.
pub fn sew_constructed_surfaces(
    name: &str,
    faces: &[crate::Face],
    sewing: ToleranceOptions,
    options: n::NumericalOptions,
    dec: &mut Decider,
) -> Result<SewnSurfaceBody, ImportError> {
    sew_with_reconstruction(name, faces, sewing, options, |_, _| Ok(()), dec)
}

/// Run bounded geometric reconstruction after seam binding, before admission.
pub(crate) fn sew_with_reconstruction(
    name: &str,
    faces: &[crate::Face],
    sewing: ToleranceOptions,
    options: n::NumericalOptions,
    reconstruct: impl FnOnce(&mut [crate::Face], &mut Decider) -> Result<(), ImportError>,
    dec: &mut Decider,
) -> Result<SewnSurfaceBody, ImportError> {
    if !operation_budget_active() {
        return with_operation_budget(n::DEFAULT_WORK, || {
            sew_with_reconstruction(name, faces, sewing, options, reconstruct, dec)
        })
        .0;
    }
    let healed = crate::tolerant::healing::constructed::sew(name, faces, sewing, dec)?;
    dec.charge_work(healed.report().issues.len() as u64)?;
    if healed.bodies().len() != 1
        || healed.bodies()[0].shells().len() != 1
        || healed.report().issues.iter().any(|issue| {
            !issue.code.ends_with("-pending")
                && issue.code != "solid.unresolved-prerequisites"
                && issue.code != "wire.small-line-requires-curved-or-empty-boundary-analysis"
        })
    {
        return Err(Refusal::new(
            "surface.sew-not-solid",
            "constructed faces did not sew into one connected closed candidate",
        )
        .into());
    }
    let mut file = healed.definition(dec)?;
    let mut model = file.surface_models.pop().ok_or_else(|| {
        Refusal::new(
            "surface.sew-kind",
            "constructed input must retain its surface model",
        )
    })?;
    let shell = model
        .shells
        .pop()
        .ok_or_else(|| Refusal::new("surface.sew-kind", "constructed input lost its shell"))?;
    file.solids.push(Solid {
        name: model.name,
        faces: shell.faces,
    });
    reconstruct(&mut file.solids[0].faces, dec)?;
    let source_checks = intersection_checks(&file, OPTIONAL_EXACT_WORK, dec)?;
    let mut topology_repairs = Vec::new();
    let body = n::build(
        definition(&file.solids[0], &mut topology_repairs, dec)?,
        options,
        dec,
    )?;
    dec.check_deadline()?;
    let faces = file
        .solids
        .pop()
        .expect("one constructed solid")
        .faces
        .into();
    Ok(SewnSurfaceBody {
        body,
        faces,
        healing: healed,
        pair_checks: source_checks.pairs,
        topology_repairs,
    })
}
