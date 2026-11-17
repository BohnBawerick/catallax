//! Per-body outcomes. Solid refusal retains readable source surfaces and trims.
//! A retained surface representation supplies no exact-model or solid authority.

use catallax_numbers::{
    Decider, current_operation_work, operation_budget_active, with_operation_budget,
};
use catallax_recipe::StepId;
use catallax_topology::{BodyId, BodyKind, Model, numerical as n};
use serde::{Deserialize, Serialize};
use std::sync::Arc;

use super::{Error, ToleranceOptions, TolerantBody, TolerantFile, admission, healing};
use crate::{
    ImportRefusal,
    numerical::{self, StepImport},
};

/// Material/topology result, independent of exact versus numerical evidence.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Tier {
    /// A solid admitted without a recorded source interpretation/repair.
    ValidSolid,
    /// An admitted solid with a recorded derivative or numerical interpretation.
    RepairedSolid,
    /// Closed native topology with sampled fits and measured local tolerances.
    /// This is not a certified construction or an exact-model conversion.
    TolerantSolid,
    /// Retained surfaces; closed material is not asserted.
    SurfaceBody,
}

#[derive(Debug)]
enum Storage {
    Source,
    File(usize),
    Certified { model: Box<Model>, body: BodyId },
    Numerical(Box<StepImport>),
    Tolerant(usize),
}

/// One source body; all identities remain scoped to the file's immutable hash.
#[derive(Debug)]
pub struct BodyOutcome {
    tier: Tier,
    approximate: bool,
    storage: Storage,
    healed: Option<Box<healing::HealedFile>>,
    refusal: Option<Error>,
    checked_refusal: Option<Error>,
    healing_failure: Option<Box<healing::HealingFailure>>,
    upgrade_refusal: Option<Error>,
    first_blocker: Option<admission::NativeStage>,
}
impl BodyOutcome {
    /// Actual outcome, never inferred from declared source closure.
    pub fn tier(&self) -> Tier {
        self.tier
    }
    /// Numerical/repair or unverified surface representation, never exact.
    pub fn approximate(&self) -> bool {
        self.approximate
    }
    /// Actual remaining construction/admission refusal.
    pub fn refusal(&self) -> Option<&Error> {
        self.refusal.as_ref()
    }
    /// Earlier checked-tier refusal retained after successful native fallback.
    pub fn checked_refusal(&self) -> Option<&Error> {
        self.checked_refusal.as_ref()
    }
    /// Successful candidate checks and complete repair history, when attempted.
    pub fn healing(&self) -> Option<&healing::HealedFile> {
        self.healed.as_deref()
    }
    /// Typed failure and stage work; original source remains available.
    pub fn healing_failure(&self) -> Option<&healing::HealingFailure> {
        self.healing_failure.as_deref()
    }
    /// Optional stricter proof failure; an already usable native body is retained.
    pub fn upgrade_refusal(&self) -> Option<&Error> {
        self.upgrade_refusal.as_ref()
    }
    /// First native repair/check that prevents this body from being usable.
    pub fn first_blocker(&self) -> Option<admission::NativeStage> {
        self.first_blocker
    }
}

/// Evidence-specific access. Retained surfaces do not expose a kernel model.
pub enum BodyView<'a> {
    /// Checked construction; source repair quality is in the BodyOutcome.
    Certified {
        /// Existing checked model; source repair quality remains separate.
        model: &'a Model,
        /// This source body's handle in that model.
        body: BodyId,
    },
    /// Separately checked approximate B-rep with its recorded evidence.
    Numerical(&'a n::NumericalBody),
    /// Native B-rep. Its file owns the source carriers, trims and fit evidence.
    Tolerant {
        /// Whole source document, including original occurrence placement.
        import: &'a admission::TolerantImport,
        /// Source body's index in that document.
        body: usize,
    },
    /// Readable finite source declarations with unresolved boundary checks.
    Surfaces {
        /// Original supported surfaces, trims and shared source identities.
        source: &'a TolerantBody,
        /// Published repair candidate, still without solid admission.
        candidate: Option<&'a TolerantBody>,
    },
}

/// All readable bodies, including siblings of a solid which could not close.
#[derive(Debug)]
pub struct OpenedStep {
    source: Arc<TolerantFile>,
    source_graph: Option<Arc<crate::geometry_source::ParsedStep>>,
    bodies: Vec<BodyOutcome>,
    admission: Option<StepImport>,
    native_admission: Option<Arc<admission::TolerantImport>>,
    first_refusal: Option<Error>,
    read_work_used: u64,
    file_admission_work_used: u64,
    native_admission_work_used: u64,
    checked_upgrade_work_used: u64,
}
impl OpenedStep {
    /// Immutable original body definitions, before any repair.
    pub fn source(&self) -> &TolerantFile {
        &self.source
    }
    /// Original geometry, units and occurrence graph. Default opening retains
    /// this even when native fitting fails; checked-only opening may omit it.
    pub fn parsed_source(&self) -> Option<&crate::geometry_source::ParsedStep> {
        self.source_graph.as_deref()
    }
    /// Source solids followed by surface models, in original reader order.
    pub fn bodies(&self) -> &[BodyOutcome] {
        &self.bodies
    }
    /// The whole-file fast path's actual first refusal, if it did not complete.
    pub fn first_refusal(&self) -> Option<&Error> {
        self.first_refusal.as_ref()
    }
    /// Work for complete source reading/retention.
    pub fn read_work_used(&self) -> u64 {
        self.read_work_used
    }
    /// Work spent on the unchanged whole-file admission path.
    pub fn file_admission_work_used(&self) -> u64 {
        self.file_admission_work_used
    }
    /// Work for body-local native construction, boundary repair and fit checks.
    pub fn native_admission_work_used(&self) -> u64 {
        self.native_admission_work_used
    }
    /// Optional checked proofs of already-built native bodies.
    pub fn checked_upgrade_work_used(&self) -> u64 {
        self.checked_upgrade_work_used
    }
    /// Native evidence for every attempted body, including body-local refusals.
    pub fn tolerant_admission(&self) -> Option<&admission::TolerantImport> {
        self.native_admission.as_deref()
    }
    /// Existing source interpretation/repair evidence for an admitted body.
    pub fn admission(&self, index: usize) -> Option<&StepImport> {
        match &self.bodies.get(index)?.storage {
            Storage::File(_) => self.admission.as_ref(),
            Storage::Numerical(admitted) => Some(admitted),
            _ => None,
        }
    }
    /// Save immutable source identity and the existing versioned controls.
    /// Reopen reruns admission/healing; no healed-body codec is introduced.
    pub fn saved_controls(&self) -> super::ToleranceRecord {
        self.source.record()
    }
    /// Access one body's actual representation without an implicit promotion.
    pub fn body(&self, index: usize) -> Option<BodyView<'_>> {
        let outcome = self.bodies.get(index)?;
        Some(match &outcome.storage {
            Storage::File(body) => admission_view(self.admission.as_ref()?, *body),
            Storage::Certified { model, body } => BodyView::Certified { model, body: *body },
            Storage::Numerical(admitted) => admission_view(admitted, 0),
            Storage::Tolerant(body) => BodyView::Tolerant {
                import: self.native_admission.as_deref()?,
                body: *body,
            },
            Storage::Source => BodyView::Surfaces {
                source: &self.source.bodies[index],
                candidate: outcome.healed.as_ref().map(|h| &h.bodies()[0]),
            },
        })
    }
}

fn admission_view(admission: &StepImport, index: usize) -> BodyView<'_> {
    match admission {
        StepImport::Certified(a) => BodyView::Certified {
            model: a.model(),
            body: a.bodies()[index],
        },
        StepImport::Numerical(a) => BodyView::Numerical(&a.bodies()[index].1),
    }
}

fn classify(admission: &StepImport, index: usize) -> (Tier, bool) {
    match admission {
        StepImport::Certified(a) => {
            let body = a.model().body(a.bodies()[index]);
            let repaired = !body.marks.is_empty();
            (
                if body.kind == BodyKind::Surface {
                    Tier::SurfaceBody
                } else if repaired {
                    Tier::RepairedSolid
                } else {
                    Tier::ValidSolid
                },
                repaired,
            )
        }
        // Finite numerical source interpretation is explicitly recorded by
        // NumericalImport even when retained 3D geometry has not moved.
        StepImport::Numerical(_) => (Tier::RepairedSolid, true),
    }
}

fn control_refusal(error: &Error) -> bool {
    current_operation_work().exhausted
        || matches!(error, Error::Numerical(r) if r.code == "numbers.deadline")
}

fn retain_surface(outcome: &mut BodyOutcome, error: Error) {
    outcome.storage = Storage::Source;
    outcome.tier = Tier::SurfaceBody;
    outcome.approximate = true;
    outcome.refusal = Some(error);
}

// A private strict candidate must finish before the caller's deadline on both
// success and refusal. Expiry cannot publish it or start optional healing.
fn strict_retry(
    outcome: &mut BodyOutcome,
    definition: &crate::StepFile,
    step: &StepId,
    dec: &mut Decider,
) -> Result<Option<ImportRefusal>, catallax_numbers::Refusal> {
    let mut model = Model::new();
    let built =
        crate::build::build_file_with_deadline(&mut model, definition, Some(step), dec.deadline());
    dec.check_deadline()?;
    match built {
        Ok(ids) if ids.len() == 1 && !current_operation_work().exhausted => {
            let body = ids[0];
            outcome.approximate = !model.body(body).marks.is_empty();
            outcome.tier = if model.body(body).kind == BodyKind::Surface {
                Tier::SurfaceBody
            } else if outcome.approximate {
                Tier::RepairedSolid
            } else {
                Tier::ValidSolid
            };
            outcome.storage = Storage::Certified {
                model: Box::new(model),
                body,
            };
            Ok(None)
        }
        Ok(_) => Ok(Some(ImportRefusal::ResourceLimit)),
        Err(error) => Ok(Some(error)),
    }
}

/// Default opening policy: preserve the existing main admission first, then
/// build and heal only its refused or incomplete bodies from the shared source.
/// A failed candidate cannot discard already usable siblings.
///
/// Exhausted optional work retains the prepaid source representation.
/// No implicit work cap is installed for native admission.
pub fn open_step(
    bytes: &[u8],
    expected_sha256: &str,
    step: &StepId,
    options: ToleranceOptions,
    dec: &mut Decider,
) -> Result<OpenedStep, Error> {
    let before = current_operation_work().used;
    let (parsed, source, angles) = admission::prepare(bytes, expected_sha256, options, dec)?;
    let read_work_used = current_operation_work().used.saturating_sub(before);
    open_prepared(bytes, step, parsed, source, angles, read_work_used, dec)
}

fn open_prepared(
    bytes: &[u8],
    step: &StepId,
    parsed: Arc<crate::geometry_source::ParsedStep>,
    source: Arc<TolerantFile>,
    angles: crate::geometry_source::AngleUnit,
    read_work_used: u64,
    dec: &mut Decider,
) -> Result<OpenedStep, Error> {
    dec.charge_work(source.bodies.len() as u64 + 1)?;
    let mut result = OpenedStep {
        bodies: source
            .bodies
            .iter()
            .map(|_| BodyOutcome {
                tier: Tier::SurfaceBody,
                approximate: true,
                storage: Storage::Source,
                healed: None,
                refusal: None,
                checked_refusal: None,
                healing_failure: None,
                upgrade_refusal: None,
                first_blocker: None,
            })
            .collect(),
        source,
        source_graph: Some(Arc::clone(&parsed)),
        admission: None,
        native_admission: None,
        first_refusal: None,
        read_work_used,
        file_admission_work_used: 0,
        native_admission_work_used: 0,
        checked_upgrade_work_used: 0,
    };
    let before = current_operation_work().used;
    if operation_budget_active() {
        checked_admit(bytes, step, &mut result, dec);
    } else {
        let reservation = catallax_numbers::OperationAllowance::reserve(n::DEFAULT_WORK)?;
        reservation.run(|| checked_admit(bytes, step, &mut result, dec));
    }
    result.checked_upgrade_work_used = current_operation_work().used.saturating_sub(before);
    let incomplete: Vec<_> = result
        .bodies
        .iter()
        .map(|body| body.tier == Tier::SurfaceBody)
        .collect();
    if incomplete.iter().all(|needed| !needed) {
        return Ok(result);
    }
    let before = current_operation_work().used;
    let native = Arc::new(admission::admit(
        parsed,
        Arc::clone(&result.source),
        angles,
        admission::AdmissionPolicy::MainFirst(&incomplete),
        dec,
    )?);
    result.native_admission_work_used = current_operation_work().used.saturating_sub(before);
    for (source_index, failure) in native.body_failures() {
        result.bodies[source_index].refusal = Some(failure.error.clone());
        result.bodies[source_index].first_blocker = Some(failure.stage);
        if result.first_refusal.is_none() {
            result.first_refusal = Some(failure.error.clone());
        }
    }
    for (index, body) in native.bodies().iter().enumerate() {
        let outcome = &mut result.bodies[body.source_index()];
        outcome.checked_refusal = outcome.refusal.take();
        outcome.storage = Storage::Tolerant(index);
        outcome.tier = Tier::TolerantSolid;
        outcome.approximate = true;
    }
    result.native_admission = Some(native);
    if !current_operation_work().exhausted {
        let before = current_operation_work().used;
        if operation_budget_active() {
            upgrade_native(&mut result, step, dec);
        } else {
            let reservation = catallax_numbers::OperationAllowance::reserve(n::DEFAULT_WORK)?;
            reservation.run(|| upgrade_native(&mut result, step, dec));
        }
        result.checked_upgrade_work_used += current_operation_work().used.saturating_sub(before);
    }
    Ok(result)
}

fn upgrade_native(result: &mut OpenedStep, step: &StepId, dec: &mut Decider) {
    let native = result.native_admission.as_ref().expect("native outcomes");
    let options = n::NumericalOptions {
        tolerance_mm: result.source.options.maximum_mm.clone(),
    };
    for (index, body) in native.bodies().iter().enumerate() {
        let outcome = &mut result.bodies[body.source_index()];
        let attempt = (|| {
            let definition = native.checked_definition(index, dec)?;
            let strict_refusal = match strict_retry(outcome, &definition, step, dec)? {
                None => {
                    if !result.source.bodies[body.source_index()].declared_solid
                        || native.evidence().repairs.iter().any(|repair| {
                            repair.body() == body.name()
                                && !matches!(repair, healing::HealingChange::TraceRebuilt { .. })
                        })
                    {
                        outcome.tier = Tier::RepairedSolid;
                        outcome.approximate = true;
                    }
                    return Ok(());
                }
                Some(error) => error,
            };
            if current_operation_work().exhausted {
                return Err(ImportRefusal::ResourceLimit.into());
            }
            let admitted = numerical::admit_retained_definition(
                &definition,
                result.source.source_sha256(),
                strict_refusal,
                &options,
                dec,
            )?;
            outcome.tier = Tier::RepairedSolid;
            outcome.approximate = true;
            outcome.storage = Storage::Numerical(Box::new(admitted));
            Ok::<_, Error>(())
        })();
        if let Err(error) = attempt {
            let stop = control_refusal(&error);
            outcome.upgrade_refusal = Some(error);
            if stop {
                break;
            }
        }
    }
    // A refused stitch must not take away a previously usable source surface.
    // This optional proof checks the original surface, never the refused candidate.
    for (index, source) in result.source.bodies.iter().enumerate() {
        if current_operation_work().exhausted {
            break;
        }
        if source.declared_solid || !matches!(result.bodies[index].storage, Storage::Source) {
            continue;
        }
        let attempt = native
            .checked_source_definition(index, dec)
            .and_then(|definition| {
                strict_retry(&mut result.bodies[index], &definition, step, dec)
                    .map_err(Error::from)?
                    .map_or(Ok(()), |error| Err(error.into()))
            });
        if let Err(error) = attempt {
            let stop = control_refusal(&error);
            result.bodies[index].upgrade_refusal = Some(error);
            if stop {
                break;
            }
        }
    }
}

/// Explicit checked policy. Global proof and fixed-tolerance admission failures
/// remain refusals here; they are not native-tolerance defaults.
pub fn open_checked_step(
    bytes: &[u8],
    expected_sha256: &str,
    step: &StepId,
    options: ToleranceOptions,
    dec: &mut Decider,
) -> Result<OpenedStep, Error> {
    if !operation_budget_active() {
        return with_operation_budget(n::DEFAULT_WORK, || {
            open_checked_step(bytes, expected_sha256, step, options, dec)
        })
        .0;
    }
    let before = current_operation_work().used;
    let source = super::read_step(bytes, expected_sha256, options, dec)?;
    let read_work_used = current_operation_work().used.saturating_sub(before);
    // Reserve output slots before any optional admission can exhaust the query.
    dec.charge_work(source.bodies.len() as u64 + 1)?;
    let bodies = source
        .bodies
        .iter()
        .map(|_| BodyOutcome {
            tier: Tier::SurfaceBody,
            approximate: true,
            storage: Storage::Source,
            healed: None,
            refusal: None,
            checked_refusal: None,
            healing_failure: None,
            upgrade_refusal: None,
            first_blocker: None,
        })
        .collect();
    let mut result = OpenedStep {
        source: Arc::new(source),
        source_graph: None,
        bodies,
        admission: None,
        native_admission: None,
        first_refusal: None,
        read_work_used,
        file_admission_work_used: 0,
        native_admission_work_used: 0,
        checked_upgrade_work_used: 0,
    };
    checked_admit(bytes, step, &mut result, dec);
    Ok(result)
}

fn checked_admit(bytes: &[u8], step: &StepId, result: &mut OpenedStep, dec: &mut Decider) {
    let expected_sha256 = result.source.source_sha256.as_str();
    let numerical_options = n::NumericalOptions {
        tolerance_mm: result.source.options.maximum_mm.clone(),
    };
    let before = current_operation_work().used;
    let admission = {
        let _probe = numerical::CheckedStageProbe::new("checked.file_import", "file");
        numerical::import_step(bytes, expected_sha256, step, &numerical_options, dec)
    };
    result.file_admission_work_used = current_operation_work().used.saturating_sub(before);
    match admission {
        Ok(admitted) => {
            let count = match &admitted {
                StepImport::Certified(a) => a.bodies().len(),
                StepImport::Numerical(a) => a.bodies().len(),
            };
            if count == result.bodies.len() && !current_operation_work().exhausted {
                for (index, body) in result.bodies.iter_mut().enumerate() {
                    (body.tier, body.approximate) = classify(&admitted, index);
                    body.storage = Storage::File(index);
                }
                result.admission = Some(admitted);
                match dec.check_deadline() {
                    Ok(()) => return,
                    Err(error) => {
                        result.admission = None;
                        for body in &mut result.bodies {
                            retain_surface(body, error.clone().into());
                        }
                        result.first_refusal = Some(error.into());
                        return;
                    }
                }
            }
            result.first_refusal = Some(
                catallax_numbers::Refusal::new(
                    "import.opened-body-inventory",
                    "file admission did not retain the complete source body inventory",
                )
                .into(),
            );
        }
        Err(error) => result.first_refusal = Some(error),
    }
    if result.first_refusal.as_ref().is_some_and(control_refusal) {
        return;
    }
    for index in 0..result.bodies.len() {
        if current_operation_work().exhausted {
            break;
        }
        let _body_probe =
            numerical::CheckedStageProbe::new("checked.body", &result.source.bodies[index].name);
        let outcome = &mut result.bodies[index];
        let single = match {
            let _probe = numerical::CheckedStageProbe::new("checked.body_decode", "body");
            result.source.one_body(index, dec)
        } {
            Ok(single) => single,
            Err(error) => {
                let stop = control_refusal(&error);
                outcome.refusal = Some(error);
                if stop {
                    break;
                }
                continue;
            }
        };
        let definition = match {
            let _probe = numerical::CheckedStageProbe::new("checked.body_definition", "body");
            single.definition(dec)
        } {
            Ok(definition) => definition,
            Err(error) => {
                let stop = control_refusal(&error);
                outcome.refusal = Some(error);
                if stop {
                    break;
                }
                continue;
            }
        };
        let strict_refusal = match {
            let _probe = numerical::CheckedStageProbe::new("checked.body_strict", "body");
            strict_retry(outcome, &definition, step, dec)
        } {
            Ok(None) => continue,
            Ok(Some(error)) => error,
            Err(error) => {
                retain_surface(outcome, error.into());
                break;
            }
        };
        outcome.refusal = Some(strict_refusal.clone().into());
        if current_operation_work().exhausted {
            break;
        }
        let healed = match {
            let _probe = numerical::CheckedStageProbe::new("checked.body_healing", "body");
            healing::attempt(&single, dec)
        } {
            Ok(healed) => Box::new(healed),
            Err(failure) => {
                let stop = control_refusal(&failure.error);
                outcome.healing_failure = Some(failure);
                if stop {
                    break;
                }
                continue;
            }
        };
        // Known bad/open boundaries remain surfaces. Pending geometric checks
        // may be resolved by the numerical builder and its source-pair checks.
        if let Err(error) = dec.charge_work(healed.report().issues.len() as u64) {
            outcome.refusal = Some(error.into());
            outcome.healed = Some(healed);
            continue;
        }
        let may_admit = healed.bodies()[0].shells.len() == 1
            && healed.report().issues.iter().all(|issue| {
                issue.code.ends_with("-pending")
                    || issue.code == "solid.unresolved-prerequisites"
                    || issue.code == "wire.small-line-requires-curved-or-empty-boundary-analysis"
            });
        if may_admit {
            let definition = {
                let _probe = numerical::CheckedStageProbe::new("checked.healed_definition", "body");
                healed.definition(dec)
            };
            match definition.and_then(|mut definition| {
                // Only an already connected, closed candidate reaches here.
                // Source surface declarations remain unchanged in `source`.
                if let Some(surface) = definition.surface_models.pop() {
                    let mut shells = surface.shells;
                    definition.solids.push(crate::Solid {
                        name: surface.name,
                        faces: shells.remove(0).faces,
                    });
                }
                numerical::admit_retained_definition(
                    &definition,
                    expected_sha256,
                    strict_refusal,
                    &numerical_options,
                    dec,
                )
            }) {
                Ok(admitted) => {
                    outcome.tier = Tier::RepairedSolid;
                    outcome.approximate = true;
                    outcome.storage = Storage::Numerical(Box::new(admitted));
                    outcome.refusal = None;
                }
                Err(error) => outcome.refusal = Some(error),
            }
        }
        outcome.healed = Some(healed);
        if let Err(error) = dec.check_deadline() {
            retain_surface(outcome, error.into());
            break;
        }
    }
}

/// Reverify the immutable source and saved controls, then replay the ordinary
/// selector in the same request allowance. Derived topology is recomputed.
pub fn reopen(
    bytes: &[u8],
    controls: &super::ToleranceRecord,
    step: &StepId,
    dec: &mut Decider,
) -> Result<OpenedStep, Error> {
    let before = current_operation_work().used;
    super::charge_record(controls, dec)?;
    if controls.format != "catallax-import-tolerances/1" {
        return Err(super::invalid("unsupported tolerance record format"));
    }
    if controls.bodies.iter().any(|body| !body.changes.is_empty()) {
        return Err(catallax_numbers::Refusal::new(
            "import.opened-controls",
            "saved source controls must precede healing",
        )
        .into());
    }
    let parse = |text: &str| {
        let value = catallax_numbers::Rational::parse(text)
            .map_err(|_| super::invalid("invalid rational millimetre value"))?;
        if value.to_string() != text {
            return Err(super::invalid("noncanonical rational millimetre value"));
        }
        Ok(value)
    };
    let options = ToleranceOptions {
        initial_mm: parse(&controls.initial_mm)?,
        maximum_mm: parse(&controls.maximum_mm)?,
    };
    let (parsed, source, angles) =
        admission::prepare(bytes, &controls.source_sha256, options, dec)?;
    if source.bodies.len() != controls.bodies.len()
        || source
            .bodies
            .iter()
            .zip(&controls.bodies)
            .any(|(body, saved)| body.name != saved.name)
        || !super::matches_record(&source, controls, dec)?
    {
        return Err(super::invalid(
            "saved source controls differ from the immutable source inventory",
        ));
    }
    let read_work_used = current_operation_work().used.saturating_sub(before);
    open_prepared(bytes, step, parsed, source, angles, read_work_used, dec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{TestDeadline, TestDeadlinePoint};

    #[test]
    fn strict_retry_expiry_retains_prepaid_surfaces_before_healing() {
        let definition =
            crate::read_step(include_bytes!("../../tests/fixtures/spline-roof.step")).unwrap();
        let step = StepId::new("source").unwrap();
        let source_slot = || BodyOutcome {
            tier: Tier::SurfaceBody,
            approximate: true,
            storage: Storage::Source,
            healed: None,
            refusal: None,
            checked_refusal: None,
            healing_failure: None,
            upgrade_refusal: None,
            first_blocker: None,
        };
        let (_, complete) = with_operation_budget(n::DEFAULT_WORK, || {
            strict_retry(
                &mut source_slot(),
                &definition,
                &step,
                &mut Decider::default(),
            )
        });
        let clock = TestDeadline::after(TestDeadlinePoint::Work, 2);
        let mut dec = Decider::default();
        dec.set_deadline(Some(clock.instant()));
        let mut outcome = source_slot();
        let (expired, work) = with_operation_budget(n::DEFAULT_WORK, || {
            strict_retry(&mut outcome, &definition, &step, &mut dec)
        });
        assert_eq!(expired.unwrap_err().code, "numbers.deadline");
        clock.assert_expired();
        // A fresh deadline-free builder would do the entire strict replay and
        // discover expiry only at the outer publication check.
        assert!(work.used < complete.used, "{work:?} versus {complete:?}");
        assert!(!work.exhausted);
        assert!(matches!(outcome.storage, Storage::Source));
        assert_eq!(outcome.tier(), Tier::SurfaceBody);
        assert!(outcome.healed.is_none());
        assert!(outcome.healing_failure.is_none());
    }
}
