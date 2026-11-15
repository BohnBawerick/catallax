//! Complete-definition correspondence supplies mappings, never entity numbers.
//! A replacement rebuilds every affected import and downstream reference together.

use std::collections::{BTreeMap, BTreeSet};

use catallax_import::{
    SourceFormat, StepFile,
    correspondence::{EntityKind, ReplacementCorrespondence},
};
use catallax_numbers::{
    Decider, Refusal, current_operation_work, operation_budget_active, with_operation_budget,
};
use catallax_recipe::{FileRef, Operation, Recipe, StepId};
use catallax_topology::Model;
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::{BuildError, Built, ImportFiles, evaluate_with_files};

/// Replace all STEP import steps whose immutable source has the old digest.
#[derive(Clone, Debug)]
pub struct StepFileReplacement {
    /// Expected actual old source; stale source selections refuse.
    pub old_sha256: String,
    /// New immutable file, already present in `ImportFiles`.
    pub file: FileRef,
    /// New hash-bound import decisions, by affected step. Omitted steps use
    /// automatic policy only; old approvals are never copied to different bytes.
    pub decisions: BTreeMap<StepId, BTreeMap<String, String>>,
}

/// Why a current structured reference needs an explicit selection decision.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ReferenceProblem {
    /// Its finite face definition has no verified unique new partner.
    UnverifiedFace,
    /// Complete face membership does not verify one replacement source solid.
    UnverifiedSolid,
    /// A wildcard or non-imported role needs an explicit selection edit.
    NonSingletonOrigin,
}

/// One unresolved occurrence, including repeated uses of the same source face.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct UnresolvedReference {
    /// Import step that supplies the origin or solid membership.
    pub step: StepId,
    /// JSON pointer into the current canonical recipe.
    pub path: String,
    /// Old source key when the reference names one entity.
    pub key: Option<String>,
    /// Reason automatic replacement cannot choose a target.
    pub reason: ReferenceProblem,
}

/// Accepted source correspondence for one affected import interpretation.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ImportReplacementAnalysis {
    /// Origin identity retained by the import operation.
    pub step: StepId,
    /// Old selected source-solid key, if this step selects one member.
    pub old_source_solid: Option<String>,
    /// Verified new selected member, if one was established.
    pub new_source_solid: Option<String>,
    /// Producer-owned snapshots were compared before this report was created.
    pub correspondence: ReplacementCorrespondence,
}

/// Reports are information, never authority accepted back by the transaction.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FileReplacementAnalysis {
    /// Digest of the complete canonical old recipe, for an application's
    /// revision comparison before swapping its document to the returned build.
    pub old_recipe_sha256: String,
    /// Every import step that uses the replaced file.
    pub imports: Vec<ImportReplacementAnalysis>,
    /// All unresolved structured occurrences. No candidate model is returned
    /// until this list is empty and the complete recipe rebuild succeeds.
    pub unresolved_references: Vec<UnresolvedReference>,
}

/// A complete new revision; the caller retains the old recipe and bytes for undo.
#[derive(Debug)]
pub struct VerifiedFileReplacement {
    /// Rebuilt document, including the complete saved replacement recipe.
    pub built: Built,
    /// Correspondence and old revision identity used by this transaction.
    pub analysis: FileReplacementAnalysis,
}

/// Either a complete replacement or explicit choices, never a partial document.
#[derive(Debug)]
pub enum FileReplacementOutcome {
    /// All current references were verified and the whole recipe rebuilt.
    Ready(Box<VerifiedFileReplacement>),
    /// The prior document stays intact while the application requests choices.
    NeedsResolution(FileReplacementAnalysis),
}

fn error(kind: &str, detail: impl Into<String>) -> Refusal {
    let code = match kind {
        "recipe" => "import.replace-recipe",
        "missing-source" => "import.replace-missing-source",
        "hash" => "import.replace-hash",
        "file-not-found" => "import.replace-file-not-found",
        "decision-step" => "import.replace-decision-step",
        "reference" => "import.replace-reference",
        "build" => "import.replace-build",
        _ => unreachable!(),
    };
    Refusal::new(code, detail)
}

struct Binding {
    faces: BTreeMap<String, String>,
}

fn source_refusal(error_value: catallax_import::ImportRefusal, stage: &str) -> Refusal {
    let mut refusal = crate::evaluate::import_refusal(error_value);
    refusal.detail = format!("{stage}: {}", refusal.detail);
    refusal
}

fn rebuild_refusal(error_value: BuildError) -> Refusal {
    match error_value {
        BuildError::Geometry { step, mut refusal } => {
            refusal.detail = format!("replacement rebuild step {step}: {}", refusal.detail);
            refusal
        }
        error_value => error("build", error_value.to_string()),
    }
}

fn checked_result(
    result: Result<FileReplacementOutcome, Refusal>,
    exhausted: bool,
) -> Result<FileReplacementOutcome, Refusal> {
    match result {
        Ok(_) if exhausted => Err(Refusal::new(
            "numbers.operation-budget",
            "file replacement exhausted its operation allowance",
        )),
        // Keep the underlying code, stage context and undecided evidence even
        // when the failing lower layer also exhausted the shared allowance.
        result => result,
    }
}

/// Recompute correspondence from actual source bytes, update every matching
/// STEP import and every structured origin, and rebuild as one transaction.
///
/// Solid selection follows complete verified face membership, not entity numbers
/// or source order. Unresolved occurrences return reports without a new model.
/// On error the borrowed recipe, models and immutable source bytes are unchanged.
/// The caller's active allowance spans source imports, matching and full rebuild;
/// without one the complete transaction shares the unchanged 20M default.
pub fn replace_step_file(
    recipe: &Recipe,
    edit: &StepFileReplacement,
    files: &ImportFiles<'_>,
    dec: &mut Decider,
) -> Result<FileReplacementOutcome, Refusal> {
    if operation_budget_active() {
        let result = replace_inner(recipe, edit, files, dec);
        return checked_result(result, current_operation_work().exhausted);
    }
    let (result, work) =
        with_operation_budget(20_000_000, || replace_inner(recipe, edit, files, dec));
    checked_result(result, work.exhausted)
}

fn checked_bytes<'a>(
    hash: &str,
    files: &ImportFiles<'a>,
    dec: &mut Decider,
) -> Result<&'a [u8], Refusal> {
    let bytes = files
        .get(hash)
        .ok_or_else(|| error("missing-source", format!("missing bytes for sha256 {hash}")))?;
    dec.charge_work(bytes.len() as u64)?;
    if catallax_import::sha256_hex(bytes) != hash {
        return Err(error("hash", format!("bytes do not match sha256 {hash}")));
    }
    Ok(bytes)
}

fn selected_solid(
    old_key: &str,
    old: &StepFile,
    new: &StepFile,
    report: &ReplacementCorrespondence,
) -> Option<String> {
    let faces = old
        .solids
        .iter()
        .find(|s| s.name == old_key)?
        .faces
        .iter()
        .map(|f| f.name.as_str())
        .collect::<BTreeSet<_>>();
    if faces.is_empty() {
        return None;
    }
    let matches = report
        .verified
        .iter()
        .filter(|m| m.old.kind == EntityKind::Face && m.old.solid == old_key)
        .collect::<Vec<_>>();
    if matches.len() != faces.len() || !matches.iter().all(|m| faces.contains(m.old.key.as_str())) {
        return None;
    }
    let targets = matches
        .iter()
        .map(|m| m.new.solid.as_str())
        .collect::<BTreeSet<_>>();
    if targets.len() != 1 {
        return None;
    }
    let target = *targets.first()?;
    let expected = new
        .solids
        .iter()
        .find(|s| s.name == target)?
        .faces
        .iter()
        .map(|f| f.name.as_str())
        .collect::<BTreeSet<_>>();
    if matches
        .iter()
        .map(|m| m.new.key.as_str())
        .collect::<BTreeSet<_>>()
        != expected
    {
        return None;
    }
    Some(target.to_owned())
}

fn replace_inner(
    recipe: &Recipe,
    edit: &StepFileReplacement,
    files: &ImportFiles<'_>,
    dec: &mut Decider,
) -> Result<FileReplacementOutcome, Refusal> {
    dec.charge_work(1)?;
    recipe.check().map_err(|e| error("recipe", e.to_string()))?;
    let targets = recipe
        .steps
        .iter()
        .enumerate()
        .filter(|(_, step)| {
            matches!(&step.op, Operation::ImportStep { file, .. } if file.sha256 == edit.old_sha256)
        })
        .collect::<Vec<_>>();
    if targets.is_empty() {
        return Err(error(
            "file-not-found",
            "no STEP import uses the expected old source",
        ));
    }
    if edit
        .decisions
        .keys()
        .any(|id| !targets.iter().any(|(_, s)| &s.id == id))
    {
        return Err(error(
            "decision-step",
            "new decisions name an unaffected import step",
        ));
    }
    let old_bytes = checked_bytes(&edit.old_sha256, files, dec)?;
    let new_bytes = checked_bytes(&edit.file.sha256, files, dec)?;
    let old = SourceFormat::Step
        .read(old_bytes)
        .map_err(|e| source_refusal(e, "old source read"))?;
    let new = SourceFormat::Step
        .read(new_bytes)
        .map_err(|e| source_refusal(e, "new source read"))?;
    let text = recipe.to_text();
    dec.charge_work(text.len() as u64)?;
    let mut analysis = FileReplacementAnalysis {
        old_recipe_sha256: catallax_import::sha256_hex(text.as_bytes()),
        imports: Vec::new(),
        unresolved_references: Vec::new(),
    };
    let mut bindings = BTreeMap::new();
    let empty = BTreeMap::new();
    for (index, step) in &targets {
        dec.charge_work(1)?;
        let Operation::ImportStep {
            file,
            decisions,
            source_solid,
            ..
        } = &step.op
        else {
            unreachable!()
        };
        let old_choices = crate::import_edit::import_decisions(&file.sha256, decisions)?;
        let new_choices = crate::import_edit::import_decisions(
            &edit.file.sha256,
            edit.decisions.get(&step.id).unwrap_or(&empty),
        )?;
        let a = catallax_import::apply(
            &mut Model::new(),
            old_bytes,
            &file.name,
            &Default::default(),
            &old_choices,
            None,
        )
        .map_err(|e| source_refusal(e, &format!("old source import at step {}", step.id)))?;
        let b = catallax_import::apply(
            &mut Model::new(),
            new_bytes,
            &edit.file.name,
            &Default::default(),
            &new_choices,
            None,
        )
        .map_err(|e| source_refusal(e, &format!("new source import at step {}", step.id)))?;
        let correspondence = catallax_import::correspondence::compare(&a, &b, dec)?;
        let member = source_solid
            .as_deref()
            .and_then(|key| selected_solid(key, &old, &new, &correspondence));
        if source_solid.is_some() && member.is_none() {
            analysis.unresolved_references.push(UnresolvedReference {
                step: step.id.clone(),
                path: format!("/steps/{index}/source_solid"),
                key: source_solid.clone(),
                reason: ReferenceProblem::UnverifiedSolid,
            });
        }
        let mut faces = BTreeMap::new();
        let mut duplicate_keys = BTreeSet::new();
        for pair in correspondence
            .verified
            .iter()
            .filter(|m| m.old.kind == EntityKind::Face)
        {
            dec.charge_work(1)?;
            if faces
                .insert(pair.old.key.clone(), pair.new.key.clone())
                .is_some()
            {
                duplicate_keys.insert(pair.old.key.clone());
            }
        }
        // Origins contain a face key without a source-solid discriminator.
        // A repeated key cannot select one verified solid-specific relation.
        for key in duplicate_keys {
            faces.remove(&key);
        }
        bindings.insert(step.id.as_str().to_owned(), Binding { faces });
        analysis.imports.push(ImportReplacementAnalysis {
            step: step.id.clone(),
            old_source_solid: source_solid.clone(),
            new_source_solid: member,
            correspondence,
        });
    }
    let mut value = serde_json::to_value(recipe).expect("recipe serializes");
    rewrite(
        &mut value,
        "",
        &bindings,
        &mut analysis.unresolved_references,
        dec,
    )?;
    if !analysis.unresolved_references.is_empty() {
        return Ok(FileReplacementOutcome::NeedsResolution(analysis));
    }
    let mut candidate =
        Recipe::parse(&value.to_string()).map_err(|e| error("reference", e.to_string()))?;
    for result in &analysis.imports {
        let step = candidate
            .steps
            .iter_mut()
            .find(|s| s.id == result.step)
            .expect("checked step");
        let Operation::ImportStep {
            file,
            decisions,
            source_solid,
            ..
        } = &mut step.op
        else {
            unreachable!()
        };
        *file = edit.file.clone();
        *decisions = edit
            .decisions
            .get(&result.step)
            .cloned()
            .unwrap_or_default();
        *source_solid = result.new_source_solid.clone();
    }
    let built = evaluate_with_files(&candidate, files, dec).map_err(rebuild_refusal)?;
    Ok(FileReplacementOutcome::Ready(Box::new(
        VerifiedFileReplacement { built, analysis },
    )))
}

fn rewrite(
    value: &mut Value,
    path: &str,
    bindings: &BTreeMap<String, Binding>,
    unresolved: &mut Vec<UnresolvedReference>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    dec.charge_work(1)?;
    match value {
        Value::Object(object) => {
            if let Some(step) = object
                .get("step")
                .and_then(Value::as_str)
                .filter(|step| bindings.contains_key(*step))
                .map(str::to_owned)
            {
                let key = object
                    .get_mut("role")
                    .and_then(|r| r.get_mut("imported"))
                    .and_then(|r| r.get_mut("key"));
                if let Some(key) = key {
                    let old = key
                        .as_str()
                        .ok_or_else(|| error("reference", "invalid imported face key"))?;
                    if let Some(new) = bindings[&step].faces.get(old) {
                        *key = Value::String(new.clone());
                    } else {
                        unresolved.push(UnresolvedReference {
                            step: StepId::new(&step).expect("checked step"),
                            path: path.to_owned(),
                            key: Some(old.to_owned()),
                            reason: ReferenceProblem::UnverifiedFace,
                        });
                    }
                } else {
                    unresolved.push(UnresolvedReference {
                        step: StepId::new(&step).expect("checked step"),
                        path: path.to_owned(),
                        key: None,
                        reason: ReferenceProblem::NonSingletonOrigin,
                    });
                }
            }
            for (key, value) in object.iter_mut() {
                let key = key.replace('~', "~0").replace('/', "~1");
                rewrite(value, &format!("{path}/{key}"), bindings, unresolved, dec)?;
            }
            if let Some(Value::Array(pair)) = object.get("between")
                && pair.len() == 2
                && pair[0] == pair[1]
            {
                return Err(error(
                    "reference",
                    "replacement collapses an edge's two origins",
                ));
            }
        }
        Value::Array(values) => {
            for (index, value) in values.iter_mut().enumerate() {
                rewrite(value, &format!("{path}/{index}"), bindings, unresolved, dec)?;
            }
        }
        _ => {}
    }
    Ok(())
}

#[cfg(test)]
mod evidence_tests {
    use super::*;
    use catallax_numbers::{NotExact, Rational};

    fn decision() -> NotExact {
        NotExact {
            bound: Rational::from_ratio(1, 1 << 20).unwrap(),
            code: "numbers.precision-limit",
            detail: "order along a finite arc".into(),
        }
    }

    #[test]
    fn source_refusal_keeps_underlying_code_and_decision_with_stage() {
        let decision = decision();
        let refusal = Refusal::from(decision.clone());
        let error = catallax_import::ImportRefusal::SelfIntersectionNotVerified {
            solid: "#1".into(),
            faces: ["#2".into(), "#3".into()],
            reason: refusal.to_string(),
            refusal: Some(Box::new(refusal)),
        };
        let result = source_refusal(error, "old source import at step first");
        assert_eq!(result.code, "import.self-intersection-not-verified");
        assert_eq!(result.undecided(), Some(&decision));
        assert!(result.detail.contains("old source import at step first"));
    }

    #[test]
    fn rebuild_refusal_keeps_underlying_code_and_decision_with_step() {
        let decision = decision();
        let result = rebuild_refusal(BuildError::Geometry {
            step: StepId::new("downstream").unwrap(),
            refusal: Refusal::from(decision.clone()),
        });
        assert_eq!(result.code, "numbers.precision-limit");
        assert_eq!(result.undecided(), Some(&decision));
        assert!(
            result
                .detail
                .contains("replacement rebuild step downstream")
        );
    }

    #[test]
    fn exhausted_request_retains_refusal_but_discards_an_outcome() {
        let decision = decision();
        let refusal = checked_result(Err(Refusal::from(decision.clone())), true).unwrap_err();
        assert_eq!(refusal.code, "numbers.precision-limit");
        assert_eq!(refusal.undecided(), Some(&decision));
        let outcome = FileReplacementOutcome::NeedsResolution(FileReplacementAnalysis {
            old_recipe_sha256: String::new(),
            imports: Vec::new(),
            unresolved_references: Vec::new(),
        });
        assert_eq!(
            checked_result(Ok(outcome), true).unwrap_err().code,
            "numbers.operation-budget"
        );
    }
}
