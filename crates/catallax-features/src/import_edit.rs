//! M24 import edits. Correspondence is supplied editorial intent, never inferred
//! from entity numbers or geometric proximity. No input recipe is mutated.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Decider, Refusal};
use catallax_recipe::{FileRef, Operation, Recipe, StepId};
use catallax_topology::{BodyId, Model};
use serde_json::Value;

use crate::{Built, ImportFiles, evaluate_with_files};

/// One source solid's membership after applying the recipe placement.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ImportedSolid {
    /// MANIFOLD_SOLID_BREP entity key in the immutable source file.
    pub entity: String,
    /// Placed topology lumps, in source order; never an inferred union.
    pub lumps: Vec<catallax_topology::LumpId>,
}

/// An explicit, hash-bound replacement proposal for one import step.
#[derive(Clone, Debug)]
pub struct StepReplacement {
    /// Import step whose origin identity is retained.
    pub step: StepId,
    /// Expected old source digest; stale edits fail.
    pub old_sha256: String,
    /// New immutable source. Its bytes must be present in `ImportFiles`.
    pub file: FileRef,
    /// Newly resolved decisions. Old decisions are never copied implicitly.
    pub decisions: BTreeMap<String, String>,
    /// Explicit relation between face entity keys in the two hash-fixed files.
    /// This bounded entry accepts singleton targets only. Splits and deletion
    /// need an explicit selection edit, not an arbitrary chosen target.
    pub faces: BTreeMap<String, Vec<String>>,
}

fn refusal(kind: &'static str, detail: impl Into<String>) -> Refusal {
    Refusal::new(kind, detail)
}

/// Validate and rebuild a replacement as one transaction. The returned build's
/// `recipe()` is the complete new saved revision. The caller can atomically swap
/// that revision into its document while retaining the old recipe and bytes for
/// undo. On any failure no document, model or source bytes have been changed.
///
/// All structured origins, including nested copies, offsets, edge/vertex names
/// and touching selectors, are visited through the recipe's own wire schema.
/// Re-parsing restores canonical name order and rejects collapsed names. A full
/// rebuild then validates downstream selections and compatibility. This does not
/// certify geometric equivalence of the editorial correspondence.
pub fn replace_step(
    recipe: &Recipe,
    edit: &StepReplacement,
    files: &ImportFiles<'_>,
    dec: &mut Decider,
) -> Result<Built, Refusal> {
    recipe
        .check()
        .map_err(|e| refusal("import.replace-recipe", e.to_string()))?;
    let step = recipe
        .steps
        .iter()
        .find(|s| s.id == edit.step)
        .ok_or_else(|| refusal("import.replace-step", "the import step does not exist"))?;
    let (file, format) = match &step.op {
        Operation::ImportStep { file, .. } => (file, catallax_import::SourceFormat::Step),
        Operation::ImportIges { file, .. } => (file, catallax_import::SourceFormat::Iges),
        _ => {
            return Err(refusal(
                "import.replace-step",
                "the selected step is not a STEP or IGES import",
            ));
        }
    };
    if file.sha256 != edit.old_sha256 {
        return Err(refusal(
            "import.replace-stale",
            "the old source hash no longer matches the recipe",
        ));
    }
    let keys = |hash: &str| -> Result<BTreeSet<String>, Refusal> {
        let bytes = files.get(hash).ok_or_else(|| {
            refusal(
                "import.replace-missing-source",
                format!("missing bytes for sha256 {hash}"),
            )
        })?;
        if catallax_import::sha256_hex(bytes) != hash {
            return Err(refusal(
                "import.replace-hash",
                format!("bytes do not match sha256 {hash}"),
            ));
        }
        let source = format
            .read(bytes)
            .map_err(|e| refusal("import.replace-source", format!("{e:?}")))?;
        Ok(source
            .solids
            .into_iter()
            .flat_map(|s| s.faces.into_iter().map(|f| f.name))
            .collect())
    };
    let old = keys(&edit.old_sha256)?;
    let new = keys(&edit.file.sha256)?;
    for (key, targets) in &edit.faces {
        if !old.contains(key) {
            return Err(refusal(
                "import.replace-key",
                format!("{key} is not a face in the old file"),
            ));
        }
        if targets.len() != 1 {
            return Err(refusal(
                "import.replace-selection",
                format!("{key} needs an explicit split or deletion selection edit"),
            ));
        }
        if !new.contains(&targets[0]) {
            return Err(refusal(
                "import.replace-key",
                format!("{} is not a face in the new file", targets[0]),
            ));
        }
    }
    let mut value = serde_json::to_value(recipe).expect("recipe serializes");
    rewrite(&mut value, edit)?;
    let mut candidate = Recipe::parse(&value.to_string())
        .map_err(|e| refusal("import.replace-reference", e.to_string()))?;
    let step = candidate
        .steps
        .iter_mut()
        .find(|s| s.id == edit.step)
        .expect("checked step");
    let (Operation::ImportStep {
        file, decisions, ..
    }
    | Operation::ImportIges {
        file, decisions, ..
    }) = &mut step.op
    else {
        unreachable!()
    };
    *file = edit.file.clone();
    *decisions = edit.decisions.clone();
    evaluate_with_files(&candidate, files, dec)
        .map_err(|e| refusal("import.replace-build", e.to_string()))
}

fn rewrite(value: &mut Value, edit: &StepReplacement) -> Result<(), Refusal> {
    match value {
        Value::Object(object) => {
            if object.get("step").and_then(Value::as_str) == Some(edit.step.as_str()) {
                let key = object
                    .get_mut("role")
                    .and_then(|r| r.get_mut("imported"))
                    .and_then(|r| r.get_mut("key"))
                    .ok_or_else(|| {
                        refusal(
                            "import.replace-selection",
                            "an import wildcard requires an explicit selection edit",
                        )
                    })?;
                let old = key.as_str().ok_or_else(|| {
                    refusal("import.replace-reference", "invalid imported face key")
                })?;
                let targets = edit.faces.get(old).ok_or_else(|| {
                    refusal(
                        "import.replace-unmapped",
                        format!("no confirmed replacement for {old}"),
                    )
                })?;
                *key = Value::String(targets[0].clone());
            }
            for v in object.values_mut() {
                rewrite(v, edit)?;
            }
            // EdgeRef is not an EdgeName: reject merged endpoints explicitly,
            // even when parsing a selector would otherwise permit them.
            if let Some(Value::Array(pair)) = object.get("between")
                && pair.len() == 2
                && pair[0] == pair[1]
            {
                return Err(refusal(
                    "import.replace-reference",
                    "replacement collapses an edge's two origins",
                ));
            }
        }
        Value::Array(values) => {
            for v in values {
                rewrite(v, edit)?;
            }
        }
        _ => {}
    }
    Ok(())
}

/// Strict enclosure separation proves that retaining all lumps neither merges
/// solid membership nor duplicates material. No union healing is performed.
pub(crate) fn collect_solids(
    model: &mut Model,
    bodies: &[(String, BodyId)],
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let Some((_, first)) = bodies.first() else {
        return Err(Refusal::new(
            "import.not-solid",
            "no closed solid was imported",
        ));
    };
    if bodies.len() == 1 {
        return Ok(*first);
    }
    let mut boxes: Vec<catallax_measure::BoundingBox> = Vec::new();
    for (index, (name, id)) in bodies.iter().enumerate() {
        let bounds = catallax_measure::bounding_box(model, *id, &Default::default(), dec)?
            .value
            .ok_or_else(|| Refusal::new("import.not-solid", format!("{name} has no bounds")))?;
        for (other, b) in boxes.iter().enumerate() {
            if !(0..3).any(|axis| bounds.max[axis] < b.min[axis] || b.max[axis] < bounds.min[axis])
            {
                return Err(Refusal::new(
                    "import.multisolid-separation",
                    format!(
                        "{} and {} need a certified separation or union; their enclosing boxes are not strictly separated",
                        bodies[other].0, bodies[index].0
                    ),
                ));
            }
        }
        boxes.push(bounds);
    }
    let mut result = model.body(*first).clone();
    for (_, id) in &bodies[1..] {
        result.lumps.extend_from_slice(&model.body(*id).lumps);
        result.marks.extend(&model.body(*id).marks);
    }
    Ok(model.add_body(result))
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn edit() -> StepReplacement {
        StepReplacement {
            step: StepId::new("import").unwrap(),
            old_sha256: "a".repeat(64),
            file: FileRef {
                name: "new.step".into(),
                sha256: "b".repeat(64),
            },
            decisions: BTreeMap::new(),
            faces: BTreeMap::from([
                ("#1".into(), vec!["#10".into()]),
                ("#2".into(), vec!["#20".into()]),
            ]),
        }
    }

    #[test]
    fn visits_nested_origins_touching_edges_vertices_and_preserves_other_imports() {
        let origin = |key| json!({"step":"import","role":{"imported":{"key":key}}});
        let a = origin("#1");
        let b = origin("#2");
        let mut references = json!({
            "face": {"step":"copy","role":{"copy":{"instance":1,"of":a}}},
            "touching":[b], "edge":{"between":[a,b]}, "vertex":{"faces":[a,b]},
            "inner":{"step":"shell","role":{"inner":{"of":a}}},
            "other":{"step":"other","role":{"imported":{"key":"#1"}}}
        });
        rewrite(&mut references, &edit()).unwrap();
        assert_eq!(references["face"]["role"]["copy"]["of"], origin("#10"));
        assert_eq!(references["touching"][0], origin("#20"));
        assert_eq!(
            references["edge"]["between"],
            json!([origin("#10"), origin("#20")])
        );
        assert_eq!(
            references["vertex"]["faces"],
            json!([origin("#10"), origin("#20")])
        );
        assert_eq!(references["inner"]["role"]["inner"]["of"], origin("#10"));
        assert_eq!(references["other"]["role"]["imported"]["key"], "#1");
    }

    #[test]
    fn wildcard_and_merged_edge_require_explicit_edits() {
        let mut wildcard = json!({"step":"import"});
        assert_eq!(
            rewrite(&mut wildcard, &edit()).unwrap_err().code,
            "import.replace-selection"
        );
        let mut edit = edit();
        edit.faces.insert("#2".into(), vec!["#10".into()]);
        let mut edge = json!({"between":[
            {"step":"import","role":{"imported":{"key":"#1"}}},
            {"step":"import","role":{"imported":{"key":"#2"}}}
        ]});
        assert_eq!(
            rewrite(&mut edge, &edit).unwrap_err().code,
            "import.replace-reference"
        );
    }
}
