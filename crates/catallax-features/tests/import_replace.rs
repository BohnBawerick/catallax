//! Hash-bound file replacement controls use real production imports and recipes.

#[path = "../../catallax-import/tests/common/mod.rs"]
mod common;

use catallax_features::import_replace::{
    FileReplacementOutcome, ReferenceProblem, StepFileReplacement, replace_step_file,
};
use catallax_features::{
    ImportFiles, evaluate_with_files,
    import_edit::{StepReplacement, replace_step},
};
use catallax_import::{read_step, sha256_hex};
use catallax_numbers::Decider;
use catallax_recipe::{FileRef, Recipe, StepId};
use common::{BoxFixture, boxes_step, r};
use serde_json::json;
use std::collections::BTreeMap;

fn sources() -> (String, String) {
    let mut second = BoxFixture::new(["1", "1", "1"]);
    second.origin[0] = r("5");
    let old = boxes_step(&[BoxFixture::new(["2", "3", "4"]), second]);
    let new = renumber(&old);
    (old, new)
}
fn renumber(source: &str) -> String {
    let mut new = String::new();
    let mut chars = source.chars().peekable();
    while let Some(c) = chars.next() {
        new.push(c);
        if c == '#' {
            let mut digits = String::new();
            while chars.peek().is_some_and(char::is_ascii_digit) {
                digits.push(chars.next().unwrap());
            }
            new.push_str(&(digits.parse::<u64>().unwrap() + 10000).to_string());
        }
    }
    new
}
fn recipe(bytes: &[u8]) -> Recipe {
    let source = read_step(bytes).unwrap();
    Recipe::parse(&json!({"format":"catallax-recipe/5", "name":"two selected imports", "steps":[
        {"id":"first", "op":"import_step", "file":{"name":"parts.step","sha256":sha256_hex(bytes)}, "source_solid":source.solids[0].name,
            "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
        {"id":"second", "op":"import_step", "file":{"name":"parts.step","sha256":sha256_hex(bytes)}, "source_solid":source.solids[1].name,
            "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}
    ],"result":"second"}).to_string()).unwrap()
}

#[test]
fn legacy_selected_solid_replacement_reproduces_the_missing_file_transaction() {
    let (old, new) = sources();
    let recipe = recipe(old.as_bytes());
    let files = ImportFiles::from([
        (sha256_hex(old.as_bytes()), old.as_bytes()),
        (sha256_hex(new.as_bytes()), new.as_bytes()),
    ]);
    evaluate_with_files(&recipe, &files, &mut Decider::default()).unwrap();
    let keys = |bytes| {
        read_step(bytes)
            .unwrap()
            .solids
            .into_iter()
            .flat_map(|s| s.faces.into_iter().map(|f| f.name))
            .collect::<Vec<_>>()
    };
    let edit = StepReplacement {
        step: StepId::new("first").unwrap(),
        old_sha256: sha256_hex(old.as_bytes()),
        file: FileRef {
            name: "replacement.step".into(),
            sha256: sha256_hex(new.as_bytes()),
        },
        decisions: BTreeMap::new(),
        faces: keys(old.as_bytes())
            .into_iter()
            .zip(keys(new.as_bytes()).into_iter().map(|key| vec![key]))
            .collect(),
    };
    let previous = recipe.to_text();
    assert_eq!(
        replace_step(&recipe, &edit, &files, &mut Decider::default())
            .unwrap_err()
            .code,
        "import.replace-source-solid"
    );
    assert_eq!(recipe.to_text(), previous);
}

fn edit(old: &[u8], new: &[u8]) -> StepFileReplacement {
    StepFileReplacement {
        old_sha256: sha256_hex(old),
        file: FileRef {
            name: "replacement.step".into(),
            sha256: sha256_hex(new),
        },
        decisions: BTreeMap::new(),
    }
}

#[test]
fn file_replacement_updates_both_selected_imports_and_survives_reopen() {
    let (old, new) = sources();
    let recipe = recipe(old.as_bytes());
    let saved = recipe.to_text();
    let edit = edit(old.as_bytes(), new.as_bytes());
    let files = ImportFiles::from([
        (edit.old_sha256.clone(), old.as_bytes()),
        (edit.file.sha256.clone(), new.as_bytes()),
    ]);
    let FileReplacementOutcome::Ready(result) =
        replace_step_file(&recipe, &edit, &files, &mut Decider::default()).unwrap()
    else {
        panic!("unchanged definitions must replace")
    };
    assert_eq!(result.analysis.imports.len(), 2);
    assert!(result.analysis.unresolved_references.is_empty());
    assert_eq!(
        result.analysis.old_recipe_sha256,
        sha256_hex(saved.as_bytes())
    );
    let new_source = read_step(new.as_bytes()).unwrap();
    for (index, step) in result.built.recipe().steps.iter().enumerate() {
        let catallax_recipe::Operation::ImportStep {
            file, source_solid, ..
        } = &step.op
        else {
            unreachable!()
        };
        assert_eq!(file.sha256, edit.file.sha256);
        assert_eq!(
            source_solid.as_deref(),
            Some(new_source.solids[index].name.as_str())
        );
    }
    let reopened = Recipe::parse(&result.built.recipe().to_text()).unwrap();
    let again = evaluate_with_files(&reopened, &files, &mut Decider::default()).unwrap();
    assert_eq!(again.model.faces_of_body(again.result).len(), 6);
    assert_eq!(recipe.to_text(), saved);
    assert_eq!(
        evaluate_with_files(&recipe, &files, &mut Decider::default())
            .unwrap()
            .recipe(),
        &recipe
    );
}

#[test]
fn changed_geometry_behind_the_same_numbers_requires_resolution_for_the_whole_file() {
    let (old, _) = sources();
    let mut second = BoxFixture::new(["1.001", "1", "1"]);
    second.origin[0] = r("5");
    let new = boxes_step(&[BoxFixture::new(["2", "3", "4"]), second]);
    let recipe = recipe(old.as_bytes());
    let saved = recipe.to_text();
    let edit = edit(old.as_bytes(), new.as_bytes());
    let files = ImportFiles::from([
        (edit.old_sha256.clone(), old.as_bytes()),
        (edit.file.sha256.clone(), new.as_bytes()),
    ]);
    let FileReplacementOutcome::NeedsResolution(report) =
        replace_step_file(&recipe, &edit, &files, &mut Decider::default()).unwrap()
    else {
        panic!("reused numbers cannot bind changed geometry")
    };
    assert_eq!(report.imports.len(), 2);
    assert!(
        report
            .unresolved_references
            .iter()
            .any(|r| r.step.as_str() == "second" && r.reason == ReferenceProblem::UnverifiedSolid)
    );
    assert_eq!(recipe.to_text(), saved);
}

#[test]
fn corrupt_bytes_and_exhaustion_return_no_candidate() {
    let (old, new) = sources();
    let recipe = recipe(old.as_bytes());
    let saved = recipe.to_text();
    let edit = edit(old.as_bytes(), new.as_bytes());
    let mut files = ImportFiles::from([
        (edit.old_sha256.clone(), old.as_bytes()),
        (edit.file.sha256.clone(), new.as_bytes()),
    ]);
    let (result, work) = catallax_numbers::with_operation_budget(1, || {
        replace_step_file(&recipe, &edit, &files, &mut Decider::default())
    });
    assert!(work.exhausted);
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    files.insert(edit.file.sha256.clone(), old.as_bytes());
    assert_eq!(
        replace_step_file(&recipe, &edit, &files, &mut Decider::default())
            .unwrap_err()
            .code,
        "import.replace-hash"
    );
    assert_eq!(recipe.to_text(), saved);
}

#[test]
fn reordered_solid_entities_follow_verified_membership_instead_of_source_order() {
    let (old, _) = sources();
    let mut second = BoxFixture::new(["1", "1", "1"]);
    second.origin[0] = r("5");
    let new = renumber(&boxes_step(&[second, BoxFixture::new(["2", "3", "4"])]));
    let recipe = recipe(old.as_bytes());
    let edit = edit(old.as_bytes(), new.as_bytes());
    let files = ImportFiles::from([
        (edit.old_sha256.clone(), old.as_bytes()),
        (edit.file.sha256.clone(), new.as_bytes()),
    ]);
    let FileReplacementOutcome::Ready(result) =
        replace_step_file(&recipe, &edit, &files, &mut Decider::default()).unwrap()
    else {
        panic!("solid order is not identity")
    };
    let source = read_step(new.as_bytes()).unwrap();
    assert_eq!(
        result.analysis.imports[0].new_source_solid.as_deref(),
        Some(source.solids[1].name.as_str())
    );
    assert_eq!(
        result.analysis.imports[1].new_source_solid.as_deref(),
        Some(source.solids[0].name.as_str())
    );
}

#[test]
fn nested_copy_edge_vertex_and_touching_references_are_all_repointed() {
    let (old, new) = sources();
    let source = read_step(old.as_bytes()).unwrap();
    let faces = &source.solids[0].faces;
    let origin = |index: usize| json!({"step":"mirrored","role":{"copy":{"instance":1,"of":{"step":"first","role":{"imported":{"key":faces[index].name}}}}}});
    let edge = json!({"between":[origin(0),origin(4)],"touching":[origin(2)]});
    let mut value = serde_json::to_value(recipe(old.as_bytes())).unwrap();
    value["description"] = json!(format!("literal source label {}", faces[0].name));
    value["steps"].as_array_mut().unwrap().extend([
        json!({"id":"mirrored","op":"mirror","body":"first","plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}),
        json!({"id":"round","op":"fillet","body":"mirrored","radius":"0.1","edges":[edge],"overflow":{"keep_neighbours":[
            {"edge":edge,"at":{"faces":[origin(0),origin(2),origin(4)]},"setback":"0.2"},
            {"edge":edge,"at":{"faces":[origin(0),origin(3),origin(4)]},"setback":"0.2"}]}})
    ]);
    value["result"] = json!("round");
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let saved = recipe.to_text();
    let edit = edit(old.as_bytes(), new.as_bytes());
    let files = ImportFiles::from([
        (edit.old_sha256.clone(), old.as_bytes()),
        (edit.file.sha256.clone(), new.as_bytes()),
    ]);
    evaluate_with_files(&recipe, &files, &mut Decider::default())
        .expect("old nested references resolve");
    let FileReplacementOutcome::Ready(result) =
        replace_step_file(&recipe, &edit, &files, &mut Decider::default()).unwrap()
    else {
        panic!("all nested references are verified")
    };
    let candidate = serde_json::to_value(result.built.recipe()).unwrap();
    fn collect(value: &serde_json::Value, keys: &mut Vec<String>) {
        match value {
            serde_json::Value::Object(object) => {
                if object.get("step").and_then(serde_json::Value::as_str) == Some("first") {
                    keys.push(
                        object["role"]["imported"]["key"]
                            .as_str()
                            .unwrap()
                            .to_owned(),
                    );
                }
                for value in object.values() {
                    collect(value, keys);
                }
            }
            serde_json::Value::Array(values) => {
                for value in values {
                    collect(value, keys);
                }
            }
            _ => {}
        }
    }
    let mut keys = Vec::new();
    collect(&candidate, &mut keys);
    assert!(keys.len() >= 15);
    let new_faces = read_step(new.as_bytes()).unwrap().solids[0]
        .faces
        .iter()
        .map(|f| f.name.clone())
        .collect::<Vec<_>>();
    assert!(keys.iter().all(|key| new_faces.contains(key)));
    assert_eq!(candidate["description"], value["description"]);
    assert_eq!(recipe.to_text(), saved);
}

#[test]
fn an_unverified_fixed_face_is_reported_even_without_solid_selection() {
    let old = common::box_step(BoxFixture::new(["2", "3", "4"]));
    let new = common::box_step(BoxFixture::new(["5", "3", "4"]));
    let face = read_step(old.as_bytes()).unwrap().solids[0].faces[1]
        .name
        .clone();
    let recipe=Recipe::parse(&json!({"format":"catallax-recipe/5","name":"fixed face", "steps":[
        {"id":"source","op":"import_step","file":{"name":"old.step","sha256":sha256_hex(old.as_bytes())},"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
        {"id":"hollow","op":"shell","body":"source","thickness":"0.1","open":[{"face":{"step":"source","role":{"imported":{"key":face}}}}]}
    ],"result":"hollow"}).to_string()).unwrap();
    let edit = edit(old.as_bytes(), new.as_bytes());
    let files = ImportFiles::from([
        (edit.old_sha256.clone(), old.as_bytes()),
        (edit.file.sha256.clone(), new.as_bytes()),
    ]);
    evaluate_with_files(&recipe, &files, &mut Decider::default()).unwrap();
    let FileReplacementOutcome::NeedsResolution(report) =
        replace_step_file(&recipe, &edit, &files, &mut Decider::default()).unwrap()
    else {
        panic!("changed fixed face needs a decision")
    };
    assert_eq!(report.unresolved_references.len(), 1);
    assert_eq!(
        report.unresolved_references[0].reason,
        ReferenceProblem::UnverifiedFace
    );
    assert_eq!(
        report.unresolved_references[0].key.as_deref(),
        Some(face.as_str())
    );
}

#[test]
fn fixed_face_and_its_adjacency_filters_follow_the_new_file() {
    let (old, new) = sources();
    let source = read_step(old.as_bytes()).unwrap();
    let origin = |index: usize| json!({"step":"first","role":{"imported":{"key":source.solids[0].faces[index].name}}});
    let mut value = serde_json::to_value(recipe(old.as_bytes())).unwrap();
    value["steps"].as_array_mut().unwrap().push(
        json!({"id":"hollow","op":"shell","body":"first","thickness":"0.1",
        "open":[{"face":origin(0),"touching":[origin(2)],"near":["0","1","1"]}]}),
    );
    value["result"] = json!("hollow");
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let edit = edit(old.as_bytes(), new.as_bytes());
    let files = ImportFiles::from([
        (edit.old_sha256.clone(), old.as_bytes()),
        (edit.file.sha256.clone(), new.as_bytes()),
    ]);
    evaluate_with_files(&recipe, &files, &mut Decider::default()).unwrap();
    let FileReplacementOutcome::Ready(result) =
        replace_step_file(&recipe, &edit, &files, &mut Decider::default()).unwrap()
    else {
        panic!("fixed face geometry is unchanged")
    };
    let candidate = serde_json::to_value(result.built.recipe()).unwrap();
    let faces = read_step(new.as_bytes()).unwrap().solids[0].faces.clone();
    assert_eq!(
        candidate["steps"][2]["open"][0]["face"]["role"]["imported"]["key"],
        faces[0].name
    );
    assert_eq!(
        candidate["steps"][2]["open"][0]["touching"][0]["role"]["imported"]["key"],
        faces[2].name
    );
    assert_eq!(
        candidate["steps"][2]["open"][0]["near"],
        value["steps"][2]["open"][0]["near"]
    );
}

#[test]
fn duplicate_geometry_remains_an_explicit_solid_selection_issue() {
    let old = boxes_step(&[
        BoxFixture::new(["2", "3", "4"]),
        BoxFixture::new(["2", "3", "4"]),
    ]);
    let new = renumber(&old);
    let recipe = recipe(old.as_bytes());
    let edit = edit(old.as_bytes(), new.as_bytes());
    let files = ImportFiles::from([
        (edit.old_sha256.clone(), old.as_bytes()),
        (edit.file.sha256.clone(), new.as_bytes()),
    ]);
    evaluate_with_files(&recipe, &files, &mut Decider::default()).unwrap();
    let FileReplacementOutcome::NeedsResolution(report) =
        replace_step_file(&recipe, &edit, &files, &mut Decider::default()).unwrap()
    else {
        panic!("duplicate geometry must not bind by number")
    };
    assert_eq!(report.unresolved_references.len(), 2);
    assert!(
        report
            .unresolved_references
            .iter()
            .all(|r| r.reason == ReferenceProblem::UnverifiedSolid)
    );
}

#[test]
fn old_hash_decisions_are_not_adopted_for_new_bytes() {
    let (old, new) = sources();
    let recipe = recipe(old.as_bytes());
    let mut edit = edit(old.as_bytes(), new.as_bytes());
    edit.decisions.insert(
        StepId::new("first").unwrap(),
        BTreeMap::from([(
            format!("repair:{}:g0", &edit.old_sha256[..12]),
            "decline".into(),
        )]),
    );
    let files = ImportFiles::from([
        (edit.old_sha256.clone(), old.as_bytes()),
        (edit.file.sha256.clone(), new.as_bytes()),
    ]);
    assert_eq!(
        replace_step_file(&recipe, &edit, &files, &mut Decider::default())
            .unwrap_err()
            .code,
        "import.invalid-decision"
    );
}
