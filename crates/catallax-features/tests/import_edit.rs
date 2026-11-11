//! Public, synthetic M24 fixtures. No external geometry oracle is used.
#[path = "../../catallax-import/tests/common/mod.rs"]
mod common;

use std::collections::BTreeMap;

use catallax_features::{
    ImportFiles, evaluate_with_files,
    import_edit::{StepReplacement, replace_step},
};
use catallax_import::{read_step, sha256_hex};
use catallax_measure::{Accuracy, volume};
use catallax_numbers::Decider;
use catallax_recipe::{FileRef, Recipe, StepId};
use common::*;
use serde_json::json;

fn recipe(bytes: &[u8]) -> Recipe {
    Recipe::parse(
        &json!({
            "format":"catallax-recipe/5", "name":"public import edit",
            "steps":[{"id":"source", "op":"import_step",
                "file":{"name":"boxes.step", "sha256":sha256_hex(bytes)},
                "pose":{"rotation":["1","0","0","0"],"translation":["7","-2","3"]}}],
            "result":"source"
        })
        .to_string(),
    )
    .unwrap()
}

fn renumber(text: &str) -> String {
    let mut output = String::new();
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        output.push(c);
        if c == '#' {
            let mut digits = String::new();
            while chars.peek().is_some_and(char::is_ascii_digit) {
                digits.push(chars.next().unwrap());
            }
            output.push_str(&(digits.parse::<u64>().unwrap() + 1000).to_string());
        }
    }
    output
}

fn proposal(old: &[u8], new: &[u8]) -> StepReplacement {
    let keys = |bytes| {
        read_step(bytes)
            .unwrap()
            .solids
            .into_iter()
            .flat_map(|s| s.faces.into_iter().map(|f| f.name))
            .collect::<Vec<_>>()
    };
    StepReplacement {
        step: StepId::new("source").unwrap(),
        old_sha256: sha256_hex(old),
        file: FileRef {
            name: "replacement.step".into(),
            sha256: sha256_hex(new),
        },
        decisions: BTreeMap::new(),
        faces: keys(old)
            .into_iter()
            .zip(keys(new).into_iter().map(|k| vec![k]))
            .collect(),
    }
}

fn shell_recipe(bytes: &[u8]) -> Recipe {
    let key = &read_step(bytes).unwrap().solids[0].faces[0].name;
    let mut value = serde_json::to_value(recipe(bytes)).unwrap();
    value["steps"].as_array_mut().unwrap().push(json!({
        "id":"hollow", "op":"shell", "body":"source", "thickness":"0.1",
        "open":[{"face":{"step":"source","role":{"imported":{"key":key}}}}]
    }));
    value["result"] = json!("hollow");
    Recipe::parse(&value.to_string()).unwrap()
}

#[test]
fn separated_solids_keep_lumps_names_volume_and_pose() {
    let mut second = BoxFixture::new(["2", "3", "4"]);
    second.origin = [r("5"), r("0"), r("0")];
    let text = boxes_step(&[BoxFixture::new(["1", "1", "1"]), second]);
    let files = ImportFiles::from([(sha256_hex(text.as_bytes()), text.as_bytes())]);
    let built =
        evaluate_with_files(&recipe(text.as_bytes()), &files, &mut Decider::default()).unwrap();
    assert_eq!(built.model.body(built.result).lumps.len(), 2);
    assert_eq!(built.model.faces_of_body(built.result).len(), 12);
    assert!(built.model.check(built.result).is_empty());
    let v = volume(
        &built.model,
        built.result,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert_eq!((v.lo, v.hi), (r("25"), r("25")));
    let bounds = catallax_measure::bounding_box(
        &built.model,
        built.result,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value
    .unwrap();
    assert_eq!(bounds.min, [r("7"), r("-2"), r("3")]);
    assert_eq!(bounds.max, [r("14"), r("1"), r("7")]);
    let mut rotated = serde_json::to_value(recipe(text.as_bytes())).unwrap();
    rotated["steps"][0]["pose"]["rotation"] = json!(["1", "0", "0", "1"]);
    let rotated = Recipe::parse(&rotated.to_string()).unwrap();
    let turned = evaluate_with_files(&rotated, &files, &mut Decider::default()).unwrap();
    let bounds = catallax_measure::bounding_box(
        &turned.model,
        turned.result,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value
    .unwrap();
    assert_eq!(bounds.min, [r("4"), r("-2"), r("3")]);
    assert_eq!(bounds.max, [r("7"), r("5"), r("7")]);
    let source = read_step(text.as_bytes()).unwrap();
    let membership = built
        .import_solids(&StepId::new("source").unwrap())
        .unwrap();
    assert_eq!(membership.len(), 2);
    for (entry, solid) in membership.iter().zip(&source.solids) {
        assert_eq!(entry.entity, solid.name);
        assert_eq!(entry.lumps.len(), 1);
    }
    for (lump, solid) in built
        .model
        .body(built.result)
        .lumps
        .iter()
        .zip(source.solids)
    {
        let faces = &built.model.shell(built.model.lump(*lump).shells[0]).faces;
        assert_eq!(faces.len(), solid.faces.len());
        for (face, source) in faces.iter().zip(solid.faces) {
            assert_eq!(
                built.model.face_origins(*face)[0].role,
                catallax_recipe::Role::Imported { key: source.name }
            );
        }
    }
}

#[test]
fn contact_overlap_and_duplicate_material_refuse() {
    for offset in ["0", "0.5", "1"] {
        let mut second = BoxFixture::new(["1", "1", "1"]);
        second.origin[0] = r(offset);
        let text = boxes_step(&[BoxFixture::new(["1", "1", "1"]), second]);
        let files = ImportFiles::from([(sha256_hex(text.as_bytes()), text.as_bytes())]);
        let error = evaluate_with_files(&recipe(text.as_bytes()), &files, &mut Decider::default())
            .unwrap_err();
        assert!(
            error.to_string().contains("import.multisolid-separation"),
            "{error}"
        );
    }
}

#[test]
fn explicit_replacement_rewrites_and_rebuilds_selection_without_mutating_old_revision() {
    let old = box_step(BoxFixture::new(["2", "3", "4"]));
    let new = renumber(&old);
    let recipe = shell_recipe(old.as_bytes());
    let previous = recipe.to_text();
    let edit = proposal(old.as_bytes(), new.as_bytes());
    let files = ImportFiles::from([
        (edit.old_sha256.clone(), old.as_bytes()),
        (edit.file.sha256.clone(), new.as_bytes()),
    ]);
    let built = replace_step(&recipe, &edit, &files, &mut Decider::default()).unwrap();
    assert_eq!(recipe.to_text(), previous);
    assert_eq!(built.recipe().steps[0].id, recipe.steps[0].id);
    assert_ne!(built.recipe(), &recipe);
    assert!(built.recipe().to_text().contains(&edit.file.sha256));
    let key = &read_step(new.as_bytes()).unwrap().solids[0].faces[0].name;
    assert!(built.recipe().to_text().contains(key));
    let old_built = evaluate_with_files(&recipe, &files, &mut Decider::default()).unwrap();
    let vol = |b: &catallax_features::Built| {
        volume(
            &b.model,
            b.result,
            &Accuracy::default(),
            &mut Decider::default(),
        )
        .unwrap()
        .value
    };
    assert_eq!(vol(&old_built), vol(&built));
}

#[test]
fn failures_are_typed_and_leave_recipe_and_bytes_unchanged() {
    let old = box_step(BoxFixture::new(["2", "3", "4"]));
    let new = renumber(&old);
    let recipe = shell_recipe(old.as_bytes());
    let saved = recipe.to_text();
    let base = proposal(old.as_bytes(), new.as_bytes());
    let files = ImportFiles::from([
        (base.old_sha256.clone(), old.as_bytes()),
        (base.file.sha256.clone(), new.as_bytes()),
    ]);
    let key = read_step(old.as_bytes()).unwrap().solids[0].faces[0]
        .name
        .clone();
    for (case, code) in [
        (0, "stale"),
        (1, "unmapped"),
        (2, "selection"),
        (3, "key"),
        (4, "build"),
        (5, "missing-source"),
    ] {
        let mut edit = base.clone();
        match case {
            0 => edit.old_sha256 = "0".repeat(64),
            1 => {
                edit.faces.remove(&key);
            }
            2 => {
                let target = edit.faces[&key][0].clone();
                edit.faces.get_mut(&key).unwrap().push(target);
            }
            3 => {
                edit.faces.insert(key.clone(), vec!["#999999".into()]);
            }
            4 => {
                edit.decisions
                    .insert("repair:stale:g0".into(), "decline".into());
            }
            _ => edit.file.sha256 = "1".repeat(64),
        }
        let error = replace_step(&recipe, &edit, &files, &mut Decider::default()).unwrap_err();
        assert_eq!(error.code, format!("import.replace-{code}"), "{error:?}");
        assert_eq!(recipe.to_text(), saved);
        assert_eq!(sha256_hex(old.as_bytes()), base.old_sha256);
        assert_eq!(sha256_hex(new.as_bytes()), base.file.sha256);
    }
}

#[test]
fn replacement_is_identical_at_one_two_and_four_workers() {
    let old = box_step(BoxFixture::new(["2", "3", "4"]));
    let new = renumber(&old);
    let recipe = shell_recipe(old.as_bytes());
    let edit = proposal(old.as_bytes(), new.as_bytes());
    let files = ImportFiles::from([
        (edit.old_sha256.clone(), old.as_bytes()),
        (edit.file.sha256.clone(), new.as_bytes()),
    ]);
    let run = || {
        let built = replace_step(&recipe, &edit, &files, &mut Decider::default()).unwrap();
        let mut step = Vec::new();
        catallax_step::write_step(
            &built.model,
            &[built.result],
            &Default::default(),
            &mut Decider::default(),
            &mut step,
        )
        .unwrap();
        (built.recipe().to_text(), step)
    };
    let expected = run();
    for workers in [1, 2, 4] {
        std::thread::scope(|scope| {
            let threads: Vec<_> = (0..workers).map(|_| scope.spawn(run)).collect();
            for thread in threads {
                assert_eq!(thread.join().unwrap(), expected);
            }
        });
    }
}

#[test]
fn incompatible_geometry_and_corrupt_bytes_roll_back() {
    let old = box_step(BoxFixture::new(["2", "3", "4"]));
    let new = renumber(&old);
    let old_faces = read_step(old.as_bytes()).unwrap().solids[0].faces.clone();
    let new_faces = read_step(new.as_bytes()).unwrap().solids[0].faces.clone();
    let origin = |key: &str| json!({"step":"source","role":{"imported":{"key":key}}});
    let mut value = serde_json::to_value(recipe(old.as_bytes())).unwrap();
    value["steps"].as_array_mut().unwrap().push(json!({
        "id":"round", "op":"fillet", "body":"source", "radius":"0.1",
        "edges":[{"between":[origin(&old_faces[0].name),origin(&old_faces[2].name)]}]
    }));
    value["result"] = json!("round");
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let saved = recipe.to_text();
    let mut edit = proposal(old.as_bytes(), new.as_bytes());
    edit.faces
        .insert(old_faces[2].name.clone(), vec![new_faces[1].name.clone()]);
    let mut files = ImportFiles::from([
        (edit.old_sha256.clone(), old.as_bytes()),
        (edit.file.sha256.clone(), new.as_bytes()),
    ]);
    evaluate_with_files(&recipe, &files, &mut Decider::default()).expect("old edge is valid");
    assert_eq!(
        replace_step(&recipe, &edit, &files, &mut Decider::default())
            .unwrap_err()
            .code,
        "import.replace-build"
    );
    files.insert(edit.file.sha256.clone(), old.as_bytes());
    assert_eq!(
        replace_step(&recipe, &edit, &files, &mut Decider::default())
            .unwrap_err()
            .code,
        "import.replace-hash"
    );
    assert_eq!(recipe.to_text(), saved);
}

#[test]
fn mapped_occurrences_are_not_silently_read_at_definition_locations() {
    let mut second = BoxFixture::new(["1", "1", "1"]);
    second.origin[0] = r("4");
    let text = boxes_step(&[BoxFixture::new(["1", "1", "1"]), second]);
    let text = text.replace(
        "ENDSEC;\nEND-ISO-10303-21;",
        "#90000=MAPPED_ITEM('',#1,#2);\nENDSEC;\nEND-ISO-10303-21;",
    );
    assert!(
        matches!(read_step(text.as_bytes()),Err(catallax_import::ImportRefusal::Unsupported {type_name,..}) if type_name.contains("MAPPED_ITEM"))
    );
}
