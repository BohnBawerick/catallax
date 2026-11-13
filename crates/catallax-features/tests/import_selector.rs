//! Source entity selection retains the original file and checkpoint history.
#[path = "../../catallax-import/tests/common/mod.rs"]
mod common;

use catallax_features::{
    BuildError, Built, ImportFiles, evaluate_incremental_with_files, evaluate_with_files,
};
use catallax_import::{SourceFormat, sha256_hex};
use catallax_measure::{Accuracy, volume};
use catallax_numbers::Decider;
use catallax_recipe::{Recipe, Role, StepId};
use common::*;
use serde_json::json;

fn recipe(bytes: &[u8], format: SourceFormat, entity: Option<&str>) -> Recipe {
    let mut step = json!({
        "id":"source", "op":if format == SourceFormat::Step {"import_step"} else {"import_iges"},
        "file":{"name":"reviewed source", "sha256":sha256_hex(bytes)},
        "pose":{"rotation":["1","0","0","0"],"translation":["7","-2","3"]}
    });
    if let Some(entity) = entity {
        step["source_solid"] = json!(entity);
    }
    Recipe::parse(&json!({
        "format":"catallax-recipe/5", "name":"selected source", "steps":[step], "result":"source"
    }).to_string()).unwrap()
}

fn measured(built: &Built) -> catallax_numbers::Rational {
    let value = volume(
        &built.model,
        built.result,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert_eq!(value.lo, value.hi);
    value.lo
}

fn selected_solids(bytes: &[u8], format: SourceFormat, expected: &[&str]) {
    let original = bytes.to_vec();
    let entities = format.solid_entities(bytes).unwrap();
    assert_eq!(entities.len(), expected.len());
    assert_eq!(entities, format.solid_entities(bytes).unwrap());
    let source = format.read(bytes).unwrap();
    let files = ImportFiles::from([(sha256_hex(bytes), bytes)]);
    let step = StepId::new("source").unwrap();
    for ((entity, solid), expected) in entities.iter().zip(&source.solids).zip(expected) {
        let recipe = recipe(bytes, format, Some(entity));
        let checkpoint =
            evaluate_incremental_with_files(&recipe, None, &files, &mut Decider::default())
                .unwrap();
        let built = checkpoint.built();
        assert_eq!(measured(built), r(expected));
        assert_eq!(built.model.body(built.result).lumps.len(), 1);
        assert!(built.model.check(built.result).is_empty());
        let membership = built.import_solids(&step).unwrap();
        assert_eq!(membership.len(), 1);
        assert_eq!(&membership[0].entity, entity);
        assert_eq!(membership[0].lumps, built.model.body(built.result).lumps);
        let mut keys = built
            .model
            .faces_of_body(built.result)
            .iter()
            .map(|face| {
                let origins = built.model.face_origins(*face);
                assert_eq!(origins.len(), 1);
                assert_eq!(origins[0].step, step);
                let Role::Imported { key } = &origins[0].role else {
                    panic!("source origin lost")
                };
                key.clone()
            })
            .collect::<Vec<_>>();
        let mut expected_keys = solid
            .faces
            .iter()
            .map(|face| face.name.clone())
            .collect::<Vec<_>>();
        keys.sort();
        expected_keys.sort();
        assert_eq!(keys, expected_keys);
        assert_eq!(
            built.import_record(&step).unwrap().source.sha256,
            sha256_hex(bytes)
        );
        assert_eq!(built.recipe(), &recipe);
        // Reusing the exact prefix must retain its topology IDs and repair record.
        let reused = evaluate_incremental_with_files(
            &recipe,
            Some(checkpoint.clone()),
            &files,
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(reused.built().result, built.result);
        assert_eq!(
            reused.built().import_record(&step),
            built.import_record(&step)
        );
        assert_eq!(
            reused.built().import_solids(&step),
            built.import_solids(&step)
        );
        let mut suffix = serde_json::to_value(&recipe).unwrap();
        suffix["steps"].as_array_mut().unwrap().push(json!({
            "id":"reflected", "op":"mirror", "body":"source",
            "plane":{"point":["0","0","0"],"normal":["1","0","0"]},
            "keep_original":false
        }));
        suffix["result"] = json!("reflected");
        let suffix = Recipe::parse(&suffix.to_string()).unwrap();
        let mut resumed_dec = Decider::default();
        let resumed = evaluate_incremental_with_files(
            &suffix,
            Some(checkpoint.clone()),
            &files,
            &mut resumed_dec,
        )
        .unwrap();
        let mut fresh_dec = Decider::default();
        let fresh = evaluate_with_files(&suffix, &files, &mut fresh_dec).unwrap();
        assert_eq!(measured(resumed.built()), measured(&fresh));
        assert_eq!(
            resumed.built().import_solids(&step),
            built.import_solids(&step)
        );
        assert_eq!(
            resumed.built().import_record(&step),
            built.import_record(&step)
        );
        assert_eq!(resumed.built().bodies[&step], built.result);
        // A different selector invalidates the cached prefix.
        if let Some(other) = entities.iter().find(|key| *key != entity) {
            let other_recipe = crate::recipe(bytes, format, Some(other));
            let changed = evaluate_incremental_with_files(
                &other_recipe,
                Some(checkpoint.clone()),
                &files,
                &mut Decider::default(),
            )
            .unwrap();
            assert_eq!(
                &changed.built().import_solids(&step).unwrap()[0].entity,
                other
            );
        }
    }
    assert_eq!(bytes, original);
}

#[test]
fn selects_distinct_step_solids_even_with_reversed_entity_records() {
    let mut second = BoxFixture::new(["2", "3", "4"]);
    second.origin[0] = r("5");
    let text = boxes_step(&[BoxFixture::new(["1", "1", "1"]), second]);
    selected_solids(text.as_bytes(), SourceFormat::Step, &["1", "24"]);
    let (header, rest) = text.split_once("DATA;\n").unwrap();
    let (entities, footer) = rest.split_once("ENDSEC;").unwrap();
    let reversed = format!(
        "{header}DATA;\n{}\nENDSEC;{footer}",
        entities.lines().rev().collect::<Vec<_>>().join("\n")
    );
    assert_eq!(
        SourceFormat::Step.solid_entities(text.as_bytes()).unwrap(),
        SourceFormat::Step
            .solid_entities(reversed.as_bytes())
            .unwrap()
    );
    selected_solids(reversed.as_bytes(), SourceFormat::Step, &["1", "24"]);
}

#[test]
fn shared_shell_solids_have_distinct_identities_and_select_before_collection() {
    let text = box_step(BoxFixture::new(["2", "2", "3"]));
    let solid = text
        .lines()
        .find(|line| line.contains("=MANIFOLD_SOLID_BREP"))
        .unwrap();
    let (_, definition) = solid.split_once('=').unwrap();
    let text = text.replace(
        "ENDSEC;\nEND-ISO",
        &format!("#999999={definition}\nENDSEC;\nEND-ISO"),
    );
    let entities = SourceFormat::Step.solid_entities(text.as_bytes()).unwrap();
    assert_eq!(entities.len(), 2);
    assert_ne!(entities[0], entities[1]);
    let files = ImportFiles::from([(sha256_hex(text.as_bytes()), text.as_bytes())]);
    assert!(
        evaluate_with_files(
            &recipe(text.as_bytes(), SourceFormat::Step, None),
            &files,
            &mut Decider::default()
        )
        .is_err()
    );
    selected_solids(text.as_bytes(), SourceFormat::Step, &["12", "12"]);
}

#[test]
fn selects_iges_solids_by_directory_entity() {
    let mut second = BoxFixture::new(["3", "2", "2"]);
    second.origin[0] = r("10");
    let step = boxes_step(&[BoxFixture::new(["2", "2", "2"]), second]);
    let source = SourceFormat::Step.read(step.as_bytes()).unwrap();
    let mut model = catallax_topology::Model::new();
    let bodies = catallax_import::build_file(&mut model, &source, None).unwrap();
    let mut bytes = Vec::new();
    catallax_iges::write_iges(
        &model,
        &bodies,
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    assert!(
        SourceFormat::Iges
            .solid_entities(&bytes)
            .unwrap()
            .iter()
            .all(|key| key.starts_with('D'))
    );
    selected_solids(&bytes, SourceFormat::Iges, &["8", "12"]);
}

#[test]
fn omitted_selector_keeps_recipe_bytes_and_missing_entity_refuses() {
    let text = boxes_step(&[BoxFixture::new(["1", "1", "1"])]);
    for format in [SourceFormat::Step, SourceFormat::Iges] {
        let recipe = recipe(text.as_bytes(), format, None);
        let saved = recipe.to_text();
        assert!(!saved.contains("source_solid"));
        assert_eq!(Recipe::parse(&saved).unwrap().to_text(), saved);
    }
    let files = ImportFiles::from([(sha256_hex(text.as_bytes()), text.as_bytes())]);
    for entity in ["", "#999999", "D1", "#1"] {
        let error = evaluate_with_files(
            &recipe(text.as_bytes(), SourceFormat::Step, Some(entity)),
            &files,
            &mut Decider::default(),
        )
        .unwrap_err();
        let BuildError::Geometry { refusal, .. } = error else {
            panic!("expected geometry refusal")
        };
        assert_eq!(refusal.code, "import.source-solid");
    }
}

#[test]
fn selected_repair_keeps_original_record_decisions_and_hash_checks() {
    let mut first = BoxFixture::new(["10", "10", "10"]);
    first.shift = Some((1, [r("0.004"), r("0"), r("0")]));
    let mut second = BoxFixture::new(["3", "3", "3"]);
    second.origin[0] = r("20");
    let text = boxes_step(&[first, second]);
    let bytes = text.as_bytes();
    let format = SourceFormat::Step;
    let entities = format.solid_entities(bytes).unwrap();
    let files = ImportFiles::from([(sha256_hex(bytes), bytes)]);
    let whole = evaluate_with_files(
        &recipe(bytes, format, None),
        &files,
        &mut Decider::default(),
    )
    .unwrap();
    let step = StepId::new("source").unwrap();
    let record = whole.import_record(&step).unwrap();
    assert!(!record.changes.is_empty());
    for entity in &entities {
        let selected = recipe(bytes, format, Some(entity));
        let checkpoint =
            evaluate_incremental_with_files(&selected, None, &files, &mut Decider::default())
                .unwrap();
        assert_eq!(checkpoint.built().import_record(&step).unwrap(), record);
        let mut changed = bytes.to_vec();
        changed.push(b'\n');
        let changed_files = ImportFiles::from([(sha256_hex(bytes), changed.as_slice())]);
        assert!(
            evaluate_incremental_with_files(
                &selected,
                Some(checkpoint),
                &changed_files,
                &mut Decider::default()
            )
            .is_err()
        );
    }
    let proposals = catallax_import::propose(bytes, &Default::default()).unwrap();
    let proposal = &proposals.proposals[0];
    let question = format!("repair:{}:{}", &sha256_hex(bytes)[..12], proposal.id);
    for entity in &entities {
        let mut value = serde_json::to_value(recipe(bytes, format, Some(entity))).unwrap();
        value["steps"][0]["decisions"] = json!({question.clone():"decline"});
        let declined = Recipe::parse(&value.to_string()).unwrap();
        let outcome = evaluate_with_files(&declined, &files, &mut Decider::default());
        let mut model = catallax_topology::Model::new();
        let imported = catallax_import::apply(
            &mut model,
            bytes,
            "reviewed source",
            &Default::default(),
            &[catallax_import::Decision::Decline {
                id: proposal.id.clone(),
            }],
            Some(&step),
        )
        .unwrap();
        if imported.left_open.contains(entity) {
            assert!(
                outcome
                    .unwrap_err()
                    .to_string()
                    .contains("import.left-open")
            );
        } else {
            let built = outcome.unwrap();
            assert_eq!(built.import_record(&step).unwrap(), &imported.record);
        }
    }
}

#[test]
fn selected_replacement_cannot_reuse_an_entity_key_in_different_bytes() {
    use catallax_features::import_edit::{StepReplacement, replace_step};
    use catallax_recipe::FileRef;
    use std::collections::BTreeMap;

    let larger = || {
        let mut solid = BoxFixture::new(["2", "3", "4"]);
        solid.origin[0] = r("5");
        solid
    };
    let old = boxes_step(&[BoxFixture::new(["1", "1", "1"]), larger()]);
    let new = boxes_step(&[larger(), BoxFixture::new(["1", "1", "1"])]);
    for format in [SourceFormat::Step, SourceFormat::Iges] {
        let encode = |text: &str| {
            if format == SourceFormat::Step {
                return text.as_bytes().to_vec();
            }
            let source = SourceFormat::Step.read(text.as_bytes()).unwrap();
            let mut model = catallax_topology::Model::new();
            let bodies = catallax_import::build_file(&mut model, &source, None).unwrap();
            let mut bytes = Vec::new();
            catallax_iges::write_iges(
                &model,
                &bodies,
                &Default::default(),
                &mut Decider::default(),
                &mut bytes,
            )
            .unwrap();
            bytes
        };
        let old = encode(&old);
        let new = encode(&new);
        let entities = format.solid_entities(&old).unwrap();
        // Swapping equal-topology solids preserves keys while changing their geometry.
        assert_eq!(entities, format.solid_entities(&new).unwrap());
        let selected = recipe(&old, format, Some(&entities[0]));
        let files = ImportFiles::from([
            (sha256_hex(&old), old.as_slice()),
            (sha256_hex(&new), new.as_slice()),
        ]);
        assert_eq!(
            measured(&evaluate_with_files(&selected, &files, &mut Decider::default()).unwrap()),
            r("1")
        );
        let edit = StepReplacement {
            step: StepId::new("source").unwrap(),
            old_sha256: sha256_hex(&old),
            file: FileRef {
                name: "replacement".into(),
                sha256: sha256_hex(&new),
            },
            decisions: BTreeMap::new(),
            faces: BTreeMap::new(),
        };
        let saved = selected.to_text();
        match replace_step(&selected, &edit, &files, &mut Decider::default()) {
            Err(error) => assert_eq!(error.code, "import.replace-source-solid"),
            Ok(built) => panic!(
                "replacement silently selected volume {} instead of 1",
                measured(&built)
            ),
        }
        assert_eq!(selected.to_text(), saved);
        let same_source = StepReplacement {
            file: FileRef {
                name: "same original bytes".into(),
                sha256: sha256_hex(&old),
            },
            ..edit
        };
        let rebuilt =
            replace_step(&selected, &same_source, &files, &mut Decider::default()).unwrap();
        assert_eq!(measured(&rebuilt), r("1"));
    }
}
