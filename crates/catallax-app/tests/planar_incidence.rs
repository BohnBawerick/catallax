//! Recorded planar source repair through recipe evaluation and checkpoint reuse.

#[path = "../../catallax-import/tests/common/mod.rs"]
mod common;

use catallax_app::build_incremental_with_files;
use catallax_features::ImportFiles;
use catallax_import::sha256_hex;
use catallax_mesh::MeshOptions;
use catallax_numbers::Decider;
use catallax_recipe::{Recipe, StepId};
use serde_json::{Value, json};

fn recipe(hash: &str, decisions: Value) -> Recipe {
    Recipe::parse(
        &json!({
            "format": "catallax-recipe/3", "name": "planar source control",
            "steps": [{"id": "part", "op": "import_step",
                "file": {"name": "source.step", "sha256": hash},
                "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]},
                "decisions": decisions}], "result": "part"
        })
        .to_string(),
    )
    .unwrap()
}

#[test]
fn recipe_build_and_checkpoint_reuse_retain_the_source_definition() {
    let bytes = common::box_with_placement(5, Some(["0.000001", "0", "1"]), Some(["1", "0", "0"]));
    let hash = sha256_hex(bytes.as_bytes());
    let files = ImportFiles::from([(hash.clone(), bytes.as_bytes())]);
    let recipe = recipe(&hash, json!({}));
    let cold = build_incremental_with_files(
        &recipe,
        None,
        &files,
        &MeshOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let step = StepId::new("part").unwrap();
    let record = cold
        .checkpoint
        .built()
        .import_record(&step)
        .unwrap()
        .clone();
    assert_eq!(record.format, "catallax-import-repair/3");
    assert_eq!(record.source_incidence.len(), 1);
    assert_eq!(
        record.source_incidence[0].source_definition,
        "projected_vertex_boundary/1"
    );
    let warm = build_incremental_with_files(
        &recipe,
        Some(cold.checkpoint),
        &files,
        &MeshOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(warm.checkpoint.built().import_record(&step), Some(&record));
    let question = format!("repair:{}:p0", &hash[..12]);
    let declined = recipe_for_decline(&hash, &question);
    let error = build_incremental_with_files(
        &declined,
        None,
        &files,
        &MeshOptions::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code(), "import.left-open");
}

fn recipe_for_decline(hash: &str, question: &str) -> Recipe {
    recipe(hash, json!({question: "decline"}))
}

#[test]
fn an_unbounded_planar_repair_keeps_its_typed_recipe_refusal() {
    let bytes = common::box_with_placement(5, Some(["1", "0", "1"]), Some(["1", "0", "0"]));
    let hash = sha256_hex(bytes.as_bytes());
    let files = ImportFiles::from([(hash.clone(), bytes.as_bytes())]);
    let error = build_incremental_with_files(
        &recipe(&hash, json!({})),
        None,
        &files,
        &MeshOptions::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code(), "import.planar-repair-not-verified");
}
