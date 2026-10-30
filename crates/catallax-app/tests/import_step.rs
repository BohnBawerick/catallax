//! V9-K acceptance K1-K3, from the independent V9 batch 1 design.
use catallax_app::{Error, build_incremental_with_files};
use catallax_features::{BuildError, ImportFiles};
use catallax_measure::{Accuracy, volume};
use catallax_mesh::MeshOptions;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Origin, Recipe, Role, StepId};
use serde_json::json;

const BLOCK: &[u8] = include_bytes!("fixtures/block-seam-4um.step");
const TWO_BLOCKS: &[u8] = include_bytes!("fixtures/two-blocks-5um.step");
const H1: &str = "1a602addaa4b534ef9427427d5231ba689a811a11159a5f9e7ebc17b3365c32b";
const H2: &str = "0a16eb4caa0073701b15903a503d8b75e629936806506918cae1814ac162a0bd";
const QUESTION: &str = "repair:1a602addaa4b:g0";

fn recipe(hash: &str, decisions: serde_json::Value) -> Recipe {
    Recipe::parse(
        &json!({
            "format": "catallax-recipe/3",
            "name": "reviewed import",
            "steps": [{
                "id": "block", "op": "import_step",
                "file": {"name": "source.step", "sha256": hash},
                "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]},
                "decisions": decisions
            }],
            "result": "block"
        })
        .to_string(),
    )
    .unwrap()
}

fn build(
    recipe: &Recipe,
    files: &ImportFiles<'_>,
) -> Result<catallax_app::IncrementalApplicationBody, Error> {
    build_incremental_with_files(
        recipe,
        None,
        files,
        &MeshOptions::default(),
        &mut Decider::default(),
    )
}

fn geometry(error: Error, code: &str, names: &[&str]) {
    assert_eq!(error.code(), code);
    let Error::Build(BuildError::Geometry { step, refusal }) = error else {
        panic!("expected a step-local typed import refusal")
    };
    assert_eq!(step.as_str(), "block");
    for name in names {
        assert!(
            refusal.detail.contains(name),
            "{} lacks {name}",
            refusal.detail
        );
    }
}

#[test]
fn k1_repaired_import_builds_exact_volume_and_file_face_origins() {
    let files = ImportFiles::from([(H1.to_string(), BLOCK)]);
    let body = build(&recipe(H1, json!({})), &files).unwrap();
    let built = body.checkpoint.built();
    let measured = volume(
        &built.model,
        built.result,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert_eq!(measured.lo, Rational::from_i64(10004));
    assert_eq!(measured.hi, Rational::from_i64(10004));
    let mut origins = body
        .names
        .faces
        .values()
        .map(|origins| {
            assert_eq!(origins.len(), 1);
            origins[0].clone()
        })
        .collect::<Vec<_>>();
    origins.sort();
    let mut expected = ["#40", "#80", "#102", "#124", "#136", "#176"].map(|key| Origin {
        step: StepId::new("block").unwrap(),
        role: Role::Imported { key: key.into() },
    });
    expected.sort();
    assert_eq!(origins, expected);
}

#[test]
fn k2_decline_leaves_no_body_and_accept_is_not_a_saved_automatic_answer() {
    let files = ImportFiles::from([(H1.to_string(), BLOCK)]);
    geometry(
        build(&recipe(H1, json!({QUESTION: "decline"})), &files).unwrap_err(),
        "import.left-open",
        &["#178"],
    );
    geometry(
        build(&recipe(H1, json!({QUESTION: "accept:0.01"})), &files).unwrap_err(),
        "import.invalid-decision",
        &[QUESTION],
    );
}

#[test]
fn k3_multiple_solids_and_unavailable_or_changed_source_refuse() {
    let files = ImportFiles::from([(H2.to_string(), TWO_BLOCKS)]);
    geometry(
        build(&recipe(H2, json!({})), &files).unwrap_err(),
        "import.multibody",
        &["2 solids", "#150", "#300"],
    );
    let recipe = recipe(H1, json!({}));
    let mut changed = BLOCK.to_vec();
    changed.push(b'\n');
    for files in [
        ImportFiles::new(),
        ImportFiles::from([(H1.to_string(), changed.as_slice())]),
    ] {
        let error = build(&recipe, &files).unwrap_err();
        assert!(
            matches!(error, Error::Build(BuildError::Recipe(_))),
            "{error:?}"
        );
        assert!(error.code().starts_with("recipe."));
        assert!(error.to_string().contains("sha256"));
    }
}
