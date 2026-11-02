//! V9-K acceptance K1-K3, from the independent V9 batch 1 design.
use catallax_app::{Error, build_incremental_with_files};
use catallax_features::{BuildError, ImportFiles};
use catallax_import::{
    Approval, Decision, Disposition, ImportOptions, RepairRecord, apply, propose, rebuild,
    sha256_hex,
};
use catallax_measure::{Accuracy, volume};
use catallax_mesh::MeshOptions;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Origin, Recipe, Role, StepId};
use catallax_step::{StepOptions, write_step};
use catallax_topology::Model;
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

fn seam_source(height: &str, gap: &str) -> Vec<u8> {
    let height = Rational::parse(height).unwrap();
    let top = &height + &Rational::parse(gap).unwrap();
    std::str::from_utf8(BLOCK)
        .unwrap()
        .replace("10.004", &top.to_decimal_exact().unwrap())
        .replace(
            ",10.)",
            &format!(",{})", height.to_decimal_exact().unwrap()),
        )
        .into_bytes()
}

fn qualify(source: &[u8], tolerance: &str, approved: bool) {
    let hash = sha256_hex(source);
    let options = ImportOptions::default();
    let proposals = propose(source, &options).unwrap();
    let [proposal] = proposals.proposals.as_slice() else {
        panic!("one seam")
    };
    assert_eq!(
        proposal.disposition,
        if approved {
            Disposition::Ask
        } else {
            Disposition::Auto
        }
    );
    let decisions = if approved {
        vec![Decision::Accept {
            id: proposal.id.clone(),
            tolerance_mm: tolerance.into(),
            approver: "designer".into(),
        }]
    } else {
        vec![]
    };
    let imported = apply(
        &mut Model::new(),
        source,
        "source.step",
        &options,
        &decisions,
        None,
    )
    .unwrap();
    assert_eq!(sha256_hex(source), hash);
    let record: RepairRecord =
        serde_json::from_slice(&serde_json::to_vec(&imported.record).unwrap()).unwrap();
    let reopened = rebuild(&mut Model::new(), source, &options, &record, None).unwrap();
    assert_eq!(reopened.derivative, imported.derivative);
    assert_eq!(reopened.record, record);
    assert_eq!(
        record.changes[0].approval,
        if approved {
            Approval::User
        } else {
            Approval::Automatic
        }
    );
    assert_eq!(record.changes[0].max_movement_mm, tolerance);
    let mut keys = record
        .face_tolerances
        .iter()
        .map(|f| {
            assert_eq!(f.deviation_lower_mm, tolerance);
            assert_eq!(f.deviation_upper_mm, tolerance);
            let upper = Rational::parse(&f.deviation_upper_mm).unwrap();
            assert!(
                &upper * &upper * Rational::from_i64(100)
                    <= Rational::parse(&f.feature_size_squared_mm2).unwrap()
            );
            f.face.as_str()
        })
        .collect::<Vec<_>>();
    keys.sort();
    assert_eq!(keys, ["#102", "#124", "#40", "#80"]);
    let answers = if approved {
        json!({&proposal.question: format!("accept:{tolerance}:designer")})
    } else {
        json!({})
    };
    let definition = recipe(&hash, answers);
    let saved = Recipe::parse(&definition.to_text()).unwrap();
    let files = ImportFiles::from([(hash, source)]);
    let app = build(&saved, &files).unwrap();
    let step = StepId::new("block").unwrap();
    let built = app.checkpoint.built();
    assert_eq!(built.import_record(&step).unwrap(), &record);
    assert!(built.model.check(built.result).is_empty());
    for face in &record.face_tolerances {
        assert!(
            app.names
                .faces
                .values()
                .any(|origins| origins.iter().any(|o| o.step == step
                    && o.role
                        == Role::Imported {
                            key: face.face.clone()
                        }))
        );
    }
    // A display-only checkpoint reuse must retain the source-relative evidence.
    let remeshed = build_incremental_with_files(
        &saved,
        Some(app.checkpoint.clone()),
        &files,
        &MeshOptions {
            chord_tolerance: Rational::parse("0.005").unwrap(),
            ..Default::default()
        },
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(
        remeshed.checkpoint.built().import_record(&step).unwrap(),
        &record
    );
    let cold = build(&saved, &files).unwrap();
    assert_eq!(
        cold.checkpoint.built().import_record(&step).unwrap(),
        &record
    );
    let mut placed: serde_json::Value = serde_json::from_str(&saved.to_text()).unwrap();
    placed["steps"][0]["pose"] = json!({
        "rotation": ["0", "0", "0", "1"], "translation": ["7", "-3", "2"]
    });
    let placed = build(&Recipe::parse(&placed.to_string()).unwrap(), &files).unwrap();
    let placed = placed.checkpoint.built();
    assert_eq!(placed.import_record(&step).unwrap(), &record);
    let bounds = catallax_measure::bounding_box(
        &placed.model,
        placed.result,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value
    .unwrap();
    assert_eq!(bounds.min[0], Rational::from_i64(-33));
    assert_eq!(bounds.max[1], Rational::from_i64(-3));

    let mut output = Vec::new();
    let exported = write_step(
        &built.model,
        &[built.result],
        &StepOptions::default(),
        &mut Decider::default(),
        &mut output,
    )
    .unwrap()
    .value;
    // The existing exchange promise is relative to the repaired model, not source.
    assert_eq!(exported.max_deviation, Rational::zero());
    let back = apply(
        &mut Model::new(),
        &output,
        "export.step",
        &options,
        &[],
        None,
    )
    .unwrap();
    assert_eq!(back.bodies.len(), 1);
    assert!(back.record.face_tolerances.is_empty());

    let mut forged = record.clone();
    forged.face_tolerances[0].deviation_upper_mm = "0".into();
    assert!(matches!(
        rebuild(&mut Model::new(), source, &options, &forged, None),
        Err(catallax_import::ImportRefusal::ValidationFailed { .. })
    ));
    let declined = recipe(
        &record.source.sha256,
        json!({&proposal.question: "decline"}),
    );
    let declined = Recipe::parse(&declined.to_text()).unwrap();
    geometry(
        build(&declined, &files).unwrap_err(),
        "import.left-open",
        &["#178"],
    );
}

#[test]
fn planar_repair_evidence_survives_application_save_reopen_display_and_export() {
    qualify(BLOCK, "0.004", false);
    qualify(&seam_source("10", "0.01"), "0.01", false);
    qualify(&seam_source("10", "0.010001"), "0.010001", true);
    qualify(&seam_source("10", "0.02"), "0.02", true);
}

#[test]
fn nearest_feature_limit_cannot_be_overridden_by_saved_approval() {
    qualify(&seam_source("0.05", "0.005"), "0.005", false);
    let source = seam_source("0.05", "0.005001");
    let hash = sha256_hex(&source);
    let proposal = propose(&source, &ImportOptions::default())
        .unwrap()
        .proposals
        .remove(0);
    assert_eq!(proposal.disposition, Disposition::Ask);
    assert!(proposal.reason.contains("tenth"));
    let files = ImportFiles::from([(hash.clone(), source.as_slice())]);
    geometry(
        build(
            &recipe(
                &hash,
                json!({&proposal.question: "accept:0.005001:designer"}),
            ),
            &files,
        )
        .unwrap_err(),
        "import.validation-failed",
        &["tenth"],
    );
}
