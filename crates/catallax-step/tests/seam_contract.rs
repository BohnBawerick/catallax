//! The seam-curve contract: `Seams::Required` writes the legacy bytes when every
//! curved face already carries its pcurves and seams, and otherwise refuses
//! `step.seam-unsupported` before any byte leaves.

use catallax_numbers::Decider;
use catallax_recipe::Recipe;
use catallax_step::{Seams, StepOptions, check_seams, write_step};
use catallax_topology::{BodyId, Model};
use serde_json::{Value, json};

fn pose(translation: [&str; 3]) -> Value {
    json!({"rotation": ["1", "0", "0", "0"], "translation": translation})
}

fn tilted() -> Value {
    json!({"rotation": ["1", "2", "3", "4"], "translation": ["7/3", "-11/7", "13/5"]})
}

fn primitive(op: &str, fields: Value, pose: Value) -> Value {
    let mut step = json!({"id": "body", "op": op, "pose": pose});
    step.as_object_mut()
        .unwrap()
        .extend(fields.as_object().unwrap().clone());
    json!({"format": "catallax-recipe/3", "name": op, "steps": [step], "result": "body"})
}

fn through_hole() -> Value {
    json!({"format": "catallax-recipe/3", "name": "through hole", "steps": [
        {"id": "base", "op": "block", "size": ["4", "4", "2"], "pose": pose(["0", "0", "0"])},
        {"id": "tool", "op": "cylinder", "radius": "1", "height": "4", "pose": pose(["2", "2", "-1"])},
        {"id": "hole", "op": "boolean", "kind": "subtract", "target": "base", "tool": "tool"}],
        "result": "hole"})
}

fn fixture(text: &str) -> Value {
    serde_json::from_str(text).unwrap()
}

fn built(recipe: &Value) -> (Model, BodyId) {
    let recipe = Recipe::parse(&recipe.to_string()).unwrap();
    let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
    (built.model, built.result)
}

fn options(seams: Seams) -> StepOptions {
    StepOptions {
        seams,
        ..StepOptions::default()
    }
}

fn write(model: &Model, body: BodyId, seams: Seams) -> Result<Vec<u8>, catallax_numbers::Refusal> {
    let mut bytes = Vec::new();
    write_step(
        model,
        &[body],
        &options(seams),
        &mut Decider::default(),
        &mut bytes,
    )?;
    Ok(bytes)
}

/// Seamed parts: (recipe, whether the file must hold a `SEAM_CURVE`).
fn seamed() -> Vec<(&'static str, Value, bool)> {
    vec![
        (
            "cylinder",
            primitive("cylinder", json!({"radius": "1", "height": "2"}), tilted()),
            true,
        ),
        (
            "frustum",
            primitive(
                "cone",
                json!({"base_radius": "3", "top_radius": "1", "height": "4"}),
                tilted(),
            ),
            true,
        ),
        (
            "pointed cone",
            primitive(
                "cone",
                json!({"base_radius": "3", "top_radius": "0", "height": "4"}),
                pose(["0", "0", "0"]),
            ),
            true,
        ),
        // A complete sphere closes at the reader's own pole with one vertex loop.
        (
            "sphere",
            primitive("sphere", json!({"radius": "2"}), tilted()),
            false,
        ),
        (
            "torus",
            primitive(
                "torus",
                json!({"major": "3", "minor": "1"}),
                pose(["0", "0", "0"]),
            ),
            true,
        ),
        ("through hole", through_hole(), true),
        (
            "spline sweep",
            fixture(include_str!(
                "../../catallax-features/tests/fixtures/spline-sweep-planar.json"
            )),
            false,
        ),
    ]
}

#[test]
fn seamed_parts_keep_their_bytes_under_the_seam_contract() {
    for (name, recipe, seam) in seamed() {
        let (model, body) = built(&recipe);
        let legacy = write(&model, body, Seams::Legacy).unwrap();
        let required = write(&model, body, Seams::Required)
            .unwrap_or_else(|refusal| panic!("{name}: {} {}", refusal.code, refusal.detail));
        assert_eq!(legacy, required, "{name}");
        assert!(check_seams(&required).is_ok(), "{name}");
        let text = String::from_utf8(required).unwrap();
        assert_eq!(text.contains("SEAM_CURVE"), seam, "{name}");
    }
}

/// Faces the current encoders leave without pcurves: an elliptic extrusion
/// and torus rounds beside a cylinder.
fn unseamed() -> Vec<(&'static str, Value, &'static str)> {
    vec![
        (
            "oblique circle extrusion",
            fixture(include_str!(
                "../../catallax-features/tests/fixtures/algebraic-extrusions/circle.json"
            )),
            "SURFACE_OF_LINEAR_EXTRUSION 1",
        ),
        (
            "offset coaxial cylinders",
            fixture(include_str!(
                "../../catallax-features/tests/fixtures/cylinder-distance/coaxial-cylinders.json"
            )),
            "TOROIDAL_SURFACE 2",
        ),
    ]
}

#[test]
fn unseamed_files_refuse_and_write_nothing() {
    for (name, recipe, kinds) in unseamed() {
        let (model, body) = built(&recipe);
        let legacy = write(&model, body, Seams::Legacy).unwrap();
        let before = b"existing output".to_vec();
        let mut bytes = before.clone();
        let refusal = write_step(
            &model,
            &[body],
            &options(Seams::Required),
            &mut Decider::default(),
            &mut bytes,
        )
        .unwrap_err();
        assert_eq!(refusal.code, "step.seam-unsupported", "{name}");
        assert!(refusal.detail.contains(kinds), "{name}: {}", refusal.detail);
        assert_eq!(bytes, before, "{name}");
        // The checker reads the legacy file the same way.
        assert_eq!(check_seams(&legacy).unwrap_err().detail, refusal.detail);
    }
}

#[test]
fn the_contract_answers_alike_at_one_two_and_four_threads() {
    let (model, body) = built(&through_hole());
    let (oblique, cut) = built(&unseamed().swap_remove(0).1);
    let mut first = None;
    for threads in [1, 2, 4] {
        catallax_numbers::set_threads(std::num::NonZeroUsize::new(threads));
        let answer = (
            write(&model, body, Seams::Required).unwrap(),
            write(&oblique, cut, Seams::Required).unwrap_err(),
        );
        let first = first.get_or_insert_with(|| answer.clone());
        assert_eq!(first.0, answer.0, "{threads} threads");
        assert_eq!(first.1.detail, answer.1.detail, "{threads} threads");
    }
    catallax_numbers::set_threads(None);
}
