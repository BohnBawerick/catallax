//! STEP files written under the seam contract import again: the reader drops
//! each `SEAM_CURVE` used twice by one face and rebuilds the same body.

use catallax_features::{ImportFiles, evaluate, evaluate_with_files};
use catallax_import::sha256_hex;
use catallax_measure::{Accuracy, volume};
use catallax_numbers::Decider;
use catallax_recipe::Recipe;
use catallax_step::{Seams, StepOptions, write_step};
use catallax_topology::{BodyId, Model};
use serde_json::{Value, json};

fn pose(translation: [&str; 3]) -> Value {
    json!({"rotation": ["1", "0", "0", "0"], "translation": translation})
}

fn recipe(name: &str, steps: Value) -> Recipe {
    Recipe::parse(
        &json!({"format": "catallax-recipe/3", "name": name, "steps": steps, "result": "body"})
            .to_string(),
    )
    .unwrap()
}

/// Seamed parts: name, recipe and whether the file holds a `SEAM_CURVE`.
fn parts() -> Vec<(&'static str, Recipe, bool)> {
    let origin = pose(["0", "0", "0"]);
    vec![
        (
            "cylinder",
            recipe(
                "cylinder",
                json!([{"id": "body", "op": "cylinder", "radius": "1", "height": "2",
                    "pose": pose(["1", "2", "3"])}]),
            ),
            true,
        ),
        (
            "frustum",
            recipe(
                "frustum",
                json!([{"id": "body", "op": "cone", "base_radius": "3", "top_radius": "1",
                    "height": "4", "pose": origin}]),
            ),
            true,
        ),
        (
            "pointed cone",
            recipe(
                "pointed cone",
                json!([{"id": "body", "op": "cone", "base_radius": "3", "top_radius": "0",
                    "height": "4", "pose": origin}]),
            ),
            true,
        ),
        (
            "sphere",
            recipe(
                "sphere",
                json!([{"id": "body", "op": "sphere", "radius": "2", "pose": origin}]),
            ),
            false,
        ),
        (
            "torus",
            recipe(
                "torus",
                json!([{"id": "body", "op": "torus", "major": "3", "minor": "1", "pose": origin}]),
            ),
            true,
        ),
        (
            "through hole",
            recipe(
                "through hole",
                json!([
                    {"id": "base", "op": "block", "size": ["4", "4", "2"], "pose": origin},
                    {"id": "tool", "op": "cylinder", "radius": "1", "height": "4",
                        "pose": pose(["2", "2", "-1"])},
                    {"id": "body", "op": "boolean", "kind": "subtract", "target": "base",
                        "tool": "tool"}]),
            ),
            true,
        ),
        (
            "blind drill tip",
            recipe(
                "blind drill tip",
                json!([
                    {"id": "plate", "op": "block", "size": ["60", "40", "10"],
                        "pose": pose(["-30", "-20", "0"])},
                    {"id": "body", "op": "hole", "target": "plate",
                        "at": [{"id": "h0", "point": ["10", "0", "10"], "direction": ["0", "0", "-1"]}],
                        "diameter": "4", "depth": {"blind": "5"}, "tip": {"cone": "90"}}]),
            ),
            true,
        ),
    ]
}

fn seamed(model: &Model, body: BodyId) -> Vec<u8> {
    let mut bytes = Vec::new();
    let options = StepOptions {
        seams: Seams::Required,
        ..StepOptions::default()
    };
    write_step(
        model,
        &[body],
        &options,
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap_or_else(|refusal| panic!("{}: {}", refusal.code, refusal.detail));
    bytes
}

fn kinds(model: &Model, body: BodyId) -> Vec<String> {
    let mut kinds: Vec<_> = model
        .faces_of_body(body)
        .into_iter()
        .map(|face| {
            format!("{:?}", model.geom().surface(model.face(face).surface))
                .chars()
                .take_while(|c| c.is_alphanumeric())
                .collect()
        })
        .collect();
    kinds.sort();
    kinds
}

#[test]
fn seamed_files_import_as_the_same_body() {
    for (name, source, seam) in parts() {
        let built = evaluate(&source, &mut Decider::default()).unwrap();
        let bytes = seamed(&built.model, built.result);
        assert_eq!(
            String::from_utf8_lossy(&bytes).contains("SEAM_CURVE"),
            seam,
            "{name}"
        );
        let hash = sha256_hex(&bytes);
        let imported = Recipe::parse(
            &json!({"format": "catallax-recipe/5", "name": name,
                "steps": [{"id": "body", "op": "import_step",
                    "file": {"name": "part.step", "sha256": hash}, "pose": pose(["0", "0", "0"])}],
                "result": "body"})
            .to_string(),
        )
        .unwrap();
        let files = ImportFiles::from([(hash, &bytes[..])]);
        let again = evaluate_with_files(&imported, &files, &mut Decider::default())
            .unwrap_or_else(|error| panic!("{name}: {error:?}"));
        assert!(again.model.check(again.result).is_empty(), "{name}");
        assert_eq!(
            kinds(&built.model, built.result),
            kinds(&again.model, again.result),
            "{name}"
        );
        let measure = |model: &Model, body| {
            volume(model, body, &Accuracy::default(), &mut Decider::default())
                .unwrap()
                .value
        };
        let (a, b) = (
            measure(&built.model, built.result),
            measure(&again.model, again.result),
        );
        assert!(a.lo <= b.hi && b.lo <= a.hi, "{name}: {a:?} {b:?}");
        // The imported body writes a seamed file again.
        seamed(&again.model, again.result);
    }
}
