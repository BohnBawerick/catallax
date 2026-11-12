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

/// A pointed cone halved through its axis keeps its apex on the cut, which no
/// cone chart can carry.
fn halved_cone() -> Value {
    json!({"format": "catallax-recipe/3", "name": "halved cone", "steps": [
        {"id": "stock", "op": "cone", "base_radius": "3", "top_radius": "0", "height": "4",
         "pose": pose(["0", "0", "0"])},
        {"id": "cut", "op": "block", "size": ["10", "20", "20"], "pose": pose(["0", "-10", "-10"])},
        {"id": "half", "op": "boolean", "kind": "intersect", "target": "stock", "tool": "cut"}],
        "result": "half"})
}

/// Faces the current encoders leave without pcurves: a cone with its apex on
/// a trim, and torus rounds beside a cylinder.
fn unseamed() -> Vec<(&'static str, Value, &'static str)> {
    vec![
        ("halved cone", halved_cone(), "CONICAL_SURFACE 1"),
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
    let (halved, cut) = built(&unseamed().swap_remove(0).1);
    let mut first = None;
    for threads in [1, 2, 4] {
        catallax_numbers::set_threads(std::num::NonZeroUsize::new(threads));
        let answer = (
            write(&model, body, Seams::Required).unwrap(),
            write(&halved, cut, Seams::Required).unwrap_err(),
        );
        let first = first.get_or_insert_with(|| answer.clone());
        assert_eq!(first.0, answer.0, "{threads} threads");
        assert_eq!(first.1.detail, answer.1.detail, "{threads} threads");
    }
    catallax_numbers::set_threads(None);
}

/// Cones trimmed by planes, curved limits and an offset of either: the
/// fixtures whose only unseamed faces were cones, or cones beside sphere and
/// torus faces that the cone charts carry.
fn cone_trims() -> Vec<(&'static str, Value)> {
    vec![
        (
            "hyperbola",
            fixture(include_str!(
                "../../catallax-app/tests/fixtures/open-conics/hyperbola.json"
            )),
        ),
        (
            "parabola",
            fixture(include_str!(
                "../../catallax-app/tests/fixtures/open-conics/parabola.json"
            )),
        ),
        (
            "oblique inward",
            fixture(include_str!(
                "../../catallax-features/tests/fixtures/cone-distance/oblique-inward.json"
            )),
        ),
        (
            "tip outward",
            fixture(include_str!(
                "../../catallax-features/tests/fixtures/cone-distance/tip-outward.json"
            )),
        ),
        (
            "tip outward rational",
            fixture(include_str!(
                "../../catallax-features/tests/fixtures/cone-distance/tip-outward-rational.json"
            )),
        ),
        (
            "cone entry next",
            fixture(include_str!(
                "../../catallax-features/tests/fixtures/curved-limits/cone-entry-next.json"
            )),
        ),
        (
            "cone entry to",
            fixture(include_str!(
                "../../catallax-features/tests/fixtures/curved-limits/cone-entry-to.json"
            )),
        ),
        (
            "pentagon",
            fixture(include_str!(
                "../../catallax-features/tests/fixtures/polygon-limits/cone-pentagon.json"
            )),
        ),
        (
            "rotated",
            fixture(include_str!(
                "../../catallax-features/tests/fixtures/polygon-limits/cone-rotated.json"
            )),
        ),
        (
            "triangle",
            fixture(include_str!(
                "../../catallax-features/tests/fixtures/polygon-limits/cone-triangle.json"
            )),
        ),
    ]
}

#[test]
fn trimmed_cones_write_their_seams_under_the_contract_and_keep_the_legacy_file() {
    for (name, recipe) in cone_trims() {
        let (model, body) = built(&recipe);
        let legacy = write(&model, body, Seams::Legacy).unwrap();
        // Version 2 keeps its encoding, which leaves the cone face without pcurves.
        assert!(check_seams(&legacy).is_err(), "{name}");
        let required = write(&model, body, Seams::Required)
            .unwrap_or_else(|refusal| panic!("{name}: {} {}", refusal.code, refusal.detail));
        assert!(check_seams(&required).is_ok(), "{name}");
        let text = String::from_utf8(required).unwrap();
        assert!(text.contains("CONICAL_SURFACE"), "{name}");
        assert!(text.contains("PCURVE"), "{name}");
    }
}

#[test]
fn trimmed_cones_answer_alike_at_one_two_and_four_threads() {
    let parts: Vec<_> = cone_trims()
        .into_iter()
        .map(|(_, recipe)| built(&recipe))
        .collect();
    let mut first: Option<Vec<Vec<u8>>> = None;
    for threads in [1, 2, 4] {
        catallax_numbers::set_threads(std::num::NonZeroUsize::new(threads));
        let answer: Vec<Vec<u8>> = parts
            .iter()
            .map(|(model, body)| write(model, *body, Seams::Required).unwrap())
            .collect();
        let first = first.get_or_insert_with(|| answer.clone());
        assert_eq!(*first, answer, "{threads} threads");
    }
    catallax_numbers::set_threads(None);
}

/// Oblique extrusions of a circle and of a slot: their walls are elliptic
/// cylinders, the first with two whole rings and the second with half walls
/// whose end lines lie along the chart's seam.
fn oblique_extrusions() -> Vec<(&'static str, Value)> {
    vec![
        (
            "circle",
            fixture(include_str!(
                "../../catallax-features/tests/fixtures/algebraic-extrusions/circle.json"
            )),
        ),
        (
            "slot",
            fixture(include_str!(
                "../../catallax-features/tests/fixtures/algebraic-extrusions/slot.json"
            )),
        ),
    ]
}

#[test]
fn oblique_extrusions_write_their_seams_under_the_contract_and_keep_the_legacy_file() {
    for (name, recipe) in oblique_extrusions() {
        let (model, body) = built(&recipe);
        let legacy = write(&model, body, Seams::Legacy).unwrap();
        assert!(check_seams(&legacy).is_err(), "{name}");
        let required = write(&model, body, Seams::Required)
            .unwrap_or_else(|refusal| panic!("{name}: {} {}", refusal.code, refusal.detail));
        assert!(check_seams(&required).is_ok(), "{name}");
        let text = String::from_utf8(required).unwrap();
        assert!(text.contains("SURFACE_OF_LINEAR_EXTRUSION"), "{name}");
        assert_eq!(text.contains("SEAM_CURVE"), name == "circle", "{name}");
    }
}

#[test]
fn oblique_extrusions_answer_alike_at_one_two_and_four_threads() {
    let parts: Vec<_> = oblique_extrusions()
        .into_iter()
        .map(|(_, recipe)| built(&recipe))
        .collect();
    let mut first: Option<Vec<Vec<u8>>> = None;
    for threads in [1, 2, 4] {
        catallax_numbers::set_threads(std::num::NonZeroUsize::new(threads));
        let answer: Vec<Vec<u8>> = parts
            .iter()
            .map(|(model, body)| write(model, *body, Seams::Required).unwrap())
            .collect();
        let first = first.get_or_insert_with(|| answer.clone());
        assert_eq!(*first, answer, "{threads} threads");
    }
    catallax_numbers::set_threads(None);
}
