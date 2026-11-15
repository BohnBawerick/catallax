//! Recipe shortcuts must preserve the distinction between surfaces and material.
#[path = "../../catallax-import/tests/common/mod.rs"]
mod common;

use catallax_features::{BuildError, ImportFiles, evaluate_with_files};
use catallax_import::sha256_hex;
use catallax_numbers::Decider;
use catallax_recipe::Recipe;
use serde_json::{Value, json};

fn source() -> String {
    common::box_step(common::BoxFixture::new(["1", "2", "3"]))
        .lines()
        .map(|line| {
            if line.contains("MANIFOLD_SOLID_BREP") {
                let (id, _) = line.split_once('=').unwrap();
                let shell = line.split(",#").nth(1).unwrap().trim_end_matches(");");
                format!("{id}=SHELL_BASED_SURFACE_MODEL('',(#{shell}));")
            } else {
                line.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn evaluate(operation: Value) -> Result<catallax_features::Built, BuildError> {
    let source = source();
    let bytes = source.as_bytes();
    let hash = sha256_hex(bytes);
    let recipe = Recipe::parse(
        &json!({
            "format":"catallax-recipe/5", "name":"surface operation",
            "steps":[{"id":"source", "op":"import_step",
                "file":{"name":"surface.step", "sha256":hash},
                "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
                operation], "result":"result"
        })
        .to_string(),
    )
    .unwrap();
    let files = ImportFiles::from([(hash, bytes)]);
    evaluate_with_files(&recipe, &files, &mut Decider::default())
}

fn refuses(operation: Value) {
    let error = evaluate(operation)
        .err()
        .expect("surface operand must refuse");
    let BuildError::Geometry { refusal, .. } = error else {
        panic!("expected geometry refusal, got {error:?}");
    };
    assert_eq!(refusal.code, "boolean.not-solid");
}

#[test]
fn separated_pattern_cannot_turn_closed_surfaces_into_material() {
    refuses(json!({"id":"result", "op":"pattern", "body":"source",
        "layout":{"linear":{"vector":["10","0","0"], "count":2}}}));
}

#[test]
fn separated_mirror_cannot_turn_closed_surfaces_into_material() {
    refuses(json!({"id":"result", "op":"mirror", "body":"source",
        "plane":{"point":["10","0","0"],"normal":["1","0","0"]},
        "keep_original":true}));
}

#[test]
fn self_boolean_identities_refuse_surface_operands() {
    for kind in ["union", "subtract", "intersect"] {
        refuses(json!({"id":"result", "op":"boolean", "kind":kind,
            "target":"source", "tool":"source"}));
    }
}

#[test]
fn placement_without_union_keeps_the_surface_kind() {
    for operation in [
        json!({"id":"result", "op":"pattern", "body":"source",
            "layout":{"linear":{"vector":["10","0","0"], "count":1}}}),
        json!({"id":"result", "op":"mirror", "body":"source",
            "plane":{"point":["10","0","0"],"normal":["1","0","0"]},
            "keep_original":false}),
    ] {
        let built = evaluate(operation).unwrap();
        assert_eq!(
            built.model.body(built.result).kind,
            catallax_topology::BodyKind::Surface
        );
        assert_eq!(
            catallax_measure::volume(
                &built.model,
                built.result,
                &Default::default(),
                &mut Decider::default()
            )
            .unwrap_err()
            .code,
            "measure.not-solid"
        );
    }
}
