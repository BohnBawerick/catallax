//! `scan_import`, RECIPE.md section 4.19.1: a reference value that is never a
//! body, never the result and never the input of a body operation.

use catallax_recipe::{Operation, Recipe, RecipeError, ValueKind};

const HASH: &str = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855";

fn recipe(format: &str, hash: &str, last: &str, result: &str) -> String {
    format!(
        r#"{{
  "format": "{format}",
  "name": "scan reference",
  "steps": [
    {{ "id": "scan", "op": "scan_import", "file": {{ "name": "part.obj", "sha256": "{hash}" }},
       "pose": {{ "rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"] }}, "scale": "1" }},
    {{ "id": "block", "op": "block", "size": ["1", "2", "3"],
       "pose": {{ "rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"] }} }}{last}
  ],
  "result": "{result}"
}}"#
    )
}

fn code(result: Result<Recipe, RecipeError>) -> &'static str {
    match result {
        Ok(_) => "ok",
        Err(e) => e.code(),
    }
}

#[test]
fn a_scan_is_reference_data_beside_a_body_result() {
    let text = recipe("catallax-recipe/4", HASH, "", "block");
    let parsed = Recipe::parse(&text).expect("a scan beside a body result is valid");
    let scan = &parsed.steps[0].op;
    assert!(matches!(scan, Operation::ScanImport { .. }));
    assert_eq!(scan.name(), "scan_import");
    assert_eq!(scan.kind(), ValueKind::Scan);
    assert_eq!(ValueKind::Scan.name(), "scan");
    assert!(scan.inputs().is_empty());
    assert_eq!(Recipe::parse(&parsed.to_text()).unwrap(), parsed);
    assert_eq!(
        code(Recipe::parse(&recipe(
            "catallax-recipe/5",
            HASH,
            "",
            "block"
        ))),
        "ok"
    );
}

#[test]
fn a_scan_is_never_a_result_a_body_input_or_a_version_3_step() {
    assert_eq!(
        code(Recipe::parse(&recipe(
            "catallax-recipe/4",
            HASH,
            "",
            "scan"
        ))),
        "recipe.result"
    );
    let boolean = r#",
    { "id": "joined", "op": "boolean", "kind": "union", "target": "block", "tool": "scan" }"#;
    assert_eq!(
        code(Recipe::parse(&recipe(
            "catallax-recipe/4",
            HASH,
            boolean,
            "joined"
        ))),
        "recipe.wrong-input"
    );
    let offset = r#",
    { "id": "grown", "op": "offset", "body": "scan", "distance": "1" }"#;
    assert_eq!(
        code(Recipe::parse(&recipe(
            "catallax-recipe/4",
            HASH,
            offset,
            "grown"
        ))),
        "recipe.wrong-input"
    );
    assert_eq!(
        code(Recipe::parse(&recipe(
            "catallax-recipe/3",
            HASH,
            "",
            "block"
        ))),
        "recipe.parse"
    );
    assert_eq!(
        code(Recipe::parse(&recipe(
            "catallax-recipe/4",
            "ABC",
            "",
            "block"
        ))),
        "recipe.file-hash"
    );
}

#[test]
fn a_scan_scale_that_is_not_positive_is_a_recipe_error() {
    for scale in ["0", "-1"] {
        let text = recipe("catallax-recipe/4", HASH, "", "block")
            .replace(r#""scale": "1""#, &format!(r#""scale": "{scale}""#));
        assert_eq!(code(Recipe::parse(&text)), "recipe.parse", "{scale}");
    }
}
