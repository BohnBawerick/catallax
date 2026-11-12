//! A `scan_import` step is reference data beside the bodies (RECIPE.md section
//! 4.19.1): both kernels build the recipe's body result without the scan's
//! bytes, because no body operation may read a scan.

use catallax_field::Plan;
use catallax_numbers::Decider;
use catallax_recipe::Recipe;

const RECIPE: &str = r#"{
  "format": "catallax-recipe/4",
  "name": "scan beside a block",
  "steps": [
    { "id": "scan", "op": "scan_import",
      "file": { "name": "scan.obj", "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855" },
      "pose": { "rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"] }, "scale": "1" },
    { "id": "block", "op": "block", "size": ["4", "3", "2"],
      "pose": { "rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"] } }
  ],
  "result": "block"
}"#;

#[test]
fn both_kernels_build_the_body_beside_a_scan_reference() {
    let recipe = Recipe::parse(RECIPE).unwrap();
    let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    assert_eq!(built.model.body(built.result).lumps.len(), 1);
    Plan::compile_checked(&recipe, &mut Decider::default()).unwrap();
}
