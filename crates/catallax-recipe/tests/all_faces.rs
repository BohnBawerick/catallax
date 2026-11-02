//! Version 5 broadens face selection without rewriting saved origin references.

use catallax_recipe::{FORMAT, FORMAT_V4, FORMAT_V5, Recipe, RecipeError};
use serde_json::{Value, json};

fn source(format: &str, reference: Value) -> Value {
    json!({"format":format,"name":"all faces","steps":[
        {"id":"stock","op":"block","size":["10","8","6"],
         "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
        {"id":"result","op":"shell","body":"stock","open":[reference],"thickness":"1"}
    ],"result":"result"})
}

#[test]
fn optional_face_requires_v5_and_preserves_filters_on_save() {
    for reference in [
        json!({}),
        json!({"near":["5","4","7"]}),
        json!({"touching":[{"step":"stock","role":"end"}]}),
    ] {
        for old in [FORMAT, FORMAT_V4] {
            assert!(matches!(
                Recipe::parse(&source(old, reference.clone()).to_string()),
                Err(RecipeError::Parse { .. })
            ));
        }
        let recipe = Recipe::parse(&source(FORMAT_V5, reference).to_string()).unwrap();
        assert_eq!(Recipe::parse(&recipe.to_text()).unwrap(), recipe);
    }
}

#[test]
fn explicit_v4_to_v5_migration_preserves_operations_names_and_edit() {
    let old = Recipe::parse(
        &source(
            FORMAT_V4,
            json!({"face":{"step":"stock","role":"end"},"near":["5","4","7"]}),
        )
        .to_string(),
    )
    .unwrap();
    let mut migrated = old.clone();
    migrated.format = FORMAT_V5.into();
    migrated.check().unwrap();
    assert_eq!(migrated.steps, old.steps);
    assert_eq!(migrated.result, old.result);
    let mut saved: Value = serde_json::from_str(&migrated.to_text()).unwrap();
    saved["steps"][0]["size"][0] = json!("12");
    let edited = Recipe::parse(&saved.to_string()).unwrap();
    assert_eq!(Recipe::parse(&edited.to_text()).unwrap(), edited);
    assert_eq!(edited.steps[1], old.steps[1]);
}

#[test]
fn all_faces_does_not_skip_touching_origin_validation() {
    let bad = source(FORMAT_V5, json!({"touching":[{"step":"missing"}]}));
    assert!(matches!(
        Recipe::parse(&bad.to_string()),
        Err(RecipeError::UnknownOrigin { .. })
    ));
}

#[test]
fn all_faces_does_not_relax_origin_pattern_or_reference_keys() {
    for reference in [
        json!({"face":{}}),
        json!({"all":true}),
        json!({"face":{"step":"stock","role":"missing"}}),
    ] {
        assert!(Recipe::parse(&source(FORMAT_V5, reference).to_string()).is_err());
    }
}
