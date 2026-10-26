//! Real application display, picking and suffix rebuild of circular first hits.

use catallax_mesh::MeshOptions;
use catallax_numbers::Decider;
use catallax_recipe::{Recipe, Role};
use serde_json::{Value, json};

const SOURCE: &str = include_str!(
    "../../catallax-features/tests/fixtures/circular-sphere-limits/sphere-entry-next.json"
);

#[test]
fn application_circle_click_and_suffix_edits_keep_selectable_faces_and_rims() {
    for to in [false, true] {
        for exit in [false, true] {
            let mut value: Value = serde_json::from_str(SOURCE).unwrap();
            if exit {
                value["steps"][0]["plane"]["translation"][2] = json!("10");
            }
            if to {
                value["steps"][2]["op"] = json!("extrude_to");
                value["steps"][2]["limit"] =
                    json!({"face":{"step":"target","role":{"side":{"segment":"wall"}}}});
            }
            let recipe = Recipe::parse(&value.to_string()).unwrap();
            let options = MeshOptions::default();
            let before = catallax_app::build(&recipe, &options, &mut Decider::default()).unwrap();
            assert_eq!(before.names.faces.len(), 3);
            assert_eq!(before.names.edges.len(), 2);
            assert!(before.display.marks.is_empty());
            for triangle in &before.display.value.triangles {
                let origins = &before.names.faces[&triangle.face];
                assert_eq!(origins.len(), 1);
                assert_eq!(origins[0].step.as_str(), "cap");
            }
            assert!(
                before
                    .names
                    .faces
                    .values()
                    .flatten()
                    .any(|o| o.role == Role::End)
            );
            assert!(
                before
                    .names
                    .faces
                    .values()
                    .flatten()
                    .any(|o| matches!(&o.role, Role::Side {segment,..} if segment.as_str()=="rim"))
            );
            let prefix =
                catallax_app::build_incremental(&recipe, None, &options, &mut Decider::default())
                    .unwrap();
            // The target is the second step: changing its radius reuses only
            // the sketch and rebuilds the real first-hit operation and display.
            value["steps"][1]["radius"] = json!("6");
            let edited = Recipe::parse(&value.to_string()).unwrap();
            let after = catallax_app::build_incremental(
                &edited,
                Some(prefix.checkpoint),
                &options,
                &mut Decider::default(),
            )
            .unwrap();
            assert_eq!(after.names.faces.len(), 3);
            assert_eq!(after.names.edges.len(), 2);
            assert!(after.display.marks.is_empty());
            let built = after.checkpoint.built();
            for origins in before.names.faces.values() {
                for origin in origins {
                    assert_eq!(built.model.faces_with_origin(built.result, origin).len(), 1);
                }
            }
            for name in before.names.edges.values() {
                assert_eq!(
                    built
                        .model
                        .edges_with_name(built.result, name.as_ref().unwrap())
                        .len(),
                    1
                );
            }
            assert!(before.built.model.check(before.built.result).is_empty());
        }
    }
}

#[test]
fn application_circle_cannot_replace_an_exhausted_caller_allowance() {
    let recipe = Recipe::parse(SOURCE).unwrap();
    let (result, work) = catallax_numbers::with_operation_budget(0, || {
        catallax_app::build(&recipe, &MeshOptions::default(), &mut Decider::default())
    });
    assert_eq!(result.unwrap_err().code(), "numbers.operation-budget");
    assert!(work.exhausted);
    assert_eq!(work.used, 0);
}
