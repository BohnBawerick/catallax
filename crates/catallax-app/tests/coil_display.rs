use catallax_app::{NumericalOptions, numerical_sweeps};
use catallax_mesh::{MeshOptions, PreviewOptions};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::json;

const DISPLAY_WORK: u64 = 20_000_000;

#[test]
fn three_turn_coil_display_and_preview_fit_the_target() {
    let recipe: Recipe = serde_json::from_value(json!({
        "format": "catallax-recipe/5", "name": "Three-turn coil display", "steps": [
            {"id": "profile", "op": "sketch", "plane": {
                "rotation": ["1", "1", "0", "0"], "translation": ["0", "0", "0"]
            }, "loops": [[{"kind": "circle", "id": "wire", "centre": ["17", "2.5"], "radius": "1"}]]},
            {"id": "path", "op": "path", "start": ["17", "0", "2.5"], "segments": [
                {"kind": "helix", "id": "coil", "centre": ["0", "0", "2.5"],
                 "axis": ["0", "0", "1"], "angle": "1080", "rise": "9"}
            ]},
            {"id": "coil", "op": "sweep", "profile": "profile", "path": "path",
             "orientation": "rotation_minimising"}
        ], "result": "coil"
    }))
    .unwrap();
    let (body, build_work) = with_operation_budget(100_000_000, || {
        numerical_sweeps::build(
            &recipe,
            NumericalOptions::default(),
            &mut Decider::default(),
        )
    });
    let body = body.unwrap_or_else(|error| panic!("coil build {build_work:?}: {error:?}"));
    let (display, display_work) = with_operation_budget(DISPLAY_WORK, || {
        body.display(&MeshOptions::default(), &mut Decider::default())
    });
    let (preview, preview_work) = with_operation_budget(DISPLAY_WORK, || {
        body.preview(&PreviewOptions::default(), &mut Decider::default())
    });
    let (angular, angular_work) = with_operation_budget(DISPLAY_WORK, || {
        body.preview(
            &PreviewOptions {
                angular_deviation: Some(4.into()),
                ..PreviewOptions::default()
            },
            &mut Decider::default(),
        )
    });
    eprintln!(
        "coil display {display_work:?}, preview {preview_work:?}, angular preview {angular_work:?}"
    );

    let display = display
        .unwrap_or_else(|error| panic!("coil display {display_work:?}: {error:?}"))
        .fitted_mesh;
    let preview =
        preview.unwrap_or_else(|error| panic!("coil preview {preview_work:?}: {error:?}"));
    let angular =
        angular.unwrap_or_else(|error| panic!("coil angular preview {angular_work:?}: {error:?}"));
    assert!(!display_work.exhausted && !preview_work.exhausted && !angular_work.exhausted);
    assert!(!display.triangles.is_empty() && !preview.triangles.is_empty());
    assert!(preview.faces.iter().all(|face| !face.triangles.is_empty()));
    assert_eq!(
        preview, angular,
        "a target above pi keeps the same checked plan"
    );
}
