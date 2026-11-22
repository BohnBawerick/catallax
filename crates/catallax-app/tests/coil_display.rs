//! Display of the retained three-turn round-wire sweep.

use std::collections::BTreeMap;

use catallax_app::{NumericalOptions, numerical_sweeps};
use catallax_mesh::{MeshOptions, PreviewOptions, Triangle};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::json;

const WORK_TARGET: u64 = 20_000_000;
const REQUEST_ALLOWANCE: u64 = 100_000_000;

fn assert_closed(triangles: &[Triangle]) {
    assert!(!triangles.is_empty());
    let mut incidence = BTreeMap::<_, (u32, i32)>::new();
    for triangle in triangles {
        for i in 0..3 {
            let a = triangle.vertices[i];
            let b = triangle.vertices[(i + 1) % 3];
            assert_ne!(a, b);
            let (count, winding) = incidence.entry((a.min(b), a.max(b))).or_default();
            *count += 1;
            *winding += if a < b { 1 } else { -1 };
        }
    }
    assert!(incidence.values().all(|&uses| uses == (2, 0)));
}

#[test]
fn three_turn_coil_display_and_preview_fit_the_request_allowance() {
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
    let (body, build_work) = with_operation_budget(REQUEST_ALLOWANCE, || {
        numerical_sweeps::build(
            &recipe,
            NumericalOptions::default(),
            &mut Decider::default(),
        )
    });
    let body = body.unwrap_or_else(|error| panic!("coil build {build_work:?}: {error:?}"));
    let options = MeshOptions::default();
    let preview_options = PreviewOptions::from(&options);
    let (display, display_work) = with_operation_budget(REQUEST_ALLOWANCE, || {
        body.display(&options, &mut Decider::default())
    });
    let (preview, preview_work) = with_operation_budget(REQUEST_ALLOWANCE, || {
        body.preview(&preview_options, &mut Decider::default())
    });
    let (angular, angular_work) = with_operation_budget(REQUEST_ALLOWANCE, || {
        body.preview(
            &PreviewOptions {
                angular_deviation: Some(4.into()),
                ..preview_options
            },
            &mut Decider::default(),
        )
    });
    for (name, work) in [
        ("display", display_work),
        ("preview", preview_work),
        ("angular preview", angular_work),
    ] {
        eprintln!(
            "coil {name} used {} work against target {WORK_TARGET} and allowance {REQUEST_ALLOWANCE} (exhausted: {})",
            work.used, work.exhausted,
        );
        assert!(
            !work.exhausted && work.used <= REQUEST_ALLOWANCE,
            "{name}: {work:?}"
        );
    }

    let display = display
        .unwrap_or_else(|error| panic!("coil display {display_work:?}: {error:?}"))
        .fitted_mesh;
    let preview =
        preview.unwrap_or_else(|error| panic!("coil preview {preview_work:?}: {error:?}"));
    let angular =
        angular.unwrap_or_else(|error| panic!("coil angular preview {angular_work:?}: {error:?}"));
    assert!(
        display
            .faces
            .iter()
            .all(|face| !face.triangles.is_empty()
                && face.max_chord_error <= options.chord_tolerance)
    );
    assert!(
        display
            .edges
            .iter()
            .all(|edge| edge.max_chord_error <= options.chord_tolerance)
    );
    assert_closed(&display.triangles);
    assert!(
        display.triangles == preview.triangles && display.positions == preview.positions,
        "preview must retain the checked display facets"
    );
    assert!(
        preview == angular,
        "a target above pi keeps the same checked plan"
    );
}
