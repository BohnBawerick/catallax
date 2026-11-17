//! Paired full/reused display meshing after a bracket extrusion-depth edit.
use catallax_mesh::{Mesh, MeshOptions, mesh_body, mesh_body_reusing};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::{Value, json};
use std::{collections::BTreeMap, time::Instant};

fn closed(mesh: &Mesh) -> bool {
    let mut edges = BTreeMap::<_, (usize, i32)>::new();
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let entry = edges.entry((a.min(b), a.max(b))).or_default();
            entry.0 += 1;
            entry.1 += if a < b { 1 } else { -1 };
        }
    }
    edges
        .values()
        .all(|&(count, winding)| count == 2 && winding == 0)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let repetitions: usize = args.next().map_or(Ok(3), |value| value.parse())?;
    let chord = args.next().unwrap_or_else(|| "0.01".into());
    if repetitions == 0 || args.next().is_some() {
        return Err("usage: mesh_reuse_profile [REPETITIONS] [CHORD_MM]".into());
    }
    let mut input: Value = serde_json::from_str(include_str!(
        "../../catallax-cli/tests/fixtures/soliddraft-gui-bracket.json"
    ))?;
    let before_recipe = Recipe::parse(&input["problem"].to_string())?;
    let step = input["problem"]["steps"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|step| step["id"] == "op_33")
        .unwrap();
    step["distance"] = "14".into();
    let after_recipe = Recipe::parse(&input["problem"].to_string())?;
    let before = catallax_features::evaluate(&before_recipe, &mut Decider::default())?;
    let after = catallax_features::evaluate(&after_recipe, &mut Decider::default())?;
    let options = MeshOptions {
        chord_tolerance: Rational::parse(&chord)?,
        ..Default::default()
    };
    let previous = mesh_body(
        &before.model,
        before.result,
        &options,
        &mut Decider::default(),
    )?;
    let unchanged = mesh_body_reusing(
        &before.model,
        before.result,
        &options,
        &before.model,
        &options,
        &previous.value,
        &mut Decider::default(),
    )?;
    assert_eq!(unchanged, previous);
    let mut pairs = Vec::new();
    // Warm up both paths, then alternate their order. Geometry construction,
    // equality/closure checks and JSON encoding are outside the measured calls.
    for iteration in 0..=repetitions {
        let mut row = json!({"iteration": iteration, "warmup": iteration == 0});
        for reuse in if iteration % 2 == 0 {
            [false, true]
        } else {
            [true, false]
        } {
            let start = Instant::now();
            let (answer, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
                if reuse {
                    mesh_body_reusing(
                        &after.model,
                        after.result,
                        &options,
                        &before.model,
                        &options,
                        &previous.value,
                        &mut Decider::default(),
                    )
                } else {
                    mesh_body(
                        &after.model,
                        after.result,
                        &options,
                        &mut Decider::default(),
                    )
                }
            });
            let elapsed = start.elapsed().as_micros();
            let mesh = answer?.value;
            assert!(closed(&mesh));
            row[if reuse { "reused" } else { "full" }] = json!({
                "mesh_us": elapsed, "work": work.used, "closed": true,
                "faces": mesh.faces.len(), "triangles": mesh.triangles.len(),
                "vertices": mesh.positions.len(),
            });
        }
        pairs.push(row);
    }
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "format": "catallax-mesh-reuse-profile/1",
            "fixture": "soliddraft-gui-bracket.json", "step": "op_33",
            "depth_before_mm": "10", "depth_after_mm": "14", "chord_mm": chord,
            "scope": "display triangulation entry, including edge sampling, reuse matching and closure",
            "unchanged_equals_fresh": true, "pairs": pairs,
        }))?
    );
    Ok(())
}
