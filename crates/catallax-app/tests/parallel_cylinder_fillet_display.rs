//! Parallel-cylinder strips through the application entry at SolidDraft's
//! automatic viewport levels: 15 degrees and `0.002 * diagonal * factor`,
//! halved by the adapter for binary32 rounding, under the default work budget.

use std::{collections::BTreeMap, time::Instant};

use catallax_mesh::{Mesh, MeshOptions, PreviewOptions, preview_mesh_body};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::json;

fn recipe(kind: &str, radius: &str) -> Recipe {
    let cylinder = |id: &str, x: &str| {
        json!({"id":id,"op":"cylinder","radius":"5","height":"12",
            "pose":{"rotation":["1","0","0","0"],"translation":[x,"0","0"]}})
    };
    Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"parallel cylinder round","steps":[
            cylinder("a", "-3"),
            cylinder("b", "3"),
            {"id":"joined","op":"boolean","kind":kind,"target":"a","tool":"b"},
            {"id":"round","op":"fillet","body":"joined","radius":radius,"edges":[
                {"between":[{"step":"a","role":{"side":{"segment":"wall"}}},
                    {"step":"b","role":{"side":{"segment":"wall"}}}]}]}
        ],"result":"round"})
        .to_string(),
    )
    .unwrap()
}

fn assert_closed(case: &str, mesh: &Mesh, options: &MeshOptions) {
    let mut incidence = BTreeMap::<_, (usize, i32)>::new();
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let uses = incidence.entry((a.min(b), a.max(b))).or_default();
            uses.0 += 1;
            uses.1 += if a < b { 1 } else { -1 };
        }
    }
    assert!(
        incidence.values().all(|v| *v == (2, 0)),
        "{case}: open or misoriented mesh"
    );
    for edge in &mesh.edges {
        for pair in edge.vertices.windows(2) {
            let key = (pair[0].min(pair[1]), pair[0].max(pair[1]));
            assert_eq!(
                incidence.get(&key),
                Some(&(2, 0)),
                "{case}: rail not shared"
            );
        }
    }
    assert!(
        mesh.faces
            .iter()
            .all(|f| f.max_chord_error <= options.chord_tolerance),
        "{case}: chord bound"
    );
    if let Some(limit) = &options.angular_deviation {
        assert!(
            mesh.faces.iter().all(|f| f
                .max_angular_deviation
                .as_ref()
                .is_some_and(|actual| actual <= limit)),
            "{case}: angular bound"
        );
    }
}

fn diagonal(mesh: &Mesh) -> f64 {
    let extent = |axis: usize| {
        let values = mesh.positions.iter().map(|p| p[axis]);
        values.clone().fold(f64::MIN, f64::max) - values.fold(f64::MAX, f64::min)
    };
    (0..3).map(|axis| extent(axis).powi(2)).sum::<f64>().sqrt()
}

#[test]
fn parallel_cylinder_fillets_stay_closed_at_every_viewport_level() {
    for (kind, radius) in [("intersect", "5/4"), ("union", "11/8"), ("subtract", "1")] {
        let recipe = recipe(kind, radius);
        let default = MeshOptions::default();
        let body = catallax_app::build(&recipe, &default, &mut Decider::default()).unwrap();
        assert_closed(kind, &body.display.value, &default);
        let diagonal = diagonal(&body.display.value);
        for factor in [4.0, 1.0, 0.25] {
            let chord = Rational::parse(&format!("{:.12}", 0.002 * diagonal * factor)).unwrap();
            let options = MeshOptions {
                chord_tolerance: chord / Rational::from_i64(2),
                angular_deviation: Some(Rational::parse("0.2617993877991494").unwrap()),
                ..MeshOptions::default()
            };
            let case = format!("{kind} factor {factor}");
            let started = Instant::now();
            let body = catallax_app::build(&recipe, &options, &mut Decider::default()).unwrap();
            let display_ms = started.elapsed().as_secs_f64() * 1000.0;
            assert_closed(&case, &body.display.value, &options);
            let (preview, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
                preview_mesh_body(
                    &body.built.model,
                    body.built.result,
                    &PreviewOptions::from(&options),
                    &mut Decider::default(),
                )
            });
            assert!(!work.exhausted, "{case} preview: {work:?}");
            let preview = preview.unwrap();
            assert!(preview.faces.iter().all(|f| !f.triangles.is_empty()));
            eprintln!(
                "{case}: build_and_display_ms={display_ms:.3} triangles={}",
                body.display.value.triangles.len()
            );
        }
    }
}
