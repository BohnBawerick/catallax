//! Fillets rolling over earlier fillets through the application entry at
//! SolidDraft's automatic viewport levels: 15 degrees and
//! `0.002 * diagonal * factor`, halved by the adapter for binary32 rounding,
//! under the default work budget.

use std::{collections::BTreeMap, time::Instant};

use catallax_mesh::{Mesh, MeshOptions, PreviewOptions, preview_mesh_body};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::{Value, json};

fn pose(at: [&str; 3]) -> Value {
    json!({"rotation":["1","0","0","0"],"translation":at})
}
fn side(step: &str, segment: &str) -> Value {
    json!({"step":step,"role":{"side":{"segment":segment}}})
}
fn cap(step: &str, role: &str) -> Value {
    json!({"step":step,"role":role})
}
fn edge(a: Value, b: Value) -> Value {
    json!({"between":[a,b]})
}
fn fillet(id: &str, body: &str, radius: &str, edges: Value) -> Value {
    json!({"id":id,"op":"fillet","body":body,"radius":radius,"edges":edges})
}
fn block(id: &str, size: [&str; 3], at: [&str; 3]) -> Value {
    json!({"id":id,"op":"block","size":size,"pose":pose(at)})
}
const SIDES: [&str; 4] = ["x_min", "y_min", "x_max", "y_max"];
fn uprights(step: &str) -> Value {
    Value::Array(
        (0..4)
            .map(|i| edge(side(step, SIDES[i]), side(step, SIDES[(i + 1) % 4])))
            .collect(),
    )
}
fn recipe(steps: Vec<Value>) -> Recipe {
    let result = steps.last().unwrap()["id"].clone();
    Recipe::parse(
        &json!({"format":"catallax-recipe/4","name":"fillet over fillet",
            "steps":steps,"result":result})
        .to_string(),
    )
    .unwrap()
}
/// A block with its uprights rounded, then one top edge rolling over them.
fn boss(upright: &str, top: &str) -> Recipe {
    recipe(vec![
        block("block", ["20", "20", "10"], ["0", "0", "0"]),
        fillet("upright", "block", upright, uprights("block")),
        fillet(
            "top",
            "upright",
            top,
            json!([edge(cap("block", "end"), side("block", "x_max"))]),
        ),
    ])
}
/// A pocket with its uprights rounded, then its floor.
fn pocket(upright: &str, floor: &str) -> Recipe {
    recipe(vec![
        block("block", ["20", "20", "10"], ["0", "0", "0"]),
        block("tool", ["10", "10", "10"], ["5", "5", "4"]),
        json!({"id":"pocket","op":"boolean","kind":"subtract","target":"block","tool":"tool"}),
        fillet("upright", "pocket", upright, uprights("tool")),
        fillet(
            "floor",
            "upright",
            floor,
            json!([edge(cap("tool", "start"), side("tool", "x_max"))]),
        ),
    ])
}
/// A boss on a plate with its uprights rounded, then its base.
fn based(upright: &str, base: &str) -> Recipe {
    recipe(vec![
        block("plate", ["30", "30", "5"], ["0", "0", "0"]),
        block("boss", ["10", "10", "10"], ["10", "10", "5"]),
        json!({"id":"joined","op":"boolean","kind":"union","target":"plate","tool":"boss"}),
        fillet("upright", "joined", upright, uprights("boss")),
        fillet(
            "base",
            "upright",
            base,
            json!([edge(cap("plate", "end"), side("boss", "x_max"))]),
        ),
    ])
}
/// A block with its top front edge rounded, then every other edge.
fn every(front: &str, rest: &str) -> Recipe {
    recipe(vec![
        block("block", ["20", "20", "10"], ["0", "0", "0"]),
        fillet(
            "front",
            "block",
            front,
            json!([edge(cap("block", "end"), side("block", "y_min"))]),
        ),
        fillet("rest", "front", rest, json!("all")),
    ])
}
/// A block's top loop rounded with mitred corners, then its uprights.
fn mitred(radius: &str) -> Recipe {
    recipe(vec![
        block("block", ["20", "20", "10"], ["0", "0", "0"]),
        fillet(
            "top",
            "block",
            radius,
            Value::Array(
                SIDES
                    .iter()
                    .map(|s| edge(cap("block", "end"), side("block", s)))
                    .collect(),
            ),
        ),
        fillet("upright", "top", radius, uprights("block")),
    ])
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
                "{case}: edge not shared"
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
            mesh.faces
                .iter()
                .all(|f| f.max_angular_deviation.as_ref().is_some_and(|a| a <= limit)),
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

fn assert_viewport_levels(kind: &str, recipe: &Recipe) {
    let default = MeshOptions::default();
    let body = catallax_app::build(recipe, &default, &mut Decider::default()).unwrap();
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
        let body = catallax_app::build(recipe, &options, &mut Decider::default()).unwrap();
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

#[test]
fn fillets_over_fillets_stay_closed_at_every_viewport_level() {
    for (kind, recipe) in [
        ("boss torus", boss("3", "1")),
        ("boss torus sheet", boss("3", "2")),
        ("boss sphere", boss("3", "3")),
        ("pocket torus", pocket("3", "1")),
        ("pocket sphere", pocket("3", "3")),
        ("every edge", every("3", "1")),
        ("every edge sphere", every("2", "2")),
        ("mitred top", mitred("2")),
        ("boss base", based("2", "1")),
    ] {
        assert_viewport_levels(kind, &recipe);
    }
}
