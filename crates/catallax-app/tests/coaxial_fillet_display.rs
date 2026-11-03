//! Torus rounds between coaxial surfaces of revolution through the application
//! entry at SolidDraft's automatic viewport levels: 15 degrees and
//! `0.002 * diagonal * factor`, halved by the adapter for binary32 rounding,
//! under the default work budget.

use std::{collections::BTreeMap, time::Instant};

use catallax_mesh::{Mesh, MeshOptions, PreviewOptions, preview_mesh_body};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::{Value, json};

const UP: [&str; 4] = ["1", "0", "0", "0"];

fn pose(z: &str) -> Value {
    json!({"rotation": UP, "translation": ["0", "0", z]})
}
fn face(step: &str, role: &str) -> Value {
    json!({"step": step, "role": role})
}
fn side(step: &str, segment: &str) -> Value {
    json!({"step": step, "role": {"side": {"segment": segment}}})
}
fn block(id: &str, z: &str, height: &str) -> Vec<Value> {
    vec![
        json!({"id": format!("{id}-section"), "op": "sketch",
        "plane": {"rotation": UP, "translation": ["0", "0", z]},
        "loops": [[
            {"kind": "line", "id": "s", "to": ["20", "-20"]},
            {"kind": "line", "id": "e", "to": ["20", "20"]},
            {"kind": "line", "id": "n", "to": ["-20", "20"]},
            {"kind": "line", "id": "w", "to": ["-20", "-20"]}
        ]]}),
        json!({"id": id, "op": "extrude", "profile": format!("{id}-section"), "distance": height}),
    ]
}
fn round(mut steps: Vec<Value>, body: &str, radius: &str, edges: Vec<[Value; 2]>) -> Recipe {
    let edges: Vec<_> = edges.into_iter().map(|e| json!({"between": e})).collect();
    steps.push(
        json!({"id": "round", "op": "fillet", "body": body, "radius": radius,
        "edges": edges}),
    );
    Recipe::parse(
        &json!({"format": "catallax-recipe/3", "name": "coaxial round", "steps": steps,
            "result": "round"})
        .to_string(),
    )
    .unwrap()
}
fn cone(base: &str, top: &str) -> Vec<Value> {
    vec![
        json!({"id": "c", "op": "cone", "base_radius": base, "top_radius": top,
        "height": "6", "pose": pose("0")}),
    ]
}
fn revolved(segments: Value) -> Vec<Value> {
    vec![
        json!({"id": "section", "op": "sketch",
            "plane": {"rotation": UP, "translation": ["0", "0", "0"]}, "loops": [segments]}),
        json!({"id": "r", "op": "revolve", "profile": "section",
            "axis": {"point": ["0", "0", "0"], "direction": ["0", "1", "0"]}, "angle": "360"}),
    ]
}

fn cases() -> Vec<(&'static str, Recipe)> {
    let chamfered = vec![
        json!({"id": "c", "op": "cylinder", "radius": "5", "height": "10", "pose": pose("0")}),
        json!({"id": "bevel", "op": "chamfer", "body": "c", "distance": "1",
            "edges": [{"between": [face("c", "end"), side("c", "wall")]}]}),
    ];
    let mut dome = vec![json!({"id": "s", "op": "sphere", "radius": "5", "pose": pose("0")})];
    dome.extend(block("b", "2", "10"));
    dome.push(json!({"id": "j", "op": "boolean", "kind": "intersect", "target": "s", "tool": "b"}));
    let mut seat = block("b", "0", "20");
    seat.push(json!({"id": "s", "op": "sphere", "radius": "5", "pose": pose("23")}));
    seat.push(json!({"id": "j", "op": "boolean", "kind": "union", "target": "b", "tool": "s"}));
    let mut sink = block("b", "0", "20");
    sink.push(
        json!({"id": "h", "op": "hole", "target": "b", "at": [{"id": "p",
        "point": ["0", "0", "20"], "direction": ["0", "0", "-1"]}],
        "diameter": "6", "depth": "through", "tip": "flat",
        "head": {"countersink": {"diameter": "12", "angle": "90"}}}),
    );
    let hole = |part: &str| json!({"step": "h", "role": {"hole": {"at": "p", "part": part}}});
    let groove = revolved(json!([
        {"kind": "line", "id": "floor", "to": ["6", "0"]},
        {"kind": "arc", "id": "groove", "centre": ["10", "0"], "to": ["10", "4"], "turn": "cw"},
        {"kind": "line", "id": "rim", "to": ["10", "6"]},
        {"kind": "line", "id": "top", "to": ["0", "6"]},
        {"kind": "line", "id": "axis", "to": ["0", "0"]}
    ]));
    let rod = revolved(json!([
        {"kind": "line", "id": "base", "to": ["6", "0"]},
        {"kind": "line", "id": "rod", "to": ["6", "8"]},
        {"kind": "arc", "id": "cap", "centre": ["0", "0"], "to": ["0", "10"], "turn": "ccw"},
        {"kind": "line", "id": "axis", "to": ["0", "0"]}
    ]));
    vec![
        (
            "frustum base",
            round(
                cone("5", "3"),
                "c",
                "1",
                vec![[face("c", "start"), side("c", "wall")]],
            ),
        ),
        (
            "narrow top",
            round(
                cone("4", "1"),
                "c",
                "1",
                vec![[face("c", "end"), side("c", "wall")]],
            ),
        ),
        (
            "chamfer both",
            round(
                chamfered,
                "bevel",
                "1/4",
                vec![
                    [face("c", "end"), json!({"step": "bevel"})],
                    [side("c", "wall"), json!({"step": "bevel"})],
                ],
            ),
        ),
        (
            "dome rim",
            round(
                dome,
                "j",
                "1",
                vec![[face("b", "start"), side("s", "wall")]],
            ),
        ),
        (
            "sphere seat",
            round(seat, "j", "1", vec![[face("b", "end"), side("s", "wall")]]),
        ),
        (
            "countersink both",
            round(
                sink,
                "h",
                "1",
                vec![
                    [face("b", "end"), hole("countersink")],
                    [hole("countersink"), hole("wall")],
                ],
            ),
        ),
        (
            "groove",
            round(
                groove,
                "r",
                "1",
                vec![
                    [side("r", "floor"), side("r", "groove")],
                    [side("r", "groove"), side("r", "rim")],
                ],
            ),
        ),
        (
            "rod cap",
            round(rod, "r", "1", vec![[side("r", "rod"), side("r", "cap")]]),
        ),
    ]
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

#[test]
fn coaxial_rounds_stay_closed_at_every_viewport_level() {
    for (kind, recipe) in cases() {
        let default = MeshOptions::default();
        let body = catallax_app::build(&recipe, &default, &mut Decider::default())
            .unwrap_or_else(|e| panic!("{kind}: {e:?}"));
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
