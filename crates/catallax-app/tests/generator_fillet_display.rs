//! Generator strips through the application entry at SolidDraft's automatic
//! viewport levels: 15 degrees and `0.002 * diagonal * factor`,
//! halved by the adapter for binary32 rounding, under the default work budget.

use std::{collections::BTreeMap, time::Instant};

use catallax_mesh::{Mesh, MeshOptions, PreviewOptions, preview_mesh_body};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::json;

fn cylinder(id: &str, x: &str) -> serde_json::Value {
    json!({"id":id,"op":"cylinder","radius":"5","height":"12",
        "pose":{"rotation":["1","0","0","0"],"translation":[x,"0","0"]}})
}
fn wall(step: &str) -> serde_json::Value {
    json!({"step":step,"role":{"side":{"segment":"wall"}}})
}
fn finish(mut steps: Vec<serde_json::Value>, radius: &str, edge: serde_json::Value) -> Recipe {
    steps.push(
        json!({"id":"round","op":"fillet","body":"joined","radius":radius,
        "edges":[{"between":edge}]}),
    );
    Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"generator round",
            "steps":steps,"result":"round"})
        .to_string(),
    )
    .unwrap()
}
fn pair(kind: &str, radius: &str) -> Recipe {
    finish(
        vec![
            cylinder("a", "-3"),
            cylinder("b", "3"),
            json!({"id":"joined","op":"boolean","kind":kind,"target":"a","tool":"b"}),
        ],
        radius,
        json!([wall("a"), wall("b")]),
    )
}
/// A cylinder about the z axis against a block whose north side is y = 3.
fn flat(kind: &str, radius: &str) -> Recipe {
    let block = [
        json!({"id":"section","op":"sketch",
            "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
            "loops":[[{"kind":"line","id":"s","to":["12","-10"]},
                {"kind":"line","id":"e","to":["12","3"]},
                {"kind":"line","id":"n","to":["-12","3"]},
                {"kind":"line","id":"w","to":["-12","-10"]}]]}),
        json!({"id":"block","op":"extrude","profile":"section","distance":"12"}),
    ];
    let (target, tool) = if kind == "subtract" {
        ("block", "a")
    } else {
        ("a", "block")
    };
    let mut steps = vec![cylinder("a", "0")];
    steps.extend(block);
    steps.push(json!({"id":"joined","op":"boolean","kind":kind,"target":target,"tool":tool}));
    finish(
        steps,
        radius,
        json!([wall("a"), {"step":"block","role":{"side":{"segment":"n"}}}]),
    )
}

/// A cylinder against a slab whose top is tilted, on an ellipse.
fn oblique(kind: &str, radius: &str) -> Recipe {
    let tilted = |at: &str, distance: &str| {
        [
            json!({"id":"section","op":"sketch",
                "plane":{"rotation":["24/25","0","7/25","0"],"translation":["0","0",at]},
                "loops":[[{"kind":"line","id":"s","to":["30","-30"]},
                    {"kind":"line","id":"e","to":["30","30"]},
                    {"kind":"line","id":"n","to":["-30","30"]},
                    {"kind":"line","id":"w","to":["-30","-30"]}]]}),
            json!({"id":"cut","op":"extrude","profile":"section","distance":distance}),
        ]
    };
    let (cylinder, slab, role, target, tool, operation) = match kind {
        "subtract" => (
            ("16", "0"),
            tilted("8", "30"),
            "start",
            "a",
            "cut",
            "subtract",
        ),
        "union" => (
            ("22", "-10"),
            tilted("-20", "20"),
            "end",
            "a",
            "cut",
            "union",
        ),
        _ => (
            ("100", "-60"),
            tilted("-20", "20"),
            "end",
            "cut",
            "a",
            "subtract",
        ),
    };
    let mut steps = vec![
        json!({"id":"a","op":"cylinder","radius":"5","height":cylinder.0,
        "pose":{"rotation":["1","0","0","0"],"translation":["0","0",cylinder.1]}}),
    ];
    steps.extend(slab);
    steps.push(json!({"id":"joined","op":"boolean","kind":operation,"target":target,"tool":tool}));
    finish(
        steps,
        radius,
        json!([wall("a"), {"step":"cut","role":role}]),
    )
}

/// A tee of two cylinders with nonparallel axes, off centre.
fn ring(kind: &str, radius: &str) -> Recipe {
    finish(
        vec![
            json!({"id":"wide","op":"cylinder","radius":"10","height":"24",
                "pose":{"rotation":["1","1","0","0"],"translation":["0","20","0"]}}),
            json!({"id":"wall","op":"cylinder","radius":"4","height":"20",
                "pose":{"rotation":["3","-1","0","0"],"translation":["1","0","0"]}}),
            json!({"id":"joined","op":"boolean","kind":kind,"target":"wide","tool":"wall"}),
        ],
        radius,
        json!([wall("wide"), wall("wall")]),
    )
}
/// Two square cylinders whose axes are 8 apart, so they cross in one loop.
fn looped(kind: &str, radius: &str) -> Recipe {
    finish(
        vec![
            json!({"id":"wide","op":"cylinder","radius":"10","height":"24",
                "pose":{"rotation":["1","1","0","0"],"translation":["0","20","0"]}}),
            json!({"id":"wall","op":"cylinder","radius":"4","height":"40",
                "pose":{"rotation":["1","0","0","0"],"translation":["8","8","-20"]}}),
            json!({"id":"joined","op":"boolean","kind":kind,"target":"wide","tool":"wall"}),
        ],
        radius,
        json!([wall("wide"), wall("wall")]),
    )
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
            mesh.faces.iter().all(|f| match &f.max_angular_deviation {
                Some(actual) => actual <= limit,
                // A band's grid certifies position, and it says so.
                None => mesh.unavailable_controls.iter().any(|u| u.face == f.face),
            }),
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

// The trim on a large cylinder used to send its refinement round the waist until
// the face ran out of vertices, with or without a fillet.
#[test]
fn a_branch_on_a_pipe_meshes_with_flat_facets() {
    let tee = Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"tee","steps":[
            {"id":"wide","op":"cylinder","radius":"10","height":"24",
                "pose":{"rotation":["1","1","0","0"],"translation":["0","20","0"]}},
            {"id":"wall","op":"cylinder","radius":"4","height":"20",
                "pose":{"rotation":["3","-1","0","0"],"translation":["1","0","0"]}},
            {"id":"joined","op":"boolean","kind":"union","target":"wide","tool":"wall"}
        ],"result":"joined"})
        .to_string(),
    )
    .unwrap();
    for chord in ["1/8", "1/32", "1/100"] {
        let options = MeshOptions {
            chord_tolerance: Rational::parse(chord).unwrap(),
            angular_deviation: Some(Rational::parse("0.2617993877991494").unwrap()),
            ..MeshOptions::default()
        };
        let body = catallax_app::build(&tee, &options, &mut Decider::default()).unwrap();
        let mesh = &body.display.value;
        assert_closed(chord, mesh, &options);
        eprintln!("tee {chord}: triangles={}", mesh.triangles.len());
        assert!(mesh.triangles.len() < 100_000, "{}", mesh.triangles.len());
        // Slivers stand up off the surface, and their area gives them away.
        let area: f64 = mesh
            .triangles
            .iter()
            .map(|triangle| {
                let [a, b, c] = triangle.vertices.map(|v| mesh.positions[v]);
                let [u, v]: [[f64; 3]; 2] = [
                    std::array::from_fn(|i| b[i] - a[i]),
                    std::array::from_fn(|i| c[i] - a[i]),
                ];
                let n = [
                    u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0],
                ];
                n.iter().map(|x| x * x).sum::<f64>().sqrt() / 2.0
            })
            .sum();
        // The two cylinders' whole walls and three end discs bound the surface.
        let bound = std::f64::consts::TAU * (10.0 * 24.0 + 4.0 * 20.0 + 100.0 + 8.0);
        assert!(area < bound, "{area}");
    }
}

#[test]
fn generator_fillets_stay_closed_at_every_viewport_level() {
    for (kind, recipe) in [
        ("lens", pair("intersect", "5/4")),
        ("notch", pair("union", "11/8")),
        ("crescent", pair("subtract", "1")),
        ("flat", flat("intersect", "7/4")),
        ("irrational flat", flat("intersect", "3/4")),
        ("boss", flat("union", "9/4")),
        ("channel", flat("subtract", "9/16")),
        ("oblique rim", oblique("subtract", "1")),
        ("oblique boss", oblique("union", "1")),
        ("oblique bore", oblique("bore", "1")),
        ("ring boss", ring("union", "1")),
        ("ring bore", ring("subtract", "1")),
    ] {
        assert_viewport_levels(kind, &recipe);
    }
}

// A loop rebuilds its fit at each level, so each case is its own test.
#[test]
fn a_loop_boss_stays_closed_at_every_viewport_level() {
    assert_viewport_levels("loop boss", &looped("union", "1"));
}

#[test]
fn a_loop_notch_stays_closed_at_every_viewport_level() {
    assert_viewport_levels("loop notch", &looped("subtract", "1"));
}

#[test]
fn a_loop_lens_stays_closed_at_every_viewport_level() {
    assert_viewport_levels("loop lens", &looped("intersect", "1"));
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
