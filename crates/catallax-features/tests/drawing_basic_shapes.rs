//! Full revolutions and tilted cylinder cuts at the application's drawing controls.

use std::{collections::BTreeMap, time::Instant};

use catallax_features::evaluate;
use catallax_mesh::{Mesh, MeshOptions, mesh_body};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::{Value, json};

fn recipe(steps: Value, result: &str) -> Recipe {
    Recipe::parse(
        &json!({"format":"catallax-recipe/4", "name":"drawing shape",
        "steps":steps, "result":result})
        .to_string(),
    )
    .unwrap()
}

fn revolved(sphere: bool) -> Recipe {
    let profile = if sphere {
        json!([
            {"kind":"arc", "id":"round", "centre":["0","0"], "to":["0","1"], "turn":"ccw"},
            {"kind":"line", "id":"axis", "to":["0","-1"]}
        ])
    } else {
        json!([{"kind":"circle", "id":"tube", "centre":["3","0"], "radius":"1"}])
    };
    recipe(
        json!([
            {"id":"profile", "op":"sketch", "plane":{"rotation":["1","0","0","0"], "translation":["0","0","0"]}, "loops":[profile]},
            {"id":"shape", "op":"revolve", "profile":"profile", "axis":{"point":["0","0","0"], "direction":["0","1","0"]}, "angle":"360"}
        ]),
        "shape",
    )
}

fn tilted_rod(angle: u32) -> Recipe {
    let [left, right] = if angle == 28 {
        [
            "7.126837726645915".to_owned(),
            "2.873162273354085".to_owned(),
        ]
    } else {
        let rise = 4.0 * f64::from(angle).to_radians().tan();
        [(5.0 + rise).to_string(), (5.0 - rise).to_string()]
    };
    recipe(
        json!([
            {"id":"section", "op":"sketch", "plane":{"rotation":["1","0","0","0"], "translation":["0","0","0"]},
             "loops":[[{"kind":"circle", "id":"rim", "centre":["0","0"], "radius":"1"}]]},
            {"id":"stock", "op":"extrude", "profile":"section", "distance":"10"},
            {"id":"plane", "op":"sketch", "plane":{"rotation":["1","1","0","0"], "translation":["0","0","0"]}, "loops":[[
                {"kind":"line", "id":"a", "to":["4",right]},
                {"kind":"line", "id":"b", "to":["4","12"]},
                {"kind":"line", "id":"c", "to":["-4","12"]},
                {"kind":"line", "id":"d", "to":["-4",left]}
            ]]},
            {"id":"cut", "op":"extrude_through", "target":"stock", "profile":"plane", "side":"both"}
        ]),
        "cut",
    )
}

fn app_tilted_rod(angle: u32) -> Recipe {
    let mut value = serde_json::to_value(tilted_rod(angle)).unwrap();
    let steps = value["steps"].as_array_mut().unwrap();
    steps.pop();
    if angle != 28 {
        // Keep the other planes exactly through z=5 too. Adding the two
        // rounded endpoint doubles independently can move that plane.
        let rise: Rational = (4.0 * f64::from(angle).to_radians().tan())
            .to_string()
            .parse()
            .unwrap();
        let five = Rational::from_i64(5);
        steps[2]["loops"][0][0]["to"][1] = json!((&five - &rise).to_string());
        steps[2]["loops"][0][3]["to"][1] = json!((&five + &rise).to_string());
    }
    // SolidDraft's ThroughAll sweeps from the front plane in the flipped
    // direction, rather than cutting in both directions like extrude_through.
    steps.push(json!({"id":"tool", "op":"extrude", "profile":"plane",
        "distance":"2.1039230484541326", "direction":["0","1","0"]}));
    steps.push(json!({"id":"cut", "op":"boolean", "kind":"subtract",
        "target":"stock", "tool":"tool"}));
    Recipe::parse(&value.to_string()).unwrap()
}

fn app_display(angle: u32) {
    let mut dec = Decider::default();
    let built = evaluate(&app_tilted_rod(angle), &mut dec).unwrap();
    let native = MeshOptions::default();
    check(
        &mesh_body(&built.model, built.result, &native, &mut dec)
            .unwrap()
            .value,
        &native,
        7.5 * std::f64::consts::PI,
    );
    // Levels::new(Auto, bbox.diagonal()) requests factors 4, 1 and 1/4
    // of 0.002 * diagonal, all at 15 degrees. The adapter halves the
    // shortest-decimal chord for f32 rounding and uses its lower bound of pi.
    let diagonal = 108.0_f64.sqrt();
    for factor in [4.0, 1.0, 0.25] {
        let chord: Rational = (0.002 * diagonal * factor).to_string().parse().unwrap();
        let options = MeshOptions {
            chord_tolerance: chord / Rational::from_i64(2),
            angular_deviation: Some("314159265/1200000000".parse().unwrap()),
            ..MeshOptions::default()
        };
        let start = Instant::now();
        let (answer, work) = with_operation_budget(20_000_000, || {
            mesh_body(&built.model, built.result, &options, &mut dec)
        });
        eprintln!(
            "app rod {angle} level {factor}: {:?}, {:?}",
            start.elapsed(),
            answer
                .as_ref()
                .map(|m| (m.value.positions.len(), m.value.triangles.len()))
        );
        check(&answer.unwrap().value, &options, 7.5 * std::f64::consts::PI);
        assert!(!work.exhausted, "{work:?}");
        assert!(work.used < 500_000, "{work:?}");
    }
}

#[test]
fn app_rod_cut_at_15_degrees() {
    app_display(15);
}

#[test]
fn app_rod_cut_at_28_degrees() {
    app_display(28);
}

#[test]
fn app_rod_cut_at_45_degrees() {
    app_display(45);
}

#[test]
fn app_rod_cut_at_53_degrees() {
    app_display(53);
}

#[test]
fn app_rod_vertex_cap_refuses_without_surface_refinement() {
    let built = evaluate(&app_tilted_rod(28), &mut Decider::default()).unwrap();
    let options = MeshOptions {
        angular_deviation: Some("314159265/1200000000".parse().unwrap()),
        max_vertices: 10,
        ..MeshOptions::default()
    };
    let (answer, work) = with_operation_budget(50_000, || {
        mesh_body(
            &built.model,
            built.result,
            &options,
            &mut Decider::default(),
        )
    });
    assert_eq!(answer.unwrap_err().code, "mesh.precision-or-budget");
    assert!(!work.exhausted, "{work:?}");
}

fn check(mesh: &Mesh, options: &MeshOptions, volume: f64) {
    let mut uses = BTreeMap::<_, (usize, i32)>::new();
    let mut enclosed = 0.0;
    for triangle in &mesh.triangles {
        let [a, b, c] = triangle.vertices.map(|i| mesh.positions[i]);
        enclosed += (a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))
            / 6.0;
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let edge = uses.entry((a.min(b), a.max(b))).or_default();
            edge.0 += 1;
            edge.1 += if a < b { 1 } else { -1 };
        }
    }
    assert!(uses.values().all(|edge| *edge == (2, 0)));
    assert!(
        (enclosed - volume).abs() < volume * 0.02,
        "{enclosed} versus {volume}"
    );
    assert!(mesh.positions.len() <= options.max_vertices);
    for face in &mesh.faces {
        assert!(face.max_chord_error <= options.chord_tolerance);
        assert_eq!(face.max_angular_deviation, options.angular_deviation);
    }
}

fn drawing(recipe: Recipe, volume: f64) {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let built = evaluate(&recipe, &mut dec).unwrap();
        let display = MeshOptions::default();
        check(
            &mesh_body(&built.model, built.result, &display, &mut dec)
                .unwrap()
                .value,
            &display,
            volume,
        );
        let bbox = catallax_measure::bounding_box(
            &built.model,
            built.result,
            &Default::default(),
            &mut dec,
        )
        .unwrap()
        .value
        .unwrap();
        let diagonal = (0..3)
            .map(|i| (bbox.max[i].to_f64_nearest() - bbox.min[i].to_f64_nearest()).powi(2))
            .sum::<f64>()
            .sqrt()
            .max(1.0);
        let viewport_chord: Rational = (0.002 * diagonal).to_string().parse().unwrap();
        let viewport = MeshOptions {
            chord_tolerance: viewport_chord / Rational::from_i64(2),
            angular_deviation: Some("314159265/1200000000".parse().unwrap()),
            ..MeshOptions::default()
        };
        check(
            &mesh_body(&built.model, built.result, &viewport, &mut dec)
                .unwrap()
                .value,
            &viewport,
            volume,
        );
        // mesh_hlr requests 2e-5 times the proven body diagonal; the adapter halves it.
        let chord: Rational = (2e-5 * diagonal).to_string().parse().unwrap();
        let options = MeshOptions {
            chord_tolerance: chord / Rational::from_i64(2),
            angular_deviation: Some("314159265/9000000000".parse().unwrap()),
            ..MeshOptions::default()
        };
        let start = Instant::now();
        let mesh = mesh_body(&built.model, built.result, &options, &mut dec);
        eprintln!(
            "chord={} mesh {:?}: {:?}",
            options.chord_tolerance,
            start.elapsed(),
            mesh.as_ref()
                .map(|m| (m.value.positions.len(), m.value.triangles.len()))
        );
        check(&mesh.unwrap().value, &options, volume);
    });
    assert!(!work.exhausted, "{work:?}");
}

#[test]
fn full_revolve_sphere_at_drawing_precision() {
    drawing(revolved(true), 4.0 * std::f64::consts::PI / 3.0);
}

#[test]
fn full_revolve_torus_at_drawing_precision() {
    drawing(revolved(false), 6.0 * std::f64::consts::PI.powi(2));
}

#[test]
fn rod_cut_at_28_degrees_at_drawing_precision() {
    drawing(tilted_rod(28), 5.0 * std::f64::consts::PI);
}

#[test]
fn rod_cut_at_53_degrees_at_drawing_precision() {
    drawing(tilted_rod(53), 5.0 * std::f64::consts::PI);
}
