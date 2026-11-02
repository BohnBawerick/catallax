//! Curved profiles on straight paths use shared, bounded display and STEP geometry.

use std::collections::BTreeMap;

use catallax_mesh::MeshOptions;
use catallax_numbers::{Ball, Decider, Rational};
use catallax_recipe::{Operation, Recipe, Role};
use serde_json::json;

fn recipe(loops: serde_json::Value, end: [&str; 3]) -> Recipe {
    let mut value: serde_json::Value = serde_json::from_str(include_str!(
        "../../catallax-features/tests/fixtures/straight-sweep.json"
    ))
    .unwrap();
    value["steps"][0]["loops"] = loops;
    value["steps"][1]["segments"][0]["to"] = json!(end);
    Recipe::parse(&value.to_string()).unwrap()
}
fn b(n: i64) -> Ball {
    Ball::from_rational(&Rational::from_i64(n), 192)
}

#[test]
fn straight_sweeps_export_and_display_at_kernel_and_soliddraft_settings() {
    let disk = json!([[{"kind":"circle","id":"outer","centre":["0","0"],"radius":"2"}]]);
    let ring = json!([
        [{"kind":"circle","id":"hole","centre":["0","0"],"radius":"1"}],
        [{"kind":"circle","id":"outer","centre":["0","0"],"radius":"2"}]
    ]);
    let slot = json!([[
        {"kind":"line","id":"bottom","to":["4","-1"]},
        {"kind":"arc","id":"right","centre":["4","0"],"to":["4","1"],"turn":"ccw"},
        {"kind":"line","id":"top","to":["0","1"]},
        {"kind":"arc","id":"left","centre":["0","0"],"to":["0","-1"],"turn":"ccw"}
    ]]);
    for (loops, area) in [
        (disk, b(4) * Ball::pi(192)),
        (ring, b(3) * Ball::pi(192)),
        (slot, b(8) + Ball::pi(192)),
    ] {
        for end in [["100", "200", "303"], ["102", "201", "297"]] {
            let source = recipe(loops.clone(), end);
            let mut step_bytes = None;
            for options in [
                MeshOptions::default(),
                MeshOptions {
                    chord_tolerance: Rational::from_ratio(1, 50).unwrap(),
                    angular_deviation: Some(Rational::from_ratio(1, 2).unwrap()),
                    ..Default::default()
                },
            ] {
                let app = catallax_app::build(&source, &options, &mut Decider::default()).unwrap();
                assert!(app.display.marks.is_empty());
                let mesh = &app.display.value;
                assert!(mesh.unavailable_controls.is_empty());
                let mut incidence = BTreeMap::<_, (usize, i32)>::new();
                for triangle in &mesh.triangles {
                    for i in 0..3 {
                        let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                        let uses = incidence.entry((a.min(b), a.max(b))).or_default();
                        uses.0 += 1;
                        uses.1 += if a < b { 1 } else { -1 };
                    }
                }
                assert!(incidence.values().all(|&uses| uses == (2, 0)));
                assert_eq!(mesh.faces.len(), app.names.faces.len());
                for face in &mesh.faces {
                    assert!(face.max_chord_error <= options.chord_tolerance);
                    assert_eq!(face.max_angular_deviation, options.angular_deviation);
                    assert!(!face.triangles.is_empty());
                    for origin in &app.names.faces[&face.face] {
                        if let Role::Side { path, .. } = &origin.role {
                            assert_eq!(path.as_ref().unwrap().as_str(), "travel");
                        }
                    }
                }
                assert!(app.names.edges.values().all(Option::is_some));
                let mass = app
                    .mass_properties(&Default::default(), &mut Decider::default())
                    .unwrap();
                let expected = &area * &b(3);
                assert!(mass.marks.is_empty());
                assert!(
                    mass.value.volume.lo <= expected.lo() && mass.value.volume.hi >= expected.hi()
                );
                assert!(mass.value.volume.width() <= Rational::parse("1e-9").unwrap());
                let mut bytes = Vec::new();
                let options = catallax_step::StepOptions::default();
                let report = catallax_step::write_step(
                    &app.built.model,
                    &[app.built.result],
                    &options,
                    &mut Decider::default(),
                    &mut bytes,
                )
                .unwrap();
                assert!(report.marks.is_empty());
                assert!(report.value.max_deviation <= options.curve_tolerance);
                assert!(bytes.starts_with(b"ISO-10303-21;"));
                if let Some(prior) = &step_bytes {
                    assert_eq!(prior, &bytes);
                }
                step_bytes = Some(bytes);
            }
        }
    }
}

#[test]
fn changing_path_location_does_not_move_the_carried_profile() {
    let source = recipe(
        json!([[{"kind":"circle","id":"outer","centre":["0","0"],"radius":"2"}]]),
        ["102", "201", "303"],
    );
    let original =
        catallax_app::build(&source, &Default::default(), &mut Decider::default()).unwrap();
    let mut translated = source;
    let Operation::Path { start, segments } = &mut translated.steps[1].op else {
        unreachable!()
    };
    start.0 = [-5, -7, -11].map(Into::into);
    let catallax_recipe::PathSegment::Line { to, .. } = &mut segments[0] else {
        unreachable!()
    };
    to.0 = [-3, -6, -8].map(Into::into);
    let moved =
        catallax_app::build(&translated, &Default::default(), &mut Decider::default()).unwrap();
    assert_eq!(original.display.value, moved.display.value);
    assert_eq!(original.display.marks, moved.display.marks);
    assert_eq!(original.names, moved.names);
}
