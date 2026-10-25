//! Public recipe and consumer regressions for exact algebraic revolutions.

use catallax_features::evaluate;
use catallax_numbers::{Ball, Decider, Rational};
use catallax_recipe::Recipe;

const CASES: &[(&str, &str)] = &[
    (
        "line-30",
        include_str!("fixtures/algebraic-revolutions/line-30.json"),
    ),
    (
        "line-45",
        include_str!("fixtures/algebraic-revolutions/line-45.json"),
    ),
    (
        "line-60",
        include_str!("fixtures/algebraic-revolutions/line-60.json"),
    ),
    (
        "arc-30",
        include_str!("fixtures/algebraic-revolutions/arc-30.json"),
    ),
    (
        "arc-45",
        include_str!("fixtures/algebraic-revolutions/arc-45.json"),
    ),
    (
        "arc-60",
        include_str!("fixtures/algebraic-revolutions/arc-60.json"),
    ),
    (
        "irrational-axis",
        include_str!("fixtures/algebraic-revolutions/irrational-axis.json"),
    ),
    (
        "sphere-radius",
        include_str!("fixtures/algebraic-revolutions/sphere-radius.json"),
    ),
    (
        "torus-minor-radius",
        include_str!("fixtures/algebraic-revolutions/torus-minor-radius.json"),
    ),
    (
        "torus-major-radius",
        include_str!("fixtures/algebraic-revolutions/torus-major-radius.json"),
    ),
];

#[test]
fn public_recipes_keep_valid_topology_origins_and_step() {
    for (name, source) in CASES {
        let recipe = Recipe::parse(source).unwrap();
        let built =
            evaluate(&recipe, &mut Decider::default()).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert!(built.model.check(built.result).is_empty(), "{name}");
        for face in built.model.faces_of_body(built.result) {
            let origins = built.model.face_origins(face);
            assert!(!origins.is_empty(), "{name}");
            assert!(origins.iter().all(|o| o.step.as_str() == "body"), "{name}");
        }
        let mut step = Vec::new();
        let answer = catallax_step::write_step(
            &built.model,
            &[built.result],
            &Default::default(),
            &mut Decider::default(),
            &mut step,
        )
        .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert!(answer.marks.is_empty(), "{name}: {:?}", answer.marks);
        assert!(!step.is_empty(), "{name}");
    }
}

#[test]
fn public_recipes_answer_requested_measurements() {
    for (name, source) in CASES {
        let recipe = Recipe::parse(source).unwrap();
        let built = evaluate(&recipe, &mut Decider::default()).unwrap();
        let answer = catallax_measure::mass_properties(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert!(answer.marks.is_empty(), "{name}: {:?}", answer.marks);
        assert!(answer.value.centroid.is_some(), "{name}");
        assert!(answer.value.inertia.is_some(), "{name}");
    }
}

#[test]
fn public_recipes_display_at_the_existing_accuracy() {
    for (name, source) in CASES {
        let recipe = Recipe::parse(source).unwrap();
        let built = evaluate(&recipe, &mut Decider::default()).unwrap();
        let options = catallax_mesh::MeshOptions::default();
        let answer = catallax_mesh::mesh_body(
            &built.model,
            built.result,
            &options,
            &mut Decider::default(),
        )
        .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert!(!answer.value.triangles.is_empty(), "{name}");
        assert_eq!(
            answer.value.faces.len(),
            built.model.faces_of_body(built.result).len(),
            "{name}"
        );
        assert!(
            answer
                .value
                .faces
                .iter()
                .all(|f| f.max_chord_error <= options.chord_tolerance),
            "{name}"
        );
    }
}

#[test]
fn parameter_edits_preserve_all_physical_face_and_edge_names() {
    for (name, source) in CASES {
        let original = Recipe::parse(source).unwrap();
        let mut edited: serde_json::Value = serde_json::from_str(source).unwrap();
        for edge in edited["steps"][0]["loops"][0].as_array_mut().unwrap() {
            for key in ["to", "centre"] {
                if let Some(point) = edge.get_mut(key) {
                    for coordinate in point.as_array_mut().unwrap() {
                        let value = catallax_numbers::Rational::parse(coordinate.as_str().unwrap())
                            .unwrap();
                        *coordinate = serde_json::Value::String(
                            (value * catallax_numbers::Rational::from_i64(2)).to_string(),
                        );
                    }
                }
            }
            if let Some(radius) = edge.get_mut("radius") {
                let value = catallax_numbers::Rational::parse(radius.as_str().unwrap()).unwrap();
                *radius = serde_json::Value::String(
                    (value * catallax_numbers::Rational::from_i64(2)).to_string(),
                );
            }
        }
        let edited = Recipe::parse(&edited.to_string()).unwrap();
        let before = evaluate(&original, &mut Decider::default()).unwrap();
        let after = evaluate(&edited, &mut Decider::default()).unwrap();
        for face in before.model.faces_of_body(before.result) {
            for origin in before.model.face_origins(face) {
                assert!(
                    !after
                        .model
                        .faces_with_origin(after.result, origin)
                        .is_empty(),
                    "{name}: {origin:?}"
                );
            }
        }
        for edge in before.model.edges_of_body(before.result) {
            let edge_name = before.model.edge_name(before.result, edge).unwrap();
            assert!(
                !after
                    .model
                    .edges_with_name(after.result, &edge_name)
                    .is_empty(),
                "{name}: {edge_name:?}"
            );
        }
    }
}

#[test]
fn angular_trims_keep_points_on_each_side_of_the_end_cap() {
    use catallax_boolean::{Containment, classify_point};
    for angle in [30, 45, 60] {
        let source = CASES
            .iter()
            .find(|(name, _)| *name == format!("line-{angle}"))
            .unwrap()
            .1;
        for reflected in [false, true] {
            let mut recipe: serde_json::Value = serde_json::from_str(source).unwrap();
            if reflected {
                for edge in recipe["steps"][0]["loops"][0].as_array_mut().unwrap() {
                    let x = edge["to"][0].as_str().unwrap();
                    edge["to"][0] = serde_json::Value::String(format!("-{x}"));
                }
            }
            let recipe = Recipe::parse(&recipe.to_string()).unwrap();
            let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
            for (x, z, want) in [
                ("2", "-1/2", Containment::Inside),
                ("1", "-2", Containment::Outside),
                ("2", "1/2", Containment::Outside),
                ("1/2", "-1/4", Containment::Outside),
            ] {
                let direction = Rational::from_i64(if reflected { -1 } else { 1 });
                let point = catallax_shapes::RPoint3::new(
                    Rational::parse(x).unwrap() * direction.clone(),
                    Rational::one(),
                    Rational::parse(z).unwrap() * direction,
                );
                let query = classify_point(
                    &mut built.model,
                    built.result,
                    &point,
                    &mut Decider::default(),
                );
                assert_eq!(
                    query.proven(),
                    Some(&want),
                    "{angle} reflected={reflected} point={point:?}: {query:?}"
                );
            }
        }
    }
}

#[test]
fn irrational_axis_and_major_radius_match_independent_cylinder_and_torus_controls() {
    let b = |n| Ball::from_rational(&Rational::from_i64(n), 256);
    let root_two = b(2).sqrt().unwrap();
    let pi = Ball::pi(256);
    for (name, volume, area) in [
        (
            "irrational-axis",
            (pi.clone() * root_two.clone()).checked_div(&b(2)).unwrap(),
            b(2) * pi.clone() + b(4),
        ),
        (
            "torus-major-radius",
            b(4) * root_two.clone() * pi.square(),
            b(8) * root_two * pi.square(),
        ),
    ] {
        let source = CASES.iter().find(|(case, _)| *case == name).unwrap().1;
        let recipe = Recipe::parse(source).unwrap();
        let built = evaluate(&recipe, &mut Decider::default()).unwrap();
        let answer = catallax_measure::volume_area_centroid(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        assert!(answer.marks.is_empty(), "{name}: {:?}", answer.marks);
        for (actual, expected) in [(&answer.value.volume, volume), (&answer.value.area, area)] {
            assert!(
                actual.lo <= expected.lo() && actual.hi >= expected.hi(),
                "{name}: {actual:?} versus {expected:?}"
            );
        }
    }
}

#[test]
fn trimmed_horn_and_spindle_sheets_keep_their_coverage_refusal() {
    use serde_json::json;
    for (centre, top, bottom, direction) in [
        (["1", "0"], ["1", "1"], ["1", "-1"], ["0", "1", "0"]),
        (["1", "0"], ["1", "2"], ["1", "-2"], ["0", "1", "0"]),
        (["2", "0"], ["3", "1"], ["3", "-1"], ["1", "1", "0"]),
    ] {
        let recipe = Recipe::parse(&json!({
            "format":"catallax-recipe/3", "name":"trimmed non-ring sheet",
            "steps":[
                {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[
                    {"kind":"arc","id":"round","centre":centre,"to":top,"turn":"ccw"},
                    {"kind":"line","id":"chord","to":bottom}
                ]]},
                {"id":"body","op":"revolve","profile":"profile","axis":{"point":["0","0","0"],"direction":direction},"angle":"360"}
            ],"result":"body"
        }).to_string()).unwrap();
        let catallax_features::BuildError::Geometry { refusal, .. } =
            evaluate(&recipe, &mut Decider::default()).unwrap_err()
        else {
            panic!("geometric coverage refusal")
        };
        assert_eq!(refusal.code, "features.unsupported-spindle-torus");
    }
}
