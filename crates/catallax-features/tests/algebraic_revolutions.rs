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

#[test]
fn field_revolutions_compile_independently_and_resolve_origins() {
    for (name, source) in CASES {
        let recipe = Recipe::parse(source).unwrap();
        let plan =
            catallax_field::Plan::compile(&recipe).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        let built = evaluate(&recipe, &mut Decider::default()).unwrap();
        for face in built.model.faces_of_body(built.result) {
            for origin in built.model.face_origins(face) {
                assert!(plan.select(origin).is_some(), "{name}: {origin:?}");
            }
        }
    }
}

#[test]
fn field_revolutions_keep_exact_angular_trims_on_both_profile_sides() {
    use catallax_field::{Location, Plan, Workspace};
    for angle in [30, 45, 60, 120, 210, 225, 300, 330] {
        for reflected in [false, true] {
            let source = CASES.iter().find(|(name, _)| *name == "line-30").unwrap().1;
            let mut value: serde_json::Value = serde_json::from_str(source).unwrap();
            value["steps"][1]["angle"] = angle.to_string().into();
            if reflected {
                for edge in value["steps"][0]["loops"][0].as_array_mut().unwrap() {
                    let x = edge["to"][0].as_str().unwrap();
                    edge["to"][0] = serde_json::Value::String(format!("-{x}"));
                }
            }
            let recipe = Recipe::parse(&value.to_string()).unwrap();
            let plan = Plan::compile(&recipe).unwrap();
            let mut workspace = Workspace::default();
            for (x, z, want) in [
                ("2", "-1/2", Location::Inside),
                (
                    "1",
                    "-2",
                    if angle <= 60 {
                        Location::Outside
                    } else {
                        Location::Inside
                    },
                ),
                ("2", "1/2", Location::Outside),
                (
                    "-2",
                    "1/2",
                    if angle >= 210 {
                        Location::Inside
                    } else {
                        Location::Outside
                    },
                ),
                (
                    "-1",
                    "2",
                    if angle >= 300 {
                        Location::Inside
                    } else {
                        Location::Outside
                    },
                ),
                (
                    "3/2",
                    "-3/2",
                    if angle == 30 {
                        Location::Outside
                    } else if angle == 45 {
                        Location::Boundary
                    } else {
                        Location::Inside
                    },
                ),
            ] {
                let f = Rational::from_i64(if reflected { -1 } else { 1 });
                let point = catallax_shapes::RPoint3::new(
                    Rational::parse(x).unwrap() * f.clone(),
                    Rational::one(),
                    Rational::parse(z).unwrap() * f,
                );
                assert_eq!(
                    plan.proven_checked(&point, &mut workspace, &mut Decider::default())
                        .unwrap()
                        .location,
                    want,
                    "{angle} reflected={reflected}: {point:?}"
                );
                let fast = plan
                    .fast(
                        [&point.x, &point.y, &point.z].map(Rational::to_f64_nearest),
                        &mut workspace,
                    )
                    .unwrap();
                assert!(!fast.true_distance_formula);
                if want == Location::Inside {
                    assert!(
                        fast.value <= 0.0,
                        "{angle} reflected={reflected}: {point:?}: {fast:?}"
                    );
                } else if want == Location::Outside {
                    assert!(
                        fast.value >= 0.0,
                        "{angle} reflected={reflected}: {point:?}: {fast:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn irrational_revolution_frames_and_radii_have_exact_field_membership() {
    use catallax_field::{Location, Plan, Workspace};
    for (name, cases) in [
        (
            "irrational-axis",
            vec![
                (["2", "1", "-1/2"], Location::Inside),
                (["2", "1", "1/2"], Location::Outside),
                (["2", "1", "0"], Location::Boundary),
                (["3/2", "3/2", "-1/2"], Location::Boundary),
            ],
        ),
        (
            "sphere-radius",
            vec![
                (["5/4", "0", "0"], Location::Inside),
                (["3/2", "0", "0"], Location::Outside),
                (["1", "0", "1"], Location::Boundary),
                (["1", "0", "0"], Location::Boundary),
                (["0", "0", "0"], Location::Outside),
            ],
        ),
        (
            "torus-minor-radius",
            vec![
                (["17/4", "0", "0"], Location::Inside),
                (["9/2", "0", "0"], Location::Outside),
                (["4", "1", "0"], Location::Boundary),
                (["4", "0", "0"], Location::Boundary),
            ],
        ),
        (
            "torus-major-radius",
            vec![
                (["4", "0", "0"], Location::Inside),
                (["5", "0", "0"], Location::Boundary),
                (["6", "0", "0"], Location::Outside),
                (["4", "0", "1/2"], Location::Inside),
                (["4", "0", "3"], Location::Outside),
                (["2", "2", "0"], Location::Outside),
            ],
        ),
    ] {
        let recipe =
            Recipe::parse(CASES.iter().find(|(case, _)| *case == name).unwrap().1).unwrap();
        let plan = Plan::compile(&recipe).unwrap();
        let mut workspace = Workspace::default();
        for (point, expected) in cases {
            let [x, y, z] = point.map(|s| Rational::parse(s).unwrap());
            let point = catallax_shapes::RPoint3::new(x, y, z);
            assert_eq!(
                plan.proven_checked(&point, &mut workspace, &mut Decider::default())
                    .unwrap()
                    .location,
                expected,
                "{name}: {point:?}"
            );
        }
    }
}

#[test]
fn native_field_display_keeps_the_algebraic_cap_coverage_refusal() {
    let source = CASES.iter().find(|(name, _)| *name == "line-30").unwrap().1;
    let recipe = Recipe::parse(source).unwrap();
    let plan = catallax_field::Plan::compile(&recipe).unwrap();
    let refusal =
        catallax_field::display::mesh_plan(&plan, &Default::default(), &mut Decider::default())
            .unwrap_err();
    assert_eq!(refusal.code, "field.display-unsupported");
}
