//! Selected-edge material changes and fillet provenance from the independent field.

use catallax_field::{FaceState, Location, Plan, Workspace};
use catallax_numbers::{BudgetPolicy, Decider, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;

fn point(p: [&str; 3]) -> RPoint3 {
    let [x, y, z] = p.map(|s| Rational::parse(s).unwrap());
    RPoint3::new(x, y, z)
}

#[test]
fn p3_convex_removal_and_concave_addition_are_independent_and_named() {
    let recipe = Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    });
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let mut workspace = Workspace::default();
    for (p, want) in [
        (["1/10", "15", "99/10"], Location::Outside),
        (["3/5", "15", "49/5"], Location::Inside),
        (["2/5", "15", "49/5"], Location::Boundary),
        (["147/5", "15", "51/5"], Location::Outside),
        (["142/5", "15", "51/5"], Location::Inside),
        (["144/5", "15", "52/5"], Location::Boundary),
    ] {
        let p = point(p);
        let query = plan.proven_checked(&p, &mut workspace, &mut dec).unwrap();
        assert_eq!(query.location, want, "at {p:?}: {query:?}");
    }
    let junction = plan
        .faces()
        .iter()
        .find(|f| {
            f.origin.step.as_str() == "junction" && matches!(f.origin.role, Role::Blend { .. })
        })
        .unwrap();
    let face = plan.select(&junction.origin).unwrap();
    let query = plan
        .face_query_checked(
            face,
            &point(["144/5", "15", "52/5"]),
            &mut workspace,
            &mut dec,
        )
        .unwrap();
    assert_eq!(query.state, FaceState::OnFace);
}

#[test]
fn convex_ring_keeps_the_unselected_end_sharp() {
    let recipe = Recipe::parse(r#"{
        "format":"catallax-recipe/3","name":"convex cylinder",
        "steps":[
            {"id":"base","op":"cylinder","radius":"8","height":"12",
             "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
            {"id":"round","op":"fillet","body":"base","radius":"2","edges":[
                {"between":[{"step":"base","role":"end"},{"step":"base","role":{"side":{"segment":"wall"}}}]}
            ]}
        ],"result":"round"
    }"#).unwrap();
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    });
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let mut workspace = Workspace::default();
    for (p, expected) in [
        (["38/5", "0", "58/5"], Location::Outside),
        (["36/5", "0", "58/5"], Location::Boundary),
        (["7", "0", "11"], Location::Inside),
        (["8", "0", "0"], Location::Boundary),
        (["39/5", "0", "1/5"], Location::Inside),
    ] {
        assert_eq!(
            plan.proven_checked(&point(p), &mut workspace, &mut dec)
                .unwrap()
                .location,
            expected,
            "at {p:?}"
        );
    }
}

#[test]
fn reflected_blend_keeps_its_named_boundary_and_material_sector() {
    let recipe = Recipe::parse(r#"{
        "format":"catallax-recipe/3","name":"placed fillet",
        "steps":[
            {"id":"base","op":"cylinder","radius":"8","height":"12",
             "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
            {"id":"round","op":"fillet","body":"base","radius":"2","edges":[
                {"between":[{"step":"base","role":"end"},{"step":"base","role":{"side":{"segment":"wall"}}}]}
            ]},
            {"id":"copies","op":"pattern","body":"round",
             "layout":{"linear":{"vector":["30","0","0"],"count":2}}},
            {"id":"reflected","op":"mirror","body":"copies",
             "plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}
        ],"result":"reflected"
    }"#).unwrap();
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    });
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let mut workspace = Workspace::default();
    for (p, expected) in [
        (["-37", "0", "11"], Location::Inside),
        (["-188/5", "0", "58/5"], Location::Outside),
        (["-186/5", "0", "58/5"], Location::Boundary),
        (["7", "0", "11"], Location::Inside),
        (["37", "0", "11"], Location::Outside),
    ] {
        assert_eq!(
            plan.proven_checked(&point(p), &mut workspace, &mut dec)
                .unwrap()
                .location,
            expected,
            "at {p:?}"
        );
    }
    let mut named_boundary = 0;
    for face in plan.faces().iter().filter(|face| {
        face.origin.step.as_str() == "reflected"
            && matches!(face.surface, catallax_shapes::Surface::Torus(_))
    }) {
        let selected = plan.select(&face.origin).unwrap();
        if plan
            .face_query_checked(
                selected,
                &point(["-186/5", "0", "58/5"]),
                &mut workspace,
                &mut dec,
            )
            .unwrap()
            .state
            == FaceState::OnFace
        {
            named_boundary += 1;
        }
    }
    assert_eq!(named_boundary, 1);
}

#[test]
fn trihedral_corner_removes_the_cube_beyond_its_ball_and_names_the_sector() {
    let recipe = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/trihedral-corner.json"
    ))
    .unwrap();
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    });
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let mut workspace = Workspace::default();
    // The ball is centred at 8,8,8 with radius 2. The 93/10 diagonal point
    // lies in all three strip cylinders but outside the ball.
    for (p, expected) in [
        (["91/10", "91/10", "91/10"], Location::Inside),
        (["93/10", "93/10", "93/10"], Location::Outside),
        (["93/10", "93/10", "5"], Location::Inside),
        (["8", "8", "6"], Location::Inside),
        (["26/3", "28/3", "28/3"], Location::Boundary),
        (["10", "10", "10"], Location::Outside),
    ] {
        assert_eq!(
            plan.proven_checked(&point(p), &mut workspace, &mut dec)
                .unwrap()
                .location,
            expected,
            "at {p:?}"
        );
    }
    let corners: Vec<_> = plan
        .faces()
        .iter()
        .filter(|f| matches!(f.origin.role, Role::Corner { .. }))
        .collect();
    let [corner] = corners.as_slice() else {
        panic!("one corner origin expected: {corners:?}");
    };
    let Role::Corner { vertex } = &corner.origin.role else {
        unreachable!()
    };
    assert_eq!(vertex.faces().len(), 3);
    let face = plan.select(&corner.origin).unwrap();
    for (p, expected) in [
        (["26/3", "28/3", "28/3"], FaceState::OnFace),
        (["8", "8", "6"], FaceState::NotOnFace),
    ] {
        assert_eq!(
            plan.face_query_checked(face, &point(p), &mut workspace, &mut dec)
                .unwrap()
                .state,
            expected,
            "at {p:?}"
        );
    }
}

#[test]
fn smooth_four_face_corner_uses_the_full_radius_sheet_without_rounding_other_edges() {
    let row: serde_json::Value = include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<serde_json::Value>(line).unwrap())
        .find(|row| row["id"] == "vertex-4-two")
        .unwrap();
    for (wall, shift, moving, extra) in [
        ("right", [0, 0, 0], false, false),
        ("left", [3, -2, 7], false, false),
        ("right", [0, 0, 0], true, false),
        ("right", [0, 0, 0], false, true),
    ] {
        let mut value = row["recipe"].clone();
        for step in value["steps"].as_array_mut().unwrap() {
            if moving {
                if step["id"] == "a" {
                    step["distance"] = serde_json::json!("24");
                }
                if step["id"] == "bs" {
                    step["loops"][0][0]["to"] = serde_json::json!(["-14", "0"]);
                    step["loops"][0][2]["to"] = serde_json::json!(["-2", "10"]);
                }
            }
            if step["op"] == "fillet" {
                step["corner"] = serde_json::json!("smooth");
                for edge in step["edges"].as_array_mut().unwrap() {
                    for face in edge["between"].as_array_mut().unwrap() {
                        if face["step"] == "b" {
                            face["role"]["side"]["segment"] = serde_json::json!(wall);
                        }
                    }
                }
                if extra {
                    step["edges"]
                        .as_array_mut()
                        .unwrap()
                        .push(serde_json::json!({
                            "between": [
                                {"step": "b", "role": {"side": {"segment": "left"}}},
                                {"step": "a", "role": {"side": {"segment": "base"}}}
                            ]
                        }));
                }
            }
            if step["op"] == "sketch" {
                for (coordinate, delta) in step["plane"]["translation"]
                    .as_array_mut()
                    .unwrap()
                    .iter_mut()
                    .zip(shift)
                {
                    let old = coordinate.as_str().unwrap().parse::<i64>().unwrap();
                    *coordinate = serde_json::json!((old + delta).to_string());
                }
            }
        }
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let mut dec = Decider::new(BudgetPolicy {
            max_work: 20_000_000,
            ..Default::default()
        });
        let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
        let mut storage = Workspace::default();
        let witnesses: &[([&str; 3], Location)] = if moving {
            &[
                (["0", "-11/5", "1959/200"], Location::Outside),
                (["0", "-11/5", "979/100"], Location::Inside),
                (["5", "4", "5"], Location::Outside),
                (["0", "-2", "5"], Location::Inside),
                (["0", "-2", "0"], Location::Boundary),
            ]
        } else if extra {
            &[
                (["0", "1/100", "99/10"], Location::Outside),
                (["0", "-1/5", "1959/200"], Location::Outside),
                (["0", "-1/5", "979/100"], Location::Inside),
                (["0", "-10", "0"], Location::Outside),
                (["0", "-9", "1"], Location::Boundary),
                (["5", "-5", "5"], Location::Boundary),
                (["0", "0", "5"], Location::Inside),
                (["0", "0", "0"], Location::Boundary),
            ]
        } else {
            &[
                (["0", "1/100", "99/10"], Location::Outside),
                (["0", "-1/100", "99/10"], Location::Outside),
                // The sphere/cylinder cell partition is not a material boundary.
                (["0", "0", "99/10"], Location::Outside),
                (["0", "0", "39/4"], Location::Inside),
                (["0", "-1/5", "1959/200"], Location::Outside),
                (["0", "-1/5", "979/100"], Location::Inside),
                (["1/3", "-1/3", "29/3"], Location::Outside),
                (["0", "0", "5"], Location::Inside),
                (["5", "-5", "5"], Location::Boundary),
                (["5", "5", "5"], Location::Outside),
                (["0", "0", "0"], Location::Boundary),
                (["0", "0", "-1/10"], Location::Outside),
            ]
        };
        for &(coordinates, expected) in witnesses {
            let [x, y, z] = coordinates.map(|v| Rational::parse(v).unwrap());
            let mirror = Rational::from_i64(if wall == "right" { 1 } else { -1 });
            let p = RPoint3::new(
                x + Rational::from_i64(shift[0]),
                y * mirror + Rational::from_i64(shift[1]),
                z + Rational::from_i64(shift[2]),
            );
            assert_eq!(
                plan.proven_checked(&p, &mut storage, &mut dec)
                    .unwrap()
                    .location,
                expected,
                "at {p:?}"
            );
        }
    }
}

#[test]
fn four_contact_field_uses_internal_sheets_without_exposing_cell_caps() {
    let row: serde_json::Value = include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<serde_json::Value>(line).unwrap())
        .find(|row| row["id"] == "vertex-4-all-unequal-slopes")
        .unwrap();
    for (mode, rotated) in [
        ("primary", false),
        ("pythagorean", false),
        ("pythagorean", true),
        ("opposite-diagonal", false),
        ("concurrent", false),
    ] {
        let mut value = row["recipe"].clone();
        for step in value["steps"].as_array_mut().unwrap() {
            if mode == "pythagorean" {
                if step["id"] == "as" {
                    step["loops"][0][0]["to"] = serde_json::json!(["-3", "0"]);
                    step["loops"][0][1]["to"] = serde_json::json!(["3", "0"]);
                    step["loops"][0][2]["to"] = serde_json::json!(["0", "4"]);
                    step["plane"]["translation"] = serde_json::json!(["0", "3", "0"]);
                }
                if step["id"] == "bs" {
                    step["loops"][0][0]["to"] = serde_json::json!(["-5/3", "0"]);
                    step["loops"][0][1]["to"] = serde_json::json!(["3", "0"]);
                    step["loops"][0][2]["to"] = serde_json::json!(["0", "4"]);
                    step["plane"]["translation"] = serde_json::json!(["-3", "0", "0"]);
                }
                if step["op"] == "extrude" {
                    step["distance"] = serde_json::json!("6");
                }
            }
            if mode == "opposite-diagonal" {
                if step["id"] == "bs" {
                    step["loops"][0][0]["to"][0] = serde_json::json!("-15");
                }
                if step["id"] == "a" {
                    step["distance"] = serde_json::json!("40");
                }
            }
            if mode == "concurrent" && step["id"] == "bs" {
                step["loops"][0][0]["to"][0] = serde_json::json!("-10");
            }
            if rotated && step["op"] == "sketch" {
                if step["id"] == "as" {
                    step["plane"]["rotation"] = serde_json::json!(["1", "1", "1", "-1"]);
                } else {
                    step["plane"]["rotation"] = serde_json::json!(["0", "1", "1", "0"]);
                    step["plane"]["translation"] = serde_json::json!(["0", "0", "3"]);
                }
            }
        }
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let mut dec = Decider::new(BudgetPolicy {
            max_work: 20_000_000,
            ..Default::default()
        });
        let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
        let corners: Vec<_> = plan.faces().iter().filter(|face| {
            matches!(&face.origin.role, Role::Corner { vertex } if vertex.faces().len() == 4)
        }).collect();
        assert_eq!(corners.len(), if mode == "concurrent" { 1 } else { 3 });
        assert!(corners.iter().all(|face| face.origin == corners[0].origin));
        let queries: &[([&str; 3], Location)] = match mode {
            "pythagorean" => &[
                // Each plane x=+/-1/8 is shared by a sphere and the internal
                // cylinder. It is not a material boundary above or below r.
                (["1/8", "1/8", "18/5"], Location::Outside),
                (["1/8", "1/8", "69/20"], Location::Inside),
                (["1/8", "1/8", "7/2"], Location::Boundary),
                (["-1/8", "1/8", "18/5"], Location::Outside),
                (["-1/8", "1/8", "69/20"], Location::Inside),
                (["-1/8", "1/8", "7/2"], Location::Boundary),
                (["0", "0", "1"], Location::Inside),
                (["0", "0", "0"], Location::Boundary),
            ],
            "opposite-diagonal" => &[
                (["0", "-0.08", "9.79"], Location::Inside),
                (["0", "-0.08", "9.80"], Location::Outside),
                (["0", "0", "5"], Location::Inside),
                (["0", "0", "0"], Location::Boundary),
            ],
            "concurrent" => &[
                (["0", "0", "9.7"], Location::Inside),
                (["0", "0", "9.9"], Location::Outside),
                (["0", "0", "5"], Location::Inside),
                (["0", "0", "0"], Location::Boundary),
            ],
            _ => &[
                (["0", "0.14", "9.65"], Location::Inside),
                (["0", "0.14", "9.66"], Location::Outside),
                (["0.31", "0.09", "9.61"], Location::Inside),
                (["0.31", "0.09", "9.63"], Location::Outside),
                (["-0.31", "0.09", "9.61"], Location::Inside),
                (["-0.31", "0.09", "9.63"], Location::Outside),
                (["0", "0", "5"], Location::Inside),
                (["0", "0", "0"], Location::Boundary),
                (["5", "0", "5"], Location::Boundary),
                (["-5", "0", "5"], Location::Boundary),
                (["0", "-2.5", "5"], Location::Boundary),
            ],
        };
        let mut workspace = Workspace::default();
        for &(coordinates, expected) in queries {
            let [x, y, z] = coordinates.map(|coordinate| Rational::parse(coordinate).unwrap());
            let p = if rotated {
                RPoint3::new(z, y, -x)
            } else {
                RPoint3::new(x, y, z)
            };
            assert_eq!(
                plan.proven_checked(&p, &mut workspace, &mut dec)
                    .unwrap()
                    .location,
                expected,
                "{mode} at {p:?}",
            );
        }
    }
}
