//! Independent material and boundary queries for planar and circular bevels.

use catallax_field::{FaceState, Location, Plan, Workspace};
use catallax_numbers::{BudgetPolicy, Decider, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;
use serde_json::json;

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    })
}
fn point(p: [&str; 3]) -> RPoint3 {
    let [x, y, z] = p.map(|x| Rational::parse(x).unwrap());
    RPoint3::new(x, y, z)
}
fn check(recipe: &Recipe, cases: &[([&str; 3], Location)]) -> Plan {
    let plan = Plan::compile_checked(recipe, &mut decider()).unwrap();
    let mut workspace = Workspace::default();
    for (p, expected) in cases {
        assert_eq!(
            plan.proven_checked(&point(*p), &mut workspace, &mut decider())
                .unwrap()
                .location,
            *expected,
            "at {p:?}"
        );
    }
    plan
}

#[test]
fn p2_has_independent_named_conical_boundaries() {
    let recipe = Recipe::parse(include_str!("../../../parts/P2-collar.json")).unwrap();
    let plan = check(
        &recipe,
        &[
            (["23/2", "1/2", "0"], Location::Boundary),
            (["47/4", "1/2", "0"], Location::Outside),
            (["11", "1/2", "0"], Location::Inside),
            (["23/2", "23/2", "0"], Location::Boundary),
            (["6", "6", "0"], Location::Boundary),
            (["5", "6", "0"], Location::Outside),
        ],
    );
    let mut named = 0;
    for face in plan
        .faces()
        .iter()
        .filter(|f| matches!(f.origin.role, Role::Bevel { .. }))
    {
        let face = plan.select(&face.origin).unwrap();
        if plan
            .face_query_checked(
                face,
                &point(["23/2", "1/2", "0"]),
                &mut Workspace::default(),
                &mut decider(),
            )
            .unwrap()
            .state
            == FaceState::OnFace
        {
            named += 1;
        }
    }
    assert_eq!(named, 1);
}

#[test]
fn all_collar_edges_include_the_inward_facing_cones() {
    let mut value: serde_json::Value =
        serde_json::from_str(include_str!("../../../parts/P2-collar.json")).unwrap();
    value["steps"][2]["edges"] = json!("all");
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let plan = check(
        &recipe,
        &[
            (["13/2", "1/2", "0"], Location::Boundary),
            (["25/4", "1/2", "0"], Location::Outside),
            (["27/4", "1/2", "0"], Location::Inside),
            (["23/2", "1/2", "0"], Location::Boundary),
        ],
    );
    assert_eq!(
        plan.faces()
            .iter()
            .filter(|f| matches!(f.origin.role, Role::Bevel { .. }))
            .count(),
        4
    );
}

#[test]
fn unequal_line_setbacks_follow_input_face_order() {
    let recipe = Recipe::parse(&json!({"format":"catallax-recipe/3", "name":"two setbacks", "steps":[
        {"id":"base", "op":"block", "size":["10","12","8"], "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
        {"id":"bevel", "op":"chamfer", "body":"base", "distance":"2", "other_distance":"3", "edges":[{"between":[
            {"step":"base","role":{"side":{"segment":"x_min"}}}, {"step":"base","role":{"side":{"segment":"y_min"}}}
        ]}]}
    ], "result":"bevel"}).to_string()).unwrap();
    check(
        &recipe,
        &[
            (["1/2", "1", "4"], Location::Outside),
            (["3/2", "1", "4"], Location::Boundary),
            (["1/2", "2", "4"], Location::Inside),
            (["1/2", "2", "9"], Location::Outside),
        ],
    );
}

#[test]
fn concave_boss_bevel_adds_material_only_at_the_selected_junction() {
    let mut value: serde_json::Value =
        serde_json::from_str(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    value["steps"].as_array_mut().unwrap().truncate(4);
    value["steps"][3]["op"] = json!("chamfer");
    value["steps"][3].as_object_mut().unwrap().remove("radius");
    value["steps"][3]["distance"] = json!("2");
    value["result"] = json!("junction");
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    check(
        &recipe,
        &[
            (["29", "15", "11"], Location::Boundary),
            (["57/2", "15", "11"], Location::Inside),
            (["59/2", "15", "11"], Location::Outside),
            (["29", "15", "21"], Location::Outside),
        ],
    );
}

#[test]
fn convex_revolved_sections_keep_axis_interior_winding_and_pose_invariant() {
    for reversed in [false, true] {
        for axis in ["2", "-2"] {
            for placed in [false, true] {
                let mut ends = [["4", "0"], ["0", "6"], ["0", "0"]];
                if reversed {
                    ends.reverse();
                }
                let segments: Vec<_> = ends
                    .iter()
                    .enumerate()
                    .map(|(i, to)| json!({"kind":"line", "id":format!("s{i}"), "to":to}))
                    .collect();
                let recipe = Recipe::parse(
                    &json!({"format":"catallax-recipe/3", "name":"convex cone", "steps":[
                    {"id":"profile", "op":"sketch", "plane":{
                        "rotation":if placed { ["1","1","1","1"] } else { ["1","0","0","0"] },
                        "translation":if placed { ["7","-8","9"] } else { ["0","0","0"] }
                    }, "loops":[segments]},
                    {"id":"solid", "op":"revolve", "profile":"profile", "angle":"360", "axis":{
                        "point":if placed { ["7","-8","9"] } else { ["0","0","0"] },
                        "direction":if placed { ["0","0",axis] } else { ["0",axis,"0"] }
                    }}
                ], "result":"solid"})
                    .to_string(),
                )
                .unwrap();
                let plan = Plan::compile_checked(&recipe, &mut decider()).unwrap();
                for (p, want) in [
                    (["0", "3", "0"], Location::Inside),
                    (["0", "6", "0"], Location::Boundary),
                    (["0", "7", "0"], Location::Outside),
                    (["1", "3", "1"], Location::Inside),
                    (["6/5", "3", "8/5"], Location::Boundary),
                    (["2", "3", "1"], Location::Outside),
                    (["-2", "3", "0"], Location::Boundary),
                ] {
                    let p = point(p);
                    let p = if placed {
                        RPoint3::new(
                            p.z + Rational::from_i64(7),
                            p.x - Rational::from_i64(8),
                            p.y + Rational::from_i64(9),
                        )
                    } else {
                        p
                    };
                    assert_eq!(
                        plan.proven_checked(&p, &mut Workspace::default(), &mut decider())
                            .unwrap()
                            .location,
                        want,
                        "reverse={reversed}, axis={axis}, placed={placed}, point={p:?}"
                    );
                }
            }
        }
    }
}
