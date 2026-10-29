//! Countersinks, drill points and plain holes beside other cylinders and cones.
use catallax_boolean::{Containment, classify_point};
use catallax_features::evaluate;
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

// The identity, and (x,y,z) -> (x,-z,y) + (7,11,-13) as an exact rigid move.
const PLACEMENTS: [bool; 2] = [false, true];

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn place(moved: bool, [x, y, z]: [&str; 3]) -> [Rational; 3] {
    if moved {
        [q(x) + q("7"), -q(z) + q("11"), q(y) - q("13")]
    } else {
        [q(x), q(y), q(z)]
    }
}
fn direction(moved: bool, [x, y, z]: [&str; 3]) -> Value {
    let [x, y, z] = [q(x), q(y), q(z)];
    let v = if moved { [x, -z, y] } else { [x, y, z] };
    json!(v.map(|c| c.to_string()))
}
fn pose(moved: bool, translation: [&str; 3]) -> Value {
    let rotation = if moved {
        ["1", "1", "0", "0"]
    } else {
        ["1", "0", "0", "0"]
    };
    let t = place(moved, translation);
    json!({"rotation":rotation,"translation":t.map(|c| c.to_string())})
}
fn hole(
    moved: bool,
    target: &str,
    at: &[[&str; 2]],
    depth: Value,
    tip: Value,
    head: Option<Value>,
) -> Value {
    let at: Vec<_> = at
        .iter()
        .enumerate()
        .map(|(i, [x, y])| {
            json!({"id":format!("h{i}"),"point":place(moved, [x, y, "10"]).map(|c| c.to_string()),
                "direction":direction(moved, ["0", "0", "-1"])})
        })
        .collect();
    let mut h = json!({"id":"hole","op":"hole","target":target,"at":at,"diameter":"4","depth":depth,"tip":tip});
    if let Some(head) = head {
        h["head"] = head;
    }
    h
}
fn boss(moved: bool) -> [Value; 2] {
    [
        json!({"id":"boss","op":"cylinder","radius":"5","height":"10","pose":pose(moved, ["15","0","10"])}),
        json!({"id":"body","op":"boolean","kind":"union","target":"plate","tool":"boss"}),
    ]
}

// Builds a 60 x 40 x 10 plate, then the steps, under the application's
// operation allowance. Each probe is checked against the independent field.
fn check(name: &str, moved: bool, mut steps: Vec<Value>, probes: &[([&str; 3], bool)]) {
    let result = steps.last().unwrap()["id"].clone();
    steps.insert(
        0,
        json!({"id":"plate","op":"block","size":["60","40","10"],"pose":pose(moved, ["-30","-20","0"])}),
    );
    let recipe = Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":name,"steps":steps,"result":result})
            .to_string(),
    )
    .unwrap();
    let (built, work) = catallax_numbers::with_operation_budget(100_000_000, || {
        evaluate(&recipe, &mut Decider::default())
    });
    let mut built = built.unwrap_or_else(|e| panic!("{name}, moved {moved}: {e}: {work:?}"));
    assert!(!work.exhausted, "{name}");
    assert!(built.model.check(built.result).is_empty(), "{name}");
    assert!(built.model.body(built.result).marks.is_empty(), "{name}");
    let field = Plan::compile(&recipe).unwrap();
    let mut workspace = Workspace::default();
    for (point, material) in probes {
        let [x, y, z] = place(moved, *point);
        let point = RPoint3::new(x, y, z);
        let expected = if *material {
            (Containment::Inside, Location::Inside)
        } else {
            (Containment::Outside, Location::Outside)
        };
        assert_eq!(
            classify_point(
                &mut built.model,
                built.result,
                &point,
                &mut Decider::default()
            )
            .proven(),
            Some(&expected.0),
            "{name}, moved {moved}, {point:?}"
        );
        assert_eq!(
            field.proven(&point, &mut workspace).location,
            expected.1,
            "{name}, moved {moved}, {point:?}"
        );
    }
}

#[test]
fn two_blind_drill_points_in_one_feature_build_within_the_allowance() {
    // The 118 degree tip cones have irrational slopes. Each hole's entry ring
    // lies far from the other hole's whole double cone.
    for moved in PLACEMENTS {
        for spacing in [["-10", "10"], ["-3", "3"]] {
            let at = [[spacing[0], "0"], [spacing[1], "0"]];
            check(
                "two drill points",
                moved,
                vec![hole(
                    moved,
                    "plate",
                    &at,
                    json!({"blind":"5"}),
                    json!({"cone":"118"}),
                    None,
                )],
                &[
                    ([spacing[0], "0", "7"], false),
                    ([spacing[1], "0", "7"], false),
                    ([spacing[1], "0", "4"], false),
                    ([spacing[1], "0", "3"], true),
                    ([spacing[1], "5/2", "7"], true),
                    (["0", "0", "1"], true),
                ],
            );
        }
    }
}

#[test]
fn countersinks_and_drill_points_build_beside_a_boss() {
    for moved in PLACEMENTS {
        for (tip, depth, head) in [
            (
                json!("flat"),
                json!("through"),
                Some(json!({"countersink":{"diameter":"8","angle":"90"}})),
            ),
            (
                json!("flat"),
                json!("through"),
                Some(json!({"countersink":{"diameter":"8","angle":"82"}})),
            ),
            (json!({"cone":"118"}), json!({"blind":"5"}), None),
        ] {
            let [boss, body] = boss(moved);
            check(
                "beside a boss",
                moved,
                vec![
                    boss,
                    body,
                    hole(moved, "body", &[["-15", "0"]], depth, tip, head),
                ],
                &[
                    (["-15", "0", "6"], false),
                    (["-15", "3", "1"], true),
                    (["15", "0", "15"], true),
                    (["15", "6", "15"], false),
                ],
            );
        }
    }
}

#[test]
fn a_countersink_builds_beside_a_bore() {
    for moved in PLACEMENTS {
        check(
            "beside a bore",
            moved,
            vec![
                json!({"id":"bore","op":"cylinder","radius":"5","height":"10","pose":pose(moved, ["15","0","0"])}),
                json!({"id":"body","op":"boolean","kind":"subtract","target":"plate","tool":"bore"}),
                hole(
                    moved,
                    "body",
                    &[["-15", "0"]],
                    json!("through"),
                    json!("flat"),
                    Some(json!({"countersink":{"diameter":"8","angle":"90"}})),
                ),
            ],
            &[
                (["-15", "3", "19/2"], false),
                (["-15", "3", "1"], true),
                (["15", "0", "5"], false),
                (["15", "6", "5"], true),
            ],
        );
    }
}

#[test]
fn a_plain_hole_builds_beside_a_bevel_cone() {
    for moved in PLACEMENTS {
        let [boss, body] = boss(moved);
        check(
            "beside a bevel",
            moved,
            vec![
                boss,
                body,
                json!({"id":"bevel","op":"chamfer","body":"body","distance":"1","edges":[{"between":[
                    {"step":"boss","role":"end"},{"step":"boss","role":{"side":{"segment":"wall"}}}]}]}),
                hole(
                    moved,
                    "bevel",
                    &[["-15", "0"]],
                    json!("through"),
                    json!("flat"),
                    None,
                ),
            ],
            &[
                (["-15", "0", "5"], false),
                (["-15", "3", "5"], true),
                (["15", "0", "39/2"], true),
                (["15", "19/4", "79/4"], false),
            ],
        );
    }
}
