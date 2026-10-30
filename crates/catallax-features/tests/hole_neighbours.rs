//! Countersinks, drill points and plain holes beside other cylinders and cones.
use catallax_boolean::{Containment, classify_point};
use catallax_features::evaluate;
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Ball, Decider, Rational, with_operation_budget};
use catallax_recipe::{Origin, Recipe, Role};
use catallax_shapes::{RPoint3, Surface};
use serde_json::{Value, json};

// Identity, exact axis swap plus translation, and the app's corner placement.
const PLACEMENTS: [u8; 3] = [0, 1, 2];

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn place(moved: u8, [x, y, z]: [&str; 3]) -> [Rational; 3] {
    match moved {
        0 => [q(x), q(y), q(z)],
        1 => [q(x) + q("7"), -q(z) + q("11"), q(y) - q("13")],
        2 => [q(x) + q("30"), q(y) + q("20"), q(z)],
        _ => unreachable!(),
    }
}
fn direction(moved: u8, [x, y, z]: [&str; 3]) -> Value {
    let [x, y, z] = [q(x), q(y), q(z)];
    let v = if moved == 1 { [x, -z, y] } else { [x, y, z] };
    json!(v.map(|c| c.to_string()))
}
fn pose(moved: u8, translation: [&str; 3]) -> Value {
    let rotation = if moved == 1 {
        ["1", "1", "0", "0"]
    } else {
        ["1", "0", "0", "0"]
    };
    let t = place(moved, translation);
    json!({"rotation":rotation,"translation":t.map(|c| c.to_string())})
}
fn hole(
    moved: u8,
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
fn boss(moved: u8) -> [Value; 2] {
    [
        json!({"id":"boss","op":"cylinder","radius":"5","height":"10","pose":pose(moved, ["15","0","10"])}),
        json!({"id":"body","op":"boolean","kind":"union","target":"plate","tool":"boss"}),
    ]
}

// Builds a 60 x 40 x 10 plate, then the steps, under the application's
// operation allowance. Each probe is checked against the independent field.
fn check(name: &str, moved: u8, mut steps: Vec<Value>, probes: &[([&str; 3], bool)]) {
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

// Quaternion (1, 2, 3, 4) as an exact matrix acting on columns, and its offset.
const COLLAR_ROTATION: [[&str; 3]; 3] = [
    ["-2/3", "2/15", "11/15"],
    ["2/3", "-1/3", "2/3"],
    ["1/3", "14/15", "2/15"],
];
const COLLAR_OFFSET: [&str; 3] = ["7", "11", "-13"];

fn collar_vector(moved: bool, v: [&str; 3]) -> [Rational; 3] {
    let v = v.map(q);
    if !moved {
        return v;
    }
    COLLAR_ROTATION
        .map(|row| (0..3).fold(Rational::zero(), |sum, i| sum + q(row[i]) * v[i].clone()))
}
fn collar_point(moved: bool, p: [&str; 3]) -> [Rational; 3] {
    let [x, y, z] = collar_vector(moved, p);
    if !moved {
        return [x, y, z];
    }
    let [a, b, c] = COLLAR_OFFSET.map(q);
    [x + a, y + b, z + c]
}
fn strings(v: [Rational; 3]) -> Value {
    json!(v.map(|c| c.to_string()))
}

// P2's turned collar with one hole of the given centre and direction.
fn collar(moved: bool, point: [&str; 3], direction: [&str; 3]) -> Recipe {
    let plane = if moved {
        json!({"rotation":["1","2","3","4"],"translation":COLLAR_OFFSET})
    } else {
        json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]})
    };
    let side = |segment: &str| json!({"step":"ring","role":{"side":{"segment":segment}}});
    Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"collar plain hole","steps":[
            {"id":"section","op":"sketch","plane":plane,"loops":[[
                {"kind":"line","id":"bottom","to":["12","0"]},
                {"kind":"line","id":"right","to":["12","12"]},
                {"kind":"line","id":"top","to":["6","12"]},
                {"kind":"line","id":"left","to":["6","0"]}
            ]]},
            {"id":"ring","op":"revolve","profile":"section","axis":{
                "point":strings(collar_point(moved, ["0", "0", "0"])),
                "direction":strings(collar_vector(moved, ["0", "1", "0"]))},"angle":"360"},
            {"id":"collar","op":"chamfer","body":"ring","distance":"1","edges":[
                {"between":[side("right"), side("bottom")]},
                {"between":[side("right"), side("top")]}
            ]},
            {"id":"drilled","op":"hole","target":"collar","at":[{"id":"radial-0",
                "point":strings(collar_point(moved, point)),
                "direction":strings(collar_vector(moved, direction))}],
                "diameter":"2","depth":"through","tip":"flat"}
        ],"result":"drilled"})
        .to_string(),
    )
    .unwrap()
}

#[test]
fn a_plain_hole_strictly_inside_a_bevel_ring_builds() {
    // The hole's radial extent is [a - 1, a + 1] and each bevel's is [11, 12].
    for moved in [false, true] {
        for a in ["9", "9999/1000"] {
            let recipe = collar(moved, [a, "12", "0"], ["0", "-1", "0"]);
            let (built, work) =
                with_operation_budget(20_000_000, || evaluate(&recipe, &mut Decider::default()));
            let mut built = built.unwrap_or_else(|e| panic!("{a}, moved {moved}: {e}: {work:?}"));
            assert!(!work.exhausted);
            assert!(built.model.check(built.result).is_empty());
            assert!(built.model.body(built.result).marks.is_empty());

            let field = Plan::compile(&recipe).unwrap();
            let mut workspace = Workspace::default();
            for (point, material) in [
                ([a, "6", "0"], false),
                ([a, "6", "3/2"], true),
                (["0", "6", "0"], false),
                (["47/4", "47/4", "0"], false),
                (["11", "47/4", "0"], true),
            ] {
                let [x, y, z] = collar_point(moved, point);
                let point = RPoint3::new(x, y, z);
                let (solid, located) = if material {
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
                    Some(&solid),
                    "{a}, moved {moved}, {point:?}"
                );
                assert_eq!(
                    field.proven(&point, &mut workspace).location,
                    located,
                    "{a}, moved {moved}, {point:?}"
                );
            }

            let origin = |value: Value| serde_json::from_value::<Origin>(value).unwrap();
            let wall =
                origin(json!({"step":"drilled","role":{"hole":{"at":"radial-0","part":"wall"}}}));
            let faces = built.model.faces_of_body(built.result);
            let named = |origin: &Origin| {
                faces
                    .iter()
                    .copied()
                    .filter(|f| built.model.face_origins(*f).contains(origin))
                    .collect::<Vec<_>>()
            };
            let walls = named(&wall);
            assert_eq!(walls.len(), 1);
            assert!(matches!(
                built
                    .model
                    .geom()
                    .surface(built.model.face(walls[0]).surface),
                Surface::Cylinder(_)
            ));
            let origins: Vec<_> = faces
                .iter()
                .flat_map(|f| built.model.face_origins(*f))
                .collect();
            assert!(
                origins
                    .iter()
                    .all(|o| !matches!(o.role, Role::Hole { .. }) || **o == wall),
                "{origins:?}"
            );
            let mut bevels: Vec<_> = origins
                .iter()
                .filter(|o| matches!(o.role, Role::Bevel { .. }))
                .collect();
            bevels.sort();
            bevels.dedup();
            assert_eq!(bevels.len(), 2, "{origins:?}");
            for segment in ["right", "left"] {
                let side = origin(json!({"step":"ring","role":{"side":{"segment":segment}}}));
                assert!(!named(&side).is_empty(), "{segment}: {origins:?}");
            }

            let volume = catallax_measure::volume(
                &built.model,
                built.result,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap()
            .value;
            let expected = Ball::from_rational(&q("3782/3"), 128) * Ball::pi(128);
            assert!(volume.lo <= expected.hi() && expected.lo() <= volume.hi);
            assert!(volume.width() <= q("1/200"));
        }
    }
}

#[test]
fn a_hole_meeting_a_bevel_keeps_its_typed_refusal() {
    let mut rows = Vec::new();
    for moved in [false, true] {
        // Exact contact at (11, 12, 0), then a genuine crossing curve.
        rows.push(collar(moved, ["10", "12", "0"], ["0", "-1", "0"]));
        rows.push(collar(moved, ["21/2", "12", "0"], ["0", "-1", "0"]));
    }
    // A skew axis crosses the retained upper cone, so no parallel tube applies.
    rows.push(collar(false, ["21/2", "12", "0"], ["1", "-1", "1"]));
    for recipe in rows {
        let (built, work) =
            with_operation_budget(20_000_000, || evaluate(&recipe, &mut Decider::default()));
        assert!(!work.exhausted, "{work:?}");
        match built {
            Err(catallax_features::BuildError::Geometry { step, refusal }) => {
                assert_eq!(step.as_str(), "drilled");
                assert_eq!(refusal.code, "intersect.traced-pair-required");
            }
            other => panic!("{:?}", other.map(|_| "built")),
        }
    }
}

#[test]
fn a_boss_builds_between_the_apps_existing_countersinks() {
    let recipe = Recipe::parse(
        &json!({
            "format":"catallax-recipe/3", "name":"app countersinks then boss",
            "steps":[
                {"id":"plate","op":"block","size":["60","40","5"],
                 "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
                {"id":"drilled","op":"hole","target":"plate","at":[
                    {"id":"left","point":["15","20","5"],"direction":["0","0","-1"]},
                    {"id":"right","point":["45","20","5"],"direction":["0","0","-1"]}
                 ],"diameter":"5","depth":"through","tip":"flat",
                 "head":{"countersink":{"diameter":"9","angle":"90"}}},
                {"id":"boss","op":"cylinder","radius":"8","height":"10",
                 "pose":{"rotation":["1","0","0","0"],"translation":["30","20","5"]}},
                {"id":"body","op":"boolean","kind":"union","target":"drilled","tool":"boss"}
            ], "result":"body"
        })
        .to_string(),
    )
    .unwrap();
    let (built, work) = catallax_numbers::with_operation_budget(100_000_000, || {
        evaluate(&recipe, &mut Decider::default())
    });
    let mut built = built.unwrap_or_else(|e| panic!("{e}: {work:?}"));
    assert!(!work.exhausted);
    assert!(built.model.check(built.result).is_empty());
    assert!(built.model.body(built.result).marks.is_empty());
    let field = Plan::compile(&recipe).unwrap();
    let mut workspace = Workspace::default();
    for (point, material) in [
        (["15", "20", "2"], false),
        (["45", "20", "2"], false),
        (["30", "20", "10"], true),
        (["39", "20", "10"], false),
        (["20", "20", "2"], true),
        (["19", "20", "19/4"], false),
    ] {
        let [x, y, z] = point.map(q);
        let p = RPoint3::new(x, y, z);
        assert_eq!(
            classify_point(&mut built.model, built.result, &p, &mut Decider::default()).proven(),
            Some(&if material {
                Containment::Inside
            } else {
                Containment::Outside
            })
        );
        assert_eq!(
            field.proven(&p, &mut workspace).location,
            if material {
                Location::Inside
            } else {
                Location::Outside
            }
        );
    }
    let volume = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let expected = catallax_numbers::Ball::from_rational(&q("12000"), 128)
        + catallax_numbers::Ball::from_rational(&q("3313/6"), 128)
            * catallax_numbers::Ball::pi(128);
    assert!(volume.value.lo <= expected.lo() && expected.hi() <= volume.value.hi);
}
