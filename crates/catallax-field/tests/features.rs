//! Definition-level signs, trims, origin sets and reproducibility for profile sweeps.
use catallax_field::{Direction, FaceState, Location, Plan, Query, Value, Workspace};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Origin, Recipe, Role, StepId};
use catallax_shapes::RPoint3;
use serde_json::{Value as Json, json};
use std::collections::BTreeSet;

fn pose() -> Json {
    json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]})
}
fn rectangle() -> Json {
    json!([[{"kind":"line","id":"bottom","to":["4","0"]},{"kind":"line","id":"right","to":["4","3"]},{"kind":"line","id":"top","to":["0","3"]},{"kind":"line","id":"left","to":["0","0"]}]])
}
fn circles(radii: &[&str], centre: [&str; 2]) -> Json {
    Json::Array(radii.iter().enumerate().map(|(i,r)|json!([{"kind":"circle","id":format!("circle{i}"),"centre":centre,"radius":r}])).collect())
}
fn slot() -> Json {
    json!([[{"kind":"line","id":"bottom","to":["2","-1"]},{"kind":"arc","id":"right","centre":["2","0"],"to":["2","1"],"turn":"ccw"},{"kind":"line","id":"top","to":["-2","1"]},{"kind":"arc","id":"left","centre":["-2","0"],"to":["-2","-1"],"turn":"ccw"}]])
}
fn recipe(loops: Json, operation: Json) -> Recipe {
    let value = json!({"format":"catallax-recipe/3","name":"feature","steps":[{"id":"profile","op":"sketch","plane":pose(),"loops":loops},operation],"result":"feature"});
    Recipe::parse(&value.to_string()).unwrap()
}
fn extrude(loops: Json, distance: &str) -> Recipe {
    recipe(
        loops,
        json!({"id":"feature","op":"extrude","profile":"profile","distance":distance}),
    )
}
fn revolve(loops: Json, angle: &str) -> Recipe {
    recipe(
        loops,
        json!({"id":"feature","op":"revolve","profile":"profile","axis":{"point":["0","0","0"],"direction":["0","1","0"]},"angle":angle}),
    )
}
fn point(p: [&str; 3]) -> RPoint3 {
    let [x, y, z] = p.map(|s| Rational::parse(s).unwrap());
    RPoint3::new(x, y, z)
}
fn query(plan: &Plan, p: [&str; 3]) -> Query {
    plan.proven(&point(p), &mut Workspace::default())
}
fn check(plan: &Plan, cases: &[([&str; 3], Location)]) {
    for &(p, location) in cases {
        let q = query(plan, p);
        assert_eq!(q.location, location, "at {p:?}: {q:?}");
        assert!(
            matches!(q.value, Value::Residual(_)),
            "sweep residual must not claim true distance"
        );
        assert!(q.distance.millimetres.lo() >= 0.0);
        if location == Location::Boundary {
            assert_eq!(q.distance.millimetres.hi().to_bits(), 0.0_f64.to_bits());
        }
        let floats = p.map(|v| Rational::parse(v).unwrap().to_f64_nearest());
        let fast = plan.fast(floats, &mut Workspace::default()).unwrap();
        assert!(!fast.true_distance_formula);
        if q.value.enclosure().hi() < 0.0 {
            assert!(fast.value <= 0.0, "fast sign {p:?}: {fast:?}");
        }
        if q.value.enclosure().lo() > 0.0 {
            assert!(fast.value >= 0.0, "fast sign {p:?}: {fast:?}");
        }
    }
}
fn normal(direction: Direction, want: [i64; 3]) {
    let Direction::Unique(n) = direction else {
        panic!("expected unique normal, got {direction:?}");
    };
    for (n, w) in n.into_iter().zip(want) {
        assert!(n.contains(w as f64), "{n:?} excludes {w}");
    }
}
fn origins(plan: &Plan, p: [&str; 3]) -> BTreeSet<Origin> {
    plan.faces()
        .iter()
        .filter_map(|f| {
            let q = plan
                .face_query(
                    plan.select(&f.origin).unwrap(),
                    &point(p),
                    &mut Workspace::default(),
                )
                .unwrap();
            (q.state == FaceState::OnFace).then(|| f.origin.clone())
        })
        .collect()
}
fn role_names(set: &BTreeSet<Origin>) -> BTreeSet<String> {
    set.iter()
        .map(|o| match &o.role {
            Role::Start => "start".into(),
            Role::End => "end".into(),
            Role::Side { segment, .. } => segment.to_string(),
            _ => panic!("unexpected role"),
        })
        .collect()
}

#[test]
fn narrow_planar_chambers_are_not_missed_by_fixed_coordinate_rays() {
    let recipe = Recipe::parse(&json!({
        "format":"catallax-recipe/3", "name":"narrow-chamber",
        "steps":[
            {"id":"profile","op":"sketch","plane":pose(),"loops":[[
                {"kind":"line","id":"lower","to":["1","2"]},
                {"kind":"line","id":"end","to":["1","2000000000001/1000000000000"]},
                {"kind":"line","id":"upper","to":["0","0"]}
            ]]},
            {"id":"wedge","op":"extrude","profile":"profile","distance":"2"},
            {"id":"remote","op":"block","size":["1","1","1"],"pose":{"rotation":["1","0","0","0"],"translation":["10","0","0"]}},
            {"id":"result","op":"boolean","kind":"union","target":"wedge","tool":"remote"}
        ],
        "result":"result"
    }).to_string()).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    assert_eq!(query(&plan, ["0", "0", "1"]).location, Location::Boundary);
    assert_eq!(
        query(&plan, ["1/2", "4000000000001/4000000000000", "1"]).location,
        Location::Inside
    );
    assert_eq!(query(&plan, ["1/2", "1", "1"]).location, Location::Boundary);
    assert_eq!(
        query(&plan, ["1/2", "999/1000", "1"]).location,
        Location::Outside
    );
}

#[test]
fn flush_extrusion_caps_regularise_at_face_interiors_edges_and_corners() {
    for kind in ["union", "subtract", "intersect"] {
        let mut upper_pose = pose();
        upper_pose["translation"] = json!(["0", "0", "2"]);
        let recipe = Recipe::parse(
            &json!({
                "format":"catallax-recipe/3", "name":"flush-extrusions",
                "steps":[
                    {"id":"lower_profile","op":"sketch","plane":pose(),"loops":rectangle()},
                    {"id":"upper_profile","op":"sketch","plane":upper_pose,"loops":rectangle()},
                    {"id":"lower","op":"extrude","profile":"lower_profile","distance":"2"},
                    {"id":"upper","op":"extrude","profile":"upper_profile","distance":"2"},
                    {"id":"result","op":"boolean","kind":kind,"target":"lower","tool":"upper"}
                ],
                "result":"result"
            })
            .to_string(),
        )
        .unwrap();
        let plan = Plan::compile(&recipe).unwrap();
        for p in [["2", "1", "2"], ["0", "1", "2"], ["0", "0", "2"]] {
            let expected = match kind {
                "intersect" => Location::Outside,
                "union" if p[0] == "2" => Location::Inside,
                _ => Location::Boundary,
            };
            assert_eq!(query(&plan, p).location, expected, "{kind} {p:?}");
        }
        if kind != "subtract" {
            assert!(origins(&plan, ["2", "1", "2"]).is_empty(), "{kind}");
        }
    }
}

#[test]
fn signed_rectangles_oblique_lines_and_cap_normals() {
    for (distance, z, inside_normal, end_normal) in [
        ("2", "1", [0, 0, -1], [0, 0, 1]),
        ("-2", "-1", [0, 0, 1], [0, 0, -1]),
    ] {
        let plan = Plan::compile(&extrude(rectangle(), distance)).unwrap();
        check(
            &plan,
            &[
                (["2", "1", z], Location::Inside),
                (["5", "1", z], Location::Outside),
                (["4", "1", z], Location::Boundary),
                (["2", "1", "0"], Location::Boundary),
                (["2", "1", distance], Location::Boundary),
            ],
        );
        normal(query(&plan, ["2", "1", "0"]).direction, inside_normal);
        normal(query(&plan, ["2", "1", distance]).direction, end_normal);
        normal(query(&plan, ["4", "1", z]).direction, [1, 0, 0]);
        assert_eq!(query(&plan, ["4", "0", z]).direction, Direction::Ambiguous);
        let q = query(&plan, ["2", "1", z]);
        assert!(q.distance.millimetres.lo() <= 1.0 && q.distance.millimetres.hi() >= 1.0);
    }
    let mut r = extrude(rectangle(), "5");
    if let catallax_recipe::Operation::Extrude { direction, .. } = &mut r.steps[1].op {
        *direction = Some(catallax_recipe::Vector3([3.into(), 0.into(), 4.into()]));
    }
    let plan = Plan::compile(&r).unwrap();
    check(
        &plan,
        &[
            (["7/2", "1", "2"], Location::Inside),
            (["3/2", "1", "2"], Location::Boundary),
            (["0", "1", "2"], Location::Outside),
            (["5", "1", "4"], Location::Boundary),
        ],
    );
}

#[test]
fn circles_slots_holes_and_islands_are_exact_not_polygonal() {
    let circle = Plan::compile(&extrude(circles(&["2"], ["0", "0"]), "4")).unwrap();
    check(
        &circle,
        &[
            (["0", "0", "2"], Location::Inside),
            (["2", "0", "2"], Location::Boundary),
            (["6/5", "8/5", "2"], Location::Boundary),
            (["3", "0", "2"], Location::Outside),
            (["0", "0", "5"], Location::Outside),
        ],
    );
    normal(query(&circle, ["2", "0", "2"]).direction, [1, 0, 0]);
    let slot = Plan::compile(&extrude(slot(), "2")).unwrap();
    check(
        &slot,
        &[
            (["0", "0", "1"], Location::Inside),
            (["5/2", "0", "1"], Location::Inside),
            (["3", "0", "1"], Location::Boundary),
            (["2", "1", "1"], Location::Boundary),
            (["0", "2", "1"], Location::Outside),
            (["13/5", "4/5", "1"], Location::Boundary),
        ],
    );
    normal(query(&slot, ["2", "1", "1"]).direction, [0, 1, 0]);
    let ring = Plan::compile(&extrude(circles(&["3", "1", "1/2"], ["0", "0"]), "2")).unwrap();
    check(
        &ring,
        &[
            (["2", "0", "1"], Location::Inside),
            (["3/4", "0", "1"], Location::Outside),
            (["0", "0", "1"], Location::Inside),
            (["1", "0", "1"], Location::Boundary),
            (["1/2", "0", "1"], Location::Boundary),
        ],
    );
    normal(query(&ring, ["1", "0", "1"]).direction, [-1, 0, 0]);
    normal(query(&ring, ["1/2", "0", "1"]).direction, [1, 0, 0]);
}

#[test]
fn revolved_circles_hollow_rings_and_cones() {
    let ring = Plan::compile(&revolve(circles(&["1", "1/2"], ["4", "0"]), "360")).unwrap();
    check(
        &ring,
        &[
            (["4", "0", "0"], Location::Outside),
            (["19/4", "0", "0"], Location::Inside),
            (["5", "0", "0"], Location::Boundary),
            (["9/2", "0", "0"], Location::Boundary),
            (["0", "0", "4"], Location::Outside),
            (["0", "0", "5"], Location::Boundary),
            (["0", "2", "4"], Location::Outside),
        ],
    );
    normal(query(&ring, ["5", "0", "0"]).direction, [1, 0, 0]);
    normal(query(&ring, ["9/2", "0", "0"]).direction, [-1, 0, 0]);
    let triangle = json!([[{"kind":"line","id":"base","to":["3","0"]},{"kind":"line","id":"slope","to":["0","4"]},{"kind":"line","id":"axis","to":["0","0"]}]]);
    let cone = Plan::compile(&revolve(triangle, "360")).unwrap();
    check(
        &cone,
        &[
            (["0", "1", "0"], Location::Inside),
            (["3/2", "2", "0"], Location::Boundary),
            (["0", "2", "3/2"], Location::Boundary),
            (["2", "2", "0"], Location::Outside),
            (["0", "4", "0"], Location::Boundary),
            (["0", "0", "0"], Location::Boundary),
        ],
    );
    assert_eq!(
        query(&cone, ["0", "4", "0"]).direction,
        Direction::Ambiguous
    );
    assert!(
        cone.faces()
            .iter()
            .all(|f| !matches!(&f.origin.role, Role::Start | Role::End)
                && !matches!(&f.origin.role,Role::Side{segment,..}if segment.to_string()=="axis"))
    );
}

#[test]
fn partial_turns_caps_axis_merges_and_angular_complements() {
    for angle in ["90", "180", "270"] {
        let plan = Plan::compile(&revolve(rectangle(), angle)).unwrap();
        check(
            &plan,
            &[
                (["1", "1", "-1"], Location::Inside),
                (["1", "1", "0"], Location::Boundary),
                (["5", "1", "0"], Location::Outside),
            ],
        );
        let expected = if angle == "90" {
            Location::Outside
        } else {
            Location::Inside
        };
        assert_eq!(query(&plan, ["-1", "1", "-1"]).location, expected);
        assert_eq!(query(&plan, ["1", "1", "1"]).location, Location::Outside);
        let end = if angle == "90" {
            ["0", "1", "-1"]
        } else if angle == "180" {
            ["-1", "1", "0"]
        } else {
            ["0", "1", "1"]
        };
        assert_eq!(query(&plan, end).location, Location::Boundary);
        let names = role_names(&origins(&plan, end));
        if angle == "180" {
            assert_eq!(names, BTreeSet::from(["start".into(), "end".into()]));
        } else {
            assert_eq!(names, BTreeSet::from(["end".into()]));
        }
        if angle == "180" {
            normal(query(&plan, ["0", "1", "0"]).direction, [0, 0, 1]);
        } else {
            assert_eq!(
                query(&plan, ["0", "1", "0"]).direction,
                Direction::Ambiguous
            );
        }
    }
    let plan = Plan::compile(&revolve(circles(&["1"], ["4", "0"]), "180")).unwrap();
    assert_eq!(
        role_names(&origins(&plan, ["4", "0", "0"])),
        BTreeSet::from(["start".into()])
    );
    assert_eq!(
        role_names(&origins(&plan, ["-4", "0", "0"])),
        BTreeSet::from(["end".into()])
    );
}

fn through_recipe(direction: Option<[&str; 3]>) -> Recipe {
    let mut through =
        json!({"id":"cut","op":"extrude_through","target":"stock","profile":"profile"});
    if let Some(d) = direction {
        through["direction"] = json!(d);
    }
    Recipe::parse(&json!({"format":"catallax-recipe/3","name":"cut","steps":[{"id":"stock","op":"block","size":["10","10","10"],"pose":pose()},{"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["5","5","1000000"]},"loops":circles(&["1"],["0","0"])},through],"result":"cut"}).to_string()).unwrap()
}
#[test]
fn through_all_is_infinite_and_has_no_construction_caps() {
    for direction in [None, Some(["0", "0", "-7"])] {
        let r = through_recipe(direction);
        let plan = Plan::compile(&r).unwrap();
        check(
            &plan,
            &[
                (["5", "5", "1"], Location::Outside),
                (["5", "5", "9"], Location::Outside),
                (["6", "5", "5"], Location::Boundary),
                (["7", "5", "5"], Location::Inside),
                (["7", "5", "0"], Location::Boundary),
            ],
        );
        normal(query(&plan, ["6", "5", "5"]).direction, [-1, 0, 0]);
        assert!(
            plan.faces()
                .iter()
                .filter(|f| f.origin.step.to_string() == "cut")
                .all(|f| matches!(f.origin.role, Role::Side { .. }))
        );
    }
}

#[test]
fn source_names_match_brep_aliases_and_survive_splitting() {
    let loops = json!([[{"kind":"line","id":"bottom_a","to":["2","0"]},{"kind":"line","id":"bottom_b","to":["4","0"]},{"kind":"line","id":"right","to":["4","3"]},{"kind":"line","id":"top","to":["0","3"]},{"kind":"line","id":"left","to":["0","0"]}]]);
    for r in [
        extrude(loops, "2"),
        extrude(slot(), "2"),
        revolve(rectangle(), "180"),
        revolve(circles(&["1"], ["4", "0"]), "360"),
        through_recipe(None),
    ] {
        let built = catallax_features::evaluate(&r, &mut Decider::default()).unwrap();
        let plan = Plan::compile(&r).unwrap();
        let brep: BTreeSet<_> = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .flat_map(|f| built.model.face_origins(f).to_vec())
            .collect();
        for origin in brep {
            assert!(plan.select(&origin).is_some(), "missing {origin:?}");
        }
    }
    let plan=Plan::compile(&extrude(json!([[{"kind":"line","id":"bottom_a","to":["2","0"]},{"kind":"line","id":"bottom_b","to":["4","0"]},{"kind":"line","id":"right","to":["4","3"]},{"kind":"line","id":"top","to":["0","3"]},{"kind":"line","id":"left","to":["0","0"]}]]),"2")).unwrap();
    assert_eq!(
        role_names(&origins(&plan, ["3", "0", "1"])),
        BTreeSet::from(["bottom_a".into(), "bottom_b".into()])
    );
    let mut r = extrude(rectangle(), "2");
    r.steps.push(catallax_recipe::Step {
        id: StepId::new("tool").unwrap(),
        op: catallax_recipe::Operation::Block {
            size: [1.into(), 3.into(), 2.into()],
            pose: catallax_recipe::Pose {
                rotation: [1.into(), 0.into(), 0.into(), 0.into()],
                translation: catallax_recipe::Vector3([1.into(), 0.into(), 0.into()]),
            },
        },
    });
    r.steps.push(catallax_recipe::Step {
        id: StepId::new("split").unwrap(),
        op: catallax_recipe::Operation::Boolean {
            kind: catallax_recipe::BooleanKind::Subtract,
            target: StepId::new("feature").unwrap(),
            tool: StepId::new("tool").unwrap(),
        },
    });
    r.result = StepId::new("split").unwrap();
    let plan = Plan::compile(&r).unwrap();
    for p in [["1/2", "0", "1"], ["3", "0", "1"]] {
        assert!(role_names(&origins(&plan, p)).contains("bottom"));
    }
    assert!(!role_names(&origins(&plan, ["3/2", "0", "1"])).contains("bottom"));
}

#[test]
fn invalid_profiles_and_unfinished_algebraic_cases_share_specific_codes() {
    let mut irrational = extrude(rectangle(), "2");
    if let catallax_recipe::Operation::Extrude { direction, .. } = &mut irrational.steps[1].op {
        *direction = Some(catallax_recipe::Vector3([1.into(), 0.into(), 1.into()]));
    }
    for r in [
        irrational,
        revolve(rectangle(), "45"),
        extrude(circles(&["0"], ["0", "0"]), "2"),
        extrude(rectangle(), "0"),
        revolve(rectangle(), "0"),
        revolve(circles(&["2"], ["0", "0"]), "360"),
    ] {
        let field = Plan::compile(&r).unwrap_err();
        let catallax_features::BuildError::Geometry { refusal, .. } =
            catallax_features::evaluate(&r, &mut Decider::default()).unwrap_err()
        else {
            panic!("expected geometric refusal")
        };
        assert_eq!(field.code, refusal.code);
        assert_ne!(field.code, "field.impossible-primitive");
    }
}

fn bits(plan: &Plan, p: &RPoint3, workspace: &mut Workspace) -> Vec<u64> {
    let q = plan.proven(p, workspace);
    let mut result = vec![match q.location {
        Location::Inside => 0,
        Location::Outside => 1,
        Location::Boundary => 2,
        Location::Indeterminate => 3,
        _ => unreachable!(),
    }];
    for i in [q.value.enclosure(), q.distance.millimetres] {
        result.extend([i.lo().to_bits(), i.hi().to_bits()]);
    }
    match q.direction {
        Direction::Unique(n) => {
            result.push(0);
            for i in n {
                result.extend([i.lo().to_bits(), i.hi().to_bits()]);
            }
        }
        Direction::Ambiguous => result.push(1),
        Direction::Undefined => result.push(2),
        _ => unreachable!(),
    }
    let fast = plan
        .fast(
            [
                p.x.to_f64_nearest(),
                p.y.to_f64_nearest(),
                p.z.to_f64_nearest(),
            ],
            workspace,
        )
        .unwrap();
    result.push(fast.value.to_bits());
    result.push(u64::from(fast.true_distance_formula));
    result
}
#[test]
fn sweeps_have_identical_complete_answers_at_one_two_and_four_threads() {
    let points: Vec<_> = (-2..=2)
        .flat_map(|x| {
            (-1..=1).flat_map(move |y| {
                (-1..=1).map(move |z| [x.to_string(), y.to_string(), z.to_string()])
            })
        })
        .map(|p| point(p.each_ref().map(String::as_str)))
        .collect();
    for r in [
        extrude(rectangle(), "-2"),
        extrude(slot(), "2"),
        extrude(circles(&["3", "1"], ["0", "0"]), "2"),
        revolve(rectangle(), "270"),
        revolve(circles(&["1"], ["4", "0"]), "360"),
        through_recipe(None),
    ] {
        let plan = Plan::compile(&r).unwrap();
        let expected: Vec<_> = points
            .iter()
            .map(|p| bits(&plan, p, &mut Workspace::default()))
            .collect();
        for count in [1, 2, 4] {
            let got = std::thread::scope(|scope| {
                let handles: Vec<_> = points
                    .chunks(points.len().div_ceil(count))
                    .map(|chunk| {
                        let plan = &plan;
                        scope.spawn(move || {
                            let mut workspace = Workspace::default();
                            chunk
                                .iter()
                                .map(|p| bits(plan, p, &mut workspace))
                                .collect::<Vec<_>>()
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    .flat_map(|h| h.join().unwrap())
                    .collect::<Vec<_>>()
            });
            assert_eq!(got, expected);
        }
    }
}

#[test]
fn all_boolean_kinds_combine_profile_features() {
    use catallax_recipe::{BooleanKind, Operation, Step};
    for kind in [
        BooleanKind::Union,
        BooleanKind::Subtract,
        BooleanKind::Intersect,
    ] {
        let mut r = extrude(rectangle(), "2");
        r.steps.push(Step {
            id: StepId::new("other_profile").unwrap(),
            op: Operation::Sketch {
                plane: catallax_recipe::Pose {
                    rotation: [1.into(), 0.into(), 0.into(), 0.into()],
                    translation: catallax_recipe::Vector3([2.into(), 0.into(), 0.into()]),
                },
                loops: match &r.steps[0].op {
                    Operation::Sketch { loops, .. } => loops.clone(),
                    _ => unreachable!(),
                },
            },
        });
        r.steps.push(Step {
            id: StepId::new("other").unwrap(),
            op: Operation::Extrude {
                profile: StepId::new("other_profile").unwrap(),
                distance: 2.into(),
                direction: None,
            },
        });
        r.steps.push(Step {
            id: StepId::new("result").unwrap(),
            op: Operation::Boolean {
                kind,
                target: StepId::new("feature").unwrap(),
                tool: StepId::new("other").unwrap(),
            },
        });
        r.result = StepId::new("result").unwrap();
        let plan = Plan::compile(&r).unwrap();
        let expected = match kind {
            BooleanKind::Union => [Location::Inside, Location::Inside, Location::Inside],
            BooleanKind::Subtract => [Location::Inside, Location::Outside, Location::Outside],
            BooleanKind::Intersect => [Location::Outside, Location::Inside, Location::Outside],
        };
        for (x, location) in ["1", "3", "5"].into_iter().zip(expected) {
            assert_eq!(query(&plan, [x, "1", "1"]).location, location);
        }
        let seam = query(&plan, ["2", "1", "1"]);
        assert_eq!(
            seam.location,
            if kind == BooleanKind::Union {
                Location::Inside
            } else {
                Location::Boundary
            }
        );
        if kind == BooleanKind::Subtract {
            normal(seam.direction, [1, 0, 0]);
        }
    }
}

#[test]
fn placements_empty_material_and_nonunit_axis_definitions() {
    let mut r = extrude(rectangle(), "2");
    if let catallax_recipe::Operation::Sketch { plane, .. } = &mut r.steps[0].op {
        plane.rotation = [1.into(), 1.into(), 1.into(), 1.into()];
        plane.translation.0 = [10.into(), 20.into(), 30.into()];
    }
    let plan = Plan::compile(&r).unwrap();
    check(
        &plan,
        &[
            (["11", "22", "31"], Location::Inside),
            (["11", "24", "31"], Location::Boundary),
            (["11", "25", "31"], Location::Outside),
        ],
    );
    normal(query(&plan, ["11", "24", "31"]).direction, [0, 1, 0]);
    let empty = Plan::compile(&extrude(json!([]), "2")).unwrap();
    assert_eq!(query(&empty, ["0", "0", "0"]).location, Location::Outside);
    assert_eq!(
        query(&empty, ["0", "0", "0"])
            .distance
            .millimetres
            .lo()
            .to_bits(),
        f64::INFINITY.to_bits()
    );
    let shifted = json!([[{"kind":"line","id":"a","to":["3","0"]},{"kind":"line","id":"b","to":["3","2"]},{"kind":"line","id":"c","to":["1","2"]},{"kind":"line","id":"d","to":["1","0"]}]]);
    let mut r = revolve(shifted, "360");
    if let catallax_recipe::Operation::Revolve { axis, .. } = &mut r.steps[1].op {
        axis.direction.0 = [0.into(), 7.into(), 0.into()];
    }
    let plan = Plan::compile(&r).unwrap();
    check(
        &plan,
        &[
            (["2", "1", "0"], Location::Inside),
            (["0", "1", "0"], Location::Outside),
            (["1", "1", "0"], Location::Boundary),
            (["3", "1", "0"], Location::Boundary),
            (["0", "1", "2"], Location::Inside),
        ],
    );
    normal(query(&plan, ["1", "1", "0"]).direction, [-1, 0, 0]);
}

#[test]
fn rational_bisector_arc_keeps_its_irrational_radius_and_trims() {
    let part = extrude(
        json!([[
            {"kind":"arc","id":"round","centre":["1","1/2"],"to":["0","0"],"turn":"ccw"},
            {"kind":"line","id":"chord","to":["2","0"]}
        ]]),
        "3",
    );
    let plan = Plan::compile(&part).unwrap();
    check(
        &plan,
        &[
            (["1", "1", "1"], Location::Inside),
            (["1", "-1/4", "1"], Location::Outside),
            (["0", "1", "1"], Location::Boundary),
            (["2", "1", "1"], Location::Boundary),
            (["1", "1", "0"], Location::Boundary),
            (["1", "1", "3"], Location::Boundary),
            (["1", "1", "4"], Location::Outside),
        ],
    );
    let built = catallax_features::evaluate(&part, &mut Decider::default()).unwrap();
    for face in built.model.faces_of_body(built.result) {
        if let catallax_shapes::Surface::Cylinder(c) =
            built.model.geom().surface(built.model.face(face).surface)
        {
            assert_eq!(
                c.radius_sq.as_rational(),
                Some(Rational::parse("5/4").unwrap())
            );
        }
    }
}

#[test]
fn generated_algebraic_radii_and_spindle_sheets_keep_the_shared_refusals() {
    let spherical = json!([[{"kind":"arc","id":"round","centre":["0","0"],"to":["1","1"],"turn":"ccw"},{"kind":"line","id":"chord","to":["1","-1"]}]]);
    let toroidal = json!([[{"kind":"arc","id":"round","centre":["3","0"],"to":["4","1"],"turn":"ccw"},{"kind":"line","id":"chord","to":["4","-1"]}]]);
    let spindle = json!([[{"kind":"arc","id":"round","centre":["1","0"],"to":["1","2"],"turn":"ccw"},{"kind":"line","id":"chord","to":["1","-2"]}]]);
    for (r, code) in [
        (
            revolve(spherical, "360"),
            "features.unsupported-algebraic-sphere-radius",
        ),
        (
            revolve(toroidal, "360"),
            "features.unsupported-algebraic-torus-radius",
        ),
        (
            revolve(spindle, "360"),
            "features.unsupported-spindle-torus",
        ),
    ] {
        assert_eq!(Plan::compile(&r).unwrap_err().code, code);
        let catallax_features::BuildError::Geometry { refusal, .. } =
            catallax_features::evaluate(&r, &mut Decider::default()).unwrap_err()
        else {
            panic!("expected geometry refusal")
        };
        assert_eq!(refusal.code, code);
    }
}

#[test]
fn algebraic_rigid_frames_preserve_sweep_trims_germs_and_complete_answers() {
    use catallax_numbers::{Real, sin_cos_degrees};
    use catallax_shapes::{Pose, Rotation, Vector3};
    let mut dec = Decider::default();
    let [sine, cosine] = sin_cos_degrees(&Rational::from_ratio(45, 2).unwrap(), &mut dec).unwrap();
    let y = Pose::new(
        Rotation::from_real_quaternion(
            cosine.clone(),
            Real::zero(),
            sine.clone(),
            Real::zero(),
            &mut dec,
        )
        .unwrap(),
        Vector3::zero(),
        &mut dec,
    )
    .unwrap();
    let z = Pose::new(
        Rotation::from_real_quaternion(cosine, Real::zero(), Real::zero(), sine, &mut dec).unwrap(),
        Vector3::zero(),
        &mut dec,
    )
    .unwrap();
    let fixtures = [
        (
            extrude(rectangle(), "2"),
            y,
            vec![
                (["2", "1", "0"], Location::Inside),
                (["1", "1", "1"], Location::Boundary),
                (["0", "4", "0"], Location::Outside),
            ],
        ),
        (
            revolve(rectangle(), "360"),
            z,
            vec![
                (["-1", "1", "0"], Location::Inside),
                (["-1", "1", "4"], Location::Boundary),
                (["-3", "3", "0"], Location::Outside),
                (["0", "0", "1"], Location::Boundary),
            ],
        ),
    ];
    for (index, (recipe, pose, cases)) in fixtures.into_iter().enumerate() {
        let plan = Plan::compile_checked(&recipe, &mut dec)
            .unwrap()
            .transformed(&pose, &mut dec)
            .unwrap();
        check(&plan, &cases);
        let refusal = plan.visit_candidates(2, |_, _| true).unwrap_err();
        assert_eq!(refusal.code, "field.sample-algebraic-placement");
        if index == 0 {
            assert_eq!(
                role_names(&origins(&plan, ["1", "1", "1"])),
                BTreeSet::from(["left".to_string()]),
            );
            let Direction::Unique(n) = query(&plan, ["1", "1", "1"]).direction else {
                panic!("the interior of the rotated extrusion wall has a unique normal");
            };
            let s = std::f64::consts::FRAC_1_SQRT_2;
            for (component, expected) in n.into_iter().zip([-s, 0.0, s]) {
                assert!(
                    component.contains(expected),
                    "{component:?} excludes {expected}"
                );
            }
        } else {
            normal(query(&plan, ["-1", "1", "4"]).direction, [0, 0, 1]);
            assert_eq!(
                role_names(&origins(&plan, ["0", "0", "1"])),
                BTreeSet::from(["bottom".to_string()]),
            );
        }
        let points: Vec<_> = cases.iter().map(|(p, _)| point(*p)).collect();
        let expected: Vec<_> = points
            .iter()
            .map(|p| bits(&plan, p, &mut Workspace::default()))
            .collect();
        for threads in [1, 2, 4] {
            let actual = std::thread::scope(|scope| {
                let handles: Vec<_> = points
                    .chunks(points.len().div_ceil(threads))
                    .map(|chunk| {
                        let plan = &plan;
                        scope.spawn(move || {
                            let mut workspace = Workspace::default();
                            chunk
                                .iter()
                                .map(|p| bits(plan, p, &mut workspace))
                                .collect::<Vec<_>>()
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    .flat_map(|h| h.join().unwrap())
                    .collect::<Vec<_>>()
            });
            assert_eq!(actual, expected);
        }
    }
}

#[test]
fn rational_oblique_circle_has_exact_membership_caps_and_conservative_metric() {
    let mut source = extrude(circles(&["2"], ["0", "0"]), "5");
    if let catallax_recipe::Operation::Extrude { direction, .. } = &mut source.steps[1].op {
        *direction = Some(catallax_recipe::Vector3([3.into(), 0.into(), 4.into()]));
    }
    let plan = Plan::compile(&source).unwrap();
    check(
        &plan,
        &[
            (["3/2", "0", "2"], Location::Inside),
            (["3/2", "2", "2"], Location::Boundary),
            (["3/2", "21/10", "2"], Location::Outside),
            (["0", "0", "0"], Location::Boundary),
            (["3", "0", "4"], Location::Boundary),
            (["3", "0", "5"], Location::Outside),
        ],
    );
    normal(query(&plan, ["3/2", "2", "2"]).direction, [0, 1, 0]);
    normal(query(&plan, ["0", "0", "0"]).direction, [0, 0, -1]);
    normal(query(&plan, ["3", "0", "4"]).direction, [0, 0, 1]);
    let outside = query(&plan, ["3/2", "21/10", "2"]);
    // The 1-Lipschitz elliptic residual is 0.08 mm; true distance here is
    // 0.1 mm. Calling the residual a nearest distance would be incorrect.
    assert!(outside.value.enclosure().contains(0.08));
    assert!(outside.distance.millimetres.lo() > 0.079);
    assert!(outside.distance.millimetres.lo() < 0.081);
    assert!(outside.distance.millimetres.hi().is_infinite());
}
