//! Drill-tool point sets, surviving names and the short-head ledge decision.
use catallax_features::evaluate;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{HolePart, Recipe, Role};
use serde_json::{Value, json};

fn recipe(depth: Value, tip: Value, head: Option<Value>) -> Recipe {
    let mut hole = json!({"id":"hole","op":"hole","target":"stock","at":[{"id":"centre","point":["0","0","0"],"direction":["0","0","1"]}],"diameter":"2","depth":depth,"tip":tip});
    if let Some(head) = head {
        hole["head"] = head;
    }
    Recipe::parse(&json!({"format":"catallax-recipe/3","name":"hole","steps":[
        {"id":"stock","op":"block","size":["10","10","10"],"pose":{"rotation":["1","0","0","0"],"translation":["-5","-5","0"]}},hole],"result":"hole"}).to_string()).unwrap()
}

#[test]
fn blind_flat_and_drill_point_keep_their_distinct_floors() {
    for (tip, expected, absent) in [
        (json!("flat"), HolePart::Floor, HolePart::Tip),
        (json!({"cone":"90"}), HolePart::Tip, HolePart::Floor),
    ] {
        let built = evaluate(
            &recipe(json!({"blind":"4"}), tip, None),
            &mut Decider::default(),
        )
        .unwrap();
        assert!(built.model.check(built.result).is_empty());
        let parts: Vec<_> = built
            .model
            .faces_of_body(built.result)
            .iter()
            .flat_map(|f| built.model.face_origins(*f))
            .filter_map(|o| match o.role {
                Role::Hole { part, .. } => Some(part),
                _ => None,
            })
            .collect();
        assert!(parts.contains(&expected), "{parts:?}");
        assert!(parts.contains(&HolePart::Wall), "{parts:?}");
        assert!(!parts.contains(&absent), "{parts:?}");
        assert!(
            !parts.contains(&HolePart::Entry),
            "a placement on the entry plane must not invent an entry wall: {parts:?}"
        );
    }
}

#[test]
fn short_countersink_with_cone_tip_has_a_named_planar_annular_floor() {
    let recipe = recipe(
        json!({"blind":"1"}),
        json!({"cone":"90"}),
        Some(json!({"countersink":{"diameter":"8","angle":"90"}})),
    );
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let floors: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|f| {
            built.model.face_origins(*f).iter().any(|o| {
                matches!(
                    o.role,
                    Role::Hole {
                        part: HolePart::Floor,
                        ..
                    }
                )
            })
        })
        .collect();
    assert_eq!(floors.len(), 1);
    let surface = built
        .model
        .geom()
        .surface(built.model.face(floors[0]).surface);
    let catallax_shapes::Surface::Plane(plane) = surface else {
        panic!("annular floor must be planar")
    };
    assert_eq!(plane.origin.z.as_rational(), Some(Rational::one()));
    assert_eq!(
        built.model.face(floors[0]).loops.len(),
        2,
        "the centre continues into the point cone"
    );
    let mass = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let pi = catallax_numbers::Ball::pi(128);
    let expected = catallax_numbers::Ball::from_rational(&Rational::from_i64(1000), 128)
        - catallax_numbers::Ball::from_rational(&Rational::from_ratio(38, 3).unwrap(), 128) * pi;
    assert!(mass.value.lo <= expected.hi() && expected.lo() <= mass.value.hi);
}

#[test]
fn a_through_counterbore_keeps_its_ledge_but_has_no_blind_floor_or_tip() {
    let built = evaluate(
        &recipe(
            json!("through"),
            json!({"cone":"118"}),
            Some(json!({"counterbore":{"diameter":"4","depth":"2"}})),
        ),
        &mut Decider::default(),
    )
    .unwrap();
    let parts: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .iter()
        .flat_map(|f| built.model.face_origins(*f))
        .filter_map(|o| match o.role {
            Role::Hole { part, .. } => Some(part),
            _ => None,
        })
        .collect();
    assert!(parts.contains(&HolePart::CounterboreFloor));
    assert!(parts.contains(&HolePart::CounterboreWall));
    assert!(parts.contains(&HolePart::Wall));
    assert!(!parts.contains(&HolePart::Floor));
    assert!(!parts.contains(&HolePart::Tip));
}

#[test]
fn p5_countersunk_hole_keeps_the_revolved_housing_and_both_drill_surfaces() {
    use catallax_numbers::{Ball, with_operation_budget};
    use catallax_shapes::{RVector3, Rotation};

    for (quaternion, translation) in [
        ([1, 0, 0, 0], [0, 0, 0]),
        ([1, 0, 0, 0], [7, 11, -13]),
        ([2, 0, 1, 0], [0, 0, 0]),
        ([1, 2, 3, 4], [7, 11, -13]),
    ] {
        let rotation = Rotation::from_quaternion(
            Rational::from_i64(quaternion[0]),
            Rational::from_i64(quaternion[1]),
            Rational::from_i64(quaternion[2]),
            Rational::from_i64(quaternion[3]),
        )
        .unwrap();
        let transform = |coordinates: &Value, point: bool| {
            let coordinates: [Rational; 3] =
                std::array::from_fn(|i| Rational::parse(coordinates[i].as_str().unwrap()).unwrap());
            let v = rotation.apply(
                &RVector3::new(
                    coordinates[0].clone(),
                    coordinates[1].clone(),
                    coordinates[2].clone(),
                )
                .to_real(),
            );
            let v = [v.x, v.y, v.z];
            json!(std::array::from_fn::<_, 3, _>(|i| {
                (v[i].as_rational().unwrap()
                    + Rational::from_i64(if point { translation[i] } else { 0 }))
                .to_string()
            }))
        };
        let mut value: Value =
            serde_json::from_str(include_str!("../../../parts/holes/P5-countersinks.json"))
                .unwrap();
        value["steps"][0]["plane"] = json!({
            "rotation": quaternion.map(|q| q.to_string()),
            "translation": translation.map(|q| q.to_string())
        });
        for key in ["point", "direction"] {
            value["steps"][1]["axis"][key] =
                transform(&value["steps"][1]["axis"][key], key == "point");
            for at in value["steps"][2]["at"].as_array_mut().unwrap() {
                at[key] = transform(&at[key], key == "point");
            }
        }
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let (built, work) =
            with_operation_budget(100_000_000, || evaluate(&recipe, &mut Decider::default()));
        let built = built.unwrap_or_else(|e| panic!("{quaternion:?}, {translation:?}: {e}"));
        assert!(!work.exhausted, "{work:?}");
        assert!(built.model.check(built.result).is_empty());
        assert!(built.model.body(built.result).marks.is_empty());
        let origins: Vec<_> = built
            .model
            .faces_of_body(built.result)
            .iter()
            .flat_map(|f| built.model.face_origins(*f))
            .collect();
        for k in 0..4 {
            for (part, count) in [
                (HolePart::Countersink, 1),
                (HolePart::Wall, 1),
                (HolePart::Entry, 0),
            ] {
                assert_eq!(
                    origins
                        .iter()
                        .filter(|o| {
                            o.step.as_str() == "drilled"
                                && matches!(&o.role, Role::Hole { at, part: p }
                            if at.as_str() == format!("bolt-{k}") && *p == part)
                        })
                        .count(),
                    count,
                    "bolt-{k}, {part:?}"
                );
            }
        }
        for segment in [
            "flange-bottom",
            "flange-rim",
            "flange-top",
            "body-wall",
            "body-top",
            "inner-wall",
        ] {
            assert!(
                origins.iter().any(|o| o.step.as_str() == "housing"
                    && matches!(&o.role, Role::Side { segment: s, .. } if s.as_str() == segment)),
                "lost housing origin {segment}"
            );
        }
        let volume = catallax_measure::volume(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        let expected =
            Ball::from_rational(&Rational::from_ratio(11932, 3).unwrap(), 128) * Ball::pi(128);
        assert!(volume.value.lo <= expected.lo() && expected.hi() <= volume.value.hi);
    }
}

#[test]
fn entry_and_wall_merge_into_one_physical_cylinder_with_both_origins() {
    let original = recipe(json!("through"), json!("flat"), None);
    let mut value: Value = serde_json::from_str(&original.to_text()).unwrap();
    value["steps"][1]["at"][0]["point"][2] = json!("5");
    let built = evaluate(
        &Recipe::parse(&value.to_string()).unwrap(),
        &mut Decider::default(),
    )
    .unwrap();
    let drilled: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|f| {
            built
                .model
                .face_origins(*f)
                .iter()
                .any(|o| matches!(o.role, Role::Hole { .. }))
        })
        .collect();
    assert_eq!(drilled.len(), 1);
    let parts: Vec<_> = built
        .model
        .face_origins(drilled[0])
        .iter()
        .filter_map(|o| match o.role {
            Role::Hole { part, .. } => Some(part),
            _ => None,
        })
        .collect();
    assert!(
        parts.contains(&HolePart::Entry) && parts.contains(&HolePart::Wall),
        "{parts:?}"
    );
}

#[test]
fn standard_angles_preserve_exact_names_and_independent_point_sets() {
    use catallax_boolean::{Containment, classify_point};
    use catallax_field::{Location, Plan, Workspace};
    use catallax_numbers::with_operation_budget;
    use catallax_shapes::RPoint3;

    for (tip, head, probes) in [
        (
            json!({"cone":"118"}),
            None,
            vec![
                ("1/2", "17/4", false),
                ("3/4", "17/4", true),
                ("0", "19/4", true),
            ],
        ),
        (
            json!("flat"),
            Some(json!({"countersink":{"diameter":"8","angle":"82"}})),
            vec![("3", "1", false), ("7/2", "1", true), ("0", "17/4", true)],
        ),
        (
            json!("flat"),
            Some(json!({"countersink":{"diameter":"8","angle":"100"}})),
            vec![("5/2", "1", false), ("3", "1", true), ("0", "17/4", true)],
        ),
    ] {
        let recipe = recipe(json!({"blind":"4"}), tip, head);
        let (built, work) =
            with_operation_budget(100_000_000, || evaluate(&recipe, &mut Decider::default()));
        let mut built = built.unwrap_or_else(|e| panic!("{e}: {work:?}"));
        assert!(!work.exhausted);
        assert!(built.model.check(built.result).is_empty());
        assert!(built.model.body(built.result).marks.is_empty());
        let plan = Plan::compile(&recipe).unwrap();
        let mut workspace = Workspace::default();
        // Independent meridian inequalities: tan(59°) is between 3/2 and 2,
        // tan(41°) between 4/5 and 1, and tan(50°) between 1 and 5/4.
        for (x, z, material) in probes
            .into_iter()
            .chain([("1/2", "2", false), ("9/2", "2", true)])
        {
            let p = RPoint3::new(
                Rational::parse(x).unwrap(),
                Rational::zero(),
                Rational::parse(z).unwrap(),
            );
            let answer =
                classify_point(&mut built.model, built.result, &p, &mut Decider::default());
            assert_eq!(
                answer.proven(),
                Some(&if material {
                    Containment::Inside
                } else {
                    Containment::Outside
                }),
                "{x}, {z}: {answer:?}"
            );
            assert_eq!(
                plan.proven(&p, &mut workspace).location,
                if material {
                    Location::Inside
                } else {
                    Location::Outside
                }
            );
        }
        assert!(
            built
                .model
                .faces_of_body(built.result)
                .iter()
                .all(|f| !built.model.face_origins(*f).is_empty())
        );
    }
}

#[test]
fn algebraic_short_head_retains_its_real_annular_floor() {
    for angle in ["82", "100"] {
        let recipe = recipe(
            json!({"blind":"1"}),
            json!({"cone":"118"}),
            Some(json!({"countersink":{"diameter":"8","angle":angle}})),
        );
        let (built, work) = catallax_numbers::with_operation_budget(100_000_000, || {
            evaluate(&recipe, &mut Decider::default())
        });
        let built = built.unwrap_or_else(|e| panic!("{angle} degree short head: {e}: {work:?}"));
        assert!(!work.exhausted);
        assert!(built.model.check(built.result).is_empty());
        let floors: Vec<_> = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .filter(|f| {
                built.model.face_origins(*f).iter().any(|o| {
                    matches!(
                        o.role,
                        Role::Hole {
                            part: HolePart::Floor,
                            ..
                        }
                    )
                })
            })
            .collect();
        assert_eq!(floors.len(), 1);
        assert_eq!(built.model.face(floors[0]).loops.len(), 2);
    }
}

#[test]
fn standard_angles_survive_scaled_rotated_nonunit_placements() {
    use catallax_boolean::{Containment, classify_point};
    use catallax_field::{Location, Plan, Workspace};
    use catallax_shapes::RPoint3;
    let q = |s: &str| Rational::parse(s).unwrap();
    for scale in ["1/1000", "1", "1000"] {
        let scale = q(scale);
        let scaled = |n: i64| (&scale * &Rational::from_i64(n)).to_string();
        for (tip, head) in [
            (json!({"cone":"118"}), None),
            (
                json!("flat"),
                Some(json!({"countersink":{"diameter":scaled(8),"angle":"82"}})),
            ),
            (
                json!("flat"),
                Some(json!({"countersink":{"diameter":scaled(8),"angle":"100"}})),
            ),
        ] {
            let original = recipe(json!({"blind":scaled(4)}), tip, head);
            let mut value: Value = serde_json::from_str(&original.to_text()).unwrap();
            value["steps"][0]["size"] = json!([scaled(10), scaled(10), scaled(10)]);
            // (x,y,z) -> (x,-z,y) + (7,11,-13), an exact rigid move.
            value["steps"][0]["pose"] = json!({"rotation":["1","1","0","0"],"translation":[(q("7")-&scale*&q("5")).to_string(),"11",(q("-13")-&scale*&q("5")).to_string()]});
            value["steps"][1]["diameter"] = json!(scaled(2));
            value["steps"][1]["at"][0]["point"] = json!(["7", "11", "-13"]);
            value["steps"][1]["at"][0]["direction"] = json!(["0", "-2", "0"]);
            let case = value["steps"][1].to_string();
            let definition = Recipe::parse(&value.to_string()).unwrap();
            let (built, work) = catallax_numbers::with_operation_budget(100_000_000, || {
                evaluate(&definition, &mut Decider::default())
            });
            let mut built =
                built.unwrap_or_else(|e| panic!("scale {scale}, {case}: {e}: {work:?}"));
            assert!(!work.exhausted);
            assert!(built.model.check(built.result).is_empty());
            assert!(built.model.body(built.result).marks.is_empty());
            let field = Plan::compile(&definition).unwrap();
            let mut workspace = Workspace::default();
            for (radius, material) in [("1/2", false), ("9/2", true)] {
                let point = RPoint3::new(
                    q("7") + &scale * &q(radius),
                    q("11") - &scale * &q("2"),
                    q("-13"),
                );
                assert_eq!(
                    classify_point(
                        &mut built.model,
                        built.result,
                        &point,
                        &mut Decider::default()
                    )
                    .proven(),
                    Some(&if material {
                        Containment::Inside
                    } else {
                        Containment::Outside
                    })
                );
                assert_eq!(
                    field.proven(&point, &mut workspace).location,
                    if material {
                        Location::Inside
                    } else {
                        Location::Outside
                    }
                );
            }
        }
    }
}
