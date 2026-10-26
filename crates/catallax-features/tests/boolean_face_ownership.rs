//! Distinct finite planar patches must not inherit each other's coplanar names.

use catallax_field::{FaceState, Location, Plan, Workspace};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{EdgeName, ElementId, Origin, Recipe, Role, StepId};
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation};
use serde_json::{Value, json};

const TWO_HOLES: &str = include_str!("fixtures/boolean-face-thicken/two-holes.json");

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn point(x: &str, y: &str) -> RPoint3 {
    RPoint3::new(q(x), q(y), q("7/2"))
}
fn rim(hole: &str, segment: &str) -> Origin {
    Origin {
        step: StepId::new("result").unwrap(),
        role: Role::Rim {
            edge: EdgeName::new(
                Origin {
                    step: StepId::new("stock").unwrap(),
                    role: Role::End,
                },
                Origin {
                    step: StepId::new(hole).unwrap(),
                    role: Role::Side {
                        segment: ElementId::new(segment).unwrap(),
                        path: None,
                    },
                },
            )
            .unwrap(),
        },
    }
}

fn check(recipe: &Recipe, map: &Pose, second_x: &str) {
    let (built, work) = with_operation_budget(100_000_000, || {
        catallax_features::evaluate(recipe, &mut Decider::default())
    });
    let mut built = built.unwrap();
    assert!(!work.exhausted);
    assert!(built.model.check(built.result).is_empty());
    assert_eq!(built.model.faces_of_body(built.result).len(), 14);
    let volume = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value
    .volume;
    assert_eq!(volume.lo, q("52"));
    assert_eq!(volume.lo, volume.hi);
    let (plan, work) = with_operation_budget(100_000_000, || Plan::compile(recipe));
    let plan = plan.unwrap();
    assert!(!work.exhausted);
    for (segment, y) in [("y_min", "2"), ("y_max", "4")] {
        let a = map.apply_rational_point(&point("3", y)).unwrap();
        let b = map.apply_rational_point(&point(second_x, y)).unwrap();
        for (hole, correct, other) in [("hole-a", &a, &b), ("hole-b", &b, &a)] {
            let origin = rim(hole, segment);
            let faces: Vec<_> = built
                .model
                .faces_of_body(built.result)
                .into_iter()
                .filter(|f| built.model.face_origins(*f).contains(&origin))
                .collect();
            assert_eq!(faces.len(), 1, "B-rep keeps one finite rim face per hole");
            assert!(
                catallax_boolean::face_contains(
                    &mut built.model,
                    faces[0],
                    correct,
                    &mut Decider::default()
                )
                .unwrap()
            );
            assert!(
                !catallax_boolean::face_contains(
                    &mut built.model,
                    faces[0],
                    other,
                    &mut Decider::default()
                )
                .unwrap()
            );
            let id = plan
                .select(&origin)
                .expect("the explicit rim origin must stay selectable");
            let positive = plan
                .face_query_checked(
                    id,
                    correct,
                    &mut Workspace::default(),
                    &mut Decider::default(),
                )
                .unwrap();
            let negative = plan
                .face_query_checked(
                    id,
                    other,
                    &mut Workspace::default(),
                    &mut Decider::default(),
                )
                .unwrap();
            println!(
                "{hole}.{segment}: correct={:?}, cross-hole={:?}",
                positive.state, negative.state
            );
            assert_eq!(
                positive.state,
                FaceState::OnFace,
                "correct-hole membership for {origin:?}"
            );
            assert_eq!(
                negative.state,
                FaceState::NotOnFace,
                "cross-hole origin leak for {origin:?} at {other:?}"
            );
        }
    }
}

#[test]
fn two_coplanar_hole_rims_keep_distinct_finite_domains() {
    check(&Recipe::parse(TWO_HOLES).unwrap(), &Pose::identity(), "7");
}

#[test]
fn finite_origin_restriction_preserves_two_hole_material_and_field_volume() {
    let recipe = Recipe::parse(TWO_HOLES).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    let (mass, work) = with_operation_budget(100_000_000, || {
        catallax_field::mass::properties(
            &plan,
            &Default::default(),
            &Default::default(),
            &mut Decider::default(),
        )
    });
    let mass = mass.unwrap();
    assert!(!work.exhausted);
    // Independent rectangle subtraction: (10 * 6 - 2 * 2 - 2 * 2) * 1.
    assert_eq!(mass.volume.lo, q("52"));
    assert_eq!(mass.volume.lo, mass.volume.hi);
    for (x, expected) in [
        ("1", Location::Inside),
        ("3", Location::Outside),
        ("5", Location::Inside),
        ("7", Location::Outside),
        ("9", Location::Inside),
    ] {
        assert_eq!(
            plan.proven_checked(
                &point(x, "3"),
                &mut Workspace::default(),
                &mut Decider::default(),
            )
            .unwrap()
            .location,
            expected,
        );
    }
}

#[test]
fn finite_rim_ownership_survives_a_rigid_source_transform() {
    let map = Pose::from_rational(
        Rotation::from_quaternion(2.into(), 1.into(), 0.into(), 0.into()).unwrap(),
        RVector3::new(5.into(), (-2).into(), 1.into()),
    );
    let mut value: Value = serde_json::from_str(TWO_HOLES).unwrap();
    for step in value["steps"].as_array_mut().unwrap() {
        if step["op"] != "block" {
            continue;
        }
        let t = step["pose"]["translation"].as_array().unwrap();
        let p = RPoint3::new(
            q(t[0].as_str().unwrap()),
            q(t[1].as_str().unwrap()),
            q(t[2].as_str().unwrap()),
        );
        let p = map.apply_rational_point(&p).unwrap();
        step["pose"] = json!({"rotation":["2","1","0","0"],"translation":[p.x.to_string(),p.y.to_string(),p.z.to_string()]});
    }
    check(&Recipe::parse(&value.to_string()).unwrap(), &map, "7");
}

#[test]
fn a_hole_position_edit_retains_origin_identity_without_cross_selection() {
    let original = Recipe::parse(TWO_HOLES).unwrap();
    let checkpoint =
        catallax_features::evaluate_incremental(&original, None, &mut Decider::default()).unwrap();
    let mut value: Value = serde_json::from_str(TWO_HOLES).unwrap();
    value["steps"][3]["pose"]["translation"][0] = json!("5");
    let edited = Recipe::parse(&value.to_string()).unwrap();
    let resumed =
        catallax_features::evaluate_incremental(&edited, Some(checkpoint), &mut Decider::default())
            .unwrap();
    assert!(
        resumed
            .built()
            .model
            .check(resumed.built().result)
            .is_empty()
    );
    check(&edited, &Pose::identity(), "6");
}
