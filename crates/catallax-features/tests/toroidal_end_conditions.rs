//! First-hit torus caps must stop at the curved branch, not its centre plane.
use catallax_features::{
    BuildError,
    end_condition::{curved::Branch, toroidal::AnnularCap},
};
use catallax_numbers::{Decider, Outcome, Rational};
use catallax_recipe::{Recipe, Role, StepId};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

fn value() -> Value {
    serde_json::from_str(include_str!("fixtures/curved-limits/torus-entry-next.json")).unwrap()
}
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn recipe(value: Value) -> Recipe {
    Recipe::parse(&value.to_string()).unwrap()
}

#[test]
fn entry_exit_selected_wall_and_reversed_scaled_rays_stop_at_exact_torus_sections() {
    for (op, start, direction, branch, inside, outside) in [
        ("extrude_next", "0", "1", Branch::Entry, "3", "9/2"),
        ("extrude_to", "0", "2", Branch::Entry, "3", "9/2"),
        ("extrude_next", "5", "3", Branch::Exit, "11/2", "7"),
        ("extrude_to", "5", "1", Branch::Exit, "11/2", "7"),
        ("extrude_next", "10", "-2", Branch::Entry, "7", "11/2"),
        ("extrude_next", "5", "-1", Branch::Exit, "9/2", "3"),
    ] {
        let mut value = value();
        value["steps"][0]["plane"]["translation"][2] = json!(start);
        value["steps"][2]["op"] = json!(op);
        value["steps"][2]["direction"] = json!(["0", "0", direction]);
        if op == "extrude_to" {
            value["steps"][2]["limit"] =
                json!({"face":{"step":"target","role":{"side":{"segment":"wall"}}}});
        }
        let recipe = recipe(value);
        let mut dec = Decider::default();
        let cap = AnnularCap::compile(
            &recipe,
            &StepId::new("target").unwrap(),
            &StepId::new("profile").unwrap(),
            match &recipe.steps[2].op {
                catallax_recipe::Operation::ExtrudeTo { direction, .. }
                | catallax_recipe::Operation::ExtrudeNext { direction, .. } => direction.as_ref(),
                _ => unreachable!(),
            },
            None,
            &mut dec,
        )
        .unwrap();
        assert_eq!(cap.branch(), branch);
        let mut built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
        assert!(built.model.check(built.result).is_empty());
        assert!(built.model.body(built.result).marks.is_empty());
        let (off_inside, off_outside) = match (start, direction) {
            ("0", _) => ("401/100", "41/10"),
            ("10", _) => ("599/100", "59/10"),
            ("5", "-1") => ("41/10", "401/100"),
            _ => ("59/10", "599/100"),
        };
        for (radial, z, expected) in [
            ("3", inside, catallax_boolean::Containment::Inside),
            ("3", outside, catallax_boolean::Containment::Outside),
            ("0", inside, catallax_boolean::Containment::Outside),
            ("4", inside, catallax_boolean::Containment::Outside),
            ("13/4", off_inside, catallax_boolean::Containment::Inside),
            ("13/4", off_outside, catallax_boolean::Containment::Outside),
        ] {
            let point = RPoint3::new(q(radial), q("0"), q(z));
            assert!(
                matches!(catallax_boolean::classify_point(&mut built.model, built.result, &point, &mut dec), Outcome::Proven(found) if found == expected),
                "{op} {start} {direction} at {radial},{z}"
            );
        }
        let names: Vec<_> = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .flat_map(|f| built.model.face_origins(f))
            .collect();
        assert!(names.iter().all(|o| o.step.as_str() == "cap"));
        for role in [Role::Start, Role::End] {
            assert!(names.iter().any(|o| o.role == role));
        }
        for segment in ["inner", "outer"] {
            assert!(names.iter().any(
                |o| matches!(&o.role, Role::Side { segment: id, .. } if id.as_str() == segment)
            ));
        }
    }
}

#[test]
fn unsupported_torus_pairs_refuse_without_claiming_a_missing_hit() {
    for change in [0, 1, 2, 3] {
        let mut value = value();
        match change {
            0 => value["steps"][2]["direction"] = json!(["1", "0", "1"]),
            1 => value["steps"][0]["plane"]["translation"][0] = json!("1/4"),
            2 => value["steps"][0]["loops"][0][0]["radius"] = json!("5"),
            _ => value["steps"][0]["plane"]["translation"][2] = json!("41/10"),
        }
        let BuildError::Geometry { refusal, .. } =
            catallax_features::evaluate(&recipe(value), &mut Decider::default()).unwrap_err()
        else {
            panic!("geometric coverage refusal")
        };
        assert_eq!(refusal.code, "features.unsupported-curved-limit");
    }
}

#[test]
fn zero_work_cannot_publish_a_toroidal_branch() {
    let recipe = recipe(value());
    let (result, work) = catallax_numbers::with_operation_budget(0, || {
        AnnularCap::compile(
            &recipe,
            &StepId::new("target").unwrap(),
            &StepId::new("profile").unwrap(),
            None,
            None,
            &mut Decider::default(),
        )
    });
    assert!(work.exhausted);
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
}

#[test]
fn rigidly_placed_torus_and_profile_keep_the_same_physical_cap() {
    let mut value = value();
    value["steps"][0]["plane"] =
        json!({"rotation":["1","1","0","0"],"translation":["7","-11","3"]});
    value["steps"][1]["pose"] = json!({"rotation":["1","1","0","0"],"translation":["7","-16","3"]});
    let mut built = catallax_features::evaluate(&recipe(value), &mut Decider::default()).unwrap();
    for (y, expected) in [
        ("-14", catallax_boolean::Containment::Inside),
        ("-31/2", catallax_boolean::Containment::Outside),
    ] {
        assert!(
            matches!(catallax_boolean::classify_point(&mut built.model, built.result,
            &RPoint3::new(q("10"),q(y),q("3")), &mut Decider::default()), Outcome::Proven(found) if found == expected)
        );
    }
    assert!(built.model.check(built.result).is_empty());
}

#[test]
fn selected_torus_wall_cannot_touch_a_nonexistent_shared_edge() {
    let mut value = value();
    value["steps"][2]["op"] = json!("extrude_to");
    value["steps"][2]["limit"] = json!({
        "face":{"step":"target","role":{"side":{"segment":"wall"}}},
        "touching":[{"step":"target","role":{"side":{"segment":"wall"}}}]
    });
    let BuildError::Geometry { refusal, .. } =
        catallax_features::evaluate(&recipe(value), &mut Decider::default()).unwrap_err()
    else {
        panic!("missing reference")
    };
    assert_eq!(refusal.code, "recipe.missing-reference");
}
