//! A roll-over fillet whose contact reaches the far edge of its neighbour
//! consumes that face and reports it.

use catallax_boolean::{Containment, classify_point};
use catallax_features::{BuildError, evaluate};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

fn input() -> Value {
    include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|case| case["id"] == "at-face-limit")
        .unwrap()["recipe"]
        .clone()
}

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    })
}

fn side(segment: &str) -> Value {
    json!({"step": "b", "role": {"side": {"segment": segment}}})
}

#[test]
fn a_ball_as_tall_as_the_side_consumes_it_and_keeps_every_other_face() {
    let recipe = Recipe::parse(&input().to_string()).unwrap();
    let mut built = evaluate(&recipe, &mut decider()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    // The cylinder's axis is the bottom-face line x = 10, z = 0.
    for (coordinates, expected) in [
        (["15", "10", "5"], Some(Containment::Inside)),
        (["5", "10", "99/10"], Some(Containment::Inside)),
        (["11", "10", "99/10"], Some(Containment::Inside)),
        (["19", "10", "9"], Some(Containment::Outside)),
        (["20", "10", "5"], Some(Containment::Outside)),
        (["16", "10", "8"], None),
        (["20", "10", "0"], None),
        (["5", "0", "5"], None),
        (["10", "20", "5"], None),
    ] {
        let [x, y, z] = coordinates.map(|n| Rational::parse(n).unwrap());
        let point = RPoint3::new(x, y, z);
        let actual = classify_point(&mut built.model, built.result, &point, &mut decider());
        assert!(
            match (&actual, expected) {
                (Outcome::Proven(found), Some(expected)) => *found == expected,
                (Outcome::ProvenCoincident(_), None) => true,
                _ => false,
            },
            "at {coordinates:?}: {actual:?}"
        );
    }
    let faces = built.model.faces_of_body(built.result);
    assert_eq!(faces.len(), 6);
    let blends: Vec<_> = faces
        .iter()
        .flat_map(|face| built.model.face_origins(*face).iter())
        .filter(|origin| matches!(origin.role, Role::Blend { .. }))
        .cloned()
        .collect();
    assert_eq!(blends.len(), 1);
    let consumed = serde_json::to_value(
        built
            .consumed
            .iter()
            .map(|record| (record.step.as_str(), &record.source, &record.by))
            .collect::<Vec<_>>(),
    )
    .unwrap();
    assert_eq!(
        consumed,
        json!([[
            "round",
            [side("x_max")],
            serde_json::to_value(&blends).unwrap()
        ]])
    );
}

#[test]
fn a_radius_past_the_side_keeps_its_limit_refusal() {
    let mut value = input();
    value["steps"][1]["radius"] = json!("10001/1000");
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let BuildError::Geometry { refusal, .. } = evaluate(&recipe, &mut decider()).unwrap_err()
    else {
        panic!("expected a construction refusal");
    };
    assert_eq!(refusal.code, "fillet.radius-too-large");
}

#[test]
fn a_shorter_side_reached_only_part_way_keeps_its_refusal() {
    // A step face below the side cuts the contact band, so the foot leaves the
    // side before the full radius and needs continuation not built here.
    let mut value = input();
    let steps = value["steps"].as_array_mut().unwrap();
    steps.insert(
        1,
        json!({"id": "cut", "op": "block", "size": ["4", "20", "4"],
            "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["18", "0", "0"]}}),
    );
    steps.insert(
        2,
        json!({"id": "notched", "op": "boolean", "kind": "subtract", "target": "b", "tool": "cut"}),
    );
    steps[3]["body"] = json!("notched");
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let error = evaluate(&recipe, &mut decider()).unwrap_err();
    let BuildError::Geometry { refusal, .. } = error else {
        panic!("expected a construction refusal, got {error:?}");
    };
    assert!(refusal.code.starts_with("fillet."), "{refusal:?}");
}

#[test]
fn a_concave_ball_as_tall_as_the_wall_side_fills_up_to_the_wall_top() {
    let mut value: Value = include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|case| case["id"] == "l-concave")
        .unwrap()["recipe"]
        .clone();
    // The wall side spans z from 5 to 10, so radius 5 reaches its top edge.
    value["steps"][3]["radius"] = json!("5");
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let mut built = evaluate(&recipe, &mut decider()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    // The blend circle is centred on y = z = 10.
    for (coordinates, expected) in [
        (["10", "6", "11/2"], Some(Containment::Inside)),
        (["10", "4", "9"], Some(Containment::Inside)),
        (["10", "7", "8"], Some(Containment::Outside)),
        (["10", "11/2", "19/2"], Some(Containment::Outside)),
        (["10", "10", "5"], None),
        (["10", "5", "10"], None),
        (["0", "6", "11/2"], None),
        (["20", "6", "11/2"], None),
    ] {
        let [x, y, z] = coordinates.map(|n| Rational::parse(n).unwrap());
        let point = RPoint3::new(x, y, z);
        let actual = classify_point(&mut built.model, built.result, &point, &mut decider());
        assert!(
            match (&actual, expected) {
                (Outcome::Proven(found), Some(expected)) => *found == expected,
                (Outcome::ProvenCoincident(_), None) => true,
                _ => false,
            },
            "at {coordinates:?}: {actual:?}"
        );
    }
    let [record] = built.consumed.as_slice() else {
        panic!("expected one consumed face, got {:?}", built.consumed);
    };
    assert_eq!(
        serde_json::to_value(&record.source).unwrap(),
        json!([{"step": "wall", "role": {"side": {"segment": "y_max"}}}])
    );
}
