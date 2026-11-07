//! Station closures preserve neighbours and change material only inside the authored run.

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
        .find(|case| case["id"] == "explicit-setbacks")
        .unwrap()["recipe"]
        .clone()
}

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    })
}

#[test]
fn stations_leave_both_end_walls_and_the_unrounded_edge_ends_intact() {
    let recipe = Recipe::parse(&input().to_string()).unwrap();
    let mut built = evaluate(&recipe, &mut decider()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    for (coordinates, expected) in [
        (["199/10", "1", "99/10"], Some(Containment::Inside)),
        (["199/10", "3", "99/10"], Some(Containment::Outside)),
        (["199/10", "19", "99/10"], Some(Containment::Inside)),
        (["96/5", "10", "46/5"], Some(Containment::Inside)),
        (["199/10", "2", "99/10"], None),
        (["199/10", "18", "99/10"], None),
        (["199/10", "0", "99/10"], None),
        (["199/10", "20", "99/10"], None),
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
    let mut corners: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|face| built.model.face_origins(face).iter())
        .filter_map(|origin| match &origin.role {
            Role::Corner { vertex } if origin.step.as_str() == "round" => Some(vertex.clone()),
            _ => None,
        })
        .collect();
    corners.sort_by(|a, b| a.faces().cmp(b.faces()));
    let mut expected: Vec<catallax_recipe::VertexName> =
        input()["steps"][1]["overflow"]["keep_neighbours"]
            .as_array()
            .unwrap()
            .iter()
            .map(|limit| serde_json::from_value(limit["at"].clone()).unwrap())
            .collect();
    expected.sort_by(|a, b| a.faces().cmp(b.faces()));
    assert_eq!(corners, expected);
}

#[test]
fn incompatible_stations_never_shrink_or_move_the_authored_run() {
    for (first, last) in [("10", "10"), ("11", "10")] {
        let mut value = input();
        let limits = &mut value["steps"][1]["overflow"]["keep_neighbours"];
        limits[0]["setback"] = json!(first);
        limits[1]["setback"] = json!(last);
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let BuildError::Geometry { refusal, .. } = evaluate(&recipe, &mut decider()).unwrap_err()
        else {
            panic!("expected an incompatible station refusal");
        };
        assert_eq!(refusal.code, "fillet.incompatible-setbacks");
    }
}

#[test]
fn conflicting_and_nonendpoint_limits_refuse_without_fallback() {
    for nonendpoint in [false, true] {
        let mut value = input();
        let limits = value["steps"][1]["overflow"]["keep_neighbours"]
            .as_array_mut()
            .unwrap();
        if nonendpoint {
            limits[0]["at"]["faces"][0]["role"] = json!("start");
        } else {
            let mut duplicate = limits[0].clone();
            duplicate["setback"] = json!("3");
            limits.push(duplicate);
        }
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let BuildError::Geometry { refusal, .. } = evaluate(&recipe, &mut decider()).unwrap_err()
        else {
            panic!("expected an incompatible station refusal");
        };
        assert_eq!(refusal.code, "fillet.incompatible-setbacks");
    }
}

#[test]
fn asymmetric_stations_follow_named_endpoints_in_either_record_order() {
    for reversed in [false, true] {
        let mut value = input();
        let limits = value["steps"][1]["overflow"]["keep_neighbours"]
            .as_array_mut()
            .unwrap();
        limits[0]["setback"] = json!("3");
        limits[1]["setback"] = json!("4");
        if reversed {
            limits.reverse();
        }
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let mut built = evaluate(&recipe, &mut decider()).unwrap();
        for (y, expected) in [
            ("2", Some(Containment::Inside)),
            ("3", None),
            ("4", Some(Containment::Outside)),
            ("15", Some(Containment::Outside)),
            ("16", None),
            ("17", Some(Containment::Inside)),
        ] {
            let point = RPoint3::new(
                Rational::parse("199/10").unwrap(),
                Rational::parse(y).unwrap(),
                Rational::parse("99/10").unwrap(),
            );
            let actual = classify_point(&mut built.model, built.result, &point, &mut decider());
            assert!(
                match (&actual, expected) {
                    (Outcome::Proven(found), Some(expected)) => *found == expected,
                    (Outcome::ProvenCoincident(_), None) => true,
                    _ => false,
                },
                "reversed={reversed}, y={y}: {actual:?}"
            );
        }
    }
}

#[test]
fn concave_station_regions_add_only_the_authored_interval() {
    let mut value: Value = include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|case| case["id"] == "l-concave")
        .unwrap()["recipe"]
        .clone();
    let edge = value["steps"][3]["edges"][0].clone();
    let adjacent = edge["between"].as_array().unwrap();
    let limits: Vec<_> = ["x_min", "x_max"]
        .map(|segment| {
            json!({
                "edge": edge,
                "at": {"faces": [adjacent[0], adjacent[1],
                    {"step": "base", "role": {"side": {"segment": segment}}}]},
                "setback": "2"
            })
        })
        .into_iter()
        .collect();
    value["steps"][3]["overflow"] = json!({"keep_neighbours": limits});
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let mut built = evaluate(&recipe, &mut decider()).unwrap();
    for (x, expected) in [
        ("1", Some(Containment::Outside)),
        ("2", None),
        ("3", Some(Containment::Inside)),
        ("18", None),
        ("19", Some(Containment::Outside)),
    ] {
        let point = RPoint3::new(
            Rational::parse(x).unwrap(),
            Rational::parse("51/10").unwrap(),
            Rational::parse("51/10").unwrap(),
        );
        let actual = classify_point(&mut built.model, built.result, &point, &mut decider());
        assert!(
            match (&actual, expected) {
                (Outcome::Proven(found), Some(expected)) => *found == expected,
                (Outcome::ProvenCoincident(_), None) => true,
                _ => false,
            },
            "x={x}: {actual:?}"
        );
    }
}
