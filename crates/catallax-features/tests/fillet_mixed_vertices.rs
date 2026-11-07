//! Convex fillets roll over the concave blends they meet at a mixed vertex.

use catallax_boolean::{Containment, classify_point};
use catallax_features::{BuildError, evaluate};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{RPoint3, Surface};
use serde_json::{Value, json};

fn input() -> Value {
    include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|case| case["id"] == "l-mixed-all")
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
fn every_edge_of_the_l_rounds_with_a_torus_at_each_mixed_vertex() {
    let recipe = Recipe::parse(&input().to_string()).unwrap();
    let mut built = evaluate(&recipe, &mut decider()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    // The concave blend's axis is y = z = 6; the end-face rounds have
    // centres at x = 1, and the torus centre circle has radius 2 about it.
    for (coordinates, expected) in [
        (["10", "26/5", "26/5"], Some(Containment::Inside)),
        (["10", "59/10", "59/10"], Some(Containment::Outside)),
        (["10", "6", "5"], None),
        (["10", "5", "6"], None),
        (["1/2", "24/5", "24/5"], Some(Containment::Inside)),
        (["1/100", "24/5", "24/5"], Some(Containment::Outside)),
        (["0", "9/2", "9/2"], None),
        (["39/2", "24/5", "24/5"], Some(Containment::Inside)),
        (["1999/100", "24/5", "24/5"], Some(Containment::Outside)),
        (["10", "1", "1"], Some(Containment::Inside)),
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
    assert_eq!(faces.len(), 38);
    let mut corners = Vec::new();
    let mut tori = 0;
    for face in &faces {
        let origins = built.model.face_origins(*face);
        let torus = matches!(
            built.model.geom().surface(built.model.face(*face).surface),
            Surface::Torus(_)
        );
        for origin in origins.iter().filter(|o| o.step.as_str() == "round") {
            match &origin.role {
                Role::Corner { vertex } => corners.push(vertex.clone()),
                // A blend over a same-step blend must be renamed as its corner.
                Role::Blend { edge } => assert!(
                    edge.between().iter().all(|o| o.step.as_str() != "round"),
                    "{edge:?}"
                ),
                role => panic!("unexpected role {role:?}"),
            }
        }
        if torus {
            tori += 1;
            assert!(
                origins
                    .iter()
                    .all(|o| matches!(o.role, Role::Corner { .. })),
                "{origins:?}"
            );
        }
    }
    assert_eq!(tori, 2);
    // Ten spherical corners and two mixed junctions: one per input vertex.
    corners.sort_by(|a, b| a.faces().cmp(b.faces()));
    corners.dedup();
    assert_eq!(corners.len(), 12);
}

#[test]
fn a_mixed_vertex_with_an_unselected_edge_keeps_its_refusal() {
    let mut value = input();
    let between = |a: Value, b: Value| json!({"between": [a, b]});
    let end = json!({"step": "base", "role": "end"});
    let side = json!({"step": "wall", "role": {"side": {"segment": "y_max"}}});
    let x_min = json!({"step": "base", "role": {"side": {"segment": "x_min"}}});
    value["steps"][3]["edges"] = json!([between(end.clone(), side), between(end, x_min)]);
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let error = evaluate(&recipe, &mut decider()).unwrap_err();
    let BuildError::Geometry { refusal, .. } = error else {
        panic!("expected a construction refusal, got {error:?}");
    };
    assert!(
        refusal.code.starts_with("fillet.unsupported"),
        "{refusal:?}"
    );
}
