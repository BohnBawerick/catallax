//! Closed rim collapse preserves material, source history and general placement.

use catallax_boolean::{Containment, classify_point};
use catallax_features::{evaluate, evaluate_incremental, fillet::FilletDefinition};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, Real};
use catallax_recipe::{Origin, Recipe};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    })
}

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn assert_real(actual: &Real, expected: &Rational) {
    let difference = actual - &Real::from(expected.clone());
    assert!(
        matches!(decider().sign(&difference), Outcome::ProvenCoincident(_)),
        "actual={actual:?}, expected={expected}"
    );
}

fn assert_declared_request(definition: &FilletDefinition, requested: &Rational, limit: &Rational) {
    if requested > limit {
        assert_real(
            definition
                .requested_radius
                .as_ref()
                .expect("terminal request"),
            requested,
        );
    } else {
        assert!(definition.requested_radius.is_none());
    }
}

#[test]
fn a_closed_rim_becomes_the_correct_hemisphere_and_records_only_consumed_faces() {
    for (radius, height, requested, cap, placed) in [
        ("5", "10", "5", "end", false),
        ("5", "10", "5", "start", true),
        ("5", "5", "5", "end", false),
        ("3/2", "3/2", "3/2", "start", true),
        ("5", "10", "501/100", "end", false),
        ("5", "5", "100", "start", true),
    ] {
        let mut value: Value = include_str!("../../../parts/fillets/cases.jsonl")
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).unwrap())
            .find(|case| case["id"] == "cylinder-collapse")
            .unwrap()["recipe"]
            .clone();
        value["steps"][0]["radius"] = json!(radius);
        value["steps"][0]["height"] = json!(height);
        value["steps"][1]["radius"] = json!(requested);
        value["steps"][1]["edges"][0]["between"][0]["role"] = json!(cap);
        if placed {
            value["steps"][0]["pose"] =
                json!({"rotation":["1","1","0","0"],"translation":["3","4","5"]});
        }
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let mut built = evaluate(&recipe, &mut decider()).unwrap();
        assert!(built.model.check(built.result).is_empty());
        let r = q(radius);
        let h = q(height);
        let centre = &h - &r;
        assert_real(&built.fillets[0].radius, &r);
        assert_declared_request(&built.fillets[0], &q(requested), &r);
        let probes = [
            (
                Rational::zero(),
                &centre + &(&r * &q("49/50")),
                Some(Containment::Inside),
            ),
            (
                &r * &q("49/50"),
                &centre + &(&r * &q("49/50")),
                Some(Containment::Outside),
            ),
            (&r * &q("3/5"), &centre + &(&r * &q("4/5")), None),
            (&r * &q("3/5"), h.clone(), Some(Containment::Outside)),
            (&r * &q("49/50"), &r * &q("1/10"), Some(Containment::Inside)),
            (Rational::zero(), Rational::zero(), None),
            (r.clone(), centre.clone(), None),
            (Rational::zero(), h.clone(), None),
            (
                Rational::zero(),
                -&(&r * &q("1/10")),
                Some(Containment::Outside),
            ),
        ];
        for (mut x, mut z, expected) in probes {
            if cap == "start" {
                z = &h - &z;
            }
            let (y, z) = if placed {
                x = &x + &q("3");
                (&q("4") - &z, q("5"))
            } else {
                (Rational::zero(), z)
            };
            let point = RPoint3::new(x, y, z);
            let actual = classify_point(&mut built.model, built.result, &point, &mut decider());
            assert!(
                match (&actual, expected) {
                    (Outcome::Proven(found), Some(expected)) => *found == expected,
                    (Outcome::ProvenCoincident(_), None) => true,
                    _ => false,
                },
                "r={radius}, h={height}, cap={cap}, placed={placed}, at {point:?}: {actual:?}"
            );
        }
        let cap_origin: Origin = serde_json::from_value(json!({"step":"c","role":cap})).unwrap();
        let wall_origin: Origin =
            serde_json::from_value(json!({"step":"c","role":{"side":{"segment":"wall"}}})).unwrap();
        let blend: Origin = serde_json::from_value(json!({"step":"round","role":{"blend":{"edge":
            {"between":value["steps"][1]["edges"][0]["between"]}}}}))
        .unwrap();
        let mut expected = vec![vec![cap_origin.clone()]];
        if h == r {
            expected.push(vec![wall_origin.clone()]);
        }
        expected.sort();
        let mut actual: Vec<_> = built
            .consumed
            .iter()
            .filter(|record| record.step.as_str() == "round")
            .map(|record| {
                assert_eq!(record.by, vec![blend.clone()]);
                record.source.clone()
            })
            .collect();
        actual.sort();
        assert_eq!(actual, expected);
        let surviving: Vec<_> = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .flat_map(|face| built.model.face_origins(face).iter())
            .collect();
        assert!(!surviving.contains(&&cap_origin));
        assert_eq!(surviving.contains(&&wall_origin), h != r);
    }
}

#[test]
fn radius_edits_restore_the_cap_and_clear_terminal_history() {
    let mut value: Value = include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|case| case["id"] == "cylinder-collapse")
        .unwrap()["recipe"]
        .clone();
    let cap: Origin = serde_json::from_value(json!({"step":"c","role":"end"})).unwrap();
    let mut prefix = None;
    for requested in ["9/2", "5", "501/100", "9/2"] {
        value["steps"][1]["radius"] = json!(requested);
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let checkpoint = evaluate_incremental(&recipe, prefix.take(), &mut decider()).unwrap();
        let mut built = checkpoint.built().clone();
        let collapsed = q(requested) >= q("5");
        assert_eq!(
            built
                .model
                .faces_of_body(built.result)
                .iter()
                .any(|face| built.model.face_origins(*face).contains(&cap)),
            !collapsed,
        );
        assert_eq!(built.consumed.is_empty(), !collapsed);
        assert_declared_request(&built.fillets[0], &q(requested), &q("5"));
        let point = RPoint3::new(q("1/4"), Rational::zero(), q("10"));
        let actual = classify_point(&mut built.model, built.result, &point, &mut decider());
        assert!(
            matches!(
                (&actual, collapsed),
                (Outcome::Proven(Containment::Outside), true)
                    | (Outcome::ProvenCoincident(_), false)
            ),
            "radius={requested}: {actual:?}"
        );
        prefix = Some(checkpoint);
    }
}

#[test]
fn thin_stock_uses_its_first_wall_event_and_keeps_the_original_cap() {
    for (height, requested, cap, placed) in [
        ("3", "3", "end", false),
        ("3", "4", "start", true),
        ("3", "100", "end", false),
        ("2", "3", "end", false),
        ("3", "5/2", "end", false),
        ("3", "3/2", "start", true),
    ] {
        let mut value: Value = include_str!("../../../parts/fillets/cases.jsonl")
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).unwrap())
            .find(|case| case["id"] == "cylinder-collapse")
            .unwrap()["recipe"]
            .clone();
        value["steps"][0]["height"] = json!(height);
        value["steps"][1]["radius"] = json!(requested);
        value["steps"][1]["edges"][0]["between"][0]["role"] = json!(cap);
        if placed {
            value["steps"][0]["pose"] =
                json!({"rotation":["1","1","0","0"],"translation":["3","4","5"]});
        }
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let mut built = evaluate(&recipe, &mut decider()).unwrap_or_else(|error| {
            panic!("height={height}, request={requested}, cap={cap}, placed={placed}: {error:?}")
        });
        let r = q("5");
        let h = q(height);
        let request = q(requested);
        let rho = std::cmp::min(&request, &h).clone();
        let major = &r - &rho;
        assert_real(&built.fillets[0].radius, &rho);
        assert_declared_request(&built.fillets[0], &request, &h);
        let cap_origin: Origin = serde_json::from_value(json!({"step":"c","role":cap})).unwrap();
        let wall_origin: Origin =
            serde_json::from_value(json!({"step":"c","role":{"side":{"segment":"wall"}}})).unwrap();
        let faces = built.model.faces_of_body(built.result);
        assert!(
            faces
                .iter()
                .any(|face| built.model.face_origins(*face).contains(&cap_origin))
        );
        let wall_gone = rho == h;
        assert_eq!(
            faces
                .iter()
                .any(|face| built.model.face_origins(*face).contains(&wall_origin)),
            !wall_gone
        );
        if wall_gone {
            assert_eq!(built.consumed.len(), 1);
            assert_eq!(built.consumed[0].source, vec![wall_origin]);
            assert_eq!(built.consumed[0].by, vec![built.fillets[0].origin.clone()]);
        } else {
            assert!(built.consumed.is_empty());
        }
        for (x, mut z, expected) in [
            (q("0"), &h / &q("2"), Some(Containment::Inside)),
            (q("0"), &h - &(&rho / &q("10")), Some(Containment::Inside)),
            (
                &major / &q("2"),
                &h - &(&rho / &q("2")),
                Some(Containment::Inside),
            ),
            (r.clone(), h.clone(), Some(Containment::Outside)),
            (
                r.clone(),
                &h - &(&rho / &q("2")),
                Some(Containment::Outside),
            ),
            (&major + &(&rho * &q("3/5")), &h - &(&rho / &q("5")), None),
            (q("0"), h.clone(), None),
            (r.clone(), q("0"), None),
        ] {
            if cap == "start" {
                z = &h - &z;
            }
            let point = if placed {
                RPoint3::new(&x + &q("3"), &q("4") - &z, q("5"))
            } else {
                RPoint3::new(x, q("0"), z)
            };
            let actual = classify_point(&mut built.model, built.result, &point, &mut decider());
            assert!(
                match &actual {
                    Outcome::Proven(value) => Some(*value) == expected,
                    Outcome::ProvenCoincident(_) => expected.is_none(),
                    _ => false,
                },
                "height={height}, request={requested}, cap={cap}, point={point:?}: {actual:?}"
            );
        }
    }
}

#[test]
fn a_disconnected_component_keeps_its_material_and_names() {
    let mut value: Value = include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|case| case["id"] == "cylinder-collapse")
        .unwrap()["recipe"]
        .clone();
    let steps = value["steps"].as_array_mut().unwrap();
    let mut round = steps.pop().unwrap();
    steps.push(json!({"id":"other","op":"cylinder","radius":"2","height":"4","pose":{"rotation":["1","0","0","0"],"translation":["20","0","0"]}}));
    steps.push(json!({"id":"j","op":"boolean","kind":"union","target":"c","tool":"other"}));
    round["body"] = json!("j");
    steps.push(round);
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let mut built = evaluate(&recipe, &mut decider()).unwrap();
    assert_eq!(built.model.body(built.result).lumps.len(), 2);
    let other_cap: Origin = serde_json::from_value(json!({"step":"other","role":"end"})).unwrap();
    assert!(
        built
            .model
            .faces_of_body(built.result)
            .iter()
            .any(|face| built.model.face_origins(*face).contains(&other_cap))
    );
    assert_eq!(built.consumed.len(), 1);
    assert_eq!(
        built.consumed[0].source,
        vec![serde_json::from_value::<Origin>(json!({"step":"c","role":"end"})).unwrap()]
    );
    for (point, expected) in [
        (["20", "0", "2"], Some(Containment::Inside)),
        (["10", "0", "2"], Some(Containment::Outside)),
        (["20", "0", "4"], None),
        (["22", "0", "2"], None),
    ] {
        let point = RPoint3::new(q(point[0]), q(point[1]), q(point[2]));
        let actual = classify_point(&mut built.model, built.result, &point, &mut decider());
        assert!(
            match &actual {
                Outcome::Proven(value) => Some(*value) == expected,
                Outcome::ProvenCoincident(_) => expected.is_none(),
                _ => false,
            },
            "{point:?}: {actual:?}"
        );
    }
}

#[test]
fn an_algebraic_axis_node_of_the_carrier_is_inside_the_trimmed_solid() {
    let mut value: Value = include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|case| case["id"] == "cylinder-collapse")
        .unwrap()["recipe"]
        .clone();
    value["steps"][0]["height"] = json!("3");
    value["steps"][1]["radius"] = json!("3");
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let mut built = evaluate(&recipe, &mut decider()).unwrap();
    let point = catallax_shapes::Point3::new(
        Real::from_i64(0),
        Real::from_i64(0),
        Real::from_i64(5).sqrt(),
    );
    let actual = classify_point(&mut built.model, built.result, &point, &mut decider());
    assert!(
        matches!(actual, Outcome::Proven(Containment::Inside)),
        "{actual:?}"
    );
}
