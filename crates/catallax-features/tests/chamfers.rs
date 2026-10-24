//! Material, distance order and persistent names for exact bevels.

use catallax_features::{Built, evaluate};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{RPoint3, Surface};
use serde_json::{Value, json};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    })
}
fn build(recipe: &Recipe) -> Built {
    let built = evaluate(recipe, &mut decider()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    built
}
fn point(p: [&str; 3]) -> RPoint3 {
    let [x, y, z] = p.map(|x| Rational::parse(x).unwrap());
    RPoint3::new(x, y, z)
}
fn side(s: &str) -> Value {
    json!({"step":"base", "role":{"side":{"segment":s}}})
}
fn block_recipe(edges: Value, distance: &str, other: Option<&str>) -> Recipe {
    let mut bevel =
        json!({"id":"bevel", "op":"chamfer", "body":"base", "edges":edges, "distance":distance});
    if let Some(other) = other {
        bevel["other_distance"] = json!(other);
    }
    Recipe::parse(&json!({"format":"catallax-recipe/3", "name":"plane bevel", "steps":[
        {"id":"base", "op":"block", "size":["10","12","8"], "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}, bevel
    ], "result":"bevel"}).to_string()).unwrap()
}

#[test]
fn unequal_setbacks_follow_reference_order_and_have_exact_volume() {
    for reverse in [false, true] {
        let mut between = vec![side("x_min"), side("y_min")];
        if reverse {
            between.reverse();
        }
        let recipe = block_recipe(json!([{"between":between}]), "2", Some("3"));
        let mut built = build(&recipe);
        let v = catallax_measure::volume(
            &built.model,
            built.result,
            &Default::default(),
            &mut decider(),
        )
        .unwrap()
        .value;
        assert_eq!(v.lo, Rational::from_i64(936));
        assert_eq!(v.hi, v.lo);
        let (outside, inside) = if reverse {
            (["3/2", "1/2", "4"], ["5/2", "1/2", "4"])
        } else {
            (["1/2", "1", "4"], ["1/2", "2", "4"])
        };
        for (p, expected) in [
            (outside, catallax_boolean::Containment::Outside),
            (inside, catallax_boolean::Containment::Inside),
        ] {
            assert!(
                matches!(catallax_boolean::classify_point(&mut built.model, built.result, &point(p), &mut decider()), Outcome::Proven(actual) if actual == expected)
            );
        }
        assert_eq!(built.chamfers.len(), 1);
        assert!(matches!(built.chamfers[0].surface, Surface::Plane(_)));
    }
}

#[test]
fn p2_builds_cones_and_setback_edits_keep_all_origins() {
    let mut recipe = Recipe::parse(include_str!("../../../parts/P2-collar.json")).unwrap();
    let mut previous = None;
    for distance in ["1", "3/4"] {
        if let catallax_recipe::Operation::Chamfer { distance: d, .. } = &mut recipe.steps[2].op {
            *d = serde_json::from_value(json!(distance)).unwrap();
        }
        let built = build(&recipe);
        let faces = built.model.faces_of_body(built.result);
        assert_eq!(faces.len(), 6);
        assert_eq!(built.chamfers.len(), 2);
        let mut names = Vec::new();
        for face in faces {
            let origins = built.model.face_origins(face).to_vec();
            if origins.iter().any(|o| matches!(o.role, Role::Bevel { .. })) {
                assert!(matches!(
                    built.model.geom().surface(built.model.face(face).surface),
                    Surface::Cone(_)
                ));
            }
            names.push(origins);
        }
        names.sort();
        if let Some(previous) = previous {
            assert_eq!(names, previous);
        }
        previous = Some(names);
    }
}

#[test]
fn conflicting_or_ambiguous_orders_refuse_without_publishing() {
    let a = side("x_min");
    let b = side("y_min");
    for (edges, expected) in [
        (
            json!([{"between":[a,b]}, {"between":[b,a]}]),
            "chamfer.conflicting-distances",
        ),
        (
            json!([{"between":[{"step":"base"},{"step":"base"}]}]),
            "chamfer.ambiguous-order",
        ),
    ] {
        let recipe = block_recipe(edges, "2", Some("3"));
        let err = evaluate(&recipe, &mut decider()).unwrap_err();
        assert!(
            matches!(err, catallax_features::BuildError::Geometry { refusal, .. } if refusal.code == expected)
        );
    }
}

#[test]
fn nonpositive_distances_and_unproved_overflow_are_distinct() {
    for (distance, expected) in [
        ("0", "chamfer.non-positive-distance"),
        ("-1", "chamfer.non-positive-distance"),
        ("20", "chamfer.unsupported-overflow"),
    ] {
        let recipe = block_recipe(
            json!([{"between":[side("x_min"),side("y_min")]}]),
            distance,
            None,
        );
        let err = evaluate(&recipe, &mut decider()).unwrap_err();
        assert!(
            matches!(err, catallax_features::BuildError::Geometry { refusal, .. } if refusal.code == expected),
            "{distance}"
        );
    }
}

#[test]
fn four_top_strips_sew_their_shared_transversal_ends() {
    let top = json!({"step":"base","role":"end"});
    let edges: Vec<_> = ["x_min", "x_max", "y_min", "y_max"]
        .map(|s| json!({"between":[top,side(s)]}))
        .into();
    let recipe = block_recipe(json!(edges), "1", None);
    let built = build(&recipe);
    assert_eq!(built.chamfers.len(), 4);
    let v = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut decider(),
    )
    .unwrap()
    .value;
    assert_eq!(v.lo, Rational::parse("2818/3").unwrap());
    assert_eq!(v.hi, v.lo);
}

#[test]
fn all_four_collar_edges_have_disjoint_annular_regions() {
    let mut value: Value =
        serde_json::from_str(include_str!("../../../parts/P2-collar.json")).unwrap();
    value["steps"][2]["edges"] = json!("all");
    let built = build(&Recipe::parse(&value.to_string()).unwrap());
    assert_eq!(built.chamfers.len(), 4);
    assert_eq!(built.model.faces_of_body(built.result).len(), 8);
    let v = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut decider(),
    )
    .unwrap()
    .value;
    let expected = catallax_numbers::Ball::pi(256)
        * catallax_numbers::Ball::from_rational(&Rational::from_i64(1260), 256);
    assert!(v.lo <= expected.hi() && v.hi >= expected.lo());
}

#[test]
fn an_unequal_inner_rim_removes_material_and_preserves_the_opposite_end() {
    let pose = json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]});
    let recipe = Recipe::parse(&json!({"format":"catallax-recipe/3", "name":"inner rim", "steps":[
        {"id":"outer", "op":"cylinder", "radius":"8", "height":"9", "pose":pose},
        {"id":"inner", "op":"cylinder", "radius":"3", "height":"9", "pose":pose},
        {"id":"ring", "op":"boolean", "kind":"subtract", "target":"outer", "tool":"inner"},
        {"id":"bevel", "op":"chamfer", "body":"ring", "distance":"2", "other_distance":"1", "edges":[{"between":[
            {"step":"outer","role":"start"}, {"step":"inner","role":{"side":{"segment":"wall"}}}
        ]}]}
    ], "result":"bevel"}).to_string()).unwrap();
    let mut built = build(&recipe);
    let v = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut decider(),
    )
    .unwrap()
    .value;
    let expected = catallax_numbers::Ball::pi(256)
        * catallax_numbers::Ball::from_rational(&Rational::parse("1463/3").unwrap(), 256);
    assert!(v.lo <= expected.hi() && v.hi >= expected.lo());
    for (p, expected) in [
        (["7/2", "0", "1/2"], catallax_boolean::Containment::Outside),
        (["9/2", "0", "1/2"], catallax_boolean::Containment::Inside),
        (["7/2", "0", "17/2"], catallax_boolean::Containment::Inside),
    ] {
        assert!(
            matches!(catallax_boolean::classify_point(&mut built.model, built.result, &point(p), &mut decider()), Outcome::Proven(actual) if actual == expected)
        );
    }
}

#[test]
fn copied_and_reflected_bevels_keep_their_generating_definitions() {
    let recipe = Recipe::parse(&json!({"format":"catallax-recipe/3", "name":"placed bevels", "steps":[
        {"id":"base", "op":"cylinder", "radius":"8", "height":"12", "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
        {"id":"bevel", "op":"chamfer", "body":"base", "distance":"2", "edges":[{"between":[{"step":"base","role":"end"},side("wall")]}]},
        {"id":"copies", "op":"pattern", "body":"bevel", "layout":{"linear":{"vector":["30","0","0"], "count":2}}},
        {"id":"reflected", "op":"mirror", "body":"copies", "plane":{"point":["0","0","0"],"normal":["1","0","0"]}, "keep_original":false}
    ], "result":"reflected"}).to_string()).unwrap();
    let mut built = build(&recipe);
    let definitions: Vec<_> = built
        .chamfers
        .iter()
        .filter(|d| d.origin.step.as_str() == "reflected")
        .collect();
    assert_eq!(definitions.len(), 2);
    for definition in definitions {
        assert!(
            built
                .model
                .faces_of_body(built.result)
                .iter()
                .any(|f| built.model.face_origins(*f).contains(&definition.origin))
        );
    }
    assert!(matches!(
        catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &point(["-37", "0", "11"]),
            &mut decider()
        ),
        Outcome::ProvenCoincident(_)
    ));
    for (p, expected) in [
        (["-75/2", "0", "11"], catallax_boolean::Containment::Outside),
        (["-36", "0", "11"], catallax_boolean::Containment::Inside),
    ] {
        assert!(
            matches!(catallax_boolean::classify_point(&mut built.model, built.result, &point(p), &mut decider()), Outcome::Proven(actual) if actual == expected)
        );
    }
}
