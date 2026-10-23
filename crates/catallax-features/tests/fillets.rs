//! Elementary strip topology and radius-independent naming.

use catallax_features::evaluate;
use catallax_numbers::{BudgetPolicy, Decider};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::Surface;

#[test]
fn p3_builds_exact_named_strips_and_keeps_names_after_radius_edit() {
    let mut recipe = Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    let mut previous = None;
    for radius in ["1", "3/2"] {
        if let catallax_recipe::Operation::Fillet { radius: r, .. } = &mut recipe.steps[4].op {
            *r = serde_json::from_value(serde_json::json!(radius)).unwrap();
        }
        let mut dec = Decider::new(BudgetPolicy {
            max_work: 1_000_000,
            ..Default::default()
        });
        let built = evaluate(&recipe, &mut dec).unwrap();
        assert_eq!(built.model.check(built.result), vec![]);
        let mut names = Vec::new();
        let mut cylinders = 0;
        let mut tori = 0;
        for face in built.model.faces_of_body(built.result) {
            let origins = built.model.face_origins(face);
            for origin in origins {
                if matches!(origin.role, Role::Blend { .. }) {
                    names.push(origin.clone());
                    match built.model.geom().surface(built.model.face(face).surface) {
                        Surface::Cylinder(_) => cylinders += 1,
                        Surface::Torus(_) => tori += 1,
                        other => panic!("non-elementary fillet {other:?}"),
                    }
                }
            }
        }
        assert_eq!((cylinders, tori), (4, 1));
        names.sort();
        if let Some(previous) = previous {
            assert_eq!(names, previous);
        }
        previous = Some(names);
    }
}

#[test]
fn convex_circular_fillet_removes_only_the_selected_end_corner() {
    let recipe = Recipe::parse(r#"{
        "format":"catallax-recipe/3","name":"convex cylinder",
        "steps":[
            {"id":"base","op":"cylinder","radius":"8","height":"12",
             "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
            {"id":"round","op":"fillet","body":"base","radius":"2","edges":[
                {"between":[{"step":"base","role":"end"},{"step":"base","role":{"side":{"segment":"wall"}}}]}
            ]}
        ],"result":"round"
    }"#).unwrap();
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    });
    let mut built = evaluate(&recipe, &mut dec).unwrap();
    for (p, inside) in [
        (["38/5", "0", "58/5"], false),
        (["7", "0", "11"], true),
        (["39/5", "0", "1/5"], true),
    ] {
        let [x, y, z] = p.map(|v| catallax_numbers::Rational::parse(v).unwrap());
        let point = catallax_shapes::RPoint3::new(x, y, z);
        let expected = if inside {
            catallax_boolean::Containment::Inside
        } else {
            catallax_boolean::Containment::Outside
        };
        assert!(
            matches!(catallax_boolean::classify_point(&mut built.model,built.result,&point,&mut dec),
            catallax_numbers::Outcome::Proven(actual) if actual == expected),
            "at {point:?}"
        );
    }
}

#[test]
fn retained_blend_names_locate_the_translated_and_reflected_boundary() {
    let recipe = Recipe::parse(r#"{
        "format":"catallax-recipe/3","name":"placed fillet",
        "steps":[
            {"id":"base","op":"cylinder","radius":"8","height":"12",
             "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
            {"id":"round","op":"fillet","body":"base","radius":"2","edges":[
                {"between":[{"step":"base","role":"end"},{"step":"base","role":{"side":{"segment":"wall"}}}]}
            ]},
            {"id":"copies","op":"pattern","body":"round",
             "layout":{"linear":{"vector":["30","0","0"],"count":2}}},
            {"id":"reflected","op":"mirror","body":"copies",
             "plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}
        ],"result":"reflected"
    }"#).unwrap();
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    });
    let mut built = evaluate(&recipe, &mut dec).unwrap();
    let definition = built
        .fillets
        .iter()
        .find(|definition| {
            definition.origin.step.as_str() == "reflected"
                && matches!(&definition.surface, Surface::Torus(t)
                if t.center.x.as_rational() == Some(catallax_numbers::Rational::from_i64(-30)))
        })
        .unwrap();
    let face = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .find(|face| built.model.face_origins(*face).contains(&definition.origin))
        .unwrap();
    let point = |x: &str, z: &str| {
        catallax_shapes::RPoint3::new(
            catallax_numbers::Rational::parse(x).unwrap(),
            catallax_numbers::Rational::zero(),
            catallax_numbers::Rational::parse(z).unwrap(),
        )
    };
    assert!(catallax_boolean::face_contains(
        &mut built.model, face, &point("-186/5", "58/5"), &mut dec,
    ).unwrap());
    for (p, want) in [
        (point("-37", "11"), catallax_boolean::Containment::Inside),
        (
            point("-188/5", "58/5"),
            catallax_boolean::Containment::Outside,
        ),
        (point("7", "11"), catallax_boolean::Containment::Inside),
        (point("37", "11"), catallax_boolean::Containment::Outside),
    ] {
        assert!(
            matches!(
                catallax_boolean::classify_point(&mut built.model, built.result, &p, &mut dec),
                catallax_numbers::Outcome::Proven(actual) if actual == want
            ),
            "at {p:?}"
        );
    }
}
