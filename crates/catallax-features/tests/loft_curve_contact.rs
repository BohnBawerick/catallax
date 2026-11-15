//! Finite line contact and nearby crossings on the recorded ruled loft.

use catallax_boolean::{Containment, classify_point};
use catallax_numbers::{Decider, Outcome, Rational, Real, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::Point3;
use catallax_topology::BodyKind;
use serde_json::{Value, json};

#[test]
fn a_finite_tangent_interval_preserves_material_and_a_nearby_crossing_cuts_it() {
    for (height, removed) in [
        ("2", false),
        ("2000001/1000000", false),
        ("1999999/1000000", true),
    ] {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut input: Value =
                serde_json::from_str(include_str!("fixtures/loft-spline-prism.json")).unwrap();
            let steps = input["steps"].as_array_mut().unwrap();
            steps.push(json!({"id":"box","op":"block","size":["10","10","1"],
                "pose":{"rotation":["1","0","0","0"],"translation":["-3",height,"1"]}}));
            steps.push(json!({"id":"cut","op":"boolean","kind":"subtract",
                "target":"loft","tool":"box"}));
            input["result"] = json!("cut");
            let recipe = Recipe::parse(&input.to_string()).unwrap();
            let mut dec = Decider::default();
            let mut built = catallax_features::evaluate(&recipe, &mut dec)
                .unwrap_or_else(|error| panic!("height={height}: {error:?}"));
            assert_eq!(built.model.body(built.result).kind, BodyKind::Solid);
            assert_eq!(built.model.body(built.result).lumps.len(), 1);
            assert!(built.model.check(built.result).is_empty());
            for (y, expected) in [
                ("1", Containment::Inside),
                (
                    "19999995/10000000",
                    if removed {
                        Containment::Outside
                    } else {
                        Containment::Inside
                    },
                ),
            ] {
                let point = Point3::new(
                    Real::from_i64(2),
                    Rational::parse(y).unwrap().into(),
                    Rational::parse("3/2").unwrap().into(),
                );
                assert!(
                    matches!(classify_point(&mut built.model, built.result, &point, &mut dec),
                    Outcome::Proven(actual) if actual == expected),
                    "height={height}, y={y}"
                );
            }
        });
        assert!(!work.exhausted, "height={height}: {} units", work.used);
    }
}
