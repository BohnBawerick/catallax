//! Physical openings, names and material after a planar side draft.
use catallax_boolean::{Containment, classify_point};
use catallax_features::evaluate;
use catallax_numbers::{Decider, Outcome, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::Point3;
use serde_json::{Value, json};

#[test]
fn side_draft_preserves_circular_and_polygonal_through_openings() {
    let source: Value =
        serde_json::from_str(include_str!("fixtures/draft-mixed-openings.json")).unwrap();
    for angle in ["3", "5", "-3"] {
        let mut source = source.clone();
        source["steps"][5]["angle"] = json!(angle);
        let recipe = Recipe::parse(&source.to_string()).unwrap();
        let (mut built, work) = with_operation_budget(20_000_000, || {
            evaluate(&recipe, &mut Decider::default()).unwrap()
        });
        eprintln!("angle={angle} {work:?}");
        assert!(!work.exhausted);
        assert!(built.model.check(built.result).is_empty());
        assert_eq!(built.model.faces_of_body(built.result).len(), 13);
        for (point, wanted) in [
            (["10", "15", "6"], Containment::Outside),
            (["28", "15", "6"], Containment::Outside),
            (["20", "15", "6"], Containment::Inside),
            (["1", "15", "6"], Containment::Inside),
        ] {
            let [x, y, z] = point.map(|v| Rational::parse(v).unwrap().into());
            assert_eq!(
                classify_point(
                    &mut built.model,
                    built.result,
                    &Point3::new(x, y, z),
                    &mut Decider::default()
                ),
                Outcome::Proven(wanted)
            );
        }
        let mut source = source.clone();
        source["steps"].as_array_mut().unwrap().pop();
        source["result"] = json!("pierced");
        let original = evaluate(
            &Recipe::parse(&source.to_string()).unwrap(),
            &mut Decider::default(),
        )
        .unwrap();
        let names = |built: &catallax_features::Built| {
            built
                .model
                .faces_of_body(built.result)
                .iter()
                .map(|face| built.model.face_origins(*face).to_vec())
                .collect::<std::collections::BTreeSet<_>>()
        };
        assert_eq!(names(&built), names(&original));
    }
}
