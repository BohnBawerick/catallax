//! Numerical source qualification for the recorded loft top intersection.

use catallax_boolean::{BoolOp, numerical::loft::loft_box_selection};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{Recipe, StepId};
use catallax_topology::{BodyKind, numerical};
use serde_json::{Value, json};

#[test]
#[ignore = "slow loft-top admission, over 20M; phase-5 speed list"]
fn recorded_top_intersection_admits_a_numerical_source_body() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut input: Value = serde_json::from_str(include_str!(
            "../../catallax-features/tests/fixtures/m10-smooth-loft.json"
        ))
        .unwrap();
        input["steps"].as_array_mut().unwrap().push(json!({
            "id":"box","op":"block","size":["10","10","1"],
            "pose":{"rotation":["1","0","0","0"],"translation":["-3","5/2","1/2"]}
        }));
        let recipe = Recipe::parse(&input.to_string()).unwrap();
        let mut dec = Decider::default();
        let built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
        let loft = built.bodies[&StepId::new("loft").unwrap()];
        let stock = built.bodies[&StepId::new("box").unwrap()];
        let tolerance = Rational::parse("1/10000000").unwrap();
        let selected = loft_box_selection(
            &built.model,
            loft,
            stock,
            BoolOp::Intersect,
            tolerance.clone(),
            &mut dec,
        )
        .unwrap();
        assert_eq!(selected.classifications.len(), 15);
        for class in &selected.classifications {
            let [x, y, z] = class
                .witness_bounds
                .each_ref()
                .map(|bound| bound.mid().to_f64_nearest());
            let expected = if class.from_loft {
                x > -3. && x < 7. && y > 2.5 && y < 12.5 && z > 0.5 && z < 1.5
            } else {
                // Independent natural-cubic loft law from the recorded profiles.
                let v = z.min(2. - z);
                let h = 4. + 3. * v - v * v * v;
                let roof = 2. * h * (x / 4.) * (1. - x / 4.);
                x > 0. && x < 4. && y > 0. && z > 0. && z < 2. && y < roof
            };
            assert_eq!(class.inside, expected, "material witness {:?}", class.face);
            assert!(!class.evidence.is_empty());
        }
        let body = selected
            .admit_with(&mut dec, |model, collection, options, dec| {
                assert_eq!(model.body(collection).kind, BodyKind::Surface);
                let counts = model.counts(collection);
                assert_eq!((counts.faces, counts.edges, counts.vertices), (4, 6, 4));
                let (definition, trims) =
                    catallax_features::numerical_embedding::definition_from_exact_model(
                        model,
                        collection,
                        recipe.name.clone(),
                        &options.tolerance_mm,
                        dec,
                    )?;
                numerical::build_with_native_trims(definition, trims, options, dec)
            })
            .unwrap();
        let counts = body.counts();
        assert_eq!(
            (
                counts.lumps,
                counts.shells,
                counts.faces,
                counts.edges,
                counts.vertices
            ),
            (1, 1, 4, 6, 4)
        );
        assert_eq!(body.evidence().tolerance_mm(), &tolerance);
        assert!(body.evidence().repairs().is_empty());
    });
    assert!(work.used > 0 && work.used <= 20_000_000);
    assert!(!work.exhausted);
}
