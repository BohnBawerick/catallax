//! Owner pocket geometry retained through numerical source embedding.

use catallax_numbers::{BudgetPolicy, Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::Surface;
use catallax_topology::numerical::{NumericalOptions, build_with_native_trims};

#[test]
fn owner_source_preserves_open_pocket_in_numerical_embedding() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let recipe = Recipe::parse(include_str!("fixtures/owner-curved-pocket.json")).unwrap();
        let mut dec = Decider::new(BudgetPolicy {
            max_work: 100_000_000,
            ..Default::default()
        });
        let built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
        assert_eq!(built.model.counts(built.result).faces, 13);
        let tolerance = Rational::parse("0.0000001").unwrap();
        let (definition, native) =
            catallax_features::numerical_embedding::definition_from_exact_model(
                &built.model,
                built.result,
                recipe.name.clone(),
                &tolerance,
                &mut dec,
            )
            .unwrap();
        let body = build_with_native_trims(
            definition,
            native,
            NumericalOptions {
                tolerance_mm: tolerance,
            },
            &mut dec,
        )
        .unwrap();
        assert_eq!(body.definition().faces.len(), 13);
        let top: Vec<_> = body
            .definition()
            .faces
            .iter()
            .filter(|face| {
                matches!(&face.surface, Surface::Plane(plane)
                if plane.origin.z.as_rational() == Some(Rational::from_i64(60))
                    && plane.normal.x.as_rational().is_some_and(|q| q.is_zero())
                    && plane.normal.y.as_rational().is_some_and(|q| q.is_zero()))
            })
            .collect();
        assert_eq!(
            top.len(),
            1,
            "one flat top, without a coincident downward lid"
        );
        // The plane chart is orthonormal; signed loop area therefore has mm² units.
        // The opening removes 17 * 22.67 from the 30 * 76 flat top.
        let twice_area = top[0]
            .loops
            .iter()
            .flatten()
            .fold(Rational::zero(), |sum, trim| {
                sum + &trim.uv[0][0] * &trim.uv[1][1] - &trim.uv[1][0] * &trim.uv[0][1]
            });
        let area = twice_area.abs().mul_pow2(-1);
        assert!(
            (&area - &Rational::parse("1894.61").unwrap()).abs()
                <= Rational::parse("0.000001").unwrap(),
            "flat top area={area} mm²"
        );
    });
    println!(
        "correct owner source and native embedding request work={}",
        work.used
    );
}
