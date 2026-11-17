//! Numerical admission keeps an atlas face's trims on that face's own charts.

use catallax_features::numerical_embedding::definition_from_exact_model;
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_shapes::{RationalAtlas, SplineSupport, Surface};
use catallax_topology::numerical::{NumericalOptions, build_with_native_trims};

#[test]
fn a_trim_on_another_component_of_the_atlas_source_refuses() {
    let recipe = catallax_recipe::Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/m10-smooth-loft.json"
    ))
    .unwrap();
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
        let tolerance = Rational::from_ratio(1, 10_000_000).unwrap();
        let (mut definition, native_trims) = definition_from_exact_model(
            &built.model,
            built.result,
            "smooth loft".into(),
            &tolerance,
            &mut dec,
        )
        .unwrap();
        let options = NumericalOptions {
            tolerance_mm: tolerance,
        };
        build_with_native_trims(
            definition.clone(),
            native_trims.clone(),
            options.clone(),
            &mut dec,
        )
        .unwrap();
        // The closed loft profile gives one source with two components: the
        // roof face and the planar floor. Give the roof's trims the floor's
        // charts; every sample still lies on the shared source surface.
        let face = definition
            .faces
            .iter_mut()
            .find(|face| {
                matches!(
                    face.surface,
                    Surface::BSpline(SplineSupport::RationalAtlas(_))
                )
            })
            .unwrap();
        let Surface::BSpline(SplineSupport::RationalAtlas(roof)) = &face.surface else {
            unreachable!()
        };
        let roof_domains: Vec<_> = roof
            .charts()
            .iter()
            .map(|c| c.patch().domain().clone())
            .collect();
        let other = RationalAtlas::partition(
            roof.source().clone(),
            [roof.periodic_u(), roof.periodic_v()],
            &mut dec,
        )
        .unwrap()
        .into_iter()
        .find(|atlas| {
            atlas
                .charts()
                .iter()
                .map(|c| c.patch().domain())
                .ne(roof_domains.iter())
        })
        .unwrap();
        face.surface = Surface::BSpline(other.into());
        let error =
            build_with_native_trims(definition, native_trims, options, &mut dec).unwrap_err();
        assert_eq!(error.code, "numerical.invalid-trim");
        assert_eq!(
            error.detail,
            "parameter trim leaves the charts of its atlas face"
        );
    });
    assert!(!work.exhausted, "{work:?}");
}
