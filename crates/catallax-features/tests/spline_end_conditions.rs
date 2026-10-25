//! A represented spline is not an implicit quadric or an exact intersection pair.
use catallax_features::end_condition::curved::AuthoredCap;
use catallax_numbers::{Decider, Real};
use catallax_recipe::{Recipe, StepId};
use catallax_shapes::{GeomStore, Plane, Point3, Surface, Vector3};

#[test]
fn represented_spline_limit_refuses_until_an_exact_pair_is_available() {
    let recipe = Recipe::parse(include_str!("fixtures/curved-limits/spline-offset.json")).unwrap();
    let mut dec = Decider::default();
    let built = catallax_features::evaluate_certified(&recipe, &mut dec).unwrap();
    let spline = built
        .result_certificate()
        .spline()
        .expect("genuine represented spline test part");
    let mut geom = GeomStore::new();
    let surface = geom.add_surface(Surface::BSpline(spline.clone().into()));
    let plane = geom.add_surface(Surface::Plane(
        Plane::new(
            Point3::new(Real::zero(), Real::zero(), Real::from_i64(4)),
            Vector3::new(Real::zero(), Real::zero(), Real::one()),
            Vector3::new(Real::one(), Real::zero(), Real::zero()),
            &mut dec,
        )
        .unwrap(),
    ));
    let refusal =
        catallax_intersect::intersect_surfaces(&mut geom, plane, surface, &mut dec).unwrap_err();
    assert_eq!(refusal.code, "intersect.unsupported-pair");
    let refusal = AuthoredCap::compile(
        &recipe,
        &StepId::new("target").unwrap(),
        &StepId::new("profile").unwrap(),
        None,
        None,
        &mut dec,
    )
    .unwrap_err();
    assert_eq!(refusal.code, "features.unsupported-curved-limit");
    assert!(built.built().model.check(built.built().result).is_empty());
}
