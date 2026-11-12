//! Complete authored collar: independent meridian moments and one reported mass.
use catallax_app::assembly;
use catallax_features::evaluate;
use catallax_field::{
    Plan,
    assembly::{self as field, PlacedPart},
    mass,
};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::Pose;
use std::num::NonZeroUsize;

#[test]
fn p2_has_independent_mass_and_centroid_agreement_at_default_accuracy() {
    let recipe = Recipe::parse(include_str!("../../../parts/P2-collar.json")).unwrap();
    let (built, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
        evaluate(&recipe, &mut Decider::default())
    });
    assert!(!work.exhausted);
    let built = built.unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    let part = PlacedPart::new(&plan, &Pose::identity(), &mut Decider::default()).unwrap();
    let density = Rational::parse("0.00785").unwrap();
    let accuracy = mass::Accuracy::default();
    assert_eq!(
        field::uniform_mass(
            &part,
            &density,
            &accuracy,
            &mass::Options::default(),
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "field.print-coverage"
    );
    let mut baseline = None;
    for threads in [1, 2, 4] {
        catallax_numbers::set_threads(NonZeroUsize::new(threads));
        let (independent, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
            field::revolved_uniform_mass(
                &part,
                &density,
                &accuracy,
                &mass::Options::default(),
                &mut Decider::default(),
            )
        });
        let independent = independent.unwrap();
        assert!(!work.exhausted);
        let brep =
            assembly::uniform_mass(&built, &density, &accuracy, &mut Decider::default()).unwrap();
        let reported = assembly::mass_agreement(&brep, &independent, &density, &accuracy).unwrap();
        assert!(reported.value.grams.lo > Rational::parse("31.3858").unwrap());
        assert!(reported.value.grams.hi < Rational::parse("31.3859").unwrap());
        assert_eq!(
            independent
                .centre_of_gravity
                .as_ref()
                .unwrap()
                .clone()
                .map(|v| v.lo),
            [Rational::zero(), Rational::from_i64(6), Rational::zero()]
        );
        let snapshot = (
            independent.grams,
            independent.centre_of_gravity,
            reported.value.grams,
            reported.value.centre_of_gravity,
            work.used,
        );
        if let Some(expected) = &baseline {
            assert_eq!(&snapshot, expected);
        } else {
            baseline = Some(snapshot);
        }
        println!(
            "P2 threads={threads} independent field mass work={}",
            work.used
        );
    }
    catallax_numbers::set_threads(None);
}
