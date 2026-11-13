//! Complete P1 mass from the independent field when B-rep construction refuses.
use catallax_field::{
    Plan,
    assembly::{self, PlacedPart},
    mass,
};
use catallax_numbers::{Ball, Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::Pose;
use std::num::NonZeroUsize;

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn affine(pi: &Ball, constant: &str, coefficient: &str) -> Ball {
    &Ball::from_rational(&q(constant), 512) + &(pi * &Ball::from_rational(&q(coefficient), 512))
}
#[test]
fn authored_p1_mass_and_centroid_match_the_exact_component_formula() {
    let recipe = Recipe::parse(include_str!("../../../parts/P1-bracket.json")).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    let part = PlacedPart::new(&plan, &Pose::identity(), &mut Decider::default()).unwrap();
    let density = q("0.00785");
    assert_eq!(
        assembly::uniform_mass(
            &part,
            &density,
            &mass::Accuracy::default(),
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
        let (result, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
            assembly::prismatic_uniform_mass(
                &part,
                &density,
                &mass::Accuracy::default(),
                &mass::Options::default(),
                &mut Decider::default(),
            )
        });
        let m = result.unwrap();
        assert!(!work.exhausted);
        let pi = Ball::pi(512);
        let volume = affine(&pi, "19080", "275");
        let grams = &volume * &Ball::from_rational(&density, 512);
        assert!(m.grams.lo <= grams.hi() && grams.lo() <= m.grams.hi);
        let centre = m.centre_of_gravity.unwrap();
        assert!(centre[0].contains(&q("30")));
        for (i, numerator) in [
            (1, affine(&pi, "257700", "8127.5")),
            (2, affine(&pi, "153690", "3590")),
        ] {
            let expected = numerator.checked_div(&volume).unwrap();
            assert!(centre[i].lo <= expected.hi() && expected.lo() <= centre[i].hi);
        }
        let snapshot = (m.grams.clone(), centre.clone(), work.used);
        if let Some(expected) = &baseline {
            assert_eq!(&snapshot, expected);
        } else {
            baseline = Some(snapshot);
        }
        println!(
            "P1 native mass work={}; grams=[{},{}]; centre={:?}",
            work.used,
            m.grams.lo.to_decimal_bounds(16).0,
            m.grams.hi.to_decimal_bounds(16).1,
            centre.map(|v| (v.lo.to_decimal_bounds(16).0, v.hi.to_decimal_bounds(16).1))
        );
    }
    catallax_numbers::set_threads(None);
}
