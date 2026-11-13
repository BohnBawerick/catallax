//! Complete authored tray checked against a separate circular-segment formula.
use catallax_field::{
    Plan,
    assembly::{self, PlacedPart},
    mass,
};
use catallax_numbers::{Ball, Decider, Rational, set_threads, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::Pose;
use std::num::NonZeroUsize;

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn b(s: &str) -> Ball {
    Ball::from_rational(&q(s), 512)
}

#[test]
fn authored_p4_uses_its_complete_shell_definition() {
    let recipe = Recipe::parse(include_str!("../../../parts/P4-shelled-tray.json")).unwrap();
    // Outer rectangle clipped by r23; the 36x26 inner rectangle is wholly inside.
    let root129 = b("129").sqrt().unwrap();
    let root304 = b("304").sqrt().unwrap();
    let dec = Decider::default();
    let angle = &b("20")
        .checked_div(&root129)
        .unwrap()
        .atan_checked(&dec)
        .unwrap()
        + &b("15")
            .checked_div(&root304)
            .unwrap()
            .atan_checked(&dec)
            .unwrap();
    let angle = &angle - &(&Ball::pi(512) * &b("0.5"));
    let area = &b("2") * &(&(&b("15") * &root304) + &(&b("20") * &root129) + (&b("529") * &angle));
    let volume = &(&b("12") * &area) - &b("9360");
    let first_z = &(&b("72") * &area) - &b("65520");
    let expected_z = first_z.checked_div(&volume).unwrap();
    let expected_mass = &volume * &b("0.00785");
    let mut baseline = None;
    for threads in [1, 2, 4] {
        set_threads(NonZeroUsize::new(threads));
        let (result, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
            let plan = Plan::compile(&recipe)?;
            let part = PlacedPart::new(&plan, &Pose::identity(), &mut Decider::default())?;
            assembly::prismatic_uniform_mass(
                &part,
                &q("0.00785"),
                &mass::Accuracy::default(),
                &mass::Options::default(),
                &mut Decider::default(),
            )
        });
        let result = result.unwrap();
        assert!(!work.exhausted);
        assert!(result.grams.lo <= expected_mass.hi() && expected_mass.lo() <= result.grams.hi);
        let centre = result.centre_of_gravity.unwrap();
        assert!(centre[0].contains(&q("0")) && centre[1].contains(&q("0")));
        assert!(centre[2].lo <= expected_z.hi() && expected_z.lo() <= centre[2].hi);
        let snapshot = (result.grams.clone(), centre.clone(), work.used);
        if let Some(expected) = &baseline {
            assert_eq!(&snapshot, expected);
        } else {
            baseline = Some(snapshot);
        }
        println!(
            "P4 field work={}; grams=[{},{}]; centre={:?}",
            work.used,
            result.grams.lo.to_decimal_bounds(16).0,
            result.grams.hi.to_decimal_bounds(16).1,
            centre.map(|v| (v.lo.to_decimal_bounds(16).0, v.hi.to_decimal_bounds(16).1))
        );
    }
    set_threads(None);
}
