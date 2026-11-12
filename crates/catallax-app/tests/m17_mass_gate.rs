//! Full authored targets at unchanged default accuracies and application work.
//! Run explicitly in release mode; these are verdicts, not substitute recipes.

use catallax_app::assembly;
use catallax_features::evaluate;
use catallax_field::{
    Plan,
    assembly::{self as field, Mass, PlacedPart},
    mass,
};
use catallax_measure::Enclosure;
use catallax_numbers::{Answer, Decider, Rational, set_threads, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::Pose;
use std::num::NonZeroUsize;

const PARTS: [(&str, &str); 6] = [
    ("P1", include_str!("../../../parts/P1-bracket.json")),
    ("P2", include_str!("../../../parts/P2-collar.json")),
    ("P3", include_str!("../../../parts/P3-filleted-boss.json")),
    ("P4", include_str!("../../../parts/P4-shelled-tray.json")),
    ("P5", include_str!("../../../parts/P5-flanged-housing.json")),
    ("P6", include_str!("../../../parts/P6-lattice-bracket.json")),
];

fn bounds(e: &Enclosure) -> String {
    format!(
        "[{}, {}]",
        e.lo.to_decimal_bounds(10).0,
        e.hi.to_decimal_bounds(10).1
    )
}
fn properties(m: &Mass) -> String {
    format!(
        "grams={}; centre={:?}",
        bounds(&m.grams),
        m.centre_of_gravity
            .as_ref()
            .map(|c| c.each_ref().map(bounds))
    )
}

fn verdict(name: &str, text: &str) -> String {
    let recipe = Recipe::parse(text).unwrap();
    let density = Rational::parse("0.00785").unwrap();
    let accuracy = mass::Accuracy::default();
    let (built, construction) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
        evaluate(&recipe, &mut Decider::default())
    });
    let (brep, measurement) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
        built.as_ref().map_err(|e| e.to_string()).and_then(|built| {
            assembly::uniform_mass(built, &density, &accuracy, &mut Decider::default())
                .map_err(|e| e.to_string())
        })
    });
    let (independent, field_work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
        Plan::compile(&recipe)
            .and_then(|p| PlacedPart::new(&p, &Pose::identity(), &mut Decider::default()))
            .and_then(|part| {
                field::uniform_mass(
                    &part,
                    &density,
                    &accuracy,
                    &mass::Options::default(),
                    &mut Decider::default(),
                )
            })
            .map_err(|e| e.to_string())
    });
    let brep_record = match &brep {
        Ok(m) => properties(&m.value),
        Err(r) => r.clone(),
    };
    let field_record = match &independent {
        Ok(m) => properties(m),
        Err(r) => r.clone(),
    };
    let agreed = match (&brep, &independent) {
        (Ok(a), Ok(b)) => assembly::mass_agreement(a, b, &density, &accuracy),
        _ => Err(catallax_numbers::Refusal::new(
            "assembly.mass-unresolved",
            "at least one independent mass route refused",
        )),
    };
    let agreement = match agreed {
        Ok(Answer { value, .. }) => properties(&value),
        Err(r) => r.to_string(),
    };
    format!(
        "{name}: build-work={}; brep-work={}; field-work={}\nB-rep: {brep_record}\nfield: {field_record}\nagreement: {agreement}",
        construction.used, measurement.used, field_work.used
    )
}

#[test]
#[ignore = "full target mass verdicts; run in release mode on the work PC"]
fn target_mass_verdicts_are_identical_at_one_two_four_threads() {
    let mut baseline = None;
    for threads in [1, 2, 4] {
        set_threads(NonZeroUsize::new(threads));
        let rows: Vec<_> = PARTS
            .iter()
            .map(|(name, text)| {
                let row = verdict(name, text);
                println!("threads={threads}\n{row}");
                row
            })
            .collect();
        if let Some(baseline) = &baseline {
            assert_eq!(&rows, baseline);
        } else {
            baseline = Some(rows);
        }
    }
    set_threads(None);
}
