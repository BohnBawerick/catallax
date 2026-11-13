//! Full authored targets at unchanged default accuracies and application work.
//! Run explicitly in release mode; these are verdicts, not substitute recipes.

use catallax_app::assembly;
use catallax_features::{BuildError, evaluate};
use catallax_field::{
    Plan,
    assembly::{self as field, Mass, PlacedPart},
    mass,
};
use catallax_measure::Enclosure;
use catallax_numbers::{Decider, Rational, Refusal, set_threads, with_operation_budget};
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
        e.lo.to_decimal_bounds(16).0,
        e.hi.to_decimal_bounds(16).1
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

fn build_refusal(error: &BuildError) -> Refusal {
    match error {
        BuildError::Geometry { step, refusal } => {
            let mut refusal = refusal.clone();
            refusal.detail = format!("step {step}: {}", refusal.detail);
            refusal
        }
        BuildError::Recipe(error) => Refusal::new("recipe.invalid", error.to_string()),
    }
}

fn verdict(name: &str, text: &str) -> (String, String) {
    let recipe = Recipe::parse(text).unwrap();
    let density = Rational::parse("0.00785").unwrap();
    let accuracy = mass::Accuracy::default();
    let (built, construction) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
        evaluate(&recipe, &mut Decider::default())
    });
    let (brep, measurement) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
        built.as_ref().map_err(build_refusal).and_then(|built| {
            assembly::uniform_mass(built, &density, &accuracy, &mut Decider::default()).map_err(
                |e| {
                    let mut refusal = e.refusal;
                    refusal.detail = format!("{}; operands: {:?}", refusal.detail, e.operands);
                    refusal
                },
            )
        })
    });
    let (independent, field_work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
        Plan::compile(&recipe)
            .and_then(|p| PlacedPart::new(&p, &Pose::identity(), &mut Decider::default()))
            .and_then(|part| {
                let revolved = field::revolved_uniform_mass(
                    &part,
                    &density,
                    &accuracy,
                    &mass::Options::default(),
                    &mut Decider::default(),
                );
                if !matches!(&revolved, Err(r) if r.code == "field.mass-revolution-coverage") {
                    return revolved;
                }
                let prismatic = field::prismatic_uniform_mass(
                    &part,
                    &density,
                    &accuracy,
                    &mass::Options::default(),
                    &mut Decider::default(),
                );
                if !matches!(&prismatic, Err(r) if r.code == "field.mass-prismatic-coverage") {
                    return prismatic;
                }
                field::uniform_mass(
                    &part,
                    &density,
                    &accuracy,
                    &mass::Options::default(),
                    &mut Decider::default(),
                )
            })
    });
    for (stage, succeeded, exhausted) in [
        ("construction", built.is_ok(), construction.exhausted),
        ("B-rep mass", brep.is_ok(), measurement.exhausted),
        ("field mass", independent.is_ok(), field_work.exhausted),
    ] {
        assert!(
            !succeeded || !exhausted,
            "{name}: successful {stage} exhausted work"
        );
    }
    let brep_record = match &brep {
        Ok(m) => properties(&m.value),
        Err(r) => r.to_string(),
    };
    let field_record = match &independent {
        Ok(m) => properties(m),
        Err(r) => r.to_string(),
    };
    let resolved = assembly::resolve_mass(brep.as_ref(), independent.as_ref(), &density, &accuracy);
    if name == "P5" {
        assert_eq!(
            resolved.as_ref().unwrap_err().code,
            "assembly.mass-unresolved"
        );
        assert_eq!(
            independent.as_ref().unwrap_err().code,
            "fillet.unsupported-overflow"
        );
        assert_eq!(brep.as_ref().unwrap_err().code, "assembly.mass-accuracy");
    } else {
        let resolved = resolved.as_ref().unwrap_or_else(|r| panic!("{name}: {r}"));
        let expected_source = match name {
            "P1" | "P6" => assembly::MassSource::Field,
            "P2" | "P4" => assembly::MassSource::Agreed,
            "P3" => assembly::MassSource::Brep,
            _ => unreachable!(),
        };
        assert_eq!(resolved.value.source, expected_source, "{name}");
    }
    let canonical = format!(
        "{brep:?}; {independent:?}; {resolved:?}; undecided={:?}/{:?}; work={}/{}/{}; exhausted={}/{}/{}",
        brep.as_ref().err().and_then(Refusal::undecided),
        independent.as_ref().err().and_then(Refusal::undecided),
        construction.used,
        measurement.used,
        field_work.used,
        construction.exhausted,
        measurement.exhausted,
        field_work.exhausted,
    );
    let reported = match resolved {
        Ok(answer) => format!(
            "{:?}: {}",
            answer.value.source,
            properties(&answer.value.mass)
        ),
        Err(r) => r.to_string(),
    };
    let report = format!(
        "{name}: build-work={}; brep-work={}; field-work={}; exhausted={}/{}/{}\nB-rep: {brep_record}\nfield: {field_record}\nreported: {reported}",
        construction.used,
        measurement.used,
        field_work.used,
        construction.exhausted,
        measurement.exhausted,
        field_work.exhausted
    );
    (canonical, report)
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
                let (canonical, report) = verdict(name, text);
                println!("threads={threads}\n{report}");
                canonical
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
