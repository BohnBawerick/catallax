//! A failed independent query does not erase another proven mass.
use catallax_app::assembly::{self, Mass, MassSource};
use catallax_field::mass;
use catallax_measure::Enclosure;
use catallax_numbers::{Answer, NotExact, Rational, Refusal};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn value(grams: &str) -> Mass {
    Mass {
        grams: Enclosure {
            lo: q(grams),
            hi: q(grams),
        },
        centre_of_gravity: Some(std::array::from_fn(|_| Enclosure {
            lo: q("1"),
            hi: q("1"),
        })),
    }
}
#[test]
fn one_source_keeps_the_other_refusal_and_two_sources_require_agreement() {
    let brep = Answer::exact(value("6"));
    let field = value("6");
    let failed = Refusal::new("numbers.operation-budget", "19999999 construction work");
    let a = mass::Accuracy::default();
    let result = assembly::resolve_mass(Err(&failed), Ok(&field), &q("1"), &a).unwrap();
    assert_eq!(result.value.source, MassSource::Field);
    assert_eq!(result.value.mass.grams, field.grams);
    assert_eq!(result.value.unavailable.unwrap().detail, failed.detail);
    let result = assembly::resolve_mass(Ok(&brep), Err(&failed), &q("1"), &a).unwrap();
    assert_eq!(result.value.source, MassSource::Brep);
    assert_eq!(result.value.unavailable.unwrap().code, failed.code);
    let result = assembly::resolve_mass(Ok(&brep), Ok(&field), &q("1"), &a).unwrap();
    assert_eq!(result.value.source, MassSource::Agreed);
    assert!(result.value.unavailable.is_none());
    assert_eq!(
        assembly::resolve_mass(Ok(&brep), Ok(&value("7")), &q("1"), &a)
            .unwrap_err()
            .code,
        "assembly.disagreement"
    );
    let r = assembly::resolve_mass(Err(&failed), Err(&failed), &q("1"), &a).unwrap_err();
    assert_eq!(r.code, "assembly.mass-unresolved");
    assert!(r.detail.contains("19999999"));
}
#[test]
fn invalid_density_and_unqualified_single_source_bounds_still_refuse() {
    let a = mass::Accuracy::default();
    let failed = Refusal::new("field.print-coverage", "not integrated");
    let mut field = value("1");
    field.grams.hi = q("2");
    assert_eq!(
        assembly::resolve_mass(Err(&failed), Ok(&field), &q("1"), &a)
            .unwrap_err()
            .code,
        "assembly.agreement-width"
    );
    assert_eq!(
        assembly::resolve_mass(Err(&failed), Ok(&value("1")), &q("0"), &a)
            .unwrap_err()
            .code,
        "assembly.density"
    );
}

#[test]
fn numeric_undecided_bounds_survive_reporting() {
    let decision = NotExact {
        bound: q("1/1000"),
        code: "numbers.budget-exhausted",
        detail: "bounded undecided predicate".into(),
    };
    let failed = Refusal::new("app.example", "could not finish").keeping(decision.clone());
    let other = Refusal::new("field.example", "uncovered definition");
    let a = mass::Accuracy::default();
    let result = assembly::resolve_mass(Err(&failed), Ok(&value("1")), &q("1"), &a).unwrap();
    assert_eq!(
        result.value.unavailable.unwrap().undecided(),
        Some(&decision)
    );
    let result = assembly::resolve_mass(Err(&failed), Err(&other), &q("1"), &a).unwrap_err();
    assert_eq!(result.undecided(), Some(&decision));
}
