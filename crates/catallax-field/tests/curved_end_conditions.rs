//! The field's support-sign query is independent of height evaluation and B-rep.

use catallax_features::end_condition::curved::{AuthoredCap, Location};
use catallax_field::CurvedField;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Recipe, Role, StepId};
use catallax_shapes::RPoint3;

fn id(s: &str) -> StepId {
    StepId::new(s).unwrap()
}
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn compile() -> (AuthoredCap, CurvedField) {
    let recipe = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/curved-limits/sphere-source.json"
    ))
    .unwrap();
    let cap = AuthoredCap::compile(
        &recipe,
        &id("target"),
        &id("profile"),
        None,
        None,
        &mut Decider::default(),
    )
    .unwrap();
    let field = CurvedField::new(cap.clone(), id("capped"));
    (cap, field)
}

#[test]
fn support_signs_match_the_closed_first_hit_set_and_authored_names() {
    let (cap, field) = compile();
    let mut dec = Decider::default();
    for (p, expected) in [
        (["0", "0", "2"], Location::Inside),
        (["0", "0", "4"], Location::Outside),
        (["0", "0", "3"], Location::Boundary),
        (["1/2", "0", "2"], Location::Boundary),
        (["3/4", "0", "2"], Location::Outside),
        (["0", "0", "0"], Location::Boundary),
        (["0", "0", "7"], Location::Outside),
    ] {
        let [x, y, z] = p.map(q);
        let point = RPoint3::new(x.clone(), y.clone(), z.clone());
        assert_eq!(field.location(&point, &mut dec).unwrap(), expected);
        assert_eq!(
            cap.definition()
                .location(&[x.into(), y.into()], &z.into(), &mut dec)
                .unwrap(),
            expected
        );
    }
    let start = field
        .origins(&RPoint3::new(q("0"), q("0"), q("0")), &mut dec)
        .unwrap();
    assert_eq!(start.len(), 1);
    assert_eq!(start[0].role, Role::Start);
    let end = field
        .origins(&RPoint3::new(q("0"), q("0"), q("3")), &mut dec)
        .unwrap();
    assert_eq!(end.len(), 1);
    assert_eq!(end[0].role, Role::End);
    let side = field
        .origins(&RPoint3::new(q("1/2"), q("0"), q("2")), &mut dec)
        .unwrap();
    assert_eq!(side.len(), 1);
    assert!(matches!(&side[0].role,Role::Side { segment,.. } if segment.as_str()=="right"));
    assert!(
        start
            .iter()
            .chain(&end)
            .chain(&side)
            .all(|o| o.step == id("capped"))
    );
}

#[test]
fn query_cannot_reuse_a_definition_to_escape_request_exhaustion() {
    let (_, field) = compile();
    let mut dec = Decider::default();
    let (result, work) = catallax_numbers::with_operation_budget(0, || {
        field.location(&RPoint3::new(q("0"), q("0"), q("2")), &mut dec)
    });
    assert!(work.exhausted);
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
}
