//! Construction remains usable when a later measurement exhausts its budget.
use catallax_cli::{build_question, parse_question};
use catallax_measure::{Accuracy, mass_properties};
use catallax_numbers::{Ball, Decider, Rational, with_operation_budget};
use serde_json::json;

#[test]
fn measurement_is_separate_and_cannot_discard_a_built_body() {
    let question = parse_question(
        &json!({
            "format":"exact-brep-question/2", "id":"deferred-union",
            "operands": {
                "block":{"shape":"block","size":["2","2","2"]},
                "sphere":{"shape":"sphere","radius":"1"}
            },
            "steps":[{"id":"union","op":"unite","args":["block","sphere"]}],
            "result":"union", "limits":{"deadline_ms":0}
        })
        .to_string(),
    )
    .unwrap();
    let built = build_question(&question, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.body).is_empty());
    let counts = built.model.counts(built.body);
    let (failed, work) = with_operation_budget(0, || {
        mass_properties(
            &built.model,
            built.body,
            &Accuracy::default(),
            &mut Decider::default(),
        )
    });
    assert!(failed.is_err());
    assert!(work.exhausted);
    let measured = mass_properties(
        &built.model,
        built.body,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap();
    // The sphere adds seven octants to the cube. Its eighth octant is inside.
    let pi = Ball::pi(256);
    let r = |s| Rational::parse(s).unwrap();
    let expected_volume =
        &Ball::from_rational(&r("8"), 256) + &(&pi * &Ball::from_rational(&r("7/6"), 256));
    let expected_area =
        &Ball::from_rational(&r("24"), 256) + &(&pi * &Ball::from_rational(&r("11/4"), 256));
    assert!(measured.value.volume.lo <= expected_volume.lo());
    assert!(measured.value.volume.hi >= expected_volume.hi());
    assert!(measured.value.area.lo <= expected_area.lo());
    assert!(measured.value.area.hi >= expected_area.hi());
    // The three cube first moments are 8 each; subtract each sphere-octant
    // first moment pi/16 because the complete sphere's first moment is zero.
    let moment =
        &Ball::from_rational(&r("8"), 256) - &(&pi * &Ball::from_rational(&r("1/16"), 256));
    let center = moment.checked_div(&expected_volume).unwrap();
    for enclosure in measured.value.centroid.unwrap() {
        assert!(enclosure.lo <= center.lo());
        assert!(enclosure.hi >= center.hi());
    }
    assert_eq!(built.model.counts(built.body), counts);
    assert!(built.model.check(built.body).is_empty());
    assert!(built.marks.is_empty());
}
