//! SolidDraft wire conversion and application display for polynomial profiles.

use catallax_app::soliddraft::{ConversionError, convert};
use catallax_numbers::{Decider, Rational};
use serde_json::{Value, json};

fn source(controls: Value) -> Value {
    json!({"format":"soliddraft-recipe/1","name":"spline-profile","steps":[
        {"id":"profile","op":"profile","plane":"XY","loops":{"outer":[
            {"type":"bezier","from":[0,0],"controls":controls,"to":[4,0]},
            {"type":"line","from":[4,0],"to":[0,0]}
        ]}},
        {"id":"body","op":"extrude","profile":"profile","distance":2}
    ],"result":"body"})
}

#[test]
fn quadratic_and_cubic_profiles_reach_display_and_exact_volume() {
    for (controls, expected) in [(json!([[2, 4]]), "32/3"), (json!([[0, 4], [4, 4]]), "96/5")] {
        let app = catallax_app::build_soliddraft(
            &source(controls).to_string(),
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap()
        .body;
        assert!(!app.display.value.triangles.is_empty());
        assert!(app.built.model.check(app.built.result).is_empty());
        let mass = app
            .mass_properties(&Default::default(), &mut Decider::default())
            .unwrap();
        assert!(mass.marks.is_empty());
        assert!(
            mass.value
                .volume
                .contains(&Rational::parse(expected).unwrap())
        );
        assert!(mass.value.volume.width() < Rational::parse("1e-8").unwrap());
    }
}

#[test]
fn endpoint_canonicalization_keeps_the_documented_join_tolerance() {
    let mut accepted = source(json!([[0, 4], [4, 4]]));
    accepted["steps"][0]["loops"]["outer"][1]["from"] = json!([4.00000001, 0]);
    assert!(convert(&accepted.to_string()).is_ok());
    accepted["steps"][0]["loops"]["outer"][1]["from"] = json!([4.001, 0]);
    assert!(
        matches!(convert(&accepted.to_string()),Err(ConversionError::Invalid { detail,.. }) if detail=="open profile loop")
    );
}

#[test]
fn unsupported_degrees_and_mixed_curve_families_have_typed_conversion_errors() {
    for controls in [json!([]), json!([[0, 4], [2, 4], [4, 4]])] {
        assert!(
            matches!(convert(&source(controls).to_string()),Err(ConversionError::Invalid { detail,.. }) if detail=="a Bezier span needs one or two interior controls")
        );
    }
    let mut mixed = source(json!([[0, 4], [4, 4]]));
    mixed["steps"][0]["loops"]["outer"][1] =
        json!({"type":"arc","center":[2,0],"from":[4,0],"to":[0,0],"clockwise":true});
    assert!(matches!(
        convert(&mixed.to_string()),
        Err(ConversionError::NotYet {
            code: "soliddraft.bezier-mix",
            ..
        })
    ));
}
