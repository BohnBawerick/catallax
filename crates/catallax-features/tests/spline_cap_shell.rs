#![allow(missing_docs)]

use catallax_features::{Built, evaluate};
use catallax_numbers::{Decider, Outcome, Rational, with_operation_budget};
use catallax_recipe::{Origin, Recipe, Role, StepId};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

fn document() -> Value {
    serde_json::from_str(include_str!("fixtures/spline-cap-shell.json")).unwrap()
}

fn build(source: &Value) -> Built {
    let recipe = Recipe::parse(&source.to_string()).unwrap();
    let (answer, work) =
        with_operation_budget(100_000_000, || evaluate(&recipe, &mut Decider::default()));
    let built = answer.unwrap();
    assert!(!work.exhausted);
    assert!(built.model.check(built.result).is_empty());
    built
}

fn volume(built: &Built, expected: &str) {
    let measured = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    let expected = Rational::parse(expected).unwrap();
    assert!(
        measured.lo <= expected && measured.hi >= expected,
        "{measured:?}"
    );
    assert!(measured.width() < Rational::parse("1e-9").unwrap());
}

#[test]
fn one_kept_spline_prism_cap_has_exact_axial_material() {
    let source = document();
    let mut built = build(&source);
    volume(&built, "48/5");
    assert!(
        catallax_topology::BezierExtrusionBand::in_body(
            &built.model,
            built.result,
            &mut Decider::default()
        )
        .unwrap()
    );
    let inner = Origin {
        step: StepId::new("shell").unwrap(),
        role: Role::Inner {
            of: Box::new(Origin {
                step: StepId::new("stock").unwrap(),
                role: Role::Start,
            }),
        },
    };
    assert!(
        built
            .model
            .faces_of_body(built.result)
            .iter()
            .any(|f| built.model.face_origins(*f).contains(&inner))
    );
    for (z, expected) in [("1/2", true), ("3/2", false), ("7/2", false)] {
        let point = RPoint3::new(2.into(), 1.into(), Rational::parse(z).unwrap());
        let answer = catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &point,
            &mut Decider::default(),
        );
        let expected = if expected {
            catallax_boolean::Containment::Inside
        } else {
            catallax_boolean::Containment::Outside
        };
        assert!(
            matches!(answer, Outcome::Proven(actual) if actual == expected),
            "{z}: {answer:?}"
        );
    }
    let boundary = RPoint3::new(2.into(), 1.into(), 1.into());
    assert!(matches!(
        catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &boundary,
            &mut Decider::default()
        ),
        Outcome::ProvenCoincident(_)
    ));
}

#[test]
fn both_caps_produce_two_slabs_and_coincident_cavity_limits_keep_the_source() {
    let mut source: Value =
        serde_json::from_str(include_str!("fixtures/spline-two-cap-shell.json")).unwrap();
    let built = build(&source);
    volume(&built, "96/5");
    assert_eq!(built.model.body(built.result).lumps.len(), 2);
    for thickness in ["2", "3", "4", "5"] {
        source["steps"][2]["thickness"] = json!(thickness);
        let built = build(&source);
        volume(&built, "192/5");
    }
}

#[test]
fn quadratic_reversed_and_rigidly_placed_profiles_keep_exact_volume() {
    for (controls, expected) in [
        (json!([["2", "4"]]), "16/3"),
        (json!([["0", "4"], ["4", "4"]]), "48/5"),
    ] {
        for reverse in [false, true] {
            let mut source = document();
            source["steps"][0]["loops"][0][0]["controls"] = controls.clone();
            if reverse {
                source["steps"][1]["distance"] = json!("-4");
            }
            source["steps"][0]["plane"] =
                json!({"rotation":["2","1","3","4"],"translation":["17","-4","9"]});
            volume(&build(&source), expected);
        }
    }
}

#[test]
fn a_spline_through_cut_keeps_the_hole_in_the_exact_cap_domain() {
    let source: Value =
        serde_json::from_str(include_str!("fixtures/spline-cut-cap-shell.json")).unwrap();
    let mut built = build(&source);
    volume(&built, "272/5");
    for (x, inside) in [("-1", true), ("2", false)] {
        let point = RPoint3::new(
            Rational::parse(x).unwrap(),
            1.into(),
            Rational::parse("1/2").unwrap(),
        );
        let actual = catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &point,
            &mut Decider::default(),
        );
        let expected = if inside {
            catallax_boolean::Containment::Inside
        } else {
            catallax_boolean::Containment::Outside
        };
        assert!(
            matches!(actual, Outcome::Proven(value) if value == expected),
            "{actual:?}"
        );
    }
}

#[test]
fn unsupported_normal_envelopes_and_invalid_selection_keep_named_refusals() {
    for (open, thickness, expected) in [
        (json!([]), "1", "distance.unsupported-envelope"),
        (
            json!([{"face":{"step":"stock","role":"start"},"touching":[{"step":"stock","role":"end"}]}]),
            "1",
            "reference.lost-face",
        ),
        (
            document()["steps"][2]["open"].clone(),
            "0",
            "distance.non-positive-thickness",
        ),
    ] {
        let mut source = document();
        source["steps"][2]["open"] = open;
        source["steps"][2]["thickness"] = json!(thickness);
        let recipe = Recipe::parse(&source.to_string()).unwrap();
        let error = evaluate(&recipe, &mut Decider::default()).unwrap_err();
        let catallax_features::BuildError::Geometry { refusal, .. } = error else {
            panic!("{error:?}")
        };
        assert_eq!(refusal.code, expected);
    }
}

#[test]
fn a_spline_pocket_is_not_mistaken_for_a_complete_prism() {
    let mut source: Value =
        serde_json::from_str(include_str!("fixtures/spline-cut-cap-shell.json")).unwrap();
    source["steps"][2]["distance"] = json!("3");
    let recipe = Recipe::parse(&source.to_string()).unwrap();
    let error = evaluate(&recipe, &mut Decider::default()).unwrap_err();
    let catallax_features::BuildError::Geometry { step, refusal } = error else {
        panic!("{error:?}")
    };
    assert_eq!(step, StepId::new("shell").unwrap());
    assert_eq!(refusal.code, "distance.unsupported-envelope");
}
