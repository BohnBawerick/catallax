//! Analytic path chains retain round-turn material and maximal physical names.

use catallax_features::{BuildError, Built, evaluate};
use catallax_numbers::{Ball, Decider, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_topology::BodyKind;
use serde_json::{Value, json};

fn source() -> Value {
    serde_json::from_str(include_str!("fixtures/round-turn-sweep.json")).unwrap()
}
fn build(source: Value) -> Built {
    let recipe = Recipe::parse(&source.to_string()).unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(
        built.model.check(built.result).is_empty(),
        "{:?}",
        built.model.check(built.result)
    );
    built
}
fn b(n: i64) -> Ball {
    Ball::from_rational(&Rational::from_i64(n), 192)
}
fn volume(built: &Built, expected: Ball) {
    let mass = catallax_measure::mass_properties(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(mass.marks.is_empty());
    assert!(
        mass.value.volume.lo <= expected.lo() && mass.value.volume.hi >= expected.hi(),
        "{:?} versus {expected:?}",
        mass.value.volume
    );
    assert!(mass.value.volume.width() <= Rational::parse("1e-9").unwrap());
}

#[test]
fn stationary_quarter_turn_adds_the_annular_sector_and_merges_planar_sides() {
    let built = build(source());
    // Two 2x1x4 prisms and a quarter annulus of radii 1 and 2, width 2.
    volume(&built, b(16) + b(3) * Ball::pi(192) / b(2));
    let faces = built.model.faces_of_body(built.result);
    let merged = faces
        .iter()
        .filter(|&&f| {
            let origins = built.model.face_origins(f);
            ["rise", "run"].iter().all(|path| {
                origins.iter().any(
                    |o| matches!(&o.role, Role::Side { path: Some(p), .. } if p.as_str() == *path),
                )
            })
        })
        .count();
    assert_eq!(
        merged, 2,
        "coplanar width sides must be maximal across the turn"
    );
    let origins: Vec<_> = faces
        .iter()
        .flat_map(|&f| built.model.face_origins(f))
        .collect();
    assert_eq!(origins.iter().filter(|o| o.role == Role::Start).count(), 1);
    assert_eq!(origins.iter().filter(|o| o.role == Role::End).count(), 1);
    assert!(
        !origins
            .iter()
            .any(|o| matches!(o.role, Role::Envelope { .. }))
    );
}

#[test]
fn tangent_line_arc_line_chain_preserves_authored_motion() {
    let mut source = source();
    source["steps"][1]["segments"] = json!([
        {"kind":"line","id":"rise","to":["0","0","4"]},
        {"kind":"arc","id":"bend","centre":["0","5","4"],"to":["0","5","9"],"axis":["-1","0","0"]},
        {"kind":"line","id":"run","to":["0","9","9"]}
    ]);
    // The profile lies at radii 6..7 about the arc centre, not about the path.
    let built = build(source);
    volume(&built, b(16) + b(13) * Ball::pi(192) / b(2));
}

#[test]
fn collinear_retrace_is_a_regularised_union_not_a_sum_of_piece_volumes() {
    let mut source = source();
    source["steps"][2]["orientation"] = json!("fixed");
    source["steps"][1]["segments"] = json!([
        {"kind":"line","id":"rise","to":["0","0","4"]},
        {"kind":"line","id":"run","to":["0","0","2"]}
    ]);
    let built = build(source);
    volume(&built, b(8));
    assert_eq!(built.model.faces_of_body(built.result).len(), 6);
    let origins: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|f| built.model.face_origins(f))
        .collect();
    assert!(!origins.iter().any(|o| o.role == Role::End));
    assert!(
        origins
            .iter()
            .any(|o| matches!(&o.role, Role::Envelope { path } if path.as_str() == "run"))
    );
}

#[test]
fn reversing_turn_uses_saved_axis_and_missing_axis_is_a_recipe_refusal() {
    let mut source = source();
    source["steps"][1]["segments"][1] = json!({"kind":"line","id":"run","to":["0","0","0"]});
    let recipe = Recipe::parse(&source.to_string()).unwrap();
    let BuildError::Geometry { step, refusal } =
        evaluate(&recipe, &mut Decider::default()).unwrap_err()
    else {
        panic!("expected missing-axis refusal")
    };
    assert_eq!(step.as_str(), "body");
    assert_eq!(refusal.code, "sweep.missing-corner-axis");
    source["steps"][1]["segments"][1]["corner_axis"] = json!(["-1", "0", "0"]);
    let built = build(source);
    // The half turn joins two disjoint prisms with a half annulus.
    volume(&built, b(16) + b(3) * Ball::pi(192));
}

#[test]
fn coplanar_translation_chain_has_no_three_dimensional_remnant() {
    let mut source = source();
    source["steps"][2]["orientation"] = json!("fixed");
    source["steps"][1]["segments"] = json!([
        {"kind":"line","id":"rise","to":["4","0","0"]},
        {"kind":"line","id":"run","to":["4","4","0"]}
    ]);
    let built = build(source);
    assert_eq!(built.model.body(built.result).kind, BodyKind::Empty);
}

#[test]
fn cardinal_arc_spans_follow_the_axis_and_sweep_angle() {
    for (end, quarters) in [
        (["0", "5", "5"], 1),
        (["0", "10", "0"], 2),
        (["0", "5", "-5"], 3),
    ] {
        let mut source = source();
        source["steps"][1]["segments"] = json!([{
            "kind":"arc","id":"bend","centre":["0","5","0"],
            "to":end,"axis":["-2","0","0"]
        }]);
        let built = build(source);
        volume(&built, b(13 * quarters) * Ball::pi(192) / b(2));
    }
}

#[test]
fn invalid_arc_and_corner_axes_fail_at_the_path_not_as_solid_geometry() {
    for segment in [
        json!({"kind":"arc","id":"bend","centre":["0","5","0"],"to":["0","0","0"],"axis":["-1","0","0"]}),
        json!({"kind":"arc","id":"bend","centre":["0","5","0"],"to":["0","5","6"],"axis":["-1","0","0"]}),
        json!({"kind":"line","id":"rise","to":["0","0","4"],"corner_axis":["1","0","0"]}),
    ] {
        let mut source = source();
        source["steps"][1]["segments"] = json!([segment]);
        let recipe = Recipe::parse(&source.to_string()).unwrap();
        let BuildError::Geometry { step, refusal } =
            evaluate(&recipe, &mut Decider::default()).unwrap_err()
        else {
            panic!("expected path failure")
        };
        assert_eq!(step.as_str(), "rail");
        assert!(matches!(refusal.code, "path.arc" | "path.corner-axis"));
    }
}
