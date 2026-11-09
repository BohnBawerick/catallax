#![allow(missing_docs)]

use catallax_features::{Built, evaluate};
use catallax_numbers::{Ball, Decider, Outcome, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

fn document() -> Value {
    serde_json::from_str(include_str!("fixtures/trimmed-planar-offset.json")).unwrap()
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
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn location(built: &mut Built, coordinates: [&str; 3], inside: bool) {
    let [x, y, z] = coordinates.map(q);
    let expected = if inside {
        catallax_boolean::Containment::Inside
    } else {
        catallax_boolean::Containment::Outside
    };
    let actual = catallax_boolean::classify_point(
        &mut built.model,
        built.result,
        &RPoint3::new(x, y, z),
        &mut Decider::default(),
    );
    assert!(
        matches!(actual, Outcome::Proven(value) if value == expected),
        "{coordinates:?}: {actual:?}"
    );
}
fn volume(built: &Built, constant: &str, pi: &str) {
    let actual = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    let expected =
        Ball::from_rational(&q(constant), 512) + Ball::pi(512) * Ball::from_rational(&q(pi), 512);
    assert!(
        actual.lo <= expected.lo() && actual.hi >= expected.hi(),
        "{actual:?}"
    );
    assert!(actual.width() < q("1e-6"));
}
#[test]
fn nonconvex_trimmed_erosion_uses_finite_edge_distance_instead_of_support_mitres() {
    let mut built = build(&document());
    for (point, inside) in [
        (["2.1", "2.1", "2"], true),
        (["2.5", "2.5", "2"], false),
        (["1.5", "5", "2"], true),
        (["0.5", "1.5", "2"], false),
        (["1.5", "1.5", "0.5"], false),
    ] {
        location(&mut built, point, inside);
    }
    let boundary = RPoint3::new(q("12/5"), q("11/5"), 2.into());
    assert!(matches!(
        catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &boundary,
            &mut Decider::default()
        ),
        Outcome::ProvenCoincident(_)
    ));
    volume(&built, "24", "-1/2");
}

#[test]
fn a_polygonal_hole_expands_with_round_corners_inside_the_shrunken_outer_loop() {
    let source: Value =
        serde_json::from_str(include_str!("fixtures/trimmed-planar-hole-offset.json")).unwrap();
    let mut built = build(&source);
    volume(&built, "64", "-2");
    for (point, inside) in [
        (["1.5", "5", "2"], true),
        (["2.5", "5", "2"], false),
        (["2.2", "2.2", "2"], true),
        (["2.4", "2.4", "2"], false),
    ] {
        location(&mut built, point, inside);
    }
}

#[test]
fn boolean_trimmed_faces_use_the_same_physical_envelope_as_an_authored_profile() {
    let source: Value =
        serde_json::from_str(include_str!("fixtures/trimmed-planar-boolean-offset.json")).unwrap();
    let mut built = build(&source);
    volume(&built, "24", "-1/2");
    location(&mut built, ["2.1", "2.1", "2"], true);
    location(&mut built, ["2.5", "2.5", "2"], false);
}

#[test]
fn consumed_cavities_regularise_to_empty_instead_of_self_intersecting() {
    for distance in ["-9/5", "-2", "-3"] {
        let mut source = document();
        source["steps"][2]["distance"] = json!(distance);
        let built = build(&source);
        assert_eq!(
            built.model.body(built.result).kind,
            catallax_topology::BodyKind::Empty
        );
    }
}

#[test]
fn rigid_placement_preserves_the_finite_trim_distance() {
    let mut source = document();
    source["steps"][0]["plane"] =
        json!({"rotation":["1","1","1","1"],"translation":["17","-4","9"]});
    let mut built = build(&source);
    volume(&built, "24", "-1/2");
    location(&mut built, ["19", "-1.9", "11.1"], true);
    location(&mut built, ["19", "-1.5", "11.5"], false);
}

#[test]
fn a_spatial_reentrant_vertex_uses_a_spherical_distance_envelope() {
    let source: Value =
        serde_json::from_str(include_str!("fixtures/trimmed-planar-corner-offset.json")).unwrap();
    let mut built = build(&source);
    location(&mut built, ["2.4", "2.4", "2.4"], true);
    location(&mut built, ["2.5", "2.5", "2.5"], false);
    let faces = built.model.faces_of_body(built.result);
    assert!(faces.iter().any(|f| matches!(
        built.model.geom().surface(built.model.face(*f).surface),
        catallax_shapes::Surface::Sphere(_)
    )));
    let error = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code, "measure.budget-exhausted");
    assert!(error.detail.contains("finite curved distance-feature trim"));
}
