//! Application conversion, exact bodies, display, measures and cuts for ellipses.

use catallax_numbers::{Ball, Decider, Rational};
use serde_json::json;

fn source(segments: serde_json::Value) -> serde_json::Value {
    json!({"format":"soliddraft-recipe/1", "name":"ellipse", "steps":[
        {"id":"profile","op":"profile","plane":"XY","loops":{"outer":segments}},
        {"id":"body","op":"extrude","profile":"profile","distance":2.0}
    ], "result":"body"})
}

#[test]
fn ellipse_and_partial_ellipse_reach_application_display_and_mass() {
    for (segments, coefficient) in [
        (
            json!([{"type":"ellipse","center":[1,2],"axis":[4,0],"ratio":0.5}]),
            16,
        ),
        (
            json!([
                {"type":"elliptic_arc","center":[0,0],"axis":[4,0],"ratio":0.5,"from":[4,0],"to":[-4,0]},
                {"type":"line","from":[-4,0],"to":[4,0]}
            ]),
            8,
        ),
    ] {
        eprintln!("ellipse {coefficient}: building application display");
        let app = catallax_app::build_soliddraft(
            &source(segments).to_string(),
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap()
        .body;
        assert!(!app.display.value.triangles.is_empty());
        assert!(app.built.model.check(app.built.result).is_empty());
        eprintln!("ellipse {coefficient}: display complete, measuring");
        let mass = app
            .mass_properties(&Default::default(), &mut Decider::default())
            .unwrap();
        let expected = Ball::from_rational(&Rational::from_i64(coefficient), 192) * Ball::pi(192);
        assert!(mass.value.volume.lo <= expected.lo() && mass.value.volume.hi >= expected.hi());
        assert!(mass.value.volume.width() < Rational::parse("1e-8").unwrap());
    }
}

#[test]
fn ellipse_tool_cuts_a_block_through_the_application_entry() {
    let mut source = source(json!([{"type":"ellipse","center":[0,0],"axis":[4,0],"ratio":0.5}]));
    let steps = source["steps"].as_array_mut().unwrap();
    steps.insert(0, json!({"id":"base_profile","op":"profile","plane":"XY","rectangle":{"center":[0,0],"width":12,"height":8}}));
    steps.insert(
        1,
        json!({"id":"base","op":"extrude","profile":"base_profile","distance":2}),
    );
    steps.push(json!({"id":"cut","op":"boolean","kind":"subtract","target":"base","tool":"body"}));
    source["result"] = json!("cut");
    eprintln!("ellipse cut: building application display");
    let app = catallax_app::build_soliddraft(
        &source.to_string(),
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    eprintln!("ellipse cut: display complete, measuring");
    let mass = app
        .mass_properties(&Default::default(), &mut Decider::default())
        .unwrap();
    let b = |n| Ball::from_rational(&Rational::from_i64(n), 192);
    let expected = b(192) - b(16) * Ball::pi(192);
    assert!(mass.value.volume.lo <= expected.lo() && mass.value.volume.hi >= expected.hi());
}
