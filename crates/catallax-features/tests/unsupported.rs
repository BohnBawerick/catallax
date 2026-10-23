//! Legal geometry gaps remain distinct from impossible recipe inputs.

use catallax_features::{BuildError, evaluate};
use catallax_numbers::Decider;
use catallax_recipe::Recipe;
use serde_json::{Value, json};

fn refusal(loops: Value, operation: Value) -> String {
    let recipe = Recipe::parse(&json!({"format":catallax_recipe::FORMAT, "name":"legal gap", "steps":[
        {"id":"profile", "op":"sketch", "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}, "loops":loops}, operation
    ], "result":"body"}).to_string()).unwrap();
    match evaluate(&recipe, &mut Decider::default()).unwrap_err() {
        BuildError::Geometry { step, refusal } => {
            assert_eq!(step.as_str(), "body");
            refusal.code.to_string()
        }
        error => panic!("valid geometry must not fail recipe structure: {error}"),
    }
}
fn extrude(direction: [&str; 3]) -> Value {
    json!({"id":"body", "op":"extrude", "profile":"profile", "distance":"5", "direction":direction})
}
fn revolve(angle: &str) -> Value {
    json!({"id":"body", "op":"revolve", "profile":"profile", "axis":{"point":["0","0","0"],"direction":["0","1","0"]},"angle":angle})
}
fn rectangle() -> Value {
    json!([[{"kind":"line","id":"a","to":["4","0"]},{"kind":"line","id":"b","to":["4","2"]},{"kind":"line","id":"c","to":["0","2"]},{"kind":"line","id":"d","to":["0","0"]}]])
}

#[test]
fn unavailable_placements_and_elliptic_cylinders_have_specific_codes() {
    assert_eq!(
        refusal(rectangle(), extrude(["1", "0", "1"])),
        "features.unsupported-algebraic-placement"
    );
    assert_eq!(
        refusal(rectangle(), revolve("45")),
        "features.unsupported-algebraic-revolve-frame"
    );
    let circle = json!([[{"kind":"circle","id":"wall","centre":["0","0"],"radius":"1"}]]);
    assert_eq!(
        refusal(circle, extrude(["3", "0", "4"])),
        "features.unsupported-elliptic-cylinder"
    );
}

#[test]
fn irrational_generated_radii_name_the_missing_surface_family() {
    let spherical_band = json!([[{"kind":"arc","id":"round","centre":["0","0"],"to":["1","1"],"turn":"ccw"},{"kind":"line","id":"chord","to":["1","-1"]}]]);
    assert_eq!(
        refusal(spherical_band, revolve("360")),
        "features.unsupported-algebraic-sphere-radius"
    );
    let toroidal_band = json!([[{"kind":"arc","id":"round","centre":["3","0"],"to":["4","1"],"turn":"ccw"},{"kind":"line","id":"chord","to":["4","-1"]}]]);
    assert_eq!(
        refusal(toroidal_band, revolve("360")),
        "features.unsupported-algebraic-torus-radius"
    );
}

#[test]
fn valid_trimmed_spindle_torus_is_not_called_impossible() {
    let section = json!([[{"kind":"arc","id":"round","centre":["1","0"],"to":["1","2"],"turn":"ccw"},{"kind":"line","id":"chord","to":["1","-2"]}]]);
    assert_eq!(
        refusal(section, revolve("360")),
        "features.unsupported-spindle-torus"
    );
}
