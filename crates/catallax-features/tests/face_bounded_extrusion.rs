//! Exact face starts and ends through the public recipe evaluator.
use catallax_features::{Built, evaluate};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::{Recipe, Role};
use serde_json::{Value, json};

fn pose(translation: [&str; 3]) -> Value {
    json!({"rotation":["1","0","0","0"],"translation":translation})
}

fn rectangle(x: [&str; 2], y: [&str; 2]) -> Value {
    json!([[
        {"kind":"line","id":"bottom","to":[x[1],y[0]]},
        {"kind":"line","id":"right","to":[x[1],y[1]]},
        {"kind":"line","id":"top","to":[x[0],y[1]]},
        {"kind":"line","id":"left","to":[x[0],y[0]]}
    ]])
}

fn run(value: &Value) -> Built {
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let (built, work) =
        with_operation_budget(100_000_000, || evaluate(&recipe, &mut Decider::default()));
    let built = built.unwrap_or_else(|error| panic!("{}: {error} ({work:?})", recipe.name));
    let (faults, validation) = with_operation_budget(100_000_000, || {
        catallax_boolean::check_body(&built.model, built.result, &mut Decider::default())
    });
    let faults = faults.unwrap_or_else(|error| panic!("{}: {error} ({validation:?})", recipe.name));
    assert!(faults.is_empty(), "{}: {faults:?}", recipe.name);
    assert!(
        built.model.body(built.result).marks.is_empty(),
        "{}",
        recipe.name
    );
    eprintln!("{}: {} build work", recipe.name, work.used);
    built
}

fn volume(built: &Built, expected: f64) {
    let (actual, work) = with_operation_budget(100_000_000, || {
        catallax_measure::volume(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
    });
    let actual = actual.unwrap_or_else(|error| panic!("{error} ({work:?})"));
    let low = actual.value.lo.to_f64_nearest();
    let high = actual.value.hi.to_f64_nearest();
    assert!(
        (low - expected).abs() < 1e-6 && (high - expected).abs() < 1e-6,
        "{low}..{high} vs {expected}"
    );
}

fn cap_names(built: &Built) {
    let names: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|face| built.model.face_origins(face))
        .collect();
    assert!(names.iter().all(|origin| origin.step.as_str() == "cap"));
    assert!(names.iter().any(|origin| origin.role == Role::Start));
    assert!(names.iter().any(|origin| origin.role == Role::End));
    assert!(
        names
            .iter()
            .any(|origin| matches!(origin.role, Role::Side { .. }))
    );
}

fn wall(segment: &str) -> Value {
    json!({"face":{"step":"target","role":{"side":{"segment":segment}}}})
}

fn circle_sweep() -> Value {
    json!({"format":"catallax-recipe/5","name":"circle sweep target","steps":[
        {"id":"section","op":"sketch",
         "plane":{"rotation":["1","1","1","1"],"translation":["-5","0","8"]},
         "loops":[[{"kind":"circle","id":"wall","centre":["0","0"],"radius":"5"}]]},
        {"id":"target","op":"extrude","profile":"section","distance":"10"},
        {"id":"profile","op":"sketch","plane":pose(["0","0","0"]),
         "loops":rectangle(["-1","1"],["-3","3"])},
        {"id":"cap","op":"extrude_to","target":"target","profile":"profile","limit":wall("wall")}
    ],"result":"cap"})
}

#[test]
fn the_wall_of_a_circle_sweep_is_a_trimmed_end() {
    let built = run(&circle_sweep());
    cap_names(&built);
    volume(&built, 72.0 - 50.0 * 0.6_f64.asin());
}

#[test]
fn the_planar_cap_of_a_circle_sweep_is_a_trimmed_end() {
    let mut value = circle_sweep();
    value["steps"][2]["plane"] = json!({"rotation":["1","1","1","1"],"translation":["-8","0","8"]});
    value["steps"][2]["loops"] = rectangle(["-1", "1"], ["-1", "1"]);
    value["steps"][3]["limit"] = json!({"face":{"step":"target","role":"end"}});
    let built = run(&value);
    cap_names(&built);
    volume(&built, 52.0);
}

#[test]
fn blind_from_an_annular_torus_has_two_exact_curved_caps() {
    let mut value: Value =
        serde_json::from_str(include_str!("fixtures/curved-limits/torus-entry-next.json")).unwrap();
    value["format"] = json!("catallax-recipe/5");
    value["steps"][2] = json!({
        "id":"cap","op":"extrude_from","target":"target","profile":"profile",
        "start":wall("wall"),"distance":"2"
    });
    let built = run(&value);
    cap_names(&built);
    volume(&built, 12.0 * std::f64::consts::PI);
}

fn from(mut value: Value, distance: &str) -> Value {
    value["format"] = json!("catallax-recipe/5");
    let steps = value["steps"].as_array_mut().unwrap();
    let cap = steps.last_mut().unwrap().as_object_mut().unwrap();
    let limit = cap.remove("limit").unwrap();
    cap.insert("op".into(), json!("extrude_from"));
    cap.insert("start".into(), limit);
    cap.insert("distance".into(), json!(distance));
    value
}

#[test]
fn blind_start_copies_quadric_caps_with_the_authored_signed_depth() {
    for (label, input, area) in [
        ("cylinder", circle_sweep(), 12.0),
        (
            "sphere",
            serde_json::from_str(include_str!(
                "fixtures/circular-sphere-limits/sphere-entry-to.json"
            ))
            .unwrap(),
            9.0 * std::f64::consts::PI,
        ),
        (
            "cone",
            serde_json::from_str(include_str!("fixtures/curved-limits/cone-entry-to.json"))
                .unwrap(),
            1.0,
        ),
    ] {
        for depth in ["2", "-2"] {
            let mut value = from(input.clone(), depth);
            value["name"] = json!(format!("{label} start {depth}"));
            let built = run(&value);
            cap_names(&built);
            volume(&built, 2.0 * area);
        }
    }
}

fn spline_roof() -> Value {
    json!({"format":"catallax-recipe/5","name":"Bezier roof limit","steps":[
        {"id":"section","op":"sketch",
         "plane":{"rotation":["1","1","0","0"],"translation":["0","0","0"]},
         "loops":[[
            {"kind":"bezier","id":"roof","controls":[["2","4"]],"to":["4","0"]},
            {"kind":"line","id":"floor","to":["0","0"]}
         ]]},
        {"id":"target","op":"extrude","profile":"section","direction":["0","1","0"],"distance":"4"},
        {"id":"profile","op":"sketch","plane":pose(["0","0","3"]),
         "loops":rectangle(["1","3"],["1","3"])},
        {"id":"cap","op":"extrude_to","target":"target","profile":"profile","direction":["0","0","-1"],"limit":wall("roof")}
    ],"result":"cap"})
}

#[test]
fn a_bezier_extrusion_wall_can_be_both_cap_types() {
    let to = spline_roof();
    let built = run(&to);
    cap_names(&built);
    volume(&built, 14.0 / 3.0);
    let built = run(&from(to, "2"));
    cap_names(&built);
    volume(&built, 8.0);
}

#[test]
fn stock_sphere_cone_and_torus_ends_match_geometric_integrals() {
    let mut torus: Value =
        serde_json::from_str(include_str!("fixtures/curved-limits/torus-entry-next.json")).unwrap();
    torus["steps"][2]["op"] = json!("extrude_to");
    torus["steps"][2]["limit"] = wall("wall");
    let pi = std::f64::consts::PI;
    for (value, expected) in [
        (
            serde_json::from_str(include_str!(
                "fixtures/circular-sphere-limits/sphere-entry-to.json"
            ))
            .unwrap(),
            pi * (90.0 - (2.0 / 3.0) * (125.0 - 64.0)),
        ),
        // Integrate 5 - sqrt((4 - z/5)^2 - y^2) over |y| <= 1/2 and |z - 5| <= 1/2.
        (
            serde_json::from_str(include_str!("fixtures/curved-limits/cone-entry-to.json"))
                .unwrap(),
            2.013_952_620_81,
        ),
        (
            torus,
            30.0 * pi - 6.0 * pi * (0.5 * 0.75_f64.sqrt() + 0.5_f64.asin()),
        ),
    ] {
        let built = run(&value);
        cap_names(&built);
        volume(&built, expected);
    }
}

#[test]
fn a_curved_start_can_merge_with_and_cut_another_solid() {
    for (kind, expected) in [("union", 40.0), ("subtract", 20.0)] {
        let mut value = from(circle_sweep(), "2");
        value["name"] = json!(format!("curved start {kind}"));
        let steps = value["steps"].as_array_mut().unwrap();
        steps.push(
            json!({"id":"block","op":"block","size":["1","2","10"],"pose":pose(["-1/2","-1","0"])}),
        );
        steps.push(json!({"id":"result","op":"boolean","kind":kind,"target":"cap","tool":"block"}));
        value["result"] = json!("result");
        let built = run(&value);
        volume(&built, expected);
    }
}

#[test]
fn a_round_profile_keeps_its_rim_on_a_transverse_cylinder() {
    let mut value = circle_sweep();
    value["steps"][2]["loops"] = json!([[
        {"kind":"circle","id":"rim","centre":["0","0"],"radius":"2"}
    ]]);
    let built = run(&value);
    cap_names(&built);
    // Integrate (8 - sqrt(25 - y²)) over the disk of radius 2.
    volume(&built, 38.982_234_412_964);
    let built = run(&from(value, "2"));
    cap_names(&built);
    volume(&built, 8.0 * std::f64::consts::PI);
}
