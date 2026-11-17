//! Finite fillet regions may be clear when their infinite carriers intersect.

use catallax_boolean::Containment;
use catallax_features::{Built, evaluate};
use catallax_numbers::{Decider, Outcome, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
use serde_json::{Value, json};
use std::f64::consts::PI;

fn pose(x: &str, y: &str, z: &str) -> Value {
    json!({"rotation":["1","0","0","0"],"translation":[x,y,z]})
}

fn build(name: &str, steps: Vec<Value>) -> Built {
    let recipe = Recipe::parse(
        &json!({
            "format":"catallax-recipe/3", "name":name, "steps":steps, "result":"round",
        })
        .to_string(),
    )
    .unwrap();
    let (answer, work) =
        with_operation_budget(100_000_000, || evaluate(&recipe, &mut Decider::default()));
    println!("{name} build {work:?}");
    assert!(!work.exhausted, "{name}: {work:?}");
    let built = answer.unwrap();
    assert_eq!(built.model.check(built.result), vec![]);
    if let Some(directory) = std::env::var_os("CATALLAX_TRIMMED_FILLET_STEP_DIR") {
        let mut bytes = Vec::new();
        let (answer, work) = with_operation_budget(100_000_000, || {
            catallax_step::write_step(
                &built.model,
                &[built.result],
                &Default::default(),
                &mut Decider::default(),
                &mut bytes,
            )
        });
        println!("{name} STEP {work:?}");
        assert!(!work.exhausted, "{name}: {work:?}");
        assert!(answer.unwrap().marks.is_empty());
        std::fs::create_dir_all(&directory).unwrap();
        std::fs::write(
            std::path::Path::new(&directory).join(format!("{name}.step")),
            bytes,
        )
        .unwrap();
    }
    built
}

fn check_volume(built: &Built, expected: f64) {
    let (answer, work) = with_operation_budget(100_000_000, || {
        catallax_measure::volume_area_centroid(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
    });
    println!("independent volume {work:?}");
    assert!(!work.exhausted, "{work:?}");
    let measured = answer.unwrap().value.volume;
    let (lo, hi) = (measured.lo.to_f64_nearest(), measured.hi.to_f64_nearest());
    assert!(
        lo - 1e-8 <= expected && expected <= hi + 1e-8 && hi - lo < 1e-6,
        "volume [{lo}, {hi}] excludes {expected}"
    );
}

fn ring_volume(wall: f64, radius: f64) -> f64 {
    2.0 * PI
        * (wall * radius * radius * (1.0 - PI / 4.0) + radius.powi(3) * (10.0 - 3.0 * PI) / 12.0)
}

fn ring_selection(body: &str, wall: &str, radius: &str) -> Value {
    json!({"id":"round", "op":"fillet", "body":body, "radius":radius,
        "edges":[{"between":[{"step":"block","role":"end"},
            {"step":wall,"role":{"side":{"segment":"wall"}}}]}]})
}

#[test]
fn trimmed_corner_sector_stays_above_the_opposite_face_at_large_radii() {
    for radius in ["10", "15", "19"] {
        let mut input: Value =
            serde_json::from_str(include_str!("fixtures/trihedral-corner.json")).unwrap();
        input["steps"][0]["size"] = json!(["60", "40", "20"]);
        input["steps"][1]["radius"] = json!(radius);
        let built = build(
            &format!("corner-r{radius}"),
            input["steps"].as_array().unwrap().clone(),
        );
        let radius: f64 = radius.parse().unwrap();
        let strips = (120.0 - 3.0 * radius) * radius * radius;
        let corner = radius.powi(3);
        check_volume(
            &built,
            48_000.0 - strips - corner + PI * (strips / 4.0 + corner / 6.0),
        );
    }
}

#[test]
fn trimmed_boss_base_rings_stop_at_the_contact_radius_not_the_full_torus() {
    for radius in ["4", "5", "6"] {
        let built = build(
            &format!("boss-base-r{radius}"),
            vec![
                json!({"id":"block","op":"block","size":["40","40","20"],"pose":pose("-20","-20","0")}),
                json!({"id":"boss","op":"cylinder","radius":"12","height":"25","pose":pose("0","0","10")}),
                json!({"id":"joined","op":"boolean","kind":"union","target":"block","tool":"boss"}),
                ring_selection("joined", "boss", radius),
            ],
        );
        let radius: f64 = radius.parse().unwrap();
        check_volume(
            &built,
            32_000.0 + PI * 144.0 * 15.0 + ring_volume(12.0, radius),
        );
    }
}

#[test]
fn trimmed_hole_rim_stays_clear_of_a_wall_twelve_mm_from_its_axis() {
    let built = build(
        "near-wall-hole-r3",
        vec![
            json!({"id":"block","op":"block","size":["40","40","20"],"pose":pose("-28","-20","0")}),
            json!({"id":"hole","op":"cylinder","radius":"6","height":"30","pose":pose("0","0","-5")}),
            json!({"id":"cut","op":"boolean","kind":"subtract","target":"block","tool":"hole"}),
            ring_selection("cut", "hole", "3"),
        ],
    );
    check_volume(&built, 32_000.0 - PI * 36.0 * 20.0 - ring_volume(6.0, 3.0));
}

#[test]
fn trimmed_cone_ring_does_not_touch_the_cylinder_below_its_chamfer() {
    let mut built = build(
        "chamfered-r2",
        vec![
            json!({"id":"c","op":"cylinder","radius":"5","height":"10","pose":pose("0","0","0")}),
            json!({"id":"bevel","op":"chamfer","body":"c","distance":"1", "edges":[{"between":[
            {"step":"c","role":"end"}, {"step":"c","role":{"side":{"segment":"wall"}}}]}]}),
            json!({"id":"round","op":"fillet","body":"bevel","radius":"2", "edges":[{"between":[
            {"step":"c","role":"end"}, {"step":"bevel"}]}]}),
        ],
    );
    for (x, z, expected) in [
        ("49/10", "17/2", Containment::Inside),
        ("4", "19/2", Containment::Inside),
        ("4", "999/100", Containment::Outside),
    ] {
        let point = RPoint3::new(
            Rational::parse(x).unwrap(),
            Rational::zero(),
            Rational::parse(z).unwrap(),
        );
        assert!(
            matches!(catallax_boolean::classify_point(
            &mut built.model, built.result, &point, &mut Decider::default(),
        ), Outcome::Proven(actual) if actual == expected),
            "{point:?}"
        );
    }
}

#[test]
fn trimmed_four_contact_cone_stays_clear_of_the_base_when_its_ball_does_not() {
    let row: Value = include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|row| row["id"] == "vertex-4-all-unequal-slopes")
        .unwrap();
    let mut steps = row["recipe"]["steps"].as_array().unwrap().clone();
    for step in &mut steps {
        if step["id"] == "bs" {
            step["loops"][0][0]["to"][0] = json!("-10");
        }
        if step["op"] == "fillet" {
            step["radius"] = json!("5");
        }
    }
    let mut built = build("four-contact-r5", steps);
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let plan = catallax_field::Plan::compile_checked(built.recipe(), &mut dec).unwrap();
        let mut storage = catallax_field::Workspace::default();
        // The common centre is at z=10-5*sqrt(2). Its retained cap ends at
        // z=15-5*sqrt(2), between 7.9 and 8, while the complete ball crosses z=0.
        for (coordinates, inside) in [
            (["0", "0", "79/10"], true),
            (["0", "0", "8"], false),
            (["2", "0", "7"], true),
            (["2", "0", "39/5"], false),
        ] {
            let [x, y, z] = coordinates.map(|value| Rational::parse(value).unwrap());
            let point = RPoint3::new(x, y, z);
            let expected = if inside {
                Containment::Inside
            } else {
                Containment::Outside
            };
            let answer = catallax_boolean::classify_point(
                &mut built.model, built.result, &point, &mut dec,
            );
            println!("four-contact B-rep {point:?}: {answer:?}; {:?}", catallax_numbers::current_operation_work());
            assert!(
                matches!(&answer, Outcome::Proven(actual) if *actual == expected),
                "B-rep at {point:?}: {answer:?}"
            );
            let expected = if inside {
                catallax_field::Location::Inside
            } else {
                catallax_field::Location::Outside
            };
            assert_eq!(
                plan.proven_checked(&point, &mut storage, &mut dec)
                    .unwrap()
                    .location,
                expected,
                "field at {point:?}"
            );
        }
    });
    assert!(!work.exhausted, "{work:?}");
    println!("four-contact field and material {work:?}");
}
