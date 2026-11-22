//! Application measurement and numerical STEP of a multi-turn round-wire sweep.

use catallax_app::{NumericalOptions, numerical_sweeps};
use catallax_measure::{Accuracy, Enclosure};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_step::numerical as step;
use std::f64::consts::PI;

const APP_WORK: u64 = 100_000_000;
const TARGET_WORK: u64 = 20_000_000;

fn spring(turns: u32, pitch_mm: u32) -> Recipe {
    serde_json::from_value(serde_json::json!({
        "format": "catallax-recipe/5", "name": "spring", "steps": [
            {"id": "profile", "op": "sketch", "plane": {
                "rotation": ["1", "1", "0", "0"], "translation": ["0", "0", "0"]
            }, "loops": [[{"kind": "circle", "id": "wire", "centre": ["17", "2.5"], "radius": "1"}]]},
            {"id": "path", "op": "path", "start": ["17", "0", "2.5"], "segments": [
                {"kind": "helix", "id": "coil", "centre": ["0", "0", "2.5"],
                 "axis": ["0", "0", "1"], "angle": (360 * turns).to_string(),
                 "rise": (pitch_mm * turns).to_string()}
            ]},
            {"id": "spring", "op": "sweep", "profile": "profile", "path": "path",
             "orientation": "rotation_minimising"}
        ], "result": "spring"
    }))
    .unwrap()
}

fn near(enclosure: &Enclosure, expected: f64, difference: f64) {
    assert!(
        enclosure.lo.to_f64_nearest() >= expected - difference
            && enclosure.hi.to_f64_nearest() <= expected + difference,
        "{enclosure:?} should lie within {difference} of {expected}"
    );
}

#[test]
fn three_turn_spring_measures_and_exports() {
    measure_and_export(3, 3, TARGET_WORK);
}

fn measure_and_export(turns: u32, pitch_mm: u32, work_target: u64) {
    let recipe = spring(turns, pitch_mm);
    let (body, build_work) = with_operation_budget(APP_WORK, || {
        numerical_sweeps::build(
            &recipe,
            NumericalOptions::default(),
            &mut Decider::default(),
        )
    });
    let body = body.unwrap();
    assert!(!build_work.exhausted && build_work.used <= work_target);

    let (properties, property_work) = with_operation_budget(APP_WORK, || {
        body.properties(&Accuracy::default(), &mut Decider::default())
    });
    let properties = properties.unwrap();
    assert!(
        !property_work.exhausted && property_work.used <= work_target,
        "{property_work:?}"
    );
    assert!(properties.marks.is_empty(), "{:?}", properties.marks);

    let rise = f64::from(turns * pitch_mm);
    let length = f64::from(turns) * (2. * PI * 17.).hypot(f64::from(pitch_mm));
    let ideal_volume = PI * length;
    let ideal_area = 2. * PI * (length + 1.);
    near(&properties.fitted.volume, ideal_volume, ideal_volume / 100.);
    near(&properties.fitted.area, ideal_area, ideal_area / 100.);
    let centroid = properties.fitted.centroid.as_ref().unwrap();
    for (value, expected) in centroid.iter().zip([0., 0., 2.5 + rise / 2.]) {
        near(value, expected, 0.01);
    }

    let (file, export_work) = with_operation_budget(APP_WORK, || {
        let dec = &mut Decider::default();
        let file = step::build(body.body(), &Default::default(), dec)?;
        if let Some(path) = std::env::var_os("CATALLAX_SPRING_STEP") {
            std::fs::write(&path, file.bytes()).unwrap();
            file.check_written(&std::fs::read(path).unwrap(), dec)?;
        } else {
            file.check_written(file.bytes(), dec)?;
        }
        Ok::<_, catallax_numbers::Refusal>(file)
    });
    let file = file.unwrap();
    assert!(!export_work.exhausted && export_work.used <= work_target);
    let text = std::str::from_utf8(file.bytes()).unwrap();
    assert!(text.contains("catallax numerical approximate"));
    assert_eq!(text.matches("MANIFOLD_SOLID_BREP(").count(), 1);
    eprintln!("sweep properties {:?}", properties.fitted);
    eprintln!(
        "spring work: build {}, properties {}, STEP and reread {}",
        build_work.used, property_work.used, export_work.used
    );
}

#[test]
#[ignore = "near-allowance spring qualification"]
fn near_limit_spring_measures_and_excess_turn_refuses() {
    measure_and_export(20, 3, APP_WORK);

    let (result, work) = with_operation_budget(APP_WORK, || {
        numerical_sweeps::build(
            &spring(21, 3),
            NumericalOptions::default(),
            &mut Decider::default(),
        )
    });
    let Err(catallax_app::Error::Build(catallax_features::BuildError::Geometry {
        refusal, ..
    })) = result
    else {
        panic!("an over-allowance spring must return a geometry refusal");
    };
    assert_eq!(refusal.code, "numbers.operation-budget");
    assert!(work.exhausted);
}
