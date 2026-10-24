//! Resource regressions run outside the test runner, including stack-aborting bugs.
#![cfg(target_os = "linux")]

use std::process::Command;

use catallax_app::{Error, build_soliddraft};
use catallax_mesh::MeshOptions;
use catallax_numbers::Decider;
use serde_json::{Value, json};

const CASE_ENV: &str = "CATALLAX_RESOURCE_CASE";

fn run_case(case: &str) {
    let output = Command::new("timeout")
        .args([
            "--kill-after=1",
            "15s",
            "prlimit",
            "--as=536870912",
            "--stack=8388608",
            "--cpu=12",
            "--core=0",
            "--",
        ])
        .arg(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "resource_child",
            "--nocapture",
            "--test-threads=1",
        ])
        .env(CASE_ENV, case)
        .output()
        .expect("Linux resource regressions require timeout and prlimit");
    assert!(
        output.status.success(),
        "{case}: {}\nstdout: {}\nstderr: {}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

#[test]
fn tiny_plane_coordinates_finish_under_process_caps() {
    run_case("ordinary");
    run_case("boss:5e-324");
    run_case("hole:1e-20");
    run_case("hole2:1e-20");
    run_case("hole:0");
    for exponent in 3..=12 {
        run_case(&format!("hole:1e-{exponent}"));
    }
}

#[test]
fn huge_circle_finishes_under_process_caps() {
    run_case("huge");
}

#[test]
fn resource_child() {
    let Ok(case) = std::env::var(CASE_ENV) else {
        return;
    };
    let mut input: Value = serde_json::from_str(include_str!("fixtures/bracket.json")).unwrap();
    if let Some((kind, number)) = case.split_once(':') {
        let step = match kind {
            "boss" => 5,
            "hole2" => 14,
            _ => 11,
        };
        input["steps"][step]["plane"]["origin"][0] = number.parse::<f64>().unwrap().into();
    } else if case == "huge" {
        input = json!({"format":"soliddraft-recipe/1","name":"large-circle","steps":[
            {"id":"circle","op":"profile","plane":"XY","circle":{"center":[0,0],"radius":500000}},
            {"id":"body","op":"extrude","profile":"circle","distance":10}
        ],"result":"body"});
    }
    let result = build_soliddraft(
        &input.to_string(),
        &MeshOptions::default(),
        &mut Decider::default(),
    );
    match result {
        Ok(body) => {
            assert_ne!(
                case, "huge",
                "huge circle must respect the display resource allowance"
            );
            let mesh = &body.body.display.value;
            assert!(
                mesh.faces
                    .iter()
                    .all(|face| face.max_chord_error <= MeshOptions::default().chord_tolerance)
            );
            println!(
                "answered {case}: {} vertices, {} triangles",
                mesh.positions.len(),
                mesh.triangles.len()
            );
        }
        Err(error @ Error::Build(_)) if case != "ordinary" && case != "huge" => {
            assert_eq!(error.code(), "boolean.sampling-budget", "{error:?}");
            println!("refused {case}: {error:?}");
        }
        Err(Error::Display(error)) if case == "huge" => {
            assert_eq!(error.code, "mesh.resource-limit", "{error:?}");
            println!("refused {case}: {error:?}");
        }
        other => panic!("unexpected outcome for {case}: {other:?}"),
    }
}
