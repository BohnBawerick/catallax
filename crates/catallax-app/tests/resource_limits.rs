//! Resource regressions run outside the test runner, including stack-aborting bugs.
#![cfg(target_os = "linux")]

use std::{
    collections::{BTreeMap, BTreeSet},
    process::Command,
};

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
fn fine_chord_drawings_finish_under_process_caps() {
    run_case("annulus-drawing");
    run_case("bracket-drawing");
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
    let mut options = MeshOptions::default();
    if case == "annulus-drawing" {
        input = json!({"format":"soliddraft-recipe/1","name":"annulus","steps":[
            {"id":"outer_profile","op":"profile","plane":"XY","circle":{"center":[0,0],"radius":10}},
            {"id":"outer","op":"extrude","profile":"outer_profile","distance":20},
            {"id":"inner_profile","op":"profile","plane":{"base":"XY","origin":[0,0,-1]},"circle":{"center":[0,0],"radius":4}},
            {"id":"inner","op":"extrude","profile":"inner_profile","distance":22},
            {"id":"annulus","op":"boolean","kind":"subtract","target":"outer","tool":"inner"}
        ],"result":"annulus"});
        options.chord_tolerance = "3467565716752893/10000000000000000000".parse().unwrap();
    } else if case == "bracket-drawing" {
        options.chord_tolerance = "16158836963098553/10000000000000000000".parse().unwrap();
    }
    let result = build_soliddraft(&input.to_string(), &options, &mut Decider::default());
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
                    .all(|face| face.max_chord_error <= options.chord_tolerance)
            );
            assert!(
                mesh.edges
                    .iter()
                    .all(|edge| edge.max_chord_error <= options.chord_tolerance)
            );
            if case.ends_with("-drawing") {
                let expected: BTreeSet<_> = body
                    .body
                    .built
                    .model
                    .faces_of_body(body.body.built.result)
                    .into_iter()
                    .collect();
                let actual: BTreeSet<_> = mesh.triangles.iter().map(|t| t.face).collect();
                assert_eq!(actual, expected);
                let mut incidence = BTreeMap::<_, (usize, i32)>::new();
                for triangle in &mesh.triangles {
                    for i in 0..3 {
                        let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                        let uses = incidence.entry((a.min(b), a.max(b))).or_default();
                        uses.0 += 1;
                        uses.1 += if a < b { 1 } else { -1 };
                    }
                }
                assert!(incidence.values().all(|uses| *uses == (2, 0)));
                for edge in &mesh.edges {
                    for pair in edge.vertices.windows(2) {
                        assert_eq!(
                            incidence.get(&(pair[0].min(pair[1]), pair[0].max(pair[1]))),
                            Some(&(2, 0))
                        );
                    }
                }
            }
            println!(
                "answered {case}: {} vertices, {} triangles",
                mesh.positions.len(),
                mesh.triangles.len()
            );
        }
        Err(error @ Error::Build(_)) if case.contains(':') => {
            assert_eq!(error.code(), "boolean.sampling-budget", "{error:?}");
            println!("refused {case}: {error:?}");
        }
        Err(Error::ResourceLimit) if case.contains(':') || case == "huge" => {
            println!("refused {case}: application.work-limit");
        }
        other => panic!("unexpected outcome for {case}: {other:?}"),
    }
}
