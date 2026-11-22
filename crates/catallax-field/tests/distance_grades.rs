//! Terminal grades share one definition for preview, export and mass.

use catallax_field::{
    Location, Workspace,
    lattice_stage::{GyroidSpec, LatticeSpec, Proof, Settings, Stage, TerminalField, Verdict},
    mass,
};
use catallax_numbers::{Decider, Rational, Refusal, with_operation_budget};
use catallax_recipe::{Grade, Recipe};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

const WORK: u64 = 19_999_999;

fn q(text: &str) -> Rational {
    Rational::parse(text).unwrap()
}

fn block(size: [&str; 3], origin: [&str; 3]) -> Recipe {
    Recipe::parse(
        &json!({"format":"catallax-recipe/5","name":"distance grade","steps":[
            {"id":"region","op":"block","size":size,
             "pose":{"rotation":["1","0","0","0"],"translation":origin}}
        ],"result":"region"})
        .to_string(),
    )
    .unwrap()
}

fn lattice(radius: Grade) -> LatticeSpec {
    LatticeSpec {
        cell: serde_json::from_value(json!("cubic")).unwrap(),
        size: 100.into(),
        frame: frame(),
        radius,
    }
}

fn frame() -> catallax_recipe::Pose {
    serde_json::from_value(json!({
        "rotation":["1","0","0","0"],"translation":["0","0","0"]
    }))
    .unwrap()
}

fn gyroid(level: Grade) -> GyroidSpec {
    GyroidSpec {
        period: 10.into(),
        kind: catallax_recipe::GyroidKind::Network,
        frame: frame(),
        level,
    }
}

fn distance(faces: Value, near: &str, far: &str, span: &str) -> Grade {
    serde_json::from_value(json!({"distance":{"faces":faces,"near":near,"far":far,"span":span}}))
        .unwrap()
}

fn linear(z0: &str, z1: &str, start: &str, end: &str) -> Grade {
    linear_between(["0", "0", z0], ["0", "0", z1], start, end)
}

fn linear_between(from: [&str; 3], to: [&str; 3], start: &str, end: &str) -> Grade {
    serde_json::from_value(json!({
        "linear":{"from":from,"to":to,"start":start,"end":end}
    }))
    .unwrap()
}

fn point(centre: [&str; 3], r0: &str, r1: &str, start: &str, end: &str) -> Grade {
    serde_json::from_value(
        json!({"point":{"centre":centre,"r0":r0,"r1":r1,"start":start,"end":end}}),
    )
    .unwrap()
}

fn radial(direction: [&str; 3], r0: &str, r1: &str, start: &str, end: &str) -> Grade {
    serde_json::from_value(json!({
        "radial":{"origin":["0","0","0"],"direction":direction,"r0":r0,"r1":r1,"start":start,"end":end}
    }))
    .unwrap()
}

fn settings() -> Settings {
    Settings {
        accuracy: q("1/10"),
        minimum_detail: q("1/10"),
        max_cells: 2_000_000,
        max_vertices: 1_000_000,
    }
}

fn accuracy() -> mass::Accuracy {
    let absolute = |width| catallax_measure::Accuracy {
        relative: q("0"),
        absolute: q(width),
    };
    mass::Accuracy {
        volume: absolute("1/100"),
        centroid: absolute("1/1000"),
        inertia: absolute("1/10"),
    }
}

#[derive(Debug, PartialEq, Eq)]
struct Receipts {
    preview: Vec<u8>,
    stl: Vec<u8>,
    verdict: Verdict,
    mass: String,
}

fn outputs<S: TerminalField>(body: &Recipe, spec: S) -> Receipts {
    let (result, work) = with_operation_budget(WORK, || -> Result<_, Refusal> {
        let mut dec = Decider::default();
        let stage = Stage::new(body, spec, settings(), &mut dec)?;
        let options =
            catallax_field::lattice_stage::default_preview_options(stage.plan(), &mut dec)?;
        let preview = stage.preview(&options, &mut dec)?;
        assert_eq!(preview.proof, Proof::Unproved);
        assert!(!preview.mesh.triangles.is_empty());
        let mut obj = Vec::new();
        preview.write_obj(&mut obj).unwrap();
        let (bytes, output) = stage.export(&mut dec)?.labelled_stl()?;
        Ok((stage, obj, bytes, output.verdict().clone()))
    });
    assert!(!work.exhausted && work.used < WORK, "{work:?}");
    let (stage, preview, stl, verdict) = result.unwrap();
    assert!(verdict.closed && verdict.non_self_intersecting && verdict.connected);
    assert!(verdict.max_deviation <= settings().accuracy);
    let (properties, work) = with_operation_budget(WORK, || {
        mass::properties(
            stage.plan(),
            &accuracy(),
            &mass::Options::default(),
            &mut Decider::default(),
        )
    });
    assert!(!work.exhausted && work.used < WORK, "{work:?}");
    let properties = properties.unwrap();
    assert!(properties.volume.lo > q("0"));
    Receipts {
        preview,
        stl,
        verdict,
        mass: format!("{properties:?}"),
    }
}

#[test]
#[cfg_attr(debug_assertions, ignore = "certification runs in release mode")]
fn one_load_face_has_the_same_outputs_as_the_equivalent_linear_grade() {
    let body = block(["6", "6", "10"], ["-3", "-3", "20"]);
    let selected = json!([{"face":{"step":"region","role":"start"}}]);
    let a = outputs(&body, lattice(distance(selected, "1", "2", "10")));
    let b = outputs(&body, lattice(linear("20", "30", "1", "2")));
    assert!(a == b, "load-face lattice differs from its linear control");

    let body = block(["2", "2", "2"], ["-1", "-1", "-1"]);
    let selected = json!([{"face":{"step":"region"},"near":["0","0","-1"]}]);
    let a = outputs(&body, gyroid(distance(selected, "-1/10", "1/10", "2")));
    let b = outputs(&body, gyroid(linear("-1", "1", "-1/10", "1/10")));
    assert!(a == b, "load-face gyroid differs from its linear control");
}

#[test]
#[cfg_attr(debug_assertions, ignore = "certification runs in release mode")]
fn all_region_faces_drive_the_preview_export_and_mass() {
    let body = block(["6", "6", "10"], ["-3", "-3", "20"]);
    let spec = lattice(distance(json!([{}]), "9/8", "1", "1/2"));
    let stage = Stage::new(&body, spec.clone(), settings(), &mut Decider::default()).unwrap();
    let location = |z: &str| {
        stage
            .plan()
            .proven(
                &RPoint3::new(q("17/16"), q("0"), q(z)),
                &mut Workspace::default(),
            )
            .location
    };
    assert_eq!(location("161/8"), Location::Inside);
    assert_eq!(location("25"), Location::Outside);
    outputs(&body, spec);

    let body = block(["2", "2", "2"], ["-1", "-1", "-1"]);
    let spec = gyroid(distance(json!([{}]), "-1/100", "1/100", "1/4"));
    outputs(&body, spec);
}

#[test]
fn distance_refusals_name_the_unproved_input() {
    let body = block(["2", "2", "2"], ["-1", "-1", "-1"]);
    for (faces, near, far, span, code) in [
        (json!([{}]), "1", "2", "0", "lattice.impossible"),
        (
            json!([{}]),
            "1",
            "2",
            "1",
            "lattice.grade-range-unsupported",
        ),
        (
            json!([{"near":["0","0","0"]}]),
            "1",
            "2",
            "10",
            "field.grade-face-ambiguous",
        ),
        (
            json!([{"face":{"step":"region","role":"start"},"touching":[{"step":"region","role":"end"}]}]),
            "1",
            "2",
            "10",
            "field.grade-face-missing",
        ),
    ] {
        let error = Stage::new(
            &body,
            lattice(distance(faces, near, far, span)),
            settings(),
            &mut Decider::default(),
        )
        .unwrap_err();
        assert_eq!(error.code, code, "{}", error.detail);
    }
}

#[test]
fn extrusion_faces_use_the_current_rectangular_trim() {
    let profile = |depth: &str| {
        Recipe::parse(
            &json!({
                "format":"catallax-recipe/4","name":"extruded region","steps":[
                    {"id":"profile","op":"sketch",
                     "plane":{"rotation":["1","0","0","0"],"translation":["0","0","20"]},
                     "loops":[[
                        {"kind":"line","id":"a","to":["3","-3"]},
                        {"kind":"line","id":"b","to":["3","3"]},
                        {"kind":"line","id":"c","to":["-3","3"]},
                        {"kind":"line","id":"d","to":["-3","-3"]}
                     ]]},
                    {"id":"region","op":"extrude","profile":"profile","distance":depth}
                ],"result":"region"
            })
            .to_string(),
        )
        .unwrap()
    };
    let spec = lattice(distance(
        json!([{"face":{"step":"region","role":"end"}}]),
        "9/8",
        "1",
        "1/2",
    ));
    let point = RPoint3::new(q("17/16"), q("0"), q("239/8"));
    for (depth, location) in [("10", Location::Inside), ("12", Location::Outside)] {
        let stage = Stage::new(
            &profile(depth),
            spec.clone(),
            settings(),
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(
            stage
                .plan()
                .proven(&point, &mut Workspace::default())
                .location,
            location
        );
    }
}

#[test]
#[cfg_attr(debug_assertions, ignore = "certification runs in release mode")]
fn point_and_radial_grades_export_and_measure_both_families() {
    let body = block(["6", "6", "10"], ["-3", "-3", "20"]);
    let a = outputs(
        &body,
        lattice(radial(["0", "0", "1"], "1/2", "3/2", "1", "5/4")),
    );
    let b = outputs(
        &body,
        lattice(radial(["0", "0", "-2"], "1/2", "3/2", "1", "5/4")),
    );
    assert!(
        a == b,
        "reversing and scaling the radial axis changed the outputs"
    );
    outputs(
        &body,
        lattice(point(["0", "0", "25"], "1", "2", "1", "9/8")),
    );

    let body = block(["2", "2", "2"], ["-1", "-1", "-1"]);
    outputs(
        &body,
        gyroid(radial(["1", "2", "3"], "1/2", "3/2", "-1/50", "1/50")),
    );
    outputs(
        &body,
        gyroid(point(["0", "0", "0"], "1/2", "3/2", "-1/50", "1/50")),
    );
}

#[test]
fn point_grade_proves_an_irrational_contact() {
    let body = block(["6", "6", "10"], ["-3", "-3", "20"]);
    let spec = lattice(point(["-1", "-1", "25"], "2", "4", "1", "2"));
    let stage = Stage::new(&body, spec, settings(), &mut Decider::default()).unwrap();
    // At (1, 1, 25), the distance to the centre is sqrt(8). The radius
    // is sqrt(8)/2, exactly the beam distance sqrt(2).
    for (x, expected) in [
        ("99/100", Location::Inside),
        ("1", Location::Boundary),
        ("101/100", Location::Outside),
    ] {
        let p = RPoint3::new(q(x), q("1"), q("25"));
        assert_eq!(
            stage.plan().proven(&p, &mut Workspace::default()).location,
            expected
        );
    }
}

#[test]
fn point_and_radial_grade_refusals_keep_family_limits() {
    let body = block(["2", "2", "2"], ["-1", "-1", "-1"]);
    for (grade, expected) in [
        (
            point(["0", "0", "0"], "1", "1", "1", "2"),
            "lattice.impossible",
        ),
        (
            radial(["0", "0", "0"], "0", "1", "1", "2"),
            "lattice.impossible",
        ),
        (
            radial(["0", "0", "2"], "0", "1", "1", "2"),
            "lattice.grade-range-unsupported",
        ),
    ] {
        let code = Stage::new(&body, lattice(grade), settings(), &mut Decider::default())
            .unwrap_err()
            .code;
        assert_eq!(code, expected);
    }
    for (grade, expected) in [
        (
            point(["0", "0", "0"], "-1", "1", "0", "1/10"),
            "gyroid.impossible",
        ),
        (
            radial(["1", "0", "0"], "0", "1", "0", "10001/100000"),
            "gyroid.grade-range-unsupported",
        ),
        (
            point(["0", "0", "0"], "0", "100", "1/4", "25001/100000"),
            "gyroid.grade-range-unsupported",
        ),
    ] {
        let code = Stage::new(&body, gyroid(grade), settings(), &mut Decider::default())
            .unwrap_err()
            .code;
        assert_eq!(code, expected);
    }
}

#[test]
#[cfg_attr(debug_assertions, ignore = "certification runs in release mode")]
fn oblique_linear_grades_export_and_measure_both_families() {
    let body = block(["6", "6", "10"], ["-3", "-3", "20"]);
    let spec = lattice(linear_between(
        ["-1", "0", "18"],
        ["1", "0", "32"],
        "1",
        "2",
    ));
    let stage = Stage::new(&body, spec.clone(), settings(), &mut Decider::default()).unwrap();
    for (x, expected) in [("3/2", Location::Inside), ("-3/2", Location::Outside)] {
        let point = RPoint3::new(q(x), q("0"), q("25"));
        assert_eq!(
            stage
                .plan()
                .proven(&point, &mut Workspace::default())
                .location,
            expected
        );
    }
    outputs(&body, spec);

    let body = block(["2", "2", "2"], ["-1", "-1", "-1"]);
    for (from, to) in [
        (["-4", "0", "-4"], ["4", "0", "4"]),
        (["-1/2", "0", "-1/2"], ["1/2", "0", "1/2"]),
    ] {
        outputs(&body, gyroid(linear_between(from, to, "-1/50", "1/50")));
    }
}

#[test]
fn oblique_slope_limits_use_the_exact_anchor_length() {
    let body = block(["2", "2", "2"], ["-1", "-1", "-1"]);
    let origin = ["0", "0", "0"];
    for (end, accepted) in [("7/2", true), ("3500001/1000000", false)] {
        let grade = linear_between(origin, ["3", "4", "0"], "1", end);
        let result = Stage::new(&body, lattice(grade), settings(), &mut Decider::default());
        if accepted {
            assert!(result.is_ok());
        } else {
            assert_eq!(result.unwrap_err().code, "lattice.grade-range-unsupported");
        }
    }
    for (end, accepted) in [("1/10", true), ("100001/1000000", false)] {
        let grade = linear_between(origin, ["3/5", "4/5", "0"], "0", end);
        let result = Stage::new(&body, gyroid(grade), settings(), &mut Decider::default());
        if accepted {
            assert!(result.is_ok());
        } else {
            assert_eq!(result.unwrap_err().code, "gyroid.grade-range-unsupported");
        }
    }
}
