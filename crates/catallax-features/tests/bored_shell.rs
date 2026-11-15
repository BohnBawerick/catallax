//! Material, source names and display of a tray pierced by through bores.

use std::collections::{BTreeMap, BTreeSet};

use catallax_features::{BuildError, Built, evaluate, evaluate_incremental};
use catallax_mesh::{Mesh, MeshOptions, PreviewOptions, mesh_body, preview_mesh_body};
use catallax_numbers::{Ball, Decider, Outcome, Rational, with_operation_budget};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{RPoint3, Surface};
use serde_json::{Value, json};

const WORK: u64 = 20_000_000;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn document() -> Value {
    serde_json::from_str(include_str!("fixtures/bored-shell.json")).unwrap()
}

fn parse(document: &Value) -> Recipe {
    Recipe::parse(&document.to_string()).unwrap()
}

fn build(document: &Value) -> Built {
    let (built, work) = with_operation_budget(WORK, || {
        evaluate(&parse(document), &mut Decider::default()).unwrap()
    });
    assert!(!work.exhausted, "{work:?}");
    assert!(built.model.check(built.result).is_empty());
    built
}

fn location(built: &mut Built, point: [&str; 3], inside: bool) {
    let [x, y, z] = point.map(q);
    let result = catallax_boolean::classify_point(
        &mut built.model,
        built.result,
        &RPoint3::new(x, y, z),
        &mut Decider::default(),
    );
    let want = if inside {
        catallax_boolean::Containment::Inside
    } else {
        catallax_boolean::Containment::Outside
    };
    assert!(
        matches!(result, Outcome::Proven(value) if value == want),
        "{point:?}: {result:?}"
    );
}

fn volume(built: &Built, constant: Rational, pi: Rational) {
    let (answer, work) = with_operation_budget(WORK, || {
        catallax_measure::mass_properties(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap()
    });
    assert!(!work.exhausted, "{work:?}");
    let expected =
        Ball::from_rational(&constant, 512) + Ball::pi(512) * Ball::from_rational(&pi, 512);
    assert!(answer.value.volume.lo <= expected.lo());
    assert!(answer.value.volume.hi >= expected.hi());
    assert!(answer.value.volume.width() < q("1e-6"));
}

#[test]
fn top_open_bore_shell_preserves_the_hole_and_constant_thickness_after_edit() {
    let mut previous = None;
    for thickness in ["2", "3/2"] {
        let mut source = document();
        source["steps"][2]["thickness"] = json!(thickness);
        let mut built = build(&source);
        location(&mut built, ["20", "15", "1"], false);
        location(&mut built, ["24", "15", "6"], true);
        location(&mut built, ["26", "15", "6"], false);
        location(&mut built, ["10", "10", "1"], true);
        location(&mut built, ["10", "10", "6"], false);
        location(&mut built, ["1", "10", "6"], true);
        location(&mut built, ["10", "10", "13"], false);
        let mut names = Vec::new();
        let mut inner_bore = false;
        for face in built.model.faces_of_body(built.result) {
            for origin in built.model.face_origins(face) {
                if let Role::Inner { of } = &origin.role {
                    if of.step.as_str() == "bored" {
                        let Surface::Cylinder(cylinder) =
                            built.model.geom().surface(built.model.face(face).surface)
                        else {
                            panic!("the bore's inner wall must remain cylindrical");
                        };
                        let radius = q("3") + q(thickness);
                        assert_eq!(
                            Decider::default()
                                .rational_value(&cylinder.radius_sq)
                                .unwrap(),
                            Some(&radius * &radius)
                        );
                        inner_bore = true;
                    }
                    names.push(origin.clone());
                }
            }
        }
        assert!(
            inner_bore,
            "the expanded bore wall needs its source's inner origin"
        );
        names.sort();
        if let Some(previous) = previous {
            assert_eq!(names, previous, "thickness edits retain inner-wall names");
        }
        previous = Some(names);
        let t = q(thickness);
        let height = &q("12") - &t;
        let twice = &t + &t;
        let width = &q("40") - &twice;
        let depth = &q("30") - &twice;
        let cavity = &(&width * &depth) * &height;
        let radius = &q("3") + &t;
        volume(
            &built,
            q("14400") - cavity,
            &radius * &radius * height - q("108"),
        );
    }
}

#[test]
fn closed_two_open_caps_and_consumed_cavities_keep_the_distance_definition() {
    let mut source = document();
    source["steps"][2]["open"] = json!([]);
    let mut closed = build(&source);
    location(&mut closed, ["10", "10", "11"], true);
    location(&mut closed, ["10", "10", "6"], false);
    volume(&closed, q("6912"), q("92"));

    source["steps"][2]["open"] = json!([
        {"face":{"step":"stock","role":"start"}},
        {"face":{"step":"stock","role":"end"}}
    ]);
    let mut open = build(&source);
    location(&mut open, ["10", "10", "1"], false);
    location(&mut open, ["24", "15", "1"], true);
    volume(&open, q("3168"), q("192"));

    source["steps"][2]["open"] = json!([]);
    source["steps"][2]["thickness"] = json!("6");
    let mut consumed = build(&source);
    location(&mut consumed, ["10", "10", "6"], true);
    location(&mut consumed, ["20", "15", "6"], false);
    volume(&consumed, q("14400"), q("-108"));
}

#[test]
fn multiple_bores_keep_separate_inner_names_and_material_columns() {
    let mut source = document();
    source["steps"][1]["at"] = json!([
        {"id":"left","point":["12","15","12"],"direction":["0","0","-1"]},
        {"id":"right","point":["28","15","12"],"direction":["0","0","-1"]}
    ]);
    let mut built = build(&source);
    location(&mut built, ["12", "15", "6"], false);
    location(&mut built, ["28", "15", "6"], false);
    location(&mut built, ["16", "15", "6"], true);
    location(&mut built, ["24", "15", "6"], true);
    location(&mut built, ["20", "15", "6"], false);
    let names: BTreeSet<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|face| {
            built
                .model
                .face_origins(face)
                .iter()
                .filter_map(|origin| match &origin.role {
                    Role::Inner { of } if of.step.as_str() == "bored" => Some(origin.clone()),
                    _ => None,
                })
        })
        .collect();
    assert_eq!(names.len(), 2);
    volume(&built, q("5040"), q("284"));
}

#[test]
fn rational_rigid_placement_preserves_material_and_the_bore() {
    let mut source = document();
    source["steps"][0]["pose"] = json!({
        "rotation":["3/5","4/5","0","0"], "translation":["7","-5","11"]
    });
    source["steps"][1]["at"][0]["point"] = json!(["27", "-518/25", "551/25"]);
    source["steps"][1]["at"][0]["direction"] = json!(["0", "24/25", "7/25"]);
    let mut built = build(&source);
    location(&mut built, ["27", "-374/25", "593/25"], false);
    location(&mut built, ["31", "-374/25", "593/25"], true);
    location(&mut built, ["33", "-374/25", "593/25"], false);
}

fn unsupported(source: &Value) {
    let error = evaluate(&parse(source), &mut Decider::default()).unwrap_err();
    assert!(matches!(error, BuildError::Geometry { refusal, .. }
        if refusal.code == "distance.unsupported-envelope"));
}

#[test]
fn touching_envelopes_and_unsupported_openings_refuse_without_changing_the_source() {
    let mut source = document();
    source["steps"][2]["open"] = json!([{"face":{"step":"bored"}}]);
    unsupported(&source);
    source["steps"][2]["open"] = json!([{"face":{"step":"stock"}}]);
    unsupported(&source);
    source = document();
    source["steps"][1]["at"][0]["point"] = json!(["7", "15", "12"]);
    unsupported(&source);
    source["result"] = json!("bored");
    source["steps"].as_array_mut().unwrap().pop();
    let mut retained = build(&source);
    location(&mut retained, ["11", "15", "6"], true);
    location(&mut retained, ["7", "15", "6"], false);

    source = document();
    source["steps"][1]["at"] = json!([
        {"id":"left","point":["15","15","12"],"direction":["0","0","-1"]},
        {"id":"right","point":["25","15","12"],"direction":["0","0","-1"]}
    ]);
    unsupported(&source);
}

fn mesh_check(mesh: &Mesh, options: &MeshOptions, expected_volume: f64) {
    assert!(mesh.unavailable_controls.is_empty());
    assert!(
        mesh.faces
            .iter()
            .all(|face| face.max_chord_error <= options.chord_tolerance)
    );
    let mut uses = BTreeMap::<_, (usize, i32)>::new();
    let mut faces = BTreeSet::new();
    let mut volume = 0.0;
    for triangle in &mesh.triangles {
        faces.insert(triangle.face);
        let [a, b, c] = triangle.vertices.map(|i| mesh.positions[i]);
        let cross = [
            b[1] * c[2] - b[2] * c[1],
            b[2] * c[0] - b[0] * c[2],
            b[0] * c[1] - b[1] * c[0],
        ];
        volume += (0..3).map(|i| a[i] * cross[i]).sum::<f64>() / 6.0;
        let center: [f64; 3] = std::array::from_fn(|i| (a[i] + b[i] + c[i]) / 3.0);
        assert!(
            (center[0] - 20.0).hypot(center[1] - 15.0)
                >= 3.0 - options.chord_tolerance.to_f64_nearest()
        );
        assert!(
            !(center[2] > 2.1
                && center[0] > 2.1
                && center[0] < 37.9
                && center[1] > 2.1
                && center[1] < 27.9
                && (center[0] - 20.0).hypot(center[1] - 15.0) > 5.1),
            "filled cavity at {center:?}"
        );
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let count = uses.entry((a.min(b), a.max(b))).or_default();
            count.0 += 1;
            count.1 += if a < b { 1 } else { -1 };
        }
    }
    assert_eq!(faces.len(), mesh.faces.len());
    assert!(
        uses.values().all(|v| *v == (2, 0)),
        "display must be closed and oriented"
    );
    assert!((volume - expected_volume).abs() < expected_volume * 0.02);
}

#[test]
fn kernel_and_soliddraft_display_levels_keep_the_cavity_open_and_selectable() {
    let built = build(&document());
    let expected = 5040.0 + 142.0 * std::f64::consts::PI;
    let mut levels = vec![MeshOptions::default()];
    // SolidDraft auto quality uses 0.002 times the diagonal, scaled by 4, 1, 1/4.
    // Its adapter halves the requested chord and uses a rational lower bound of pi.
    let diagonal = (40.0_f64.powi(2) + 30.0_f64.powi(2) + 12.0_f64.powi(2)).sqrt();
    for factor in [4.0, 1.0, 0.25] {
        levels.push(MeshOptions {
            chord_tolerance: q(&(0.002 * diagonal * factor).to_string()) / q("2"),
            angular_deviation: Some(q("314159265/1200000000")),
            ..Default::default()
        });
    }
    for options in levels {
        let (mesh, work) = with_operation_budget(WORK, || {
            mesh_body(
                &built.model,
                built.result,
                &options,
                &mut Decider::default(),
            )
            .unwrap()
            .value
        });
        assert!(!work.exhausted, "{options:?}: {work:?}");
        mesh_check(&mesh, &options, expected);
        let (preview, work) = with_operation_budget(WORK, || {
            preview_mesh_body(
                &built.model,
                built.result,
                &PreviewOptions::from(&options),
                &mut Decider::default(),
            )
            .unwrap()
        });
        assert!(!work.exhausted, "{work:?}");
        assert_eq!(
            preview
                .faces
                .iter()
                .map(|f| f.face)
                .collect::<BTreeSet<_>>(),
            mesh.faces.iter().map(|f| f.face).collect()
        );
    }
    let mut bytes = Vec::new();
    let (report, work) = with_operation_budget(WORK, || {
        catallax_step::write_step(
            &built.model,
            &[built.result],
            &Default::default(),
            &mut Decider::default(),
            &mut bytes,
        )
        .unwrap()
    });
    assert!(!work.exhausted, "{work:?}");
    assert!(report.marks.is_empty());
    assert!(bytes.starts_with(b"ISO-10303-21;"));
}

fn side_document() -> Value {
    serde_json::from_str(include_str!("fixtures/bored-side-shell.json")).unwrap()
}

fn step_bytes(built: &Built) -> Vec<u8> {
    let mut bytes = Vec::new();
    let report = catallax_step::write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    assert!(report.marks.is_empty());
    assert!(bytes.starts_with(b"ISO-10303-21;"));
    bytes
}

fn side_rim(built: &Built) {
    let recipe = parse(&side_document());
    let reference = &recipe.steps[2].op;
    let catallax_recipe::Operation::Shell { open, .. } = reference else {
        unreachable!()
    };
    let faces: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|face| open[0].matches_origins(built.model.face_origins(*face)))
        .collect();
    assert_eq!(faces.len(), 1, "one finite face rings the side opening");
    let face = built.model.face(faces[0]);
    assert_eq!(
        face.loops.len(),
        2,
        "the side rim has outer and cavity wires"
    );
    for wire in &face.loops {
        let catallax_topology::Loop::Coedges(wire) = built.model.loop_(*wire) else {
            panic!("a rectangular rim needs edges")
        };
        assert_eq!(wire.len(), 4);
        for coedge in wire {
            let edge = built.model.edge(built.model.coedge(*coedge).edge);
            assert!(
                edge.start.is_some() && edge.end.is_some(),
                "rim edges are finite"
            );
        }
    }
    assert!(built.model.check(built.result).is_empty());
}

#[test]
fn side_open_bore_shell_has_a_finite_rim_and_thickness_after_edit() {
    let mut source = side_document();
    let mut prefix = source.clone();
    prefix["steps"].as_array_mut().unwrap().pop();
    prefix["result"] = json!("bored");
    let checkpoint = evaluate_incremental(&parse(&prefix), None, &mut Decider::default()).unwrap();
    let saved_source = step_bytes(checkpoint.built());
    let mut previous_names = None;
    let mut previous_build = None;
    for thickness in ["2", "3/2"] {
        source["steps"][2]["thickness"] = json!(thickness);
        // First append to a cached source, then edit the complete prior build.
        let (result, work) = with_operation_budget(WORK, || {
            evaluate_incremental(
                &parse(&source),
                previous_build.take().or_else(|| Some(checkpoint.clone())),
                &mut Decider::default(),
            )
            .unwrap()
        });
        assert!(!work.exhausted, "{work:?}");
        let mut built = result.built().clone();
        side_rim(&built);
        location(&mut built, ["1", "10", "6"], false);
        location(&mut built, ["10", "10", "1"], true);
        location(&mut built, ["10", "10", "3"], false);
        location(&mut built, ["10", "10", "11"], true);
        location(&mut built, ["10", "1", "6"], true);
        location(&mut built, ["39", "10", "6"], true);
        location(&mut built, ["20", "15", "6"], false);
        location(&mut built, ["24", "15", "6"], true);
        location(&mut built, ["26", "15", "6"], false);
        let mut names = Vec::new();
        let mut inner_bore = 0;
        for face in built.model.faces_of_body(built.result) {
            for origin in built.model.face_origins(face) {
                names.push(origin.clone());
                if let Role::Inner { of } = &origin.role
                    && of.step.as_str() == "bored"
                {
                    let Surface::Cylinder(cylinder) =
                        built.model.geom().surface(built.model.face(face).surface)
                    else {
                        unreachable!()
                    };
                    let radius = q("3") + q(thickness);
                    assert_eq!(
                        Decider::default()
                            .rational_value(&cylinder.radius_sq)
                            .unwrap(),
                        Some(&radius * &radius)
                    );
                    inner_bore += 1;
                }
            }
        }
        assert_eq!(inner_bore, 1);
        names.sort();
        if let Some(previous) = previous_names {
            assert_eq!(names, previous);
        }
        previous_names = Some(names);
        let t = q(thickness);
        let twice = &t + &t;
        let height = &q("12") - &twice;
        let width = &q("40") - &t;
        let depth = &q("30") - &twice;
        let cavity = &(&width * &depth) * &height;
        let radius = &q("3") + &t;
        volume(
            &built,
            q("14400") - cavity,
            &radius * &radius * height - q("108"),
        );
        assert_eq!(
            step_bytes(&built),
            step_bytes(&build(&source)),
            "edit equals a fresh rebuild"
        );
        assert_eq!(
            step_bytes(checkpoint.built()),
            saved_source,
            "the cached input is immutable"
        );
        previous_build = Some(result);
    }
}

#[test]
fn side_open_multiple_bores_and_rigid_placement_keep_their_holes() {
    let mut source = side_document();
    source["steps"][1]["at"] = json!([
        {"id":"left","point":["12","15","12"],"direction":["0","0","-1"]},
        {"id":"right","point":["28","15","12"],"direction":["0","0","-1"]}
    ]);
    let mut built = build(&source);
    side_rim(&built);
    location(&mut built, ["12", "15", "6"], false);
    location(&mut built, ["28", "15", "6"], false);
    location(&mut built, ["16", "15", "6"], true);
    location(&mut built, ["24", "15", "6"], true);
    volume(&built, q("6496"), q("184"));

    source = side_document();
    source["steps"][0]["pose"] = json!({
        "rotation":["3/5","4/5","0","0"], "translation":["7","-5","11"]
    });
    source["steps"][1]["at"][0]["point"] = json!(["27", "-518/25", "551/25"]);
    source["steps"][1]["at"][0]["direction"] = json!(["0", "24/25", "7/25"]);
    let mut built = build(&source);
    side_rim(&built);
    location(&mut built, ["8", "-339/25", "473/25"], false);
    location(&mut built, ["27", "-374/25", "593/25"], false);
    location(&mut built, ["31", "-374/25", "593/25"], true);
    location(&mut built, ["33", "-374/25", "593/25"], false);
    volume(&built, q("6496"), q("92"));
}

#[test]
fn side_opening_refuses_multiple_openings_and_nonprismatic_side_walls() {
    let mut source = side_document();
    source["steps"][2]["open"]
        .as_array_mut()
        .unwrap()
        .push(json!({"face":{"step":"stock","role":"end"}}));
    unsupported(&source);

    source = side_document();
    let mut shell = source["steps"].as_array_mut().unwrap().pop().unwrap();
    source["steps"].as_array_mut().unwrap().push(json!({
        "id":"drafted", "op":"draft", "body":"bored",
        "faces":[{"face":{"step":"stock","role":{"side":{"segment":"x_max"}}}}],
        "pull":["0","0","1"],
        "neutral":{"point":["0","0","0"],"normal":["0","0","1"]},
        "angle":"45"
    }));
    source["result"] = json!("drafted");
    let _drafted = build(&source);
    shell["body"] = json!("drafted");
    source["steps"].as_array_mut().unwrap().push(shell);
    source["result"] = json!("shell");
    let error = evaluate(&parse(&source), &mut Decider::default()).unwrap_err();
    assert!(matches!(error, BuildError::Geometry { refusal, .. }
        if refusal.code == "distance.unsupported-envelope"
        && refusal.detail.contains("every side plane parallel")));
}
