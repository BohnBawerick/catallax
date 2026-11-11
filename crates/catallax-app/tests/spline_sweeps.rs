//! Fixed spline sweeps use the ordinary named display, property and STEP paths.

use catallax_features::BuildError;
use catallax_measure::Accuracy;
use catallax_mesh::MeshOptions;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Recipe, Role};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet};

fn definition(spatial: bool) -> Value {
    serde_json::from_str(if spatial {
        include_str!("../../catallax-features/tests/fixtures/spline-sweep-spatial.json")
    } else {
        include_str!("../../catallax-features/tests/fixtures/spline-sweep-planar.json")
    })
    .unwrap()
}
fn recipe(value: &Value) -> Recipe {
    Recipe::parse(&value.to_string()).unwrap()
}
fn build(value: &Value) -> catallax_app::ApplicationBody {
    catallax_app::build(
        &recipe(value),
        &MeshOptions::default(),
        &mut Decider::default(),
    )
    .unwrap()
}
fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

#[test]
fn planar_and_spatial_paths_have_named_watertight_display_requested_mass_and_step() {
    for spatial in [false, true] {
        let app = build(&definition(spatial));
        let model = &app.built.model;
        assert!(model.check(app.built.result).is_empty());
        assert_eq!(
            (
                app.names.faces.len(),
                app.names.edges.len(),
                app.names.vertices.len()
            ),
            (6, 12, 8)
        );
        assert!(app.display.marks.is_empty());
        let mesh = &app.display.value;
        assert_eq!(
            mesh.faces.iter().map(|f| f.face).collect::<BTreeSet<_>>(),
            app.names.faces.keys().copied().collect()
        );
        assert_eq!(
            mesh.edges.iter().map(|e| e.edge).collect::<BTreeSet<_>>(),
            app.names.edges.keys().copied().collect()
        );
        let mut uses = BTreeMap::<_, (usize, i32)>::new();
        for triangle in &mesh.triangles {
            assert!(app.names.faces.contains_key(&triangle.face));
            for i in 0..3 {
                let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                let entry = uses.entry((a.min(b), a.max(b))).or_default();
                entry.0 += 1;
                entry.1 += if a < b { 1 } else { -1 };
            }
        }
        assert!(uses.values().all(|&entry| entry == (2, 0)));
        assert!(
            mesh.faces
                .iter()
                .all(|face| face.max_chord_error <= MeshOptions::default().chord_tolerance)
        );
        let sides: BTreeSet<_> = app
            .names
            .faces
            .values()
            .flatten()
            .filter_map(|o| match &o.role {
                Role::Side { segment, path } => {
                    assert_eq!(path.as_ref().unwrap().as_str(), "travel");
                    Some(segment.as_str())
                }
                _ => None,
            })
            .collect();
        assert_eq!(sides, BTreeSet::from(["a", "b", "c", "d"]));
        let accuracy = Accuracy {
            relative: Rational::zero(),
            absolute: q("1e-6"),
        };
        let mass = app
            .mass_properties(&accuracy, &mut Decider::default())
            .unwrap();
        assert!(mass.value.volume.contains(&q("24")));
        assert!(mass.value.volume.width() <= accuracy.absolute);
        let centroid = mass.value.centroid.unwrap();
        for (actual, expected) in
            centroid
                .iter()
                .zip(["5/2", if spatial { "3/2" } else { "1" }, "3/2"])
        {
            assert!(
                actual.contains(&q(expected)),
                "{actual:?} excludes {expected}"
            );
            assert!(actual.width() <= accuracy.absolute);
        }
        let mut bytes = Vec::new();
        let options = catallax_step::StepOptions::default();
        let step = catallax_step::write_step(
            model,
            &[app.built.result],
            &options,
            &mut Decider::default(),
            &mut bytes,
        )
        .unwrap();
        assert!(step.marks.is_empty());
        assert_eq!(step.value.faces, 6);
        assert_eq!(step.value.face_map.len(), 6);
        assert!(step.value.max_deviation <= options.curve_tolerance);
        let text = std::str::from_utf8(&bytes).unwrap();
        assert!(text.contains("MANIFOLD_SOLID_BREP"));
        assert_eq!(text.matches("ADVANCED_FACE").count(), 6);
    }
}

#[test]
fn path_edits_keep_names_and_path_start_is_not_a_profile_anchor() {
    let value = definition(true);
    let before = build(&value);
    let mut edited = value.clone();
    edited["steps"][1]["segments"][0]["control_points"][1][0] = json!("2");
    let after = build(&edited);
    assert_eq!(
        before.names.faces.values().collect::<BTreeSet<_>>(),
        after.names.faces.values().collect::<BTreeSet<_>>()
    );
    assert_eq!(before.names.edges.len(), after.names.edges.len());
    for name in before.names.edges.values() {
        assert!(
            after
                .names
                .edges
                .values()
                .any(|candidate| candidate == name)
        );
    }
    assert_eq!(before.names.vertices.len(), after.names.vertices.len());
    for name in before.names.vertices.values() {
        assert!(
            after
                .names
                .vertices
                .values()
                .any(|candidate| candidate == name)
        );
    }
    let mut shifted = value;
    shifted["steps"][1]["start"] = json!(["100", "200", "300"]);
    for point in shifted["steps"][1]["segments"][0]["control_points"]
        .as_array_mut()
        .unwrap()
    {
        for (axis, delta) in [100, 200, 300].into_iter().enumerate() {
            point[axis] =
                json!((q(point[axis].as_str().unwrap()) + Rational::from_i64(delta)).to_string());
        }
    }
    let translated = build(&shifted);
    assert_eq!(before.names.faces, translated.names.faces);
    assert_eq!(
        before.display.value.positions,
        translated.display.value.positions
    );
}

#[test]
fn self_crossing_and_tight_return_paths_need_a_typed_envelope_proof() {
    for controls in [
        json!([
            ["0", "0", "0"],
            ["1", "0", "2"],
            ["-1", "0", "2"],
            ["0", "0", "0"]
        ]),
        json!([
            ["0", "0", "0"],
            ["0", "0", "1/10"],
            ["1/10", "0", "1/10"],
            ["1/10", "0", "0"]
        ]),
    ] {
        let mut value = definition(false);
        value["steps"][1]["segments"][0]["control_points"] = controls;
        let error = catallax_app::build(
            &recipe(&value),
            &MeshOptions::default(),
            &mut Decider::default(),
        )
        .unwrap_err();
        assert!(
            matches!(error, catallax_app::Error::Build(BuildError::Geometry { refusal, .. }) if refusal.code == "sweep.spline-admissibility-unproven")
        );
    }
}

#[test]
fn positive_rational_weights_internal_knots_and_reverse_travel_are_supported() {
    for reverse in [false, true] {
        let mut value = definition(true);
        let path = &mut value["steps"][1]["segments"][0];
        path["degree"] = json!(2);
        path["control_points"] = json!([
            ["0", "0", "0"],
            ["1", "0", "1"],
            ["1", "1", "2"],
            ["0", "1", "3"]
        ]);
        path["weights"] = json!(["1", "2", "2", "1"]);
        path["knots"] = json!(["0", "0", "0", "1/2", "1", "1", "1"]);
        if reverse {
            for point in path["control_points"].as_array_mut().unwrap() {
                point[2] = json!((-q(point[2].as_str().unwrap())).to_string());
            }
        }
        let app = build(&value);
        assert!(app.built.model.check(app.built.result).is_empty());
        assert_eq!(app.names.faces.len(), 6);
        let volume = catallax_measure::volume(
            &app.built.model,
            app.built.result,
            &Accuracy::default(),
            &mut Decider::default(),
        )
        .unwrap();
        assert!(volume.value.contains(&q("24")));
    }
}

#[test]
fn crossed_linear_spline_is_refused_before_publishing_a_body() {
    let mut value = definition(false);
    let path = &mut value["steps"][1]["segments"][0];
    path["degree"] = json!(1);
    path["control_points"] = json!([
        ["0", "0", "0"],
        ["1", "0", "1"],
        ["0", "0", "1"],
        ["1", "0", "0"]
    ]);
    path["weights"] = json!(["1", "1", "1", "1"]);
    path["knots"] = json!(["0", "0", "1/3", "2/3", "1", "1"]);
    let error = catallax_features::evaluate(&recipe(&value), &mut Decider::default()).unwrap_err();
    assert!(
        matches!(error, BuildError::Geometry { refusal, .. } if refusal.code == "sweep.spline-admissibility-unproven")
    );
}

#[test]
fn polygon_holes_and_collinear_origins_survive_the_atlas() {
    let mut value = definition(true);
    value["steps"][0]["loops"][0] = json!([
        {"kind":"line","id":"a","to":["2","0"]},
        {"kind":"line","id":"a_split","to":["4","0"]},
        {"kind":"line","id":"b","to":["4","2"]},
        {"kind":"line","id":"c","to":["0","2"]},
        {"kind":"line","id":"d","to":["0","0"]}
    ]);
    value["steps"][0]["loops"]
        .as_array_mut()
        .unwrap()
        .push(json!([
            {"kind":"line","id":"h0","to":["1","1/2"]},
            {"kind":"line","id":"h1","to":["2","1/2"]},
            {"kind":"line","id":"h2","to":["2","3/2"]},
            {"kind":"line","id":"h3","to":["1","3/2"]}
        ]));
    let app = build(&value);
    assert_eq!(app.names.faces.len(), 10);
    assert!(app.names.faces.values().any(|origins| {
        let ids: BTreeSet<_> = origins
            .iter()
            .filter_map(|origin| match &origin.role {
                Role::Side { segment, .. } => Some(segment.as_str()),
                _ => None,
            })
            .collect();
        ids == BTreeSet::from(["a", "a_split"])
    }));
    let volume = catallax_measure::volume(
        &app.built.model,
        app.built.result,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(volume.value.contains(&q("21")));
}
