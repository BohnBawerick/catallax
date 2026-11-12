//! Complete revolved shells prove their finite walls, corners and opening rims.
use catallax_boolean::{Containment, classify_point};
use catallax_features::{Built, evaluate};
use catallax_numbers::{Decider, Rational, with_operation_budget, with_operation_log};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
use catallax_topology::BodyKind;
use serde_json::{Value, json};

fn source(name: &str) -> Value {
    serde_json::from_str(match name {
        "frustum" => include_str!("fixtures/revolved-shell/revolved-frustum.json"),
        "joined" => include_str!("fixtures/revolved-shell/revolved-cone-cylinder.json"),
        "cylinder" => include_str!("fixtures/revolved-shell/revolved-cylinder.json"),
        _ => unreachable!(),
    })
    .unwrap()
}
fn build(source: &Value, path: Option<&str>) -> Built {
    let recipe = Recipe::parse(&source.to_string()).unwrap();
    let ((built, log), work) = with_operation_budget(1_000_000, || {
        with_operation_log(|| evaluate(&recipe, &mut Decider::default()))
    });
    let built = built.unwrap();
    assert!(!work.exhausted);
    assert!(
        built.model.check(built.result).is_empty(),
        "{:?}",
        built.model.check(built.result)
    );
    if let Some(path) = path {
        assert_eq!(
            log.iter().find(|entry| entry.id == "result").unwrap().path,
            Some(path)
        );
    }
    built
}
fn point(coordinates: [&str; 3], posed: bool) -> RPoint3 {
    let [x, y, z] = coordinates.map(|value| Rational::parse(value).unwrap());
    if posed {
        RPoint3::new(
            Rational::from_i64(3) - y,
            x - Rational::from_i64(2),
            z + Rational::from_i64(5),
        )
    } else {
        RPoint3::new(x, y, z)
    }
}
fn locations(built: &mut Built, posed: bool, points: &[([&str; 3], Option<Containment>)]) {
    for (coordinates, expected) in points {
        let actual = classify_point(
            &mut built.model,
            built.result,
            &point(*coordinates, posed),
            &mut Decider::default(),
        );
        match expected {
            Some(expected) => assert_eq!(
                actual.proven(),
                Some(expected),
                "at {coordinates:?}: {actual:?}"
            ),
            None => assert!(actual.is_coincident(), "at {coordinates:?}: {actual:?}"),
        }
    }
}
fn pose(source: &mut Value) {
    source["steps"][0]["plane"] =
        json!({"rotation":["1","0","0","1"],"translation":["3","-2","5"]});
    source["steps"][1]["axis"] = json!({"point":["3","-2","5"],"direction":["-3","0","0"]});
}

#[test]
fn open_revolved_frustum_sews_true_inset_wall_to_the_closed_planar_rim() {
    for posed in [false, true] {
        let mut source = source("frustum");
        if posed {
            pose(&mut source);
        }
        let mut built = build(&source, Some("distance.shell.cone-sewn"));
        locations(
            &mut built,
            posed,
            &[
                (["11/2", "5", "0"], Some(Containment::Inside)),
                (["27/5", "5", "0"], Some(Containment::Outside)),
                (["0", "1/4", "0"], Some(Containment::Inside)),
                (["0", "1", "0"], Some(Containment::Outside)),
                (["31/4", "10", "0"], None),
                (["7", "10", "0"], Some(Containment::Outside)),
                (["8", "10", "0"], None),
                (["0", "-1/4", "0"], Some(Containment::Outside)),
            ],
        );
    }
}

#[test]
fn cone_cylinder_join_rebuilds_the_exact_inset_meeting() {
    let mut built = build(&source("joined"), Some("distance.shell.cone-sewn"));
    locations(
        &mut built,
        false,
        &[
            (["19/4", "7", "0"], Some(Containment::Inside)),
            (["17/4", "7", "0"], Some(Containment::Outside)),
            (["449/100", "101/20", "0"], Some(Containment::Inside)),
            (["449/100", "51/10", "0"], Some(Containment::Outside)),
            (["19/4", "10", "0"], None),
            (["4", "10", "0"], Some(Containment::Outside)),
        ],
    );
}

#[test]
fn cylindrical_rims_share_the_inset_edges_at_both_openings() {
    let mut source = source("cylinder");
    source["steps"][2]["open"]
        .as_array_mut()
        .unwrap()
        .push(json!({"face":{"step":"turned","role":{"side":{"segment":"base"}}}}));
    let mut built = build(&source, Some("distance.shell.cylinder-sewn"));
    assert_eq!(built.model.faces_of_body(built.result).len(), 4);
    locations(
        &mut built,
        false,
        &[
            (["15/4", "5", "0"], Some(Containment::Inside)),
            (["13/4", "5", "0"], Some(Containment::Outside)),
            (["15/4", "0", "0"], None),
            (["0", "0", "0"], Some(Containment::Outside)),
            (["15/4", "10", "0"], None),
        ],
    );
    assert_eq!(built.model.body(built.result).lumps.len(), 1);
    let lump = built.model.body(built.result).lumps[0];
    assert_eq!(built.model.lump(lump).shells.len(), 1);
}

#[test]
fn closed_revolved_shell_keeps_a_reversed_void_shell() {
    let mut source = source("frustum");
    source["steps"][2]["open"] = json!([]);
    let mut built = build(&source, Some("distance.shell.cone-sewn"));
    let lump = built.model.body(built.result).lumps[0];
    assert_eq!(built.model.lump(lump).shells.len(), 2);
    locations(
        &mut built,
        false,
        &[
            (["0", "1/4", "0"], Some(Containment::Inside)),
            (["0", "5", "0"], Some(Containment::Outside)),
            (["0", "39/4", "0"], Some(Containment::Inside)),
        ],
    );
}

#[test]
fn consumed_cavity_and_no_kept_faces_follow_the_distance_definition() {
    let mut source = source("frustum");
    source["steps"][2]["thickness"] = json!("20");
    let mut built = build(&source, None);
    locations(
        &mut built,
        false,
        &[(["0", "5", "0"], Some(Containment::Inside))],
    );
    source["steps"][2]["open"] = json!([{"face":{"step":"turned"}}]);
    let empty = build(&source, None);
    assert_eq!(empty.model.body(empty.result).kind, BodyKind::Empty);
}

#[test]
fn opening_the_narrow_conical_end_keeps_the_finite_projection_refusal() {
    let mut source = source("frustum");
    source["steps"][2]["open"] =
        json!([{"face":{"step":"turned","role":{"side":{"segment":"base"}}}}]);
    let recipe = Recipe::parse(&source.to_string()).unwrap();
    let catallax_features::BuildError::Geometry { refusal, .. } =
        evaluate(&recipe, &mut Decider::default()).unwrap_err()
    else {
        panic!("expected geometric refusal");
    };
    assert_eq!(refusal.code, "distance.cone-domain");
}

fn radial_paths(kind: &str) {
    let stock = if kind == "sphere" {
        json!({"id":"stock","op":"sphere","radius":"2","pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}})
    } else {
        json!({"id":"stock","op":"torus","major":"5","minor":"1","pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}})
    };
    let operations = [
        json!({"id":"result","op":"offset","body":"stock","distance":"1/4"}),
        json!({"id":"result","op":"shell","body":"stock","open":[],"thickness":"1/4"}),
        json!({"id":"result","op":"thicken","body":"stock","faces":[{"face":{"step":"stock"}}],"side":"front","thickness":"1/2"}),
    ];
    for (index, operation) in operations.into_iter().enumerate() {
        let name = operation["op"].as_str().unwrap();
        let path = format!(
            "distance.{name}.{}",
            if kind == "sphere" {
                "sphere-radius"
            } else {
                "torus-minor-radius"
            }
        );
        let source = json!({"format":"catallax-recipe/4","name":"radial distance proof","steps":[stock.clone(),operation],"result":"result"});
        let mut built = build(&source, Some(&path));
        let mut dec = Decider::default();
        let mut radii = Vec::new();
        for face in built.model.faces_of_body(built.result) {
            let surface = built.model.geom().surface(built.model.face(face).surface);
            let radius = match surface {
                catallax_shapes::Surface::Sphere(sphere) => {
                    assert!(
                        dec.sign(
                            &sphere
                                .center
                                .minus(&catallax_shapes::Point3::origin())
                                .norm_squared()
                        )
                        .is_coincident()
                    );
                    &sphere.radius
                }
                catallax_shapes::Surface::Torus(torus) => {
                    assert!(
                        dec.sign(&(&torus.major - &catallax_numbers::Real::from_i64(5)))
                            .is_coincident()
                    );
                    assert!(
                        dec.sign(
                            &torus
                                .center
                                .minus(&catallax_shapes::Point3::origin())
                                .norm_squared()
                        )
                        .is_coincident()
                    );
                    &torus.minor
                }
                _ => panic!("a radial distance path changed carrier family"),
            };
            radii.push(dec.rational_value(radius).unwrap().unwrap());
        }
        radii.sort();
        let expected: &[&str] = match (kind, index) {
            ("sphere", 0) => &["9/4"],
            ("sphere", 1) => &["7/4", "2"],
            ("sphere", _) => &["2", "5/2"],
            (_, 0) => &["5/4"],
            (_, 1) => &["3/4", "1"],
            _ => &["1", "3/2"],
        };
        assert_eq!(
            radii,
            expected
                .iter()
                .map(|value| Rational::parse(value).unwrap())
                .collect::<Vec<_>>()
        );
        let points = match (kind, index) {
            ("sphere", 0) => [
                ("11/5", Some(Containment::Inside)),
                ("23/10", Some(Containment::Outside)),
                ("9/4", None),
            ],
            ("sphere", 1) => [
                ("19/10", Some(Containment::Inside)),
                ("3/2", Some(Containment::Outside)),
                ("7/4", None),
            ],
            ("sphere", _) => [
                ("9/4", Some(Containment::Inside)),
                ("7/4", Some(Containment::Outside)),
                ("5/2", None),
            ],
            (_, 0) => [
                ("31/5", Some(Containment::Inside)),
                ("63/10", Some(Containment::Outside)),
                ("25/4", None),
            ],
            (_, 1) => [
                ("59/10", Some(Containment::Inside)),
                ("11/2", Some(Containment::Outside)),
                ("23/4", None),
            ],
            _ => [
                ("25/4", Some(Containment::Inside)),
                ("23/4", Some(Containment::Outside)),
                ("13/2", None),
            ],
        };
        for (x, expected) in points {
            locations(&mut built, false, &[([x, "0", "0"], expected)]);
        }
    }
}

#[test]
fn sphere_radius_paths_have_zero_analytic_deviation_and_the_defined_material() {
    radial_paths("sphere");
}

#[test]
fn torus_radius_paths_preserve_the_core_and_have_zero_analytic_deviation() {
    radial_paths("torus");
}

#[test]
fn pointed_revolved_shell_keeps_distinct_outer_and_inset_apices() {
    let source: Value = serde_json::from_str(include_str!(
        "fixtures/revolved-shell/pointed-open-cone.json"
    ))
    .unwrap();
    let mut built = build(&source, Some("distance.shell.cone-sewn"));
    assert_eq!(built.model.vertices_of_body(built.result).len(), 2);
    locations(
        &mut built,
        false,
        &[
            (["7/4", "5", "0"], Some(Containment::Inside)),
            (["5/4", "5", "0"], Some(Containment::Outside)),
            (["0", "9", "0"], Some(Containment::Inside)),
            (["0", "8", "0"], Some(Containment::Outside)),
            (["15/4", "0", "0"], None),
            (["3", "0", "0"], Some(Containment::Outside)),
        ],
    );
}

#[test]
fn revolved_sphere_cylinder_shell_uses_the_finite_projection_proof() {
    let source: Value = serde_json::from_str(include_str!(
        "fixtures/revolved-shell/revolved-spherical-cap.json"
    ))
    .unwrap();
    let mut built = build(&source, Some("distance.shell.circular-sewn"));
    locations(
        &mut built,
        false,
        &[
            (["21/10", "4", "0"], Some(Containment::Inside)),
            (["2", "4", "0"], Some(Containment::Outside)),
            (["249/100", "15/4", "0"], Some(Containment::Inside)),
            (["249/100", "373/100", "0"], Some(Containment::Outside)),
            (["11/4", "0", "0"], None),
            (["2", "0", "0"], Some(Containment::Outside)),
            (["0", "19/4", "0"], Some(Containment::Inside)),
            (["0", "17/4", "0"], Some(Containment::Outside)),
        ],
    );
}

#[test]
fn nonplanar_revolved_opening_keeps_a_typed_assembly_refusal() {
    let mut source = source("cylinder");
    source["steps"][2]["open"] =
        json!([{"face":{"step":"turned","role":{"side":{"segment":"wall"}}}}]);
    let recipe = Recipe::parse(&source.to_string()).unwrap();
    let catallax_features::BuildError::Geometry { refusal, .. } =
        evaluate(&recipe, &mut Decider::default()).unwrap_err()
    else {
        panic!("expected geometric refusal");
    };
    assert_eq!(refusal.code, "distance.analytic-rim-domain");
}
