//! Oblique neutral generators, complete finite material and independent volume.
use catallax_boolean::{Containment, classify_point};
use catallax_features::evaluate;
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Ball, Decider, Outcome, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
use catallax_topology::DraftedCylinder;

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn recipe() -> Recipe {
    Recipe::parse(include_str!("fixtures/draft-oblique.json")).unwrap()
}

#[test]
fn both_evaluators_keep_the_authored_generators_and_finite_material() {
    let recipe = recipe();
    let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
    DraftedCylinder::check(&built.model, built.result, &mut Decider::default()).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    for (xyz, expected) in [
        (["0", "0", "1"], Location::Inside),
        (["6", "0", "1"], Location::Inside),
        (["7", "0", "1"], Location::Boundary),
        (["8", "0", "1"], Location::Outside),
        (["-9", "0", "1"], Location::Boundary),
        (["0", "6", "1"], Location::Boundary),
        (["0", "7", "1"], Location::Outside),
        (["0", "-10", "1"], Location::Boundary),
        (["0", "0", "0"], Location::Boundary),
        (["0", "0", "2"], Location::Boundary),
        (["0", "0", "3"], Location::Outside),
    ] {
        let [x, y, z] = xyz.map(q);
        let p = RPoint3::new(x, y, z);
        let field = plan.proven(&p, &mut Workspace::default());
        assert_eq!(field.location, expected, "field {xyz:?}: {field:?}");
        let brep = classify_point(&mut built.model, built.result, &p, &mut Decider::default());
        assert!(
            match expected {
                Location::Inside => matches!(brep, Outcome::Proven(Containment::Inside)),
                Location::Outside => matches!(brep, Outcome::Proven(Containment::Outside)),
                Location::Boundary => matches!(brep, Outcome::ProvenCoincident(_)),
                _ => false,
            },
            "B-rep {xyz:?}: {brep:?}"
        );
    }
    let volume = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let oracle = Ball::pi(192) * Ball::from_rational(&q("401/3"), 192);
    assert!(
        volume.value.lo <= oracle.lo() && volume.value.hi >= oracle.hi(),
        "{volume:?}, {oracle:?}"
    );
    assert!(volume.value.width() < q("1e-9"));
    assert!(built.consumed.is_empty());
}

fn check_closed(positions: &[[f64; 3]], triangles: impl Iterator<Item = [usize; 3]>) {
    let mut edges = std::collections::BTreeMap::<(usize, usize), (usize, i32)>::new();
    let mut volume = 0.0;
    for t in triangles {
        for i in 0..3 {
            let a = t[i];
            let b = t[(i + 1) % 3];
            let entry = edges.entry((a.min(b), a.max(b))).or_default();
            entry.0 += 1;
            entry.1 += if a < b { 1 } else { -1 };
        }
        let [a, b, c] = t.map(|i| positions[i]);
        volume += (a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))
            / 6.0;
    }
    assert!(edges.values().all(|entry| *entry == (2, 0)));
    let oracle = 401.0 * std::f64::consts::PI / 3.0;
    assert!(
        (volume - oracle).abs() < 0.5,
        "closed display volume {volume}, oracle {oracle}"
    );
}

#[test]
fn both_displays_close_the_periodic_seam_and_export_the_real_ruled_surface() {
    let recipe = recipe();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let mesh = catallax_mesh::mesh_body(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    check_closed(
        &mesh.value.positions,
        mesh.value.triangles.iter().map(|t| t.vertices),
    );
    assert!(
        mesh.value
            .faces
            .iter()
            .all(|f| f.max_chord_error <= q("1/100"))
    );
    let plan = Plan::compile(&recipe).unwrap();
    let field =
        catallax_field::display::mesh_plan(&plan, &Default::default(), &mut Decider::default())
            .unwrap();
    check_closed(&field.positions, field.triangles.iter().map(|t| t.vertices));
    assert!(field.max_deviation <= q("1/100"));
    let mut bytes = Vec::new();
    let exported = catallax_step::write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    assert!(exported.value.max_deviation <= q("1e-7"));
    assert!(exported.value.regularizations.is_empty());
}

fn document() -> serde_json::Value {
    serde_json::from_str(include_str!("fixtures/draft-oblique.json")).unwrap()
}

fn from_document(source: &serde_json::Value) -> Recipe {
    Recipe::parse(&source.to_string()).unwrap()
}

fn refuses(source: &serde_json::Value, code: &str) {
    let recipe = from_document(source);
    let error = evaluate(&recipe, &mut Decider::default()).unwrap_err();
    assert!(
        matches!(error, catallax_features::BuildError::Geometry { refusal, .. } if refusal.code == code)
    );
    assert_eq!(Plan::compile(&recipe).unwrap_err().code, code);
}

#[test]
fn signed_angles_pull_reversal_and_scaled_neutral_normals_keep_the_polar_law() {
    use serde_json::json;
    for (angle, pull, factor, boundary) in [
        ("45", "1", "311/3", "6"),
        ("-45", "1", "503/3", "10"),
        ("45", "-1", "503/3", "10"),
    ] {
        for normal in [json!(["1", "2", "8"]), json!(["-2", "-4", "-16"])] {
            let mut source = document();
            source["steps"][2]["neutral"]["point"][2] = json!("0");
            source["steps"][2]["neutral"]["normal"] = normal;
            source["steps"][2]["angle"] = json!(angle);
            source["steps"][2]["pull"][2] = json!(pull);
            let recipe = from_document(&source);
            let built = evaluate(&recipe, &mut Decider::default()).unwrap();
            let plan = Plan::compile(&recipe).unwrap();
            let p = RPoint3::new(q(boundary), q("0"), q("1"));
            assert_eq!(
                plan.proven(&p, &mut Workspace::default()).location,
                Location::Boundary
            );
            let volume = catallax_measure::volume(
                &built.model,
                built.result,
                &Default::default(),
                &mut Decider::default(),
            )
            .unwrap();
            let oracle = Ball::pi(192) * Ball::from_rational(&q(factor), 192);
            assert!(volume.value.lo <= oracle.lo() && volume.value.hi >= oracle.hi());
        }
    }
}

#[test]
fn a_rigidly_placed_cylinder_keeps_the_same_finite_material_and_volume() {
    use serde_json::json;
    let mut source = document();
    source["steps"][0]["plane"] =
        json!({"rotation":["1","0","1","0"],"translation":["3","-2","1"]});
    source["steps"][2]["pull"] = json!(["1", "0", "0"]);
    source["steps"][2]["neutral"] = json!({"point":["4","-2","1"],"normal":["8","2","-1"]});
    let recipe = from_document(&source);
    let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    for (xyz, expected) in [
        (["4", "-2", "1"], Location::Inside),
        (["4", "-2", "-6"], Location::Boundary),
        (["4", "-2", "-7"], Location::Outside),
    ] {
        let [x, y, z] = xyz.map(q);
        let p = RPoint3::new(x, y, z);
        assert_eq!(
            plan.proven(&p, &mut Workspace::default()).location,
            expected
        );
        let brep = classify_point(&mut built.model, built.result, &p, &mut Decider::default());
        assert!(match expected {
            Location::Inside => matches!(brep, Outcome::Proven(Containment::Inside)),
            Location::Boundary => matches!(brep, Outcome::ProvenCoincident(_)),
            _ => matches!(brep, Outcome::Proven(Containment::Outside)),
        });
    }
    let volume = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let oracle = Ball::pi(192) * Ball::from_rational(&q("401/3"), 192);
    assert!(volume.value.lo <= oracle.lo() && volume.value.hi >= oracle.hi());
}

#[test]
fn retention_uses_strict_material_picks_and_preserves_consumed_face_history() {
    use serde_json::json;
    for (height, keep, empty) in [
        ("1", Some(json!([])), true),
        ("20", None, true),
        ("20", Some(json!([["0", "0", "1"]])), false),
    ] {
        let mut source = document();
        source["steps"][2]["neutral"]["point"][2] = json!(height);
        if let Some(keep) = keep {
            source["steps"][2]["keep_regions"] = keep;
        }
        let recipe = from_document(&source);
        let built = evaluate(&recipe, &mut Decider::default()).unwrap();
        let field = Plan::compile(&recipe).unwrap();
        assert_eq!(
            built.model.body(built.result).kind == catallax_topology::BodyKind::Empty,
            empty
        );
        assert_eq!(field.is_proven_empty(), empty);
        if empty {
            assert!(
                built
                    .consumed
                    .iter()
                    .any(|f| f.source.iter().any(|o| o.step.as_str() == "base"))
            );
        }
    }
    for point in [
        json!(["7", "0", "1"]),
        json!(["0", "0", "0"]),
        json!(["20", "0", "1"]),
    ] {
        let mut source = document();
        source["steps"][2]["keep_regions"] = json!([point]);
        refuses(&source, "draft.invalid-keep-region");
    }
}

#[test]
fn folds_hollow_sections_and_unproved_remote_neutral_cuts_keep_typed_refusals() {
    use serde_json::json;
    let mut source = document();
    source["steps"][2]["neutral"]["normal"] = json!(["8", "0", "8"]);
    refuses(&source, "draft.unsupported-oblique-arrangement");
    let mut source = document();
    source["steps"][2]["neutral"]["point"][2] = json!("3");
    refuses(&source, "draft.unsupported-oblique-retention");
    let mut source = document();
    source["steps"][0]["loops"]
        .as_array_mut()
        .unwrap()
        .push(json!([{"kind":"circle","id":"bore","centre":["0","0"],"radius":"2"}]));
    refuses(&source, "draft.unsupported-oblique-annulus");
}

#[test]
fn hollow_oblique_refusals_name_the_unproved_gap_with_source_numbers() {
    use serde_json::json;
    for selected in [vec!["wall"], vec!["bore"], vec!["wall", "bore"]] {
        let mut source = document();
        source["steps"][0]["loops"]
            .as_array_mut()
            .unwrap()
            .push(json!([{"kind":"circle","id":"bore","centre":["0","0"],"radius":"2"}]));
        source["steps"][2]["faces"] = json!(
            selected
                .iter()
                .map(|segment| {
                    json!({"face":{"step":"base","role":{"side":{"segment":segment}}}})
                })
                .collect::<Vec<_>>()
        );
        let recipe = from_document(&source);
        let error = evaluate(&recipe, &mut Decider::default()).unwrap_err();
        let catallax_features::BuildError::Geometry { refusal, .. } = error else {
            panic!("expected geometric coverage refusal");
        };
        let field = Plan::compile(&recipe).unwrap_err();
        assert_eq!(refusal, field);
        assert_eq!(refusal.code, "draft.unsupported-oblique-annulus");
        for text in [
            "outer_radius_squared_mm2=[64,64]".to_owned(),
            "inner_radius_squared_mm2=[4,4]".to_owned(),
            "cap_span_mm=[2,2]".to_owned(),
            "tan_angle=[1,1]".to_owned(),
            format!("selected_outer={}", selected.contains(&"wall")),
            format!("selected_inner={}", selected.contains(&"bore")),
            "certified_radial_gap=unknown".to_owned(),
        ] {
            assert!(
                refusal.detail.contains(&text),
                "{refusal:?}: missing {text}"
            );
        }
    }
}

#[test]
fn the_quartic_axis_is_not_a_selectable_wall_or_a_false_cap_edge() {
    let plan = Plan::compile(&recipe()).unwrap();
    let wall = plan
        .faces()
        .iter()
        .enumerate()
        .find_map(|(index, face)| {
            matches!(face.surface, catallax_shapes::Surface::BSpline(_))
                .then(|| plan.result_face(index))
                .flatten()
        })
        .unwrap();
    for z in ["0", "2"] {
        let point = RPoint3::new(q("0"), q("0"), q(z));
        let query = plan
            .face_query(wall, &point, &mut Workspace::default())
            .unwrap();
        assert_eq!(query.state, catallax_field::FaceState::NotOnFace);
        let point_query = plan.proven(&point, &mut Workspace::default());
        assert!(
            matches!(point_query.direction, catallax_field::Direction::Unique(_)),
            "{point_query:?}"
        );
    }
}
