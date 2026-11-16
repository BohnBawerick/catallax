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
fn folds_and_unproved_remote_neutral_cuts_keep_typed_refusals() {
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
    evaluate(&from_document(&source), &mut Decider::default()).unwrap();
}

#[test]
fn hollow_oblique_builds_keep_the_checked_field_gap_diagnostic() {
    use serde_json::json;
    for (selected, gap) in [
        (vec!["wall"], "[3,7]"),
        (vec!["bore"], "[29/8,51/8]"),
        (vec!["wall", "bore"], "[13/8,67/8]"),
    ] {
        let mut source = document();
        source["steps"][0]["loops"]
            .as_array_mut()
            .unwrap()
            .push(json!([{"kind":"circle","id":"bore","centre":["0","0"],"radius":"3"}]));
        source["steps"][2]["neutral"]["normal"] = json!(["1", "0", "8"]);
        source["steps"][2]["faces"] = json!(
            selected
                .iter()
                .map(|segment| {
                    json!({"face":{"step":"base","role":{"side":{"segment":segment}}}})
                })
                .collect::<Vec<_>>()
        );
        let recipe = from_document(&source);
        let built = evaluate(&recipe, &mut Decider::default()).unwrap();
        assert_eq!(built.model.faces_of_body(built.result).len(), 4);
        let refusal = Plan::compile(&recipe).unwrap_err();
        assert_eq!(refusal.code, "draft.unsupported-oblique-annulus");
        for text in [
            "outer_radius_squared_mm2=[64,64]".to_owned(),
            "inner_radius_squared_mm2=[9,9]".to_owned(),
            "cap_span_mm=[2,2]".to_owned(),
            "tan_angle=[1,1]".to_owned(),
            format!("selected_outer={}", selected.contains(&"wall")),
            format!("selected_inner={}", selected.contains(&"bore")),
            format!("certified_radial_gap_mm={gap}"),
            "separation=Separated".to_owned(),
            "cap_intersections=4 exact closed planar carrier rows".to_owned(),
            "shared_trimmed_cap_rings=unbuilt".to_owned(),
            "body_published=false".to_owned(),
        ] {
            assert!(
                refusal.detail.contains(&text),
                "{refusal:?}: missing {text}"
            );
        }
    }
}

#[test]
fn hollow_axial_neutral_names_coverage_without_dividing_by_zero() {
    use serde_json::json;
    let mut source = document();
    source["steps"][0]["loops"]
        .as_array_mut()
        .unwrap()
        .push(json!([{"kind":"circle","id":"bore","centre":["0","0"],"radius":"2"}]));
    source["steps"][2]["neutral"]["normal"] = json!(["1", "0", "0"]);
    let recipe = from_document(&source);
    let catallax_features::BuildError::Geometry { refusal, .. } =
        evaluate(&recipe, &mut Decider::default()).unwrap_err()
    else {
        panic!("expected geometric coverage refusal");
    };
    assert_eq!(refusal, Plan::compile(&recipe).unwrap_err());
    assert_eq!(refusal.code, "draft.oblique-annulus-zero-neutral-rate");
    assert!(refusal.detail.contains("neutral_axis_rate=0"));
    assert!(refusal.detail.contains("certified_radial_gap=undefined"));
    assert!(refusal.detail.contains("body_published=false"));
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

#[test]
fn annular_contacts_crossings_and_zero_neutral_rate_remain_distinct_in_both_evaluators() {
    use serde_json::json;
    for (rate, code, relation, gap) in [
        ("6", "draft.oblique-annulus-touching", "Touching", "[0,8]"),
        ("4", "draft.oblique-annulus-crossing", "Crossing", "[-1,9]"),
        ("0", "draft.oblique-annulus-zero-neutral-rate", "", ""),
    ] {
        let mut source: serde_json::Value =
            serde_json::from_str(include_str!("fixtures/draft-annular/oblique-bound.json"))
                .unwrap();
        source["steps"][2]["neutral"]["normal"] = json!(["1", "0", rate]);
        let recipe = from_document(&source);
        let catallax_features::BuildError::Geometry { refusal, .. } =
            evaluate(&recipe, &mut Decider::default()).unwrap_err()
        else {
            panic!("expected geometric refusal, never a body");
        };
        assert_eq!(refusal, Plan::compile(&recipe).unwrap_err());
        assert_eq!(refusal.code, code);
        assert!(refusal.detail.contains("body_published=false"));
        if rate != "0" {
            assert!(
                refusal.detail.contains(&format!("separation={relation}")),
                "{refusal:?}"
            );
            assert!(
                refusal
                    .detail
                    .contains(&format!("certified_radial_gap_mm={gap}")),
                "{refusal:?}"
            );
        }
    }
}

#[test]
fn separated_oblique_annulus_builds_shared_caps_and_preserves_the_hole() {
    let recipe = Recipe::parse(include_str!("fixtures/draft-annular/oblique-bound.json")).unwrap();
    let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    let faces = built.model.faces_of_body(built.result);
    assert_eq!(faces.len(), 4);
    assert!(faces.iter().all(|f| built.model.face(*f).loops.len() == 2));
    let volume = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let oracle = Ball::pi(192) * Ball::from_rational(&q("387/4"), 192);
    assert!(
        volume.value.lo <= oracle.lo() && volume.value.hi >= oracle.hi(),
        "{volume:?}"
    );
    let mesh = catallax_mesh::mesh_body(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(mesh.value.faces.len(), 4);
    assert_eq!(mesh.value.edges.len(), 4);
    let mut bytes = Vec::new();
    let step = catallax_step::write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    assert_eq!(
        (step.value.faces, step.value.edges, step.value.vertices),
        (6, 12, 8)
    );
    assert!(bytes.starts_with(b"ISO-10303-21;"));
    for (xyz, expected) in [
        (["0", "0", "1"], Containment::Outside),
        (["5", "0", "1"], Containment::Inside),
        (["0", "0", "0"], Containment::Outside),
        (["5", "0", "3"], Containment::Outside),
    ] {
        let [x, y, z] = xyz.map(q);
        let point = RPoint3::new(x, y, z);
        assert_eq!(
            classify_point(
                &mut built.model,
                built.result,
                &point,
                &mut Decider::default()
            ),
            Outcome::Proven(expected),
            "{xyz:?}"
        );
    }
}

#[test]
fn oblique_annular_radius_edit_and_recipe_reopen_keep_face_and_edge_names() {
    let mut source: serde_json::Value =
        serde_json::from_str(include_str!("fixtures/draft-annular/oblique-bound.json")).unwrap();
    let before = evaluate(&from_document(&source), &mut Decider::default()).unwrap();
    let origins: Vec<_> = before
        .model
        .faces_of_body(before.result)
        .iter()
        .flat_map(|f| before.model.face_origins(*f))
        .cloned()
        .collect();
    let edges: Vec<_> = before
        .model
        .edges_of_body(before.result)
        .iter()
        .map(|e| {
            before
                .model
                .edge_name(before.result, *e)
                .expect("named physical ring")
        })
        .collect();
    source["steps"][0]["loops"][0][0]["radius"] = serde_json::json!("9");
    let edited = from_document(&source);
    let saved = edited.to_text();
    let reopened = Recipe::parse(&saved).unwrap();
    assert_eq!(reopened, edited);
    let mut exports = Vec::new();
    for recipe in [&edited, &reopened] {
        let built = evaluate(recipe, &mut Decider::default()).unwrap();
        for origin in &origins {
            assert_eq!(built.model.faces_with_origin(built.result, origin).len(), 1);
        }
        for edge in &edges {
            assert_eq!(built.model.edges_with_name(built.result, edge).len(), 1);
        }
        assert!(built.consumed.is_empty());
        let mut bytes = Vec::new();
        catallax_step::write_step(
            &built.model,
            &[built.result],
            &Default::default(),
            &mut Decider::default(),
            &mut bytes,
        )
        .unwrap();
        exports.push(bytes);
    }
    assert_eq!(exports[0], exports[1]);
}

#[test]
fn oblique_annular_checker_rejects_a_filled_cap_reversed_hole_or_duplicate_ring() {
    use catallax_topology::{Body, BodyKind, Lump, Sense, Shell};
    let recipe = Recipe::parse(include_str!("fixtures/draft-annular/oblique-bound.json")).unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    for fault in 0..3 {
        let mut model = built.model.clone();
        let checked =
            DraftedCylinder::check(&model, built.result, &mut Decider::default()).unwrap();
        let face_id = if fault == 0 {
            checked.faces()[0]
        } else {
            checked.hole().unwrap().1
        };
        let mut face = model.face(face_id).clone();
        match fault {
            0 => {
                face.loops.pop();
            }
            1 => face.sense = Sense::Same,
            _ => face.loops[1] = face.loops[0],
        }
        let mut faces = model.faces_of_body(built.result);
        let slot = faces.iter().position(|f| *f == face_id).unwrap();
        faces[slot] = model.add_face(face);
        let shell = model.add_shell(Shell { faces });
        let lump = model.add_lump(Lump {
            shells: vec![shell],
        });
        let body = model.add_body(Body {
            kind: BodyKind::Solid,
            lumps: vec![lump],
            wires: Vec::new(),
            marks: Default::default(),
        });
        assert!(
            DraftedCylinder::check(&model, body, &mut Decider::default()).is_err(),
            "fault {fault}"
        );
    }
}
