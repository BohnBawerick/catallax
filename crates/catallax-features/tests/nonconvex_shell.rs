//! Polygonal shell material, boundaries and exact volume oracles.

use std::collections::{BTreeMap, BTreeSet};

use catallax_features::{Built, evaluate};
use catallax_numbers::{Ball, Decider, Outcome, Rational, Real, with_operation_budget};
use catallax_recipe::{EdgeName, Origin, Recipe, Role};
use catallax_shapes::{RPoint3, Surface};
use catallax_topology::Sense;
use serde_json::{Value, json};

const WORK: u64 = 100_000_000;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn point(xyz: [&str; 3]) -> RPoint3 {
    let [x, y, z] = xyz.map(q);
    RPoint3::new(x, y, z)
}
fn document() -> Value {
    serde_json::from_str(include_str!("fixtures/nonconvex-shell.json")).unwrap()
}
fn build(source: &Value) -> Built {
    let recipe = Recipe::parse(&source.to_string()).unwrap();
    let (built, work) = with_operation_budget(WORK, || {
        evaluate(&recipe, &mut Decider::default()).unwrap()
    });
    assert!(!work.exhausted, "{work:?}");
    assert!(built.model.check(built.result).is_empty());
    built
}
fn location(built: &mut Built, xyz: [&str; 3], inside: bool) {
    let (answer, work) = with_operation_budget(WORK, || {
        catallax_boolean::classify_point(
            &mut built.model, built.result, &point(xyz), &mut Decider::default(),
        )
    });
    assert!(!work.exhausted, "{work:?}");
    let expected = if inside {
        catallax_boolean::Containment::Inside
    } else {
        catallax_boolean::Containment::Outside
    };
    assert!(matches!(answer, Outcome::Proven(value) if value == expected),
        "{xyz:?}: {answer:?}");
}
fn boundary(built: &mut Built, xyz: [&str; 3]) {
    let (answer, work) = with_operation_budget(WORK, || {
        catallax_boolean::classify_point(
            &mut built.model, built.result, &point(xyz), &mut Decider::default(),
        )
    });
    assert!(!work.exhausted, "{work:?}");
    assert!(matches!(answer, Outcome::ProvenCoincident(_)), "{xyz:?}: {answer:?}");
}
fn volume(built: &Built, constant: &str, pi: &str) {
    let (mass, work) = with_operation_budget(WORK, || {
        catallax_measure::mass_properties(
            &built.model, built.result, &Default::default(), &mut Decider::default(),
        ).unwrap()
    });
    assert!(!work.exhausted, "{work:?}");
    let ball = |s| Ball::from_rational(&q(s), 512);
    let expected = ball(constant) + Ball::pi(512) * ball(pi);
    let actual = mass.value.volume;
    assert!(actual.lo <= expected.lo() && actual.hi >= expected.hi(), "{actual:?}");
    assert!(actual.width() < q("1e-6"), "{actual:?}");
}
fn chain(points: &[[i64; 2]], prefix: &str) -> Value {
    json!(points.iter().enumerate().map(|(i, _)| json!({
        "kind":"line", "id":format!("{prefix}-{i}"),
        "to":points[(i + 1) % points.len()].map(|n| n.to_string())
    })).collect::<Vec<_>>())
}
fn open_except(source: &mut Value, segment: &str) {
    let mut open = vec![json!({"face":{"step":"stock","role":"start"}}),
                        json!({"face":{"step":"stock","role":"end"}})];
    for edge in source["steps"][0]["loops"].as_array().unwrap().iter()
        .flat_map(|cycle| cycle.as_array().unwrap())
    {
        if edge["id"] != segment {
            open.push(json!({"face":{"step":"stock","role":{"side":{"segment":edge["id"]}}}}));
        }
    }
    source["steps"][2]["open"] = json!(open);
}
fn holed_document() -> Value {
    let mut source = document();
    source["steps"][0]["loops"] = json!([
        chain(&[[0,0],[10,0],[10,10],[0,10]], "outer"),
        chain(&[[3,3],[3,7],[7,7],[7,3]], "hole")
    ]);
    source["steps"][2]["open"] = json!([
        {"face":{"step":"stock","role":"start"}},
        {"face":{"step":"stock","role":"end"}}
    ]);
    source
}

#[test]
fn reentrant_tray_keeps_finite_corner_material_and_the_original_volume() {
    let mut built = build(&document());
    for (point, inside) in [
        (["2.1", "2.1", "2"], false),
        (["2.5", "2.5", "2"], true),
        (["4", "1.5", "2"], false),
        (["4", "0.5", "2"], true),
        (["4", "1.5", "0.5"], true),
        (["4", "4", "2"], false),
    ] {
        location(&mut built, point, inside);
    }
    boundary(&mut built, ["12/5", "11/5", "2"]);
    volume(&built, "120", "3/4");
}

#[test]
fn one_kept_face_retains_material_on_both_normal_sides() {
    let mut source = document();
    open_except(&mut source, "notch-bottom");
    let mut built = build(&source);
    for (point, inside) in [
        (["2.5", "3.5", "2"], true),
        (["2", "3.5", "2"], false),
        (["4", "2.5", "2"], true),
        (["1", "2.5", "2"], false),
    ] {
        location(&mut built, point, inside);
    }
    boundary(&mut built, ["12/5", "19/5", "2"]);
    volume(&built, "20", "2");
}

#[test]
fn material_across_a_gap_has_the_positive_parallel_inner_face() {
    let mut source = document();
    source["steps"][0]["loops"] = json!([chain(
        &[[0,0],[8,0],[8,8],[5,8],[5,2],[3,2],[3,8],[0,8]], "u",
    )]);
    open_except(&mut source, "u-5");
    source["steps"][2]["thickness"] = json!("3");
    let mut built = build(&source);
    for (point, inside) in [
        (["5.5", "5", "2"], true),
        (["6.5", "5", "2"], false),
        (["1", "5", "2"], true),
        (["4", "5", "2"], false),
    ] {
        location(&mut built, point, inside);
    }
    boundary(&mut built, ["6", "5", "2"]);
    let inner: Origin = serde_json::from_value(json!({
        "step":"shell","role":{"inner":{"of":{
            "step":"stock","role":{"side":{"segment":"u-5"}}
        }}}
    })).unwrap();
    let target = point(["6", "5", "2"]).to_real();
    let mut dec = Decider::default();
    let mut positive = false;
    for id in built.model.faces_of_body(built.result) {
        if !built.model.face_origins(id).contains(&inner) {
            continue;
        }
        let face = built.model.face(id);
        let Surface::Plane(plane) = built.model.geom().surface(face.surface) else {
            panic!("a planar inner image must remain planar");
        };
        let height = target.minus(&plane.origin).dot(&plane.normal);
        if dec.rational_value(&height).unwrap() != Some(Rational::zero()) {
            continue;
        }
        let normal = if face.sense == Sense::Same {
            plane.normal.clone()
        } else {
            plane.normal.scaled(&Real::from_i64(-1))
        };
        assert!(dec.enclose_at(&normal.x, 128).unwrap().lo() > Rational::zero());
        assert_eq!(dec.rational_value(&normal.y).unwrap(), Some(Rational::zero()));
        assert_eq!(dec.rational_value(&normal.z).unwrap(), Some(Rational::zero()));
        positive = true;
    }
    assert!(positive, "the positive parallel image must retain its source name");
}

#[test]
fn holed_and_consumed_shells_preserve_the_void_and_original_volumes() {
    let mut source = holed_document();
    source["steps"][2]["thickness"] = json!("3/2");
    let mut holed = build(&source);
    location(&mut holed, ["1.6", "1.6", "2"], false);
    location(&mut holed, ["1.4", "5", "2"], true);
    location(&mut holed, ["5", "5", "2"], false);
    volume(&holed, "300", "9");

    source["steps"][2]["thickness"] = json!("2");
    let mut consumed = build(&source);
    location(&mut consumed, ["1.6", "1.6", "2"], true);
    location(&mut consumed, ["5", "5", "2"], false);
    boundary(&mut consumed, ["0", "5", "2"]);
    volume(&consumed, "336", "0");
    for face in consumed.model.faces_of_body(consumed.result) {
        for origin in consumed.model.face_origins(face) {
            assert!(!(origin.step.as_str() == "shell" && matches!(
                &origin.role, Role::Inner { .. } | Role::Round { .. } | Role::RoundCorner { .. }
            )), "a consumed envelope cannot retain a final generated face: {origin:?}");
        }
    }
}

fn final_names(built: &Built) -> (BTreeSet<Origin>, BTreeSet<[Origin; 2]>) {
    let mut faces = BTreeSet::new();
    for face in built.model.faces_of_body(built.result) {
        let origins = built.model.face_origins(face);
        origins.first().expect("a retained physical face needs source provenance");
        faces.extend(origins.iter().cloned());
    }
    let edges = built.model.edges_of_body(built.result).into_iter().map(|edge| {
        built.model.edge_name(built.result, edge)
            .unwrap_or_else(|| panic!("a retained physical edge needs incident names: {edge:?}"))
            .between().clone()
    }).collect();
    (faces, edges)
}

#[test]
fn final_face_and_edge_names_resolve_after_a_thickness_edit() {
    let mut source = holed_document();
    source["steps"][2]["thickness"] = json!("1");
    let before = build(&source);
    let saved = final_names(&before);
    let corner: Origin = serde_json::from_value(json!({
        "step":"shell","role":{"round":{"edge":{"between":[
            {"step":"stock","role":{"side":{"segment":"hole-0"}}},
            {"step":"stock","role":{"side":{"segment":"hole-3"}}}
        ]}}}
    })).unwrap();
    let start: Origin = serde_json::from_value(json!({"step":"stock","role":"start"})).unwrap();
    let edge = EdgeName::new(start, corner.clone()).unwrap();
    assert!(saved.1.contains(edge.between()));
    volume(&before, "208", "4");

    source["steps"][2]["thickness"] = json!("5/4");
    let mut after = build(&source);
    assert_eq!(final_names(&after), saved, "saved physical names survive the edit");
    boundary(&mut after, ["9/4", "2", "2"]);
    volume(&after, "255", "25/4");
    let inner: Origin = serde_json::from_value(json!({
        "step":"shell","role":{"inner":{"of":{
            "step":"stock","role":{"side":{"segment":"hole-3"}}
        }}}
    })).unwrap();
    let mut dec = Decider::default();
    let inner_faces = after.model.faces_with_origin(after.result, &inner);
    inner_faces.first().expect("the saved inner-face reference must resolve");
    for face in inner_faces {
        let Surface::Plane(plane) = after.model.geom().surface(after.model.face(face).surface) else {
            panic!("the edited inner wall must remain planar");
        };
        assert_eq!(dec.rational_value(&plane.origin.y).unwrap(), Some(q("7/4")));
    }
    let corner_faces = after.model.faces_with_origin(after.result, &corner);
    corner_faces.first().expect("the saved round-face reference must resolve");
    for face in corner_faces {
        let Surface::Cylinder(cylinder) =
            after.model.geom().surface(after.model.face(face).surface)
        else {
            panic!("the edited reentrant round must remain cylindrical");
        };
        assert_eq!(dec.rational_value(&cylinder.radius_sq).unwrap(), Some(q("25/16")));
        assert_eq!(dec.rational_value(&cylinder.origin.x).unwrap(), Some(q("3")));
        assert_eq!(dec.rational_value(&cylinder.origin.y).unwrap(), Some(q("3")));
    }
    let edges = after.model.edges_with_name(after.result, &edge);
    edges.first().expect("the saved cap-to-round edge reference must resolve");
    for edge in edges {
        let curve = after.model.edge(edge).curve;
        let catallax_intersect::CurveForm::Circle { center, radius_sq, .. } =
            catallax_intersect::curve_form(after.model.geom(), curve, &mut dec).unwrap()
        else {
            panic!("the saved cap-to-round edge must resolve to its edited circular arc");
        };
        assert_eq!(dec.rational_value(&radius_sq).unwrap(), Some(q("25/16")));
        assert_eq!(dec.rational_value(&center.x).unwrap(), Some(q("3")));
        assert_eq!(dec.rational_value(&center.y).unwrap(), Some(q("3")));
        assert_eq!(dec.rational_value(&center.z).unwrap(), Some(q("0")));
    }
}

#[test]
fn display_mesh_keeps_the_cavity_open_with_closed_oriented_edges() {
    let built = build(&document());
    let options = catallax_mesh::MeshOptions::default();
    let (mesh, work) = with_operation_budget(WORK, || {
        catallax_mesh::mesh_body(
            &built.model, built.result, &options, &mut Decider::default(),
        ).unwrap().value
    });
    eprintln!("SHELL_MESH work={work:?}");
    assert!(!work.exhausted, "{work:?}");
    assert!(mesh.unavailable_controls.is_empty());
    assert!(mesh.faces.iter().all(|face| face.max_chord_error <= options.chord_tolerance));
    assert_eq!(
        mesh.faces.iter().map(|face| face.face).collect::<BTreeSet<_>>(),
        built.model.faces_of_body(built.result).into_iter().collect(),
    );
    let mut uses = BTreeMap::<_, (usize, i32)>::new();
    let mut volume = 0.0;
    for triangle in &mesh.triangles {
        let [a, b, c] = triangle.vertices.map(|i| mesh.positions[i]);
        let cross = [
            b[1] * c[2] - b[2] * c[1],
            b[2] * c[0] - b[0] * c[2],
            b[0] * c[1] - b[1] * c[0],
        ];
        volume += (0..3).map(|i| a[i] * cross[i]).sum::<f64>() / 6.0;
        let center: [f64; 3] = std::array::from_fn(|i| (a[i] + b[i] + c[i]) / 3.0);
        assert!(!(center[0] > 3.01 && center[1] > 3.01),
            "display filled the missing source quadrant: {center:?}");
        let left_void = center[0] > 1.1 && center[0] < 1.9
            && center[1] > 1.1 && center[1] < 6.9;
        let bottom_void = center[1] > 1.1 && center[1] < 1.9
            && center[0] > 1.1 && center[0] < 6.9;
        assert!(!(center[2] > 1.1 && (left_void || bottom_void)),
            "display filled the open cavity: {center:?}");
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let count = uses.entry((a.min(b), a.max(b))).or_default();
            count.0 += 1;
            count.1 += if a < b { 1 } else { -1 };
        }
    }
    assert!(uses.values().all(|uses| *uses == (2, 0)),
        "display edges must have two opposite oriented uses");
    let expected = 120.0 + 0.75 * std::f64::consts::PI;
    assert!((volume - expected).abs() < expected * 0.02, "{volume}");
}
