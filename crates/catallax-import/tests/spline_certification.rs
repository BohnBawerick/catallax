//! Synthetic STEP inputs exercise the production spline solid entry.

mod common;

use catallax_import::{ImportRefusal, import_step_bodies, read_step, sha256_hex};
use catallax_numbers::with_operation_budget;
use catallax_recipe::StepId;
use catallax_topology::{BodyId, Model};

fn spline_box(bow: &str) -> Vec<u8> {
    spline_box_dimensions(bow, ["2", "3", "4"])
}

fn spline_box_dimensions(bow: &str, dimensions: [&str; 3]) -> Vec<u8> {
    let text = common::box_step(common::BoxFixture::new(dimensions));
    let file = read_step(text.as_bytes()).unwrap();
    let edge = &file.solids[0].faces[0].loops[0].edges[0];
    let edge_line = text
        .lines()
        .find(|line| line.starts_with(&format!("{}=", edge.edge)))
        .unwrap();
    let curve_id = edge_line.split(',').nth(3).unwrap();
    let line = text
        .lines()
        .find(|line| line.starts_with(&format!("{curve_id}=")))
        .unwrap();
    let a = &edge.start.point;
    let b = &edge.end.point;
    let mut mid = std::array::from_fn::<_, 3, _>(|i| (&a[i] + &b[i]).mul_pow2(-1));
    mid[0] = &mid[0] + &common::r(bow);
    let mut points = String::new();
    for (id, point) in [a, &mid, b].into_iter().enumerate() {
        points.push_str(&format!(
            "#{}=CARTESIAN_POINT('',({},{},{}));\n",
            900001 + id,
            common::dec(&point[0]),
            common::dec(&point[1]),
            common::dec(&point[2])
        ));
    }
    text.replace(line, &format!("{curve_id}=B_SPLINE_CURVE_WITH_KNOTS('',2,(#900001,#900002,#900003),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,1.),.UNSPECIFIED.);"))
        .replace("DATA;", &format!("DATA;\n{points}")).into_bytes()
}

fn import(bytes: &[u8]) -> Result<usize, ImportRefusal> {
    let (result, work) = with_operation_budget(20_000_000, || {
        import_step_bodies(
            &mut Model::new(),
            bytes,
            &sha256_hex(bytes),
            &StepId::new("spline-import").unwrap(),
        )
        .map(|bodies| bodies.len())
    });
    assert!(!work.exhausted);
    println!(
        "spline-import {}",
        serde_json::json!({
            "input_sha256": sha256_hex(bytes), "result": &result, "work_used": work.used,
        })
    );
    result
}

fn import_marked(bytes: &[u8], step: &str) -> (Model, Vec<BodyId>) {
    let mut model = Model::new();
    let (result, work) = with_operation_budget(20_000_000, || {
        import_step_bodies(
            &mut model,
            bytes,
            &sha256_hex(bytes),
            &StepId::new(step).unwrap(),
        )
    });
    assert!(!work.exhausted);
    let bodies = result.unwrap();
    let marks=bodies.iter().flat_map(|b|model.body(*b).marks.iter()).map(|m|serde_json::json!({"code":m.code,"bound":m.bound.to_string(),"definition":m.detail})).collect::<Vec<_>>();
    println!(
        "marked-spline-import {}",
        serde_json::json!({"source_sha256":sha256_hex(bytes),"work_used":work.used,"marks":marks})
    );
    (model, bodies)
}

#[test]
fn admits_an_exact_straight_spline_edge_through_the_step_entry() {
    assert_eq!(import(&spline_box("0")).unwrap(), 1);
}

#[test]
fn records_the_measured_deviation_of_a_repaired_planar_spline_edge() {
    let bytes = spline_box("0.001");
    let (model, bodies) = import_marked(&bytes, "bounded-spline-import");
    assert!(model.check(bodies[0]).is_empty());
    let marks = &model.body(bodies[0]).marks;
    assert!(!marks.is_empty());
    assert!(marks.iter().any(|mark|mark.code=="import.incidence-repair" && mark.bound==common::r("0.0005")));
}

#[test]
fn refuses_a_spline_edge_above_the_import_cap_with_its_measured_size() {
    let Err(ImportRefusal::IncidenceGap { gap, .. }) = import(&spline_box("0.1")) else {
        panic!("expected measured cap refusal")
    };
    assert!(common::r(&gap.distance_lower_mm) > common::r("0.01"));
    assert!(common::r(&gap.distance_lower_mm) <= common::r("0.05"));
    assert!(common::r(&gap.distance_upper_mm) >= common::r("0.05"));
}

#[test]
fn accepts_the_exact_cap_and_refuses_the_next_larger_gap() {
    let bytes = spline_box("0.02");
    let (model, bodies) = import_marked(&bytes, "cap-spline-import");
    assert_eq!(
        model
            .body(bodies[0])
            .marks
            .worst(catallax_numbers::BoundKind::Distance),
        Some(&common::r("0.01"))
    );
    assert!(matches!(
        import(&spline_box("0.020000002")),
        Err(ImportRefusal::IncidenceGap { .. })
    ));
}

#[test]
fn a_small_gap_cannot_consume_a_neighbouring_feature() {
    let result = import(&spline_box_dimensions("0.001", ["0.001", "3", "4"]));
    let Err(ImportRefusal::IncidenceGap { reason, .. }) = result else {
        panic!("expected feature refusal")
    };
    assert!(reason.contains("feature"));
}

#[test]
fn a_vertex_only_repair_records_the_affected_source_entity() {
    let text = String::from_utf8(spline_box("0")).unwrap();
    let source = read_step(text.as_bytes()).unwrap();
    let vertex = &source.solids[0].faces[0].loops[0].edges[0].start;
    let line = text
        .lines()
        .find(|line| line.starts_with(&format!("{}=", vertex.name)))
        .unwrap();
    let point_id = line.split(',').next_back().unwrap().trim_end_matches(");");
    let original = text
        .lines()
        .find(|line| line.starts_with(&format!("{point_id}=")))
        .unwrap();
    let mut position = vertex.point.clone();
    position[0] = &position[0] + &common::r("0.000001");
    let bytes = text
        .replace(
            original,
            &format!(
                "{point_id}=CARTESIAN_POINT('',({},{},{}));",
                common::dec(&position[0]),
                common::dec(&position[1]),
                common::dec(&position[2])
            ),
        )
        .into_bytes();
    let (model, bodies) = import_marked(&bytes, "vertex-spline-import");
    assert!(model.check(bodies[0]).is_empty());
    let mark = model
        .body(bodies[0])
        .marks
        .iter()
        .find(|m| m.code == "import.incidence-vertex-repair")
        .unwrap();
    assert_eq!(mark.bound, common::r("0.000001"));
    let record: serde_json::Value = serde_json::from_str(&mark.detail).unwrap();
    assert_eq!(record["entity"], vertex.name);
    assert_eq!(record["entity_kind"], "vertex");
    assert_eq!(record["approved_tolerance_mm"], "0.01");
}

fn spline_face_box() -> Vec<u8> {
    let text = common::box_step(common::BoxFixture::new(["2", "3", "4"]));
    let file = read_step(text.as_bytes()).unwrap();
    let face = &file.solids[0].faces[0];
    let bound = &face.loops[0];
    let mut corners: Vec<_> = bound
        .edges
        .iter()
        .map(|edge| {
            if edge.forward {
                edge.start.point.clone()
            } else {
                edge.end.point.clone()
            }
        })
        .collect();
    if !face.same_sense {
        corners.swap(1, 3);
    }
    let mut points = String::new();
    for (id, point) in [0, 3, 1, 2].map(|i| &corners[i]).into_iter().enumerate() {
        points.push_str(&format!(
            "#{}=CARTESIAN_POINT('',({},{},{}));\n",
            900001 + id,
            common::dec(&point[0]),
            common::dec(&point[1]),
            common::dec(&point[2])
        ));
    }
    let face_line = text
        .lines()
        .find(|line| line.starts_with(&format!("{}=", face.name)))
        .unwrap();
    let surface_id = face_line.split(',').nth(2).unwrap();
    let line = text
        .lines()
        .find(|line| line.starts_with(&format!("{surface_id}=")))
        .unwrap();
    text.replace(line, &format!("{surface_id}=B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#900001,#900002),(#900003,#900004)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.);"))
        .replace("DATA;", &format!("DATA;\n{points}")).into_bytes()
}

#[test]
fn admits_a_complete_embedded_patch_face_through_the_step_entry() {
    assert_eq!(import(&spline_face_box()).unwrap(), 1);
}
