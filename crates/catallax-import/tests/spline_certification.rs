//! Synthetic STEP inputs exercise the production spline solid entry.

mod common;

use catallax_import::{ImportRefusal, import_step_bodies, read_step, sha256_hex};
use catallax_numbers::with_operation_budget;
use catallax_recipe::StepId;
use catallax_topology::Model;

fn spline_box(bow: &str) -> Vec<u8> {
    let text = common::box_step(common::BoxFixture::new(["2", "3", "4"]));
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

#[test]
fn admits_an_exact_straight_spline_edge_through_the_step_entry() {
    assert_eq!(import(&spline_box("0")).unwrap(), 1);
}

#[test]
fn refuses_a_spline_edge_that_leaves_its_face() {
    assert!(matches!(
        import(&spline_box("0.001")),
        Err(ImportRefusal::NotSolid { .. })
    ));
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
