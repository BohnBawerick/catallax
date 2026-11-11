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
    spline_shape(text, bow, 0)
}

fn spline_shape(text: String, bow: &str, axis: usize) -> Vec<u8> {
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
    mid[axis] = &mid[axis] + &common::r(bow);
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
    let bytes = spline_box("0.1");
    let result = import(&bytes);
    let Err(refusal @ ImportRefusal::IncidenceGap { .. }) = result else {
        panic!("expected measured cap refusal")
    };
    let repeated = catallax_import::measure_incidence(&bytes, &refusal).unwrap();
    let ImportRefusal::IncidenceGap { gap, .. } = refusal else {
        unreachable!()
    };
    assert_eq!(*gap, repeated);
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

#[test]
fn the_source_gap_cap_applies_when_both_vertex_and_curve_move() {
    let mut text = String::from_utf8(spline_box("0")).unwrap();
    let source = read_step(text.as_bytes()).unwrap();
    let edge = &source.solids[0].faces[0].loops[0].edges[0];
    let vertex_line = text
        .lines()
        .find(|line| line.starts_with(&format!("{}=", edge.start.name)))
        .unwrap();
    let vertex_point = vertex_line
        .split(',')
        .next_back()
        .unwrap()
        .trim_end_matches(");")
        .to_owned();
    for (id, x, y, z) in [
        (
            "#900001",
            "0.006",
            common::dec(&edge.start.point[1]),
            common::dec(&edge.start.point[2]),
        ),
        (
            "#900002",
            "0.004",
            common::dec(&((&edge.start.point[1] + &edge.end.point[1]).mul_pow2(-1))),
            common::dec(&edge.start.point[2]),
        ),
        (
            "#900003",
            "0.002",
            common::dec(&edge.end.point[1]),
            common::dec(&edge.end.point[2]),
        ),
        (
            vertex_point.as_str(),
            "-0.0040000088888849382751175941786105555773972114115377466928292691521294620970443500018202888",
            common::dec(&edge.start.point[1]),
            common::dec(&edge.start.point[2]),
        ),
    ] {
        let old = text
            .lines()
            .find(|line| line.starts_with(&format!("{id}=")))
            .unwrap();
        text = text.replace(old, &format!("{id}=CARTESIAN_POINT('',({x},{y},{z}));"));
    }
    let source = read_step(text.as_bytes()).unwrap();
    let edge = &source.solids[0].faces[0].loops[0].edges[0];
    let catallax_import::Curve::BSpline(s) = &edge.curve else {
        panic!("spline fixture")
    };
    let point = |p: &[catallax_numbers::Rational; 3]| {
        catallax_shapes::RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone())
    };
    let a = point(&s.geometry.controls()[0]);
    let b = point(&s.geometry.controls()[2]);
    assert_eq!(
        s.geometry.controls()[1],
        std::array::from_fn(
            |i| (&s.geometry.controls()[0][i] + &s.geometry.controls()[2][i]).mul_pow2(-1)
        )
    );
    let p = point(&edge.start.point);
    let d = b.minus(&a);
    let v = p.minus(&a);
    let length2 = d.norm_squared();
    let parameter = v.dot(&d) / length2.clone();
    assert!(parameter > common::r("0") && parameter < common::r("1"));
    let distance2 = v.cross(&d).norm_squared() / length2;
    assert!(
        distance2 > common::r("0.0001"),
        "the exact source gap exceeds 0.01 mm"
    );
    let bytes = text.into_bytes();
    let result = import(&bytes);
    assert!(
        matches!(result, Err(ImportRefusal::IncidenceGap { .. })),
        "{result:?}"
    );
}

#[test]
fn repaired_line_carriers_preserve_reversed_step_sense() {
    let mut text = String::from_utf8(spline_box("0.001")).unwrap();
    let source = read_step(text.as_bytes()).unwrap();
    let edge = &source.solids[0].faces[0].loops[0].edges[1];
    let catallax_import::Curve::Line(_, direction) = &edge.curve else {
        panic!("line fixture")
    };
    let edge_line = text
        .lines()
        .find(|l| l.starts_with(&format!("{}=", edge.edge)))
        .unwrap()
        .to_owned();
    let curve_id = edge_line.split(',').nth(3).unwrap();
    let curve_line = text
        .lines()
        .find(|l| l.starts_with(&format!("{curve_id}=")))
        .unwrap();
    let vector_id = curve_line.split(',').nth(2).unwrap().trim_end_matches(");");
    let vector_line = text
        .lines()
        .find(|l| l.starts_with(&format!("{vector_id}=")))
        .unwrap();
    let direction_id = vector_line.split(',').nth(1).unwrap().to_owned();
    let direction_line = text
        .lines()
        .find(|l| l.starts_with(&format!("{direction_id}=")))
        .unwrap()
        .to_owned();
    let reversed = direction.clone().map(|x| -x);
    text = text.replace(
        &direction_line,
        &format!(
            "{direction_id}=DIRECTION('',({},{},{}));",
            common::dec(&reversed[0]),
            common::dec(&reversed[1]),
            common::dec(&reversed[2])
        ),
    );
    text = text.replace(&edge_line, &edge_line.replace(",.T.);", ",.F.);"));
    let (model, bodies) = import_marked(text.as_bytes(), "reversed-line-import");
    let point = |vertex| match model.geom().point(model.vertex(vertex).point) {
        catallax_shapes::PointDef::Explicit(p) => p.clone(),
        _ => panic!("explicit imported vertex"),
    };
    let mut reversed_count = 0;
    for id in model.edges_of_body(bodies[0]) {
        let edge = model.edge(id);
        let catallax_shapes::Curve::Line(line) = model.geom().curve(edge.curve) else {
            continue;
        };
        let direction = point(edge.end.unwrap()).minus(&point(edge.start.unwrap()));
        let dot = line.dir.dot(&direction).as_rational().unwrap();
        if edge.sense == catallax_topology::Sense::Same {
            assert!(dot > common::r("0"));
        } else {
            reversed_count += 1;
            assert!(
                dot < common::r("0"),
                "a reversed carrier must run against its endpoint order"
            );
        }
    }
    assert!(reversed_count > 0);
}

#[test]
fn bounded_planar_repair_keeps_a_concave_prism_closed() {
    let base = [
        ["0", "0"],
        ["3", "0"],
        ["3", "1"],
        ["1", "1"],
        ["1", "3"],
        ["0", "3"],
    ];
    let mut faces = vec![
        base.iter().rev().map(|p| [p[0], p[1], "0"]).collect(),
        base.iter().map(|p| [p[0], p[1], "2"]).collect(),
    ];
    for i in 0..base.len() {
        let a = base[i];
        let b = base[(i + 1) % base.len()];
        faces.push(vec![
            [a[0], a[1], "0"],
            [b[0], b[1], "0"],
            [b[0], b[1], "2"],
            [a[0], a[1], "2"],
        ]);
    }
    let bytes = spline_shape(common::polyhedra_step(&[faces]), "0.001", 2);
    let (model, bodies) = import_marked(&bytes, "concave-spline-import");
    assert!(model.check(bodies[0]).is_empty());
    let volume = catallax_measure::volume(
        &model,
        bodies[0],
        &catallax_measure::Accuracy::default(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap()
    .value;
    assert_eq!(volume.lo, common::r("10"));
    assert_eq!(volume.hi, common::r("10"));
}

#[test]
fn a_spline_file_repairs_analytic_sibling_solids_independently() {
    let a = common::BoxFixture::new(["2", "3", "4"]);
    let mut b = common::BoxFixture::new(["2", "3", "4"]);
    b.origin = [common::r("5"), common::r("0"), common::r("0")];
    let mut text =
        String::from_utf8(spline_shape(common::boxes_step(&[a, b]), "0.001", 0)).unwrap();
    let source = read_step(text.as_bytes()).unwrap();
    let vertex = &source.solids[1].faces[0].loops[0].edges[0].start;
    let line = text
        .lines()
        .find(|line| line.starts_with(&format!("{}=", vertex.name)))
        .unwrap();
    let id = line
        .split(',')
        .next_back()
        .unwrap()
        .trim_end_matches(");")
        .to_owned();
    let old = text
        .lines()
        .find(|line| line.starts_with(&format!("{id}=")))
        .unwrap();
    let mut p = vertex.point.clone();
    p[0] = &p[0] + &common::r("0.000001");
    text = text.replace(
        old,
        &format!(
            "{id}=CARTESIAN_POINT('',({},{},{}));",
            common::dec(&p[0]),
            common::dec(&p[1]),
            common::dec(&p[2])
        ),
    );
    let (model, bodies) = import_marked(text.as_bytes(), "sibling-spline-import");
    assert_eq!(bodies.len(), 2);
    for body in bodies {
        assert!(model.check(body).is_empty());
        assert!(!model.body(body).marks.is_empty());
    }
}

#[test]
fn a_straight_repair_can_have_a_zero_source_parameter_derivative_at_an_endpoint() {
    let mut text = String::from_utf8(spline_box("0")).unwrap();
    let source = read_step(text.as_bytes()).unwrap();
    let mut p = source.solids[0].faces[0].loops[0].edges[0]
        .start
        .point
        .clone();
    p[0] = &p[0] + &common::r("0.001");
    for id in ["#900001", "#900002"] {
        let old = text
            .lines()
            .find(|line| line.starts_with(&format!("{id}=")))
            .unwrap();
        text = text.replace(
            old,
            &format!(
                "{id}=CARTESIAN_POINT('',({},{},{}));",
                common::dec(&p[0]),
                common::dec(&p[1]),
                common::dec(&p[2])
            ),
        );
    }
    let (model, bodies) = import_marked(text.as_bytes(), "flat-parameter-import");
    assert!(model.check(bodies[0]).is_empty());
    assert_eq!(
        model
            .body(bodies[0])
            .marks
            .worst(catallax_numbers::BoundKind::Distance),
        Some(&common::r("0.001"))
    );
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

#[test]
fn bounded_planar_repair_preserves_a_through_hole() {
    use catallax_numbers::Decider;
    use catallax_shapes::{Pose, RVector3, Rotation};
    let mut source = Model::new();
    let outer = catallax_topology::block(
        &mut source,
        [4.into(), 4.into(), 2.into()],
        &Pose::identity(),
    )
    .unwrap();
    let pose = Pose::from_rational(
        Rotation::identity(),
        RVector3::new(1.into(), 1.into(), (-1).into()),
    );
    let inner =
        catallax_topology::block(&mut source, [2.into(), 2.into(), 4.into()], &pose).unwrap();
    let mut dec = Decider::default();
    let body = catallax_boolean::boolean(
        &mut source,
        outer,
        inner,
        catallax_boolean::BoolOp::Subtract,
        &mut dec,
    )
    .unwrap()
    .body;
    let mut bytes = Vec::new();
    catallax_step::write_step(
        &source,
        &[body],
        &catallax_step::StepOptions::default(),
        &mut dec,
        &mut bytes,
    )
    .unwrap();
    let file = read_step(&bytes).unwrap();
    assert!(file.solids[0].faces.iter().any(|f| f.loops.len() > 1));
    let catallax_import::Surface::Plane(plane) = &file.solids[0].faces[0].surface else {
        panic!("planar fixture")
    };
    let axis = (0..3).max_by_key(|&i| plane.axis[i].abs()).unwrap();
    let bytes = spline_shape(String::from_utf8(bytes).unwrap(), "0.001", axis);
    let (model, bodies) = import_marked(&bytes, "hole-spline-import");
    assert!(model.check(bodies[0]).is_empty());
    let volume = catallax_measure::volume(
        &model,
        bodies[0],
        &catallax_measure::Accuracy::default(),
        &mut dec,
    )
    .unwrap()
    .value;
    assert_eq!(volume.lo, common::r("24"));
    assert_eq!(volume.hi, common::r("24"));
}
