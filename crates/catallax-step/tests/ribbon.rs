//! A planar step narrower than a reader's shortest edge leaves the file, not the model.
mod common;
use catallax_boolean::{BoolOp, boolean};
use common::*;

// A 12 by 9 by 1 plate whose left wall stands `width` proud of a 1 by 9 by 6
// flange's left wall, multiplied through by `scale`.
fn flush_step(width: &str, scale: &str) -> (Model, BodyId) {
    let k = r(scale);
    let w = &r(width) * &k;
    let mut model = Model::new();
    let plate = Primitive::Block([&r("12") * &k + w.clone(), &r("9") * &k, k.clone()]);
    let flange = Primitive::Block([k.clone(), &r("9") * &k, &r("6") * &k]);
    let plate = plate.build(&mut model, &translated(&(-w).to_string(), "0", "0"));
    let flange = flange.build(&mut model, &Pose::identity());
    let body = boolean(
        &mut model,
        plate,
        flange,
        BoolOp::Unite,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    (model, body)
}

fn written(model: &Model, body: BodyId) -> (Step, StepReport) {
    let mut bytes = Vec::new();
    let answer = write_step(
        model,
        &[body],
        &StepOptions::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    (Step::parse(&bytes, "catallax"), answer.value)
}

fn assert_reconciled(model: &Model, body: BodyId, report: &StepReport) {
    let exact = model.counts(body);
    let faces: usize = report.regularizations.iter().map(|x| x.removed_faces).sum();
    let edges: usize = report.regularizations.iter().map(|x| x.removed_edges).sum();
    assert_eq!(exact.faces, report.faces + faces);
    assert_eq!(exact.edges + report.added_edges, report.edges + edges);
}

#[test]
fn a_step_too_narrow_to_read_is_written_as_one_wall() {
    for width in ["1e-30", "2.8e-15", "9.9e-8", "1e-7"] {
        let (model, body) = flush_step(width, "1");
        let exact = model.counts(body);
        assert_eq!(exact.faces, 10, "{width}");
        let (file, report) = written(&model, body);
        assert_eq!(report.faces, 8, "{width}");
        assert_eq!(file.count("ADVANCED_FACE"), 8);
        assert_eq!(file.count("EDGE_CURVE"), report.edges);
        let [record] = report.regularizations.as_slice() else {
            panic!("{width}: {:?}", report.regularizations);
        };
        assert_eq!((record.removed_faces, record.removed_edges), (2, 4));
        assert_eq!(record.max_deviation, r(width));
        assert!(report.max_deviation >= record.max_deviation);
        assert_reconciled(&model, body, &report);
        // The exact model is untouched.
        assert_eq!(model.counts(body), exact);
        assert!(model.check(body).is_empty());
    }
}

#[test]
fn a_step_rounded_to_nothing_is_repaired_too() {
    // At x = 50 the written walls coincide, so the end edges have no length.
    let mut model = Model::new();
    let plate = Primitive::Block([&r("12") + &r("1e-30"), n(9), n(1)]);
    let plate = plate.build(
        &mut model,
        &translated("49.999999999999999999999999999999", "0", "0"),
    );
    let flange = Primitive::Block([n(1), n(9), n(6)]);
    let flange = flange.build(&mut model, &translated("50", "0", "0"));
    let body = boolean(
        &mut model,
        plate,
        flange,
        BoolOp::Unite,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let (_, report) = written(&model, body);
    assert_eq!(report.faces, 8);
    assert_eq!(report.regularizations[0].max_deviation, r("1e-30"));
    assert_reconciled(&model, body, &report);
}

#[test]
fn a_shell_with_a_vertex_loop_keeps_its_topology() {
    let (mut model, body) = flush_step("1e-8", "1");
    let cone = Primitive::Cone(r("1/2"), n(0), n(2));
    let cone = cone.build(&mut model, &translated("6", "4", "1"));
    let body = boolean(
        &mut model,
        body,
        cone,
        BoolOp::Unite,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let (_, report) = written(&model, body);
    assert_eq!(report.faces, model.counts(body).faces);
    assert!(report.regularizations.is_empty());
}

#[test]
fn a_step_above_the_floor_is_not_repaired() {
    let (model, body) = flush_step("1.00000000001e-7", "1");
    let (_, report) = written(&model, body);
    assert_eq!(report.faces, 10);
    assert!(report.regularizations.is_empty());
}

#[test]
fn a_uniformly_tiny_body_has_no_short_ended_ribbon() {
    // Scaling the whole part scales the width and the length together, so the
    // step stays wide against its length or every edge falls below the floor.
    for scale in ["1e-30", "1e-20", "1e-8"] {
        for width in ["1", "1/2", "1e-3"] {
            let (model, body) = flush_step(width, scale);
            let (_, report) = written(&model, body);
            assert!(
                report.regularizations.is_empty(),
                "{scale} {width}: {:?}",
                report.regularizations
            );
            assert_eq!(report.faces, 10);
        }
    }
}

#[test]
fn tiny_primitives_and_crossing_cylinders_are_left_alone() {
    for scale in ["1e-30", "1e-20"] {
        let k = r(scale);
        for primitive in primitives() {
            let mut model = Model::new();
            let body = primitive.scale(&k).build(&mut model, &Pose::identity());
            assert!(written(&model, body).1.regularizations.is_empty());
        }
        let mut model = Model::new();
        let a = Primitive::Cylinder(n(1), n(4)).scale(&k);
        let a = a.build(
            &mut model,
            &translated("0", "0", &(&r("-2") * &k).to_string()),
        );
        let shift = (&r("-2") * &k).to_string();
        let b = Primitive::Cylinder(r("1/2"), n(4)).scale(&k);
        let b = b.build(
            &mut model,
            &xpose(&shift, &(&r("1/4") * &k).to_string(), "0"),
        );
        let body = boolean(&mut model, a, b, BoolOp::Unite, &mut Decider::default())
            .unwrap()
            .body;
        assert!(written(&model, body).1.regularizations.is_empty());
    }
}

#[test]
fn a_width_over_a_tenth_of_the_length_stays() {
    // 9e-8 wide and 8.1e-7 long: under the floor, but not thin against its length.
    let (model, body) = flush_step("1", "9e-8");
    let (_, report) = written(&model, body);
    assert!(report.regularizations.is_empty());
}

#[test]
fn a_thin_plate_is_not_a_step() {
    // Its narrow faces lie between opposing faces, so no wall can be moved.
    let mut model = Model::new();
    let plate = Primitive::Block([r("1e-8"), r("9"), r("5")]).build(&mut model, &Pose::identity());
    let (_, report) = written(&model, plate);
    assert_eq!(report.faces, 6);
    assert!(report.regularizations.is_empty());
}

#[test]
fn a_wall_with_a_hole_stays_and_its_neighbour_moves() {
    let (mut model, body) = flush_step("1e-8", "1");
    let tool =
        Primitive::Cylinder(r("1/4"), r("1")).build(&mut model, &xpose("-1/2", "9/2", "1/2"));
    let cut = boolean(
        &mut model,
        body,
        tool,
        BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let exact = model.counts(cut);
    let (_, report) = written(&model, cut);
    assert_eq!(report.faces + 2, exact.faces);
    let [record] = report.regularizations.as_slice() else {
        panic!("{:?}", report.regularizations);
    };
    let [wall, target] = record.merged_faces;
    let loops = |face| model.face(face).loops.len();
    assert_eq!((loops(wall), loops(target)), (1, 2));
    assert_eq!(record.max_deviation, r("1e-8"));
    assert_reconciled(&model, cut, &report);
}

#[test]
fn the_smaller_wall_moves() {
    let (model, body) = flush_step("1e-8", "1");
    let (_, report) = written(&model, body);
    let [wall, target] = report.regularizations[0].merged_faces;
    let area = |face| {
        catallax_measure::face_area(
            &model,
            face,
            &catallax_measure::Accuracy::default(),
            &mut Decider::default(),
        )
        .unwrap()
        .value
    };
    assert!(area(wall).hi < area(target).lo);
}
