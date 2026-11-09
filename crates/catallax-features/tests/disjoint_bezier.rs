//! Disconnected text-like profiles retain their counters and construction names.
use catallax_boolean::{BoolOp, Stage, boolean};
use catallax_features::{Profile, extrude};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{ElementId, Loop, Point2, Pose, Segment, StepId};
use catallax_topology::{BodyId, Model};
use std::collections::BTreeSet;

fn point(x: i64, y: i64) -> Point2 {
    Point2([x.into(), y.into()])
}

fn arch(x: i64, hole: bool) -> Loop {
    let (left, right, bottom, top, prefix) = if hole {
        (x + 1, x + 3, 1, 3, "inner")
    } else {
        (x, x + 4, 0, 4, "outer")
    };
    Loop(vec![
        Segment::Bezier {
            id: ElementId::new(format!("{prefix}-arch-{x}")).unwrap(),
            controls: vec![point(left, top), point(right, top)],
            to: point(right, bottom),
        },
        Segment::Line {
            id: ElementId::new(format!("{prefix}-base-{x}")).unwrap(),
            to: point(left, bottom),
        },
    ])
}

fn prism(model: &mut Model, x: i64, counter: bool, name: &str) -> BodyId {
    let mut loops = vec![arch(x, false)];
    if counter {
        loops.push(arch(x, true));
    }
    let profile = Profile::new(&Pose::identity(), &loops, &mut Decider::default()).unwrap();
    extrude(
        model,
        &profile,
        &StepId::new(name).unwrap(),
        &2.into(),
        None,
    )
    .unwrap()
}

fn check_volume(
    model: &Model,
    body: BodyId,
    expected: &str,
) -> catallax_measure::VolumeAreaCentroid {
    let measured = catallax_measure::volume_area_centroid(
        model,
        body,
        &catallax_measure::Accuracy {
            relative: Rational::zero(),
            absolute: Rational::parse("1e-9").unwrap(),
        },
        &mut Decider::default(),
    )
    .unwrap();
    assert!(measured.marks.is_empty());
    assert!(
        measured
            .value
            .volume
            .contains(&Rational::parse(expected).unwrap())
    );
    assert!(measured.value.volume.width() <= Rational::parse("1e-9").unwrap());
    measured.value
}

#[test]
fn separate_glyph_prisms_keep_counters_topology_and_names() {
    let mut model = Model::new();
    let a = prism(&mut model, 0, true, "o");
    let b = prism(&mut model, 6, true, "a");
    let c = prism(&mut model, 12, false, "e");
    let faces: BTreeSet<_> = [a, b, c]
        .into_iter()
        .flat_map(|body| model.faces_of_body(body))
        .collect();
    let origins: BTreeSet<_> = faces
        .iter()
        .map(|f| model.face_origins(*f).to_vec())
        .collect();
    // The accumulated outer glyphs straddle the middle glyph's entire box.
    let first = boolean(&mut model, a, c, BoolOp::Unite, &mut Decider::default()).unwrap();
    let answer = boolean(
        &mut model,
        first.body,
        b,
        BoolOp::Unite,
        &mut Decider::default(),
    )
    .unwrap();
    assert!(answer.marks.is_empty());
    assert!(model.check(answer.body).is_empty());
    assert_eq!(model.body(answer.body).lumps.len(), 3);
    let result_faces = model.faces_of_body(answer.body);
    assert_eq!(result_faces.len(), faces.len());
    assert_eq!(
        result_faces
            .iter()
            .map(|f| model.face_origins(*f).to_vec())
            .collect::<BTreeSet<_>>(),
        origins
    );
    let measured = check_volume(&model, answer.body, "48");
    assert!(measured.area.contains(&Rational::from_i64(144)));
    for (actual, expected) in measured
        .centroid
        .unwrap()
        .iter()
        .zip(["43/5", "17/14", "1"])
    {
        assert!(actual.contains(&Rational::parse(expected).unwrap()));
        assert!(actual.width() <= Rational::parse("1e-9").unwrap());
    }
    let mesh = catallax_mesh::mesh_body(
        &model,
        answer.body,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(mesh.marks.is_empty());
    assert_eq!(mesh.value.faces.len(), faces.len());
}

#[test]
fn separate_analytic_and_spline_prisms_work_in_both_orders() {
    for reverse in [false, true] {
        let mut model = Model::new();
        let spline = prism(&mut model, 0, true, "glyph");
        let rectangle = Loop(
            [(8, 0), (8, 2), (6, 2), (6, 0)]
                .into_iter()
                .enumerate()
                .map(|(i, (x, y))| Segment::Line {
                    id: ElementId::new(format!("line{i}")).unwrap(),
                    to: point(x, y),
                })
                .collect(),
        );
        let profile =
            Profile::new(&Pose::identity(), &[rectangle], &mut Decider::default()).unwrap();
        let analytic = extrude(
            &mut model,
            &profile,
            &StepId::new("analytic").unwrap(),
            &2.into(),
            None,
        )
        .unwrap();
        let (a, b) = if reverse {
            (analytic, spline)
        } else {
            (spline, analytic)
        };
        let answer = boolean(&mut model, a, b, BoolOp::Unite, &mut Decider::default()).unwrap();
        assert!(model.check(answer.body).is_empty());
        assert_eq!(model.body(answer.body).lumps.len(), 2);
        check_volume(&model, answer.body, "112/5");
    }
}

#[test]
fn meeting_or_overlapping_spline_prisms_keep_the_intersection_refusal() {
    for x in [3, 4] {
        let mut model = Model::new();
        let a = prism(&mut model, 0, false, "left");
        let b = prism(&mut model, x, false, "right");
        let error = boolean(&mut model, a, b, BoolOp::Unite, &mut Decider::default()).unwrap_err();
        assert_eq!(error.stage, Stage::Intersect);
        assert!(matches!(
            error.refusal.code,
            "intersect.unsupported-pair" | "intersect.bezier-plane-oblique"
        ));
    }
}

#[test]
fn near_touching_spline_prisms_never_use_a_guessed_gap() {
    use catallax_shapes::{Isometry, RVector3};
    for gap in ["1e-20", "1e-80", "-1e-20"] {
        let mut model = Model::new();
        let a = prism(&mut model, 0, false, "left");
        let right = prism(&mut model, 4, false, "right");
        let b = catallax_topology::transformed_body(
            &mut model,
            right,
            &Isometry::translation(RVector3::new(
                Rational::parse(gap).unwrap(),
                Rational::zero(),
                Rational::zero(),
            )),
            Clone::clone,
        );
        let result = boolean(&mut model, a, b, BoolOp::Unite, &mut Decider::default());
        match result {
            Ok(answer) => {
                assert!(!gap.starts_with('-'));
                assert!(model.check(answer.body).is_empty());
                assert_eq!(model.body(answer.body).lumps.len(), 2);
                check_volume(&model, answer.body, "192/5");
            }
            Err(error) => assert!(matches!(
                error.refusal.code,
                "intersect.unsupported-pair" | "intersect.bezier-plane-oblique"
            )),
        }
    }
}

#[test]
fn one_extrude_of_separate_closed_regions_preserves_each_counter() {
    let loops = [
        arch(0, false),
        arch(0, true),
        arch(6, false),
        arch(6, true),
        arch(12, false),
    ];
    let profile = Profile::new(&Pose::identity(), &loops, &mut Decider::default()).unwrap();
    let mut model = Model::new();
    let body = extrude(
        &mut model,
        &profile,
        &StepId::new("word").unwrap(),
        &2.into(),
        None,
    )
    .unwrap();
    assert!(model.check(body).is_empty());
    assert_eq!(model.body(body).lumps.len(), 3);
    check_volume(&model, body, "48");
}
