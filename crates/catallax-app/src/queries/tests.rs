use super::*;
use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{BudgetPolicy, MarkLog, with_operation_budget};
use catallax_topology::{Body, BodyKind, sphere};

fn vector(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(x.into(), y.into(), z.into())
}
fn point(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(x.into(), y.into(), z.into())
}
fn cube(model: &mut Model, xyz: [i64; 3], size: i64) -> BodyId {
    block(
        model,
        [size.into(), size.into(), size.into()],
        &Pose::from_rational(Rotation::identity(), vector(xyz[0], xyz[1], xyz[2])),
    )
    .unwrap()
}
fn plane(origin: RPoint3, normal: RVector3) -> SectionPlane {
    SectionPlane { origin, normal }
}
fn exact(value: &Real, expected: i64) {
    assert_eq!(
        sign(&(value - &Real::from(expected)), &mut Decider::default()).unwrap(),
        0
    );
}

#[test]
fn sections_have_exact_trims_for_axis_oblique_and_reversed_planes() {
    let mut model = Model::new();
    let body = cube(&mut model, [0, 0, 0], 2);
    let counts = model.counts(body);
    for (request, count) in [
        (plane(point(0, 0, 1), vector(0, 0, 1)), 4),
        (plane(point(0, 0, 1), vector(0, 0, -7)), 4),
        (plane(point(1, 1, 1), vector(1, 1, 1)), 6),
        (plane(point(0, 0, 3), vector(0, 0, 1)), 0),
    ] {
        let answer = section_curves(
            &model,
            body,
            &request,
            &Accuracy::default(),
            &mut Decider::default(),
        )
        .unwrap();
        assert!(answer.marks.is_empty());
        assert_eq!(answer.value.curves.len(), count);
        assert!(answer.value.points.is_empty());
        for curve in &answer.value.curves {
            assert!(!curve.faces_of_first.is_empty());
            assert!(matches!(
                curve_form(
                    answer.value.model.geom(),
                    curve.curve,
                    &mut Decider::default()
                )
                .unwrap(),
                CurveForm::Line { .. }
            ));
            for endpoint in [curve.start.unwrap(), curve.end.unwrap()] {
                let position = answer.value.model.geom().position(endpoint);
                exact(
                    &position
                        .minus(&request.origin.to_real())
                        .dot(&request.normal.to_real()),
                    0,
                );
            }
        }
    }
    assert_eq!(model.counts(body), counts);
}

#[test]
fn sphere_section_and_silhouette_are_exact_circles() {
    let mut model = Model::new();
    let body = sphere(&mut model, 5.into(), &Pose::identity()).unwrap();
    let answer = section_curves(
        &model,
        body,
        &plane(point(0, 0, 3), vector(0, 0, 1)),
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(answer.value.curves.len(), 1);
    let curve = &answer.value.curves[0];
    assert_eq!((curve.start, curve.end), (None, None));
    let CurveForm::Circle {
        center, radius_sq, ..
    } = curve_form(
        answer.value.model.geom(),
        curve.curve,
        &mut Decider::default(),
    )
    .unwrap()
    else {
        panic!("circle")
    };
    exact(&center.z, 3);
    exact(&radius_sq, 16);
    let outline = silhouette(&model, body, &vector(1, 2, 3), &mut Decider::default()).unwrap();
    assert!(outline.value.edges.is_empty());
    assert_eq!(outline.value.circles.len(), 1);
    let (_, circle) = &outline.value.circles[0];
    exact(&circle.radius_sq, 25);
    exact(&circle.normal.dot(&circle.ref_dir), 0);
}

#[test]
fn convex_outline_omits_hidden_edges_and_keeps_edge_on_views() {
    let mut model = Model::new();
    let body = cube(&mut model, [0, 0, 0], 2);
    for (direction, count) in [
        (vector(0, 0, 1), 4),
        (vector(1, 1, 1), 6),
        (vector(0, 0, -1), 4),
    ] {
        let outline = silhouette(&model, body, &direction, &mut Decider::default()).unwrap();
        assert_eq!(outline.value.edges.len(), count);
        assert!(outline.value.circles.is_empty());
        if direction.x.is_zero() && direction.y.is_zero() {
            let expected = if direction.z.sign() == catallax_numbers::Sign::Positive {
                2
            } else {
                0
            };
            for edge in &outline.value.edges {
                let crate::curves::CurveTrim::Finite { start, end } = &edge.trim else {
                    panic!()
                };
                exact(&start.position.z, expected);
                exact(&end.position.z, expected);
            }
        }
    }
}

#[test]
fn overlap_empty_touching_contained_and_crossing_have_proved_volumes() {
    let mut model = Model::new();
    let a = cube(&mut model, [0, 0, 0], 2);
    let counts = model.counts(a);
    for (xyz, size, expected) in [
        ([3, 0, 0], 2, 0),
        ([2, 0, 0], 2, 0),
        ([0, 0, 0], 1, 1),
        ([1, 1, 1], 2, 1),
        ([0, 0, 0], 2, 8),
    ] {
        let b = cube(&mut model, xyz, size);
        let before = model.counts(a);
        let volume =
            overlap_volume(&model, a, b, &Accuracy::default(), &mut Decider::default()).unwrap();
        assert!(volume.value.contains(&Rational::from(expected)));
        if expected == 0 {
            assert!(volume.value.is_exact());
        } else {
            assert!(
                volume.value.width()
                    <= Rational::parse("1e-12").unwrap() * Rational::from(expected)
            );
        }
        assert!(volume.marks.is_empty());
        assert_eq!(model.counts(a), before);
    }
    assert_eq!(model.counts(a), counts);
    let volume =
        overlap_volume(&model, a, a, &Accuracy::default(), &mut Decider::default()).unwrap();
    assert!(volume.value.contains(&Rational::from(8)));
}

#[test]
fn contact_patches_are_clipped_and_edge_touches_have_no_area_contact() {
    let mut model = Model::new();
    let a = cube(&mut model, [0, 0, 0], 2);
    for (xyz, size, count, area) in [
        ([2, 0, 0], 1, 1, 1),
        ([2, 2, 0], 2, 0, 0),
        ([3, 0, 0], 1, 0, 0),
    ] {
        let b = cube(&mut model, xyz, size);
        let contacts = contact_faces(&model, a, b, &mut Decider::default()).unwrap();
        assert_eq!(contacts.value.first.len(), count);
        assert_eq!(contacts.value.second.len(), count);
        for patch in contacts.value.first.iter().chain(&contacts.value.second) {
            assert_eq!(patch.orientation, FaceClass::OnOpposite);
            let measured = catallax_measure::face_area(
                &contacts.value.model,
                patch.face,
                &Accuracy::default(),
                &mut Decider::default(),
            )
            .unwrap();
            assert_eq!(measured.value.lo, Rational::from(area));
            assert_eq!(measured.value.hi, Rational::from(area));
        }
    }
    let same = contact_faces(&model, a, a, &mut Decider::default()).unwrap();
    assert_eq!(same.value.first.len(), 6);
    assert!(
        same.value
            .first
            .iter()
            .all(|p| p.orientation == FaceClass::OnSame)
    );
}

#[test]
fn invalid_coplanar_and_starved_requests_keep_typed_evidence() {
    let mut model = Model::new();
    let body = cube(&mut model, [0, 0, 0], 2);
    let coplanar = section_curves(
        &model,
        body,
        &plane(point(0, 0, 0), vector(0, 0, 1)),
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(coplanar.refusal.code, "app.section-coplanar");
    assert_eq!(
        (
            coplanar.operands[0].faces,
            coplanar.operands[0].edges,
            coplanar.operands[0].vertices
        ),
        (6, 12, 8)
    );
    assert_eq!(
        silhouette(&model, body, &vector(0, 0, 0), &mut Decider::default())
            .unwrap_err()
            .refusal
            .code,
        "app.silhouette-direction"
    );
    assert_eq!(
        section_curves(
            &model,
            body,
            &plane(point(0, 0, 1), vector(0, 0, 0)),
            &Accuracy::default(),
            &mut Decider::default()
        )
        .unwrap_err()
        .refusal
        .code,
        "app.section-normal"
    );
    let mut starved = Decider::new(BudgetPolicy {
        max_work: 0,
        max_precision_bits: 64,
    });
    let failed = silhouette(&model, body, &vector(1, 2, 3), &mut starved).unwrap_err();
    assert!(failed.refusal.undecided().is_some(), "{failed:?}");
    let (failed, work) = with_operation_budget(0, || {
        silhouette(&model, body, &vector(1, 2, 3), &mut Decider::default())
    });
    assert!(failed.is_err());
    assert!(work.exhausted);
}

#[test]
fn empty_bodies_have_empty_curves_zero_overlap_and_no_contacts() {
    let mut model = Model::new();
    let body = model.add_body(Body {
        kind: BodyKind::Empty,
        lumps: vec![],
        marks: MarkLog::new(),
    });
    let solid = cube(&mut model, [0, 0, 0], 2);
    assert!(
        silhouette(&model, body, &vector(0, 0, 1), &mut Decider::default())
            .unwrap()
            .value
            .edges
            .is_empty()
    );
    assert!(
        section_curves(
            &model,
            body,
            &plane(point(0, 0, 0), vector(0, 0, 1)),
            &Accuracy::default(),
            &mut Decider::default()
        )
        .unwrap()
        .value
        .curves
        .is_empty()
    );
    assert_eq!(
        overlap_volume(
            &model,
            solid,
            body,
            &Accuracy::default(),
            &mut Decider::default()
        )
        .unwrap()
        .value
        .lo,
        Rational::zero()
    );
    let contact = contact_faces(&model, solid, body, &mut Decider::default()).unwrap();
    assert!(contact.value.first.is_empty() && contact.value.second.is_empty());
}

#[test]
fn nonconvex_outline_refuses_with_operand_counts() {
    let mut model = Model::new();
    let a = cube(&mut model, [0, 0, 0], 2);
    let b = cube(&mut model, [1, 1, 0], 2);
    let body = boolean(&mut model, a, b, BoolOp::Unite, &mut Decider::default())
        .unwrap()
        .body;
    let failure = silhouette(&model, body, &vector(1, 2, 3), &mut Decider::default()).unwrap_err();
    assert_eq!(failure.refusal.code, "app.silhouette-nonconvex");
    assert!(failure.operands[0].faces > 6);
}

#[test]
fn tangent_section_keeps_the_isolated_contact_point() {
    let mut model = Model::new();
    let body = sphere(&mut model, 5.into(), &Pose::identity()).unwrap();
    let answer = section_curves(
        &model,
        body,
        &plane(point(0, 0, 5), vector(0, 0, 1)),
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(answer.value.curves.is_empty());
    assert_eq!(answer.value.points.len(), 1);
    let position = answer.value.model.geom().position(answer.value.points[0]);
    exact(&position.x, 0);
    exact(&position.y, 0);
    exact(&position.z, 5);
}

#[test]
fn construction_marks_survive_each_query() {
    use catallax_numbers::{BoundKind, Mark};
    let mut model = Model::new();
    let source = cube(&mut model, [0, 0, 0], 2);
    let mut marked = model.body(source).clone();
    let mark = Mark::new(
        "test.source-bound",
        BoundKind::Distance,
        Rational::parse("1e-9").unwrap(),
        "retained construction bound",
    );
    marked.marks.push(mark.clone());
    let body = model.add_body(marked);
    let accuracy = Accuracy::default();
    let mut dec = Decider::default();
    let section = section_curves(
        &model,
        body,
        &plane(point(0, 0, 1), vector(0, 0, 1)),
        &accuracy,
        &mut dec,
    )
    .unwrap();
    let outline = silhouette(&model, body, &vector(0, 0, 1), &mut dec).unwrap();
    let overlap = overlap_volume(&model, body, body, &accuracy, &mut dec).unwrap();
    let contacts = contact_faces(&model, body, body, &mut dec).unwrap();
    for marks in [
        &section.marks,
        &outline.marks,
        &overlap.marks,
        &contacts.marks,
    ] {
        assert!(marks.iter().any(|m| m == &mark));
    }
}
