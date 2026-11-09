#![allow(missing_docs)]

use catallax_features::{Profile, extrude};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{ElementId, Loop, Point2, Pose, Segment, StepId};
use catallax_topology::Model;

fn p(x: i64, y: i64) -> Point2 {
    Point2([x.into(), y.into()])
}
fn id(s: &str) -> ElementId {
    ElementId::new(s).unwrap()
}
fn arch(reverse: bool) -> Loop {
    let (controls, end, start) = if reverse {
        (vec![p(4, 4), p(0, 4)], p(0, 0), p(4, 0))
    } else {
        (vec![p(0, 4), p(4, 4)], p(4, 0), p(0, 0))
    };
    Loop(vec![
        Segment::Bezier {
            id: id("arch"),
            controls,
            to: end,
        },
        Segment::Line {
            id: id("base"),
            to: start,
        },
    ])
}
fn validate(loops: &[Loop]) -> Result<Profile, catallax_numbers::Refusal> {
    Profile::new(&Pose::identity(), loops, &mut Decider::default())
}

#[test]
fn quadratic_and_cubic_closed_profiles_validate_in_both_orders() {
    for reverse in [false, true] {
        let mut lp = arch(reverse);
        assert!(validate(&[lp.clone()]).is_ok());
        let Segment::Bezier { controls, .. } = &mut lp.0[0] else {
            unreachable!()
        };
        *controls = vec![p(2, 4)];
        assert!(validate(&[lp]).is_ok());
    }
}

#[test]
fn nested_and_disjoint_line_loops_use_the_curved_boundary() {
    let rectangle = |x: i64, y: i64| {
        Loop(
            [(x + 1, y), (x + 1, y + 1), (x, y + 1), (x, y)]
                .into_iter()
                .enumerate()
                .map(|(i, (x, y))| Segment::Line {
                    id: id(&format!("line{i}")),
                    to: p(x, y),
                })
                .collect(),
        )
    };
    assert!(validate(&[arch(false), rectangle(1, 1), rectangle(7, 7)]).is_ok());
    assert!(validate(&[arch(true), rectangle(1, 1)]).is_ok());
}

#[test]
fn crossing_or_touching_loops_never_pass_from_sampled_chords() {
    for y in [2, 3] {
        let crossing = Loop(vec![
            Segment::Line {
                id: id("a"),
                to: p(3, y),
            },
            Segment::Line {
                id: id("b"),
                to: p(3, y + 1),
            },
            Segment::Line {
                id: id("c"),
                to: p(1, y + 1),
            },
            Segment::Line {
                id: id("d"),
                to: p(1, y),
            },
        ]);
        let error = validate(&[arch(false), crossing]).unwrap_err();
        assert!(matches!(
            error.code,
            "sketch.loop-intersection" | "sketch.bezier-proof"
        ));
    }
}

#[test]
fn invalid_or_unproved_spans_keep_named_refusals() {
    for (controls, expected) in [
        (vec![], "sketch.bezier-controls"),
        (vec![p(-2, 4), p(4, 4)], "shapes.bezier-extrusion-embedding"),
        (vec![p(0, 0), p(4, 4)], "shapes.bezier-extrusion-singular"),
    ] {
        let mut lp = arch(false);
        let Segment::Bezier {
            controls: original, ..
        } = &mut lp.0[0]
        else {
            unreachable!()
        };
        *original = controls;
        assert_eq!(validate(&[lp]).unwrap_err().code, expected);
    }
    let circle = Loop(vec![Segment::Circle {
        id: id("circle"),
        centre: p(9, 9),
        radius: 1.into(),
    }]);
    assert_eq!(
        validate(&[arch(false), circle]).unwrap_err().code,
        "sketch.bezier-curve-mix"
    );
}

#[test]
fn normal_extrusions_share_exact_spline_rims_and_keep_sweep_origins() {
    use catallax_shapes::{Curve, SplineBoundary, SplineSupport, Surface};
    use catallax_topology::Sense;
    for height in [-2, 2] {
        for reverse in [false, true] {
            let profile = validate(&[arch(reverse)]).unwrap();
            let mut model = Model::new();
            let body = extrude(
                &mut model,
                &profile,
                &StepId::new("sweep").unwrap(),
                &Rational::from_i64(height),
                None,
            )
            .unwrap();
            assert!(
                catallax_topology::BezierExtrusionBand::in_body(
                    &model,
                    body,
                    &mut Decider::default()
                )
                .unwrap()
            );
            assert!(model.check(body).is_empty(), "{:?}", model.check(body));
            assert_eq!(model.faces_of_body(body).len(), 4);
            assert_eq!(model.edges_of_body(body).len(), 6);
            for face in model.faces_of_body(body) {
                assert!(!model.face_origins(face).is_empty());
                if let Surface::BSpline(SplineSupport::BezierExtrusion(span)) =
                    model.geom().surface(model.face(face).surface)
                {
                    assert_eq!(
                        model.face(face).sense,
                        if height > 0 {
                            Sense::Same
                        } else {
                            Sense::Reversed
                        }
                    );
                    for boundary in [SplineBoundary::VMin, SplineBoundary::VMax] {
                        assert_eq!(model.edges_of_body(body).iter().filter(|edge| matches!(model.geom().curve(model.edge(**edge).curve), Curve::BSpline(curve) if curve.is_boundary_of(span.carrier(), boundary))).count(), 1);
                    }
                }
            }
        }
    }
}

#[test]
fn oblique_spline_extrusions_keep_the_named_refusal() {
    let profile = validate(&[arch(false)]).unwrap();
    let direction =
        catallax_shapes::RVector3::new(Rational::one(), Rational::zero(), Rational::one());
    assert_eq!(
        extrude(
            &mut Model::new(),
            &profile,
            &StepId::new("sweep").unwrap(),
            &Rational::from_i64(2),
            Some(&direction)
        )
        .unwrap_err()
        .code,
        "extrude.bezier-oblique"
    );
}

#[test]
fn band_checks_reject_false_curve_incidence_and_wrong_or_missing_trims() {
    use catallax_shapes::{BSplineCurve, Curve, Surface};
    use catallax_topology::{BezierExtrusionBand, Coedge, Loop as TopoLoop, Sense};
    let profile = validate(&[arch(false)]).unwrap();
    let mut model = Model::new();
    let body = extrude(
        &mut model,
        &profile,
        &StepId::new("sweep").unwrap(),
        &Rational::from_i64(2),
        None,
    )
    .unwrap();
    let face = model
        .faces_of_body(body)
        .into_iter()
        .find(|face| {
            matches!(
                model.geom().surface(model.face(*face).surface),
                Surface::BSpline(_)
            )
        })
        .unwrap();
    let definition = model.face(face).clone();
    let TopoLoop::Coedges(coedges) = model.loop_(definition.loops[0]) else {
        unreachable!()
    };
    for alteration in 0..4 {
        let mut damaged = model.clone();
        let mut changed = definition.clone();
        let mut coedges = coedges.clone();
        match alteration {
            0 => changed.sense = Sense::Reversed,
            1 => {
                coedges.pop();
            }
            2 => coedges.push(coedges[0]),
            _ => {
                let (index, coedge) = coedges
                    .iter()
                    .enumerate()
                    .find(|(_, c)| {
                        matches!(
                            damaged
                                .geom()
                                .curve(damaged.edge(damaged.coedge(**c).edge).curve),
                            Curve::BSpline(_)
                        )
                    })
                    .unwrap();
                let coedge = damaged.coedge(*coedge).clone();
                let mut edge = damaged.edge(coedge.edge).clone();
                let Curve::BSpline(curve) = damaged.geom().curve(edge.curve) else {
                    unreachable!()
                };
                let mut controls = curve.controls().cloned().collect::<Vec<_>>();
                controls[1].x = &controls[1].x + &catallax_numbers::Real::one();
                let false_curve = BSplineCurve::new(
                    curve.degree(),
                    curve.knots().to_vec(),
                    controls,
                    curve.weights().cloned().collect(),
                    &mut Decider::default(),
                )
                .unwrap();
                edge.curve = damaged.geom_mut().add_curve(Curve::BSpline(false_curve));
                let edge = damaged.add_edge(edge);
                coedges[index] = damaged.add_coedge(Coedge {
                    edge,
                    sense: coedge.sense,
                });
            }
        }
        changed.loops = vec![damaged.add_loop(TopoLoop::Coedges(coedges))];
        let changed = damaged.add_face(changed);
        assert!(BezierExtrusionBand::of_face(&damaged, changed, &mut Decider::default()).is_err());
    }
}

#[test]
fn signed_spline_prisms_preview_and_mesh_with_shared_pickable_edges() {
    for height in [-2, 2] {
        let profile = validate(&[arch(false)]).unwrap();
        let mut model = Model::new();
        let body = extrude(
            &mut model,
            &profile,
            &StepId::new("sweep").unwrap(),
            &Rational::from_i64(height),
            None,
        )
        .unwrap();
        let mesh =
            catallax_mesh::mesh_body(&model, body, &Default::default(), &mut Decider::default())
                .unwrap();
        assert!(mesh.marks.is_empty());
        assert_eq!(mesh.value.faces.len(), 4);
        assert_eq!(mesh.value.edges.len(), 6);
        assert!(
            mesh.value
                .faces
                .iter()
                .all(|face| !face.triangles.is_empty())
        );
        assert!(mesh.value.edges.iter().all(|edge| edge.vertices.len() >= 2));
        let preview = catallax_mesh::preview_mesh_body(
            &model,
            body,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(preview.faces.len(), 4);
        assert_eq!(preview.edges.len(), 6);
    }
}

#[test]
fn exact_copies_keep_band_domains_and_outward_mesh_orientation() {
    use catallax_shapes::{Isometry, RPoint3, RVector3, Rotation};
    let r = Rational::from_i64;
    let profile = validate(&[arch(false)]).unwrap();
    let mut model = Model::new();
    let body = extrude(
        &mut model,
        &profile,
        &StepId::new("sweep").unwrap(),
        &r(2),
        None,
    )
    .unwrap();
    let pose = catallax_shapes::Pose::from_rational(
        Rotation::from_quaternion(r(1), r(0), r(0), r(1)).unwrap(),
        RVector3::new(r(7), r(11), r(13)),
    );
    for (index, map) in [
        Isometry::translation(RVector3::new(r(7), r(11), r(13))),
        Isometry::reflection(&RPoint3::origin(), &RVector3::new(r(1), r(0), r(0))).unwrap(),
        Isometry::from(&pose),
    ]
    .into_iter()
    .enumerate()
    {
        let copy = catallax_topology::transformed_body(&mut model, body, &map, Clone::clone);
        assert!(model.check(copy).is_empty());
        assert!(
            catallax_topology::BezierExtrusionBand::in_body(&model, copy, &mut Decider::default())
                .unwrap()
        );
        let mesh =
            catallax_mesh::mesh_body(&model, copy, &Default::default(), &mut Decider::default())
                .unwrap();
        assert!(mesh.marks.is_empty());
        assert_eq!(mesh.value.faces.len(), 4);
        assert_eq!(mesh.value.edges.len(), 6);
        let bounds = catallax_measure::bounding_box(
            &model,
            copy,
            &catallax_measure::Accuracy {
                relative: Rational::zero(),
                absolute: Rational::parse("1e-9").unwrap(),
            },
            &mut Decider::default(),
        )
        .unwrap();
        assert!(bounds.marks.is_empty());
        let bounds = bounds.value.expect("a solid has finite bounds");
        let (low, high) = [
            (["7", "11", "13"], ["11", "14", "15"]),
            (["-4", "0", "0"], ["0", "3", "2"]),
            (["4", "11", "13"], ["7", "15", "15"]),
        ][index];
        for axis in 0..3 {
            let low = Rational::parse(low[axis]).unwrap();
            let high = Rational::parse(high[axis]).unwrap();
            let tolerance = Rational::parse("1e-9").unwrap();
            assert!(bounds.min[axis] <= low && high <= bounds.max[axis]);
            assert!(&low - &bounds.min[axis] <= tolerance);
            assert!(&bounds.max[axis] - &high <= tolerance);
        }
    }
}

fn encloses(actual: &catallax_measure::Enclosure, expected: &str) {
    let expected = Rational::parse(expected).unwrap();
    assert!(
        actual.contains(&expected),
        "{actual:?} does not enclose {expected}"
    );
    assert!(actual.width() <= Rational::parse("1e-9").unwrap());
}

#[test]
fn signed_curved_prisms_enclose_independent_volume_area_and_centroid_values() {
    let accuracy = catallax_measure::Accuracy {
        relative: Rational::zero(),
        absolute: Rational::parse("1e-9").unwrap(),
    };
    for height in [-2, 2] {
        let profile = validate(&[arch(false)]).unwrap();
        let mut model = Model::new();
        let body = extrude(
            &mut model,
            &profile,
            &StepId::new("sweep").unwrap(),
            &Rational::from_i64(height),
            None,
        )
        .unwrap();
        let measured = catallax_measure::volume_area_centroid(
            &model,
            body,
            &accuracy,
            &mut Decider::default(),
        )
        .unwrap();
        assert!(measured.marks.is_empty());
        encloses(&measured.value.volume, "96/5");
        // The cubic has polynomial speed and exact arc length eight.
        encloses(&measured.value.area, "216/5");
        for (actual, expected) in measured.value.centroid.unwrap().iter().zip([
            "2",
            "9/7",
            if height > 0 { "1" } else { "-1" },
        ]) {
            encloses(actual, expected);
        }
    }
}

#[test]
fn cubic_prism_second_moments_include_the_curved_cap() {
    let profile = validate(&[arch(false)]).unwrap();
    let mut model = Model::new();
    let body = extrude(
        &mut model,
        &profile,
        &StepId::new("sweep").unwrap(),
        &Rational::from_i64(2),
        None,
    )
    .unwrap();
    let measured = catallax_measure::mass_properties(
        &model,
        body,
        &catallax_measure::Accuracy {
            relative: Rational::zero(),
            absolute: Rational::parse("1e-9").unwrap(),
        },
        &mut Decider::default(),
    )
    .unwrap();
    assert!(measured.marks.is_empty());
    let inertia = measured.value.inertia.unwrap();
    for (i, row) in inertia.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            encloses(
                value,
                if i == j {
                    ["4544/245", "10016/385", "17120/539"][i]
                } else {
                    "0"
                },
            );
        }
    }
}

#[test]
fn a_curved_hole_subtracts_its_exact_volume_and_adds_its_wall_area() {
    let hole = Loop(vec![
        Segment::Bezier {
            id: id("inner-arch"),
            controls: vec![p(1, 3), p(3, 3)],
            to: p(3, 1),
        },
        Segment::Line {
            id: id("inner-base"),
            to: p(1, 1),
        },
    ]);
    let profile = validate(&[arch(false), hole]).unwrap();
    let mut model = Model::new();
    let body = extrude(
        &mut model,
        &profile,
        &StepId::new("sweep").unwrap(),
        &Rational::from_i64(2),
        None,
    )
    .unwrap();
    let measured = catallax_measure::volume_area_centroid(
        &model,
        body,
        &catallax_measure::Accuracy {
            relative: Rational::zero(),
            absolute: Rational::parse("1e-9").unwrap(),
        },
        &mut Decider::default(),
    )
    .unwrap();
    assert!(measured.marks.is_empty());
    encloses(&measured.value.volume, "72/5");
    encloses(&measured.value.area, "252/5");
}
