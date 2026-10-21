use super::*;
use catallax_numbers::BudgetPolicy;
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation};

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn vector(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(r(x), r(y), r(z))
}
fn ball(n: i64) -> Ball {
    Ball::from_rational(&r(n), 256)
}
fn exact(value: Rational) -> Enclosure {
    Enclosure {
        lo: value.clone(),
        hi: value,
    }
}
fn actual_width(value: &Enclosure) -> Rational {
    &value.hi - &value.lo
}
fn meets_width(value: &Enclosure, accuracy: &Accuracy) -> bool {
    let magnitude = if value.lo >= r(0) {
        value.lo.clone()
    } else if value.hi <= r(0) {
        -&value.hi
    } else {
        r(0)
    };
    actual_width(value) <= (&accuracy.relative * &magnitude).max(accuracy.absolute.clone())
}
fn contains(enclosure: &Enclosure, expected: &Ball) {
    assert!(
        enclosure.lo <= expected.lo() && enclosure.hi >= expected.hi(),
        "{enclosure:?} does not contain {expected:?}"
    );
    assert!(
        meets_width(enclosure, &Accuracy::default()),
        "accuracy not reached: {enclosure:?}"
    );
}

#[test]
fn rational_blocks_and_empty_need_no_budget() {
    for dimensions in [
        [r(2), r(3), r(5)],
        [
            Rational::parse("1/3").unwrap(),
            r(7),
            Rational::parse("1e-3").unwrap(),
        ],
    ] {
        let pose = Pose::new(
            Rotation::from_quaternion(r(2), r(-3), r(5), r(7)).unwrap(),
            vector(1000000, -2000000, 123),
        );
        let mut model = Model::new();
        let body = catallax_topology::block(&mut model, dimensions.clone(), &pose).unwrap();
        let mut dec = Decider::new(BudgetPolicy {
            max_precision_bits: 0,
            max_work: 0,
        });
        let answer = mass_properties(&model, body, &Accuracy::default(), &mut dec).unwrap();
        assert!(answer.marks.is_empty());
        let [a, b, c] = dimensions;
        assert_eq!(answer.value.volume, exact((&a * &b) * c.clone()));
        assert_eq!(
            answer.value.area,
            exact(r(2) * (&a * &b + &b * &c + &c * &a))
        );
        let expected = pose.apply_point(&RPoint3::new(a / r(2), b / r(2), c / r(2)));
        assert_eq!(
            answer.value.centroid.unwrap(),
            [&expected.x, &expected.y, &expected.z].map(|c| exact(c.clone()))
        );
    }
    let mut model = Model::new();
    let body = catallax_topology::empty(&mut model);
    let result =
        mass_properties(&model, body, &Accuracy::default(), &mut Decider::default()).unwrap();
    assert_eq!(result.value.volume, exact(r(0)));
    assert_eq!(result.value.area, exact(r(0)));
    assert_eq!(result.value.centroid, None);
}

#[test]
fn primitive_closed_forms() {
    let pi = Ball::pi(256);
    for kind in 0..5 {
        let mut model = Model::new();
        let (body, volume, area, centroid) = match kind {
            0 => (
                catallax_topology::cylinder(&mut model, r(2), r(5), &Pose::identity()).unwrap(),
                &pi * &ball(20),
                &pi * &ball(28),
                Rational::parse("5/2").unwrap(),
            ),
            1 => (
                catallax_topology::sphere(&mut model, r(5), &Pose::identity()).unwrap(),
                (&pi * &ball(500)).checked_div(&ball(3)).unwrap(),
                &pi * &ball(100),
                r(0),
            ),
            2 => (
                catallax_topology::cone(&mut model, r(2), r(1), r(2), &Pose::identity()).unwrap(),
                (&pi * &ball(14)).checked_div(&ball(3)).unwrap(),
                &pi * &(ball(5) + ball(3) * ball(5).sqrt().unwrap()),
                Rational::parse("11/14").unwrap(),
            ),
            3 => (
                catallax_topology::cone(&mut model, r(3), r(0), r(4), &Pose::identity()).unwrap(),
                &pi * &ball(12),
                &pi * &ball(24),
                r(1),
            ),
            _ => (
                catallax_topology::torus(&mut model, r(3), r(1), &Pose::identity()).unwrap(),
                pi.square() * ball(6),
                pi.square() * ball(12),
                r(0),
            ),
        };
        assert!(model.check(body).is_empty());
        let result =
            mass_properties(&model, body, &Accuracy::default(), &mut Decider::default()).unwrap();
        contains(&result.value.volume, &volume);
        contains(&result.value.area, &area);
        let center = result.value.centroid.unwrap();
        assert!(
            center[0].contains(&r(0)) && center[1].contains(&r(0)) && center[2].contains(&centroid)
        );
        assert!(
            result.marks.is_empty(),
            "primitive {kind}: {:?}",
            result.marks
        );
    }
}

#[test]
fn hemisphere_and_cap_include_the_pole() {
    let pi = Ball::pi(256);
    for (radius, height, volume, area, center) in [
        (
            2,
            0,
            Rational::parse("16/3").unwrap(),
            r(12),
            Rational::parse("3/4").unwrap(),
        ),
        (
            5,
            3,
            Rational::parse("52/3").unwrap(),
            r(36),
            Rational::parse("48/13").unwrap(),
        ),
    ] {
        let (model, body) = crate::fixtures::cap(radius, height);
        let result =
            mass_properties(&model, body, &Accuracy::default(), &mut Decider::default()).unwrap();
        contains(
            &result.value.volume,
            &(&pi * &Ball::from_rational(&volume, 256)),
        );
        contains(
            &result.value.area,
            &(&pi * &Ball::from_rational(&area, 256)),
        );
        let centroid = result.value.centroid.unwrap();
        assert!(
            centroid[0].contains(&r(0))
                && centroid[1].contains(&r(0))
                && centroid[2].contains(&center)
        );
        assert!(result.marks.is_empty());
    }
}

#[test]
fn zero_work_refuses_transcendental_properties() {
    let mut model = Model::new();
    let body = catallax_topology::sphere(&mut model, r(5), &Pose::identity()).unwrap();
    let mut dec = Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        max_work: 0,
    });
    assert_eq!(
        volume(&model, body, &Accuracy::default(), &mut dec)
            .unwrap_err()
            .code,
        "measure.budget-exhausted"
    );
}

#[test]
fn holes_and_void_shells_are_measured_as_built() {
    let (model, body) = crate::fixtures::void_block();
    let result = mass_properties(
        &model,
        body,
        &Accuracy::default(),
        &mut Decider::new(BudgetPolicy {
            max_precision_bits: 0,
            max_work: 0,
        }),
    )
    .unwrap();
    assert_eq!(result.value.volume, exact(r(992)));
    assert_eq!(result.value.area, exact(r(624)));
    assert_eq!(
        result.value.centroid,
        Some(std::array::from_fn(|_| exact(r(5))))
    );
    let pi = Ball::pi(256);
    let pose = Pose::new(
        Rotation::from_quaternion(r(7), r(-3), r(2), r(5)).unwrap(),
        vector(1000000, -12345, 987654),
    );
    for scale in ["1", "1/3", "1e-30", "1e30"] {
        let scale = Rational::parse(scale).unwrap();
        for pose in [&Pose::identity(), &pose] {
            let (model, body) = crate::fixtures::through_hole(&scale, pose);
            let result =
                mass_properties(&model, body, &Accuracy::default(), &mut Decider::default())
                    .unwrap();
            contains(
                &result.value.volume,
                &((ball(64) - &ball(4) * &pi) * Ball::from_rational(&scale.pow(3).unwrap(), 256)),
            );
            contains(
                &result.value.area,
                &((ball(96) + &ball(6) * &pi) * Ball::from_rational(&scale.pow(2).unwrap(), 256)),
            );
            let center = pose.apply_point(&RPoint3::new(
                &scale * &r(2),
                &scale * &r(2),
                &scale * &r(2),
            ));
            for (actual, expected) in result
                .value
                .centroid
                .unwrap()
                .iter()
                .zip([&center.x, &center.y, &center.z])
            {
                assert!(
                    actual.contains(expected),
                    "centroid {actual:?} excludes {expected:?}"
                );
            }
            assert!(result.marks.is_empty());
        }
    }
}

#[test]
fn quarter_arcs_seam_boundaries_and_ellipses() {
    let pi = Ball::pi(256);
    for kind in 0..3 {
        let (model, body, expected_volume, expected_area, centroid, face_areas) = match kind {
            0 => {
                let (model, body) = crate::fixtures::octant(&Pose::identity());
                (
                    model,
                    body,
                    (&pi * &ball(9)).checked_div(&ball(2)).unwrap(),
                    (&pi * &ball(45)).checked_div(&ball(4)).unwrap(),
                    std::array::from_fn(|_| Rational::parse("9/8").unwrap()),
                    vec![
                        (&pi * &ball(9)).checked_div(&ball(4)).unwrap(),
                        (&pi * &ball(9)).checked_div(&ball(4)).unwrap(),
                        (&pi * &ball(9)).checked_div(&ball(4)).unwrap(),
                        (&pi * &ball(9)).checked_div(&ball(2)).unwrap(),
                    ],
                )
            }
            1 => {
                let (model, body) = crate::fixtures::half_cylinder(&Pose::identity());
                (
                    model,
                    body,
                    &pi * &ball(6),
                    &pi * &ball(10) + ball(12),
                    [r(0), r(0), Rational::parse("3/2").unwrap()],
                    vec![&pi * &ball(2), &pi * &ball(2), ball(12), &pi * &ball(6)],
                )
            }
            _ => {
                let (model, body) = crate::fixtures::steinmetz(&Pose::identity());
                (
                    model,
                    body,
                    ball(16).checked_div(&ball(3)).unwrap(),
                    ball(16),
                    [r(0), r(0), r(0)],
                    vec![ball(4); 4],
                )
            }
        };
        let result =
            mass_properties(&model, body, &Accuracy::default(), &mut Decider::default()).unwrap();
        contains(&result.value.volume, &expected_volume);
        contains(&result.value.area, &expected_area);
        let center = result.value.centroid.unwrap();
        for i in 0..3 {
            if kind != 1 || i != 1 {
                assert!(center[i].contains(&centroid[i]));
            }
            assert!(meets_width(&center[i], &Accuracy::default()));
        }
        if kind == 1 {
            contains(&center[1], &ball(8).checked_div(&(&ball(3) * &pi)).unwrap());
        }
        assert!(
            result.marks.is_empty(),
            "hand case {kind}: {:?}",
            result.marks
        );
        for (face, expected) in model.faces_of_body(body).into_iter().zip(face_areas) {
            contains(
                &face_area(&model, face, &Accuracy::default(), &mut Decider::default())
                    .unwrap()
                    .value,
                &expected,
            );
        }
    }
}

#[test]
fn tight_accuracy_enclosures_reach_requested_width() {
    let accuracy = Accuracy {
        relative: Rational::parse("1e-30").unwrap(),
        absolute: Rational::parse("1e-40").unwrap(),
    };
    let mut model = Model::new();
    let body = catallax_topology::sphere(&mut model, r(5), &Pose::identity()).unwrap();
    let expected = mass_properties(&model, body, &accuracy, &mut Decider::default())
        .unwrap()
        .value;
    assert!(meets_width(&expected.volume, &accuracy) && meets_width(&expected.area, &accuracy));
}

#[test]
fn torus_meridians_need_the_noncontractible_loop_term() {
    let (model, body) = crate::fixtures::half_torus();
    let pi = Ball::pi(256);
    let result =
        mass_properties(&model, body, &Accuracy::default(), &mut Decider::default()).unwrap();
    contains(&result.value.volume, &(pi.square() * ball(3)));
    contains(
        &result.value.area,
        &(pi.square() * ball(6) + &pi * &ball(2)),
    );
    let center = result.value.centroid.unwrap();
    assert!(center[0].contains(&r(0)) && center[2].contains(&r(0)));
    contains(
        &center[1],
        &ball(37).checked_div(&(&pi * &ball(6))).unwrap(),
    );
    assert!(result.marks.is_empty());
}

#[test]
fn limited_precision_marks_property_widths_in_order() {
    let mut model = Model::new();
    let body = catallax_topology::cylinder(&mut model, r(2), r(5), &Pose::identity()).unwrap();
    let accuracy = Accuracy {
        relative: r(0),
        absolute: r(0),
    };
    let mut dec = Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        max_work: 10_000,
    });
    let result = mass_properties(&model, body, &accuracy, &mut dec).unwrap();
    let mut widths = [&result.value.volume, &result.value.area]
        .into_iter()
        .filter(|value| !meets_width(value, &accuracy))
        .map(actual_width)
        .collect::<Vec<_>>();
    let center = result.value.centroid.as_ref().unwrap();
    if center.iter().any(|value| !meets_width(value, &accuracy)) {
        widths.push(center.iter().map(actual_width).max().unwrap());
    }
    let marks = result.marks.iter().collect::<Vec<_>>();
    assert_eq!(marks.len(), widths.len());
    for (mark, width) in marks.into_iter().zip(widths) {
        assert_eq!(mark.code, "measure.accuracy-not-reached");
        assert_eq!(mark.kind, BoundKind::Value);
        assert_eq!(mark.bound, width);
    }
    let pi = Ball::pi(256);
    assert!(result.value.volume.lo <= (&pi * &ball(20)).lo());
    assert!(result.value.volume.hi >= (&pi * &ball(20)).hi());
    assert!(result.value.area.lo <= (&pi * &ball(28)).lo());
    assert!(result.value.area.hi >= (&pi * &ball(28)).hi());
    assert!(center[0].contains(&r(0)) && center[1].contains(&r(0)));
    assert!(center[2].contains(&Rational::parse("5/2").unwrap()));
}

#[test]
fn centroid_separation_ignores_the_absolute_floor() {
    use catallax_topology::{Body, Face, Lump, Sense, Shell};
    let mut model = Model::new();
    let outer = catallax_topology::sphere(&mut model, r(1), &Pose::identity()).unwrap();
    let inner = catallax_topology::sphere(
        &mut model,
        r(1) - Rational::parse("1e-100").unwrap(),
        &Pose::identity(),
    )
    .unwrap();
    let surface = model.face(model.faces_of_body(inner)[0]).surface;
    let face = model.add_face(Face {
        surface,
        sense: Sense::Reversed,
        loops: Vec::new(),
    });
    let void = model.add_shell(Shell { faces: vec![face] });
    let outside = model.lump(model.body(outer).lumps[0]).shells[0];
    let lump = model.add_lump(Lump {
        shells: vec![outside, void],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    });
    assert!(model.check(body).is_empty());
    let accuracy = Accuracy {
        relative: r(0),
        absolute: r(1),
    };
    let mut limited = Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        max_work: 10_000,
    });
    assert_eq!(
        mass_properties(&model, body, &accuracy, &mut limited)
            .unwrap_err()
            .code,
        "measure.centroid-not-separated"
    );
    let result = mass_properties(&model, body, &accuracy, &mut Decider::default()).unwrap();
    assert_eq!(
        result.value.centroid,
        Some(std::array::from_fn(|_| exact(r(0))))
    );
}

#[test]
fn fresh_deciders_and_threads_produce_identical_enclosures() {
    let pose = Pose::new(
        Rotation::from_quaternion(r(2), r(7), r(-5), r(3)).unwrap(),
        vector(123456, -987654, 300000),
    );
    let mut model = Model::new();
    let body = catallax_topology::cone(&mut model, r(2), r(1), r(2), &pose).unwrap();
    let expected = mass_properties(&model, body, &Accuracy::default(), &mut Decider::default())
        .unwrap()
        .value;
    std::thread::scope(|scope| {
        let handles = (0..8)
            .map(|_| {
                scope.spawn(|| {
                    mass_properties(&model, body, &Accuracy::default(), &mut Decider::default())
                        .unwrap()
                        .value
                })
            })
            .collect::<Vec<_>>();
        for handle in handles {
            assert_eq!(handle.join().unwrap(), expected);
        }
    });
}

#[test]
fn coplanar_curved_subdivisions_preserve_rational_exactness() {
    let (model, body) = crate::fixtures::coplanar_split_block();
    let mut dec = Decider::new(BudgetPolicy {
        max_precision_bits: 0,
        max_work: 0,
    });
    let result = mass_properties(&model, body, &Accuracy::default(), &mut dec).unwrap();
    assert_eq!(result.value.volume, exact(r(64)));
    assert_eq!(result.value.area, exact(r(96)));
    assert_eq!(
        result.value.centroid,
        Some(std::array::from_fn(|_| exact(r(2))))
    );
    assert!(result.marks.is_empty());
}

#[test]
fn integer_periodic_lifts_keep_the_whole_surface_in_the_enclosure() {
    use catallax_numbers::{Fact, FactKind, FactOrigin};
    use catallax_shapes::{Circle, Curve, GeomKey, Surface, Torus};
    use catallax_topology::{Body, Coedge, Edge, Face, Loop, Lump, Sense, Shell};
    let mut model = Model::new();
    let surface = model.geom_mut().add_surface(Surface::Torus(
        Torus::new(
            RPoint3::origin(),
            vector(0, 0, 1),
            r(3),
            r(1),
            vector(1, 0, 0),
        )
        .unwrap(),
    ));
    let curve = model.geom_mut().add_curve(Curve::Circle(
        Circle::new(
            RPoint3::new(r(3), r(0), r(0)),
            vector(0, -1, 0),
            r(1),
            vector(1, 0, 0),
        )
        .unwrap(),
    ));
    model.geom_mut().record_fact(Fact {
        kind: FactKind::Incident,
        a: GeomKey::Curve(curve),
        b: GeomKey::Surface(surface),
        at: None,
        origin: FactOrigin::Construction("measure.test-torus-seam"),
    });
    let edge = model.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: None,
        end: None,
    });
    let mut loops = Vec::new();
    for sense in [Sense::Same, Sense::Reversed] {
        let coedge = model.add_coedge(Coedge { edge, sense });
        loops.push(model.add_loop(Loop::Coedges(vec![coedge])));
    }
    let face = model.add_face(Face {
        surface,
        sense: Sense::Same,
        loops,
    });
    let shell = model.add_shell(Shell { faces: vec![face] });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    });
    assert!(model.check(body).is_empty());
    let accuracy = Accuracy::default();
    let mut dec = Decider::new(BudgetPolicy {
        max_precision_bits: 256,
        max_work: 1_000_000,
    });
    let result = face_area(&model, face, &accuracy, &mut dec).unwrap();
    let expected = Ball::pi(256).square() * ball(12);
    assert!(result.value.lo <= expected.lo() && result.value.hi >= expected.hi());
    if !meets_width(&result.value, &accuracy) {
        let marks = result.marks.iter().collect::<Vec<_>>();
        assert_eq!(marks.len(), 1);
        assert_eq!(marks[0].code, "measure.accuracy-not-reached");
        assert_eq!(marks[0].bound, actual_width(&result.value));
    }
}
