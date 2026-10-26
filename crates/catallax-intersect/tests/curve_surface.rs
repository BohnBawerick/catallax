//! Blind curve-surface roots, contacts, ordering and refusal transactions.
mod common;
use common::curves::*;
use common::*;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Exact,
    Budget,
    Deadline,
}
fn check(case: &HitCase, scale: &Rational, pose: &Pose, mut d: Decider, mode: Mode) {
    let limited = mode != Mode::Exact;
    let mut g = GeomStore::new();
    let c = case.curve.transformed(scale, pose).add(&mut g);
    let s = g.add_surface(scaled_surface(&case.surface, scale).transformed(pose));
    let before = fingerprint(&g);
    let counts = (
        g.surface_count(),
        g.curve_count(),
        g.point_count(),
        g.facts().len(),
    );
    let answer = intersect_curve_surface(&mut g, c, s, &mut d);
    let answer = match answer {
        Ok(a) => a,
        Err(e) if limited => {
            assert!(
                [
                    "numbers.deadline",
                    "numbers.root-not-certified",
                    "numbers.budget-exhausted",
                    "numbers.precision-limit",
                    "numbers.divisor-not-separated",
                    "numbers.division-by-zero",
                    "numbers.sqrt-of-negative",
                    "numbers.sqrt-not-separated",
                    "numbers.zero-polynomial",
                ]
                .contains(&e.code),
                "{e:?}"
            );
            if mode == Mode::Deadline {
                assert_eq!(e.code, "numbers.deadline");
            }
            assert_eq!(
                fingerprint(&g),
                before,
                "refusals must roll back all definitions and facts"
            );
            return;
        }
        Err(e) => panic!("{}: unexpected refusal {e:?}", case.name),
    };
    assert_eq!((g.surface_count(), g.curve_count()), (counts.0, counts.1));
    let marked = !answer.marks.is_empty();
    assert!(limited || !marked, "default answer must be exact");
    if mode == Mode::Deadline {
        assert!(!marked, "expired decisions may refuse, never mark");
    }
    for mark in answer.marks.iter() {
        assert!(matches!(mark.kind, BoundKind::Distance | BoundKind::Value));
        assert!(mark.bound >= n(0));
        assert!(["numbers.budget-exhausted", "numbers.precision-limit"].contains(&mark.code));
    }
    match (answer.value, &case.expected) {
        (CurveSurface::Disjoint, Hits::Disjoint) | (CurveSurface::Contained, Hits::Contained) => {
            assert_eq!((g.point_count(), g.facts().len()), (counts.2, counts.3));
        }
        (CurveSurface::Meets(actual), Hits::Meets(expected)) => {
            assert!(!actual.is_empty());
            assert_eq!(g.point_count() - counts.2, actual.len());
            if !marked {
                assert_eq!(actual.len(), expected.len());
            }
            let expected: Vec<_> = expected
                .iter()
                .map(|(p, c)| {
                    let k = Real::from_rational(scale.clone());
                    (
                        pose.apply_real_point(&Point3::new(&p.x * &k, &p.y * &k, &p.z * &k)),
                        *c,
                    )
                })
                .collect();
            let mut seen = std::collections::BTreeSet::new();
            for (i, h) in actual.iter().enumerate() {
                assert!(seen.insert(h.point), "duplicate hit handle");
                let PointDef::Meet(m) = g.point(h.point) else {
                    panic!("hits must be meeting definitions")
                };
                if !marked {
                    assert_eq!(m.of, MeetOf::CurveSurface(c, s));
                    assert_eq!(h.contact, expected[i].1);
                    let mut comparison = g.clone();
                    let wanted = comparison.add_point(PointDef::Meet(Meet {
                        of: MeetOf::CurveSurface(c, s),
                        position: Position::Closed(expected[i].0.clone()),
                    }));
                    coincident(comparison.same_point(h.point, wanted, &mut Decider::default()));
                    if h.contact == Contact::Tangent || matches!(g.curve(c), Curve::Line(_)) {
                        assert!(matches!(m.position, Position::Closed(_)));
                    }
                    coincident(g.side_of(s, &g.position(h.point), &mut Decider::default()));
                    for equation in g
                        .curve(c)
                        .equations(&g)
                        .expect("analytic fixture equations")
                    {
                        coincident(
                            Decider::default()
                                .sign(&equation.eval_real(&xyz(&g.position(h.point)))),
                        );
                    }
                } else {
                    // A merged hit may choose either mathematical candidate, never an average.
                    assert!(expected.iter().any(|(p, _)| {
                        let mut comparison = g.clone();
                        let wanted = comparison.add_point(PointDef::Meet(Meet {
                            of: MeetOf::CurveSurface(c, s),
                            position: Position::Closed(p.clone()),
                        }));
                        comparison
                            .same_point(h.point, wanted, &mut Decider::default())
                            .is_coincident()
                    }));
                    if h.contact == Contact::Tangent
                        && expected.iter().all(|(_, c)| *c == Contact::Crossing)
                    {
                        assert!(
                            answer.marks.iter().any(
                                |m| m.kind == BoundKind::Value || m.kind == BoundKind::Distance
                            )
                        );
                        assert!(
                            !g.facts().iter().any(|fact| fact.kind == FactKind::Tangent
                                && ((fact.a == GeomKey::Curve(c)
                                    && fact.b == GeomKey::Surface(s))
                                    || (fact.b == GeomKey::Curve(c)
                                        && fact.a == GeomKey::Surface(s)))
                                && fact.at == Some(GeomKey::Point(h.point))),
                            "bounded contact must not become an exact tangent fact"
                        );
                    }
                    match &m.of {
                        MeetOf::CurveSurface(cc, ss) => assert_eq!((*cc, *ss), (c, s)),
                        MeetOf::Surfaces(ss) => {
                            let Curve::Intersection(branch) = g.curve(c) else {
                                panic!("only branch uncertainty changes the definition")
                            };
                            assert_eq!(ss, &vec![branch.surfaces[0], branch.surfaces[1], s]);
                            assert!(!g.facts().holds(
                                FactKind::Incident,
                                GeomKey::Point(h.point),
                                GeomKey::Curve(c)
                            ));
                            assert!(answer.marks.iter().any(|m| m.kind == BoundKind::Distance));
                        }
                        other => panic!("invalid marked point definition {other:?}"),
                    }
                }
            }
            if marked && actual.len() < expected.len() {
                // Every omitted exact root must be covered by a merge-distance
                // mark to a retained tangent root. Unmerged hits may still cross.
                for (p, _) in &expected {
                    let retained = actual.iter().any(|h| {
                        Decider::default()
                            .sign(&p.minus(&g.position(h.point)).norm_squared())
                            .is_coincident()
                    });
                    if !retained {
                        assert!(
                            actual
                                .iter()
                                .filter(|h| h.contact == Contact::Tangent)
                                .any(|h| {
                                    let distance2 = p.minus(&g.position(h.point)).norm_squared();
                                    answer
                                        .marks
                                        .iter()
                                        .filter(|m| m.kind == BoundKind::Distance)
                                        .any(|m| {
                                            matches!(
                                                Decider::default().sign(
                                                    &(Real::from_rational(&m.bound * &m.bound)
                                                        - distance2.clone())
                                                ),
                                                Outcome::Proven(StrictSign::Positive)
                                                    | Outcome::ProvenCoincident(_)
                                            )
                                        })
                                })
                        );
                    }
                }
            }
        }
        (CurveSurface::Meets(actual), Hits::Disjoint) if marked => {
            // A branch decision may conservatively retain a root from the other
            // component. It must not forge incidence on the selected component.
            let Curve::Intersection(branch) = g.curve(c) else {
                panic!("analytic disjoint answers cannot gain roots")
            };
            assert!(!actual.is_empty());
            assert_eq!(g.point_count() - counts.2, actual.len());
            let cylinder = branch
                .surfaces
                .iter()
                .find_map(|s| match g.surface(*s) {
                    Surface::Cylinder(c) => Some(c),
                    _ => None,
                })
                .unwrap();
            let sphere = branch
                .surfaces
                .iter()
                .find_map(|s| match g.surface(*s) {
                    Surface::Sphere(s) => Some(s),
                    _ => None,
                })
                .unwrap();
            let f = cylinder.origin.plus(&cylinder.axis.scaled(
                &(sphere.center.minus(&cylinder.origin).dot(&cylinder.axis)
                    / cylinder.axis.norm_squared()),
            ));
            for h in actual {
                let PointDef::Meet(m) = g.point(h.point) else {
                    panic!("meeting point required")
                };
                assert_eq!(
                    m.of,
                    MeetOf::Surfaces(vec![branch.surfaces[0], branch.surfaces[1], s])
                );
                assert!(!g.facts().holds(
                    FactKind::Incident,
                    GeomKey::Point(h.point),
                    GeomKey::Curve(c)
                ));
                let p = g.position(h.point);
                for s in [branch.surfaces[0], branch.surfaces[1], s] {
                    coincident(g.side_of(s, &p, &mut Decider::default()));
                }
                let distance2 =
                    p.minus(&f).dot(&cylinder.axis).square() / cylinder.axis.norm_squared();
                assert!(
                    answer
                        .marks
                        .iter()
                        .filter(|m| m.kind == BoundKind::Distance)
                        .any(|m| {
                            matches!(
                                Decider::default().sign(
                                    &(Real::from_rational(&m.bound * &m.bound) - distance2.clone())
                                ),
                                Outcome::Proven(StrictSign::Positive)
                                    | Outcome::ProvenCoincident(_)
                            )
                        })
                );
            }
        }
        (actual, expected) => panic!("{}: {actual:?} instead of {expected:?}", case.name),
    }
}

blind_test!(cases_28_to_32_exact_hits_contacts_order_and_definitions, {
    for case in hits() {
        check(
            &case,
            &n(1),
            &Pose::identity(),
            Decider::default(),
            Mode::Exact,
        );
    }
});
blind_test!(case_33_own_surfaces_are_fact_settled_even_after_deadline, {
    for pair in cases::pairs() {
        let oracle::Expected::Meets(parts) = oracle::surfaces(&pair.a, &pair.b) else {
            continue;
        };
        for part in parts {
            let oracle::Part::Curve { shape, .. } = part else {
                continue;
            };
            let oracle::Shape::Intersection { witness, .. } = shape else {
                continue;
            };
            let mut g = GeomStore::new();
            let c = Carrier::Branch {
                a: pair.a.clone(),
                b: Box::new(pair.b.clone()),
                witness,
            }
            .add(&mut g);
            let Curve::Intersection(branch) = g.curve(c) else {
                unreachable!()
            };
            let supports = branch.surfaces;
            for s in supports {
                let before = fingerprint(&g);
                let mut d = expired();
                assert_eq!(
                    exact(intersect_curve_surface(&mut g, c, s, &mut d)),
                    CurveSurface::Contained
                );
                assert_eq!(d.stats().decisions, 1);
                assert_eq!(d.stats().coincident_fact, 1);
                assert_eq!(d.stats().work, 0);
                assert_eq!(fingerprint(&g), before);
            }
        }
    }
});
blind_test!(cases_44_45_hits_scale_and_placement, {
    let fixtures = hits();
    for exponent in -30..=30 {
        for c in &fixtures {
            check(
                c,
                &n(10).pow(exponent).unwrap(),
                &Pose::identity(),
                Decider::default(),
                Mode::Exact,
            );
        }
    }
    let mut rng = Rng(0x043245);
    for _ in 0..20 {
        let pose = rng.pose();
        for c in &fixtures {
            check(c, &n(1), &pose, Decider::default(), Mode::Exact);
        }
    }
});
blind_test!(case_46_hit_determinism, {
    for case in hits() {
        let mut base = GeomStore::new();
        let c = case.curve.add(&mut base);
        let s = base.add_surface(case.surface);
        let outputs = std::thread::scope(|scope| {
            (0..8)
                .map(|_| {
                    let mut g = base.clone();
                    scope.spawn(move || {
                        let mut d = Decider::default();
                        let out = intersect_curve_surface(&mut g, c, s, &mut d);
                        (format!("{out:?}"), fingerprint(&g), d.stats().clone())
                    })
                })
                .collect::<Vec<_>>()
                .into_iter()
                .map(|h| h.join().unwrap())
                .collect::<Vec<_>>()
        });
        for output in &outputs[1..] {
            assert_eq!(output, &outputs[0]);
        }
    }
});
blind_test!(case_47_limited_and_expired_curve_surface_answers, {
    for c in hits() {
        check(&c, &n(1), &Pose::identity(), tiny(), Mode::Budget);
        check(&c, &n(1), &Pose::identity(), starved(), Mode::Budget);
        check(&c, &n(1), &Pose::identity(), expired(), Mode::Deadline);
    }
});

blind_test!(
    line_component_intersections_do_not_depend_on_witness_offset,
    {
        let mut geom = GeomStore::new();
        let curve = Carrier::Branch {
            a: axis_plane(1, n(0)),
            b: Box::new(axis_plane(2, n(0))),
            witness: Point3::new(n(1).mul_pow2(65_536).into(), rr(0), rr(0)),
        }
        .add(&mut geom);
        let CurveForm::Line { origin, .. } =
            curve_form(&geom, curve, &mut Decider::default()).unwrap()
        else {
            panic!("two intersecting planes have a line carrier");
        };
        coincident(
            Decider::default().sign(
                &origin
                    .minus(&geom.curve(curve).parameter_origin(&geom))
                    .norm_squared(),
            ),
        );
        let surface = geom.add_surface(sphere(pi([0, 0, 0]), n(1)));
        let mut decider = Decider::new(BudgetPolicy {
            max_precision_bits: 64,
            max_work: 2048,
        });
        let CurveSurface::Meets(hits) = exact(intersect_curve_surface(
            &mut geom,
            curve,
            surface,
            &mut decider,
        )) else {
            panic!("the plane intersection crosses the unit sphere twice");
        };
        assert_eq!(hits.len(), 2);
        for (hit, x) in hits.iter().zip([-1, 1]) {
            assert_eq!(hit.contact, Contact::Crossing);
            let expected = explicit(&mut geom, pi([x, 0, 0]));
            coincident(geom.same_point(hit.point, expected, &mut Decider::default()));
        }
    }
);

blind_test!(line_meeting_coordinate_bounds_ignore_parameter_scale, {
    let directions = [
        Vector3::new(rr(3), rr(1), rr(0)),
        Vector3::new(rr(1), Real::from_rational(n(1) / n(3)), rr(0)),
    ];
    let mut expected = None;
    for direction in directions {
        let mut geom = GeomStore::new();
        let curve = geom.add_curve(Curve::Line(Line {
            origin: Point3::origin(),
            dir: direction,
        }));
        let surface = geom.add_surface(sphere(pi([0, 0, 0]), n(1)));
        let CurveSurface::Meets(hits) = exact(intersect_curve_surface(
            &mut geom,
            curve,
            surface,
            &mut Decider::default(),
        )) else {
            panic!("the line crosses the unit sphere twice");
        };
        assert_eq!(hits.len(), 2);
        let bounds: Vec<_> = hits
            .iter()
            .map(|hit| {
                let point = geom.position(hit.point);
                [&point.x, &point.y, &point.z].map(|value| {
                    let ball = Decider::default().enclose_at(value, 128).unwrap();
                    (ball.lo(), ball.hi())
                })
            })
            .collect();
        if let Some(expected) = &expected {
            assert_eq!(&bounds, expected);
        } else {
            expected = Some(bounds);
        }
    }
});

blind_test!(
    offset_circular_sections_keep_roots_contacts_and_containment,
    {
        let crossing = Hits::Meets(vec![
            (pi([3, 4, 0]).to_real(), Contact::Crossing),
            (pi([3, -4, 0]).to_real(), Contact::Crossing),
        ]);
        let tangent = Hits::Meets(vec![(pi([5, 0, 0]).to_real(), Contact::Tangent)]);
        for (offset, radius, expected) in [
            (6, 5, crossing),
            (10, 5, tangent.clone()),
            (2, 3, tangent),
            (11, 5, Hits::Disjoint),
            (1, 3, Hits::Disjoint),
            (0, 3, Hits::Disjoint),
            (0, 5, Hits::Contained),
        ] {
            let case = HitCase {
                name: format!("parallel circular section {offset}/{radius}"),
                curve: Carrier::Analytic(k5()),
                surface: cyl(pi([offset, 0, 7]), vi([0, 0, -3]), n(radius)),
                expected,
            };
            for pose in [Pose::identity(), Rng(0x513).pose()] {
                check(&case, &n(1), &pose, Decider::default(), Mode::Exact);
            }
        }
    }
);

blind_test!(
    algebraic_lines_keep_correlated_hits_on_rational_cylinders,
    {
        let root = rr(2).sqrt();
        let surface = cyl(pi([0, 0, 0]), vi([0, 0, 1]), n(2));
        for (origin, expected) in [
            (
                Point3::new(root.clone(), rr(0), rr(0)),
                Hits::Meets(vec![
                    (Point3::new(root.clone(), -&root, rr(0)), Contact::Crossing),
                    (
                        Point3::new(root.clone(), root.clone(), rr(0)),
                        Contact::Crossing,
                    ),
                ]),
            ),
            (
                Point3::new(rr(2), rr(0), root.clone()),
                Hits::Meets(vec![(
                    Point3::new(rr(2), rr(0), root.clone()),
                    Contact::Tangent,
                )]),
            ),
            (Point3::new(&rr(2) * &root, rr(0), rr(0)), Hits::Disjoint),
        ] {
            let case = HitCase {
                name: "algebraic line on rational cylinder".into(),
                curve: Carrier::Analytic(Curve::Line(
                    Line::new(origin, vi([0, 1, 0]).into(), &mut Decider::default()).unwrap(),
                )),
                surface: surface.clone(),
                expected,
            };
            check(
                &case,
                &n(1),
                &Pose::identity(),
                Decider::default(),
                Mode::Exact,
            );
        }
    }
);

blind_test!(circle_planes_parallel_to_cylinder_axes_keep_tangent_hits, {
    let circle_at = |height: Rational| {
        Carrier::Analytic(Curve::Circle(
            catallax_shapes::Circle::new(
                Point3::new(rr(0), rr(0), height.into()),
                vi([0, 0, 7]).into(),
                rr(9),
                vi([2, 0, 0]).into(),
                &mut Decider::default(),
            )
            .unwrap(),
        ))
    };
    let surface = cyl(pi([0, 0, 0]), vi([0, 3, 0]), n(2));
    for height in [n(3), n(2) + n(1).mul_pow2(-128)] {
        let case = HitCase {
            name: "circle plane strictly beyond cylinder".into(),
            curve: circle_at(height),
            surface: surface.clone(),
            expected: Hits::Disjoint,
        };
        check(
            &case,
            &n(1),
            &Pose::identity(),
            Decider::new(BudgetPolicy {
                max_precision_bits: 128,
                max_work: 512,
            }),
            Mode::Exact,
        );
    }
    let root = rr(5).sqrt();
    for (height, expected) in [
        (
            n(0),
            Hits::Meets(vec![
                (Point3::new(rr(2), root.clone(), rr(0)), Contact::Crossing),
                (Point3::new(rr(-2), root.clone(), rr(0)), Contact::Crossing),
                (Point3::new(rr(-2), -&root, rr(0)), Contact::Crossing),
                (Point3::new(rr(2), -&root, rr(0)), Contact::Crossing),
            ]),
        ),
        (
            n(2),
            Hits::Meets(vec![
                (pi([0, 3, 2]).to_real(), Contact::Tangent),
                (pi([0, -3, 2]).to_real(), Contact::Tangent),
            ]),
        ),
    ] {
        let case = HitCase {
            name: "circle plane intersects or touches cylinder".into(),
            curve: circle_at(height),
            surface: surface.clone(),
            expected,
        };
        for pose in [Pose::identity(), Rng(0x5A13).pose()] {
            check(&case, &n(1), &pose, Decider::default(), Mode::Exact);
        }
    }
});

blind_test!(intersection_supports_preserve_tangent_plane_contacts, {
    let carrier = Carrier::Branch {
        a: cyl(pi([0, 0, 0]), vi([0, 0, -3]), n(5)),
        b: Box::new(cyl(pi([3, 0, 2]), vi([0, 5, 0]), n(1))),
        witness: pi([3, 4, 3]).to_real(),
    };
    for (height, expected) in [
        (n(4), Hits::Disjoint),
        (n(3) + n(1).mul_pow2(-128), Hits::Disjoint),
        (
            n(3),
            Hits::Meets(vec![(pi([3, 4, 3]).to_real(), Contact::Tangent)]),
        ),
    ] {
        let case = HitCase {
            name: "intersection component against a support-parallel plane".into(),
            curve: carrier.clone(),
            surface: axis_plane(2, height),
            expected,
        };
        for pose in [Pose::identity(), Rng(0xC71).pose()] {
            check(&case, &n(1), &pose, Decider::default(), Mode::Exact);
        }
    }
});
