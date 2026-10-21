//! Blind surface-pair configurations, construction facts and generated families.
mod common;
use common::oracle::{Expected, Part, Shape};
use common::*;
use std::collections::BTreeSet;

fn compare_form(a: CurveForm, b: &CurveForm) {
    match (a, b) {
        (CurveForm::Line { origin, dir }, CurveForm::Line { origin: o, dir: d }) => {
            point_eq(&origin, o);
            assert_eq!(&dir, d);
        }
        (
            CurveForm::Circle {
                center,
                normal,
                radius_sq,
                ref_dir,
            },
            CurveForm::Circle {
                center: c,
                normal: n,
                radius_sq: r,
                ref_dir: d,
            },
        ) => {
            point_eq(&center, c);
            assert_eq!((&normal, &radius_sq, &ref_dir), (n, r, d));
        }
        (
            CurveForm::Ellipse {
                center,
                normal,
                major_dir,
                major_sq,
                minor_sq,
            },
            CurveForm::Ellipse {
                center: c,
                normal: n,
                major_dir: d,
                major_sq: a,
                minor_sq: b,
            },
        ) => {
            point_eq(&center, c);
            assert_eq!((&normal, &major_dir, &major_sq, &minor_sq), (n, d, a, b));
        }
        (CurveForm::Branch, CurveForm::Branch) => {}
        (a, b) => panic!("wrong form {a:?}, expected {b:?}"),
    }
}
fn expected_form(c: &Curve) -> CurveForm {
    match c {
        Curve::Line(c) => CurveForm::Line {
            origin: c.origin.to_real(),
            dir: c.dir.clone(),
        },
        Curve::Circle(c) => CurveForm::Circle {
            center: c.center.to_real(),
            normal: c.normal.clone(),
            radius_sq: c.radius_sq.clone(),
            ref_dir: c.ref_dir.clone(),
        },
        Curve::Ellipse(c) => CurveForm::Ellipse {
            center: c.center.to_real(),
            normal: c.normal.clone(),
            major_dir: c.major_dir.clone(),
            major_sq: c.major_sq.clone(),
            minor_sq: c.minor_sq.clone(),
        },
        _ => unreachable!(),
    }
}
fn check_point(g: &GeomStore, id: PointId, expected: &Point3, supports: [SurfaceId; 2]) {
    match g.point(id) {
        PointDef::Meet(Meet {
            of,
            position: Position::Closed(_),
        }) => assert_eq!(*of, MeetOf::Surfaces(supports.to_vec())),
        other => panic!("surface result must be closed meeting point: {other:?}"),
    }
    // Expected point lives in a clone; no oracle fact can affect the subject.
    let mut clone = g.clone();
    let expected_id = closed(&mut clone, expected.clone(), &supports);
    coincident(clone.same_point(id, expected_id, &mut Decider::default()));
    for s in supports {
        coincident(g.side(s, id, &mut Decider::default()));
        coincident(g.side_of(s, &g.position(id), &mut Decider::default()));
    }
}
fn fact(
    kind: FactKind,
    a: GeomKey,
    b: GeomKey,
    at: Option<GeomKey>,
    origin: &'static str,
) -> Fact<GeomKey> {
    Fact {
        kind,
        a,
        b,
        at,
        origin: FactOrigin::Construction(origin),
    }
}
fn check(a: &Surface, b: &Surface, d: &mut Decider, deep: bool) -> (GeomStore, SurfacePair) {
    let expected = oracle::surfaces(a, b);
    let mut g = GeomStore::new();
    let a = g.add_surface(a.clone());
    let b = g.add_surface(b.clone());
    let before = fingerprint(&g);
    let stats = d.stats().clone();
    let result = intersect_surfaces(&mut g, a, b, d);
    assert_eq!(d.stats(), &stats, "configuration must perform no decisions");
    if matches!(expected, Expected::Refused) {
        assert_eq!(result.unwrap_err().code, "intersect.unsupported-pair");
        assert_eq!(fingerprint(&g), before);
        return (g, SurfacePair::Disjoint);
    }
    let result = exact(result);
    let mut expected_facts = vec![];
    let mut points = BTreeSet::new();
    let mut curves = BTreeSet::new();
    match (&result, expected) {
        (SurfacePair::Disjoint, Expected::Disjoint) => {}
        (SurfacePair::Coincident { same_sense }, Expected::Coincident(s)) => {
            assert_eq!(*same_sense, s)
        }
        (SurfacePair::Meets(parts), Expected::Meets(expected)) => {
            assert_eq!(parts.len(), expected.len());
            let mut singular_ids: Vec<(Point3, PointId)> = vec![];
            for (part, expected) in parts.iter().zip(expected) {
                let mut record_point = |id: PointId| {
                    if points.insert(id) {
                        for s in [a, b] {
                            expected_facts.push(fact(
                                FactKind::Incident,
                                GeomKey::Point(id),
                                GeomKey::Surface(s),
                                None,
                                "shapes.definition",
                            ));
                        }
                    }
                };
                match (part, expected) {
                    (Component::Point { point }, Part::Point(p)) => {
                        record_point(*point);
                        if deep {
                            check_point(&g, *point, &p, [a, b]);
                        }
                        expected_facts.push(fact(
                            FactKind::Tangent,
                            GeomKey::Surface(a),
                            GeomKey::Surface(b),
                            Some(GeomKey::Point(*point)),
                            "intersect.surfaces",
                        ));
                    }
                    (
                        Component::Curve {
                            curve,
                            contact,
                            singular,
                        },
                        Part::Curve {
                            shape,
                            contact: expected_contact,
                            singular: sp,
                        },
                    ) => {
                        assert!(curves.insert(*curve));
                        assert_eq!(*contact, expected_contact);
                        assert_eq!(singular.len(), sp.len());
                        match *shape {
                            Shape::Analytic(expected) => {
                                assert_eq!(g.curve(*curve), &expected);
                                if deep {
                                    compare_form(
                                        curve_form(&g, *curve, &mut Decider::default()).unwrap(),
                                        &expected_form(&expected),
                                    );
                                    let position = expected.parameter_origin(&g);
                                    let mut query = g.clone();
                                    let p = if let Some(p) = rational_point(&position) {
                                        explicit(&mut query, p)
                                    } else {
                                        closed(&mut query, position, &[a, b])
                                    };
                                    coincident(on_curve(
                                        &query,
                                        *curve,
                                        p,
                                        &mut Decider::default(),
                                    ));
                                }
                                for s in [a, b] {
                                    expected_facts.push(fact(
                                        FactKind::Incident,
                                        GeomKey::Curve(*curve),
                                        GeomKey::Surface(s),
                                        None,
                                        "intersect.surfaces",
                                    ));
                                }
                            }
                            Shape::Intersection { witness, form } => {
                                let Curve::Intersection(c) = g.curve(*curve) else {
                                    panic!("expected intersection curve")
                                };
                                assert_eq!(c.surfaces, [a, b]);
                                record_point(c.witness);
                                if deep {
                                    check_point(&g, c.witness, &witness, [a, b]);
                                    compare_form(
                                        curve_form(&g, *curve, &mut Decider::default()).unwrap(),
                                        &form,
                                    );
                                }
                                for s in [a, b] {
                                    expected_facts.push(fact(
                                        FactKind::Incident,
                                        GeomKey::Curve(*curve),
                                        GeomKey::Surface(s),
                                        None,
                                        "shapes.definition",
                                    ));
                                }
                                expected_facts.push(fact(
                                    FactKind::Incident,
                                    GeomKey::Point(c.witness),
                                    GeomKey::Curve(*curve),
                                    None,
                                    "intersect.surfaces",
                                ));
                            }
                        }
                        for (id, p) in singular.iter().zip(sp) {
                            if points.insert(*id) {
                                for s in [a, b] {
                                    expected_facts.push(fact(
                                        FactKind::Incident,
                                        GeomKey::Point(*id),
                                        GeomKey::Surface(s),
                                        None,
                                        "shapes.definition",
                                    ));
                                }
                            }
                            if deep {
                                check_point(&g, *id, &p, [a, b]);
                                for (earlier, prior_id) in &singular_ids {
                                    if Decider::default()
                                        .sign(&p.minus(earlier).norm_squared())
                                        .is_coincident()
                                    {
                                        assert_eq!(
                                            id, prior_id,
                                            "shared node must have one handle"
                                        );
                                    }
                                }
                                singular_ids.push((p, *id));
                                coincident(on_curve(&g, *curve, *id, &mut Decider::default()));
                            }
                            expected_facts.push(fact(
                                FactKind::Incident,
                                GeomKey::Point(*id),
                                GeomKey::Curve(*curve),
                                None,
                                "intersect.surfaces",
                            ));
                        }
                        if *contact == Contact::Tangent {
                            expected_facts.push(fact(
                                FactKind::Tangent,
                                GeomKey::Surface(a),
                                GeomKey::Surface(b),
                                Some(GeomKey::Curve(*curve)),
                                "intersect.surfaces",
                            ));
                        }
                    }
                    (actual, expected) => {
                        panic!("wrong component {actual:?}, expected {expected:?}")
                    }
                }
            }
        }
        (actual, expected) => panic!("wrong configuration {actual:?}, expected {expected:?}"),
    }
    assert_eq!(g.curve_count(), curves.len());
    assert_eq!(g.point_count(), points.len());
    assert_eq!(g.surface_count(), 2);
    let facts: Vec<_> = g.facts().iter().copied().collect();
    assert_eq!(facts.len(), expected_facts.len());
    for f in expected_facts {
        assert!(facts.contains(&f), "missing fact {f:?}");
    }
    (g, result)
}

blind_test!(cases_01_to_27_both_argument_orders, {
    for case in cases::pairs() {
        eprintln!("{}", case.name);
        check(&case.a, &case.b, &mut Decider::default(), true);
        check(&case.b, &case.a, &mut Decider::default(), true);
    }
});
blind_test!(same_handle_and_repeated_calls, {
    for s in [axis_plane(2, n(0)), zc(5), sphere(pi([0, 0, 0]), n(5))] {
        let mut g = GeomStore::new();
        let a = g.add_surface(s);
        assert_eq!(
            exact(intersect_surfaces(&mut g, a, a, &mut expired())),
            SurfacePair::Coincident { same_sense: true }
        );
        assert_eq!(
            (g.curve_count(), g.point_count(), g.facts().len()),
            (0, 0, 0)
        );
    }
    for case in cases::pairs() {
        if matches!(oracle::surfaces(&case.a, &case.b), Expected::Refused) {
            continue;
        }
        let mut g = GeomStore::new();
        let a = g.add_surface(case.a);
        let b = g.add_surface(case.b);
        let first = exact(intersect_surfaces(&mut g, a, b, &mut Decider::default()));
        let counts = (g.curve_count(), g.point_count(), g.facts().len());
        let second = exact(intersect_surfaces(&mut g, a, b, &mut Decider::default()));
        assert_eq!(
            (g.curve_count(), g.point_count(), g.facts().len()),
            (2 * counts.0, 2 * counts.1, 2 * counts.2)
        );
        if let (SurfacePair::Meets(x), SurfacePair::Meets(y)) = (first, second) {
            assert_eq!(x.len(), y.len());
            for (x, y) in x.iter().zip(y) {
                assert_ne!(x, &y, "calls must not reuse definitions");
                match (x, y) {
                    (Component::Point { point: x }, Component::Point { point: y }) => {
                        point_eq(&g.position(*x), &g.position(y))
                    }
                    (
                        Component::Curve {
                            curve: x,
                            contact: xc,
                            singular: xs,
                        },
                        Component::Curve {
                            curve: y,
                            contact: yc,
                            singular: ys,
                        },
                    ) => {
                        assert_eq!(*xc, yc);
                        assert_eq!(xs.len(), ys.len());
                        for (x, y) in xs.iter().zip(ys) {
                            assert_ne!(*x, y);
                            point_eq(&g.position(*x), &g.position(y));
                        }
                        match (g.curve(*x), g.curve(y)) {
                            (Curve::Intersection(x), Curve::Intersection(y)) => {
                                assert_eq!(x.surfaces, y.surfaces);
                                assert_ne!(x.witness, y.witness);
                                point_eq(&g.position(x.witness), &g.position(y.witness));
                            }
                            (x, y) => assert_eq!(x, y),
                        }
                    }
                    other => panic!("repeat changed component kind: {other:?}"),
                }
            }
        }
    }
});
blind_test!(cases_44_45_scale_and_rational_placement, {
    // Sweep scales and poses separately, keeping expensive radical proofs out
    // of the large configuration family. Every scale and pose still checks
    // witness positions, forms, singular points, fields and construction facts.
    let fixtures = cases::pairs();
    let mut rng = Rng(0x041445);
    let mut count = 0usize;
    for k in -30..=30 {
        let scale = n(10).pow(k).unwrap();
        for c in &fixtures {
            let a = scaled_surface(&c.a, &scale);
            let b = scaled_surface(&c.b, &scale);
            check(&a, &b, &mut Decider::default(), true);
            check(&b, &a, &mut Decider::default(), true);
            count += 2;
        }
    }
    for _ in 0..20 {
        let pose = rng.pose();
        for c in &fixtures {
            check(
                &c.a.transformed(&pose),
                &c.b.transformed(&pose),
                &mut Decider::default(),
                true,
            );
            check(
                &c.b.transformed(&pose),
                &c.a.transformed(&pose),
                &mut Decider::default(),
                true,
            );
            count += 2;
        }
    }
    // Rational homotheties randomise radii, offsets and near-boundary distances.
    // A fresh pose randomises axes and placement. Cycling fixtures visits every
    // table configuration, including both sides of tangent boundaries.
    for i in 0..30_000 {
        let case = &fixtures[i % fixtures.len()];
        let scale =
            Rational::from_ratio((rng.next() % 1000 + 1) as i64, (rng.next() % 97 + 1) as i64)
                .unwrap();
        let pose = rng.pose();
        let a = scaled_surface(&case.a, &scale).transformed(&pose);
        let b = scaled_surface(&case.b, &scale).transformed(&pose);
        if i % 2 == 0 {
            check(&a, &b, &mut Decider::default(), false);
        } else {
            check(&b, &a, &mut Decider::default(), false);
        }
        count += 1;
    }
    assert!(count >= 30_000);
    eprintln!("maths-answered surface questions: {count}");
});
blind_test!(case_46_eight_thread_determinism, {
    for c in cases::pairs() {
        let mut base = GeomStore::new();
        let a = base.add_surface(c.a);
        let b = base.add_surface(c.b);
        let runs: Vec<_> = std::thread::scope(|scope| {
            (0..8)
                .map(|_| {
                    let mut g = base.clone();
                    scope.spawn(move || {
                        let mut d = Decider::default();
                        let result = intersect_surfaces(&mut g, a, b, &mut d);
                        (format!("{result:?}"), fingerprint(&g), d.stats().clone())
                    })
                })
                .collect::<Vec<_>>()
                .into_iter()
                .map(|h| h.join().unwrap())
                .collect()
        });
        for run in &runs[1..] {
            assert_eq!(run, &runs[0]);
        }
    }
});
blind_test!(case_47_surface_configuration_ignores_budget_and_deadline, {
    for c in cases::pairs() {
        check(&c.a, &c.b, &mut tiny(), true);
        check(&c.a, &c.b, &mut starved(), true);
        check(&c.a, &c.b, &mut expired(), true);
    }
});
