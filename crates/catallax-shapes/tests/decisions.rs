#![allow(missing_docs)]
mod common;
use catallax_numbers::*;
use catallax_shapes::*;
use common::cases::*;
use common::*;

fn supports_z0(store: &mut GeomStore) -> Vec<SurfaceId> {
    // Distinct handles, same plane. Every supplied point has z=0.
    (0..2)
        .map(|_| {
            store.add_surface(Surface::Plane(
                Plane::new(
                    (p("0", "0", "0")).into(),
                    (v("0", "0", "1")).into(),
                    (v("1", "0", "0")).into(),
                    &mut catallax_numbers::Decider::default(),
                )
                .unwrap(),
            ))
        })
        .collect()
}
fn check_small_side(out: Outcome<Side>, expected: i8, lower: &Rational) {
    match out {
        Outcome::NotExact(n) => {
            assert!(n.bound >= *lower);
            assert!(n.code.starts_with("numbers."));
        }
        other => side_answer(other, expected),
    }
}

blind_test!(cases_16_17_18_21_rational_sides_and_tiny_budget, {
    let translation = Pose::new(
        Rotation::identity(),
        v("1000000", "1000000", "1000000").into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap();
    for (s, point, expected) in side_cases() {
        assert_eq!(sign(&implicit_oracle(&s, &point)), expected);
        for pose in [Pose::identity(), translation.clone()] {
            let mut store = GeomStore::new();
            let sid = store.add_surface(s.transformed(&pose));
            let point = pose.apply_point(&point.to_real());
            let pid = store.add_point(catallax_shapes::PointDef::Explicit(point.clone()));
            side_answer(store.side(sid, pid, &mut Decider::default()), expected);
            side_answer(
                store.side_of(sid, &point, &mut Decider::default()),
                expected,
            );
            check_small_side(store.side(sid, pid, &mut tiny()), expected, &r("0"));
        }
    }
    let e = r("1e-30");
    for (s, point, distance) in [
        (
            cylinder("5"),
            rp([r("3") * (&r("1") + &e), r("4") * (&r("1") + &e), r("7")]),
            &r("5") * &e,
        ),
        (
            sphere("3"),
            rp([
                &r("1") + &e,
                r("2") * (&r("1") + &e),
                r("2") * (&r("1") + &e),
            ]),
            &r("3") * &e,
        ),
        (plane(), rp([r("5"), r("7"), &r("1") + &e]), e.clone()),
    ] {
        let mut store = GeomStore::new();
        let sid = store.add_surface(s);
        let pid = explicit(&mut store, point);
        check_small_side(store.side(sid, pid, &mut tiny()), 1, &distance);
    }
});

blind_test!(cases_19_20_closed_and_isolated_points, {
    let mut store = GeomStore::new();
    let cylinder = store.add_surface(cylinder("1"));
    let supports = supports_z0(&mut store);
    let root = real("1/2").sqrt();
    for (offset, expected) in [(r("0"), 0), (r("1e-30"), 1), (-r("1e-30"), -1)] {
        let point = Point3::new(
            &root + &Real::from_rational(offset),
            root.clone(),
            real("0"),
        );
        let id = closed_on(&mut store, point.clone(), supports.clone());
        side_answer(store.side(cylinder, id, &mut Decider::default()), expected);
        check_small_side(store.side(cylinder, id, &mut tiny()), expected, &r("0"));
        point_eq(&store.position(id), &point);
    }
    let s = store.add_surface(sphere("1"));
    let support = Surface::Plane(
        Plane::new(
            p("0", "0", "0").into(),
            v("1", "-1", "0").into(),
            v("0", "0", "1").into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    );
    let supports = vec![
        store.add_surface(support.clone()),
        store.add_surface(support),
    ];
    let point = Point3::new(real("1/2"), real("1/2"), root);
    let id = closed_on(&mut store, point, supports);
    coincident(store.side(s, id, &mut Decider::default()), None);
    check_small_side(store.side(s, id, &mut tiny()), 0, &r("0"));
    for positive in [true, false] {
        let root = sphere_root(positive);
        let supports: Vec<_> = [p("0", "0", "0"), p("6", "0", "0"), p("0", "6", "0")]
            .into_iter()
            .map(|center| {
                store.add_surface(Surface::Sphere(
                    Sphere::new(
                        center.into(),
                        r("5").into(),
                        v("0", "0", "1").into(),
                        v("1", "0", "0").into(),
                        &mut catallax_numbers::Decider::default(),
                    )
                    .unwrap(),
                ))
            })
            .collect();
        let id = store.add_point(PointDef::Meet(Meet {
            of: MeetOf::Surfaces(supports.clone()),
            position: Position::Isolated(root),
        }));
        let point = store.position(id);
        // side_of deliberately bypasses the construction incidences.
        for sid in supports {
            coincident(store.side_of(sid, &point, &mut Decider::default()), None);
            check_small_side(store.side_of(sid, &point, &mut tiny()), 0, &r("0"));
        }
        let plane = store.add_surface(Surface::Plane(
            Plane::new(
                p("0", "0", "0").into(),
                v("0", "0", "1").into(),
                v("1", "0", "0").into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ));
        side_answer(
            store.side_of(plane, &point, &mut Decider::default()),
            if positive { 1 } else { -1 },
        );
        check_small_side(
            store.side_of(plane, &point, &mut tiny()),
            if positive { 1 } else { -1 },
            &r("0"),
        );
        let s = store.add_surface(Surface::Sphere(
            Sphere::new(
                p("3", "3", "0").into(),
                r("3").into(),
                v("0", "0", "1").into(),
                v("1", "0", "0").into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ));
        side_answer(store.side_of(s, &point, &mut Decider::default()), -1);
        check_small_side(store.side_of(s, &point, &mut tiny()), -1, &r("0"));
    }
});

blind_test!(cases_22_23_deadline_and_undefined_coordinate, {
    let mut store = GeomStore::new();
    let s = store.add_surface(sphere("1"));
    let point = explicit(&mut store, p("1", "0", "0"));
    refusal(store.side(s, point, &mut expired()), "numbers.deadline");
    store.record_fact(fact(
        FactKind::Incident,
        GeomKey::Point(point),
        GeomKey::Surface(s),
    ));
    coincident(
        store.side(s, point, &mut expired()),
        Some(ZeroProof::ConstructionFact),
    );
    refusal(
        store.side_of(s, &p("1", "0", "0").to_real(), &mut expired()),
        "numbers.deadline",
    );
    let s2 = real("2").sqrt();
    let bad = real("1") / (&s2 * &s2 - real("2"));
    let point = Point3::new(bad.clone(), real("0"), real("0"));
    let expected = Decider::default().sign(&(bad.square() - real("1")));
    assert_eq!(
        store.side_of(s, &point, &mut Decider::default()),
        expected.map(|sign| match sign {
            StrictSign::Positive => Side::Front,
            StrictSign::Negative => Side::Back,
        })
    );
    refusal(
        store.side_of(s, &point, &mut Decider::default()),
        "numbers.division-by-zero",
    );
});

blind_test!(cases_24_25_26_29_analytic_curve_incidence, {
    for (curve, point, on) in curve_cases() {
        let mut store = GeomStore::new();
        let c = store.add_curve(curve);
        let p = explicit(&mut store, point);
        incidence(store.point_on_curve(c, p, &mut Decider::default()), on);
        match store.point_on_curve(c, p, &mut tiny()) {
            Outcome::NotExact(n) => assert!(n.bound >= r("0")),
            other => incidence(other, on),
        }
    }
    let e = r("1e-30");
    let mut store = GeomStore::new();
    let c = store.add_curve(circle("25"));
    for point in [
        rp([r("3"), r("4"), e.clone()]),
        rp([r("3"), &r("4") + &e, r("0")]),
    ] {
        let id = explicit(&mut store, point.clone());
        if let Outcome::NotExact(n) = store.point_on_curve(c, id, &mut tiny()) {
            if point.z != r("0") {
                assert!(n.bound >= e);
            } else {
                // distance=sqrt(9+(4+e)^2)-5, hence b>=distance iff (5+b)^2>=9+(4+e)^2.
                assert!((r("5") + n.bound).pow(2).unwrap() >= r("9") + point.y.pow(2).unwrap());
            }
        }
    }
});

blind_test!(case_27_intersection_carrier_not_branch, {
    let mut store = GeomStore::new();
    let s = store.add_surface(sphere("5"));
    let plane = store.add_surface(Surface::Plane(
        Plane::new(
            p("0", "0", "3").into(),
            v("0", "0", "1").into(),
            v("1", "0", "0").into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let witness = explicit(&mut store, p("4", "0", "3"));
    let c = store.add_curve(Curve::Intersection(IntersectionCurve {
        surfaces: [s, plane],
        witness,
    }));
    assert_eq!(store.curve(c).kind(), CurveKind::Intersection);
    assert_eq!(store.curve(c).is_closed(), None);
    assert!(
        store
            .curve(c)
            .evaluate(&ball("0", 64), &catallax_numbers::Decider::default())
            .unwrap()
            .is_none()
    );
    assert!(store.curve(c).transformed(&Pose::identity()).is_none());
    point_eq(
        &store.curve(c).parameter_origin(&store),
        &p("4", "0", "3").to_real(),
    );
    for point in [p("0", "4", "3"), p("-4", "0", "3")] {
        let id = explicit(&mut store, point);
        coincident(store.point_on_curve(c, id, &mut Decider::default()), None);
    }
    let near = explicit(&mut store, rp([r("0"), r("4"), r("3") + r("1e-30")]));
    incidence(
        store.point_on_curve(c, near, &mut Decider::default()),
        false,
    );
    match store.point_on_curve(c, near, &mut tiny()) {
        Outcome::NotExact(n) => assert!(n.bound >= r("1e-30")),
        other => incidence(other, false),
    }
    let meet = closed_on(&mut store, p("0", "4", "3").to_real(), vec![s, plane]);
    coincident(
        store.point_on_curve(c, meet, &mut expired()),
        Some(ZeroProof::ConstructionFact),
    );
    let s = store.add_surface(cylinder("1"));
    let plane = store.add_surface(Surface::Plane(
        Plane::new(
            p("1/2", "0", "0").into(),
            v("1", "0", "0").into(),
            v("0", "1", "0").into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let y = real("3").sqrt() / real("2");
    let w = closed_on(
        &mut store,
        Point3::new(real("1/2"), y.clone(), real("0")),
        vec![s, plane],
    );
    let c = store.add_curve(Curve::Intersection(IntersectionCurve {
        surfaces: [s, plane],
        witness: w,
    }));
    // Supports z=5, duplicated, avoid answering this branch test through facts.
    let support = Surface::Plane(
        Plane::new(
            p("0", "0", "5").into(),
            v("0", "0", "1").into(),
            v("1", "0", "0").into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    );
    let supports = vec![
        store.add_surface(support.clone()),
        store.add_surface(support),
    ];
    let other = closed_on(
        &mut store,
        Point3::new(real("1/2"), -y, real("5")),
        supports,
    );
    let mut d = Decider::default();
    coincident(store.point_on_curve(c, other, &mut d), None);
    assert_eq!(d.stats().coincident_fact, 0);
});

blind_test!(cases_28_29_same_point_representations, {
    let mut store = GeomStore::new();
    let supports = supports_z0(&mut store);
    let a = explicit(&mut store, p("1/3", "0", "0"));
    let b = closed_on(
        &mut store,
        Point3::new(real("1") / real("3"), real("0"), real("0")),
        supports.clone(),
    );
    coincident(store.same_point(a, b, &mut Decider::default()), None);
    let radical = Point3::new(real("2").sqrt() * real("1/2").sqrt(), real("1"), real("0"));
    let a = explicit(&mut store, p("1", "1", "0"));
    let b = closed_on(&mut store, radical, supports);
    coincident(store.same_point(a, b, &mut Decider::default()), None);
    let near = explicit(&mut store, rp([r("1") + r("1e-30"), r("1"), r("0")]));
    incidence(store.same_point(near, b, &mut Decider::default()), false);
    match store.same_point(near, b, &mut tiny()) {
        Outcome::NotExact(n) => assert!(n.bound >= r("1e-30")),
        other => incidence(other, false),
    }
    let root = circle_root();
    let s = store.add_surface(sphere("1"));
    let plane = store.add_surface(Surface::Plane(
        Plane::new(
            p("0", "0", "0").into(),
            v("1", "-1", "0").into(),
            v("0", "0", "1").into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let a = store.add_point(PointDef::Meet(Meet {
        of: MeetOf::Surfaces(vec![s, plane]),
        position: Position::Isolated(root),
    }));
    let b = closed_on(
        &mut store,
        Point3::new(real("1/2").sqrt(), real("1/2").sqrt(), real("0")),
        vec![s, plane],
    );
    coincident(store.same_point(a, b, &mut Decider::default()), None);
    match store.same_point(a, b, &mut tiny()) {
        Outcome::NotExact(n) => assert!(n.bound >= r("0")),
        other => coincident(other, None),
    }
    let a = explicit(&mut store, p("1", "2", "3"));
    let b = explicit(&mut store, p("1", "2", "3"));
    assert_ne!(a, b);
    let mut d = Decider::default();
    coincident(store.same_point(a, b, &mut d), None);
    assert_eq!(d.stats().decisions, 3);
    assert_eq!(d.stats().coincident_fact, 0);
    coincident(
        store.same_point(a, a, &mut expired()),
        Some(ZeroProof::ConstructionFact),
    );
});

blind_test!(cases_30_to_36_same_surfaces_and_normals, {
    for (a, b, same) in surface_pairs() {
        let mut store = GeomStore::new();
        let a = store.add_surface(a);
        let b = store.add_surface(b);
        let mut d = expired();
        let out = store.same_surface(a, b, &mut d);
        if same {
            coincident(out, None);
        } else {
            incidence(out, false);
        }
    }
    for (a, b, at, expected) in normal_cases() {
        let mut store = GeomStore::new();
        let a = store.add_surface(a);
        let b = store.add_surface(b);
        let at = explicit(&mut store, at);
        let mut d = Decider::default();
        let out = store.normals_agree(a, b, at, &mut d);
        match expected {
            0 => coincident(out, None),
            1 => assert_eq!(out, Outcome::Proven(StrictSign::Positive)),
            -1 => assert_eq!(out, Outcome::Proven(StrictSign::Negative)),
            _ => unreachable!(),
        }
        assert_eq!(d.stats().decisions, 1);
        refusal(
            store.normals_agree(a, b, at, &mut expired()),
            "numbers.deadline",
        );
    }
});

blind_test!(
    numeric_incidence_rejects_off_geometry_and_obeys_deadlines,
    {
        for (curve, point) in [
            (line(), p("5", "0", "4")),
            (line(), p("6", "0", "3")),
            (circle("25"), p("3", "4", "1")),
            (circle("25"), p("3", "5", "0")),
            (ellipse(), p("3", "12/5", "1")),
            (ellipse(), p("6", "0", "0")),
        ] {
            let mut store = GeomStore::new();
            let c = store.add_curve(curve);
            let p = explicit(&mut store, point);
            let mut d = Decider::default();
            incidence(store.point_on_curve(c, p, &mut d), false);
            let mut d = expired();
            refusal(store.point_on_curve(c, p, &mut d), "numbers.deadline");
        }
        for point in [p("1", "0", "0"), p("0", "1", "0"), p("0", "0", "1")] {
            let mut store = GeomStore::new();
            let a = explicit(&mut store, p("0", "0", "0"));
            let b = explicit(&mut store, point);
            let mut d = Decider::default();
            incidence(store.same_point(a, b, &mut d), false);
            let mut d = expired();
            refusal(store.same_point(a, b, &mut d), "numbers.deadline");
        }
    }
);
