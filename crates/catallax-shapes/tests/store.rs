#![allow(missing_docs)]
mod common;
use catallax_numbers::*;
use catallax_shapes::*;
use common::*;

blind_test!(case_45_arenas_handles_and_no_deduplication, {
    let mut arena = Arena::new();
    assert!(arena.is_empty());
    assert_eq!(arena.len(), 0);
    let mut handles = Vec::new();
    for i in 0..100 {
        let h = arena.push(format!("value-{i}"));
        assert_eq!(h.index(), i);
        assert_eq!(format!("{h:?}"), format!("#{i}"));
        handles.push(h);
    }
    assert_eq!(arena.len(), 100);
    assert!(!arena.is_empty());
    assert_eq!(arena.handles(), handles);
    for h in handles {
        assert_eq!(arena.get(h), &format!("value-{}", h.index()));
    }
    let mut store = GeomStore::new();
    let mut clone = store.clone();
    assert_eq!(clone.surface_count(), 0);
    let s0 = store.add_surface(plane());
    let s1 = store.add_surface(plane());
    assert_ne!(s0, s1);
    assert_eq!((s0.index(), s1.index()), (0, 1));
    let c0 = store.add_curve(line());
    let c1 = store.add_curve(line());
    assert_ne!(c0, c1);
    let p0 = explicit(&mut store, p("1", "2", "3"));
    let p1 = explicit(&mut store, p("1", "2", "3"));
    assert_ne!(p0, p1);
    assert_eq!(
        (
            store.surface_count(),
            store.curve_count(),
            store.point_count()
        ),
        (2, 2, 2)
    );
    assert!(store.facts().is_empty());
    let keys = [
        GeomKey::Surface(s0),
        GeomKey::Surface(s1),
        GeomKey::Curve(c0),
        GeomKey::Curve(c1),
        GeomKey::Point(p0),
        GeomKey::Point(p1),
    ];
    for pair in keys.windows(2) {
        assert!(pair[0] < pair[1]);
    }
    clone.add_surface(sphere("1"));
    assert_eq!(clone.surface_count(), 1);
    assert_eq!(store.surface_count(), 2);
});

blind_test!(case_46_construction_incidence_order_and_origin, {
    let mut store = GeomStore::new();
    // All supports contain the origin; both intersection witnesses are regular.
    let s0 = store.add_surface(Surface::Plane(
        Plane::new(
            p("0", "0", "0").into(),
            v("1", "0", "0").into(),
            v("0", "1", "0").into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let s1 = store.add_surface(Surface::Plane(
        Plane::new(
            p("0", "0", "0").into(),
            v("0", "1", "0").into(),
            v("1", "0", "0").into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let s2 = store.add_surface(Surface::Plane(
        Plane::new(
            p("0", "0", "0").into(),
            v("0", "0", "1").into(),
            v("1", "0", "0").into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let w = explicit(&mut store, p("0", "0", "0"));
    assert!(store.facts().is_empty());
    let c0 = store.add_curve(Curve::Intersection(IntersectionCurve {
        surfaces: [s1, s0],
        witness: w,
    }));
    let c1 = store.add_curve(Curve::Intersection(IntersectionCurve {
        surfaces: [s2, s0],
        witness: w,
    }));
    let p0 = closed_on(&mut store, p("0", "0", "0").to_real(), vec![s2, s0, s1]);
    let p1 = store.add_point(PointDef::Meet(Meet {
        of: MeetOf::CurveSurface(c0, s2),
        position: Position::Closed(p("0", "0", "0").to_real()),
    }));
    let p2 = store.add_point(PointDef::Meet(Meet {
        of: MeetOf::Curves(c1, c0),
        position: Position::Closed(p("0", "0", "0").to_real()),
    }));
    let expected = [
        (GeomKey::Curve(c0), GeomKey::Surface(s1)),
        (GeomKey::Curve(c0), GeomKey::Surface(s0)),
        (GeomKey::Curve(c1), GeomKey::Surface(s2)),
        (GeomKey::Curve(c1), GeomKey::Surface(s0)),
        (GeomKey::Point(p0), GeomKey::Surface(s2)),
        (GeomKey::Point(p0), GeomKey::Surface(s0)),
        (GeomKey::Point(p0), GeomKey::Surface(s1)),
        (GeomKey::Point(p1), GeomKey::Curve(c0)),
        (GeomKey::Point(p1), GeomKey::Surface(s2)),
        (GeomKey::Point(p2), GeomKey::Curve(c1)),
        (GeomKey::Point(p2), GeomKey::Curve(c0)),
    ];
    let facts: Vec<_> = store.facts().iter().copied().collect();
    assert_eq!(facts.len(), expected.len());
    for (fact, (a, b)) in facts.iter().zip(expected) {
        assert_eq!(
            *fact,
            Fact {
                kind: FactKind::Incident,
                a,
                b,
                at: None,
                origin: FactOrigin::Construction("shapes.definition")
            }
        );
    }
});

blind_test!(
    case_47_fact_closure_precedes_deadlines_and_costs_one_decision,
    {
        let mut store = GeomStore::new();
        let s = store.add_surface(sphere("5"));
        let s2 = store.add_surface(sphere("5"));
        let plane = store.add_surface(Surface::Plane(
            Plane::new(
                p("0", "0", "3").into(),
                v("0", "0", "1").into(),
                v("1", "0", "0").into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ));
        let point = closed_on(&mut store, p("4", "0", "3").to_real(), vec![s, plane]);
        let witness = explicit(&mut store, p("4", "0", "3"));
        let c = store.add_curve(Curve::Intersection(IntersectionCurve {
            surfaces: [s, plane],
            witness,
        }));
        let analytic = store.add_curve(Curve::Circle(
            Circle::new(
                p("0", "0", "3").into(),
                v("0", "0", "1").into(),
                r("16").into(),
                v("1", "0", "0").into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ));
        let f = fact(
            FactKind::Incident,
            GeomKey::Curve(analytic),
            GeomKey::Surface(s),
        );
        assert!(store.record_fact(f));
        assert!(!store.record_fact(f));
        store.record_fact(fact(
            FactKind::Incident,
            GeomKey::Point(witness),
            GeomKey::Curve(analytic),
        ));
        store.record_fact(fact(
            FactKind::Coincident,
            GeomKey::Surface(s),
            GeomKey::Surface(s2),
        ));
        store.record_fact(fact(
            FactKind::Coincident,
            GeomKey::Point(point),
            GeomKey::Point(witness),
        ));
        let mut d = expired();
        for out in [
            store.side(s, point, &mut d).map(|_| Apart),
            store.side(s, witness, &mut d).map(|_| Apart),
            store.side(s2, witness, &mut d).map(|_| Apart),
            store.point_on_curve(c, point, &mut d),
            store.point_on_curve(analytic, witness, &mut d),
            store.same_point(point, witness, &mut d),
            store.same_surface(s, s, &mut d),
            store.same_surface(s, s2, &mut d),
        ] {
            coincident(out, Some(ZeroProof::ConstructionFact));
        }
        assert_eq!(
            d.stats(),
            &DecisionStats {
                decisions: 8,
                coincident_fact: 8,
                ..DecisionStats::default()
            }
        );
        // Neither side_of nor normals_agree has a fact shortcut.
        refusal(
            store.side_of(s, &store.position(point), &mut expired()),
            "numbers.deadline",
        );
        refusal(
            store.normals_agree(s, s, point, &mut expired()),
            "numbers.deadline",
        );
    }
);

blind_test!(cases_48_49_positions_cache_and_ordered_enclosures, {
    let mut store = GeomStore::new();
    let point = explicit(&mut store, p("1/3", "0", "0"));
    point_eq(&store.position(point), &p("1/3", "0", "0").to_real());
    let cached = store.approx_point(point);
    assert_eq!(
        cached.map(f64::to_bits),
        [r("1/3").to_f64_nearest(), 0.0, 0.0].map(f64::to_bits)
    );
    let s = store.add_surface(sphere("5"));
    let sx = store.add_surface(Surface::Sphere(
        Sphere::new(
            p("6", "0", "0").into(),
            r("5").into(),
            v("0", "0", "1").into(),
            v("1", "0", "0").into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let root = sphere_root(true);
    let id = store.add_point(PointDef::Meet(Meet {
        of: MeetOf::Surfaces(vec![s, sx]),
        position: Position::Isolated(root),
    }));
    let pos = store.position(id);
    real_eq(&pos.x, &real("3"));
    real_eq(&pos.y, &real("3"));
    real_eq(&pos.z, &real("7").sqrt());
    for (cached, coord) in store.approx_point(id).into_iter().zip(reals(&pos)) {
        let a = coord.approx();
        assert!(cached >= a.lo() && cached <= a.hi());
    }
    let width = r("1").mul_pow2(-100);
    let mut actual = Decider::default();
    let enclosure = store.enclose_point(id, &width, &mut actual).unwrap();
    assert!(enclosure[0].contains(&r("3")));
    assert!(enclosure[1].contains(&r("3")));
    // Positive enclosure contains sqrt(7) iff lo² <= 7 <= hi². No rounded decimal oracle.
    assert!(enclosure[2].lo() > r("0"));
    assert!(enclosure[2].lo().pow(2).unwrap() <= r("7"));
    assert!(enclosure[2].hi().pow(2).unwrap() >= r("7"));
    for b in &enclosure {
        assert!(b.hi() - b.lo() <= width);
    }
    vector_eq(
        &store.normal(s, point, &mut actual).unwrap(),
        &v("2/3", "0", "0"),
    );
});
