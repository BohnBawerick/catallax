#![allow(missing_docs)]
mod common;
use catallax_numbers::*;
use catallax_shapes::*;
use common::*;

fn starved() -> Decider {
    Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        max_work: 0,
    })
}
fn delta() -> Real {
    // c² < 2 < (c+10^-33)², checked independently in oracle.rs.
    // The positive gap is too small for a binary64 sign filter.
    real("2").sqrt() - real("1.414213562373095048801688724209698")
}
fn nonnegative(x: Real) {
    match Decider::default().sign(&x) {
        Outcome::Proven(StrictSign::Positive) | Outcome::ProvenCoincident(_) => (),
        out => panic!("bound is too small: {out:?}"),
    }
}
fn bounded<T: std::fmt::Debug>(out: Outcome<T>, distance: Real) -> NotExact {
    match out {
        Outcome::NotExact(n) => {
            assert!(n.bound >= r("0"));
            nonnegative(Real::from_rational(n.bound.clone()) - distance);
            n
        }
        other => panic!("zero-work ambiguous decision must remain bounded: {other:?}"),
    }
}
fn closed_z(store: &mut GeomStore, point: Point3, z: &str) -> PointId {
    let s = Surface::Plane(Plane::new(p("0", "0", z), v("0", "0", "1"), v("1", "0", "0")).unwrap());
    let supports = vec![store.add_surface(s.clone()), store.add_surface(s)];
    closed_on(store, point, supports)
}

blind_test!(not_exact_side_distance_units_for_every_surface, {
    let gap = delta();
    let cases = [
        (
            Surface::Plane(
                Plane::new(p("0", "0", "0"), v("2", "0", "0"), v("0", "1", "0")).unwrap(),
            ),
            Point3::new(gap.clone(), real("0"), real("0")),
            gap.clone(),
            &real("2") * &gap,
        ),
        (
            cylinder("5"),
            Point3::new(&real("5") + &gap, real("0"), real("0")),
            gap.clone(),
            (&real("5") + &gap).square() - real("25"),
        ),
        (
            sphere("3"),
            Point3::new(&real("3") + &gap, real("0"), real("0")),
            gap.clone(),
            (&real("3") + &gap).square() - real("9"),
        ),
        // In the meridian, distance to x=z/2 is delta/sqrt(1+1/4).
        (
            cone(),
            Point3::new(&real("1") + &gap, real("0"), real("2")),
            &gap / &real("5/4").sqrt(),
            (&real("1") + &gap).square() - real("1"),
        ),
        (
            torus(),
            Point3::new(&real("4") + &gap, real("0"), real("0")),
            gap.clone(),
            ((&real("4") + &gap).square() + real("8")).square()
                - real("36") * (&real("4") + &gap).square(),
        ),
    ];
    for (surface, point, distance, value) in cases {
        let mut store = GeomStore::new();
        let sid = store.add_surface(surface);
        let original = bounded(starved().sign(&value), real("0"));
        let mut d = starved();
        let converted = bounded(store.side_of(sid, &point, &mut d), distance);
        assert_eq!(converted.code, original.code);
        assert_eq!(converted.detail, original.detail);
        assert_eq!(d.stats().decisions, 1);
        assert_eq!(d.stats().not_exact, 1);
    }
});

blind_test!(not_exact_curve_and_point_distance_bounds, {
    let gap = delta();
    let mut store = GeomStore::new();
    let line = store.add_curve(Curve::Line(
        Line::new(p("0", "0", "0"), v("0", "2", "0")).unwrap(),
    ));
    let point = closed_z(
        &mut store,
        Point3::new(gap.clone(), real("0"), real("0")),
        "0",
    );
    bounded(
        store.point_on_curve(line, point, &mut starved()),
        gap.clone(),
    );
    for curve in [circle("25"), ellipse()] {
        let id = store.add_curve(curve);
        let point = closed_z(
            &mut store,
            Point3::new(&real("5") + &gap, real("0"), real("0")),
            "0",
        );
        // On the outward major-axis ray, the endpoint (5,0,0) is nearest.
        bounded(store.point_on_curve(id, point, &mut starved()), gap.clone());
    }
    let origin = explicit(&mut store, p("0", "0", "0"));
    let point = closed_z(
        &mut store,
        Point3::new(gap.clone(), &real("2") * &gap, real("0")),
        "0",
    );
    let mut d = starved();
    let n = bounded(
        store.same_point(point, origin, &mut d),
        &real("5").sqrt() * &gap,
    );
    let first = bounded(starved().sign(&gap), real("0"));
    assert_eq!(n.code, first.code);
    assert_eq!(n.detail, first.detail);
    assert_eq!(d.stats().decisions, 3);
    // A later nonzero coordinate overrides an earlier unresolved coordinate.
    let point = closed_z(
        &mut store,
        Point3::new(gap.clone(), real("0"), real("1")),
        "1",
    );
    let filter_budget = BudgetPolicy {
        max_precision_bits: 64,
        max_work: 1,
    };
    let mut d = Decider::new(filter_budget);
    incidence(store.same_point(point, origin, &mut d), false);
    assert_eq!(d.stats().decisions, 3);
    let point = closed_z(&mut store, Point3::new(real("1"), gap, real("0")), "0");
    let mut d = Decider::new(filter_budget);
    incidence(store.same_point(point, origin, &mut d), false);
    assert_eq!(d.stats().decisions, 1);
    refusal(
        store.same_point(point, origin, &mut expired()),
        "numbers.deadline",
    );
    refusal(
        store.point_on_curve(line, point, &mut expired()),
        "numbers.deadline",
    );
});

blind_test!(intersection_bound_is_to_both_carriers_not_the_curve, {
    let gap = delta();
    let mut store = GeomStore::new();
    // A tiny dihedral angle: y=0 and y=epsilon*x. The meeting is the z axis.
    // x=delta, y=0 is delta from the meeting but only epsilon*delta/sqrt(1+epsilon²)
    // from the second carrier. Do not demand a curve-distance bound here.
    let a = store.add_surface(Surface::Plane(
        Plane::new(p("0", "0", "0"), v("0", "1", "0"), v("1", "0", "0")).unwrap(),
    ));
    let b = store.add_surface(Surface::Plane(
        Plane::new(
            p("0", "0", "0"),
            v("-1e-30", "1", "0"),
            v("1", "1e-30", "0"),
        )
        .unwrap(),
    ));
    let witness = explicit(&mut store, p("0", "0", "0"));
    let c = store.add_curve(Curve::Intersection(IntersectionCurve {
        surfaces: [a, b],
        witness,
    }));
    let point = closed_z(
        &mut store,
        Point3::new(gap.clone(), real("0"), real("0")),
        "0",
    );
    let distance = real("1e-30") * gap / (real("1") + real("1e-60")).sqrt();
    bounded(store.point_on_curve(c, point, &mut starved()), distance);
});

blind_test!(
    normals_keep_value_bound_and_decisions_do_not_record_facts,
    {
        let gap = delta();
        let mut store = GeomStore::new();
        let a = store.add_surface(sphere("1"));
        let b = store.add_surface(Surface::Plane(
            Plane::new(p("0", "0", "0"), v("3", "0", "0"), v("0", "1", "0")).unwrap(),
        ));
        let point = closed_z(
            &mut store,
            Point3::new(gap.clone(), real("0"), real("0")),
            "0",
        );
        let before: Vec<_> = store.facts().iter().copied().collect();
        let expected = bounded(starved().sign(&(&real("6") * &gap)), real("0"));
        let actual = bounded(
            store.normals_agree(a, b, point, &mut starved()),
            real("6") * gap,
        );
        // Equivalent expression trees may reach different valid enclosures.
        assert_eq!(actual.code, expected.code);
        assert_eq!(actual.detail, expected.detail);
        let _ = store.side(a, point, &mut Decider::default());
        let _ = store.same_point(point, point, &mut Decider::default());
        let _ = store.same_surface(a, a, &mut Decider::default());
        assert_eq!(store.facts().iter().copied().collect::<Vec<_>>(), before);
    }
);
