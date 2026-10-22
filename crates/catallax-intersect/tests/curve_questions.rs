//! Blind closed forms, membership, equality, ordering, tangents and arc points.
mod common;
use common::curves::*;
use common::*;

fn membership(out: Outcome<Apart>, on: bool, limited: bool, min_bound: Rational) {
    match out {
        Outcome::Proven(Apart) => assert!(!on),
        Outcome::ProvenCoincident(_) => assert!(on),
        Outcome::NotExact(x) if limited => {
            assert!(["numbers.budget-exhausted", "numbers.precision-limit"].contains(&x.code));
            assert!(x.bound >= min_bound);
        }
        Outcome::Refused(e) if limited => assert!(
            [
                "numbers.deadline",
                "numbers.budget-exhausted",
                "numbers.root-not-certified",
                "numbers.precision-limit",
                "numbers.divisor-not-separated"
            ]
            .contains(&e.code)
        ),
        other => panic!("unexpected membership {other:?}"),
    }
}
fn ordered<T: std::fmt::Debug + PartialEq>(
    out: Outcome<T>,
    expected: Option<T>,
    limited: bool,
    bound: Rational,
) {
    match out {
        Outcome::Proven(v) => assert_eq!(Some(v), expected),
        Outcome::ProvenCoincident(_) => assert!(expected.is_none()),
        Outcome::NotExact(x) if limited => {
            assert!(["numbers.budget-exhausted", "numbers.precision-limit"].contains(&x.code));
            assert!(x.bound >= bound);
        }
        Outcome::Refused(e) if limited => {
            assert!(
                [
                    "numbers.deadline",
                    "numbers.budget-exhausted",
                    "numbers.precision-limit",
                    "numbers.divisor-not-separated"
                ]
                .contains(&e.code)
            )
        }
        other => panic!("unexpected ordering {other:?}"),
    }
}
fn queries(mut dec: impl FnMut() -> Decider, limited: bool) {
    let mut g = GeomStore::new();
    let k = g.add_curve(k5());
    let l = g.add_curve(line28());
    let on = explicit(&mut g, pi([3, 4, 0]));
    let off = explicit(&mut g, p(["3", "4", "1e-30"]));
    membership(on_curve(&g, k, on, &mut dec()), true, limited, n(0));
    membership(on_curve(&g, k, off, &mut dec()), false, limited, r("1e-30"));
    let upper = upper26().add(&mut g);
    let opposite = explicit(&mut g, pi([1, 0, -3]));
    membership(
        on_curve(&g, upper, opposite, &mut dec()),
        false,
        limited,
        n(0),
    );
    let top = viviani(true).add(&mut g);
    let bottom = viviani(false).add(&mut g);
    let node = explicit(&mut g, pi([-1, 0, 0]));
    for c in [top, bottom] {
        membership(on_curve(&g, c, node, &mut dec()), true, limited, n(0));
    }
    let same = g.add_curve(Curve::Circle(
        Circle::new(pi([0, 0, 0]), vi([0, 0, -5]), n(25), vi([0, 1, 0])).unwrap(),
    ));
    let different = g.add_curve(Curve::Circle(
        Circle::new(
            pi([0, 0, 0]),
            vi([0, 0, 1]),
            n(25) + r("1e-30"),
            vi([1, 0, 0]),
        )
        .unwrap(),
    ));
    for (other, equal) in [(same, true), (different, false), (l, false), (k, true)] {
        membership(same_curve(&g, k, other, &mut dec()), equal, limited, n(0));
        membership(same_curve(&g, other, k, &mut dec()), equal, limited, n(0));
    }
    let a = explicit(&mut g, pi([-5, 0, 0]));
    let b = explicit(&mut g, pi([5, 0, 0]));
    ordered(
        order_on_curve(&g, l, a, b, &mut dec()),
        Some(StrictOrder::Less),
        limited,
        n(10),
    );
    ordered(
        order_on_curve(&g, l, b, a, &mut dec()),
        Some(StrictOrder::Greater),
        limited,
        n(10),
    );
    ordered(order_on_curve(&g, l, a, a, &mut dec()), None, limited, n(0));
    let p = explicit(&mut g, pi([5, 0, 0]));
    let q = explicit(&mut g, pi([0, 5, 0]));
    let r = explicit(&mut g, pi([-5, 0, 0]));
    for points in [[p, q, r], [q, r, p], [r, p, q]] {
        ordered(
            cyclic_order(&g, k, points, &mut dec()),
            Some(Cyclic::Forward),
            limited,
            n(50),
        );
    }
    for points in [[p, r, q], [q, p, r], [r, q, p]] {
        ordered(
            cyclic_order(&g, k, points, &mut dec()),
            Some(Cyclic::Backward),
            limited,
            n(50),
        );
    }
    for points in [[p, p, q], [p, q, p], [p, q, q]] {
        ordered(cyclic_order(&g, k, points, &mut dec()), None, limited, n(0));
    }
    let Curve::Intersection(branch) = g.curve(top) else {
        unreachable!()
    };
    let w = branch.witness;
    let supports = branch.surfaces;
    let y = rr(3).sqrt() / rr(2);
    let q = closed(
        &mut g,
        Point3::new(-rr(1) / rr(2), -y.clone(), rr(1)),
        &supports,
    );
    let r = closed(&mut g, Point3::new(-rr(1) / rr(2), y, rr(1)), &supports);
    ordered(
        cyclic_order(&g, top, [w, q, r], &mut dec()),
        Some(Cyclic::Forward),
        limited,
        n(0),
    );
    for (c, p, expected) in [
        (k, on, vi([-4, 3, 0]).to_real()),
        (l, a, vi([1, 0, 0]).to_real()),
    ] {
        let actual = tangent_at(&g, c, p, &mut dec()).unwrap();
        for v in actual {
            vector_eq(&v, &expected);
        }
    }
}
blind_test!(cases_34_to_38_predicates_and_ordering_identities, {
    queries(Decider::default, false);
});
blind_test!(case_47_small_budgets_never_prove_wrong_queries, {
    queries(tiny, true);
    queries(starved, true);
    queries(expired, true);
});

blind_test!(case_35_closed_forms_and_alternate_witnesses, {
    let mut g = GeomStore::new();
    let a = g.add_surface(zc(5));
    let b = g.add_surface(sphere(pi([0, 0, 1]), n(6)));
    let h = rr(11).sqrt();
    let mut cs = vec![];
    for (x, z) in [
        (rr(5), &rr(1) + &h),
        (rr(-5), &rr(1) + &h),
        (rr(-5), &rr(1) - &h),
    ] {
        let w = closed(&mut g, Point3::new(x, rr(0), z), &[a, b]);
        cs.push(g.add_curve(Curve::Intersection(IntersectionCurve {
            surfaces: [a, b],
            witness: w,
        })));
    }
    coincident(same_curve(&g, cs[0], cs[1], &mut Decider::default()));
    assert_eq!(
        same_curve(&g, cs[0], cs[2], &mut Decider::default()),
        Outcome::Proven(Apart)
    );
    for (c, sign) in [(cs[0], 1), (cs[2], -1)] {
        let CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ref_dir,
        } = curve_form(&g, c, &mut Decider::default()).unwrap()
        else {
            panic!("circle form required")
        };
        point_eq(
            &center,
            &Point3::new(rr(0), rr(0), rr(1) + (&rr(sign) * &h)),
        );
        assert_eq!(normal, vi([0, 0, -sign]));
        real_eq(&radius_sq, &rr(25));
        assert_eq!(ref_dir, vi([1, 0, 0]));
    }
    let plane = g.add_surface(axis_plane(0, n(1)));
    let mut lines = vec![];
    for y in [-rr(2) * rr(6).sqrt(), rr(2) * rr(6).sqrt()] {
        let w = closed(&mut g, Point3::new(rr(1), y, rr(0)), &[plane, a]);
        lines.push(g.add_curve(Curve::Intersection(IntersectionCurve {
            surfaces: [plane, a],
            witness: w,
        })));
    }
    assert_eq!(
        same_curve(&g, lines[0], lines[1], &mut Decider::default()),
        Outcome::Proven(Apart)
    );
    let first = g.add_curve(line28());
    let same = g.add_curve(Curve::Line(
        Line::new(pi([7, 0, 0]), vi([-3, 0, 0])).unwrap(),
    ));
    coincident(same_curve(&g, first, same, &mut Decider::default()));
    let e = g.add_curve(ellipse4());
    let Curve::Ellipse(mut other) = ellipse4() else {
        unreachable!()
    };
    other.normal = other.normal.scaled(&n(-3));
    other.major_dir = other.major_dir.scaled(&n(-7));
    let same = g.add_curve(Curve::Ellipse(other));
    coincident(same_curve(&g, e, same, &mut Decider::default()));
});

blind_test!(cases_17_19_24_38_regular_and_one_sided_tangents, {
    let mut g = GeomStore::new();
    let e = g.add_curve(ellipse4());
    for (p, expected) in [
        (p(["0", "-5", "15/4"]), v(["4/5", "0", "0"])),
        (pi([5, 0, 0]), v(["0", "4/5", "-3/5"])),
    ] {
        let p = explicit(&mut g, p);
        for tangent in tangent_at(&g, e, p, &mut expired()).unwrap() {
            vector_eq(&tangent, &expected.to_real());
        }
    }
    let regular = Carrier::Branch {
        a: zc(2),
        b: Box::new(cyl(pi([0, 0, 0]), vi([1, 0, 0]), n(1))),
        witness: pi([2, 0, 1]).to_real(),
    }
    .add(&mut g);
    let Curve::Intersection(b) = g.curve(regular) else {
        unreachable!()
    };
    for t in tangent_at(&g, regular, b.witness, &mut Decider::default()).unwrap() {
        positive_multiple(&t, &vi([0, -1, 0]).to_real());
    }
    let upper = viviani(true).add(&mut g);
    let Curve::Intersection(b) = g.curve(upper) else {
        unreachable!()
    };
    for t in tangent_at(&g, upper, b.witness, &mut Decider::default()).unwrap() {
        positive_multiple(&t, &vi([0, -1, 0]).to_real());
    }
    let node = explicit(&mut g, pi([-1, 0, 0]));
    let t = tangent_at(&g, upper, node, &mut Decider::default()).unwrap();
    positive_multiple(&t[0], &vi([0, 1, -1]).to_real());
    positive_multiple(&t[1], &vi([0, 1, 1]).to_real());
    // At the CC4 node y=2-h, x²=4h-h² and z²=2h-h².
    // The x>0 branch arrives along (-sqrt(2),0,1), leaves along (sqrt(2),0,1).
    let branch = Carrier::Branch {
        a: zc(2),
        b: Box::new(cyl(pi([0, 1, 0]), vi([1, 0, 0]), n(1))),
        witness: Point3::new(rr(3).sqrt(), rr(1), rr(1)),
    }
    .add(&mut g);
    let node = explicit(&mut g, pi([0, 2, 0]));
    let t = tangent_at(&g, branch, node, &mut Decider::default()).unwrap();
    positive_multiple(&t[0], &Vector3::new(-rr(2).sqrt(), rr(0), rr(1)));
    positive_multiple(&t[1], &Vector3::new(rr(2).sqrt(), rr(0), rr(1)));
});

blind_test!(
    case_39_arc_points_are_strictly_inside_without_oracle_facts,
    {
        let mut g = GeomStore::new();
        let k = g.add_curve(k5());
        let l = g.add_curve(line28());
        let branch = viviani(true).add(&mut g);
        let z = g.add_surface(axis_plane(2, n(0)));
        let y = g.add_surface(axis_plane(1, n(0)));
        let cylinder = g.add_surface(zc(5));
        let p = explicit(&mut g, pi([5, 0, 0]));
        let q = explicit(&mut g, pi([0, 5, 0]));
        let a = explicit(&mut g, pi([-5, 0, 0]));
        for (limited, factory) in [
            (false, Decider::default as fn() -> Decider),
            (true, tiny),
            (true, starved),
            (true, expired),
        ] {
            for (c, from, to) in [
                (k, Some(p), Some(q)),
                (l, Some(a), Some(p)),
                (branch, None, None),
                (k, None, None),
                (k, Some(q), Some(q)),
            ] {
                let before = fingerprint(&g);
                let first = arc_point(&g, c, from, to, &mut factory());
                assert_eq!(fingerprint(&g), before);
                let answer = match first {
                    Ok(a) => {
                        assert!(a.marks.is_empty());
                        a.value
                    }
                    Err(e) => {
                        assert!(limited, "default arc_point refused: {e:?}");
                        assert!(
                            [
                                "numbers.deadline",
                                "numbers.budget-exhausted",
                                "numbers.precision-limit",
                                "numbers.divisor-not-separated",
                                "numbers.root-not-certified"
                            ]
                            .contains(&e.code)
                        );
                        continue;
                    }
                };
                if c == k {
                    coincident(g.side_of(z, &answer, &mut Decider::default()));
                    coincident(g.side_of(cylinder, &answer, &mut Decider::default()));
                    if from == Some(p) && to == Some(q) {
                        positive(
                            &answer
                                .minus(&g.position(p))
                                .cross(&g.position(q).minus(&g.position(p)))
                                .dot(&vi([0, 0, 1]).to_real()),
                        );
                    } else {
                        let start = from
                            .map(|p| g.position(p))
                            .unwrap_or_else(|| pi([5, 0, 0]).to_real());
                        positive(&answer.minus(&start).norm_squared());
                    }
                } else if c == l {
                    coincident(g.side_of(y, &answer, &mut Decider::default()));
                    coincident(g.side_of(z, &answer, &mut Decider::default()));
                    positive(&(&answer.x + &rr(5)));
                    positive(&(&rr(5) - &answer.x));
                } else {
                    let Curve::Intersection(b) = g.curve(c) else {
                        unreachable!()
                    };
                    for s in b.surfaces {
                        coincident(g.side_of(s, &answer, &mut Decider::default()));
                    }
                    assert!(matches!(
                        Decider::default().sign(&answer.z),
                        Outcome::Proven(StrictSign::Positive) | Outcome::ProvenCoincident(_)
                    ));
                    positive(&answer.minus(&pi([1, 0, 2]).to_real()).norm_squared());
                }
                let again = exact(arc_point(&g, c, from, to, &mut factory()));
                assert_eq!(format!("{answer:?}"), format!("{again:?}"));
            }
        }
    }
);

blind_test!(
    facts_bypass_expired_deadline_and_branch_uncertainty_has_distance_bound,
    {
        let mut g = GeomStore::new();
        let c = upper26().add(&mut g);
        let Curve::Intersection(b) = g.curve(c) else {
            unreachable!()
        };
        let surfaces = b.surfaces;
        let lower = closed(&mut g, pi([1, 0, -3]).to_real(), &surfaces);
        // Both surface incidences are known, but membership in the upper component is not.
        match on_curve(&g, c, lower, &mut starved()) {
            Outcome::NotExact(e) => {
                assert_eq!(e.code, "numbers.budget-exhausted");
                assert!(e.bound >= n(3), "distance to z=0 is 3");
            }
            other => panic!("branch sign with zero work must be undecided: {other:?}"),
        }
        let node = explicit(&mut g, pi([-1, 0, 0]));
        let c = viviani(true).add(&mut g);
        g.record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Point(node),
            b: GeomKey::Curve(c),
            at: None,
            origin: FactOrigin::Construction("test.exact-node"),
        });
        for which in [false, true] {
            let mut d = expired();
            let before = fingerprint(&g);
            if which {
                coincident(same_curve(&g, c, c, &mut d));
            } else {
                coincident(on_curve(&g, c, node, &mut d));
            }
            assert_eq!(
                (
                    d.stats().decisions,
                    d.stats().coincident_fact,
                    d.stats().work
                ),
                (1, 1, 0)
            );
            assert_eq!(fingerprint(&g), before);
        }
        let a = g.add_curve(k5());
        let b = g.add_curve(k5());
        let s = g.add_surface(axis_plane(2, n(0)));
        g.record_fact(Fact {
            kind: FactKind::Coincident,
            a: GeomKey::Curve(a),
            b: GeomKey::Curve(b),
            at: None,
            origin: FactOrigin::Construction("test.equal-circles"),
        });
        g.record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Curve(a),
            b: GeomKey::Surface(s),
            at: None,
            origin: FactOrigin::Construction("test.circle-plane"),
        });
        let before = fingerprint(&g);
        let mut d = expired();
        coincident(same_curve(&g, a, b, &mut d));
        assert_eq!(
            exact(intersect_curve_surface(&mut g, b, s, &mut d)),
            CurveSurface::Contained
        );
        assert_eq!(
            (
                d.stats().decisions,
                d.stats().coincident_fact,
                d.stats().work
            ),
            (2, 2, 0)
        );
        assert_eq!(fingerprint(&g), before);
    }
);

blind_test!(
    unsupported_intersection_definition_leaves_the_store_unchanged,
    {
        let cc6 = Carrier::Branch {
            a: zc(1),
            b: Box::new(cyl(pi([0, 0, 0]), vi([1, 0, 0]), n(1))),
            witness: Point3::new(root(&r("1/2")), root(&r("1/2")), root(&r("1/2"))),
        };
        let carrier = cc6;
        let mut g = GeomStore::new();
        let c = carrier.add(&mut g);
        let Curve::Intersection(b) = g.curve(c) else {
            unreachable!()
        };
        let p = b.witness;
        let s = b.surfaces[1];
        let good = g.add_curve(line28());
        let before = fingerprint(&g);
        assert_eq!(
            curve_form(&g, c, &mut expired()).unwrap_err().code,
            "intersect.unsupported-curve"
        );
        assert_eq!(
            curve_is_closed(&g, c, &mut expired()).unwrap_err().code,
            "intersect.unsupported-curve"
        );
        for out in [
            on_curve(&g, c, p, &mut expired()),
            same_curve(&g, c, c, &mut expired()),
            same_curve(&g, c, good, &mut expired()),
            same_curve(&g, good, c, &mut expired()),
        ] {
            match out {
                Outcome::Refused(e) => assert_eq!(e.code, "intersect.unsupported-curve"),
                other => panic!("{other:?}"),
            }
        }
        match order_on_curve(&g, c, p, p, &mut expired()) {
            Outcome::Refused(e) => assert_eq!(e.code, "intersect.unsupported-curve"),
            other => panic!("{other:?}"),
        }
        match cyclic_order(&g, c, [p, p, p], &mut expired()) {
            Outcome::Refused(e) => assert_eq!(e.code, "intersect.unsupported-curve"),
            other => panic!("{other:?}"),
        }
        assert_eq!(
            tangent_at(&g, c, p, &mut expired()).unwrap_err().code,
            "intersect.unsupported-curve"
        );
        assert_eq!(
            arc_point(&g, c, Some(p), Some(p), &mut expired())
                .unwrap_err()
                .code,
            "intersect.unsupported-curve"
        );
        assert_eq!(
            trace(&g, c, Some(p), Some(p), &mut expired())
                .unwrap_err()
                .code,
            "intersect.unsupported-curve"
        );
        assert_eq!(
            intersect_curve_surface(&mut g, c, s, &mut expired())
                .unwrap_err()
                .code,
            "intersect.unsupported-curve"
        );
        assert_eq!(fingerprint(&g), before);
    }
);
