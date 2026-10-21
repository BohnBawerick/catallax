//! Additional carrier forms, exact identities and eight-thread query replay.
mod common;
use common::curves::*;
use common::*;

fn samples() -> Vec<(Carrier, [Point3; 3], bool)> {
    let y = rr(3).sqrt() / rr(2);
    let z = &rr(1) + &rr(11).sqrt();
    let line_y = -rr(2) * rr(6).sqrt();
    vec![
        (
            Carrier::Analytic(line28()),
            [
                pi([-5, 0, 0]).to_real(),
                pi([5, 0, 0]).to_real(),
                pi([0, 0, 0]).to_real(),
            ],
            true,
        ),
        (
            Carrier::Analytic(k5()),
            [
                pi([5, 0, 0]).to_real(),
                pi([0, 5, 0]).to_real(),
                pi([-5, 0, 0]).to_real(),
            ],
            false,
        ),
        (
            Carrier::Analytic(ellipse4()),
            [
                p(["0", "-5", "15/4"]).to_real(),
                pi([5, 0, 0]).to_real(),
                p(["0", "5", "-15/4"]).to_real(),
            ],
            false,
        ),
        (
            viviani(true),
            [
                pi([1, 0, 2]).to_real(),
                Point3::new(-rr(1) / rr(2), -y.clone(), rr(1)),
                Point3::new(-rr(1) / rr(2), y, rr(1)),
            ],
            false,
        ),
        (
            loop25(),
            [
                Point3::new(rr(1), rr(0), rr(3).sqrt()),
                Point3::new(rr(1) / rr(4), rr(15).sqrt() / rr(4), rr(0)),
                Point3::new(rr(1), rr(0), -rr(3).sqrt()),
            ],
            false,
        ),
        (
            Carrier::Branch {
                a: axis_plane(0, n(1)),
                b: Box::new(zc(5)),
                witness: Point3::new(rr(1), line_y.clone(), rr(0)),
            },
            [
                Point3::new(rr(1), line_y.clone(), rr(0)),
                Point3::new(rr(1), line_y.clone(), rr(-1)),
                Point3::new(rr(1), line_y, rr(-2)),
            ],
            true,
        ),
        (
            Carrier::Branch {
                a: zc(5),
                b: Box::new(sphere(pi([0, 0, 1]), n(6))),
                witness: Point3::new(rr(5), rr(0), z.clone()),
            },
            [
                Point3::new(rr(5), rr(0), z.clone()),
                Point3::new(rr(0), rr(-5), z.clone()),
                Point3::new(rr(-5), rr(0), z),
            ],
            false,
        ),
    ]
}
blind_test!(case_46_every_query_replays_identically_on_eight_threads, {
    for (carrier, positions, line) in samples() {
        let mut base = GeomStore::new();
        let c = carrier.add(&mut base);
        let supports = match base.curve(c) {
            Curve::Intersection(b) => Some(b.surfaces),
            _ => None,
        };
        let points = positions.map(|p| match supports {
            Some(ss) => closed(&mut base, p, &ss),
            None => explicit(&mut base, rational_point(&p).unwrap()),
        });
        assert_eq!(
            curve_is_closed(&base, c, &mut Decider::default()).unwrap(),
            !line
        );
        let results = std::thread::scope(|scope| {
            (0..8)
                .map(|_| {
                    let g = base.clone();
                    scope.spawn(move || {
                        let before = fingerprint(&g);
                        let mut d = Decider::default();
                        let [p, q, r] = points;
                        let mut values = vec![];
                        values.push(format!("{:?}", curve_form(&g, c, &mut d)));
                        values.push(format!("{:?}", on_curve(&g, c, q, &mut d)));
                        values.push(format!("{:?}", same_curve(&g, c, c, &mut d)));
                        values.push(if line {
                            format!("{:?}", order_on_curve(&g, c, p, q, &mut d))
                        } else {
                            format!("{:?}", cyclic_order(&g, c, [p, q, r], &mut d))
                        });
                        values.push(format!("{:?}", tangent_at(&g, c, q, &mut d)));
                        values.push(format!("{:?}", arc_point(&g, c, Some(p), Some(q), &mut d)));
                        values.push(format!("{:?}", trace(&g, c, Some(p), Some(q), &mut d)));
                        assert_eq!(fingerprint(&g), before);
                        (values, d.stats().clone())
                    })
                })
                .collect::<Vec<_>>()
                .into_iter()
                .map(|h| h.join().unwrap())
                .collect::<Vec<_>>()
        });
        for r in &results[1..] {
            assert_eq!(r, &results[0]);
        }
    }
});
blind_test!(handmade_plane_pair_and_budgeted_component_form, {
    let mut g = GeomStore::new();
    let a = g.add_surface(axis_plane(2, n(1)));
    let b = g.add_surface(axis_plane(0, n(3)));
    let witness = closed(&mut g, pi([3, 7, 1]).to_real(), &[a, b]);
    let c = g.add_curve(Curve::Intersection(IntersectionCurve {
        surfaces: [a, b],
        witness,
    }));
    let mut d = expired();
    let CurveForm::Line { origin, dir } = curve_form(&g, c, &mut d).unwrap() else {
        panic!("line form required")
    };
    point_eq(&origin, &pi([3, 7, 1]).to_real());
    assert_eq!(dir, vi([0, 1, 0]));
    assert_eq!(d.stats().work, 0);
    let a = g.add_surface(axis_plane(0, n(1)));
    let b = g.add_surface(zc(5));
    let pos = Point3::new(rr(1), -rr(2) * rr(6).sqrt(), rr(7));
    let witness = closed(&mut g, pos.clone(), &[a, b]);
    let c = g.add_curve(Curve::Intersection(IntersectionCurve {
        surfaces: [a, b],
        witness,
    }));
    let mut d = Decider::default();
    let CurveForm::Line { origin, dir } = curve_form(&g, c, &mut d).unwrap() else {
        panic!("line form required")
    };
    point_eq(&origin, &pos);
    assert_eq!(dir, vi([0, 0, -1]));
    assert_eq!(d.stats().decisions, 1);
    assert_eq!(
        curve_form(&g, c, &mut starved()).unwrap_err().code,
        "numbers.budget-exhausted"
    );
    assert_eq!(
        curve_form(&g, c, &mut expired()).unwrap_err().code,
        "numbers.deadline"
    );
});
blind_test!(deadline_refusals_remain_deadlines_not_budget_marks, {
    let mut g = GeomStore::new();
    let l = g.add_curve(line28());
    let k = g.add_curve(k5());
    let p = explicit(&mut g, pi([5, 0, 0]));
    let q = explicit(&mut g, pi([0, 5, 0]));
    let r = explicit(&mut g, pi([-5, 0, 0]));
    match order_on_curve(&g, l, r, p, &mut expired()) {
        Outcome::Refused(e) => assert_eq!(e.code, "numbers.deadline"),
        other => panic!("{other:?}"),
    }
    match cyclic_order(&g, k, [p, q, r], &mut expired()) {
        Outcome::Refused(e) => assert_eq!(e.code, "numbers.deadline"),
        other => panic!("{other:?}"),
    }
    for outcome in [
        on_curve(&g, k, p, &mut expired()),
        same_curve(&g, k, k, &mut expired()),
    ] {
        match outcome {
            Outcome::ProvenCoincident(_) => {}
            Outcome::Refused(e) => assert_eq!(e.code, "numbers.deadline"),
            other => panic!("{other:?}"),
        }
    }
    match arc_point(&g, k, Some(p), Some(q), &mut expired()) {
        Ok(a) => assert!(a.marks.is_empty()),
        Err(e) => assert_eq!(e.code, "numbers.deadline"),
    }
    match trace(&g, k, Some(p), Some(q), &mut expired()) {
        Ok(a) => assert!(a.marks.is_empty()),
        Err(e) => assert_eq!(e.code, "numbers.deadline"),
    }
});
