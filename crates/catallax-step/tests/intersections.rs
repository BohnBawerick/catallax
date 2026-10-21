//! Analytic intersection preservation and exact spline sample checks.
mod common;
use common::*;

blind_test!(
    #[ignore = "requires catallax-boolean: boolean() panics at todo!() before STEP export"]
    case16_irrational_intersection_lines_remain_lines,
    {
        let (m, b) = boolean_fixture(13);
        let (_, a, s) = export(&m, &[b], &options(17));
        report(&a, [1, 4, 6, 4, 0, 0, 0], false);
        counts(
            &s,
            &[("LINE", 4), ("CIRCLE", 2), ("B_SPLINE_CURVE_WITH_KNOTS", 0)],
        );
        assert!(
            m.edges_of_body(b)
                .iter()
                .any(|&e| matches!(m.geom().curve(m.edge(e).curve), Curve::Intersection(_)))
        );
        // PC3: (y-5)^2=4e-30-e-60 at z=10. This is non-square rational.
        let expected_sq = r("4e-30") - r("1e-60");
        let delta = &a.value.max_deviation;
        for id in s.ids("VERTEX_POINT") {
            let p = s.xyz(s.vertex(id));
            assert!((&p[2] - &n(10)).abs() <= *delta);
            // Written distance to exact two section lines, allowing decimal error.
            let dy = (&p[1] - &n(5)).abs();
            assert!(expected_sq <= sq(&(&dy + delta)));
            if dy > *delta {
                assert!(expected_sq >= sq(&(&dy - delta)));
            }
        }
    }
);
blind_test!(
    #[ignore = "requires catallax-boolean: boolean() panics at todo!() before STEP export"]
    case16_irrational_circle_centres_remain_circles,
    {
        let (m, b) = boolean_fixture(30);
        let (_, a, s) = export(&m, &[b], &options(17));
        report(&a, [1, 3, 2, 2, 2, 0, 0], false);
        counts(
            &s,
            &[
                ("CIRCLE", 2),
                ("B_SPLINE_CURVE_WITH_KNOTS", 0),
                ("SPHERICAL_SURFACE", 1),
                ("CYLINDRICAL_SURFACE", 1),
            ],
        );
        // r_sphere^2-r_cylinder^2=36-9=27, centres z=+/-sqrt(27).
        let d = &a.value.max_deviation;
        let mut signs = BTreeSet::new();
        for id in s.ids("CIRCLE") {
            let c = &s.entities[&id];
            rounded(&c.args[2], &n(3), 17);
            let (p, _, _) = s.placement(reference(&c.args[1]));
            assert_eq!(&p[..2], [n(0), n(0)]);
            signs.insert(p[2] > n(0));
            let h = p[2].abs();
            assert!(n(27) <= sq(&(&h + d)));
            if h > *d {
                assert!(n(27) >= sq(&(&h - d)));
            }
        }
        assert_eq!(signs.len(), 2);
    }
);
blind_test!(
    #[ignore = "requires catallax-boolean: boolean() panics at todo!() before STEP export"]
    analytic_ellipses_shared_by_distinct_arcs,
    {
        let (m, b) = boolean_fixture(24);
        let (_, a, s) = export(&m, &[b], &options(17));
        report(&a, [1, 4, 4, 2, 0, 0, 0], false);
        counts(&s, &[("ELLIPSE", 2), ("B_SPLINE_CURVE_WITH_KNOTS", 0)]);
        // Each ellipse has semiaxes sqrt(2),1, because x=+/-z on x²+y²=1.
        for id in s.ids("ELLIPSE") {
            let e = &s.entities[&id];
            assert!(["1.414213562373095", "1.4142135623730951"].contains(&e.args[2].as_str()));
            assert_eq!(e.args[3], "1.");
        }
    }
);
blind_test!(
    #[ignore = "requires catallax-boolean: boolean() panics at todo!() before STEP export"]
    case17_branch_splines_at_default_and_tight_tolerance,
    {
        for case in [23, 33] {
            let (m, b) = boolean_fixture(case);
            for tolerance in ["1e-7", "1e-12"] {
                let mut o = options(17);
                o.curve_tolerance = r(tolerance);
                let (_, a, s) = export(&m, &[b], &o);
                let expected = if case == 23 { 2 } else { 1 };
                assert_eq!(a.value.approximated_curves, expected);
                assert_eq!(a.value.added_vertices, expected);
                assert_eq!(a.value.added_edges, 0);
                assert_eq!(s.count("CIRCLE"), 0);
                if tolerance == "1e-7" {
                    assert!(a.value.max_deviation <= o.curve_tolerance);
                    assert!(a.marks.is_empty());
                }
                if case == 23 {
                    check_splines(&s, &a.value.max_deviation, crossing_cylinders);
                } else {
                    check_splines(&s, &a.value.max_deviation, offset_sphere);
                }
            }
        }
    }
);
blind_test!(
    #[ignore = "requires catallax-boolean: boolean() panics at todo!() before STEP export"]
    case18_budget_limited_fit_still_has_proven_bound,
    {
        let (m, b) = boolean_fixture(23);
        let mut o = options(17);
        o.curve_tolerance = r("1e-12");
        let mut dec = Decider::new(BudgetPolicy {
            max_precision_bits: 64,
            max_work: 10_000,
        });
        assert_eq!(m.check(b), vec![], "invalid test fixture");
        let before = b"existing output".to_vec();
        let mut out = before.clone();
        match write_step(&m, &[b], &o, &mut dec, &mut out) {
            Ok(answer) => {
                assert_eq!(&out[..before.len()], before);
                let bytes = out.split_off(before.len());
                let (_, a, s) = check_export(&m, &[b], &o, bytes, answer);
                assert_eq!(a.value.approximated_curves, 2);
                check_splines(&s, &a.value.max_deviation, crossing_cylinders);
            }
            Err(refusal) => assert_budget_refusal(refusal.code, &before, &out),
        }
    }
);
blind_test!(
    #[ignore = "requires catallax-boolean: boolean() panics at todo!() before STEP export"]
    branches_through_nodes_keep_single_topological_vertex,
    {
        for case in [28, 32] {
            let (m, b) = boolean_fixture(case);
            let (_, a, s) = export(&m, &[b], &options(17));
            report(&a, [1, 3, 2, 1, 0, 0, 2], false);
            let expected = if case == 28 {
                [n(0), n(2), n(0)]
            } else {
                [n(-1), n(0), n(0)]
            };
            assert_eq!(vertex_set(&s), BTreeSet::from([expected]));
            check_splines(&s, &a.value.max_deviation, |p, d| {
                if case == 28 {
                    radial_distance(sq(&p[0]) + sq(&p[1]), &n(2), d);
                    radial_distance(sq(&(&p[1] - &n(1))) + sq(&p[2]), &n(1), d);
                } else {
                    radial_distance(sq(&p[0]) + sq(&p[1]), &n(1), d);
                    radial_distance(sq(&(&p[0] - &n(1))) + sq(&p[1]) + sq(&p[2]), &n(2), d);
                }
            });
        }
    }
);
blind_test!(
    #[ignore = "requires catallax-boolean: boolean() panics at todo!() before STEP export"]
    case19_branch_fits_repeat_on_eight_threads,
    {
        let (m, b) = boolean_fixture(23);
        let o = options(17);
        let (bytes, a, _) = export(&m, &[b], &o);
        let threads = (0..8)
            .map(|_| {
                let m = m.clone();
                let o = o.clone();
                std::thread::spawn(move || {
                    let (b, a, _) = export(&m, &[b], &o);
                    (b, a)
                })
            })
            .collect::<Vec<_>>();
        for t in threads {
            let (other, result) = t.join().unwrap();
            assert_eq!(other, bytes);
            assert_eq!(result.value, a.value);
            assert_eq!(result.marks, a.marks);
        }
    }
);
