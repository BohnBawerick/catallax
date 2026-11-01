#![allow(missing_docs)]
mod common;
use catallax_numbers::*;
use common::*;
use rug::ops::Pow;
use rug::{Integer as Z, Rational as Q};
use std::sync::Arc;

blind_test!(polynomial_normal_form_arithmetic_and_factorisation, {
    assert_eq!(poly(&[0, 0]).degree(), None);
    assert!(UPoly::constant(r("0")).coefficients().is_empty());
    assert_eq!(UPoly::from_roots(&[]), poly(&[1]));
    let p = poly(&[-6, 11, -6, 1]);
    assert_eq!(p, UPoly::from_roots(&[r("1"), r("2"), r("3")]));
    assert_eq!(p.derivative(), poly(&[11, -12, 3]));
    let (quotient, remainder) = p.div_rem(&poly(&[-1, 1])).unwrap();
    assert_eq!(quotient, poly(&[6, -5, 1]));
    assert!(remainder.is_zero());
    assert!(p.div_rem(&UPoly::zero()).is_none());
    assert_eq!(p.compose(&poly(&[1, 1])), poly(&[0, 2, -3, 1]));
    assert_eq!(
        UPoly::new(vec![r("-1/2"), r("1/3"), r("-1/6")]).primitive_integer(),
        vec![Integer::from(3), Integer::from(-2), Integer::from(1)]
    );
    assert!(UPoly::zero().primitive_integer().is_empty());
    assert!(UPoly::zero().gcd(&UPoly::zero()).is_zero());
    let a = poly(&[-1, 1]);
    let b = poly(&[2, 1]);
    let c = poly(&[1, 0, 1]);
    let repeated = &(&(&a * &a) * &a) * &(&b * &b);
    let full = &repeated * &c;
    assert_eq!(
        full.square_free_factors(),
        vec![(c, 1), (b.clone(), 2), (a.clone(), 3)]
    );
    assert_eq!(full.gcd(&full.derivative()), &(&a * &a) * &b);
    assert!(poly(&[5]).square_free_factors().is_empty());
    assert!(UPoly::zero().square_free_factors().is_empty());
    let mut rng = Rng(0xa01);
    for _ in 0..256 {
        let x = rng.rational();
        let expected = ((x.clone() - 1) * (x.clone() - 2)) * (x.clone() - 3);
        let input = eq(&x);
        assert_q(&p.eval(&input), &expected);
        interval_contains(p.eval_interval(Interval::from_rational(&input)), &expected);
        ball_contains(&p.eval_ball(&Ball::from_rational(&input, 128)), &expected);
        assert_q(
            &p.eval_real(&Real::from_rational(input))
                .as_rational()
                .unwrap(),
            &expected,
        );
    }
});

blind_test!(gcd_preserves_rational_factors_at_large_integer_scales, {
    let prime = r("2147483647");
    let common = UPoly::new(vec![r("1/2147483647"), r("1")]);
    let scale = UPoly::constant(prime.clone());
    let left = &(&common * &poly(&[-3, 1])) * &scale;
    let right = &(&common * &poly(&[5, 1])) * &scale;
    // Reducing modulo this integer loses the common factor's degree.
    assert_eq!(left.gcd(&right), common);

    let common = UPoly::new(vec![r("-1/2147483647"), r("1")]);
    let left = &common * &poly(&[-3, 1]);
    let right = &common * &poly(&[5, 1]);
    // Rational denominators need not be invertible in a chosen finite field.
    assert_eq!(left.gcd(&right), common);

    let left = poly(&[1, 0, 1]);
    let right = &left + &UPoly::constant(prime);
    // A shared finite-field image does not imply a shared rational factor.
    assert_eq!(left.gcd(&right), poly(&[1]));
});

blind_test!(
    multivariate_order_substitution_derivatives_and_resultants,
    {
        let p = MPoly::from_terms(
            2,
            vec![
                (vec![2, 0], r("1")),
                (vec![0, 1], r("2")),
                (vec![1, 0], r("4")),
                (vec![0, 0], r("3")),
                (vec![1, 0], r("-4")),
                (vec![1, 1], r("5")),
            ],
        );
        let terms = p.terms();
        let exponents: Vec<_> = terms.iter().map(|(m, _)| m.exponents().to_vec()).collect();
        assert_eq!(
            exponents,
            vec![vec![0, 0], vec![0, 1], vec![1, 1], vec![2, 0]]
        );
        let mut monomials = [
            Monomial::new(vec![2, 0]),
            Monomial::new(vec![0, 1]),
            Monomial::new(vec![1, 1]),
            Monomial::new(vec![0, 0]),
        ];
        monomials.sort();
        assert_eq!(
            monomials
                .iter()
                .map(|m| m.exponents().to_vec())
                .collect::<Vec<_>>(),
            exponents
        );
        assert_eq!(p.total_degree(), Some(2));
        assert!(MPoly::zero(2).total_degree().is_none());
        let x = MPoly::var(2, 0);
        let y = MPoly::var(2, 1);
        assert_eq!(p.partial(0), &x.scale(&r("2")) + &y.scale(&r("5")));
        assert_eq!(
            p.partial(1),
            &MPoly::constant(2, r("2")) + &x.scale(&r("5"))
        );
        assert!(p.to_upoly(0).is_none());
        assert_eq!(x.to_upoly(0), Some(poly(&[0, 1])));
        let substituted = p.substitute(0, &y);
        assert_eq!(substituted.to_upoly(1), Some(poly(&[3, 2, 6])));
        // Res_x(x^2-y,x-2)=4-y; Res_x(x-y,x-2)=y-2.
        let quadratic = &(&x * &x) - &y;
        let linear = &x - &MPoly::constant(2, r("2"));
        assert_eq!(
            quadratic.resultant(&linear, 0),
            &MPoly::constant(2, r("4")) - &y
        );
        assert_eq!(
            (&x - &y).resultant(&linear, 0),
            &y - &MPoly::constant(2, r("2"))
        );
        assert_eq!(
            linear.resultant(&(&x - &y), 0),
            &MPoly::constant(2, r("2")) - &y
        );
        assert_eq!(
            MPoly::constant(2, r("3")).resultant(&quadratic, 0),
            MPoly::constant(2, r("9"))
        );
        assert!(quadratic.resultant(&MPoly::zero(2), 0).is_zero());
        let point = [r("2"), r("3")];
        assert_q(&p.eval(&point), &q(43, 1));
        interval_contains(
            p.eval_interval(&point.each_ref().map(Interval::from_rational)),
            &q(43, 1),
        );
        ball_contains(
            &p.eval_ball(&point.each_ref().map(|x| Ball::from_rational(x, 128))),
            &q(43, 1),
        );
        assert_q(
            &p.eval_real(&point.map(Real::from_rational))
                .as_rational()
                .unwrap(),
            &q(43, 1),
        );
    }
);

fn verify_rational_roots(p: &UPoly, expected: &[(Q, u32)]) {
    let mut decider = Decider::default();
    let got = roots(&mut decider, p);
    assert_eq!(got.len(), expected.len());
    assert_eq!(decider.stats().decisions, 1);
    assert_eq!(decider.stats().exact, 1);
    stats_sum(decider.stats());
    for (root, (value, multiplicity)) in got.iter().zip(expected) {
        assert_eq!(root.multiplicity, *multiplicity);
        assert_q(
            root.value
                .as_rational()
                .expect("rational root must be exact"),
            value,
        );
        let (lo, hi) = root.value.isolating_interval();
        assert_q(&lo, value);
        assert_q(&hi, value);
        let defining = root.value.defining_polynomial();
        let want = vec![-value.numer().clone(), value.denom().clone()];
        assert_eq!(
            defining.coefficients().iter().map(oq).collect::<Vec<_>>(),
            want.into_iter().map(Q::from).collect::<Vec<_>>()
        );
        ball_contains(&root.value.to_ball(128), value);
    }
}

blind_test!(
    rational_roots_wilkinson_multiplicity_and_generated_families,
    {
        let values: Vec<_> = (1..=20).map(|i| q(i, 1)).collect();
        verify_rational_roots(
            &UPoly::from_roots(&values.iter().map(eq).collect::<Vec<_>>()),
            &values.into_iter().map(|v| (v, 1)).collect::<Vec<_>>(),
        );
        let repeated = UPoly::from_roots(&[r("1"), r("1"), r("1"), r("-2"), r("-2")]);
        verify_rational_roots(
            &(&repeated * &poly(&[1, 0, 1])),
            &[(q(-2, 1), 2), (q(1, 1), 3)],
        );
        let mut rng = Rng(0xa007);
        for _ in 0..128 {
            let mut values: Vec<_> = (0..8).map(|_| rng.rational()).collect();
            values.sort();
            let mut expected: Vec<(Q, u32)> = Vec::new();
            for v in &values {
                if let Some((last, m)) = expected.last_mut()
                    && last == v
                {
                    *m += 1;
                    continue;
                }
                expected.push((v.clone(), 1));
            }
            // x^2+c has no real root for c>0.
            let c = (rng.next() % 20 + 1) as i64;
            let p =
                &UPoly::from_roots(&values.iter().map(eq).collect::<Vec<_>>()) * &poly(&[c, 0, 1]);
            verify_rational_roots(&p, &expected);
        }
        verify_rational_roots(&poly(&[5]), &[]);
        refused(
            Decider::default().real_roots(&UPoly::zero()),
            "numbers.zero-polynomial",
        );
        refused(
            Decider::new(BudgetPolicy {
                max_precision_bits: 1,
                max_work: 0,
            })
            .real_roots(&mignotte(21)),
            "numbers.budget-exhausted",
        );
    }
);

blind_test!(irrational_roots_isolation_refinement_and_real_equality, {
    let mut d = Decider::default();
    let got = roots(&mut d, &poly(&[-2, 0, 1]));
    assert_eq!(got.len(), 2);
    for (i, root) in got.iter().enumerate() {
        assert_eq!(root.multiplicity, 1);
        assert!(root.value.as_rational().is_none());
        let refined = root.value.refined(&r("1e-60"));
        let (lo, hi) = refined.isolating_interval();
        let (a, b) = (oq(&lo), oq(&hi));
        assert!(a < b && b.clone() - &a <= Q::from((Z::from(1), Z::from(10).pow(60))));
        check_irrational_certificate(&root.value);
        check_irrational_certificate(&refined);
        let fa = a.clone() * &a - 2;
        let fb = b.clone() * &b - 2;
        assert!(fa * fb < 0);
        assert!(if i == 0 { b < 0 } else { a > 0 });
        let expected = if i == 0 {
            -real("2").sqrt()
        } else {
            real("2").sqrt()
        };
        assert!(
            d.compare(&Real::from_algebraic(refined.clone()), &expected)
                .is_coincident()
        );
        let oracle = Oracle::rational(&q(2, 1), 1024).sqrt().unwrap();
        let oracle = if i == 0 { oracle.neg() } else { oracle };
        oracle.assert_ball(&refined.to_ball(128));
    }
    let positive = Real::from_algebraic(got[1].value.clone());
    assert_eq!(
        d.compare(&positive, &real("1.41421356")),
        Outcome::Proven(StrictOrder::Greater)
    );
    for k in -30i32..=30 {
        let s = Real::from_rational(eq(&Q::from(10).pow(k)));
        assert!(
            d.compare(&(&s * &positive), &(&s * &real("2").sqrt()))
                .is_coincident()
        );
        assert_eq!(
            d.compare(&(&s * &positive), &(&s * &real("1.41421356"))),
            Outcome::Proven(StrictOrder::Greater)
        );
    }
});

blind_test!(mignotte_close_roots_are_distinct_and_complete, {
    // Three positive roots: signs alternate at 0, 1/10, 1, and +infinity;
    // Descartes' three sign changes forbid any more. For even n, p(-x)
    // has one sign change and opposite signs at 0 and -infinity. For odd
    // n every term of p(-x) is negative, so there is no negative root.
    for n in [20usize, 21] {
        let got = roots(&mut Decider::default(), &mignotte(n));
        assert_eq!(got.len(), if n % 2 == 0 { 4 } else { 3 });
        let mut boxes = Vec::new();
        for root in got {
            assert_eq!(root.multiplicity, 1);
            let (a, b) = root.value.refined(&r("1e-25")).isolating_interval();
            let (a, b) = (oq(&a), oq(&b));
            check_irrational_certificate(&root.value);
            let eval =
                |x: &Q| x.clone().pow(n as i32) - q(2, 1) * (x.clone() * 10i32 - 1i32).square();
            assert!(a < b);
            assert!(eval(&a) * eval(&b) < 0);
            if let Some((_, last)) = boxes.last() {
                assert!(last < &a);
            }
            boxes.push((a, b));
        }
        let offset = if n % 2 == 0 { 1 } else { 0 };
        let (left, right) = (&boxes[offset], &boxes[offset + 1]);
        assert!(left.1 < q(1, 10));
        assert!(right.0 > q(1, 10));
        // Squared form of the specified bound avoids any approximate sqrt:
        // gap < 2*10^(-(n+2)/2) iff gap^2 < 4*10^(-(n+2)).
        let outer_gap = right.1.clone() - &left.0;
        assert!(outer_gap.square() < Q::from((Z::from(4), Z::from(10).pow((n + 2) as u32))));
    }
});

blind_test!(chebyshev_roots_one_through_sixty_four, {
    for n in 1..=64 {
        let got = roots(&mut Decider::default(), &chebyshev(n));
        assert_eq!(got.len(), n);
        let mut previous = None;
        for (i, root) in got.iter().enumerate() {
            assert_eq!(root.multiplicity, 1);
            let value = root.value.refined(&r("1e-30"));
            let (lo, hi) = value.isolating_interval();
            let (lo, hi) = (oq(&lo), oq(&hi));
            assert!(lo > -1 && hi < 1);
            if let Some(last) = previous {
                assert!(last < lo);
            }
            previous = Some(hi.clone());
            // Sorted roots cos((2(n-k)+1)*pi/(2n)), k=1..n. cos decreases on [0,pi].
            let angle = Oracle::pi(1024).mul(&Oracle::rational(
                &q((2 * (n - i - 1) + 1) as i64, (2 * n) as i64),
                1024,
            ));
            let expected = angle.unary(rug::Float::cos_round, false);
            assert!(
                lo <= expected.hi.to_rational().unwrap()
                    && hi >= expected.lo.to_rational().unwrap()
            );
            let pi = Ball::pi(256);
            let angle = pi
                * Ball::from_rational(&eq(&q((2 * (n - i - 1) + 1) as i64, (2 * n) as i64)), 256);
            overlap(&value.to_ball(256), &angle.cos());
        }
    }
});

fn circle_system(tangent: bool) -> Arc<PolySystem> {
    // x=y on the circle gives x=y=+/-sqrt(1/2). The positive box selects
    // the plus sign and det(J)=-2(x+y) is nonzero. With y=1 instead,
    // x^2=0 is a double root and det(J)=2x=0.
    let circle = MPoly::from_terms(
        2,
        vec![
            (vec![2, 0], r("1")),
            (vec![0, 2], r("1")),
            (vec![0, 0], r("-1")),
        ],
    );
    let line = if tangent {
        MPoly::from_terms(2, vec![(vec![0, 1], r("1")), (vec![0, 0], r("-1"))])
    } else {
        MPoly::from_terms(2, vec![(vec![1, 0], r("1")), (vec![0, 1], r("-1"))])
    };
    Arc::new(PolySystem::new(vec![circle, line]).unwrap())
}
fn box2(lo: &str, hi: &str) -> RootBox {
    RootBox::new(vec![r(lo); 2], vec![r(hi); 2]).unwrap()
}

blind_test!(system_certification_unique_empty_multiple_and_singular, {
    let mut d = Decider::default();
    let system = circle_system(false);
    let candidate = box2("1/2", "1");
    let root = match d.certify_root(&system, &candidate) {
        Outcome::Proven(BoxRoots::One(root)) => root,
        out => panic!("{out:?}"),
    };
    for i in 0..2 {
        assert!(
            d.compare(&root.coordinate(i), &real("1/2").sqrt())
                .is_coincident()
        );
        assert!(root.root_box().lo()[i] >= candidate.lo()[i]);
        assert!(root.root_box().hi()[i] <= candidate.hi()[i]);
    }
    assert_eq!(root.system().dimension(), 2);
    assert_eq!(root.enclose(128).len(), 2);
    for b in root.enclose(128) {
        Oracle::rational(&q(1, 2), 1024)
            .sqrt()
            .unwrap()
            .assert_ball(&b);
    }
    assert!(matches!(
        d.certify_root(&system, &box2("-1/10", "1/10")),
        Outcome::Proven(BoxRoots::None)
    ));
    refused(
        d.certify_root(&system, &box2("-1", "1")),
        "numbers.root-not-certified",
    );
    let tangent_box = RootBox::new(vec![r("-1/2"), r("1/2")], vec![r("1/2"), r("3/2")]).unwrap();
    refused(
        d.certify_root(&circle_system(true), &tangent_box),
        "numbers.root-not-certified",
    );
    refused(
        Decider::new(BudgetPolicy {
            max_precision_bits: 1,
            max_work: 0,
        })
        .certify_root(&system, &candidate),
        "numbers.root-not-certified",
    );
    stats_sum(d.stats());
});

blind_test!(
    single_precision_root_certification_preserves_fallback_work,
    {
        let gap = Rational::one().mul_pow2(-160);
        let system = Arc::new(
            PolySystem::new(vec![MPoly::from_terms(
                1,
                vec![(vec![2], r("1")), (vec![0], -(&gap * &gap))],
            )])
            .unwrap(),
        );
        let radius = Rational::one().mul_pow2(-32);
        let both = RootBox::new(vec![-&radius], vec![radius]).unwrap();
        let mut d = Decider::default();
        let ((), work) = with_operation_budget(1_000_000, || {
            refused(
                d.certify_root_at(&system, &both, 64),
                "numbers.root-not-certified",
            );
            let low = Real::from_rational(-&gap);
            let high = Real::from_rational(gap);
            assert_eq!(d.compare(&low, &high), Outcome::Proven(StrictOrder::Less));
        });
        assert!(!work.exhausted);
        stats_sum(d.stats());
    }
);

blind_test!(
    single_precision_root_certification_keeps_closed_box_obligations,
    {
        let system = circle_system(false);
        let candidate = box2("1/2", "1");
        let mut d = Decider::default();
        let Outcome::Proven(BoxRoots::One(root)) = d.certify_root_at(&system, &candidate, 64)
        else {
            panic!("regular root not certified");
        };
        for axis in 0..2 {
            assert!(root.root_box().lo()[axis] >= candidate.lo()[axis]);
            assert!(root.root_box().hi()[axis] <= candidate.hi()[axis]);
            assert!(
                d.compare(&root.coordinate(axis), &real("1/2").sqrt())
                    .is_coincident()
            );
        }
        assert!(matches!(
            d.certify_root_at(&system, &box2("-1/10", "1/10"), 64),
            Outcome::Proven(BoxRoots::None)
        ));
        let linear = Arc::new(PolySystem::new(vec![MPoly::var(1, 0)]).unwrap());
        let boundary = RootBox::new(vec![r("0")], vec![r("1")]).unwrap();
        let Outcome::Proven(BoxRoots::One(root)) = d.certify_root_at(&linear, &boundary, 64) else {
            panic!("closed boundary root lost");
        };
        proven_sign(&mut d, &root.coordinate(0), 0);
        stats_sum(d.stats());
    }
);

blind_test!(single_precision_root_certification_obeys_caller_limits, {
    let system = circle_system(false);
    let candidate = box2("1/2", "1");
    let mut d = Decider::new(BudgetPolicy {
        max_precision_bits: 32,
        max_work: 1_000_000,
    });
    refused(
        d.certify_root_at(&system, &candidate, 64),
        "numbers.precision-limit",
    );
    stats_sum(d.stats());
    let mut d = Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        max_work: 0,
    });
    refused(
        d.certify_root_at(&system, &candidate, 64),
        "numbers.root-not-certified",
    );
    stats_sum(d.stats());
    let mut d = Decider::default();
    let (answer, work) = with_operation_budget(0, || d.certify_root_at(&system, &candidate, 64));
    assert!(answer.is_refused());
    assert!(work.exhausted);
    stats_sum(d.stats());
});

blind_test!(system_certification_with_pinned_coordinate, {
    // The positive solution is (sqrt(1/2), sqrt(1/2), 0).
    // Its Jacobian determinant is -2 sqrt(2), so it is nonsingular.
    let sphere = MPoly::from_terms(
        3,
        vec![
            (vec![2, 0, 0], r("1")),
            (vec![0, 2, 0], r("1")),
            (vec![0, 0, 2], r("1")),
            (vec![0, 0, 0], r("-1")),
        ],
    );
    let line = MPoly::from_terms(3, vec![(vec![1, 0, 0], r("1")), (vec![0, 1, 0], r("-1"))]);
    let system = Arc::new(PolySystem::new(vec![sphere, line, MPoly::var(3, 2)]).unwrap());
    for (zlo, zhi) in [("-1/2", "1/2"), ("0", "1/2"), ("-1/2", "0"), ("0", "0")] {
        let candidate = RootBox::new(
            vec![r("1/2"), r("1/2"), r(zlo)],
            vec![r("1"), r("1"), r(zhi)],
        )
        .unwrap();
        let mut d = Decider::new(BudgetPolicy::default());
        let root = match d.certify_root(&system, &candidate) {
            Outcome::Proven(BoxRoots::One(root)) => root,
            out => panic!("z in [{zlo}, {zhi}]: {out:?}"),
        };
        for axis in 0..2 {
            assert!(
                d.compare(&root.coordinate(axis), &real("1/2").sqrt())
                    .is_coincident()
            );
        }
        proven_sign(&mut d, &root.coordinate(2), 0);
        let balls = root.enclose(128);
        let expected = Oracle::rational(&q(1, 2), 1024).sqrt().unwrap();
        expected.assert_ball(&balls[0]);
        expected.assert_ball(&balls[1]);
        ball_contains(&balls[2], &q(0, 1));
        for axis in 0..3 {
            assert!(root.root_box().lo()[axis] >= candidate.lo()[axis]);
            assert!(root.root_box().hi()[axis] <= candidate.hi()[axis]);
        }
    }
});

blind_test!(system_certification_with_two_pinned_coordinates, {
    // y=1/2 and z=-1/2 leave x^2=3/2, with two simple roots.
    let sphere = MPoly::from_terms(
        3,
        vec![
            (vec![2, 0, 0], r("1")),
            (vec![0, 2, 0], r("1")),
            (vec![0, 0, 2], r("1")),
            (vec![0, 0, 0], r("-2")),
        ],
    );
    let y = MPoly::from_terms(3, vec![(vec![0, 1, 0], r("1")), (vec![0, 0, 0], r("-1/2"))]);
    let z = MPoly::from_terms(3, vec![(vec![0, 0, 1], r("1")), (vec![0, 0, 0], r("1/2"))]);
    let system = Arc::new(PolySystem::new(vec![sphere, y, z]).unwrap());
    for (ylo, yhi, zlo, zhi) in [
        ("0", "1", "-1", "0"),
        ("1/2", "1", "-1", "-1/2"),
        ("1/2", "1/2", "-1/2", "-1/2"),
    ] {
        let candidate =
            RootBox::new(vec![r("1"), r(ylo), r(zlo)], vec![r("3/2"), r(yhi), r(zhi)]).unwrap();
        let mut d = Decider::new(BudgetPolicy::default());
        let root = match d.certify_root(&system, &candidate) {
            Outcome::Proven(BoxRoots::One(root)) => root,
            out => panic!("{candidate:?}: {out:?}"),
        };
        for (axis, expected) in [real("3/2").sqrt(), real("1/2"), real("-1/2")]
            .iter()
            .enumerate()
        {
            assert!(d.compare(&root.coordinate(axis), expected).is_coincident());
            assert!(root.root_box().lo()[axis] >= candidate.lo()[axis]);
            assert!(root.root_box().hi()[axis] <= candidate.hi()[axis]);
        }
        let balls = root.enclose(128);
        Oracle::rational(&q(3, 2), 1024)
            .sqrt()
            .unwrap()
            .assert_ball(&balls[0]);
        ball_contains(&balls[1], &q(1, 2));
        ball_contains(&balls[2], &q(-1, 2));
    }
    let empty = RootBox::new(
        vec![r("1"), r("0"), r("-1")],
        vec![r("6/5"), r("1"), r("0")],
    )
    .unwrap();
    assert!(matches!(
        Decider::default().certify_root(&system, &empty),
        Outcome::Proven(BoxRoots::None)
    ));
    let both = RootBox::new(
        vec![r("-3/2"), r("1/2"), r("-1/2")],
        vec![r("3/2"), r("1/2"), r("-1/2")],
    )
    .unwrap();
    refused(
        Decider::default().certify_root(&system, &both),
        "numbers.root-not-certified",
    );
});

blind_test!(three_sphere_roots_and_closed_boundary_boxes, {
    let system = sphere_system();
    let mut d = Decider::default();
    for positive in [false, true] {
        let candidate = sphere_box(positive);
        let root = match d.certify_root(&system, &candidate) {
            Outcome::Proven(BoxRoots::One(root)) => root,
            out => panic!("{out:?}"),
        };
        assert!(d.compare(&root.coordinate(0), &real("3")).is_coincident());
        assert!(d.compare(&root.coordinate(1), &real("3")).is_coincident());
        let expected = if positive {
            real("7").sqrt()
        } else {
            -real("7").sqrt()
        };
        assert!(d.compare(&root.coordinate(2), &expected).is_coincident());
        assert!(
            d.compare(&Real::from_root_coordinate(&root, 2), &expected)
                .is_coincident()
        );
        let balls = root.enclose(256);
        assert_eq!(balls.len(), 3);
        ball_contains(&balls[0], &q(3, 1));
        ball_contains(&balls[1], &q(3, 1));
        let oracle = Oracle::rational(&q(7, 1), 1024).sqrt().unwrap();
        let oracle = if positive { oracle } else { oracle.neg() };
        oracle.assert_ball(&balls[2]);
        for i in 0..3 {
            assert!(root.root_box().lo()[i] >= candidate.lo()[i]);
            assert!(root.root_box().hi()[i] <= candidate.hi()[i]);
        }
    }
    // x=0 is on a closed candidate's boundary. Refusal is allowed, exclusion is not.
    let linear = Arc::new(PolySystem::new(vec![MPoly::var(1, 0)]).unwrap());
    let boundary = RootBox::new(vec![r("0")], vec![r("1")]).unwrap();
    match d.certify_root(&linear, &boundary) {
        Outcome::Proven(BoxRoots::One(root)) => proven_sign(&mut d, &root.coordinate(0), 0),
        Outcome::Refused(e) => assert_eq!(e.code, "numbers.root-not-certified"),
        out => panic!("closed boundary lost: {out:?}"),
    }
    assert!(PolySystem::new(vec![]).is_none());
    assert!(PolySystem::new(vec![MPoly::var(2, 0)]).is_none());
    assert!(RootBox::new(vec![], vec![]).is_none());
    assert!(RootBox::new(vec![r("1")], vec![r("0")]).is_none());
    assert!(RootBox::new(vec![r("0")], vec![r("1"), r("2")]).is_none());
});

fn enclosed(value: &Real) -> Ball {
    Decider::default().enclose_at(value, 256).unwrap()
}

#[test]
fn a_simple_root_stays_exact_in_a_quadratic_tower_without_an_eliminant() {
    let q = |text: &str| Rational::parse(text).unwrap();
    let root = Real::from_i64(2).sqrt();
    let nested = (&root + &Real::from_i64(3)).sqrt();
    // t^3 + sqrt(2) t - sqrt(3 + sqrt(2)), and s = t sqrt(2).
    let t = RealPoly::var(2, 0);
    let s = RealPoly::var(2, 1);
    let cubic = &(&t * &t) * &t + t.scale(&root) - RealPoly::constant(2, nested.clone());
    let line = s - t.scale(&root);
    let (solved, work) = with_operation_budget(1_000_000, || {
        RealPoly::simple_root_in(
            &[cubic, line],
            &[[q("9/10"), q("19/20")], [q("6/5"), q("7/5")]],
            128,
            &mut Decider::default(),
        )
    });
    assert!(!work.exhausted);
    let solved = solved.unwrap().unwrap();
    let residual = enclosed(
        &(&(&(&solved[0] * &solved[0]) * &solved[0]) + &(&root * &solved[0]) - nested.clone()),
    );
    assert!(residual.lo().abs().max(residual.hi().abs()) < Rational::one().mul_pow2(-40));
    // No root in the box, and a coefficient off every quadratic tower.
    let lone = |value: &Real| {
        RealPoly::var(1, 0) * RealPoly::var(1, 0) - RealPoly::constant(1, value.clone())
    };
    let empty = RealPoly::simple_root_in(
        &[lone(&root)],
        &[[q("3"), q("4")]],
        128,
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(empty.code, "numbers.root-not-certified");
    let x = RealPoly::var(1, 0);
    let irreducible = &(&x * &x) * &x - x.clone() - RealPoly::constant(1, Real::from_i64(2));
    let cubic_root = irreducible.real_roots(&mut Decider::default()).unwrap()[0].clone();
    let leaves = RealPoly::simple_root_in(
        &[lone(&cubic_root)],
        &[[q("0"), q("2")]],
        128,
        &mut Decider::default(),
    );
    assert!(matches!(leaves, Ok(None)), "{leaves:?}");
}
