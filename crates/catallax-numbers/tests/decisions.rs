#![allow(missing_docs)]
mod common;
use catallax_numbers::*;
use common::*;
use rug::ops::Pow;
use rug::{Integer as Z, Rational as Q};
use std::time::{Duration, Instant};

blind_test!(orientation_2d_all_65536_points_and_all_scales, {
    let mut d = Decider::default();
    for k in -30i32..=30 {
        let scale = Q::from(10).pow(k);
        let s = Real::from_rational(eq(&scale));
        let qx = real("12") * s.clone();
        let rx = real("24") * s.clone();
        for i in 0i32..256 {
            for j in 0i32..256 {
                let px =
                    Real::from_rational(eq(&(q(1, 2) + Q::from((Z::from(i), Z::from(1) << 53)))))
                        * s.clone();
                let py =
                    Real::from_rational(eq(&(q(1, 2) + Q::from((Z::from(j), Z::from(1) << 53)))))
                        * s.clone();
                let determinant = (&qx - &px) * (&rx - &py) - (&qx - &py) * (&rx - &px);
                // Expand: determinant=12*(py-px)=12*(j-i)*2^-53.
                // Scaling both axes multiplies it by the positive factor 10^(2k).
                proven_sign(&mut d, &determinant, (j - i).signum() as i8);
            }
        }
    }
});

blind_test!(orientation_3d_near_plane_exact_grid_and_scales, {
    let mut rng = Rng(0x3d0a1);
    let mut d = Decider::default();
    for _ in 0..256 {
        let u = rng.small();
        let v = rng.small();
        let w = (rng.next() % 3) as i64 - 1;
        for k in -30i32..=30 {
            let scale = Q::from(10).pow(k);
            let cv = |v: Q| Real::from_rational(eq(&(v * &scale)));
            // Columns b-a=(1,0,1), c-a=(0,1,1), d-a=(u,v,u+v+w*2^-100).
            // det = w*2^-100. A large translation tests exact cancellation.
            let shift = Q::from(Z::from(1) << 80);
            let a = [cv(shift.clone()), cv(shift.clone()), cv(shift.clone())];
            let b = [
                cv(shift.clone() + 1),
                cv(shift.clone()),
                cv(shift.clone() + 1),
            ];
            let c = [
                cv(shift.clone()),
                cv(shift.clone() + 1),
                cv(shift.clone() + 1),
            ];
            let e = [
                cv(shift.clone() + u),
                cv(shift.clone() + v),
                cv(shift + u + v + Q::from((Z::from(w), Z::from(1) << 100))),
            ];
            let b: Vec<_> = (0..3).map(|i| &b[i] - &a[i]).collect();
            let c: Vec<_> = (0..3).map(|i| &c[i] - &a[i]).collect();
            let e: Vec<_> = (0..3).map(|i| &e[i] - &a[i]).collect();
            let det = &b[0] * &(&c[1] * &e[2] - &c[2] * &e[1])
                - &b[1] * &(&c[0] * &e[2] - &c[2] * &e[0])
                + &b[2] * &(&c[0] * &e[1] - &c[1] * &e[0]);
            proven_sign(&mut d, &det, w as i8);
        }
    }
});

blind_test!(known_signs_hidden_zeros_and_scale_invariance, {
    let mut d = Decider::default();
    let exact_zero = real("0.1") + real("0.2") - real("0.3");
    assert_eq!(
        d.sign(&exact_zero),
        Outcome::ProvenCoincident(Coincidence {
            proof: ZeroProof::ExactArithmetic
        })
    );
    let mut cases: Vec<(Real, i8)> = hidden_zeros().into_iter().map(|x| (x, 0)).collect();
    cases.push((exact_zero, 0));
    cases.push((near_tie(), 1));
    // sqrt is strictly concave on positive arguments: sqrt(N+1)+sqrt(N-1)<2sqrt(N).
    for n in ["1000000", "100000000000000000000"] {
        let n = real(n);
        cases.push((
            (&n + &real("1")).sqrt() + (&n - &real("1")).sqrt() - real("2") * n.sqrt(),
            -1,
        ));
    }
    let a = Real::from_rational(Rational::one().mul_pow2(1100));
    cases.push((&a * &a - &a * &a + real("1"), 1));
    cases.push((&a * &a - &a * &a, 0));
    for k in -30i32..=30 {
        let scale = Real::from_rational(eq(&Q::from(10).pow(k)));
        for (x, sign) in &cases {
            proven_sign(&mut d, &(&scale * x), *sign);
        }
    }
    let mut rng = Rng(0x5a27);
    for _ in 0..256 {
        let x = eq(&rng.rational().abs());
        let y = eq(&rng.rational().abs());
        // Principal roots are nonnegative; both sides square to x*y.
        let e = Real::from_rational(x.clone()).sqrt() * Real::from_rational(y.clone()).sqrt()
            - Real::from_rational(x * y).sqrt();
        for k in -30i32..=30 {
            proven_sign(
                &mut d,
                &(Real::from_rational(eq(&Q::from(10).pow(k))) * e.clone()),
                0,
            );
        }
    }
});

blind_test!(tiny_budgets_never_invent_signs_or_coincidence, {
    let x = near_tie();
    let mut d = Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        max_work: 10_000,
    });
    match d.sign(&x) {
        Outcome::Proven(StrictSign::Positive) => (),
        Outcome::NotExact(ne) => {
            assert!(matches!(
                ne.code,
                "numbers.precision-limit" | "numbers.budget-exhausted"
            ));
            let b = oq(&ne.bound);
            assert!(b >= 0);
            // x=sqrt(N^2+1)-N. b>=x iff (N+b)^2>=N^2+1, since N+b>=0.
            let n = Q::from(z("100000000000000000000"));
            let radicand = n.clone() * &n + q(1, 1);
            assert!((n.clone() + &b) * (n + b) >= radicand);
        }
        out => panic!("unsafe or unexpected near-tie result: {out:?}"),
    }
    for max_work in [0, 1, 8, 64] {
        let mut zeros = hidden_zeros();
        zeros.push(real("2/3").sqrt() * real("5/7").sqrt() - real("10/21").sqrt());
        for x in zeros {
            let mut d = Decider::new(BudgetPolicy {
                max_precision_bits: 64,
                max_work,
            });
            match d.sign(&x) {
                Outcome::ProvenCoincident(c) => assert_ne!(c.proof, ZeroProof::ConstructionFact),
                Outcome::NotExact(ne) => assert!(oq(&ne.bound) >= 0),
                out => panic!("zero under tiny budget: {out:?}"),
            }
            assert!(d.stats().work <= max_work);
            stats_sum(d.stats());
        }
    }
});

blind_test!(undefined_subexpressions_and_postorder_precedence, {
    let mut d = Decider::default();
    let zero = real("2").sqrt().square() - real("2");
    refused(
        d.sign(&(real("1") / zero.clone())),
        "numbers.division-by-zero",
    );
    refused(d.sign(&real("-1").sqrt()), "numbers.sqrt-of-negative");
    proven_sign(&mut d, &zero.sqrt(), 0);
    let bad_div = real("1") / real("0");
    let bad_sqrt = real("-1").sqrt();
    for x in [
        real("0") * bad_div.clone(),
        bad_div.clone() - bad_div.clone(),
        bad_div.pow_u32(0),
    ] {
        refused(d.sign(&x), "numbers.division-by-zero");
    }
    refused(
        d.sign(&(bad_sqrt.clone() + bad_div.clone())),
        "numbers.sqrt-of-negative",
    );
    refused(
        d.sign(&(bad_div.clone() + bad_sqrt.clone())),
        "numbers.division-by-zero",
    );
    refused(d.sign(&(bad_sqrt / real("0"))), "numbers.sqrt-of-negative");
    assert_eq!(
        d.enclose(&bad_div, &r("1/100")).unwrap_err().code,
        "numbers.division-by-zero"
    );
    refused(d.compare(&bad_div, &bad_div), "numbers.division-by-zero");
    let budget = BudgetPolicy {
        max_precision_bits: 64,
        max_work: 10_000,
    };
    let n = Q::from(z("100000000000000000000"));
    let exact = Oracle::rational(&(n.clone() * &n + 1), 1024)
        .sqrt()
        .unwrap()
        .sub(&Oracle::rational(&n, 1024));
    let out = Decider::new(budget).sign(&(real("1") / near_tie()));
    match out {
        Outcome::Refused(e) => assert_eq!(e.code, "numbers.divisor-not-separated"),
        Outcome::Proven(StrictSign::Positive) => (),
        Outcome::NotExact(ne) => Oracle::rational(&q(1, 1), 1024)
            .div(&exact)
            .unwrap()
            .assert_bound(&ne.bound),
        out => panic!("{out:?}"),
    }
    let out = Decider::new(budget).sign(&near_tie().sqrt());
    match out {
        Outcome::Refused(e) => assert_eq!(e.code, "numbers.sqrt-not-separated"),
        Outcome::Proven(StrictSign::Positive) => (),
        Outcome::NotExact(ne) => exact.sqrt().unwrap().assert_bound(&ne.bound),
        out => panic!("{out:?}"),
    }
});

// Exact oracle accompanies every generated expression. Square-root branches
// take sqrt(v^2)=|v|, or a deliberately negative argument. This proves zeros
// and undefinedness, which a finite-precision floating oracle cannot do alone.
fn expression(rng: &mut Rng, depth: u32) -> (Real, Option<Q>) {
    if depth == 0 {
        let v = rng.rational();
        return (Real::from_rational(eq(&v)), Some(v));
    }
    let op = rng.next() % 5;
    if op == 4 && depth >= 2 {
        if rng.next().is_multiple_of(5) {
            if depth < 4 {
                return (real("-1").sqrt(), None);
            }
            let (x, _) = expression(rng, depth - 4);
            return ((-x.square() - real("1")).sqrt(), None);
        }
        let (x, v) = expression(rng, depth - 2);
        return (x.square().sqrt(), v.map(|v| v.abs()));
    }
    let (x, a) = expression(rng, depth - 1);
    let (y, b) = expression(rng, depth - 1);
    let value = match (a, b) {
        (Some(a), Some(b)) => match op % 4 {
            0 => Some(a + b),
            1 => Some(a - b),
            2 => Some(a * b),
            _ => {
                if b == 0 {
                    None
                } else {
                    Some(a / b)
                }
            }
        },
        _ => None,
    };
    (
        match op % 4 {
            0 => x + y,
            1 => x - y,
            2 => x * y,
            _ => x / y,
        },
        value,
    )
}

blind_test!(generated_depth_six_expressions_with_random_budgets, {
    let mut rng = Rng(0xe8a355);
    for _ in 0..4096 {
        let depth = (rng.next() % 7) as u32;
        let (x, value) = expression(&mut rng, depth);
        let policy = BudgetPolicy {
            max_precision_bits: [1, 64, 128, 512][(rng.next() % 4) as usize],
            max_work: [0, 1, 32, 10_000, 1_000_000][(rng.next() % 5) as usize],
        };
        let mut d = Decider::new(policy);
        let out = d.sign(&x);
        match (out, value) {
            (Outcome::Refused(_), _) => (),
            (Outcome::Proven(s), Some(v)) => assert!(if s == StrictSign::Positive {
                v > 0
            } else {
                v < 0
            }),
            (Outcome::ProvenCoincident(c), Some(v)) => {
                assert_eq!(v, 0);
                assert_ne!(c.proof, ZeroProof::ConstructionFact);
            }
            (Outcome::NotExact(ne), Some(v)) => {
                assert!(oq(&ne.bound) >= v.abs());
                assert!(matches!(
                    ne.code,
                    "numbers.precision-limit" | "numbers.budget-exhausted"
                ));
            }
            (out, None) => panic!("undefined expression accepted: {out:?}"),
        }
        stats_sum(d.stats());
        assert_eq!(d.stats().decisions, 1);
        assert!(d.stats().work <= policy.max_work);
    }
});

blind_test!(deadline_all_decisions_no_work_and_reset, {
    let mut d = Decider::default();
    d.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
    refused(d.sign(&real("1")), "numbers.deadline");
    refused(d.compare(&real("1"), &real("2")), "numbers.deadline");
    assert_eq!(
        d.enclose(&real("1"), &r("1")).unwrap_err().code,
        "numbers.deadline"
    );
    refused(d.real_roots(&poly(&[-2, 0, 1])), "numbers.deadline");
    refused(
        d.certify_root(&sphere_system(), &sphere_box(true)),
        "numbers.deadline",
    );
    assert_eq!(d.stats().decisions, 5);
    assert_eq!(d.stats().refused, 5);
    assert_eq!(d.stats().work, 0);
    stats_sum(d.stats());
    d.set_deadline(None);
    proven_sign(&mut d, &real("1"), 1);
    let old = d.take_stats();
    assert_eq!(old.decisions, 6);
    assert_eq!(d.stats(), &DecisionStats::default());
    d.note_fact_decision();
    assert_eq!(d.stats().coincident_fact, 1);
    assert_eq!(d.stats().decisions, 1);
    stats_sum(d.stats());
});

blind_test!(compare_matches_subtraction_and_threads_are_deterministic, {
    fn run() -> Vec<(Outcome<StrictSign>, DecisionStats)> {
        let mut d = Decider::default();
        let mut out = Vec::new();
        let mut values = hidden_zeros();
        values.extend([real("7"), near_tie(), real("-1").sqrt()]);
        for x in values {
            let first = d.sign(&x);
            let first_stats = d.take_stats();
            let second = d.sign(&x);
            let second_stats = d.take_stats();
            assert_eq!(first, second);
            assert_eq!(first_stats, second_stats);
            stats_sum(&first_stats);
            out.push((first, first_stats));
            let mut a = Decider::default();
            let mut b = Decider::default();
            let right = real("1/3");
            let compare = a.compare(&x, &right);
            let sign = b.sign(&(&x - &right)).map(|s| match s {
                StrictSign::Negative => StrictOrder::Less,
                StrictSign::Positive => StrictOrder::Greater,
            });
            assert_eq!(compare, sign);
            assert_eq!(a.stats(), b.stats());
        }
        out
    }
    let expected = run();
    let workers: Vec<_> = (0..8).map(|_| std::thread::spawn(run)).collect();
    for worker in workers {
        assert_eq!(worker.join().unwrap(), expected);
    }
});

blind_test!(enclose_width_accounting_and_rational_evaluation, {
    let mut d = Decider::default();
    let x = real("2").sqrt();
    let width = r("1e-50");
    let b = d.enclose(&x, &width).unwrap();
    Oracle::rational(&q(2, 1), b.precision().saturating_add(256).max(1024))
        .sqrt()
        .unwrap()
        .assert_ball(&b);
    assert!(b.hi() - b.lo() <= width);
    assert_eq!(d.stats().refined, 1);
    assert_eq!(d.stats().decisions, 1);
    let mut tiny = Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        max_work: 1,
    });
    let b = tiny.enclose(&x, &r("1e-100")).unwrap();
    Oracle::rational(&q(2, 1), 1024)
        .sqrt()
        .unwrap()
        .assert_ball(&b);
    assert!(b.hi() - b.lo() > r("1e-100"));
    assert_eq!(tiny.stats().not_exact, 1);
    let rational = (real("7/3") - real("1/3")) * real("3") / real("2");
    assert_q(&rational.as_rational().unwrap(), &q(3, 1));
    interval_contains(rational.approx(), &q(3, 1));
});
