#![allow(missing_docs, clippy::float_cmp)]
mod common;
use catallax_numbers::*;
use common::*;
use rug::{Float, Rational as Q};

blind_test!(interval_endpoints_infinities_and_domains, {
    assert!(Interval::new(2.0, 1.0).is_none());
    assert!(Interval::point(f64::NAN).is_none());
    assert!(Interval::new(0.0, f64::NAN).is_none());
    let zero = Interval::point(-0.0).unwrap();
    assert_eq!(zero.lo().to_bits(), 0);
    assert_eq!(zero.hi().to_bits(), 0);
    assert_eq!(zero, Interval::point(0.0).unwrap());
    let entire = Interval::entire();
    assert_eq!(entire.mid(), 0.0);
    assert_eq!(Interval::new(1.0, f64::INFINITY).unwrap().mid(), f64::MAX);
    assert_eq!(
        Interval::new(f64::NEG_INFINITY, 1.0).unwrap().mid(),
        -f64::MAX
    );
    assert_eq!(zero * entire, zero);
    assert_eq!(entire * zero, zero);
    assert_eq!(zero / zero, entire);
    assert_eq!(
        Interval::point(1.0).unwrap() / Interval::new(-1.0, 1.0).unwrap(),
        entire
    );
    assert!(Interval::new(-2.0, -1.0).unwrap().sqrt().is_none());
    assert_eq!(
        Interval::new(-1.0, 4.0).unwrap().sqrt().unwrap(),
        Interval::new(0.0, 2.0).unwrap()
    );
    assert_eq!(Interval::new(0.0, 1.0).unwrap().sign(), None);
    assert_eq!(Interval::new(-1.0, 0.0).unwrap().sign(), None);
    let overflow = Interval::point(1e300).unwrap() * Interval::point(1e300).unwrap();
    assert_eq!(overflow.sign(), Some(StrictSign::Positive));
    assert_eq!(overflow.hi(), f64::INFINITY);
    assert!(!overflow.lo().is_nan());
    let tenth = Interval::from_rational(&r("1/10"));
    interval_contains(tenth, &q(1, 10));
    assert!(Q::from_f64(tenth.lo()).unwrap() < q(1, 10));
    assert!(Q::from_f64(tenth.hi()).unwrap() > q(1, 10));
    assert_eq!(tenth.lo().to_bits() + 1, tenth.hi().to_bits());
    assert_eq!(
        Interval::from_rational(&r("1/2")),
        Interval::point(0.5).unwrap()
    );
    let mut sum = zero;
    for _ in 0..10 {
        sum = sum + tenth;
    }
    interval_contains(sum, &q(1, 1));
    for text in ["1e-320", "1e-400", "-1e-320", "1e400", "-1e400"] {
        let x = r(text);
        let i = Interval::from_rational(&x);
        interval_contains(i, &oq(&x));
        if i.sign() == Some(StrictSign::Positive) {
            assert!(i.lo() > 0.0);
        }
    }
    let x = Interval::new(-2.0, 3.0).unwrap();
    for v in [0, 4, 9] {
        interval_contains(x.square(), &q(v, 1));
    }
});

blind_test!(generated_interval_arithmetic_contains_exact_answers, {
    let mut rng = Rng(0x1a7e2);
    for _ in 0..8192 {
        let (a, b) = (rng.rational(), rng.rational());
        let (x, y) = (
            Interval::from_rational(&eq(&a)),
            Interval::from_rational(&eq(&b)),
        );
        interval_contains(x + y, &(a.clone() + &b));
        interval_contains(x - y, &(a.clone() - &b));
        interval_contains(x * y, &(a.clone() * &b));
        if b != 0 {
            interval_contains(x / y, &(a.clone() / &b));
        }
        interval_contains(-x, &(-a.clone()));
        interval_contains(x.square(), &(a.clone() * &a));
        let h = x.hull(&y);
        interval_contains(h, &a);
        interval_contains(h, &b);
        assert!(x.width() >= x.hi() - x.lo());
    }
});

blind_test!(ball_precision_hulls_and_domain_edges, {
    let b = Ball::from_rational(&r("1/3"), 128);
    ball_contains(&b, &q(1, 3));
    let reduced = b.with_precision(16);
    assert_eq!(reduced.precision(), 16);
    assert!(oq(&reduced.lo()) <= oq(&b.lo()));
    assert!(oq(&reduced.hi()) >= oq(&b.hi()));
    let raised = b.with_precision(256);
    assert_eq!(raised.mid(), b.mid());
    assert_eq!(raised.rad(), b.rad());
    assert_eq!(raised.precision(), 256);
    assert_eq!(Ball::from_rational(&r("1"), 0).precision(), 1);
    let a = Ball::from_rational(&r("-1"), 64);
    let c = Ball::from_rational(&r("4"), 128);
    let h = a.hull(&c);
    assert_eq!(h.precision(), 128);
    ball_contains(&h, &q(-1, 1));
    ball_contains(&h, &q(4, 1));
    let sqrt = h.sqrt().unwrap();
    ball_contains(&sqrt, &q(0, 1));
    ball_contains(&sqrt, &q(2, 1));
    assert!(a.sqrt().is_none());
    assert!(h.ln().is_none());
    assert!(h.asin().is_none());
    assert!(h.acos().is_none());
    assert!(h.recip().is_none());
    assert!(c.checked_div(&h).is_none());
    assert!(Ball::atan2(&h, &h).is_none());
    for x in [-1, 1] {
        let v = Ball::from_rational(&Rational::from_i64(x), 128);
        assert!(v.asin().is_some());
        assert!(v.acos().is_some());
    }
    assert_eq!((&a + &c).precision(), 128);
    assert_eq!((&a * &c).precision(), 128);
    let exact = Ball::from_rational(&r("5/8"), 3);
    assert_q(&exact.rad(), &q(0, 1));
    assert_q(&exact.mid(), &q(5, 8));
    let zero = Ball::from_rational(&r("0"), 128);
    assert!(zero.ln().is_none());
    assert!(zero.recip().is_none());
});

blind_test!(ball_arithmetic_encloses_every_interval_sign_pair, {
    let ranges = [(-5, -2), (-5, 3), (2, 5), (-3, -3), (0, 0), (3, 3)];
    for precision in [1, 3, 128] {
        let point = |n| Ball::from_rational(&Rational::from_i64(n), precision);
        for (a, b) in ranges {
            for (c, d) in ranges {
                let x = point(a).hull(&point(b));
                let y = point(c).hull(&point(d));
                let (sum, difference, product) = (&x + &y, &x - &y, &x * &y);
                for u in [a, b] {
                    for v in [c, d] {
                        ball_contains(&sum, &q(u + v, 1));
                        ball_contains(&difference, &q(u - v, 1));
                        ball_contains(&product, &q(u * v, 1));
                    }
                }
            }
        }
    }
});

#[test]
fn ball_alignment_handles_signed_exponent_gaps_and_exact_cancellation() {
    for precision in [1, 53, 1024] {
        for left_power in [-257, -1, 0, 1, 257] {
            for right_power in [-257, -1, 0, 1, 257] {
                for sign in [-1, 1] {
                    let left = Rational::from_ratio(3 * sign, 7)
                        .unwrap()
                        .mul_pow2(left_power);
                    let right = Rational::from_ratio(5, 11).unwrap().mul_pow2(right_power);
                    let a = Ball::from_rational(&left, precision);
                    let b = Ball::from_rational(&right, precision);
                    let x = oq(&left);
                    let y = oq(&right);
                    ball_contains(&(&a + &b), &(x.clone() + &y));
                    ball_contains(&(&a - &b), &(x.clone() - &y));
                    ball_contains(&(&b - &a), &(y.clone() - &x));
                    ball_contains(&(&a * &b), &(x * y));
                }
                let value = Rational::from_i64(3).mul_pow2(left_power);
                let a = Ball::from_rational(&value, precision.max(2));
                let cancelled = &a - &a;
                assert!(cancelled.is_exact_zero());
                let doubled = &a + &a;
                assert_eq!(doubled.mid(), value.mul_pow2(1));
                assert_eq!(doubled.rad(), Rational::zero());
            }
        }
    }
}

blind_test!(generated_ball_arithmetic_and_elementary_functions, {
    let mut rng = Rng(0xba11);
    for _ in 0..512 {
        // Dyadic points are represented exactly. Restrict |x|<=1/2 so sine is
        // increasing and cosine is decreasing on the nonnegative test inputs.
        let n = (rng.next() % 129) as i64;
        let x = q(n, 256);
        let y = q((rng.next() % 129 + 1) as i64, 256);
        let a = Ball::from_rational(&eq(&x), 96);
        let b = Ball::from_rational(&eq(&y), 160);
        let ox = Oracle::rational(&x, 768);
        let oy = Oracle::rational(&y, 768);
        ball_contains(&(&a + &b), &(x.clone() + &y));
        ball_contains(&(&a - &b), &(x.clone() - &y));
        ball_contains(&(&a * &b), &(x.clone() * &y));
        ball_contains(&a.checked_div(&b).unwrap(), &(x.clone() / &y));
        ball_contains(&a.square(), &(x.clone() * &x));
        ball_contains(&a.pow_u32(3), &(x.clone() * &x * &x));
        ox.sqrt().unwrap().assert_ball(&a.sqrt().unwrap());
        ox.unary(Float::exp_round, true).assert_ball(&a.exp());
        ox.unary(Float::sin_round, true).assert_ball(&a.sin());
        ox.unary(Float::cos_round, false).assert_ball(&a.cos());
        ox.unary(Float::atan_round, true).assert_ball(&a.atan());
        ox.unary(Float::asin_round, true)
            .assert_ball(&a.asin().unwrap());
        ox.unary(Float::acos_round, false)
            .assert_ball(&a.acos().unwrap());
        oy.unary(Float::ln_round, true)
            .assert_ball(&b.ln().unwrap());
        Oracle::rational(&q(1, 1), 768)
            .div(&oy)
            .unwrap()
            .assert_ball(&b.recip().unwrap());
        // x>=0,y>0: atan2(x,y)=atan(x/y), with no branch cut.
        ox.div(&oy)
            .unwrap()
            .unary(Float::atan_round, true)
            .assert_ball(&Ball::atan2(&a, &b).unwrap());
        interval_contains(a.to_interval(), &x);
        for input in [&a, &(-a.clone())] {
            let first = input.exp();
            let second = input.exp();
            assert_eq!(first.mid(), second.mid());
            assert_eq!(first.rad(), second.rad());
        }
    }
});

blind_test!(high_precision_constants_and_function_identities, {
    // MPFR's directed constants give certified reference intervals beyond the
    // requested 3000 decimal places of pi and 1000 of e and ln(2).
    let opi = Oracle::pi(12000);
    for p in [64, 1000, 10000] {
        let pi = Ball::pi(p);
        opi.assert_ball(&pi);
        ball_contains(&pi.sin(), &q(0, 1));
        ball_contains(&pi.cos(), &q(-1, 1));
        let one = Ball::from_rational(&r("1"), p);
        let two = Ball::from_rational(&r("2"), p);
        let four = Ball::from_rational(&r("4"), p);
        overlap(&(one.atan() * four), &pi);
        ball_contains(&two.ln().unwrap().exp(), &q(2, 1));
        overlap(&one.asin().unwrap(), &pi.checked_div(&two).unwrap());
        overlap(
            &Ball::atan2(&one, &(-one.clone())).unwrap(),
            &(pi.clone() * Ball::from_rational(&r("3/4"), p)),
        );
        let again = Ball::pi(p);
        assert_eq!(pi.mid(), again.mid());
        assert_eq!(pi.rad(), again.rad());
        eprintln!("quality pi bits={p} radius={}", pi.rad());
    }
    Oracle::rational(&q(1, 1), 4096)
        .unary(Float::exp_round, true)
        .assert_ball(&Ball::from_rational(&r("1"), 3500).exp());
    Oracle::rational(&q(2, 1), 4096)
        .unary(Float::ln_round, true)
        .assert_ball(&Ball::from_rational(&r("2"), 3500).ln().unwrap());
});

#[test]
fn oracle_directed_rounding_smoke() {
    // Runs even against the interface skeleton. This validates the independent
    // oracle path, not any catallax-numbers implementation.
    let sqrt = Oracle::rational(&q(2, 1), 256).sqrt().unwrap();
    let lo = sqrt.lo.to_rational().unwrap();
    let hi = sqrt.hi.to_rational().unwrap();
    assert!(lo.clone() * lo < 2);
    assert!(hi.clone() * hi > 2);
    let pi = Oracle::pi(256);
    assert!(pi.lo > 3);
    assert!(pi.hi < 4);
    let exp = Oracle::rational(&q(0, 1), 256).unary(Float::exp_round, true);
    assert_eq!(exp.lo, 1);
    assert_eq!(exp.hi, 1);
}

blind_test!(balls_cover_ranges_extrema_and_atan2_quadrants, {
    let ball = |x: &str| Ball::from_rational(&r(x), 192);
    let broad = ball("-10").hull(&ball("10"));
    for result in [broad.sin(), broad.cos()] {
        ball_contains(&result, &q(-1, 1));
        ball_contains(&result, &q(1, 1));
    }
    let neg = ball("-1/2");
    let oracle = Oracle::rational(&q(-1, 2), 1024);
    oracle.unary(Float::sin_round, true).assert_ball(&neg.sin());
    oracle.unary(Float::cos_round, true).assert_ball(&neg.cos());
    oracle
        .unary(Float::atan_round, true)
        .assert_ball(&neg.atan());
    oracle
        .unary(Float::asin_round, true)
        .assert_ball(&neg.asin().unwrap());
    oracle
        .unary(Float::acos_round, false)
        .assert_ball(&neg.acos().unwrap());
    for (y, x, pi_multiple) in [
        ("1", "1", q(1, 4)),
        ("1", "-1", q(3, 4)),
        ("-1", "-1", q(-3, 4)),
        ("-1", "1", q(-1, 4)),
        ("0", "-1", q(1, 1)),
    ] {
        Oracle::pi(1024)
            .mul(&Oracle::rational(&pi_multiple, 1024))
            .assert_ball(&Ball::atan2(&ball(y), &ball(x)).unwrap());
    }
    // Crossing the negative x-axis includes angles approaching both -pi and pi.
    let crossing = Ball::atan2(&ball("-1/8").hull(&ball("1/8")), &ball("-1")).unwrap();
    ball_contains(&crossing, &q(-3, 1));
    ball_contains(&crossing, &q(3, 1));
    assert_eq!(ball("1").sign(), Some(StrictSign::Positive));
    assert_eq!(ball("-1").sign(), Some(StrictSign::Negative));
    assert_eq!(ball("0").sign(), None);
    let range = ball("1/2").hull(&ball("2"));
    for x in [q(1, 2), q(2, 1)] {
        Oracle::rational(&x, 1024)
            .unary(Float::ln_round, true)
            .assert_ball(&range.ln().unwrap());
        Oracle::rational(&x, 1024)
            .unary(Float::exp_round, true)
            .assert_ball(&range.exp());
    }
});

#[test]
fn independent_sturm_oracle_counts_known_roots() {
    let square_free = sturm(&[q(-2, 1), q(0, 1), q(1, 1)]);
    assert_eq!(
        variations(&square_free, &q(-2, 1)) - variations(&square_free, &q(2, 1)),
        2
    );
    assert_eq!(
        variations(&square_free, &q(1, 1)) - variations(&square_free, &q(2, 1)),
        1
    );
    assert_eq!(
        variations(&square_free, &q(-1, 1)) - variations(&square_free, &q(1, 1)),
        0
    );
    let repeated = sturm(&[q(1, 1), q(-2, 1), q(1, 1)]);
    assert_eq!(repeated.last().unwrap().len(), 2);
}
