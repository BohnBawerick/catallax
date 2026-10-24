#![allow(clippy::float_cmp)] // Bit-exact conversion and endpoint tests.

use super::*;
use std::{
    collections::hash_map::DefaultHasher,
    hash::{Hash, Hasher},
    sync::Arc,
    time::{Duration, Instant},
};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn real(s: &str) -> Real {
    Real::from_rational(q(s))
}
fn ball(s: &str, p: u32) -> Ball {
    Ball::from_rational(&q(s), p)
}
fn hash(r: &Rational) -> u64 {
    let mut h = DefaultHasher::new();
    r.hash(&mut h);
    h.finish()
}
fn overlaps(a: &Ball, b: &Ball) -> bool {
    a.lo() <= b.hi() && b.lo() <= a.hi()
}

#[test]
fn exact_parsing_normalization_and_error_precedence() {
    for (s, expected) in [
        ("1.1", "11/10"),
        ("0.1", "1/10"),
        ("1e-7", "1/10000000"),
        ("2.5E3", "2500"),
        ("-6/8", "-3/4"),
        ("-0", "0"),
        ("+007", "7"),
    ] {
        assert_eq!(q(s).to_string(), expected);
    }
    for (s, expected) in [
        ("", ParseErrorKind::Empty),
        ("1/0", ParseErrorKind::ZeroDenominator),
        ("0/000", ParseErrorKind::ZeroDenominator),
        ("1e100001", ParseErrorKind::ExponentTooLarge),
        ("1e-000100001", ParseErrorKind::ExponentTooLarge),
    ] {
        assert_eq!(Rational::parse(s).unwrap_err().kind, expected);
    }
    for s in [
        ".",
        "1.",
        ".5",
        "1.1.1",
        " 1",
        "1 ",
        "1_000",
        "0x10",
        "inf",
        "nan",
        "1e",
        "1/-2",
        "-6/-4",
        "   ",
        "1/0/2",
        "1e100001x",
        "1.0/0",
    ] {
        assert_eq!(
            Rational::parse(s).unwrap_err().kind,
            ParseErrorKind::InvalidSyntax,
            "{s}"
        );
    }
    let a = Rational::new(Integer::from_i64(-6), Integer::from_i64(-4)).unwrap();
    assert_eq!(a, q("3/2"));
    assert_eq!(hash(&a), hash(&q("6/4")));
    assert_eq!(hash(&q("-0")), hash(&q("0/5")));
    assert_eq!(&q("0.1") + &q("0.2"), q("0.3"));
    assert_eq!(Rational::zero().pow(0), Some(Rational::one()));
    assert_eq!(Integer::zero().pow(0), Integer::one());
    assert_eq!(Rational::zero().pow(-1), None);
}

#[test]
fn floor_division_and_large_integer_field_laws() {
    for (a, d, quotient, remainder) in [
        (7, -3, -3, -2),
        (-7, 3, -3, 2),
        (-7, -3, 2, -1),
        (7, 3, 2, 1),
    ] {
        let (quot, rem) = Integer::from_i64(a)
            .div_rem_floor(&Integer::from_i64(d))
            .unwrap();
        assert_eq!(
            (quot, rem),
            (Integer::from_i64(quotient), Integer::from_i64(remainder))
        );
    }
    let a = Integer::parse(&"1234567890".repeat(1000)).unwrap();
    let b = Integer::parse(&"9876543210".repeat(1000)).unwrap();
    assert_eq!(
        (&a * &b).div_rem_floor(&b),
        Some((a.clone(), Integer::zero()))
    );
    assert_eq!(a.gcd(&Integer::zero()), a);
    assert_eq!(Integer::from_i64(i64::MIN).to_i64(), Some(i64::MIN));
    assert_eq!(q("-7/3").floor(), Integer::from_i64(-3));
    assert_eq!(q("-7/3").ceil(), Integer::from_i64(-2));
}

#[test]
fn double_conversion_roundtrips_and_ties() {
    assert_eq!(
        Rational::from_f64_exact(0.1),
        Some(q("3602879701896397/36028797018963968"))
    );
    assert!(Rational::from_f64_exact(f64::NAN).is_none());
    let mut seed = 0x1234_5678_90ab_cdefu64;
    for _ in 0..4000 {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let f = f64::from_bits(seed);
        if f.is_finite() {
            assert_eq!(
                Rational::from_f64_exact(f)
                    .unwrap()
                    .to_f64_nearest()
                    .to_bits(),
                if f == 0.0 { 0 } else { f.to_bits() }
            );
        }
    }
    let one = Rational::one();
    let half_ulp = one.mul_pow2(-53);
    assert_eq!(
        (&one + &half_ulp).to_f64_nearest().to_bits(),
        1.0f64.to_bits()
    );
    let three_half = &half_ulp * &Rational::from_i64(3);
    assert_eq!(
        (&one + &three_half).to_f64_nearest().to_bits(),
        1.0f64.to_bits() + 2
    );
    assert_eq!(one.mul_pow2(-1075).to_f64_nearest().to_bits(), 0);
    assert_eq!(
        (&one.mul_pow2(-1075) * &q("3")).to_f64_nearest().to_bits(),
        2
    );
    assert!(one.mul_pow2(1024).to_f64_nearest().is_infinite());
}

#[test]
fn directed_decimal_rounding_and_terminating_decimals() {
    for (value, digits, lo, hi) in [
        ("1/3", 5, "0.33333", "0.33334"),
        ("2/3", 5, "0.66666", "0.66667"),
        ("5/4", 2, "1.2", "1.3"),
        ("5/4", 3, "1.25", "1.25"),
        ("-1/3", 3, "-0.334", "-0.333"),
        ("1e30", 3, "1e+30", "1e+30"),
        ("1e-7", 2, "1e-7", "1e-7"),
        ("999/100", 1, "9", "10"),
        ("-999/100", 1, "-10", "-9"),
    ] {
        assert_eq!(
            q(value).to_decimal_bounds(digits),
            (lo.to_owned(), hi.to_owned())
        );
    }
    assert_eq!(q("-1/8").to_decimal_exact(), Some("-0.125".into()));
    assert_eq!(q("100").to_decimal_exact(), Some("100".into()));
    assert_eq!(q("1/3").to_decimal_exact(), None);
}

fn interval_contains(i: Interval, value: &Rational) -> bool {
    (i.lo() == f64::NEG_INFINITY || Rational::from_f64_exact(i.lo()).unwrap() <= *value)
        && (i.hi() == f64::INFINITY || Rational::from_f64_exact(i.hi()).unwrap() >= *value)
}
#[test]
fn interval_arithmetic_contains_exact_results_and_edges() {
    let tenth = Interval::from_rational(&q("1/10"));
    assert_eq!(tenth.lo().next_up().to_bits(), tenth.hi().to_bits());
    assert!(interval_contains(tenth, &q("1/10")));
    let zero = Interval::point(0.0).unwrap();
    assert_eq!(zero * Interval::entire(), zero);
    assert_eq!(zero / zero, Interval::entire());
    assert_eq!(
        Interval::new(-1.0, 4.0).unwrap().sqrt(),
        Interval::new(0.0, 2.0)
    );
    assert!(Interval::new(-2.0, -1.0).unwrap().sqrt().is_none());
    assert_eq!(Interval::entire().mid().to_bits(), 0.0f64.to_bits());
    assert_eq!(
        Interval::new(1.0, f64::INFINITY).unwrap().mid().to_bits(),
        f64::MAX.to_bits()
    );
    for a in -12..=12 {
        for b in 1..=12 {
            let x = Rational::from_ratio(a, 13).unwrap();
            let y = Rational::from_ratio(b, 17).unwrap();
            let ix = Interval::from_rational(&x);
            let iy = Interval::from_rational(&y);
            for (i, r) in [
                (ix + iy, &x + &y),
                (ix - iy, &x - &y),
                (ix * iy, &x * &y),
                (ix / iy, &x / &y),
            ] {
                assert!(interval_contains(i, &r), "{i:?} misses {r}");
            }
        }
    }
    let overflow = Interval::point(1e300).unwrap() * Interval::point(1e300).unwrap();
    assert_eq!(overflow.sign(), Some(StrictSign::Positive));
    assert!(overflow.hi().is_infinite());
}

#[test]
fn ball_rounding_domains_and_sqrt_containment() {
    for p in [0, 1, 10, 53, 128] {
        for n in -20..=20 {
            let r = Rational::from_ratio(n, 7).unwrap();
            let b = Ball::from_rational(&r, p);
            assert!(b.contains(&r));
            assert_eq!(b.precision(), p.max(1));
            for d in [b.mid().denominator().clone(), b.rad().denominator().clone()] {
                let bit = d.bit_length() - 1;
                assert_eq!(d, Integer::one().mul_pow2(bit));
            }
            assert!(b.with_precision(1).contains(&r));
            assert!(b.with_precision(256).contains(&r));
            if n >= 0 {
                let root = b.sqrt().unwrap();
                assert!(
                    root.lo() <= Rational::zero()
                        || &root.lo() * &root.lo() <= b.lo().max(Rational::zero())
                );
                assert!(&root.hi() * &root.hi() >= b.hi());
            }
        }
    }
    assert!(ball("-1", 64).sqrt().is_none());
    assert!(ball("0", 64).ln().is_none());
    assert!(ball("0", 64).recip().is_none());
    assert!(ball("1001/1000", 64).asin().is_none());
    assert!(Ball::atan2(&ball("0", 64), &ball("0", 64)).is_none());
}

#[test]
fn ratio_balls_keep_signed_rounding_and_common_factors() {
    for (n, d, mid, rad) in [
        (7, 2, "4", "1/2"),
        (-7, 2, "-3", "1/2"),
        (7, -2, "-3", "1/2"),
        (-7, -2, "4", "1/2"),
        (0, -7, "0", "0"),
    ] {
        let b = Ball::from_ratio(&Integer::from_i64(n), &Integer::from_i64(d), 2).unwrap();
        assert_eq!(b.mid(), q(mid));
        assert_eq!(b.rad(), q(rad));
    }
    let factor = Integer::one().mul_pow2(1024) + Integer::from_i64(37);
    for bits in [0, 1, 53, 128, 512] {
        for n in -12..=12 {
            for d in [-17, -1, 1, 17] {
                let n = Integer::from_i64(n);
                let d = Integer::from_i64(d);
                let exact = Rational::new(n.clone(), d.clone()).unwrap();
                let b = Ball::from_ratio(&n, &d, bits).unwrap();
                let scaled = Ball::from_ratio(&(&n * &factor), &(&d * &factor), bits).unwrap();
                assert!(b.contains(&exact));
                assert_eq!(b.precision(), bits.max(1));
                assert_eq!(b.mid(), scaled.mid());
                assert_eq!(b.rad(), scaled.rad());
            }
        }
    }
    assert!(Ball::from_ratio(&Integer::one(), &Integer::zero(), 128).is_none());
}

#[test]
fn enclosure_errors_survive_deep_precision_and_underflow() {
    for p in [2048, 8192, 10000] {
        let tiny = q("1/3").mul_pow2(-20000);
        let a = Ball::from_rational(&tiny, p);
        let b = Ball::from_rational(&q("5/7"), p);
        for (enclosure, exact) in [
            (a.clone(), tiny.clone()),
            (a.with_precision(17), tiny.clone()),
            (&a + &b, &tiny + &q("5/7")),
            (&a - &b, &tiny - &q("5/7")),
            (&a * &b, &tiny * &q("5/7")),
            (a.checked_div(&b).unwrap(), &tiny / &q("5/7")),
            (a.square(), &tiny * &tiny),
            (a.recip().unwrap(), tiny.recip().unwrap()),
        ] {
            assert!(enclosure.contains(&exact));
            assert_eq!(enclosure.rad().sign(), Sign::Positive);
        }
        let root = Ball::from_rational(&Rational::one().mul_pow2(-20001), p)
            .sqrt()
            .unwrap();
        assert_eq!(root.rad().sign(), Sign::Positive);
        assert!(&root.lo() * &root.lo() <= Rational::one().mul_pow2(-20001));
        assert!(&root.hi() * &root.hi() >= Rational::one().mul_pow2(-20001));
        let cancellation = &a - &a;
        assert!(cancellation.contains_zero());
        assert_eq!(cancellation.sign(), None);
        let interval = a.to_interval();
        assert_eq!(interval.lo(), 0.0);
        assert_eq!(interval.hi(), f64::from_bits(1));
        assert_eq!(interval.sign(), None);
    }
    let smallest = Interval::point(f64::from_bits(1)).unwrap();
    let half = Interval::point(0.5).unwrap();
    for interval in [smallest * half, smallest / Interval::point(2.0).unwrap()] {
        assert!(interval.lo() <= 0.0);
        assert!(interval.hi() >= f64::from_bits(1));
        assert_eq!(interval.sign(), None);
    }
}

#[test]
fn elementary_functions_match_constant_brackets_and_identities() {
    let p = 256;
    let pi = Ball::pi(p);
    let e = ball("1", p).exp();
    let ln2 = ball("2", p).ln().unwrap();
    for (b, lo, hi) in [
        (
            &pi,
            "3.14159265358979323846264338327950288419716939937510",
            "3.14159265358979323846264338327950288419716939937511",
        ),
        (
            &e,
            "2.71828182845904523536028747135266249775724709369995",
            "2.71828182845904523536028747135266249775724709369996",
        ),
        (
            &ln2,
            "0.69314718055994530941723212145817656807550013436025",
            "0.69314718055994530941723212145817656807550013436026",
        ),
    ] {
        assert!(b.lo() > q(lo) && b.hi() < q(hi));
        assert!(b.rad() < Rational::one().mul_pow2(-240));
    }
    assert!(pi.sin().contains(&q("0")));
    assert!(pi.cos().contains(&q("-1")));
    assert!(ln2.exp().contains(&q("2")));
    assert!(overlaps(&(&ball("1", p).atan() * &ball("4", p)), &pi));
    assert!(overlaps(
        &ball("1", p).asin().unwrap(),
        &(&pi * &ball("1/2", p))
    ));
    assert!(overlaps(
        &Ball::atan2(&ball("1", p), &ball("-1", p)).unwrap(),
        &(&pi * &ball("3/4", p))
    ));
    assert!(ball("0", p).exp().contains(&q("1")));
    assert!(ball("1", p).ln().unwrap().contains(&q("0")));
    assert!(ball("0", p).sin().contains(&q("0")));
    assert!(ball("0", p).cos().contains(&q("1")));
}

#[test]
fn polynomial_arithmetic_factorization_and_resultant() {
    let p = &UPoly::from_roots(&[q("1"), q("1"), q("1")]) * &UPoly::from_roots(&[q("-2"), q("-2")]);
    let factors = p.square_free_factors();
    assert_eq!(
        factors,
        vec![
            (UPoly::from_roots(&[q("-2")]), 2),
            (UPoly::from_roots(&[q("1")]), 3)
        ]
    );
    assert!(p.eval(&q("1")).is_zero());
    let divisor = UPoly::from_roots(&[q("1")]);
    let (quot, rem) = p.div_rem(&divisor).unwrap();
    assert!(rem.is_zero());
    assert_eq!(&quot * &divisor, p);
    assert_eq!(
        UPoly::new(vec![q("1/2"), q("-2/3")]).primitive_integer(),
        vec![Integer::from_i64(-3), Integer::from_i64(4)]
    );
    let x = MPoly::var(2, 0);
    let y = MPoly::var(2, 1);
    assert_eq!((&x - &y).resultant(&(&x + &y), 0), y.scale(&q("2")));
    assert_eq!((&x - &y).substitute(0, &y), MPoly::zero(2));
    let terms = MPoly::from_terms(
        2,
        vec![
            (vec![2, 0], q("1")),
            (vec![0, 1], q("2")),
            (vec![1, 1], q("3")),
            (vec![2, 0], q("-1")),
        ],
    )
    .terms();
    assert_eq!(
        terms
            .iter()
            .map(|(m, _)| m.exponents().to_vec())
            .collect::<Vec<_>>(),
        vec![vec![0, 1], vec![1, 1]]
    );
}

#[test]
fn polynomial_gcd_removes_rational_content_and_keeps_common_factors() {
    let common = UPoly::new(vec![q("1/3"), q("-7/5"), q("1")]);
    let left = &(&common * &UPoly::new(vec![q("1"), q("0"), q("1")]))
        * &UPoly::constant(q("123456789012345678901234567890/37"));
    let right = &(&common * &UPoly::new(vec![q("-2"), q("1")]))
        * &UPoly::constant(q("-98765432109876543210987654321/41"));
    assert_eq!(left.gcd(&right), common);
    assert_eq!(right.gcd(&left), common);
    assert_eq!(UPoly::zero().gcd(&common), common);
    assert_eq!(common.gcd(&UPoly::zero()), common);
    assert!(UPoly::zero().gcd(&UPoly::zero()).is_zero());
}

#[test]
fn radical_zeros_and_nonzero_cancellations_are_certified() {
    let two = real("2").sqrt();
    let three = real("3").sqrt();
    let g = (&Real::one() + &real("5").sqrt()) / real("2");
    let zeros = [
        &(&two * &two) - &real("2"),
        &(&two + &three) - &(&real("5") + &(&real("2") * &real("6").sqrt())).sqrt(),
        &(&real("3") + &(&real("2") * &two)).sqrt() - &(&Real::one() + &two),
        &(&(&g * &g) - &g) - &Real::one(),
        &(&real("2/7").sqrt() * &real("9/11").sqrt()) - &real("18/77").sqrt(),
    ];
    for x in zeros {
        assert!(Decider::default().sign(&x).is_coincident());
        for max_work in [0, 1, 100] {
            let out = Decider::new(BudgetPolicy {
                max_precision_bits: 16,
                max_work,
            })
            .sign(&x);
            assert!(out.is_coincident() || out.is_not_exact(), "{out:?}");
        }
    }
    let tiny =
        &real("10000000000000000000000000000000000000001").sqrt() - &real("100000000000000000000");
    assert!(matches!(
        Decider::default().sign(&tiny),
        Outcome::Proven(StrictSign::Positive)
    ));
    let concave = &(&real("100000000000000000001").sqrt() + &real("99999999999999999999").sqrt())
        - &real("20000000000");
    assert!(matches!(
        Decider::default().sign(&concave),
        Outcome::Proven(StrictSign::Negative)
    ));
}

#[test]
fn squared_algebraic_values_keep_correlated_conjugates() {
    let golden = (&Real::one() + &real("5").sqrt()) / real("2");
    for squared in [golden.square(), golden.pow_u32(2)] {
        // If t²-t-1=0, then s=t² satisfies s²-3s+1=0.
        let residual = &squared.square() - &(&real("3") * &squared) + Real::one();
        assert!(Decider::default().sign(&residual).is_coincident());
        for (offset, expected) in [
            ("1e-40", StrictSign::Positive),
            ("-1e-40", StrictSign::Negative),
        ] {
            assert_eq!(
                Decider::default().sign(&(&residual + &real(offset))),
                Outcome::Proven(expected)
            );
        }
    }
    let even = real("2").sqrt().square() - real("2");
    assert!(Decider::default().sign(&even).is_coincident());
}

#[test]
fn scaled_nested_radical_identity_is_proved_with_default_budget() {
    let k = real("1e-30");
    let radius = &real("2") * &k;
    let offset = (&real("2") * &k.square()).sqrt();
    let direct = (radius.square() - (&radius - &offset).square()).sqrt();
    let factored = &k * &(&real("4") * &real("2").sqrt() - real("2")).sqrt();
    let residual = &direct - &factored;
    assert!(Decider::default().sign(&residual).is_coincident());
    assert_eq!(
        Decider::default().sign(&(&residual + &real("1e-100"))),
        Outcome::Proven(StrictSign::Positive)
    );
}

#[test]
fn case_39_four_radical_expansion_is_exactly_zero() {
    let sum = real("2").sqrt() + real("3").sqrt() + real("5").sqrt() + real("7").sqrt();
    let cross = ["6", "10", "14", "15", "21", "35"]
        .into_iter()
        .fold(Real::zero(), |sum, n| sum + real(n).sqrt());
    let residual = sum.square() - (real("17") + real("2") * cross);
    let outcome = Decider::default().sign(&residual);
    assert!(outcome.is_coincident(), "{outcome:?}");
    assert!(Decider::default().sign(&residual.sqrt()).is_coincident());
    for (offset, sign) in [
        ("1e-100", StrictSign::Positive),
        ("-1e-100", StrictSign::Negative),
    ] {
        assert_eq!(
            Decider::default().sign(&(&residual + &real(offset))),
            Outcome::Proven(sign)
        );
    }
}

#[test]
fn case_39_nested_square_in_an_existing_field_is_exactly_zero() {
    let residual =
        (real("3") + real("2") * real("2").sqrt()).sqrt() - Real::one() - real("2").sqrt();
    let outcome = Decider::default().sign(&residual);
    assert!(outcome.is_coincident(), "{outcome:?}");
}

#[test]
fn case_39_three_radical_sum_has_negative_sign() {
    let residual = real("2").sqrt() + real("3").sqrt() + real("5").sqrt() - real("5.3824");
    assert_eq!(
        Decider::default().sign(&residual),
        Outcome::Proven(StrictSign::Negative)
    );
}

#[test]
fn edge_hunt_0000141_quadratic_root_and_independent_radical_cancel_exactly() {
    let p = UPoly::new(vec![q("145"), q("1116"), q("648")]);
    let Outcome::Proven(roots) = Decider::default().real_roots(&p) else {
        panic!("quadratic root isolation failed");
    };
    let independent = real("5/9").sqrt();
    for root in roots {
        let alpha = Real::from_algebraic(root.value);
        let translated = &alpha + &independent;
        // P(alpha)=0 and independent²=5/9, with no shared expression
        // for the expanded cross terms. This needs both quadratic fields.
        let expanded = real("648") * translated.square()
            - (&real("1296") * &alpha) * independent.clone()
            + &real("1116") * &alpha
            - real("215");
        let outcome = Decider::default().sign(&expanded);
        assert!(outcome.is_coincident(), "{outcome:?}");
        for (offset, sign) in [
            ("1e-100", StrictSign::Positive),
            ("-1e-100", StrictSign::Negative),
        ] {
            assert_eq!(
                Decider::default().sign(&(&expanded + &real(offset))),
                Outcome::Proven(sign)
            );
        }
    }
}

#[test]
fn shifted_irreducible_biquadratic_roots_keep_their_radical_embeddings() {
    let p = UPoly::new(
        ["7508521", "2008080", "-4890600", "-648000", "810000"]
            .map(q)
            .to_vec(),
    );
    let Outcome::Proven(roots) = Decider::default().real_roots(&p) else {
        panic!("quartic root isolation failed");
    };
    assert_eq!(roots.len(), 4);
    let inner = real("1/3").sqrt();
    let low = (real("113/36") - inner.clone()).sqrt();
    let high = (real("113/36") + inner).sqrt();
    for (root, expected) in roots
        .into_iter()
        .zip([-high.clone(), -low.clone(), low, high])
    {
        assert_eq!(root.value.defining_polynomial().degree(), Some(4));
        let actual = Real::from_algebraic(root.value);
        let residual = actual - (real("1/5") + expected);
        assert!(Decider::default().sign(&residual).is_coincident());
        for (offset, sign) in [
            ("1e-100", StrictSign::Positive),
            ("-1e-100", StrictSign::Negative),
        ] {
            assert_eq!(
                Decider::default().sign(&(&residual + &real(offset))),
                Outcome::Proven(sign)
            );
        }
    }
}

#[test]
fn irreducible_quartic_with_nonzero_linear_term_embeds_in_the_tower() {
    let p = UPoly::new(["-1", "-8", "-10", "0", "1"].map(q).to_vec());
    let Outcome::Proven(roots) = Decider::default().real_roots(&p) else {
        panic!("quartic root isolation failed");
    };
    assert_eq!(roots.len(), 4);
    let a = real("2").sqrt();
    let low = (real("3") - a.clone()).sqrt();
    let high = (real("3") + a.clone()).sqrt();
    let expected = [
        -a.clone() - low.clone(),
        a.clone() - high.clone(),
        -a.clone() + low,
        a + high,
    ];
    for (root, expected) in roots.into_iter().zip(expected) {
        let residual = Real::from_algebraic(root.value) - expected;
        assert!(Decider::default().sign(&residual).is_coincident());
        for (offset, sign) in [
            ("1e-100", StrictSign::Positive),
            ("-1e-100", StrictSign::Negative),
        ] {
            assert_eq!(
                Decider::default().sign(&(&residual + &real(offset))),
                Outcome::Proven(sign)
            );
        }
    }
}

#[test]
fn tower_square_roots_and_inverses_keep_the_selected_real_embedding() {
    let a = real("2").sqrt();
    let b = real("3").sqrt();
    for (m, n, negative) in [(-3, 2, true), (-1, -2, true), (1, 1, false), (2, -1, false)] {
        let x = Real::from_i64(m) * a.clone() + Real::from_i64(n) * b.clone();
        let expected = if negative { -x.clone() } else { x.clone() };
        assert!(
            Decider::default()
                .compare(&x.square().sqrt(), &expected)
                .is_coincident()
        );
        let inverse = &Real::one() / &x;
        let residual = &inverse * &x - Real::one();
        assert!(Decider::default().sign(&residual).is_coincident());
    }
    // The first radical creates a genuine quartic tower; a later unrelated
    // generator must not change cached values' field levels.
    let fourth_root = real("2").sqrt().sqrt();
    let sum = &fourth_root + &b;
    let expansion = fourth_root.square() + real("3") + real("2") * (&fourth_root * &b);
    assert!(
        Decider::default()
            .sign(&(sum.pow_u32(2) - expansion))
            .is_coincident()
    );
    let residual = fourth_root.pow_u32(4) - real("2");
    assert!(Decider::default().sign(&residual).is_coincident());
}

#[test]
fn oversized_tower_falls_through_and_leaves_budget_for_refinement() {
    let sum = ["2", "3", "5", "7", "11", "13", "17", "19", "23", "29"]
        .into_iter()
        .fold(real("100000"), |sum, n| sum + real(n).sqrt());
    let residual = &(&sum + &real("1e-104")) - &sum;
    let mut d = Decider::new(BudgetPolicy {
        max_precision_bits: 512,
        max_work: 1_000_000,
    });
    assert_eq!(d.sign(&residual), Outcome::Proven(StrictSign::Positive));
    assert!(d.stats().work <= d.policy().max_work);
}

#[test]
fn translated_tiny_sign_is_preserved_by_exact_cancellation() {
    let mut sum = Real::zero();
    for prime in ["2", "3", "5", "7", "11", "13"] {
        sum = sum + real(prime).sqrt();
    }
    let translated = real("100000") + sum;
    // Both terms are near 10^5; their exact difference is 10^-104.
    let residual = &(&translated + &real("1e-104")) - &translated;
    let mut d = Decider::default();
    assert_eq!(d.sign(&residual), Outcome::Proven(StrictSign::Positive));
    let stats = d.take_stats();
    assert_eq!(stats.decisions, 1);
    assert_eq!(stats.not_exact, 0);
    assert_eq!(stats.refused, 0);
    assert!(stats.work <= d.policy().max_work);
    assert_eq!(
        d.sign(&(-residual.clone())),
        Outcome::Proven(StrictSign::Negative)
    );
}

#[test]
fn quadratic_root_branches_and_positive_square_roots_keep_their_signs() {
    let p = UPoly::new(vec![q("-3"), q("2"), q("7")]);
    let Outcome::Proven(roots) = Decider::default().real_roots(&p) else {
        panic!("quadratic root isolation failed");
    };
    for (root, numerator) in roots.iter().zip([
        -real("2") - real("88").sqrt(),
        -real("2") + real("88").sqrt(),
    ]) {
        let expected = numerator / real("14");
        let actual = Real::from_algebraic(root.value.clone());
        assert!(
            Decider::default()
                .compare(&actual, &expected)
                .is_coincident()
        );
        for (offset, sign) in [
            ("1e-100", StrictSign::Positive),
            ("-1e-100", StrictSign::Negative),
        ] {
            assert_eq!(
                Decider::default().sign(&(actual.clone() - expected.clone() + real(offset))),
                Outcome::Proven(sign)
            );
        }
    }
    let negative = Real::one() - real("2").sqrt();
    let positive_root = negative.square().sqrt();
    assert!(
        Decider::default()
            .compare(&positive_root, &(-negative))
            .is_coincident()
    );
    let scaled = real("12").sqrt() / real("2");
    assert!(
        Decider::default()
            .compare(&scaled, &real("3").sqrt())
            .is_coincident()
    );
}

#[test]
fn independent_quartic_witness_and_tiny_nonzero_offsets_are_separated() {
    let p = UPoly::new(vec![q("-924"), q("0"), q("108"), q("0"), q("1")]);
    let Outcome::Proven(roots) = Decider::default().real_roots(&p) else {
        panic!("quartic root isolation failed");
    };
    let x = real("8") - real("15").sqrt();
    let expected = (real("25") - x.square()).sqrt();
    let actual = Real::from_algebraic(roots[1].value.clone());
    let residual = actual - expected;
    let mut d = Decider::default();
    assert!(d.sign(&residual).is_coincident());
    let first = d.take_stats();
    assert!(d.sign(&residual).is_coincident());
    assert_eq!(first, d.take_stats());
    for (offset, sign) in [
        ("1e-100", StrictSign::Positive),
        ("-1e-100", StrictSign::Negative),
    ] {
        assert_eq!(d.sign(&(&residual + &real(offset))), Outcome::Proven(sign));
    }
}

#[test]
fn generator_reduction_does_not_assume_a_defining_polynomial_is_minimal() {
    let a = UPoly::new(vec![q("-2"), q("0"), q("1")]);
    let b = UPoly::new(vec![q("-3"), q("0"), q("1")]);
    let polynomial = &a * &b;
    let alpha = Real::from_algebraic(Algebraic::isolated(polynomial.clone(), q("1"), q("3/2")));
    let mut d = Decider::default();
    assert!(d.sign(&polynomial.eval_real(&alpha)).is_coincident());
    // This zero has a NONZERO remainder modulo the defining polynomial.
    assert!(d.sign(&a.eval_real(&alpha)).is_coincident());
    assert_eq!(
        d.sign(&b.eval_real(&alpha)),
        Outcome::Proven(StrictSign::Negative)
    );
    // The divisor is nonzero at alpha but zero at another root of P.
    let inverse = Real::one() / b.eval_real(&alpha);
    assert!(d.sign(&(inverse + Real::one())).is_coincident());
    let divisor = alpha.square() + Real::one();
    let quotient = &alpha / &divisor;
    assert!(d.sign(&(&(&quotient * &divisor) - &alpha)).is_coincident());
    assert_eq!(
        d.sign(&(a.eval_real(&alpha) + real("1e-100"))),
        Outcome::Proven(StrictSign::Positive)
    );
}

#[test]
fn undefined_nodes_are_not_removed_by_outer_operations() {
    let zero = &(&real("2").sqrt() * &real("2").sqrt()) - &real("2");
    let bad = Real::one() / zero;
    for x in [&Real::zero() * &bad, bad.pow_u32(0)] {
        match Decider::default().sign(&x) {
            Outcome::Refused(r) => assert_eq!(r.code, "numbers.division-by-zero"),
            other => panic!("{other:?}"),
        }
    }
    let first = &real("-1").sqrt() + &(&Real::one() / &Real::zero());
    match Decider::default().sign(&first) {
        Outcome::Refused(r) => assert_eq!(r.code, "numbers.sqrt-of-negative"),
        other => panic!("{other:?}"),
    }
    assert!(
        Decider::default()
            .sign(&(&(&real("2").sqrt() * &real("2").sqrt()) - &real("2")).sqrt())
            .is_coincident()
    );
}

#[test]
fn constant_polynomials_do_not_introduce_unused_real_dependencies() {
    let undefined = Real::one() / Real::zero();
    assert_eq!(
        UPoly::constant(q("7")).eval_real(&undefined).as_rational(),
        Some(q("7"))
    );
    let p = MPoly::var(2, 0);
    let value = p.eval_real(&[real("3"), undefined]);
    assert_eq!(value.as_rational(), Some(q("3")));
}

#[test]
fn transcendental_series_abandon_in_progress_enclosures() {
    for pi in [true, false] {
        use crate::decider::test_deadline::{Deadline, Point};
        let clock = Deadline::after(Point::Enclosure, 4);
        let mut dec = Decider::default();
        dec.set_deadline(Some(clock.instant()));
        let result = if pi {
            Ball::pi_checked(128, &dec)
        } else {
            ball("1/3", 128).atan_checked(&dec)
        };
        assert_eq!(result.unwrap_err().code, "numbers.deadline");
        clock.assert_expired();
    }
}

#[test]
fn root_isolation_abandons_in_progress_decisions() {
    use crate::decider::test_deadline::{Deadline, Point};
    let polynomial = UPoly::from_roots(&(1..=8).map(Rational::from_i64).collect::<Vec<_>>());
    let mut dec = Decider::new(BudgetPolicy {
        max_precision_bits: 16_384,
        max_work: u64::MAX,
    });
    let clock = Deadline::after(Point::Work, 8);
    dec.set_deadline(Some(clock.instant()));
    assert!(matches!(
        dec.real_roots(&polynomial),
        Outcome::Refused(Refusal {
            code: "numbers.deadline",
            ..
        })
    ));
    clock.assert_expired();
}

#[test]
fn resultant_abandons_in_progress_exact_elimination() {
    let f = MPoly::from_terms(
        2,
        (0..=16)
            .map(|i| (vec![16 - i, i], Rational::from_i64(i64::from(i * i + 1))))
            .collect(),
    );
    let g = MPoly::from_terms(
        2,
        (0..=16)
            .map(|i| {
                (
                    vec![16 - i, i],
                    Rational::from_i64(i64::from(i * i * i + 3)),
                )
            })
            .collect(),
    );
    let dec = Decider::default();
    // The entry and 32 Sylvester-row checks must pass before expiry can
    // occur inside Bareiss elimination. The existing callback needs no hook.
    let checks = std::cell::Cell::new(0);
    let dec = std::cell::RefCell::new(dec);
    let result = f.resultant_checked(&g, 1, &|| {
        let count = checks.get() + 1;
        checks.set(count);
        if count == 64 {
            dec.borrow_mut().set_deadline(Some(Instant::now()));
        }
        dec.borrow().check_deadline().map_err(|r| r.code)
    });
    assert_eq!(result.unwrap_err(), "numbers.deadline");
    assert_eq!(checks.get(), 64);
}

#[test]
fn precision_enclosures_refine_rational_radical_and_algebraic_values() {
    let algebraic = Real::from_algebraic(Algebraic::isolated(
        UPoly::new(vec![q("-2"), q("0"), q("0"), q("1")]),
        q("1"),
        q("2"),
    ));
    let mut dec = Decider::default();
    assert!(dec.sign(&real("1")).is_proven());
    let stats = dec.stats().clone();
    for (value, power, expected) in [
        (real("1/3"), 1, q("1/3")),
        (real("2").sqrt(), 2, q("2")),
        (algebraic, 3, q("2")),
    ] {
        let mut previous_width = None;
        for precision in [32, 64, 128, 256] {
            let enclosure = dec.enclose_at(&value, precision).unwrap();
            assert_eq!(enclosure.precision(), precision);
            assert!(enclosure.lo().sign() == Sign::Positive);
            assert!(enclosure.lo().pow(power).unwrap() <= expected);
            assert!(enclosure.hi().pow(power).unwrap() >= expected);
            let width = enclosure.rad().mul_pow2(1);
            if let Some(previous) = previous_width {
                assert!(width < previous);
            }
            previous_width = Some(width);
            assert_eq!(dec.stats(), &stats);
        }
        dec.set_deadline(Some(Instant::now()));
        assert_eq!(
            dec.enclose_at(&value, 64).unwrap_err().code,
            "numbers.deadline"
        );
        assert_eq!(dec.stats(), &stats);
        dec.set_deadline(None);
    }
}

#[test]
fn precision_enclosures_clip_known_square_roots_and_retry_divisors() {
    let dec = Decider::new(BudgetPolicy {
        max_precision_bits: 32,
        max_work: 0,
    });
    let root = real("2").sqrt();
    let zero = &root.square() - &real("2");
    let enclosed = dec.enclose_at(&zero.sqrt(), 64).unwrap();
    assert!(enclosed.contains_zero());
    let divisor = &(&root + &real("1e-100")) - &root;
    let nested = dec.enclose_at(&zero.sqrt().sqrt(), 64).unwrap();
    assert!(nested.contains_zero());
    let positive = dec.enclose_at(&divisor.sqrt(), 64).unwrap();
    assert!(positive.contains(&q("1e-50")));
    let reciprocal = &Real::one() / &divisor;
    assert_eq!(
        dec.enclose_at(&reciprocal, 64).unwrap_err().code,
        "numbers.divisor-not-separated"
    );
    let precise = dec.enclose_at(&reciprocal, 512).unwrap();
    assert!(precise.contains(&q("1e100")));
    assert!(precise.lo().sign() == Sign::Positive);
    assert_eq!(dec.stats(), &DecisionStats::default());
}

#[test]
fn precision_enclosures_abandon_in_progress_without_accounting_decisions() {
    use crate::decider::test_deadline::{Deadline, Point};
    let root = real("2").sqrt();
    let expression = (&root + &real("3").sqrt()).square();
    let mut dec = Decider::default();
    let clock = Deadline::after(Point::Work, 4);
    dec.set_deadline(Some(clock.instant()));
    assert_eq!(
        dec.enclose_at(&expression, 128).unwrap_err().code,
        "numbers.deadline"
    );
    clock.assert_expired();
    assert_eq!(dec.stats(), &DecisionStats::default());
}

#[test]
fn decision_accounting_repeats_and_deadlines_do_no_work() {
    let a = real("2").sqrt();
    let b = real("2").sqrt();
    let mut d = Decider::default();
    let first = d.sign(&(&a - &b));
    let stats = d.take_stats();
    let second = d.sign(&(&a - &b));
    assert_eq!(first, second);
    assert_eq!(stats, d.take_stats());
    let mut compare = Decider::default();
    let out = compare.compare(&a, &b);
    assert_eq!(
        out,
        first.map(|s| if s == StrictSign::Negative {
            StrictOrder::Less
        } else {
            StrictOrder::Greater
        })
    );
    assert_eq!(compare.take_stats(), stats);
    d.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
    assert!(matches!(
        d.sign(&a),
        Outcome::Refused(Refusal {
            code: "numbers.deadline",
            ..
        })
    ));
    assert!(matches!(
        d.real_roots(&UPoly::x()),
        Outcome::Refused(Refusal {
            code: "numbers.deadline",
            ..
        })
    ));
    assert_eq!(d.stats().work, 0);
    assert_eq!(d.stats().refused, 2);
    let mut d = Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        max_work: 0,
    });
    let out = d.sign(&real("5"));
    match out {
        Outcome::NotExact(n) => assert!(n.bound >= q("5")),
        other => panic!("{other:?}"),
    }
    assert_eq!(d.stats().work, 0);
}

#[test]
fn real_roots_are_sorted_exact_when_rational_and_refinable() {
    let p = UPoly::from_roots(&[q("7/5"), q("-2"), q("1/3"), q("-2")]);
    let Outcome::Proven(roots) = Decider::default().real_roots(&p) else {
        panic!("root isolation refused");
    };
    assert_eq!(
        roots
            .iter()
            .map(|r| (r.value.as_rational().unwrap().clone(), r.multiplicity))
            .collect::<Vec<_>>(),
        vec![(q("-2"), 2), (q("1/3"), 1), (q("7/5"), 1)]
    );
    let p = UPoly::new(vec![q("-2"), q("0"), q("1")]);
    let Outcome::Proven(roots) = Decider::default().real_roots(&p) else {
        panic!("root isolation refused");
    };
    assert_eq!(roots.len(), 2);
    for root in roots {
        assert!(root.value.as_rational().is_none());
        let refined = root.value.refined(&q("1e-30"));
        let (lo, hi) = refined.isolating_interval();
        assert!(&hi - &lo <= q("1e-30"));
        assert_ne!(p.eval(&lo).sign(), p.eval(&hi).sign());
        assert!(!p.eval(&lo).is_zero());
        assert!(!p.eval(&hi).is_zero());
    }
}

#[test]
fn system_certification_separates_unique_empty_and_singular_boxes() {
    let x = MPoly::var(2, 0);
    let y = MPoly::var(2, 1);
    let one = MPoly::constant(2, q("1"));
    let circle = &(&x * &x) + &(&(&y * &y) - &one);
    let system = Arc::new(PolySystem::new(vec![circle.clone(), &x - &y]).unwrap());
    let mut d = Decider::default();
    let positive = RootBox::new(vec![q("1/2"), q("1/2")], vec![q("1"), q("1")]).unwrap();
    let Outcome::Proven(BoxRoots::One(root)) = d.certify_root(&system, &positive) else {
        panic!("unique root not certified");
    };
    assert!(
        d.compare(&root.coordinate(0), &real("1/2").sqrt())
            .is_coincident()
    );
    let mut capped = Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        max_work: 1_000_000,
    });
    assert!(
        capped
            .compare(&root.coordinate(0), &real("1/2").sqrt())
            .is_coincident()
    );
    assert!(capped.stats().max_bits <= 64);
    assert!(capped.stats().work <= 1_000_000);
    let empty = RootBox::new(vec![q("-1/10"), q("-1/10")], vec![q("1/10"), q("1/10")]).unwrap();
    assert!(matches!(
        d.certify_root(&system, &empty),
        Outcome::Proven(BoxRoots::None)
    ));
    let both = RootBox::new(vec![q("-1"), q("-1")], vec![q("1"), q("1")]).unwrap();
    assert!(d.certify_root(&system, &both).is_refused());
    let tangent = Arc::new(PolySystem::new(vec![circle, &y - &one]).unwrap());
    let bounds = RootBox::new(vec![q("-1/10"), q("9/10")], vec![q("1/10"), q("11/10")]).unwrap();
    assert!(d.certify_root(&tangent, &bounds).is_refused());
}

#[test]
fn facts_keep_origin_order_and_relation_specific_closures() {
    let mut book = FactBook::new();
    let make = |kind, a, b| Fact {
        kind,
        a,
        b,
        at: None,
        origin: FactOrigin::Construction("unit"),
    };
    assert!(book.record(make(FactKind::Coaxial, 1, 2)));
    assert!(book.record(make(FactKind::Coaxial, 2, 3)));
    assert!(book.holds(FactKind::Coaxial, 3, 1));
    assert!(book.holds(FactKind::Coaxial, 99, 99));
    assert!(!book.record(Fact {
        origin: FactOrigin::Proof(ZeroProof::Algebraic),
        ..make(FactKind::Coaxial, 2, 1)
    }));
    book.record(make(FactKind::Coincident, 10, 11));
    book.record(make(FactKind::Parallel, 11, 12));
    assert!(!book.holds(FactKind::Parallel, 10, 12));
    book.record(make(FactKind::Incident, 20, 21));
    book.record(make(FactKind::Coincident, 21, 22));
    book.record(make(FactKind::Incident, 22, 23));
    assert!(book.holds(FactKind::Incident, 20, 23));
    assert!(!book.holds(FactKind::Incident, 23, 20));
    book.record(make(FactKind::Tangent, 30, 31));
    book.record(make(FactKind::Tangent, 31, 32));
    assert!(!book.holds(FactKind::Tangent, 30, 32));
    assert!(!book.holds(FactKind::Tangent, 30, 30));
    assert_eq!(
        book.facts_about(2).iter().map(|f| f.a).collect::<Vec<_>>(),
        vec![1, 2]
    );
}

#[test]
fn negated_construction_identity_keeps_signs_and_denominators() {
    for value in ["0", "1", "-1", "12345678901234567890123456789/17"] {
        let a = real(value);
        let b = -&a;
        assert!(a.same_negated_definition(&b));
        assert!(b.same_negated_definition(&a));
    }
    for (a, b) in [
        ("1", "1"),
        ("1/2", "-1/3"),
        ("1", "-1.0000000000000000000001"),
    ] {
        assert!(!real(a).same_negated_definition(&real(b)));
    }
    let root = real("2").sqrt();
    let opposite = -&root;
    assert!(root.same_negated_definition(&opposite));
    assert!(opposite.same_negated_definition(&root));
    assert!(!root.same_negated_definition(&root));
    let mut dec = Decider::default();
    assert!(dec.sign(&(&root + &opposite)).is_coincident());
}

#[test]
fn mark_bounds_and_answer_mapping_preserve_errors() {
    let mut marks = MarkLog::new();
    marks.push(Mark::new(
        "test.first",
        BoundKind::Distance,
        q("1/7"),
        "first",
    ));
    marks.push(Mark::new("test.second", BoundKind::Angle, q("10"), "angle"));
    let mut extra = MarkLog::new();
    extra.push(Mark::new(
        "test.third",
        BoundKind::Distance,
        q("1/3"),
        "third",
    ));
    marks.extend(&extra);
    assert_eq!(marks.worst(BoundKind::Distance), Some(&q("1/3")));
    assert_eq!(marks.worst(BoundKind::Value), None);
    let answer = Answer::with_marks(7, marks).map(|x| x + 2);
    assert_eq!(answer.value, 9);
    assert!(!answer.is_exact());
    assert_eq!(
        answer.marks.iter().map(|m| m.code).collect::<Vec<_>>(),
        vec!["test.first", "test.second", "test.third"]
    );
}

#[test]
fn deep_shared_defined_rational_zero_is_proven() {
    let mut zero = Real::zero();
    for _ in 0..64 {
        zero = &zero + &zero;
    }
    assert!(Decider::default().sign(&zero).is_coincident());
    let undefined = &zero * &(Real::one() / Real::zero());
    assert!(matches!(
        Decider::default().sign(&undefined),
        Outcome::Refused(Refusal {
            code: "numbers.division-by-zero",
            ..
        })
    ));
}

#[test]
fn singleton_zero_is_distinct_from_a_ball_that_only_contains_zero() {
    let zero = ball("0", 64);
    let tiny = ball("1e-1000", 64);
    assert!(zero.is_exact_zero());
    for uncertain in [zero.hull(&tiny), tiny.hull(&(-&tiny))] {
        assert!(uncertain.contains_zero());
        assert!(!uncertain.is_exact_zero());
    }
    assert!(!tiny.is_exact_zero());
}

#[test]
fn open_root_endpoints_prove_strict_sign_without_refinement() {
    let leading = &q("1e600") + &q("1");
    for (constant, lo, hi, expected) in [
        ("-2", "0", "1", StrictSign::Positive),
        ("2", "-1", "0", StrictSign::Negative),
    ] {
        // Each cubic is strictly increasing and changes sign between these
        // endpoints. Its sole root is nonzero despite the zero endpoint.
        let value = Real::from_algebraic(Algebraic::isolated(
            UPoly::new(vec![q(constant), q("0"), q("0"), leading.clone()]),
            q(lo),
            q(hi),
        ));
        let mut decider = Decider::new(BudgetPolicy {
            max_precision_bits: 0,
            max_work: 1,
        });
        assert!(matches!(decider.sign(&value), Outcome::Proven(sign) if sign == expected));
    }
    let zero = Real::from_algebraic(Algebraic::isolated(
        UPoly::new(vec![q("0"), q("1")]),
        q("-1"),
        q("1"),
    ));
    assert!(Decider::default().sign(&zero).is_coincident());
}

fn odd_denominator_cubic_root() -> Real {
    // The cubic is strictly increasing, and these exact endpoints bracket
    // cbrt(2). A large odd denominator is certificate representation, not
    // mathematical uncertainty that requires huge exact Horner products.
    let denominator = Integer::one().mul_pow2(4096) - Integer::one();
    let epsilon = Rational::new(Integer::one(), denominator).unwrap();
    let root = Algebraic::isolated(
        UPoly::new(vec![q("-2"), q("0"), q("0"), q("1")]),
        &Rational::one() + &epsilon,
        &Rational::from_i64(2) - &epsilon,
    );
    Real::from_algebraic(root)
}

#[test]
fn high_denominator_root_brackets_refine_with_a_small_work_budget() {
    let value = odd_denominator_cubic_root();
    let mut dec = Decider::new(BudgetPolicy {
        max_precision_bits: 512,
        max_work: 1_000_000,
    });
    let width = q("1e-30");
    let enclosure = dec.enclose(&value, &width).unwrap();
    assert!(enclosure.rad().mul_pow2(1) <= width);
    let oracle = |r: Rational| {
        rug::Rational::from((
            r.numerator().to_string().parse::<rug::Integer>().unwrap(),
            r.denominator().to_string().parse::<rug::Integer>().unwrap(),
        ))
    };
    let lo = oracle(enclosure.lo());
    let hi = oracle(enclosure.hi());
    assert!(lo.clone().square() * lo < 2);
    assert!(hi.clone().square() * hi > 2);
}

#[test]
fn filtered_bisection_prefixes_preserve_cold_work_at_budget_boundaries() {
    let warmed = odd_denominator_cubic_root();
    for precision in [128, 512, 768] {
        Decider::default().enclose_at(&warmed, precision).unwrap();
    }
    let target = q("1e-220");
    let fingerprint = |result: Result<Ball, Refusal>| {
        result.map(|ball| (ball.mid(), ball.rad(), ball.precision()))
    };
    for max_work in [1_000_000, BudgetPolicy::default().max_work] {
        let policy = BudgetPolicy {
            max_precision_bits: 768,
            max_work,
        };
        let mut cold = Decider::new(policy);
        let mut warm = Decider::new(policy);
        let cold_result = cold.enclose(&odd_denominator_cubic_root(), &target);
        let warm_result = warm.enclose(&warmed, &target);
        if max_work == BudgetPolicy::default().max_work {
            assert!(cold_result.as_ref().unwrap().rad().mul_pow2(1) <= target);
        }
        assert_eq!(fingerprint(cold_result), fingerprint(warm_result));
        assert_eq!(cold.stats(), warm.stats());
    }
}

#[test]
fn singleton_zero_certificates_do_not_require_elimination() {
    let zero = Real::zero() * real("2").sqrt();
    let policy = |max_work| BudgetPolicy {
        max_precision_bits: 64,
        max_work,
    };
    match Decider::new(policy(0)).sign(&zero) {
        Outcome::NotExact(n) => assert!(n.bound.is_zero()),
        other => panic!("an empty allowance must still stop: {other:?}"),
    }
    assert!(matches!(
        Decider::new(policy(1)).sign(&zero),
        Outcome::ProvenCoincident(_)
    ));

    let tiny = real("1e-1000") * real("2").sqrt();
    match Decider::new(policy(1)).sign(&tiny) {
        Outcome::Proven(StrictSign::Positive) => {}
        Outcome::NotExact(n) => assert!(n.bound.sign() == Sign::Positive),
        other => panic!("a nonzero enclosure cannot prove coincidence: {other:?}"),
    }
    let gap = real("2").sqrt() - real("1.414213562373095048801688724209698");
    let delayed = Real::zero() * (Real::one() / gap);
    assert!(matches!(
        Decider::default().sign(&delayed),
        Outcome::ProvenCoincident(_)
    ));
    // The outer zero has a singleton coarse ball, but its hidden divisor is
    // undefined. Definedness must precede every enclosure-based zero proof.
    let undefined = &zero * &(Real::one() / Real::zero()).pow_u32(0);
    assert!(matches!(
        Decider::default().sign(&undefined),
        Outcome::Refused(Refusal {
            code: "numbers.division-by-zero",
            ..
        })
    ));
}

#[test]
fn shared_algebraic_rational_identities_preserve_domains_and_small_gaps() {
    let [sine, cosine] = sin_cos_degrees(&q("82"), &mut Decider::default()).unwrap();
    let slope = sine / (Real::one() + cosine);
    let radius = real("3");
    let depth = &radius / &slope;
    // A formal rational identity holds for every defined nonzero slope, without
    // eliminating its high-degree algebraic generator or approximating the gap.
    let reconstructed = &slope * &(&depth + &real("4") - real("4"));
    let zero = &reconstructed - &radius;
    let (answer, work) = with_operation_budget(100_000_000, || Decider::default().sign(&zero));
    assert!(answer.is_coincident(), "{answer:?} {work:?}");
    assert!(!work.exhausted);
    for epsilon in ["1e-30", "-1e-30"] {
        let result = Decider::default().sign(&(&zero + &real(epsilon)));
        assert_eq!(
            result.proven(),
            Some(&if epsilon.starts_with('-') {
                StrictSign::Negative
            } else {
                StrictSign::Positive
            }),
            "{result:?}"
        );
    }
    for undefined in [
        &zero / &zero,
        &zero * &(Real::one() / Real::zero()).pow_u32(0),
        &zero * &real("-1").sqrt(),
    ] {
        assert!(matches!(
            Decider::default().sign(&undefined),
            Outcome::Refused(_)
        ));
    }
    let (answer, work) = with_operation_budget(0, || Decider::default().sign(&zero));
    assert!(matches!(answer, Outcome::Refused(_)));
    assert!(work.exhausted);
}

#[test]
fn algebraic_coefficient_quadratics_keep_order_multiplicity_and_domains() {
    let [sine, cosine] = sin_cos_degrees(&q("82"), &mut Decider::default()).unwrap();
    let a = sine / (Real::one() + cosine);
    for scale in [Real::one(), real("-3")] {
        let polynomial = RealPoly::from_terms(
            1,
            vec![(vec![2], &scale * &a), (vec![0], -&(&scale * &a.square()))],
        );
        let (roots, work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::default();
            let roots = polynomial.real_roots(&mut dec).unwrap();
            assert_eq!(roots.len(), 2);
            assert_eq!(dec.sign(&roots[0]).proven(), Some(&StrictSign::Negative));
            assert_eq!(dec.sign(&roots[1]).proven(), Some(&StrictSign::Positive));
            for root in &roots {
                let outcome = dec.sign(&polynomial.eval_real(std::slice::from_ref(root)));
                assert!(
                    outcome.is_coincident(),
                    "{outcome:?} {:?}",
                    current_operation_work()
                );
            }
            roots
        });
        assert_eq!(roots.len(), 2);
        assert!(!work.exhausted);
    }
    for (terms, count) in [
        (
            vec![(vec![1], a.clone()), (vec![0], -&(&a * &real("3")))],
            1,
        ),
        (
            vec![
                (vec![2], a.clone()),
                (vec![1], -&(&a * &real("6"))),
                (vec![0], &a * &real("9")),
            ],
            1,
        ),
        (vec![(vec![2], a.clone()), (vec![0], a.clone())], 0),
        (vec![(vec![0], a.clone())], 0),
    ] {
        let polynomial = RealPoly::from_terms(1, terms);
        let roots = polynomial.real_roots(&mut Decider::default()).unwrap();
        assert_eq!(roots.len(), count);
        if count == 1 {
            assert!(
                Decider::default()
                    .sign(&(&roots[0] - &real("3")))
                    .is_coincident()
            );
        }
    }
    let polynomial = RealPoly::from_terms(1, vec![(vec![1], a), (vec![0], real("1"))]);
    let (result, work) =
        with_operation_budget(0, || polynomial.real_roots(&mut Decider::default()));
    assert!(result.is_err());
    assert!(work.exhausted);
    let undefined =
        RealPoly::from_terms(1, vec![(vec![1], real("-1").sqrt()), (vec![0], real("1"))]);
    assert!(undefined.real_roots(&mut Decider::default()).is_err());
}
