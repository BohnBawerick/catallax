#![allow(missing_docs, clippy::float_cmp)]
mod common;
use catallax_numbers::*;
use common::*;
use rug::ops::Pow;
use rug::{Integer as Z, Rational as Q};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

blind_test!(parsing_grammar_and_error_precedence, {
    for (text, n, d) in [
        ("1.1", 11, 10),
        ("0.1", 1, 10),
        ("1e-7", 1, 10_000_000),
        ("2.5E3", 2500, 1),
        ("-6/8", -3, 4),
        ("-0", 0, 1),
        ("007", 7, 1),
        ("+000.125e+1", 5, 4),
        ("+6/008", 3, 4),
    ] {
        assert_q(&r(text), &q(n, d));
        assert_q(&text.parse::<Rational>().unwrap(), &q(n, d));
    }
    for text in [
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
        "1e+",
        "--1",
        "1/",
        "1/0e1",
        "１",
        "1\n",
        "1\0",
    ] {
        assert_eq!(
            Rational::parse(text).unwrap_err().kind,
            ParseErrorKind::InvalidSyntax,
            "{text:?}"
        );
    }
    assert_eq!(Rational::parse("").unwrap_err().kind, ParseErrorKind::Empty);
    for text in ["1/0", "0/000", "-5/000", "+0/0"] {
        assert_eq!(
            Rational::parse(text).unwrap_err().kind,
            ParseErrorKind::ZeroDenominator
        );
    }
    for text in [
        "1e100001",
        "1e-100001",
        "0e999999999999999999999999999999",
        "1e+000100001",
    ] {
        assert_eq!(
            Rational::parse(text).unwrap_err().kind,
            ParseErrorKind::ExponentTooLarge
        );
    }
    for (text, e) in [
        ("1e100000", 100000i32),
        ("1e-100000", -100000),
        ("1e+000100000", 100000),
    ] {
        let ten = Z::from(10).pow(e.unsigned_abs());
        let expected = if e >= 0 {
            Q::from(ten)
        } else {
            Q::from((Z::from(1), ten))
        };
        assert_q(&r(text), &expected);
    }
    for text in ["1.0", "1/2", "1e2", " 1", "+", "0x1", "1_0", "1/-2"] {
        assert_eq!(
            Integer::parse(text).unwrap_err().kind,
            ParseErrorKind::InvalidSyntax
        );
    }
    assert_eq!(Integer::parse("").unwrap_err().kind, ParseErrorKind::Empty);
    assert_eq!(Integer::parse("-0007").unwrap().to_string(), "-7");
});

blind_test!(normal_form_field_and_zero_edges, {
    let hash = |x: &Rational| {
        let mut h = DefaultHasher::new();
        x.hash(&mut h);
        h.finish()
    };
    for (a, b) in [
        (r("6/4"), r("3/2")),
        (r("0/5"), r("0")),
        (Rational::from_ratio(-6, -4).unwrap(), r("3/2")),
    ] {
        assert_eq!(a, b);
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Equal);
        assert_eq!(hash(&a), hash(&b));
    }
    assert_q(&(r("0.1") + r("0.2")), &q(3, 10));
    assert_q(
        &Rational::new(Integer::from(-6), Integer::from(-4)).unwrap(),
        &q(3, 2),
    );
    for (a, b) in [
        ("-2", "-1/3"),
        ("-1/3", "0"),
        ("1/2", "2/3"),
        ("100/101", "1"),
    ] {
        assert!(r(a) < r(b));
    }
    assert_eq!(r("-6/8").to_string(), "-3/4");
    assert_eq!(r("0/5").to_string(), "0");
    assert!(Rational::from_ratio(1, 0).is_none());
    assert!(r("0").recip().is_none());
    assert!(r("1").checked_div(&r("0")).is_none());
    assert!(r("0").pow(-1).is_none());
    assert_q(&r("0").pow(0).unwrap(), &q(1, 1));
    assert_eq!(Integer::zero().pow(0), Integer::one());
    assert_q(&r("-2/3").pow(-3).unwrap(), &q(-27, 8));
    assert_eq!(r("-7/3").floor().to_string(), "-3");
    assert_eq!(r("-7/3").ceil().to_string(), "-2");
    assert!(std::panic::catch_unwind(|| r("1") / r("0")).is_err());
});

blind_test!(generated_256_bit_field_laws, {
    let mut rng = Rng(0xf1e1d);
    for _ in 0..1024 {
        let (a, b, c) = (rng.wide(), rng.wide(), rng.wide());
        let (x, y, z) = (eq(&a), eq(&b), eq(&c));
        assert_q(&(&x + &y), &(a.clone() + &b));
        assert_q(&(&x - &y), &(a.clone() - &b));
        assert_q(&(&x * &y), &(a.clone() * &b));
        assert_eq!((&x + &y) + z.clone(), x.clone() + (&y + &z));
        assert_eq!((&x * &y) * z.clone(), x.clone() * (&y * &z));
        assert_eq!(&x * &(&y + &z), &(&x * &y) + &(&x * &z));
        if b != 0 {
            assert_q(&(&x / &y), &(a.clone() / &b));
            assert_eq!((&x / &y) * y.clone(), x);
        }
        assert_q(&(&x - &x), &q(0, 1));
    }
});

#[test]
fn dyadic_normalization_and_cancellation_match_gmp() {
    for bits in [1, 53, 2048, 10000] {
        let denominator = Integer::one().mul_pow2(bits);
        for n in [-24, -3, 0, 3, 24] {
            for sign in [-1, 1] {
                let a =
                    Rational::new(Integer::from(n), &denominator * &Integer::from(sign)).unwrap();
                let expected = Q::from((Z::from(n * sign), Z::from(1) << bits as u32));
                assert_q(&a, &expected);
                for shift in [-(bits as i64) - 7, -1, 0, 1, bits as i64 + 7] {
                    let scale = Q::from((Z::from(1) << shift.unsigned_abs() as u32, 1));
                    let scaled = if shift < 0 {
                        expected.clone() / scale
                    } else {
                        expected.clone() * scale
                    };
                    assert_q(&a.mul_pow2(shift), &scaled);
                }
                let b = Rational::from_ratio(5, 8).unwrap();
                for (value, oracle) in [
                    (&a + &b, expected.clone() + q(5, 8)),
                    (&a - &b, expected.clone() - q(5, 8)),
                    (&a * &b, expected.clone() * q(5, 8)),
                    (&a / &b, expected.clone() / q(5, 8)),
                    (&a - &a, Q::from(0)),
                ] {
                    assert_q(&value, &oracle);
                }
            }
        }
    }
}

#[test]
fn shared_odd_denominators_and_reciprocals_remain_canonical() {
    for bits in [53, 607, 2048] {
        let shared = (Integer::one().mul_pow2(bits) - Integer::one()) * Integer::from(11);
        let a = Rational::new(Integer::one(), &shared * &Integer::from(3)).unwrap();
        let b = Rational::new(Integer::from(2), &shared * &Integer::from(5)).unwrap();
        for b in [b.clone(), -b, a.clone(), -&a] {
            let x = oq(&a);
            let y = oq(&b);
            assert_q(&(&a + &b), &(x.clone() + &y));
            assert_q(&(&a - &b), &(x.clone() - &y));
            assert_q(&b.recip().unwrap(), &(Q::from(1) / &y));
        }
    }
}

blind_test!(integer_floor_division_gcd_and_ten_thousand_digits, {
    // The technical review resolution: q=floor(a/d); nonzero remainder has d's sign.
    for a in -20i64..=20 {
        for d in -9i64..=9 {
            if d == 0 {
                continue;
            }
            let mut quotient = a / d;
            if a % d != 0 && (a < 0) != (d < 0) {
                quotient -= 1;
            }
            let (got_q, got_r) = Integer::from_i64(a)
                .div_rem_floor(&Integer::from_i64(d))
                .unwrap();
            assert_eq!(got_q.to_i64(), Some(quotient));
            assert_eq!(got_r.to_i64(), Some(a - quotient * d));
        }
    }
    assert!(Integer::one().div_rem_floor(&Integer::zero()).is_none());
    assert_eq!(Integer::zero().gcd(&Integer::zero()), Integer::zero());
    let sa = "9".repeat(10_000);
    let sb = format!("1{}7", "0".repeat(9_998));
    let (a, b) = (Integer::parse(&sa).unwrap(), Integer::parse(&sb).unwrap());
    let (oa, ob) = (z(&sa), z(&sb));
    let product = &a * &b;
    let (quotient, remainder) = product.div_rem_floor(&b).unwrap();
    assert_eq!(quotient, a);
    assert!(remainder.is_zero());
    assert_eq!(a.gcd(&b).to_string(), oa.clone().gcd(&ob).to_string());
    assert_eq!(product.to_string(), (oa * ob).to_string());
    assert_eq!(Integer::one().mul_pow2(1100).bit_length(), 1101);
    assert_eq!(Integer::zero().bit_length(), 0);
    assert_eq!(Integer::from_i64(i64::MIN).to_i64(), Some(i64::MIN));
    assert_eq!(Integer::from_u64(u64::MAX).to_i64(), None);
});

blind_test!(binary_conversion_roundtrips_and_ties, {
    assert_q(
        &Rational::from_f64_exact(0.1).unwrap(),
        &Q::from((z("3602879701896397"), z("36028797018963968"))),
    );
    for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(Rational::from_f64_exact(v).is_none());
    }
    let mut bits = vec![
        0,
        1,
        2,
        0x8000000000000000,
        0x8000000000000001,
        0x000fffffffffffff,
        0x0010000000000000,
        f64::MAX.to_bits(),
        (-f64::MAX).to_bits(),
    ];
    let mut rng = Rng(0xb1a4);
    bits.extend((0..10000).map(|_| rng.next()));
    for bits in bits {
        let v = f64::from_bits(bits);
        if !v.is_finite() {
            continue;
        }
        let got = Rational::from_f64_exact(v).unwrap();
        assert_q(&got, &Q::from_f64(v).unwrap());
        assert_eq!(
            got.to_f64_nearest().to_bits(),
            if v == 0.0 { 0 } else { bits }
        );
    }
    // Adjacent doubles at 1 are 2^-52 apart. Midpoints round to the even significand.
    let unit = Q::from((Z::from(1), Z::from(1) << 53));
    assert_eq!(
        eq(&(q(1, 1) + &unit)).to_f64_nearest().to_bits(),
        1.0f64.to_bits()
    );
    assert_eq!(
        eq(&(q(1, 1) + unit * 3)).to_f64_nearest().to_bits(),
        1.0f64.to_bits() + 2
    );
    let half_subnormal = Q::from((Z::from(1), Z::from(1) << 1075));
    assert_eq!(eq(&half_subnormal).to_f64_nearest().to_bits(), 0);
    assert_eq!(eq(&(half_subnormal * 3)).to_f64_nearest().to_bits(), 2);
    assert_eq!(r("1e400").to_f64_nearest(), f64::INFINITY);
    assert_eq!(r("-1e400").to_f64_nearest(), f64::NEG_INFINITY);
});

blind_test!(decimal_exact_and_tight_bounds, {
    for (x, p, lo, hi) in [
        ("1/3", 5, "0.33333", "0.33334"),
        ("2/3", 5, "0.66666", "0.66667"),
        ("5/4", 2, "1.2", "1.3"),
        ("5/4", 3, "1.25", "1.25"),
        ("-1/3", 3, "-0.334", "-0.333"),
        ("1e30", 3, "1e+30", "1e+30"),
        ("1e-7", 2, "1e-7", "1e-7"),
        ("1e-6", 1, "0.000001", "0.000001"),
        ("1e20", 1, "100000000000000000000", "100000000000000000000"),
        ("1e21", 1, "1e+21", "1e+21"),
        ("999/100", 2, "9.9", "10"),
        ("-999/100", 2, "-10", "-9.9"),
        ("1/3", 0, "0.3", "0.4"),
        ("0", 0, "0", "0"),
    ] {
        assert_eq!(r(x).to_decimal_bounds(p), (lo.into(), hi.into()));
    }
    for (x, s) in [
        ("1/8", "0.125"),
        ("-7/20", "-0.35"),
        ("5/2", "2.5"),
        ("0", "0"),
    ] {
        let decimal = r(x).to_decimal_exact().unwrap();
        assert!(!decimal.contains('/'));
        assert_eq!(r(&decimal), r(s));
    }
    for x in ["1/3", "1/7", "11/30"] {
        assert!(r(x).to_decimal_exact().is_none());
    }
    let mut rng = Rng(0xdec1a1);
    for _ in 0..2048 {
        let value = rng.rational();
        let p = (rng.next() % 8 + 1) as u32;
        let (lo, hi) = eq(&value).to_decimal_bounds(p);
        let (lo, hi) = (oq(&r(&lo)), oq(&r(&hi)));
        assert!(lo <= value && value <= hi);
        if value == 0 {
            assert_eq!(lo, hi);
            continue;
        }
        // Determine decimal exponent independently using exact powers of ten.
        let magnitude = value.clone().abs();
        let mut e = 0i32;
        let mut power = q(1, 1);
        while magnitude < power {
            e -= 1;
            power /= 10;
        }
        while magnitude >= power.clone() * 10 {
            e += 1;
            power *= 10;
        }
        let shift = e - p as i32 + 1;
        let ten = Z::from(10).pow(shift.unsigned_abs());
        let step = if shift >= 0 {
            Q::from(ten)
        } else {
            Q::from((Z::from(1), ten))
        };
        let scaled = value.clone() / &step;
        let floor = scaled.clone().floor();
        let ceil = scaled.ceil();
        assert_eq!(lo, floor * &step);
        assert_eq!(hi, ceil * step);
    }
});
