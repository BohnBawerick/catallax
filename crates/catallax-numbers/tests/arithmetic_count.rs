//! `OperationWork::uncharged`: the number types count their own arithmetic.
//!
//! Checkpoint 2 decisions, section 4.5, predictor step 1, and the owner's
//! ruling of 21 November 2026: count every ball and rational operation at the
//! decider's rates and charge none of it. The rules are in
//! `src/operation.rs`, "The count the number types keep".

use catallax_numbers::{
    Ball, Decider, OperationAllowance, OperationTask, OperationWork, Rational, Real,
    uncounted_memo, with_operation_budget,
};

fn r(text: &str) -> Rational {
    Rational::parse(text).unwrap()
}

fn work(run: impl FnOnce()) -> OperationWork {
    with_operation_budget(1_000_000, run).1
}

// What a request that only runs `run` counts. Nothing may be charged.
fn counted(run: impl FnOnce()) -> u64 {
    let work = work(run);
    assert_eq!(work.used, 0, "counting charges nothing");
    assert!(!work.exhausted);
    work.uncharged
}

// Limbs are 64 bits. A rational's size is its numerator's limbs plus its
// denominator's. Every operation counts one unit more, the unit a charge of
// the decider pays for its deadline check.
#[test]
fn a_rational_operation_counts_the_deciders_limb_price() {
    let third = r("1/3");
    let wide = Rational::from_i64(1).mul_pow2(64);
    assert_eq!(counted(|| drop(&third * &wide)), 2 * 3 + 1);
    assert_eq!(counted(|| drop(&wide / &third)), 3 * 2 + 1);
    assert_eq!(counted(|| drop(&third + &wide)), 2 + 3 + 1);
    assert_eq!(counted(|| drop(&wide - &wide)), 3 + 3 + 1);
    assert_eq!(counted(|| drop(third.checked_div(&wide))), 2 * 3 + 1);
    // An operand of zero or one is priced as any other: the decider prices
    // the operation before it knows what it will find.
    assert_eq!(counted(|| drop(&third * &Rational::zero())), 2 * 2 + 1);
    assert_eq!(counted(|| drop(&wide + &Rational::zero())), 3 + 2 + 1);
    // A power is the products of a repeated multiplication, one fewer than the
    // exponent, each priced as the power by the base.
    assert_eq!(counted(|| drop(third.pow(3))), 2 * (2 * 2) + 1);
    assert_eq!(counted(|| drop(third.pow(1))), 1);
    // A decimal conversion is the value times itself, as one operation: the
    // arithmetic it is made of is not counted again.
    assert_eq!(counted(|| drop(third.to_decimal_bounds(5))), 2 * 2 + 1);
    assert_eq!(counted(|| drop(wide.to_decimal_bounds(5))), 3 * 3 + 1);
    assert_eq!(counted(|| drop(wide.to_decimal_exact())), 3 * 3 + 1);
}

#[test]
fn what_the_decider_has_no_price_for_is_not_counted() {
    let third = r("1/3");
    let wide = Rational::from_i64(1).mul_pow2(200);
    assert_eq!(
        counted(|| {
            assert!(third < wide);
            drop(-&third);
            drop(third.abs());
            drop(third.recip());
            drop(wide.mul_pow2(-7));
            drop(wide.floor());
            let _ = wide.to_f64_nearest();
            drop(Ball::from_rational(&third, 128).mid());
        }),
        0
    );
}

// L = ceil(p / 64). A sum or a difference counts 4L and a product, a square,
// a reciprocal or a square root 8L^2, at the larger of the two working
// precisions, each as one operation.
#[test]
fn a_ball_operation_counts_the_deciders_ball_price() {
    for (precision, limbs) in [(1u32, 1u64), (64, 1), (65, 2), (200, 4)] {
        let a = Ball::from_rational(&r("7/3"), precision);
        let b = Ball::from_rational(&r("5/11"), precision);
        let sum = 4 * limbs + 1;
        let product = 8 * limbs * limbs + 1;
        assert_eq!(counted(|| drop(&a + &b)), sum, "{precision}");
        assert_eq!(counted(|| drop(&a - &b)), sum, "{precision}");
        assert_eq!(counted(|| drop(&a * &b)), product, "{precision}");
        assert_eq!(counted(|| drop(a.square())), product, "{precision}");
        assert_eq!(counted(|| drop(a.sqrt())), product, "{precision}");
        assert_eq!(counted(|| drop(a.recip())), product, "{precision}");
        // A quotient is its reciprocal and its product.
        assert_eq!(
            counted(|| drop(a.checked_div(&b))),
            2 * product,
            "{precision}"
        );
    }
    // The larger working precision prices a mixed operation.
    let narrow = Ball::from_rational(&r("7/3"), 64);
    let wide = Ball::from_rational(&r("5/11"), 200);
    assert_eq!(counted(|| drop(&narrow + &wide)), 4 * 4 + 1);
    assert_eq!(counted(|| drop(&narrow * &wide)), 8 * 4 * 4 + 1);
    // A function made of ball operations counts the operations it runs.
    let turn = Ball::from_rational(&r("1/5"), 128);
    assert!(counted(|| drop(turn.atan())) > 8 * 2 * 2);
}

#[test]
fn series_reciprocal_reuse_preserves_enclosure_and_work() {
    let input = Ball::from_rational(&r("1/5"), 128);
    let run = || {
        let (value, work) = with_operation_budget(1_000_000, || input.sin());
        let OperationWork {
            used,
            exhausted,
            refunded,
            uncharged,
        } = work;
        (
            value.mid(),
            value.rad(),
            (used, exhausted, refunded, uncharged),
        )
    };
    let uncached = [run(), run()];
    let cached = Ball::with_series_reciprocals(|| [run(), run()]);
    assert_eq!(cached, uncached);
}

// A decision charges its arithmetic to `used` as it goes. The number types
// do not count it a second time.
#[test]
fn a_decision_charges_its_arithmetic_and_counts_none_of_it() {
    let root = Real::from_rational(r("2")).sqrt();
    let value = &(&root * &Real::from_rational(r("1/3"))) + &root;
    let work = work(|| {
        let mut dec = Decider::default();
        assert!(dec.enclose_at(&value, 256).is_ok());
        assert!(matches!(
            dec.sign(&(&value - &root)),
            catallax_numbers::Outcome::Proven(_)
        ));
        assert!(dec.enclose(&value, &r("1/1000000000000")).is_ok());
    });
    assert!(work.used > 100, "the decisions charged {}", work.used);
    assert_eq!(work.uncharged, 0);
}

// `reserve_exact_arithmetic` charges `used` for a rational operation its
// caller runs next. That operation is counted once, by the charge.
#[test]
fn a_reserved_operation_is_charged_and_not_counted() {
    let a = Rational::from_i64(3).mul_pow2(100);
    let b = r("5/7");
    let price = 3 * 2 + 1;
    let reserved = work(|| {
        let mut dec = Decider::default();
        dec.reserve_exact_arithmetic(&a, &b, true).unwrap();
        drop(&a * &b);
    });
    assert_eq!((reserved.used, reserved.uncharged), (price, 0));
    let plain = work(|| drop(&a * &b));
    assert_eq!((plain.used, plain.uncharged), (0, price));
    // A second operation the reservation did not pay for is counted.
    let both = work(|| {
        let mut dec = Decider::default();
        dec.reserve_exact_arithmetic(&a, &b, true).unwrap();
        drop(&a * &b);
        drop(&a + &b);
    });
    assert_eq!((both.used, both.uncharged), (price, 3 + 2 + 1));
    // A reservation whose operation never runs leaves nothing to report, and
    // the report never goes below zero.
    let unused = work(|| {
        let mut dec = Decider::default();
        dec.reserve_exact_arithmetic(&a, &b, true).unwrap();
    });
    assert_eq!((unused.used, unused.uncharged), (price, 0));
    // A refused reservation paid for nothing.
    let (_, refused) = with_operation_budget(3, || {
        let mut dec = Decider::default();
        assert!(dec.reserve_exact_arithmetic(&a, &b, true).is_err());
        drop(&a * &b);
    });
    assert!(refused.exhausted);
    assert_eq!(refused.uncharged, price);
}

// Arithmetic that does not count itself is reserved with
// `reserve_uncounted_arithmetic`: the same charge, and nothing off the count.
// A product of two `Real`s runs when its rational form is first read, inside
// the value's memo, and a comparison cross-multiplies inside itself.
#[test]
fn a_reservation_for_arithmetic_that_does_not_count_itself_takes_nothing_off() {
    let a = Rational::from_i64(3).mul_pow2(100);
    let b = r("5/7");
    let price = 3 * 2 + 1;
    let expected = &a * &b;
    let lazy = work(|| {
        let mut dec = Decider::default();
        dec.reserve_uncounted_arithmetic(&a, &b, true).unwrap();
        let product = &Real::from_rational(a.clone()) * &Real::from_rational(b.clone());
        assert_eq!(product.as_rational(), Some(expected.clone()));
        dec.reserve_uncounted_arithmetic(&a, &b, true).unwrap();
        assert!(a > b);
        // A counted operation after them is counted in full.
        drop(&a + &b);
    });
    assert_eq!((lazy.used, lazy.uncharged), (2 * price, 3 + 2 + 1));
    // Reserved as if it counted itself, the same work leaves the count
    // short by what was charged: the sum after it is not reported.
    let short = work(|| {
        let mut dec = Decider::default();
        dec.reserve_exact_arithmetic(&a, &b, true).unwrap();
        let product = &Real::from_rational(a.clone()) * &Real::from_rational(b.clone());
        assert_eq!(product.as_rational(), Some(expected.clone()));
        drop(&a + &b);
    });
    assert_eq!((short.used, short.uncharged), (price, 0));
}

// A memo on a shared value is filled once, by whichever thread or request
// reads it first. Its fill is not counted, so a request counts the same
// whether it finds the memo filled or fills it.
#[test]
fn a_shared_memo_fill_is_not_counted() {
    let a = Rational::from_i64(3).mul_pow2(100);
    let b = r("5/7");
    assert_eq!(counted(|| drop(uncounted_memo(|| &a * &b))), 0);
    // The count is open again after the fill.
    assert_eq!(
        counted(|| {
            drop(uncounted_memo(|| &a * &b));
            drop(&a + &b);
        }),
        3 + 2 + 1
    );
    // A value's own rational form is such a memo.
    let shared = &(&Real::from_rational(a.clone()) * &Real::from_rational(b.clone()))
        + &Real::from_rational(r("1/3"));
    let first = counted(|| assert!(shared.as_rational().is_some()));
    let second = counted(|| assert!(shared.as_rational().is_some()));
    assert_eq!(first, second, "the first reader fills the memo");
    let fresh = &(&Real::from_rational(a.clone()) * &Real::from_rational(b.clone()))
        + &Real::from_rational(r("1/3"));
    std::thread::scope(|scope| {
        scope
            .spawn(|| assert!(fresh.as_rational().is_some()))
            .join()
            .unwrap();
    });
    let filled_elsewhere = counted(|| assert!(fresh.as_rational().is_some()));
    assert_eq!(first, filled_elsewhere, "another thread filled the memo");
}

// A reserved allowance run inside a request has a counter of its own. What
// it counts is its own report, and the request's count goes on where it was.
#[test]
fn a_nested_allowance_counts_for_itself() {
    let a = Rational::from_i64(3).mul_pow2(100);
    let b = r("5/7");
    let (inner, outer) = with_operation_budget(1_000, || {
        drop(&a * &b);
        let nested = OperationAllowance::reserve(100).unwrap();
        let ((), inner) = nested.run(|| {
            drop(&a + &b);
            drop(&a + &b);
        });
        drop(&a * &b);
        inner
    });
    assert_eq!(inner.uncharged, 2 * (3 + 2 + 1));
    assert_eq!(outer.uncharged, 2 * (3 * 2 + 1));
}

// An item of a numeric batch reports what it counted to its request, admitted
// or discarded, as it reports its refused attempts.
#[test]
fn a_batch_item_reports_its_count() {
    let a = Rational::from_i64(3).mul_pow2(100);
    let b = r("5/7");
    let (_, work) = with_operation_budget(1_000, || {
        drop(&a + &b);
        let admitted = OperationTask::capture().run(|| drop(&a * &b));
        let discarded = OperationTask::capture().run(|| drop(&a * &b));
        admitted.commit().unwrap();
        drop(discarded);
        drop(&a + &b);
    });
    assert_eq!(work.uncharged, 2 * (3 + 2 + 1) + 2 * (3 * 2 + 1));
}

// A request's count does not depend on what its thread ran before it.
#[test]
fn a_request_counts_the_same_after_another_request() {
    let run = || {
        let turn = Ball::from_rational(&r("3/5"), 192);
        let dec = Decider::default();
        drop(turn.atan_checked(&dec));
        drop(turn.sin());
        drop(&r("22/7") * &r("355/113"));
    };
    let first = work(run);
    let second = work(run);
    assert!(first.uncharged > 0);
    assert_eq!(
        (first.used, first.uncharged),
        (second.used, second.uncharged)
    );
    // Arithmetic outside any request is nobody's.
    run();
    assert_eq!(work(run).uncharged, first.uncharged);
}
