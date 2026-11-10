#![allow(missing_docs)]
use catallax_numbers::{Ball, Decider, OperationAllowance, with_operation_budget};
use std::time::Instant;

fn same(a: &Ball, b: &Ball) {
    assert_eq!(a.mid(), b.mid());
    assert_eq!(a.rad(), b.rad());
}

#[test]
fn pi_reuse_is_exact_precision_and_request_local() {
    let dec = Decider::default();
    for precision in [1, 64, 128, 256] {
        let expected = Ball::pi(precision);
        let (first, cold) =
            with_operation_budget(100_000, || Ball::pi_checked(precision, &dec).unwrap());
        same(&first, &expected);
        let (_, repeated) = with_operation_budget(100_000, || {
            same(&Ball::pi_checked(precision, &dec).unwrap(), &expected);
            same(&Ball::pi_checked(precision, &dec).unwrap(), &expected);
        });
        assert_eq!(repeated.used, cold.used + 1);
        let (_, next) = with_operation_budget(100_000, || {
            Ball::pi_checked(precision, &dec).unwrap();
        });
        assert_eq!(next.used, cold.used);
    }
    let (_, work) = with_operation_budget(100_000, || {
        for precision in [128, 64, 256, 128, 64] {
            same(
                &Ball::pi_checked(precision, &dec).unwrap(),
                &Ball::pi(precision),
            );
        }
    });
    assert!(!work.exhausted);
}

#[test]
fn pi_cache_hit_honors_budget_and_deadline() {
    let mut dec = Decider::default();
    let (_, cold) = with_operation_budget(100_000, || {
        Ball::pi_checked(128, &dec).unwrap();
    });
    let (_, work) = with_operation_budget(cold.used, || {
        Ball::pi_checked(128, &dec).unwrap();
        assert!(Ball::pi_checked(128, &dec).is_err());
    });
    assert!(work.exhausted);
    with_operation_budget(100_000, || {
        Ball::pi_checked(128, &dec).unwrap();
        dec.set_deadline(Some(Instant::now()));
        assert_eq!(
            Ball::pi_checked(128, &dec).unwrap_err().code,
            "numbers.deadline"
        );
    });
}

#[test]
fn deferred_allowance_keeps_pi_proofs_isolated() {
    let dec = Decider::default();
    let (_, cold) = with_operation_budget(100_000, || {
        Ball::pi_checked(128, &dec).unwrap();
    });
    with_operation_budget(100_000, || {
        Ball::pi_checked(128, &dec).unwrap();
        let deferred = OperationAllowance::reserve(cold.used - 1).unwrap();
        let (answer, work) = deferred.run(|| Ball::pi_checked(128, &dec));
        assert!(answer.is_err());
        assert!(work.exhausted);
        same(&Ball::pi_checked(128, &dec).unwrap(), &Ball::pi(128));
    });
}
