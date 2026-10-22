//! Algebraic coefficient identities, root selectors, angle sectors and budgets.

use catallax_numbers::{
    Ball, BudgetPolicy, Decider, MPoly, Outcome, Rational, Real, RealPoly, StrictSign,
    sin_cos_degrees,
};

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn zero(value: &Real, dec: &mut Decider) {
    assert!(
        matches!(dec.sign(value), Outcome::ProvenCoincident(_)),
        "{value:?}"
    );
}

#[test]
fn polynomial_degree_requires_proof_of_coefficient_cancellation() {
    let root = Real::from_i64(2).sqrt();
    let cancellation = &root.square() - &Real::from_i64(2);
    let p = RealPoly::from_terms(1, vec![(vec![3], cancellation), (vec![1], root.clone())]);
    assert!(p.as_rational().is_none());
    let mut dec = Decider::default();
    assert_eq!(p.total_degree(&mut dec).unwrap(), Some(1));
    zero(
        &(&p.eval_real(std::slice::from_ref(&root)) - &Real::from_i64(2)),
        &mut dec,
    );
    zero(
        &(&p.partial(0).eval_real(&[Real::from_i64(7)]) - &root),
        &mut dec,
    );
    let substituted = p.substitute(0, &RealPoly::constant(1, root));
    zero(
        &(&substituted.eval_real(&[Real::from_i64(9)]) - &Real::from_i64(2)),
        &mut dec,
    );
    assert_eq!(substituted.total_degree(&mut dec).unwrap(), Some(0));
}

#[test]
fn rational_polynomial_arithmetic_preserves_exact_coefficients() {
    let x = MPoly::var(2, 0);
    let y = MPoly::var(2, 1);
    let a = &x + &y;
    let b = &x - &y;
    let real_a = RealPoly::from(a.clone());
    let real_b = RealPoly::from(b.clone());
    assert_eq!((&real_a * &real_b).as_rational(), Some(&a * &b));
    assert_eq!(real_a.partial(1).as_rational(), Some(a.partial(1)));
    assert_eq!(
        real_a.substitute(0, &real_b).as_rational(),
        Some(a.substitute(0, &b))
    );
}

#[test]
fn polynomial_enclosure_includes_irrational_coefficient_error() {
    let root = Real::from_i64(2).sqrt();
    let p = RealPoly::from_terms(1, vec![(vec![2], root.clone()), (vec![0], -root)]);
    let dec = Decider::default();
    let values = [Ball::from_rational(&q("1"), 128).hull(&Ball::from_rational(&q("2"), 128))];
    let bound = p.eval_ball(&values, &dec).unwrap();
    assert!(bound.lo() <= Rational::zero());
    assert!(bound.hi() >= q("4"));
    assert!(bound.hi() < q("5"));
}

#[test]
fn degree_angles_match_closed_forms_and_exact_sectors() {
    for (angle, sine, cosine) in [
        ("0", Real::zero(), Real::one()),
        (
            "30",
            Real::from(q("1/2")),
            &Real::from_i64(3).sqrt() * &Real::from(q("1/2")),
        ),
        (
            "45",
            &Real::from_i64(2).sqrt() * &Real::from(q("1/2")),
            &Real::from_i64(2).sqrt() * &Real::from(q("1/2")),
        ),
        ("90", Real::one(), Real::zero()),
        ("180", Real::zero(), -Real::one()),
        ("-270", Real::one(), Real::zero()),
        (
            "405",
            &Real::from_i64(2).sqrt() * &Real::from(q("1/2")),
            &Real::from_i64(2).sqrt() * &Real::from(q("1/2")),
        ),
    ] {
        let mut dec = Decider::default();
        let [s, c] = sin_cos_degrees(&q(angle), &mut dec).unwrap();
        zero(&(&s - &sine), &mut dec);
        zero(&(&c - &cosine), &mut dec);
    }
}

#[test]
fn nonquadratic_angle_keeps_selected_algebraic_root() {
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 100_000_000,
        max_precision_bits: 4096,
    });
    let [s, c] = sin_cos_degrees(&q("20"), &mut dec).unwrap();
    assert!(c.as_rational().is_none());
    // cos(60 degrees) = 4*cos(20 degrees)^3 - 3*cos(20 degrees).
    let triple = &(&Real::from_i64(4) * &c.pow_u32(3)) - &(&Real::from_i64(3) * &c);
    zero(&(&triple - &Real::from(q("1/2"))), &mut dec);
    zero(&(&(&s.square() + &c.square()) - &Real::one()), &mut dec);
    assert!(matches!(
        dec.sign(&s),
        Outcome::Proven(StrictSign::Positive)
    ));
    assert!(matches!(
        dec.sign(&(&c - &Real::from(q("9/10")))),
        Outcome::Proven(StrictSign::Positive)
    ));
}

#[test]
fn angle_certificate_work_is_thread_and_cache_independent() {
    fn run() -> (Rational, Rational, catallax_numbers::DecisionStats) {
        let mut dec = Decider::default();
        let [s, c] = sin_cos_degrees(&q("72"), &mut dec).unwrap();
        let a = dec.enclose_at(&s, 96).unwrap();
        let b = dec.enclose_at(&c, 96).unwrap();
        (a.lo(), b.hi(), dec.stats().clone())
    }
    let expected = run();
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..4).map(|_| scope.spawn(run)).collect();
        for handle in handles {
            assert_eq!(handle.join().unwrap(), expected);
        }
    });
    let mut dec = Decider::default();
    sin_cos_degrees(&q("72"), &mut dec).unwrap();
    let cold = dec.take_stats();
    sin_cos_degrees(&q("72"), &mut dec).unwrap();
    assert_eq!(dec.take_stats(), cold);
}

#[test]
fn angle_construction_stops_at_the_existing_work_limit() {
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 0,
        max_precision_bits: 64,
    });
    assert_eq!(
        sin_cos_degrees(&q("20"), &mut dec).unwrap_err().code,
        "numbers.budget-exhausted"
    );
}

#[test]
fn algebraic_polynomial_roots_exclude_the_wrong_square_root_conjugate() {
    let root = Real::from_i64(2).sqrt();
    let p = RealPoly::from_terms(1, vec![(vec![1], Real::one()), (vec![0], -root.clone())]);
    let mut dec = Decider::default();
    let roots = p.real_roots(&mut dec).unwrap();
    assert_eq!(roots.len(), 1);
    zero(&(&roots[0] - &root), &mut dec);
    assert!(matches!(
        dec.sign(&roots[0]),
        Outcome::Proven(StrictSign::Positive)
    ));
}

#[test]
fn algebraic_polynomial_roots_keep_shared_coefficient_correlations() {
    let root = Real::from_i64(2).sqrt();
    let p = RealPoly::from_terms(
        1,
        vec![
            (vec![2], Real::one()),
            (vec![1], -(&Real::from_i64(2) * &root)),
            (vec![0], root.square()),
        ],
    );
    let mut dec = Decider::default();
    let roots = p.real_roots(&mut dec).unwrap();
    assert_eq!(roots.len(), 1);
    zero(&(&roots[0] - &root), &mut dec);
}

#[test]
fn rational_projection_reduces_dependent_radicals_without_guessing() {
    let a = Real::from_i64(2).sqrt();
    let b = Real::from_i64(3).sqrt();
    let product = &(&a + &b) * &(&a - &b);
    let mut dec = Decider::default();
    assert_eq!(dec.rational_value(&product).unwrap(), Some(q("-1")));
    assert_eq!(
        dec.rational_value(&(&Real::one() / &product)).unwrap(),
        Some(q("-1"))
    );
    assert!(dec.rational_value(&a).unwrap().is_none());
}

#[test]
fn rational_projection_never_erases_an_undefined_operand() {
    for (bad, reason) in [
        (&Real::one() / &Real::zero(), "numbers.division-by-zero"),
        (Real::from_i64(-1).sqrt(), "numbers.sqrt-of-negative"),
    ] {
        for erased in [&bad - &bad, bad.pow_u32(0), &Real::zero() * &bad] {
            assert_eq!(
                Decider::default().rational_value(&erased).unwrap_err().code,
                reason
            );
        }
    }
}

#[test]
fn algebraic_resultants_preserve_signs_and_ignore_proven_zero_leading_terms() {
    let root = Real::from_i64(2).sqrt();
    let x = RealPoly::var(2, 0);
    let y = RealPoly::var(2, 1);
    let constant = |n| RealPoly::constant(2, Real::from_i64(n));
    let a = &(&x - &RealPoly::constant(2, root.clone())) * &(&x + &constant(1));
    let mut dec = Decider::default();
    let expected = &(&y - &RealPoly::constant(2, root.clone())) * &(&y + &constant(1));
    assert!(
        a.resultant(&(&x - &y), 0, &mut dec)
            .unwrap()
            .equals(&expected, &mut dec)
            .unwrap()
    );
    let cancellation =
        RealPoly::from_terms(2, vec![(vec![3, 0], &root.square() - &Real::from_i64(2))]);
    let b = &(&x.pow_u32(2) - &y) + &cancellation;
    let expected = &(&constant(2) - &y) * &(&constant(1) - &y);
    assert!(
        a.resultant(&b, 0, &mut dec)
            .unwrap()
            .equals(&expected, &mut dec)
            .unwrap()
    );
}

#[test]
fn constructible_half_angles_obey_independent_double_angle_identities() {
    let half = Real::from(q("1/2"));
    let two = Real::from_i64(2);
    let root_two_half = &two.sqrt() * &half;
    let root_three_half = &Real::from_i64(3).sqrt() * &half;
    for (angle, sine_double, cosine_double) in [
        ("15", half.clone(), root_three_half.clone()),
        ("75", half, -root_three_half),
        ("45/2", root_two_half.clone(), root_two_half.clone()),
        ("135/2", root_two_half.clone(), -root_two_half),
    ] {
        let mut dec = Decider::default();
        let [sine, cosine] = sin_cos_degrees(&q(angle), &mut dec).unwrap();
        zero(&(&(&two * &sine) * &cosine - sine_double), &mut dec);
        zero(
            &(&cosine.square() - &sine.square() - cosine_double),
            &mut dec,
        );
        zero(&(&sine.square() + &cosine.square() - Real::one()), &mut dec);
        assert!(matches!(
            dec.sign(&sine),
            Outcome::Proven(StrictSign::Positive)
        ));
        assert!(matches!(
            dec.sign(&cosine),
            Outcome::Proven(StrictSign::Positive)
        ));
    }
}
