#![allow(missing_docs)]

use catallax_numbers::{Decider, Outcome, Rational, Real, RealPoly, UPoly, with_operation_budget};

fn constant(value: Real) -> RealPoly {
    RealPoly::constant(1, value)
}

#[test]
fn algebraic_common_parameter_survives_independent_zero_leading_coefficients() {
    let mut dec = Decider::default();
    let x = RealPoly::var(1, 0);
    let alpha = Real::from_i64(2).sqrt();
    let beta = Real::from_i64(3).sqrt();
    let Outcome::Proven(mut roots) =
        dec.real_roots(&UPoly::new([-2, 0, 1].map(Rational::from_i64).to_vec()))
    else {
        panic!("root certificate")
    };
    let same_alpha = Real::from_algebraic(roots.pop().unwrap().value);
    let zero = &alpha - &same_alpha;
    let factor = &x - &constant(alpha.clone());
    let f = &factor * &(&x.pow_u32(2) + &constant(beta.clone()))
        + RealPoly::from_terms(1, vec![(vec![7], zero)]);
    let g = &factor * &(&x + &constant(beta));
    let divisor = f.gcd(&g, &mut dec).unwrap();
    assert!(divisor.equals(&factor, &mut dec).unwrap());
    let root = -&divisor.eval_real(&[Real::zero()]);
    assert!(dec.compare(&root, &alpha).is_coincident());
}

#[test]
fn coprime_coefficients_and_zero_inputs_keep_the_field_gcd_semantics() {
    let mut dec = Decider::default();
    let x = RealPoly::var(1, 0);
    let a = &x.pow_u32(2) + &constant(Real::from_i64(2).sqrt());
    let b = &x - &constant(Real::from_i64(3).sqrt());
    assert!(
        a.gcd(&b, &mut dec)
            .unwrap()
            .equals(&constant(Real::one()), &mut dec)
            .unwrap()
    );
    let zero = RealPoly::zero(1);
    assert!(
        zero.gcd(&zero, &mut dec)
            .unwrap()
            .is_zero(&mut dec)
            .unwrap()
    );
    let scaled = b.scale(&Real::from_i64(7));
    assert!(
        zero.gcd(&scaled, &mut dec)
            .unwrap()
            .equals(&b, &mut dec)
            .unwrap()
    );
}

#[test]
fn cancellation_and_invalid_coefficients_cannot_short_circuit_to_a_unit() {
    let mut dec = Decider::default();
    let unit = constant(Real::one());
    let invalid = RealPoly::from_terms(
        1,
        vec![
            (vec![0], &Real::one() / &Real::zero()),
            (vec![2], Real::from_i64(2).sqrt()),
        ],
    );
    assert_eq!(
        unit.gcd(&invalid, &mut dec).unwrap_err().code,
        "numbers.division-by-zero"
    );
    let (result, work) = with_operation_budget(0, || unit.gcd(&RealPoly::var(1, 0), &mut dec));
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    assert!(work.exhausted);
}
