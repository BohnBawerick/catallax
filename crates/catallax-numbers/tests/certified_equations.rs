//! Certified root coordinates retain exact defining-equation authority.

use catallax_numbers::{Decider, Outcome, Rational, Real, RealPoly, StrictSign, ZeroProof};

fn q(value: &str) -> Rational {
    value.parse().unwrap()
}

fn cubic_root(coefficient: Real) -> Real {
    let equation = RealPoly::var(1, 0).pow_u32(3) - RealPoly::constant(1, coefficient);
    let mut coordinates = RealPoly::simple_root_in(
        &[equation],
        &[[q("1"), q("3/2")]],
        128,
        &mut Decider::default(),
    )
    .unwrap()
    .unwrap();
    coordinates.remove(0)
}

fn exact_zero(dec: &mut Decider, value: &Real) {
    match dec.sign(value) {
        Outcome::ProvenCoincident(proof) => assert_ne!(proof.proof, ZeroProof::SeparationBound),
        other => panic!("expected an exact defining-equation zero, got {other:?}"),
    }
}

#[test]
fn rational_cubic_retains_its_defining_equation() {
    let root = cubic_root(Real::from_i64(2));
    let residual = root.pow_u32(3) - Real::from_i64(2);
    let mut dec = Decider::default();
    exact_zero(&mut dec, &residual);
    exact_zero(&mut dec, &residual.square());
    exact_zero(&mut dec, &(&residual / &(&root + &root)));
    assert_eq!(dec.stats().coincident_separation, 0);
}

#[test]
fn quadratic_coefficient_definitions_survive_lowering() {
    let root = cubic_root(Real::from_i64(2).sqrt());
    let residual = root.pow_u32(3) - Real::from_i64(2).sqrt();
    let mut dec = Decider::default();
    exact_zero(&mut dec, &residual);
    exact_zero(&mut dec, &residual.square());
    assert_eq!(dec.stats().coincident_separation, 0);
    let shifted = &residual + &Real::from_rational(q("1/1000000000000"));
    assert!(matches!(
        dec.sign(&shifted),
        Outcome::Proven(StrictSign::Positive)
    ));
}

#[test]
fn four_unknowns_share_original_equation_authority() {
    let [x, y, z, w] = std::array::from_fn(|axis| RealPoly::var(4, axis));
    let equations = [
        x.pow_u32(3) - RealPoly::constant(4, Real::from_i64(2).sqrt()),
        y.clone() - x.pow_u32(2),
        z.clone() + x.clone() - y.clone(),
        w - x - y - z,
    ];
    let bounds = [
        [q("28/25"), q("113/100")],
        [q("5/4"), q("32/25")],
        [q("3/25"), q("4/25")],
        [q("5/2"), q("64/25")],
    ];
    let mut dec = Decider::default();
    let coordinates = RealPoly::simple_root_in(&equations, &bounds, 128, &mut dec)
        .unwrap()
        .unwrap();
    for equation in equations {
        exact_zero(&mut dec, &equation.eval_real(&coordinates));
    }
    exact_zero(
        &mut dec,
        &(&coordinates[3] - &(&Real::from_i64(2) * &coordinates[1])),
    );
    assert_eq!(dec.stats().coincident_separation, 0);
}

#[test]
fn nearby_certificates_do_not_exchange_equations() {
    let epsilon = Real::from_rational(q("1/1000000000000"));
    let left = cubic_root(Real::from_i64(2));
    let right = cubic_root(&Real::from_i64(2) + &epsilon);
    let difference = left.pow_u32(3) - right.pow_u32(3);
    assert!(matches!(
        Decider::default().sign(&difference),
        Outcome::Proven(StrictSign::Negative)
    ));
}

#[test]
fn defining_zero_cannot_hide_an_undefined_divisor() {
    let root = cubic_root(Real::from_i64(2));
    let residual = root.pow_u32(3) - Real::from_i64(2);
    let undefined = &residual / &Real::zero();
    assert!(matches!(
        Decider::default().sign(&undefined),
        Outcome::Refused(error) if error.code == "numbers.division-by-zero"
    ));
}

#[test]
fn a_conjugate_outside_the_original_coefficient_definition_is_not_a_root() {
    let equation = RealPoly::var(1, 0).pow_u32(3) - RealPoly::constant(1, Real::from_i64(2).sqrt());
    assert!(matches!(
        RealPoly::simple_root_in(
            &[equation],
            &[[q("-3/2"), q("-1")]],
            128,
            &mut Decider::default(),
        ),
        Err(error) if error.code == "numbers.root-not-certified"
    ));
}

fn joint_cubic_root(coefficient: Real) -> Real {
    let equation = RealPoly::var(1, 0).pow_u32(3) - RealPoly::constant(1, coefficient);
    RealPoly::simple_algebraic_root_in(
        &[equation],
        &[[q("1"), q("3/2")]],
        128,
        &mut Decider::default(),
    )
    .unwrap()
    .remove(0)
}

#[test]
fn joint_equation_preflight_preserves_undefined_dependencies() {
    let root = joint_cubic_root(Real::from_i64(2));
    let residual = root.pow_u32(3) - Real::from_i64(2);
    for (undefined, expected) in [
        (&residual / &Real::zero(), "numbers.division-by-zero"),
        (
            &residual * &Real::from_i64(-1).sqrt(),
            "numbers.sqrt-of-negative",
        ),
    ] {
        assert!(matches!(
            Decider::default().sign(&undefined),
            Outcome::Refused(error) if error.code == expected
        ));
    }
}

#[test]
fn joint_equation_preflight_keeps_nearby_certificates_distinct() {
    let epsilon = Real::from_rational(q("1/1000000000000"));
    let left = joint_cubic_root(Real::from_i64(2));
    let right = joint_cubic_root(&Real::from_i64(2) + &epsilon);
    assert!(matches!(
        Decider::default().sign(&(left.pow_u32(3) - right.pow_u32(3))),
        Outcome::Proven(StrictSign::Negative)
    ));
}

#[test]
fn joint_roots_preserve_original_coefficient_dags() {
    let scalar = (RealPoly::var(1, 0).pow_u32(3) - RealPoly::constant(1, Real::from_i64(2)))
        .real_roots(&mut Decider::default())
        .unwrap()
        .remove(0);
    let square_root = Real::from_i64(2).sqrt();
    let arithmetic = ((&square_root + &Real::from_i64(3)).sqrt() + &square_root) / &square_root;
    for coefficient in [
        square_root,
        scalar,
        cubic_root(Real::from_i64(2)),
        arithmetic,
    ] {
        let root = joint_cubic_root(coefficient.clone());
        let residual = root.pow_u32(3) - coefficient;
        let mut dec = Decider::default();
        exact_zero(&mut dec, &residual);
        exact_zero(&mut dec, &residual.square());
        assert!(matches!(
            dec.sign(&(&residual + &Real::from_rational(q("1/1000000000000")))),
            Outcome::Proven(StrictSign::Positive)
        ));
        assert_eq!(dec.stats().coincident_separation, 0);
    }
}

#[test]
fn four_unknown_joint_roots_retain_each_original_equation() {
    let [x, y, z, w] = std::array::from_fn(|axis| RealPoly::var(4, axis));
    let equations = [
        x.pow_u32(3) - RealPoly::constant(4, Real::from_i64(2).sqrt()),
        y.clone() - x.pow_u32(2),
        z.clone() + x.clone() - y.clone(),
        w - x - y - z,
    ];
    let bounds = [
        [q("28/25"), q("113/100")],
        [q("5/4"), q("32/25")],
        [q("3/25"), q("4/25")],
        [q("5/2"), q("64/25")],
    ];
    let mut dec = Decider::default();
    let coordinates =
        RealPoly::simple_algebraic_root_in(&equations, &bounds, 128, &mut dec).unwrap();
    for equation in equations {
        exact_zero(&mut dec, &equation.eval_real(&coordinates));
    }
    exact_zero(
        &mut dec,
        &(&coordinates[3] - &(&Real::from_i64(2) * &coordinates[1])),
    );
    assert_eq!(dec.stats().coincident_separation, 0);
}

#[test]
fn joint_root_refuses_a_wrong_coefficient_conjugate() {
    let equation = RealPoly::var(1, 0).pow_u32(3) - RealPoly::constant(1, Real::from_i64(2).sqrt());
    assert!(matches!(
        RealPoly::simple_algebraic_root_in(
            &[equation],
            &[[q("-3/2"), q("-1")]],
            128,
            &mut Decider::default(),
        ),
        Err(error) if error.code == "numbers.root-not-certified"
    ));
}
