//! Exact coordinate reconstruction, closed bounds and conjugate selection.

use catallax_numbers::{
    BoxRoots, Decider, IsolatedRoot, MPoly, Outcome, Rational, Real, RealPoly, RootBox, ZeroProof,
    with_operation_budget,
};
use std::sync::Arc;

fn certified(
    equations: Vec<MPoly>,
    lower: Vec<Rational>,
    upper: Vec<Rational>,
    dec: &mut Decider,
) -> IsolatedRoot {
    let system = Arc::new(catallax_numbers::PolySystem::new(equations).unwrap());
    let bounds = RootBox::new(lower, upper).unwrap();
    match dec.certify_root_at(&system, &bounds, 64) {
        Outcome::Proven(BoxRoots::One(root)) => root,
        other => panic!("the fixture must certify one nonsingular root: {other:?}"),
    }
}

#[test]
fn a_rational_coordinate_on_the_closed_box_boundary_is_reconstructed() {
    let (_, work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let x = MPoly::var(1, 0);
        let three = Rational::from_i64(3);
        let root = certified(
            vec![&x - &MPoly::constant(1, three.clone())],
            vec![three.clone()],
            vec![three.clone()],
            &mut dec,
        );
        let value = dec.reconstruct_coordinate(&root, 0).unwrap().unwrap();
        assert_eq!(dec.rational_value(&value).unwrap(), Some(three));
    });
    assert!(!work.exhausted);
}

#[test]
fn coupled_coordinates_use_the_exact_positive_conjugate() {
    let (_, work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let x = MPoly::var(2, 0);
        let y = MPoly::var(2, 1);
        let root = certified(
            vec![
                &x - &y,
                &(&y * &y) - &MPoly::constant(2, Rational::from_i64(2)),
            ],
            vec![Rational::one(); 2],
            vec![Rational::from_i64(2); 2],
            &mut dec,
        );
        let expected = Real::from_i64(2).sqrt();
        for axis in 0..2 {
            let value = dec.reconstruct_coordinate(&root, axis).unwrap().unwrap();
            match dec.sign(&(&value - &expected)) {
                Outcome::ProvenCoincident(proof) => {
                    assert_ne!(proof.proof, ZeroProof::SeparationBound)
                }
                other => {
                    panic!("reconstruction must equal the positive square root exactly: {other:?}")
                }
            }
        }
    });
    assert!(!work.exhausted);
}

#[test]
fn two_projected_values_inside_the_coordinate_box_do_not_select_a_conjugate() {
    let (_, work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let x = MPoly::var(2, 0);
        let y = MPoly::var(2, 1);
        let delta = Rational::one().mul_pow2(-200);
        let one = MPoly::constant(2, Rational::one());
        let y_squared = &y * &y;
        let root = certified(
            vec![
                &(&(&x - &one) - &y_squared) - &y.scale(&delta),
                &y_squared - &one,
            ],
            vec![
                Rational::from_ratio(9, 5).unwrap(),
                Rational::from_ratio(9, 10).unwrap(),
            ],
            vec![
                Rational::from_ratio(11, 5).unwrap(),
                Rational::from_ratio(11, 10).unwrap(),
            ],
            &mut dec,
        );
        let other = &Rational::from_i64(2) - &delta;
        let selected = &Rational::from_i64(2) + &delta;
        assert!(root.root_box().lo()[0] <= other && selected <= root.root_box().hi()[0]);
        assert!(dec.reconstruct_coordinate(&root, 0).unwrap().is_none());
    });
    assert!(!work.exhausted);
}

#[test]
fn reconstructed_algebraic_system_coordinates_keep_the_selected_coefficient_and_conjugate() {
    let (_, work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let coefficient = Real::from_i64(2).sqrt();
        let equations = [
            RealPoly::var(2, 0) - RealPoly::constant(2, coefficient.clone()),
            RealPoly::var(2, 1).pow_u32(2) - RealPoly::constant(2, Real::from_i64(3)),
        ];
        let bounds = [
            [Rational::from_i64(1), Rational::from_i64(2)],
            [Rational::from_i64(-2), Rational::from_i64(-1)],
        ];
        let [y, x]: [Real; 2] = RealPoly::simple_algebraic_root_in_reconstructed(
            &equations,
            &bounds,
            &[1, 0],
            128,
            &mut dec,
        )
        .unwrap()
        .try_into()
        .unwrap();
        for (value, expected) in [(&x, coefficient), (&y, -Real::from_i64(3).sqrt())] {
            match dec.sign(&(value - &expected)) {
                Outcome::ProvenCoincident(proof) => {
                    assert_ne!(proof.proof, ZeroProof::SeparationBound);
                }
                other => panic!("the selected coordinate must be exactly equal: {other:?}"),
            }
        }
    });
    assert!(!work.exhausted);
}

#[test]
fn compacted_affine_coefficients_keep_the_principal_root_choices() {
    let (_, work) = with_operation_budget(20_000_000, || {
        let a = Real::from_i64(2).sqrt();
        let b = Real::from_i64(3).sqrt();
        let coefficient = &(&(&a + &b) * &Real::from_i64(3)) - &(&a - &b);
        let equation = (RealPoly::var(1, 0) - RealPoly::constant(1, coefficient.clone()))
            .scale(&Real::from_i64(-3));
        let mut dec = Decider::default();
        let value = RealPoly::simple_algebraic_root_in_reconstructed(
            std::slice::from_ref(&equation),
            &[[Rational::from_i64(9), Rational::from_i64(10)]],
            &[0],
            128,
            &mut dec,
        )
        .unwrap()
        .remove(0);
        assert!(matches!(
            dec.sign(&(&value - &coefficient)),
            Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound
        ));
        assert!(matches!(
            RealPoly::simple_algebraic_root_in_reconstructed(
                &[equation],
                &[[Rational::from_i64(-10), Rational::from_i64(-9)]],
                &[0],
                128,
                &mut dec,
            ),
            Err(error) if error.code == "numbers.root-not-certified"
        ));
    });
    assert!(!work.exhausted);
}

#[test]
fn compacted_affine_coefficients_do_not_discard_cancelled_undefined_dependencies() {
    let (_, work) = with_operation_budget(20_000_000, || {
        for undefined in [&Real::one() / &Real::zero(), Real::from_i64(-1).sqrt()] {
            let coefficient = &Real::from_i64(2).sqrt() + &(&Real::zero() * &undefined);
            let equation = RealPoly::var(1, 0) - RealPoly::constant(1, coefficient);
            let result = RealPoly::simple_algebraic_root_in_reconstructed(
                &[equation],
                &[[Rational::from_i64(1), Rational::from_i64(2)]],
                &[0],
                128,
                &mut Decider::default(),
            );
            assert!(
                result.is_err(),
                "an undefined coefficient must not certify a coordinate: {result:?}"
            );
        }
    });
    assert!(!work.exhausted);
}

#[test]
fn reconstructed_joint_roots_are_not_tied_to_coordinate_units() {
    let (_, work) = with_operation_budget(20_000_000, || {
        let precision = 128u32;
        for exponent in [0, 2 * precision] {
            let scale = Rational::one().mul_pow2(i64::from(exponent));
            let y_squared = RealPoly::var(2, 1).pow_u32(2);
            let equations = [
                RealPoly::var(2, 0) - y_squared.scale(&Real::from(scale.clone())),
                y_squared - RealPoly::constant(2, Real::from_i64(2)),
            ];
            let bounds = [
                [scale.clone(), &scale * &Rational::from_i64(3)],
                [Rational::one(), Rational::from_i64(2)],
            ];
            let mut dec = Decider::default();
            let [x, y]: [Real; 2] = RealPoly::simple_algebraic_root_in_reconstructed(
                &equations,
                &bounds,
                &[0, 1],
                precision,
                &mut dec,
            )
            .unwrap()
            .try_into()
            .unwrap();
            assert_eq!(
                dec.rational_value(&x).unwrap(),
                Some(&scale * &Rational::from_i64(2))
            );
            assert!(matches!(
                dec.sign(&(y.square() - Real::from_i64(2))),
                Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound
            ));
        }
    });
    assert!(!work.exhausted);
}

#[test]
fn reconstructed_joint_roots_keep_closed_non_dyadic_coordinates_exact() {
    let (_, work) = with_operation_budget(20_000_000, || {
        for fixed in ["1/3", "-2/7"] {
            let fixed = Rational::parse(fixed).unwrap();
            let equations = [
                RealPoly::var(2, 0) - RealPoly::constant(2, Real::from(fixed.clone())),
                RealPoly::var(2, 1).pow_u32(2) - RealPoly::constant(2, Real::from_i64(2)),
            ];
            let bounds = [
                [fixed.clone(), fixed.clone()],
                [Rational::one(), Rational::from_i64(2)],
            ];
            let mut dec = Decider::default();
            let [x, y]: [Real; 2] = RealPoly::simple_algebraic_root_in_reconstructed(
                &equations,
                &bounds,
                &[0, 1],
                128,
                &mut dec,
            )
            .unwrap()
            .try_into()
            .unwrap();
            assert_eq!(dec.rational_value(&x).unwrap(), Some(fixed));
            assert!(matches!(
                dec.sign(&(y.square() - Real::from_i64(2))),
                Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound
            ));
        }
    });
    assert!(!work.exhausted);
}

#[test]
fn selected_coefficient_field_preserves_correlated_quadratic_coordinates() {
    let (_, work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let a = Real::from(2_i64).sqrt();
        let b = Real::from(3_i64).sqrt();
        let x = RealPoly::var(2, 0);
        let y = RealPoly::var(2, 1);
        let equations = [
            &x + &y - RealPoly::constant(2, &a + &b),
            x.pow_u32(2) + y.pow_u32(2).scale(&Real::from(2_i64))
                - RealPoly::constant(2, Real::from(7_i64)),
        ];
        let half_three = Rational::from_ratio(3, 2).unwrap();
        let values = RealPoly::simple_algebraic_root_in_reconstructed(
            &equations,
            &[
                [half_three.clone(), Rational::from_i64(2)],
                [Rational::one(), half_three],
            ],
            &[1, 0, 1],
            128,
            &mut dec,
        )
        .unwrap();
        for (value, expected) in values.iter().zip([&a, &b, &a]) {
            assert!(matches!(
                dec.sign(&(value - expected)),
                Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound
            ));
        }
        for equation in equations {
            assert!(matches!(
                dec.sign(&equation.eval_real(&[values[1].clone(), values[0].clone()])),
                Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound
            ));
        }
    });
    assert!(!work.exhausted);
}
