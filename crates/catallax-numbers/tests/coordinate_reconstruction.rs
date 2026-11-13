use catallax_numbers::{
    BoxRoots, Decider, IsolatedRoot, MPoly, Outcome, Rational, Real, RootBox, ZeroProof,
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
