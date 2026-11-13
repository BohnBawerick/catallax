use super::*;

fn q(value: &str) -> Rational {
    value.parse().unwrap()
}

fn square(lo: &str, hi: &str) -> Domain {
    let interval = [q(lo), q(hi)];
    [interval.clone(), interval]
}

fn diagonal_roots(polynomial: RealPoly, domain: Domain, dec: &mut Decider) -> Vec<[Real; 2]> {
    let diagonal = RealPoly::var(2, 1) - RealPoly::var(2, 0);
    linear_roots(&[polynomial, diagonal], &domain, dec)
        .unwrap()
        .unwrap()
}

fn cubic() -> RealPoly {
    RealPoly::var(2, 0).pow_u32(3) - RealPoly::constant(2, Real::from_i64(2).sqrt())
}

#[test]
fn a_monotone_algebraic_cubic_keeps_its_exact_root() {
    let mut dec = Decider::default();
    let before = dec.stats().coincident_separation;
    let roots = diagonal_roots(cubic(), square("1", "3/2"), &mut dec);
    assert_eq!(roots.len(), 1);
    assert_eq!(
        exact_sign(&(&roots[0][0] - &roots[0][1]), &mut dec).unwrap(),
        Sign::Zero
    );
    assert_eq!(
        exact_sign(
            &(&roots[0][0].pow_u32(3) - &Real::from_i64(2).sqrt()),
            &mut dec
        )
        .unwrap(),
        Sign::Zero
    );
    assert!(within(&roots[0][0], &[q("1"), q("3/2")], &mut dec).unwrap());
    assert_eq!(dec.stats().coincident_separation, before);
}

#[test]
fn a_monotone_root_on_a_closed_cell_end_is_not_discarded() {
    let u = RealPoly::var(2, 0);
    let polynomial = (u.clone() - RealPoly::constant(2, Real::one()))
        * (u.pow_u32(2) + RealPoly::constant(2, Real::from_i64(2).sqrt()));
    let mut dec = Decider::default();
    let roots = diagonal_roots(polynomial, square("1", "3/2"), &mut dec);
    assert_eq!(roots.len(), 1);
    for coordinate in &roots[0] {
        assert_eq!(
            exact_sign(&(coordinate - &Real::one()), &mut dec).unwrap(),
            Sign::Zero
        );
    }
}

#[test]
fn a_root_outside_the_recovered_coordinate_is_not_a_meeting() {
    let roots = diagonal_roots(
        cubic(),
        [[q("1"), q("3/2")], [q("0"), q("1")]],
        &mut Decider::default(),
    );
    assert!(roots.is_empty());
}

#[test]
fn nonmonotone_cells_keep_all_their_distinct_crossings() {
    let u = RealPoly::var(2, 0);
    let polynomial = u.pow_u32(3) - u.scale(&Real::from_i64(2));
    let mut dec = Decider::default();
    let roots = diagonal_roots(polynomial, square("-2", "2"), &mut dec);
    let sqrt = Real::from_i64(2).sqrt();
    let expected = [-&sqrt, Real::zero(), sqrt];
    assert_eq!(roots.len(), expected.len());
    for (root, expected) in roots.iter().zip(expected) {
        for coordinate in root {
            assert_eq!(
                exact_sign(&(coordinate - &expected), &mut dec).unwrap(),
                Sign::Zero
            );
        }
    }
}

#[test]
fn nonmonotone_algebraic_cells_keep_each_certified_crossing() {
    let u = RealPoly::var(2, 0);
    let polynomial = u.pow_u32(3) - u.scale(&Real::from_i64(2).sqrt());
    let mut dec = Decider::default();
    let roots = diagonal_roots(polynomial.clone(), square("-2", "2"), &mut dec);
    let intervals = [
        [q("-6/5"), q("-59/50")],
        [q("-1/100"), q("1/100")],
        [q("59/50"), q("6/5")],
    ];
    assert_eq!(roots.len(), intervals.len());
    for (root, interval) in roots.iter().zip(intervals) {
        for coordinate in root {
            assert!(within(coordinate, &interval, &mut dec).unwrap());
        }
        assert_eq!(
            exact_sign(&polynomial.eval_real(root), &mut dec).unwrap(),
            Sign::Zero
        );
        assert_eq!(
            exact_sign(&(&root[0] - &root[1]), &mut dec).unwrap(),
            Sign::Zero
        );
    }
    assert_eq!(dec.stats().coincident_separation, 0);
}

#[test]
fn algebraic_cells_do_not_duplicate_shared_rational_endpoints() {
    let u = RealPoly::var(2, 0);
    let polynomial = (u.pow_u32(3) - u).scale(&Real::from_i64(2).sqrt());
    let mut dec = Decider::default();
    let roots = diagonal_roots(polynomial, square("-2", "2"), &mut dec);
    let expected = [-1, 0, 1];
    assert_eq!(roots.len(), expected.len());
    for (root, expected) in roots.iter().zip(expected) {
        for coordinate in root {
            assert_eq!(
                exact_sign(&(coordinate - &Real::from_i64(expected)), &mut dec).unwrap(),
                Sign::Zero
            );
        }
    }
    assert_eq!(dec.stats().coincident_separation, 0);
}
