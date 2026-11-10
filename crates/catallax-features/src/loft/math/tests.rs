use super::*;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn half_root() -> Real {
    Real::from_rational(q("1/2")).sqrt()
}

/// The equal-angle quarter circle of RECIPE 4.12 at one station height.
fn quarter(radius: i64, height: i64) -> Vec<H> {
    let w = half_root();
    let r = Real::from_i64(radius);
    let z = Real::from_i64(height);
    vec![
        [r.clone(), Real::zero(), z.clone(), Real::one()],
        [&r * &w, &r * &w, &z * &w, w.clone()],
        [Real::zero(), r, z, Real::one()],
    ]
}

fn point(net: &[H], t: &Real) -> [Real; 3] {
    let s = &Real::one() - t;
    let basis = [s.square(), &(&s * t) * &Real::from_i64(2), t.square()];
    let h: [Real; 4] = std::array::from_fn(|axis| {
        basis
            .iter()
            .zip(net)
            .fold(Real::zero(), |sum, (b, control)| sum + b * &control[axis])
    });
    std::array::from_fn(|axis| &h[axis] / &h[3])
}

fn coincide(a: &Real, b: &Real, dec: &mut Decider) -> bool {
    matches!(dec.sign(&(a - b)), Outcome::ProvenCoincident(_))
}

#[test]
fn one_common_mobius_map_makes_quarter_circles_rational_without_moving_them() {
    let original = vec![quarter(1, 0), quarter(2, 1), quarter(1, 2)];
    let mut inputs = original.clone();
    let mut dec = Decider::default();
    rational_interval(&mut inputs, &mut dec).unwrap();
    for (old, new) in original.iter().zip(&inputs) {
        for (old, new) in old.iter().zip(new) {
            assert!(new.iter().all(|x| x.rational_ref().is_some()));
            // Euclidean controls are unchanged.
            for axis in 0..3 {
                assert!(coincide(
                    &(&old[axis] / &old[3]),
                    &(&new[axis] / &new[3]),
                    &mut dec
                ));
            }
        }
        let weights: Vec<_> = new
            .iter()
            .map(|control| control[3].rational_ref().unwrap().clone())
            .collect();
        assert_eq!(weights, [q("1"), q("1/2"), q("1/2")]);
    }
    // The new chart at t is the old chart at lambda t / ((1 - t) + lambda t),
    // with the same lambda for every section.
    let lambda = half_root();
    for t in ["0", "1/3", "1/2", "4/5", "1"] {
        let t = Real::from_rational(q(t));
        let mapped = &(&lambda * &t) / &(&(&Real::one() - &t) + &(&lambda * &t));
        for (old, new) in original.iter().zip(&inputs) {
            let expected = point(old, &mapped);
            let actual = point(new, &t);
            for (a, e) in actual.iter().zip(&expected) {
                assert!(coincide(a, e, &mut dec));
            }
        }
    }
}

#[test]
fn intervals_without_a_common_rational_map_are_unchanged() {
    // One map cannot make a quarter circle and a straight quadratic chart
    // rational together, so neither is reparameterized.
    let straight: Vec<H> = [["2", "0"], ["1", "1"], ["0", "2"]]
        .iter()
        .map(|[x, y]| {
            [
                Real::from_rational(q(x)),
                Real::from_rational(q(y)),
                Real::one(),
                Real::one(),
            ]
        })
        .collect();
    let original = vec![quarter(1, 0), straight];
    let mut inputs = original.clone();
    rational_interval(&mut inputs, &mut Decider::default()).unwrap();
    for (old, new) in original.iter().flatten().zip(inputs.iter().flatten()) {
        for (old, new) in old.iter().zip(new) {
            assert!(old.same_definition(new));
        }
    }
}
