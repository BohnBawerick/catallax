use catallax_numbers::{Decider, Outcome, Rational, Real, with_operation_budget};
use catallax_shapes::{BSplineCurve, Curve, GeomStore, Point3};

fn roof(shift: [Real; 3], weight: Real, dec: &mut Decider) -> BSplineCurve {
    BSplineCurve::new(
        2,
        [0, 0, 0, 1, 1, 2, 2, 2].map(Rational::from_i64).to_vec(),
        [[0, 0], [1, 2], [2, 2], [3, 2], [4, 0]]
            .into_iter()
            .map(|[x, y]| {
                Point3::new(
                    &Real::from_i64(x) + &shift[0],
                    &Real::from_i64(y) + &shift[1],
                    shift[2].clone(),
                )
            })
            .collect(),
        vec![
            Real::one(),
            weight.clone(),
            Real::one(),
            weight,
            Real::one(),
        ],
        dec,
    )
    .unwrap()
}

fn shift(x: i64, z: i64) -> [Real; 3] {
    [Real::from_i64(x), Real::zero(), Real::from_i64(z)]
}

#[test]
fn separated_native_curves_use_the_patch_engine_exclusion_certificate() {
    for weight in [Real::one(), Real::from_i64(2).sqrt()] {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let first = roof(shift(0, 0), weight.clone(), &mut dec);
            let second = roof(shift(0, 2), weight, &mut dec);
            assert!(
                first
                    .proves_control_hull_separation(&second, &mut dec)
                    .unwrap()
            );
            assert!(
                second
                    .proves_control_hull_separation(&first, &mut dec)
                    .unwrap()
            );
            let mut geom = GeomStore::new();
            let first = geom.add_curve(Curve::BSpline(first));
            let second = geom.add_curve(Curve::BSpline(second));
            for (a, b) in [(first, second), (second, first)] {
                assert!(matches!(
                    crate::same_curve(&geom, a, b, &mut dec),
                    Outcome::Proven(_)
                ));
            }
        });
        assert!(!work.exhausted, "{} units", work.used);
    }
}

#[test]
fn touching_hulls_and_unsupported_controls_supply_no_separation_proof() {
    let ((first, separated), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let first = roof(shift(0, 0), Real::one(), &mut dec);
        let touching = roof(shift(4, 0), Real::one(), &mut dec);
        assert!(
            !first
                .proves_control_hull_separation(&touching, &mut dec)
                .unwrap()
        );
        assert!(
            !touching
                .proves_control_hull_separation(&first, &mut dec)
                .unwrap()
        );
        assert!(
            !first
                .proves_control_hull_separation(&first, &mut dec)
                .unwrap()
        );
        let irrational = roof(
            [Real::zero(), Real::zero(), Real::from_i64(2).sqrt()],
            Real::one(),
            &mut dec,
        );
        assert!(
            !first
                .proves_control_hull_separation(&irrational, &mut dec)
                .unwrap()
        );
        let separated = roof(shift(0, 2), Real::one(), &mut dec);
        (first, separated)
    });
    assert!(!work.exhausted, "{} units", work.used);
    let (answer, stopped) = with_operation_budget(0, || {
        first.proves_control_hull_separation(&separated, &mut Decider::default())
    });
    assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
    assert!(stopped.exhausted);
}
