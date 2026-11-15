//! Complete finite point sets survive degree, weight and sense differences.

use catallax_intersect::same_curve;
use catallax_numbers::{Decider, Outcome, Rational, Real, with_operation_budget};
use catallax_shapes::{Apart, BSplineCurve, Curve, GeomStore, Point3};

fn curve(points: &[[i64; 2]], weights: &[i64], dec: &mut Decider) -> BSplineCurve {
    let degree = points.len() - 1;
    BSplineCurve::new(
        degree,
        std::iter::repeat_n(Rational::zero(), degree + 1)
            .chain(std::iter::repeat_n(Rational::one(), degree + 1))
            .collect(),
        points
            .iter()
            .map(|[x, y]| Point3::new(Real::from_i64(*x), Real::from_i64(*y), Real::zero()))
            .collect(),
        weights
            .iter()
            .map(|weight| Real::from_i64(*weight))
            .collect(),
        dec,
    )
    .unwrap()
}

#[test]
fn finite_straight_point_sets_are_equal_only_with_the_same_complete_interval() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let mut geom = GeomStore::new();
        let a = geom.add_curve(Curve::BSpline(curve(
            &[[0, 0], [1, 0], [2, 0]],
            &[1, 3, 1],
            &mut dec,
        )));
        for (points, equal) in [
            (&[[2, 0], [0, 0]][..], true),
            (&[[0, 0], [1, 0]][..], false),
        ] {
            let b = geom.add_curve(Curve::BSpline(curve(points, &[1, 1], &mut dec)));
            for (a, b) in [(a, b), (b, a)] {
                let outcome = same_curve(&geom, a, b, &mut dec);
                if equal {
                    assert!(outcome.is_coincident());
                } else {
                    assert!(matches!(outcome, Outcome::Proven(Apart)));
                }
            }
        }
        let bent = geom.add_curve(Curve::BSpline(curve(
            &[[0, 0], [1, 1], [2, 0]],
            &[1, 1, 1],
            &mut dec,
        )));
        assert!(matches!(
            same_curve(&geom, a, bent, &mut dec),
            Outcome::Proven(Apart)
        ));
        let overshoot = geom.add_curve(Curve::BSpline(curve(
            &[[0, 0], [3, 0], [2, 0]],
            &[1, 1, 1],
            &mut dec,
        )));
        assert!(matches!(
            same_curve(&geom, a, overshoot, &mut dec),
            Outcome::Refused(_)
        ));
    });
    assert!(!work.exhausted, "{} units", work.used);
}
