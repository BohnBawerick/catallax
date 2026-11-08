#![allow(missing_docs)]

use catallax_numbers::{Ball, Decider, Rational};
use catallax_shapes::{BSplineCurve, BezierExtrusion, RPoint3, RVector3, SplineBoundary};

fn r(value: i64) -> Rational {
    Rational::from_i64(value)
}
fn point(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(r(x), r(y), r(z))
}
fn vector(z: i64) -> RVector3 {
    RVector3::new(r(0), r(0), r(z))
}
fn cubic() -> Vec<RPoint3> {
    vec![
        point(0, 0, 0),
        point(0, 4, 0),
        point(4, 4, 0),
        point(4, 0, 0),
    ]
}

#[test]
fn cubic_arch_keeps_exact_shared_boundaries_and_signed_height() {
    for height in [-2, 2] {
        let mut dec = Decider::default();
        let extrusion = BezierExtrusion::new(cubic(), vector(height), &mut dec).unwrap();
        let half = Ball::from_rational(&Rational::parse("1/2").unwrap(), 128);
        let middle = extrusion.carrier().evaluate(&half, &half, &dec).unwrap();
        for (coordinate, expected) in middle.iter().zip([r(2), r(3), r(height / 2)]) {
            assert!(coordinate.contains(&expected));
        }
        let lower = BSplineCurve::boundary(extrusion.carrier(), SplineBoundary::VMin);
        let upper = BSplineCurve::boundary(extrusion.carrier(), SplineBoundary::VMax);
        assert!(lower.is_boundary_of(extrusion.carrier(), SplineBoundary::VMin));
        assert!(upper.is_boundary_of(extrusion.carrier(), SplineBoundary::VMax));
        assert!(lower.evaluate(&half, &dec).unwrap()[2].contains(&r(0)));
        assert!(upper.evaluate(&half, &dec).unwrap()[2].contains(&r(height)));
        assert_eq!(extrusion.controls(), cubic());
        assert_eq!(extrusion.direction(), &vector(height));
    }
}

#[test]
fn quadratic_and_reversed_cubic_spans_are_regular() {
    let mut reversed = cubic();
    reversed.reverse();
    for controls in [
        vec![point(0, 0, 0), point(2, 4, 0), point(4, 0, 0)],
        reversed,
    ] {
        let extrusion = BezierExtrusion::new(controls, vector(2), &mut Decider::default()).unwrap();
        let half = Ball::from_rational(&Rational::parse("1/2").unwrap(), 128);
        let jet = extrusion
            .carrier()
            .jet(&half, &half, &Decider::default())
            .unwrap();
        assert!(!jet.first[0][0].contains(&Rational::zero()));
        assert!(jet.first[1][2].contains(&r(2)));
    }
}

#[test]
fn invalid_or_unproved_spans_refuse_by_name() {
    let mut folded = cubic();
    folded[1] = point(-2, 4, 0);
    let mut singular = cubic();
    singular[1] = singular[0].clone();
    let mut nonplanar = cubic();
    nonplanar[1] = point(0, 4, 1);
    for (controls, direction, code) in [
        (vec![], vector(2), "shapes.bezier-extrusion-controls"),
        (cubic(), vector(0), "shapes.bezier-extrusion-direction"),
        (nonplanar, vector(2), "shapes.bezier-extrusion-plane"),
        (singular, vector(2), "shapes.bezier-extrusion-singular"),
        (folded, vector(2), "shapes.bezier-extrusion-embedding"),
    ] {
        let error = BezierExtrusion::new(controls, direction, &mut Decider::default()).unwrap_err();
        assert_eq!(error.code, code);
    }
}
