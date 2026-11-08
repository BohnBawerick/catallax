#![allow(missing_docs)]

use catallax_numbers::{Ball, Decider, Outcome, Rational, Real, StrictSign, ZeroProof};
use catallax_shapes::{
    BSplineCurve, BezierExtrusion, Pose, RPoint3, RVector3, Rotation, SplineBoundary,
};

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
        assert_eq!(
            extrusion
                .controls()
                .iter()
                .map(|p| p.as_rational().unwrap())
                .collect::<Vec<_>>(),
            cubic()
        );
        assert_eq!(extrusion.direction().as_rational(), Some(vector(height)));
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

#[test]
fn normals_and_shared_edges_follow_a_rigid_pose() {
    let mut dec = Decider::default();
    let original = BezierExtrusion::new(cubic(), vector(2), &mut dec).unwrap();
    let normal = original
        .normal_at(&point(2, 3, 1).to_real(), &mut dec)
        .unwrap();
    assert_eq!(
        normal.as_rational(),
        Some(RVector3::new(r(0), r(-12), r(0)))
    );
    assert!(original.shares_definition_nodes(&original.clone()));
    let pose = Pose::from_rational(
        Rotation::from_quaternion(r(1), r(0), r(0), r(1)).unwrap(),
        RVector3::new(r(7), r(11), r(13)),
    );
    let moved = original.transformed(&pose);
    let at = pose.apply_point(&point(2, 3, 1).to_real());
    let normal = moved.normal_at(&at, &mut dec).unwrap();
    assert_eq!(normal.as_rational(), Some(RVector3::new(r(12), r(0), r(0))));
    let half = Ball::from_rational(&Rational::parse("1/2").unwrap(), 128);
    let middle = moved.carrier().evaluate(&half, &half, &dec).unwrap();
    for (coordinate, expected) in middle.iter().zip([r(4), r(13), r(14)]) {
        assert!(coordinate.contains(&expected));
    }
    assert!(
        moved
            .lower_edge()
            .is_boundary_of(moved.carrier(), SplineBoundary::VMin)
    );
    assert!(
        moved
            .upper_edge()
            .is_boundary_of(moved.carrier(), SplineBoundary::VMax)
    );
    assert!(!original.shares_definition_nodes(&moved));
}

#[test]
fn normal_queries_refuse_points_off_the_finite_patch() {
    let mut dec = Decider::default();
    let extrusion = BezierExtrusion::new(cubic(), vector(2), &mut dec).unwrap();
    for at in [point(2, 3, -1), point(2, 3, 3), point(2, 1, 1)] {
        let error = extrusion.normal_at(&at.to_real(), &mut dec).unwrap_err();
        assert_eq!(error.code, "shapes.bezier-extrusion-domain");
    }
}

#[test]
fn line_crossings_are_exact_finite_and_follow_line_order() {
    let mut dec = Decider::default();
    let extrusion = BezierExtrusion::new(cubic(), vector(2), &mut dec).unwrap();
    let y = Rational::parse("3/2").unwrap();
    for dx in [1, -1] {
        let origin = RPoint3::new(r(-10 * dx), y.clone(), r(1)).to_real();
        let direction = RVector3::new(r(dx), r(0), r(0)).to_real();
        let hits = extrusion.line_hits(&origin, &direction, &mut dec).unwrap();
        assert_eq!(hits.len(), 2);
        for hit in &hits {
            let delta = &hit.x - &Real::from_i64(2);
            let identity = &delta.square() - &Real::from(Rational::parse("25/8").unwrap());
            assert!(
                matches!(dec.sign(&identity), Outcome::ProvenCoincident(proof)
                if proof.proof != ZeroProof::SeparationBound)
            );
            for identity in [&hit.y - &Real::from(y.clone()), &hit.z - &Real::one()] {
                assert!(
                    matches!(dec.sign(&identity), Outcome::ProvenCoincident(proof)
                    if proof.proof != ZeroProof::SeparationBound)
                );
            }
        }
        let order = dec.sign(&(&hits[1].x - &hits[0].x));
        assert!(
            matches!(order, Outcome::Proven(sign) if sign == if dx > 0 { StrictSign::Positive } else { StrictSign::Negative })
        );
    }
    for (origin, direction) in [
        (point(-10, 4, 1), RVector3::new(r(1), r(0), r(0))),
        (
            RPoint3::new(r(-10), y, r(3)),
            RVector3::new(r(1), r(0), r(0)),
        ),
        (point(2, 2, -10), vector(1)),
    ] {
        assert!(
            extrusion
                .line_hits(&origin.to_real(), &direction.to_real(), &mut dec)
                .unwrap()
                .is_empty()
        );
    }
}

#[test]
fn tangent_and_contained_lines_keep_named_refusals() {
    let mut dec = Decider::default();
    let extrusion = BezierExtrusion::new(cubic(), vector(2), &mut dec).unwrap();
    for (origin, direction) in [
        (point(-10, 3, 1), RVector3::new(r(1), r(0), r(0))),
        (point(2, 3, -10), vector(1)),
    ] {
        let error = extrusion
            .line_hits(&origin.to_real(), &direction.to_real(), &mut dec)
            .unwrap_err();
        assert_eq!(error.code, "intersect.parametric-surface");
    }
}
