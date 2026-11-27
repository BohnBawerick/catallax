use super::*;
use crate::Sense;
use crate::numerical::CopyTag;
use catallax_numbers::Real;
use catallax_shapes::{Plane, Point3, Vector3};

fn q(value: i64) -> Rational {
    Rational::from(value)
}

fn point(x: i64, y: i64, z: i64) -> Point3 {
    Point3 {
        x: Real::from(x),
        y: Real::from(y),
        z: Real::from(z),
    }
}

#[test]
fn a_gap_away_from_the_midpoint_is_not_admitted_as_a_certificate() {
    let mut dec = Decider::default();
    let spline = BSplineCurve::new(
        3,
        [0, 0, 0, 0, 1, 1, 1, 1].map(q).to_vec(),
        vec![
            point(0, 0, 0),
            point(1, 0, 0),
            point(2, 0, 1),
            point(3, 0, 0),
        ],
        vec![Real::from(1); 4],
        &mut dec,
    )
    .unwrap();
    let edge = NumericalEdge {
        key: "bump".into(),
        curve: Curve::BSpline(spline),
        vertices: [0, 1],
        parameters: [q(0), q(1)],
    };
    let face = NumericalFace {
        key: "plane".into(),
        surface: Surface::Plane(Plane {
            origin: point(0, 0, 0),
            normal: Vector3 {
                x: Real::from(0),
                y: Real::from(0),
                z: Real::from(1),
            },
            ref_dir: Vector3 {
                x: Real::from(1),
                y: Real::from(0),
                z: Real::from(0),
            },
        }),
        sense: Sense::Same,
        loops: Vec::new(),
    };
    let trim = NumericalTrim {
        edge: 0,
        sense: Sense::Reversed,
        uv: [[q(3), q(0)], [q(0), q(0)]],
        tag: CopyTag::Independent,
    };
    let cancel = AtomicBool::new(false);
    let prove = |limit: &Rational, left| {
        prove_use(
            &edge,
            &face,
            &trim,
            None,
            &(limit * limit),
            &mut Work {
                left,
                cancel: &cancel,
            },
            &dec,
        )
    };
    // The gap is 3t²(1-t): zero at the endpoints, 3/8 at the midpoint,
    // and 4/9 at t=2/3. The midpoint and endpoints alone miss the violation.
    assert_eq!(
        prove(&(&q(2) / &q(5)), 256),
        Err(EdgeProofFailure::AboveLimit)
    );
    let limit = &q(9) / &q(20);
    assert_eq!(prove(&limit, 3), Err(EdgeProofFailure::WorkCap));
    let upper = prove(&limit, 256).unwrap();
    assert!(upper >= &q(16) / &q(81));
    assert!(upper <= &limit * &limit);
}
