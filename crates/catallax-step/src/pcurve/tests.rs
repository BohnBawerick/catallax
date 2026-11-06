use super::*;
use catallax_shapes::{Point3, Surface as Carrier, Torus, Vector3};
use std::rc::Rc;
use uv::Parametric;

fn projection() -> Projection {
    Projection::new(
        Rc::new(Surface {
            placement: Rc::new(
                Placement::new(
                    [0.into(), 0.into(), 0.into()],
                    [0.into(), 0.into(), 1.into()],
                    [1.into(), 0.into(), 0.into()],
                )
                .unwrap(),
            ),
            parameters: SurfaceParameters::Torus {
                major: 4.into(),
                minor: 1.into(),
            },
        }),
        &Carrier::Torus(Torus {
            center: Point3::new(0.into(), 0.into(), 0.into()),
            axis: Vector3::new(0.into(), 0.into(), 1.into()),
            ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
            major: 4.into(),
            minor: 1.into(),
        }),
        128,
        17,
        &Decider::default(),
    )
    .unwrap()
}

fn circle(center: [i64; 3], axis: [i64; 3], major: i64, minor: i64) -> Curve {
    Curve::Ellipse {
        placement: Rc::new(
            Placement::new(
                center.map(Rational::from_i64),
                axis.map(Rational::from_i64),
                [1.into(), 0.into(), 0.into()],
            )
            .unwrap(),
        ),
        major: major.into(),
        minor: minor.into(),
    }
}

#[test]
fn latitude_certificate_requires_a_regular_coaxial_written_circle() {
    let projection = projection();
    for (center, axis, major, minor, expected) in [
        ([0, 0, 1], [0, 0, 1], 5, 5, Some(1)),
        ([0, 0, 1], [0, 0, 1], 4, 4, Some(1)),
        ([0, 0, 1], [0, 0, -7], 5, 5, Some(-1)),
        ([0, 0, 0], [0, 0, 1], 3, 3, Some(1)),
        ([1, 0, 1], [0, 0, 1], 5, 5, None),
        ([0, 0, 1], [1, 0, 1], 5, 5, None),
        ([0, 0, 1], [0, 0, 1], 5, 3, None),
        ([0, 0, 1], [0, 0, 1], 0, 0, None),
        ([0, 0, 0], [0, 0, 1], 4, 4, None),
    ] {
        assert_eq!(
            projection.latitude_speed(&circle(center, axis, major, minor)),
            expected
        );
    }
}

#[test]
fn latitude_jets_match_the_written_curve_in_both_directions() {
    let projection = projection();
    let dec = Decider::default();
    for (direction, reference) in [
        (-1, [1, 0, 0]),
        (1, [1, 0, 0]),
        (-1, [1, 1, 0]),
        (1, [1, 1, 0]),
    ] {
        let curve = Curve::Ellipse {
            placement: Rc::new(
                Placement::new(
                    [0.into(), 0.into(), 1.into()],
                    [0.into(), 0.into(), direction.into()],
                    reference.map(Rational::from_i64),
                )
                .unwrap(),
            ),
            major: 5.into(),
            minor: 5.into(),
        };
        let piece = |speed| Piece {
            latitude_speed: speed,
            curve: curve.local(&projection.written.placement, 128),
            guide: Guide::Fixed(scalar(0, 128)),
            projection: &projection,
            range: [(-1).into(), 1.into()],
            endpoints: [None, None],
        };
        let fast = piece(projection.latitude_speed(&curve));
        let generic = piece(None);
        let domain = ball(&(-1).into(), 128).hull(&ball(&1.into(), 128));
        assert!(fast.regular(&domain, &dec).unwrap());
        let derivatives = fast.coefficients::<9>(&domain, false, &dec).unwrap();
        assert_eq!(derivatives[1][0].lo(), direction.into());
        assert_eq!(derivatives[1][0].hi(), direction.into());
        assert!(derivatives[1][1].is_exact_zero());
        assert!(derivatives[2..].iter().flatten().all(Ball::is_exact_zero));
        for numerator in -4..=4 {
            let t = Rational::from_ratio(numerator, 4).unwrap();
            let parameter = ball(&t, 128);
            let ordinary = generic.coefficients::<9>(&parameter, true, &dec).unwrap();
            let exact = fast.coefficients::<9>(&parameter, true, &dec).unwrap();
            assert_eq!(ordinary[0][0].lo(), exact[0][0].lo());
            assert_eq!(ordinary[0][1].hi(), exact[0][1].hi());
            for (order, pair) in ordinary.iter().enumerate().skip(1) {
                for (coordinate, coefficient) in pair.iter().enumerate() {
                    let expected = if order == 1 && coordinate == 0 {
                        direction
                    } else {
                        0
                    };
                    let expected = Rational::from_i64(expected);
                    assert!(coefficient.lo() <= expected);
                    assert!(coefficient.hi() >= expected);
                }
            }
        }
    }
}
