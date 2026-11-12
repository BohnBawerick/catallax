use super::*;
use catallax_shapes::{Point3, Surface as Carrier, Torus, Vector3};
use std::sync::Arc;
use uv::Parametric;

#[test]
fn a_guide_near_half_a_turn_refines_before_selecting_its_winding() {
    let p = 96;
    let dec = Decider::default();
    let pi = Ball::pi(p);
    let guide = Guide::Line {
        ends: [[scalar(1, p), scalar(1, p)], [scalar(1, p), scalar(1, p)]],
    };
    let parameter = ball(&Rational::from_ratio(1, 2).unwrap(), p);
    let range = [Rational::zero(), Rational::one()];
    let x = scalar(-1, p);
    let y = ball(&(-Rational::one() - Rational::one().mul_pow2(-80)), p);
    let base = uv::angle(&y, &x, &pi, Ball::atan_reduced_checked, &dec).unwrap();
    let coarse = guide
        .angle(
            &parameter,
            &range,
            &pi,
            32,
            Ball::atan_reduced_checked,
            &dec,
        )
        .unwrap();
    assert!(shifted_angle(&base, &coarse, &pi).is_none());
    let target = guide
        .angle(&parameter, &range, &pi, p, Ball::atan_reduced_checked, &dec)
        .unwrap();
    let reference = lift(&y, &x, &target, &pi, Ball::atan_reduced_checked, &dec).unwrap();
    let mut projection = projection();
    projection.pi = pi.clone();
    projection.atan = Ball::atan_reduced_checked;
    projection.fillet = true;
    let answer = guide
        .lift(&parameter, &range, [&x, &y], &projection, &dec)
        .unwrap();
    assert_eq!(answer.lo(), reference.lo());
    assert_eq!(answer.hi(), reference.hi());
}

fn projection() -> Projection {
    Projection::new(
        Arc::new(Surface {
            placement: Arc::new(
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
        false,
        Seams::Legacy,
        &Decider::default(),
    )
    .unwrap()
}

fn circle(center: [i64; 3], axis: [i64; 3], major: i64, minor: i64) -> Curve {
    Curve::Ellipse {
        placement: Arc::new(
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
            placement: Arc::new(
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

fn open_conic(hyperbola: bool) -> Curve {
    let placement = Arc::new(
        Placement::new(
            [1.into(), 2.into(), 3.into()],
            [0.into(), 0.into(), 1.into()],
            [1.into(), 0.into(), 0.into()],
        )
        .unwrap(),
    );
    if hyperbola {
        Curve::Hyperbola {
            placement,
            major: 2.into(),
            minor: 3.into(),
        }
    } else {
        Curve::Parabola {
            placement,
            focal: Rational::from_ratio(3, 2).unwrap(),
        }
    }
}

#[test]
fn open_conic_jets_are_the_written_parametrizations() {
    let dec = Decider::default();
    let at = |t: f64| ball(&Rational::from_ratio((t * 8.0) as i64, 8).unwrap(), 128);
    let inside = |ball: &Ball, expected: f64| {
        let (lo, hi) = (ball.lo().to_f64_nearest(), ball.hi().to_f64_nearest());
        assert!(
            lo <= expected + 1e-12 && expected - 1e-12 <= hi,
            "{expected} outside [{lo}, {hi}]"
        );
    };
    let t = 0.625_f64;
    for hyperbola in [true, false] {
        let frame = Placement::new(
            [0.into(), 0.into(), 0.into()],
            [0.into(), 0.into(), 1.into()],
            [1.into(), 0.into(), 0.into()],
        )
        .unwrap();
        let local = open_conic(hyperbola).local(&frame, 128);
        let [x, y, z] = local.coefficients::<3>(&at(t), &dec).unwrap();
        // Value, then first derivative and half the second, per coordinate.
        let (value, first, second) = if hyperbola {
            (
                [1.0 + 2.0 * t.cosh(), 2.0 + 3.0 * t.sinh()],
                [2.0 * t.sinh(), 3.0 * t.cosh()],
                [t.cosh(), 1.5 * t.sinh()],
            )
        } else {
            (
                [1.0 + 1.5 * t * t, 2.0 + 3.0 * t],
                [3.0 * t, 3.0],
                [1.5, 0.0],
            )
        };
        for (jet, axis) in [(&x, 0), (&y, 1)] {
            inside(&jet.coefficients[0], value[axis]);
            inside(&jet.coefficients[1], first[axis]);
            inside(&jet.coefficients[2], second[axis]);
        }
        inside(&z.coefficients[0], 3.0);
        assert!(z.coefficients[1].is_exact_zero());
    }
}

#[test]
fn open_conic_ranges_run_between_the_written_ends_and_refuse_a_reversed_pair() {
    let pi = Ball::pi(96);
    let dec = Decider::default();
    let range = |hyperbola, ends: [R3; 2]| {
        open_conic(hyperbola).analytic_range(
            [&ends[0], &ends[1]],
            false,
            12,
            &pi,
            Ball::atan_checked,
            &dec,
        )
    };
    let point = |x: f64, y: f64| {
        let exact = |v: f64| Rational::from_ratio((v * 1e6).round() as i64, 1_000_000).unwrap();
        [exact(x), exact(y), Rational::from_i64(3)]
    };
    // A parabola `(1 + 1.5 t^2, 2 + 3 t)` between t = -1/2 and t = 1.
    let [lo, hi] = range(false, [point(1.375, 0.5), point(2.5, 5.0)]).unwrap();
    assert!((lo.to_f64_nearest() + 0.5).abs() < 1e-5);
    assert!((hi.to_f64_nearest() - 1.0).abs() < 1e-5);
    // A hyperbola `(1 + 2 cosh t, 2 + 3 sinh t)` between t = -1 and t = 1/2.
    let (c, s) = (1.0_f64.cosh(), 1.0_f64.sinh());
    let (c2, s2) = (0.5_f64.cosh(), 0.5_f64.sinh());
    let [lo, hi] = range(
        true,
        [
            point(1.0 + 2.0 * c, 2.0 - 3.0 * s),
            point(1.0 + 2.0 * c2, 2.0 + 3.0 * s2),
        ],
    )
    .unwrap();
    assert!((lo.to_f64_nearest() + 1.0).abs() < 1e-5);
    assert!((hi.to_f64_nearest() - 0.5).abs() < 1e-5);
    // Ends against the written direction have no parameter extent.
    assert!(range(false, [point(2.5, 5.0), point(1.375, 0.5)]).is_err());
    assert!(
        range(
            true,
            [
                point(1.0 + 2.0 * c2, 2.0 + 3.0 * s2),
                point(1.0 + 2.0 * c, 2.0 - 3.0 * s)
            ]
        )
        .is_err()
    );
}

#[test]
fn an_elliptic_chart_runs_on_the_swept_ellipses_own_parameter() {
    let p = 128;
    let dec = Decider::default();
    let pi = Ball::pi(p);
    let map = SurfaceMap::EllipticCylinder {
        inverse_major: ball(&Rational::from_ratio(1, 3).unwrap(), p),
        inverse_minor: ball(&Rational::from_ratio(1, 2).unwrap(), p),
    };
    // The point at parameter 3/4 of an ellipse with semi-axes 3 and 2, at height 5.
    let t = Rational::from_ratio(3, 4).unwrap();
    let angle = ball(&t, p);
    let coordinates = [
        Jet::<1>::constant(&angle.cos() * &scalar(3, p)),
        Jet::constant(&angle.sin() * &scalar(2, p)),
        Jet::constant(scalar(5, p)),
    ];
    assert!(map.regular(&coordinates));
    // The circular angle of this point is not the parameter.
    let [x, y] = map.angle_xy(
        &coordinates[0].coefficients[0],
        &coordinates[1].coefficients[0],
    );
    let [u, v] = map
        .inverse(
            &coordinates,
            angle.clone(),
            &pi,
            true,
            Ball::atan_checked,
            &dec,
        )
        .unwrap()[0]
        .clone();
    assert!(u.lo() <= t && t <= u.hi());
    assert!(v.lo() <= Rational::from_i64(5) && Rational::from_i64(5) <= v.hi());
    let unscaled = uv::angle(
        &coordinates[1].coefficients[0],
        &coordinates[0].coefficients[0],
        &pi,
        Ball::atan_checked,
        &dec,
    )
    .unwrap();
    assert!(number::magnitude(&(&unscaled - &angle)) > Rational::from_ratio(1, 100).unwrap());
    let scaled = uv::angle(&y, &x, &pi, Ball::atan_checked, &dec).unwrap();
    assert!(number::magnitude(&(&scaled - &angle)) < Rational::from_ratio(1, 1_000_000).unwrap());
}
