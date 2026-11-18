//! Independent native-support rolling-ball geometry, under the request work cap.

use catallax_numbers::{Ball, Decider, Rational, Real, with_operation_budget};
use catallax_shapes::{
    BSplineCurve, Cone, Cylinder, NumericalFilletSpine, NumericalFilletSupport, Plane, Point3,
    Surface, Vector3,
};

fn point(x: i64, y: i64, z: i64) -> Point3 {
    Point3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
}

fn axis(x: i64, y: i64, z: i64) -> Vector3 {
    Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
}

fn quarter_circle(height: i64, dec: &mut Decider) -> BSplineCurve {
    BSplineCurve::new(
        2,
        [0, 0, 0, 1, 1, 1].map(Rational::from_i64).to_vec(),
        vec![
            point(4, 0, height),
            point(4, 4, height),
            point(0, 4, height),
        ],
        vec![
            Real::one(),
            Real::from_rational(Rational::from_ratio(1, 2).unwrap()).sqrt(),
            Real::one(),
        ],
        dec,
    )
    .unwrap()
}

fn check_walk(
    build_surface: impl FnOnce(&mut Decider) -> Surface,
    height: i64,
    expected_radius: Real,
    seed: [[Rational; 2]; 2],
) {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let curved = build_surface(&mut dec);
        let plane = Surface::Plane(
            Plane::new(point(0, 0, height), axis(0, 0, 1), axis(1, 0, 0), &mut dec).unwrap(),
        );
        let surfaces = [plane, curved];
        let supports = surfaces.clone().map(|surface| {
            NumericalFilletSupport::new(surface, Rational::from_i64(-1), &mut dec).unwrap()
        });
        let tolerance = Rational::parse("0.0000001").unwrap();
        let walk = NumericalFilletSpine::walk(
            supports,
            quarter_circle(height, &mut dec),
            [Rational::zero(), Rational::one()],
            seed,
            tolerance.clone(),
            &mut dec,
        )
        .unwrap();
        let radius = dec.enclose_at(&expected_radius, 128).unwrap();
        let target_square = radius.square();
        for node in walk.nodes() {
            let centre = node
                .centre
                .each_ref()
                .map(|value| Ball::from_rational(value, 128));
            let radial_square = &centre[0].square() + &centre[1].square();
            let radial_error = &radial_square - &target_square;
            assert!(
                radial_error.lo().abs().max(radial_error.hi().abs())
                    <= &tolerance * &Rational::from_i64(8)
            );
            assert!((&node.centre[2] - &Rational::from_i64(height - 1)).abs() <= tolerance);
            for (surface, uv) in surfaces.iter().zip(&node.parameters) {
                let uv = uv.each_ref().map(|value| Ball::from_rational(value, 128));
                let contact = surface.evaluate(&uv[0], &uv[1], &dec).unwrap();
                let radius_square = (0..3)
                    .fold(Ball::from_rational(&Rational::zero(), 128), |sum, i| {
                        sum + (&contact[i] - &centre[i]).square()
                    });
                let error = &radius_square - &Ball::from_rational(&Rational::one(), 128);
                assert!(
                    error.lo().abs().max(error.hi().abs()) <= &tolerance * &Rational::from_i64(4)
                );
            }
        }
        let first = &walk.nodes().first().unwrap().centre;
        let last = &walk.nodes().last().unwrap().centre;
        assert!((&first[0] - &radius.mid()).abs() <= tolerance);
        assert!(first[1].abs() <= tolerance);
        assert!(last[0].abs() <= tolerance);
        assert!((&last[1] - &radius.mid()).abs() <= tolerance);
        println!("rolling-ball stations={}", walk.nodes().len());
    });
    assert!(!work.exhausted, "{work:?}");
    println!("rolling-ball request {work:?}");
}

#[test]
fn circular_cylinder_rim_keeps_the_selected_quarter_and_unit_ball_contacts() {
    check_walk(
        |dec| {
            Surface::Cylinder(
                Cylinder::new(
                    point(0, 0, 0),
                    axis(0, 0, 1),
                    Real::from_i64(4),
                    axis(1, 0, 0),
                    dec,
                )
                .unwrap(),
            )
        },
        0,
        Real::from_i64(3),
        [
            [Rational::from_i64(3), Rational::zero()],
            [Rational::zero(), Rational::from_i64(-1)],
        ],
    );
}

#[test]
fn conical_rim_uses_the_native_sloped_normal_not_a_cylinder_normal() {
    check_walk(
        |dec| {
            Surface::Cone(
                Cone::new(
                    point(0, 0, 0),
                    axis(0, 0, 1),
                    Real::one(),
                    axis(1, 0, 0),
                    dec,
                )
                .unwrap(),
            )
        },
        4,
        Real::from_i64(3) - Real::from_i64(2).sqrt(),
        [
            [Rational::from_ratio(8, 5).unwrap(), Rational::zero()],
            [Rational::zero(), Rational::from_ratio(23, 10).unwrap()],
        ],
    );
}

#[test]
fn three_scaled_plane_contacts_isolate_one_ball_but_parallel_contacts_do_not() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let planes = [
            (axis(-3, 0, 0), axis(0, 2, 0)),
            (axis(0, -4, 0), axis(0, 0, 7)),
            (axis(0, 0, -5), axis(11, 0, 0)),
        ];
        let supports = planes.map(|(normal, reference)| {
            NumericalFilletSupport::new(
                Surface::Plane(Plane::new(point(0, 0, 0), normal, reference, &mut dec).unwrap()),
                Rational::from_i64(-1),
                &mut dec,
            )
            .unwrap()
        });
        let tolerance = Rational::parse("0.0000001").unwrap();
        let seed = std::array::from_fn(|_| {
            [
                Rational::from_ratio(9, 10).unwrap(),
                Rational::from_ratio(-6, 5).unwrap(),
            ]
        });
        let corner = catallax_shapes::NumericalFilletCorner::solve(
            &supports,
            seed,
            tolerance.clone(),
            &mut dec,
        )
        .unwrap();
        for coordinate in corner.centre() {
            assert!((coordinate - &Rational::one()).abs() <= tolerance);
        }
        for (contact, expected) in corner
            .contacts()
            .iter()
            .zip([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
        {
            for (coordinate, value) in contact.iter().zip(expected) {
                assert!((coordinate - &Rational::from_i64(value)).abs() <= tolerance);
            }
        }
        let parallel = std::array::from_fn(|_| {
            NumericalFilletSupport::new(
                Surface::Plane(
                    Plane::new(point(0, 0, 0), axis(0, 0, -1), axis(1, 0, 0), &mut dec).unwrap(),
                ),
                Rational::from_i64(-1),
                &mut dec,
            )
            .unwrap()
        });
        let failure = catallax_shapes::NumericalFilletCorner::solve(
            &parallel,
            std::array::from_fn(|_| [Rational::one(), Rational::from_i64(-1)]),
            tolerance,
            &mut dec,
        )
        .unwrap_err();
        assert_eq!(failure.code, "fillet.numerical-spine");
    });
    assert!(!work.exhausted, "{work:?}");
    println!("corner request {work:?}");
}
