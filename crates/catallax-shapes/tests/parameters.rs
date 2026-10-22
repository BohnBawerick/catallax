#![allow(missing_docs)]
mod common;
use catallax_numbers::*;
use catallax_shapes::*;
use common::*;

blind_test!(cases_37_to_41_surface_values_domains_and_derivatives, {
    for bits in [64, 160, 1000] {
        let zero = ball("0", bits);
        let pi = Ball::pi(bits);
        let half = &pi * &ball("1/2", bits);
        let cylinder = Surface::Cylinder(
            Cylinder::new(
                p("0", "0", "0").into(),
                v("0", "0", "1").into(),
                r("2").into(),
                v("3", "0", "0").into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        );
        let sphere = Surface::Sphere(
            Sphere::new(
                p("1", "2", "3").into(),
                r("2").into(),
                v("0", "0", "1").into(),
                v("1", "0", "0").into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        );
        let plane = Surface::Plane(
            Plane::new(
                p("1", "1", "1").into(),
                v("0", "0", "1").into(),
                v("0", "2", "0").into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        );
        let reversed = Surface::Cone(
            Cone::new(
                p("0", "0", "0").into(),
                v("0", "0", "-1").into(),
                r("1/2").into(),
                v("1", "0", "0").into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        );
        for (s, u, v, expected) in [
            (
                cylinder.clone(),
                zero.clone(),
                ball("3", bits),
                p("2", "0", "3"),
            ),
            (
                cylinder.clone(),
                half.clone(),
                ball("3", bits),
                p("0", "2", "3"),
            ),
            (
                cylinder.clone(),
                pi.clone(),
                ball("-1", bits),
                p("-2", "0", "-1"),
            ),
            (sphere.clone(), zero.clone(), zero.clone(), p("3", "2", "3")),
            (torus(), zero.clone(), zero.clone(), p("4", "0", "0")),
            (torus(), zero.clone(), pi.clone(), p("2", "0", "0")),
            (torus(), half.clone(), half.clone(), p("0", "3", "1")),
            (cone(), zero.clone(), ball("10", bits), p("5", "0", "10")),
            (cone(), zero.clone(), ball("-10", bits), p("-5", "0", "-10")),
            (reversed, half.clone(), ball("2", bits), p("0", "-1", "-2")),
            (
                plane.clone(),
                ball("1", bits),
                ball("1", bits),
                p("0", "2", "1"),
            ),
            // Outside the advisory sphere latitude domain, cos(pi)=-1, sin(pi)=0.
            (sphere.clone(), zero.clone(), pi.clone(), p("-1", "2", "3")),
        ] {
            contains(
                &s.evaluate(&u, &v, &catallax_numbers::Decider::default())
                    .unwrap(),
                &expected,
            );
        }
        for u in [zero.clone(), half.clone(), pi.clone(), ball("10", bits)] {
            contains(
                &sphere
                    .evaluate(&u, &half, &catallax_numbers::Decider::default())
                    .unwrap(),
                &p("1", "2", "5"),
            );
            contains(
                &sphere
                    .evaluate(&u, &(-half.clone()), &catallax_numbers::Decider::default())
                    .unwrap(),
                &p("1", "2", "1"),
            );
            contains(
                &cone()
                    .evaluate(&u, &zero, &catallax_numbers::Decider::default())
                    .unwrap(),
                &p("0", "0", "0"),
            );
        }
        for (s, u, v, du, dv) in [
            (
                plane,
                zero.clone(),
                zero.clone(),
                p("0", "1", "0"),
                p("-1", "0", "0"),
            ),
            (
                cylinder,
                zero.clone(),
                ball("7", bits),
                p("0", "2", "0"),
                p("0", "0", "1"),
            ),
            (
                sphere,
                zero.clone(),
                zero.clone(),
                p("0", "2", "0"),
                p("0", "0", "2"),
            ),
            (
                cone(),
                zero.clone(),
                ball("10", bits),
                p("0", "5", "0"),
                p("1/2", "0", "1"),
            ),
            (
                cone(),
                half.clone(),
                zero.clone(),
                p("0", "0", "0"),
                p("0", "1/2", "1"),
            ),
            (
                torus(),
                zero.clone(),
                zero.clone(),
                p("0", "4", "0"),
                p("0", "0", "1"),
            ),
            (
                torus(),
                half.clone(),
                half.clone(),
                p("-3", "0", "0"),
                p("0", "-1", "0"),
            ),
        ] {
            let d = s
                .derivatives(&u, &v, &catallax_numbers::Decider::default())
                .unwrap();
            contains(&d[0], &du);
            contains(&d[1], &dv);
        }
    }
    let expected = [
        (SurfaceKind::Plane, [ParamKind::Line, ParamKind::Line]),
        (
            SurfaceKind::Cylinder,
            [ParamKind::Periodic, ParamKind::Line],
        ),
        (SurfaceKind::Cone, [ParamKind::Periodic, ParamKind::Ray]),
        (
            SurfaceKind::Sphere,
            [ParamKind::Periodic, ParamKind::Latitude],
        ),
        (
            SurfaceKind::Torus,
            [ParamKind::Periodic, ParamKind::Periodic],
        ),
    ];
    for (s, (kind, domain)) in surfaces().into_iter().zip(expected) {
        assert_eq!(s.kind(), kind);
        assert_eq!(s.domain(), domain);
        for b in s
            .evaluate(
                &ball("0", 64),
                &ball("0", 160),
                &catallax_numbers::Decider::default(),
            )
            .unwrap()
        {
            assert_eq!(b.precision(), 160);
        }
        for derivative in s
            .derivatives(
                &ball("0", 160),
                &ball("0", 64),
                &catallax_numbers::Decider::default(),
            )
            .unwrap()
        {
            for b in derivative {
                assert_eq!(b.precision(), 160);
            }
        }
    }
});

blind_test!(cases_42_44_curves_and_exact_parameter_origins, {
    let store = GeomStore::new();
    let bits = 160;
    let half = &Ball::pi(bits) * &ball("1/2", bits);
    let a = Curve::Circle(
        Circle::new(
            p("0", "0", "0").into(),
            v("0", "0", "1").into(),
            r("4").into(),
            v("0", "2", "0").into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    );
    let b = Curve::Circle(
        Circle::new(
            p("0", "0", "0").into(),
            v("0", "0", "1").into(),
            r("2").into(),
            v("1", "1", "0").into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    );
    for (c, t, expected) in [
        (a.clone(), ball("0", bits), p("0", "2", "0")),
        (a.clone(), half.clone(), p("-2", "0", "0")),
        (b.clone(), ball("0", bits), p("1", "1", "0")),
        (ellipse(), ball("0", bits), p("5", "0", "0")),
        (ellipse(), half, p("0", "3", "0")),
        (line(), ball("2", bits), p("5", "0", "3")),
    ] {
        contains(
            &c.evaluate(&t, &catallax_numbers::Decider::default())
                .unwrap()
                .unwrap(),
            &expected,
        );
    }
    for (c, expected) in [
        (a, p("0", "2", "0")),
        (b, p("1", "1", "0")),
        (ellipse(), p("5", "0", "0")),
        (line(), p("1", "2", "3")),
    ] {
        point_eq(&c.parameter_origin(&store), &expected.to_real());
    }
    // Nonsquare radius with a unit reference must remain exact, not a rational approximation.
    let point = circle("2").parameter_origin(&store);
    real_eq(&point.x, &real("2").sqrt());
    real_eq(&point.y, &real("0"));
    real_eq(&point.z, &real("0"));
});

blind_test!(case_43_parameters_seams_poles_and_off_surface, {
    let bits = 160;
    let cylinder = cylinder("5");
    for (point, u, v) in [
        (
            p("3", "4", "7"),
            Ball::atan2(&ball("4", 256), &ball("3", 256)).unwrap(),
            r("7"),
        ),
        (p("-5", "0", "0"), Ball::pi(256), r("0")),
        (
            p("0", "-5", "0"),
            &Ball::pi(256) * &ball("3/2", 256),
            r("0"),
        ),
        (
            p("6", "8", "-2"),
            Ball::atan2(&ball("4", 256), &ball("3", 256)).unwrap(),
            r("-2"),
        ),
    ] {
        let uv = cylinder
            .parameters_of(&balls(&point, bits), &catallax_numbers::Decider::default())
            .unwrap()
            .expect("zero-radius box away from seams");
        overlap(&uv[0], &u);
        assert!(uv[1].contains(&v));
    }
    if let Some(uv) = cylinder
        .parameters_of(
            &balls(&p("5", "0", "0"), bits),
            &catallax_numbers::Decider::default(),
        )
        .unwrap()
    {
        assert!(uv[0].contains(&r("0")));
        assert!(uv[1].contains(&r("0")));
    }
    let uv = torus()
        .parameters_of(
            &balls(&p("0", "2", "0"), bits),
            &catallax_numbers::Decider::default(),
        )
        .unwrap()
        .expect("inner equator is not the v seam");
    overlap(&uv[0], &(&Ball::pi(256) * &ball("1/2", 256)));
    overlap(&uv[1], &Ball::pi(256));
    let s = sphere("3");
    for pole in [p("0", "0", "3"), p("0", "0", "-3")] {
        // u is undefined here, so do not prescribe u if Some is returned.
        if let Some(uv) = s
            .parameters_of(&balls(&pole, bits), &catallax_numbers::Decider::default())
            .unwrap()
        {
            let sign = if pole.z > r("0") { "1/2" } else { "-1/2" };
            overlap(&uv[1], &(&Ball::pi(256) * &ball(sign, 256)));
        }
    }
    let plane = Surface::Plane(
        Plane::new(
            p("1", "1", "1").into(),
            v("0", "0", "1").into(),
            v("0", "2", "0").into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    );
    let box_ = [
        ball("-3", bits).hull(&ball("2", bits)),
        ball("-4", bits).hull(&ball("6", bits)),
        ball("-100", bits).hull(&ball("100", bits)),
    ];
    let uv = plane
        .parameters_of(&box_, &catallax_numbers::Decider::default())
        .unwrap()
        .expect("plane parameters are total");
    // u=y-1 in [-5,5], v=1-x in [-1,4], irrespective of z.
    for x in ["-5", "5"] {
        assert!(uv[0].contains(&r(x)));
    }
    for x in ["-1", "4"] {
        assert!(uv[1].contains(&r(x)));
    }
    // Closed boxes crossing a seam may return None, but Some must cover both sides.
    let crossing = [
        ball("5", bits),
        ball("-1", bits).hull(&ball("1", bits)),
        ball("0", bits),
    ];
    if let Some(uv) = cylinder
        .parameters_of(&crossing, &catallax_numbers::Decider::default())
        .unwrap()
    {
        let a = Ball::atan2(&ball("1", 256), &ball("5", 256)).unwrap();
        let b = &(&Ball::pi(256) * &ball("2", 256)) - &a;
        overlap(&uv[0], &a);
        overlap(&uv[0], &b);
        assert!(uv[0].contains(&r("0")));
    }
    for point in [
        p("4", "0", "0"),
        p("0", "4", "0"),
        p("3", "0", "0"),
        p("0", "0", "0"),
    ] {
        if let Some(uv) = torus()
            .parameters_of(&balls(&point, bits), &catallax_numbers::Decider::default())
            .unwrap()
        {
            if point.y == r("4") {
                overlap(&uv[0], &(&Ball::pi(256) * &ball("1/2", 256)));
                assert!(uv[1].contains(&r("0")));
            }
            if point.x == r("4") {
                assert!(uv[0].contains(&r("0")));
                assert!(uv[1].contains(&r("0")));
            }
        }
    }
});

blind_test!(wide_input_boxes_enclose_values_and_derivatives, {
    let bits = 160;
    let u = ball("0", bits).hull(&Ball::pi(bits));
    let v = ball("-2", bits).hull(&ball("3", bits));
    let cylinder = cylinder("5");
    let c = cylinder
        .evaluate(&u, &v, &catallax_numbers::Decider::default())
        .unwrap();
    let d = cylinder
        .derivatives(&u, &v, &catallax_numbers::Decider::default())
        .unwrap();
    for point in [p("5", "0", "-2"), p("-5", "0", "3"), p("0", "5", "1")] {
        contains(&c, &point);
    }
    for du in [p("0", "5", "0"), p("-5", "0", "0"), p("0", "-5", "0")] {
        contains(&d[0], &du);
    }
    contains(&d[1], &p("0", "0", "1"));
    let torus = torus();
    let ring = torus
        .evaluate(&u, &u, &catallax_numbers::Decider::default())
        .unwrap();
    let d = torus
        .derivatives(&u, &u, &catallax_numbers::Decider::default())
        .unwrap();
    for point in [
        p("4", "0", "0"),
        p("-4", "0", "0"),
        p("2", "0", "0"),
        p("0", "3", "1"),
    ] {
        contains(&ring, &point);
    }
    for du in [p("0", "4", "0"), p("-3", "0", "0"), p("0", "-2", "0")] {
        contains(&d[0], &du);
    }
    for dv in [p("0", "0", "1"), p("-1", "0", "0"), p("0", "0", "-1")] {
        contains(&d[1], &dv);
    }
    let sphere = sphere("3");
    let b = sphere
        .evaluate(&u, &u, &catallax_numbers::Decider::default())
        .unwrap();
    for point in [p("3", "0", "0"), p("0", "0", "3"), p("-3", "0", "0")] {
        contains(&b, &point);
    }
    let cone = cone();
    let b = cone
        .evaluate(&u, &v, &catallax_numbers::Decider::default())
        .unwrap();
    for point in [p("-1", "0", "-2"), p("0", "0", "0"), p("0", "3/2", "3")] {
        contains(&b, &point);
    }
    let circle = circle("25")
        .evaluate(&u, &catallax_numbers::Decider::default())
        .unwrap()
        .unwrap();
    for point in [p("5", "0", "0"), p("0", "5", "0"), p("-5", "0", "0")] {
        contains(&circle, &point);
    }
    let ellipse = ellipse()
        .evaluate(&u, &catallax_numbers::Decider::default())
        .unwrap()
        .unwrap();
    for point in [p("5", "0", "0"), p("0", "3", "0"), p("-5", "0", "0")] {
        contains(&ellipse, &point);
    }
});

blind_test!(case_43_generated_roundtrips, {
    let mut rng = Rng(0x243);
    for _ in 0..128 {
        let u = Ball::from_rational(
            &Rational::from_ratio((rng.next() % 500 + 1) as i64, 100).unwrap(),
            160,
        );
        let v = Ball::from_rational(
            &Rational::from_ratio((rng.next() % 100 + 1) as i64, 100).unwrap(),
            160,
        );
        for surface in surfaces() {
            let point = surface
                .evaluate(&u, &v, &catallax_numbers::Decider::default())
                .unwrap();
            if let Some(uv) = surface
                .parameters_of(&point, &catallax_numbers::Decider::default())
                .unwrap()
            {
                // Inputs are enclosures of rationals, so compare their exact rational midpoints.
                assert!(uv[0].contains(&u.mid()));
                assert!(uv[1].contains(&v.mid()));
            }
        }
    }
});
