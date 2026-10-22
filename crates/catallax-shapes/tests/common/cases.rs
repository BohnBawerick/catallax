use super::*;
use catallax_numbers::*;
use catallax_shapes::*;

pub fn side_cases() -> Vec<(Surface, RPoint3, i8)> {
    let mut cases = vec![
        (plane(), p("5", "7", "1"), 0),
        (cylinder("5"), p("3", "4", "7"), 0),
        (cone(), p("3", "4", "10"), 0),
        (cone(), p("0", "0", "0"), 0),
        (cone(), p("3", "4", "-10"), 0),
        (sphere("3"), p("1", "2", "2"), 0),
        (torus(), p("54/25", "72/25", "4/5"), 0),
        (torus(), p("4", "0", "0"), 0),
        (torus(), p("2", "0", "0"), 0),
        (torus(), p("3", "0", "1"), 0),
    ];
    let pose = Pose::new(
        rotation(["2", "1", "0", "0"]),
        (v("1000000", "-7", "1/3")).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap();
    cases.push((
        cylinder("5").transformed(&pose),
        p("1000003", "-51/5", "116/15"),
        0,
    ));
    for e in [r("1e-30"), r("1e-7"), r("1e-3")] {
        for sign in [-1, 1] {
            let t = r("1") + &e * &Rational::from_i64(sign);
            let opposite = r("1") - &e * &Rational::from_i64(sign);
            cases.push((
                cylinder("5"),
                rp([&r("3") * &t, &r("4") * &t, r("7")]),
                sign as i8,
            ));
            cases.push((
                sphere("3"),
                rp([t.clone(), &r("2") * &t, &r("2") * &t]),
                sign as i8,
            ));
            cases.push((plane(), rp([r("5"), r("7"), t.clone()]), sign as i8));
            cases.push((cone(), rp([r("3"), r("4"), r("10") * opposite]), sign as i8));
            // rho=18/5, so D²=9/25+16t²/25. D²-1 has sign(t-1).
            cases.push((
                torus(),
                rp([r("54/25"), r("72/25"), r("4/5") * t]),
                sign as i8,
            ));
        }
    }
    cases
}
pub fn curve_cases() -> Vec<(Curve, RPoint3, bool)> {
    let e = r("1e-30");
    vec![
        (line(), p("5", "0", "3"), true),
        (line(), rp([r("5"), r("0"), &r("3") + &e]), false),
        (line(), rp([&r("5") + &e, r("0"), r("3")]), false),
        (circle("25"), p("3", "4", "0"), true),
        (circle("25"), rp([r("3"), r("4"), e.clone()]), false),
        (circle("25"), rp([r("3"), &r("4") + &e, r("0")]), false),
        (circle("2"), p("1", "1", "0"), true),
        (circle(&(&r("2") + &e).to_string()), p("1", "1", "0"), false),
        (ellipse(), p("3", "12/5", "0"), true),
        (ellipse(), p("5", "0", "0"), true),
        (ellipse(), p("0", "3", "0"), true),
        (ellipse(), rp([r("3"), &r("12/5") + &e, r("0")]), false),
    ]
}
pub fn surface_pairs() -> Vec<(Surface, Surface, bool)> {
    let a = Surface::Plane(
        Plane::new(
            (p("0", "0", "1")).into(),
            (v("0", "0", "1")).into(),
            (v("1", "0", "0")).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    );
    let b = Surface::Plane(
        Plane::new(
            (p("5", "7", "1")).into(),
            (v("0", "0", "-4")).into(),
            (v("0", "1", "0")).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    );
    let mut near = b.clone();
    if let Surface::Plane(b) = &mut near {
        b.origin.z = &b.origin.z + &real("1e-30");
    }
    let mut pairs = vec![(a.clone(), b, true), (a, near, false)];
    let a = cylinder("2");
    let b = Surface::Cylinder(
        Cylinder::new(
            (p("0", "0", "5")).into(),
            (v("0", "0", "-3")).into(),
            (r("2")).into(),
            (v("0", "1", "0")).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    );
    pairs.push((a.clone(), b, true));
    for b in [
        Cylinder::new(
            (p("1e-30", "0", "0")).into(),
            (v("0", "0", "1")).into(),
            (r("2")).into(),
            (v("1", "0", "0")).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
        Cylinder::new(
            (p("0", "0", "0")).into(),
            (v("0", "0", "1")).into(),
            (r("2") + r("1e-30")).into(),
            (v("1", "0", "0")).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
        Cylinder::new(
            (p("0", "0", "0")).into(),
            (v("1e-30", "0", "1")).into(),
            (r("2")).into(),
            (v("1", "0", "-1e-30")).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ] {
        pairs.push((a.clone(), Surface::Cylinder(b), false));
    }
    for (axis, slope, same) in [
        ("2", r("1/2"), true),
        ("-1", r("1/2"), false),
        ("1", r("1/2") + r("1e-30"), false),
    ] {
        pairs.push((
            cone(),
            Surface::Cone(
                Cone::new(
                    (p("0", "0", "0")).into(),
                    (v("0", "0", axis)).into(),
                    (slope).into(),
                    (v("0", "1", "0")).into(),
                    &mut catallax_numbers::Decider::default(),
                )
                .unwrap(),
            ),
            same,
        ));
    }
    pairs.push((
        sphere("3"),
        Surface::Sphere(
            Sphere::new(
                (p("0", "0", "0")).into(),
                (r("3")).into(),
                (v("1", "0", "0")).into(),
                (v("0", "0", "-2")).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
        true,
    ));
    pairs.push((
        sphere("3"),
        sphere("3.000000000000000000000000000001"),
        false,
    ));
    pairs.push((
        sphere("3"),
        Surface::Sphere(
            Sphere::new(
                (p("1e-30", "0", "0")).into(),
                (r("3")).into(),
                (v("0", "0", "1")).into(),
                (v("1", "0", "0")).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
        false,
    ));
    pairs.push((
        torus(),
        Surface::Torus(
            Torus::new(
                (p("0", "0", "0")).into(),
                (v("0", "0", "-7")).into(),
                (r("3")).into(),
                (r("1")).into(),
                (v("0", "1", "0")).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
        true,
    ));
    pairs.push((
        torus(),
        Surface::Torus(
            Torus::new(
                (p("0", "0", "0")).into(),
                (v("0", "0", "1")).into(),
                (r("3")).into(),
                (r("1") + r("1e-30")).into(),
                (v("1", "0", "0")).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
        false,
    ));
    pairs.push((cylinder("1"), sphere("1"), false));
    pairs
}
pub fn normal_cases() -> Vec<(Surface, Surface, RPoint3, i8)> {
    let plane_at = |normal| {
        Surface::Plane(
            Plane::new(
                (p("0", "0", "2")).into(),
                (v("0", "0", normal)).into(),
                (v("1", "0", "0")).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        )
    };
    vec![
        (sphere("2"), plane_at("1"), p("0", "0", "2"), 1),
        (sphere("2"), plane_at("-1"), p("0", "0", "2"), -1),
        (
            Surface::Plane(
                Plane::new(
                    p("0", "0", "0").into(),
                    v("1", "0", "0").into(),
                    v("0", "1", "0").into(),
                    &mut catallax_numbers::Decider::default(),
                )
                .unwrap(),
            ),
            Surface::Plane(
                Plane::new(
                    p("0", "0", "0").into(),
                    v("0", "1", "0").into(),
                    v("1", "0", "0").into(),
                    &mut catallax_numbers::Decider::default(),
                )
                .unwrap(),
            ),
            p("0", "0", "5"),
            0,
        ),
        (cone(), plane(), p("0", "0", "0"), 0),
    ]
}
