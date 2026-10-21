#![allow(missing_docs)]
mod common;
use catallax_shapes::*;
use common::*;

blind_test!(case_11_constructor_refusals_and_precedence, {
    let o = p("0", "0", "0");
    let z = v("0", "0", "1");
    let x = v("1", "0", "0");
    let zero = v("0", "0", "0");
    let bad = v("1", "0", "1");
    for (a, b) in [
        (zero.clone(), x.clone()),
        (z.clone(), zero.clone()),
        (zero.clone(), zero.clone()),
    ] {
        rejected(
            Plane::new(o.clone(), a.clone(), b.clone()),
            "shapes.zero-vector",
        );
        rejected(
            Cylinder::new(o.clone(), a.clone(), r("-1"), b.clone()),
            "shapes.zero-vector",
        );
        rejected(
            Cone::new(o.clone(), a.clone(), r("-1"), b.clone()),
            "shapes.zero-vector",
        );
        rejected(
            Sphere::new(o.clone(), r("-1"), a.clone(), b.clone()),
            "shapes.zero-vector",
        );
        rejected(
            Torus::new(o.clone(), a.clone(), r("0"), r("1"), b.clone()),
            "shapes.zero-vector",
        );
        rejected(
            Circle::new(o.clone(), a.clone(), r("-1"), b.clone()),
            "shapes.zero-vector",
        );
        rejected(
            Ellipse::new(o.clone(), a, b, r("1"), r("1")),
            "shapes.zero-vector",
        );
    }
    rejected(Line::new(o.clone(), zero), "shapes.zero-vector");
    for radius in ["0", "-1"] {
        rejected(
            Cylinder::new(o.clone(), z.clone(), r(radius), bad.clone()),
            "shapes.non-positive-radius",
        );
        rejected(
            Sphere::new(o.clone(), r(radius), z.clone(), bad.clone()),
            "shapes.non-positive-radius",
        );
        rejected(
            Circle::new(o.clone(), z.clone(), r(radius), bad.clone()),
            "shapes.non-positive-radius",
        );
        rejected(
            Cone::new(o.clone(), z.clone(), r(radius), bad.clone()),
            "shapes.non-positive-slope",
        );
        rejected(
            Torus::new(o.clone(), z.clone(), r(radius), r("1"), bad.clone()),
            "shapes.non-positive-radius",
        );
        rejected(
            Torus::new(o.clone(), z.clone(), r("3"), r(radius), bad.clone()),
            "shapes.non-positive-radius",
        );
    }
    for minor in ["3", "4"] {
        rejected(
            Torus::new(o.clone(), z.clone(), r("3"), r(minor), bad.clone()),
            "shapes.torus-self-intersecting",
        );
    }
    for (major, minor) in [
        ("1", "1"),
        ("1", "2"),
        ("1", "0"),
        ("1", "-1"),
        ("0", "1"),
        ("-1", "-2"),
    ] {
        rejected(
            Ellipse::new(o.clone(), z.clone(), bad.clone(), r(major), r(minor)),
            "shapes.ellipse-axes",
        );
    }
    rejected(
        Plane::new(o.clone(), z.clone(), bad.clone()),
        "shapes.not-perpendicular",
    );
    rejected(
        Cylinder::new(o.clone(), z.clone(), r("1"), bad.clone()),
        "shapes.not-perpendicular",
    );
    rejected(
        Cone::new(o.clone(), z.clone(), r("1"), bad.clone()),
        "shapes.not-perpendicular",
    );
    rejected(
        Sphere::new(o.clone(), r("1"), z.clone(), bad.clone()),
        "shapes.not-perpendicular",
    );
    rejected(
        Torus::new(o.clone(), z.clone(), r("3"), r("1"), bad.clone()),
        "shapes.not-perpendicular",
    );
    rejected(
        Circle::new(o.clone(), z.clone(), r("1"), bad.clone()),
        "shapes.not-perpendicular",
    );
    rejected(
        Ellipse::new(o.clone(), z.clone(), bad, r("2"), r("1")),
        "shapes.not-perpendicular",
    );
    // Exact validation, with no epsilon around degeneracies.
    let tiny = r("1e-100");
    let almost = v("1", "0", "1e-100");
    rejected(
        Plane::new(o.clone(), z.clone(), almost),
        "shapes.not-perpendicular",
    );
    assert_eq!(
        Cylinder::new(o.clone(), z.clone(), tiny.clone(), x.clone())
            .unwrap()
            .radius,
        tiny
    );
    assert_eq!(
        Cone::new(o.clone(), z.clone(), r("1e-100"), x.clone())
            .unwrap()
            .slope,
        r("1e-100")
    );
    assert_eq!(
        Torus::new(
            o.clone(),
            z.clone(),
            r("1"),
            r("1") - r("1e-100"),
            x.clone()
        )
        .unwrap()
        .minor,
        r("1") - r("1e-100")
    );
    assert_eq!(
        Ellipse::new(o, z, x, r("1") + r("1e-100"), r("1"))
            .unwrap()
            .minor_sq,
        r("1")
    );
});

blind_test!(case_12_expanded_polynomials_without_rescaling, {
    let mut c = cylinder("3");
    if let Surface::Cylinder(c) = &mut c {
        c.axis = v("0", "0", "2");
    }
    let s = Surface::Sphere(
        Sphere::new(p("1", "2", "3"), r("2"), v("0", "0", "1"), v("1", "0", "0")).unwrap(),
    );
    let cases = [
        (plane(), poly(&[([0, 0, 1], "2"), ([0, 0, 0], "-2")])),
        (
            cylinder("3"),
            poly(&[([2, 0, 0], "1"), ([0, 2, 0], "1"), ([0, 0, 0], "-9")]),
        ),
        (
            c,
            poly(&[([2, 0, 0], "4"), ([0, 2, 0], "4"), ([0, 0, 0], "-36")]),
        ),
        (
            cone(),
            poly(&[([2, 0, 0], "1"), ([0, 2, 0], "1"), ([0, 0, 2], "-1/4")]),
        ),
        (
            s,
            poly(&[
                ([2, 0, 0], "1"),
                ([0, 2, 0], "1"),
                ([0, 0, 2], "1"),
                ([1, 0, 0], "-2"),
                ([0, 1, 0], "-4"),
                ([0, 0, 1], "-6"),
                ([0, 0, 0], "10"),
            ]),
        ),
        // Expanding (x²+y²+z²+8)²-36(x²+y²) fixes every coefficient.
        (
            torus(),
            poly(&[
                ([4, 0, 0], "1"),
                ([0, 4, 0], "1"),
                ([0, 0, 4], "1"),
                ([2, 2, 0], "2"),
                ([2, 0, 2], "2"),
                ([0, 2, 2], "2"),
                ([2, 0, 0], "-20"),
                ([0, 2, 0], "-20"),
                ([0, 0, 2], "16"),
                ([0, 0, 0], "64"),
            ]),
        ),
    ];
    for (s, expected) in cases {
        assert_eq!(s.implicit(), expected);
    }
    for (point, expected) in [
        (p("4", "0", "0"), "0"),
        (p("2", "0", "0"), "0"),
        (p("3", "0", "1"), "0"),
        (p("0", "0", "0"), "64"),
        (p("3", "0", "0"), "-35"),
    ] {
        assert_eq!(torus().implicit().eval(&xyz(&point)), r(expected));
    }
});

blind_test!(cases_13_14_15_transform_and_gradients, {
    let pose = Pose::new(rotation(["2", "1", "0", "0"]), v("1000000", "-7", "1/3"));
    let expected = Surface::Cylinder(
        Cylinder::new(
            p("1000000", "-7", "1/3"),
            v("0", "-4/5", "3/5"),
            r("5"),
            v("1", "0", "0"),
        )
        .unwrap(),
    );
    assert_eq!(cylinder("5").transformed(&pose), expected);
    assert_eq!(
        pose.apply_point(&p("3", "4", "7")),
        p("1000003", "-51/5", "116/15")
    );
    assert_eq!(
        expected
            .implicit()
            .eval(&xyz(&p("1000003", "-51/5", "116/15"))),
        r("0")
    );
    vector_eq(
        &sphere("3").gradient(&p("1", "2", "2").to_real()),
        &v("2", "4", "4"),
    );
    vector_eq(
        &plane().gradient(&p("7", "-8", "9").to_real()),
        &v("0", "0", "2"),
    );
    vector_eq(
        &cone().gradient(&p("0", "0", "0").to_real()),
        &v("0", "0", "0"),
    );
    vector_eq(
        &cylinder("5").gradient(&p("3", "4", "7").to_real()),
        &v("6", "8", "0"),
    );
    // f=(q+8)^2-36(x²+y²); grad=(4x(q-10),4y(q-10),4z(q+8)).
    vector_eq(
        &torus().gradient(&p("3", "0", "1").to_real()),
        &v("0", "0", "72"),
    );
    let mut rng = Rng(0x213);
    for _ in 0..200 {
        let pose = rng.pose();
        let point = rp(std::array::from_fn(|_| rng.rational()));
        for surface in surfaces() {
            let expected = implicit_oracle(&surface, &point);
            let moved = surface.transformed(&pose);
            assert_eq!(oq(&surface.implicit().eval(&xyz(&point))), expected);
            assert_eq!(
                oq(&moved.implicit().eval(&xyz(&pose.apply_point(&point)))),
                expected
            );
            assert_eq!(moved.transformed(&pose.inverse()), surface);
        }
        for curve in curves() {
            let moved = curve.transformed(&pose).unwrap();
            assert_eq!(moved.transformed(&pose.inverse()), Some(curve));
        }
    }
});

blind_test!(cases_24_25_26_equation_order, {
    let store = GeomStore::new();
    assert_eq!(
        line().equations(&store),
        vec![
            poly(&[([0, 0, 1], "1"), ([0, 0, 0], "-3")]),
            poly(&[([0, 0, 1], "2"), ([0, 0, 0], "-6")]),
            poly(&[([1, 0, 0], "-1"), ([0, 1, 0], "-2"), ([0, 0, 0], "5")])
        ]
    );
    let axis = Curve::Line(Line::new(p("0", "0", "0"), v("1", "0", "0")).unwrap());
    assert_eq!(
        axis.equations(&store),
        vec![poly(&[([0, 0, 1], "1")]), poly(&[([0, 1, 0], "-1")])]
    );
    assert_eq!(
        circle("25").equations(&store),
        vec![
            poly(&[([0, 0, 1], "1")]),
            poly(&[
                ([2, 0, 0], "1"),
                ([0, 2, 0], "1"),
                ([0, 0, 2], "1"),
                ([0, 0, 0], "-25")
            ])
        ]
    );
    assert_eq!(
        ellipse().equations(&store),
        vec![
            poly(&[([0, 0, 1], "1")]),
            poly(&[([2, 0, 0], "1/25"), ([0, 2, 0], "1/9"), ([0, 0, 0], "-1")])
        ]
    );
    assert_eq!(line().is_closed(), Some(false));
    assert_eq!(circle("2").is_closed(), Some(true));
    assert_eq!(ellipse().is_closed(), Some(true));
    assert_eq!(line().kind(), CurveKind::Line);
    assert_eq!(circle("2").kind(), CurveKind::Circle);
    assert_eq!(ellipse().kind(), CurveKind::Ellipse);
});

blind_test!(transforms_move_every_field_without_normalising, {
    let pose = Pose::new(rotation(["1", "0", "0", "1"]), v("10", "20", "30"));
    let expected = [
        Surface::Plane(
            Plane::new(p("10", "20", "31"), v("0", "0", "2"), v("0", "1", "0")).unwrap(),
        ),
        Surface::Cylinder(
            Cylinder::new(
                p("10", "20", "30"),
                v("0", "0", "1"),
                r("5"),
                v("0", "1", "0"),
            )
            .unwrap(),
        ),
        Surface::Cone(
            Cone::new(
                p("10", "20", "30"),
                v("0", "0", "1"),
                r("1/2"),
                v("0", "1", "0"),
            )
            .unwrap(),
        ),
        Surface::Sphere(
            Sphere::new(
                p("10", "20", "30"),
                r("3"),
                v("0", "0", "1"),
                v("0", "1", "0"),
            )
            .unwrap(),
        ),
        Surface::Torus(
            Torus::new(
                p("10", "20", "30"),
                v("0", "0", "1"),
                r("3"),
                r("1"),
                v("0", "1", "0"),
            )
            .unwrap(),
        ),
    ];
    for (surface, expected) in surfaces().into_iter().zip(expected) {
        assert_eq!(surface.transformed(&pose), expected);
    }
    let expected = [
        Curve::Line(Line::new(p("8", "21", "33"), v("1", "2", "0")).unwrap()),
        Curve::Circle(
            Circle::new(
                p("10", "20", "30"),
                v("0", "0", "1"),
                r("25"),
                v("0", "1", "0"),
            )
            .unwrap(),
        ),
        Curve::Ellipse(
            Ellipse::new(
                p("10", "20", "30"),
                v("0", "0", "1"),
                v("0", "2", "0"),
                r("25"),
                r("9"),
            )
            .unwrap(),
        ),
    ];
    for (curve, expected) in curves().into_iter().zip(expected) {
        assert_eq!(curve.transformed(&pose), Some(expected));
    }
});
