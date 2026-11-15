//! STEP revolution readers preserve the swept carrier without rounding its axis.
use catallax_import::{ImportRefusal, Surface, read_step};
use catallax_numbers::{Ball, Decider, Rational};

const SQUARE: &str = include_str!("fixtures/surface-square.step");
fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn source(axis: &str, origin: &str, points: [&str; 2], rational: bool) -> String {
    let profile = if rational {
        "#42=(BOUNDED_CURVE() B_SPLINE_CURVE(1,(#40,#41),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS((2,2),(0.,1.),.UNSPECIFIED.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE((1.,2.)) REPRESENTATION_ITEM(''));"
    } else {
        "#42=B_SPLINE_CURVE_WITH_KNOTS('',1,(#40,#41),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);"
    };
    let extra = format!(
        "#40=CARTESIAN_POINT('',({}));\n#41=CARTESIAN_POINT('',({}));\n{profile}\n#43=AXIS1_PLACEMENT('',#1,#7);\n",
        points[0], points[1]
    );
    SQUARE
        .replace("#9=PLANE('',#8);", "#9=SURFACE_OF_REVOLUTION('',#42,#43);")
        .replace(
            "#7=DIRECTION('',(0.,0.,1.));",
            &format!("#7=DIRECTION('',({axis}));"),
        )
        .replace(
            "#1=CARTESIAN_POINT('',(0.,0.,0.));",
            &format!("#1=CARTESIAN_POINT('',({origin}));"),
        )
        .replace(
            "ENDSEC;\nEND-ISO-10303-21;",
            &format!("{extra}ENDSEC;\nEND-ISO-10303-21;"),
        )
}
fn geometry(source: &str) -> catallax_import::StepSplineSurface {
    let file = read_step(source.as_bytes()).unwrap();
    let Surface::BSpline(surface) = &file.surface_models[0].shells[0].faces[0].surface else {
        panic!("spline revolution")
    };
    surface.clone()
}
fn point(surface: &catallax_import::StepSplineSurface, u: &Rational, v: &Rational) -> [Ball; 3] {
    surface
        .geometry
        .to_surface(&mut Decider::default())
        .unwrap()
        .jet_box([u, u], [v, v], 128, &Decider::default())
        .unwrap()
        .value
}
fn contains(value: &Ball, expected: &Rational) {
    assert!(value.lo() <= *expected && value.hi() >= *expected);
    assert!(value.hi() - value.lo() < q("0.000000000000000000000001"));
}

#[test]
fn a_polynomial_meridian_sweeps_a_complete_cylinder_with_the_source_orientation() {
    for (axis, sign) in [("0.,0.,1.", 1), ("0.,0.,-1.", -1)] {
        let surface = geometry(&source(axis, "0.,0.,0.", ["2.,0.,0.", "2.,0.,3."], false));
        assert_eq!(surface.geometry.degrees(), [2, 1]);
        assert_eq!(surface.geometry.control_count(), [9, 2]);
        assert_eq!(surface.closed, [Some(true), Some(false)]);
        assert_eq!(surface.self_intersect, None);
        for (u, want) in [
            (0, [q("2"), q("0"), q("1.5")]),
            (1, [q("0"), Rational::from_i64(2 * sign), q("1.5")]),
            (2, [q("-2"), q("0"), q("1.5")]),
            (3, [q("0"), Rational::from_i64(-2 * sign), q("1.5")]),
            (4, [q("2"), q("0"), q("1.5")]),
        ] {
            for (actual, expected) in point(&surface, &Rational::from_i64(u), &q("0.5"))
                .iter()
                .zip(want)
            {
                contains(actual, &expected);
            }
        }
    }
}

#[test]
fn a_rational_meridian_and_irrational_axis_length_keep_exact_circular_sections() {
    let surface = geometry(&source(
        "1.,1.,0.",
        "5.,-2.,3.",
        ["7.,-4.,3.", "10.,-1.,3."],
        true,
    ));
    let center = [q("5"), q("-2"), q("3")];
    let half = Ball::from_rational(&q("0.5"), 128);
    for u in [
        "0", "0.125", "0.5", "1", "1.25", "2", "2.5", "3", "3.875", "4",
    ] {
        let value = point(&surface, &q(u), &q("0.5"));
        let delta: [Ball; 3] =
            std::array::from_fn(|i| &value[i] - &Ball::from_rational(&center[i], 128));
        let axial = &delta[0] + &delta[1];
        contains(&axial, &q("4"));
        let transverse = &delta[0] - &delta[1];
        let radial = &(&transverse.square() * &half) + &delta[2].square();
        contains(&radial, &q("8"));
    }
}

#[test]
fn omitted_axis_defaults_to_z_and_length_units_do_not_scale_weights_or_knots() {
    let original = source("0.,0.,1.", "0.,0.,0.", ["2.,0.,0.", "2.,0.,3."], true);
    let surface = geometry(&original);
    assert_eq!(
        surface,
        geometry(&original.replace(
            "#43=AXIS1_PLACEMENT('',#1,#7)",
            "#43=AXIS1_PLACEMENT('',#1,$)"
        ))
    );
    let metres = geometry(&original.replace("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT($,.METRE.)"));
    assert_eq!(surface.geometry.weights(), metres.geometry.weights());
    assert_eq!(surface.geometry.knots(), metres.geometry.knots());
    for (before, after) in surface
        .geometry
        .controls()
        .iter()
        .zip(metres.geometry.controls())
    {
        for i in 0..3 {
            assert_eq!(&before[i] * &q("1000"), after[i]);
        }
    }
}

#[test]
fn a_zero_axis_and_wrong_axis_placement_arity_refuse() {
    let zero = source("0.,0.,0.", "0.,0.,0.", ["2.,0.,0.", "2.,0.,3."], false);
    assert!(matches!(
        read_step(zero.as_bytes()),
        Err(ImportRefusal::Malformed { .. })
    ));
    let bad = source("0.,0.,1.", "0.,0.,0.", ["2.,0.,0.", "2.,0.,3."], false).replace(
        "#43=AXIS1_PLACEMENT('',#1,#7)",
        "#43=AXIS1_PLACEMENT('',#1,#7,$)",
    );
    assert!(matches!(
        read_step(bad.as_bytes()),
        Err(ImportRefusal::Malformed { .. })
    ));
}

#[test]
fn the_closed_seam_has_matching_position_and_first_derivatives() {
    let surface = geometry(&source(
        "0.,0.,1.",
        "0.,0.,0.",
        ["2.,0.,0.", "2.,0.,3."],
        false,
    ));
    let evaluator = surface
        .geometry
        .to_surface(&mut Decider::default())
        .unwrap();
    let v = q("0.5");
    for u in [q("0"), q("4")] {
        let jet = evaluator
            .jet_box([&u, &u], [&v, &v], 128, &Decider::default())
            .unwrap();
        for (actual, expected) in jet.first[0].iter().zip([q("0"), q("4"), q("0")]) {
            contains(actual, &expected);
        }
        for (actual, expected) in jet.first[1].iter().zip([q("0"), q("0"), q("3")]) {
            contains(actual, &expected);
        }
    }
}

#[test]
fn every_quadratic_span_satisfies_the_exact_circle_polynomial() {
    let surface = geometry(&source(
        "1.,1.,0.",
        "5.,-2.,3.",
        ["7.,-4.,3.", "10.,-1.,3."],
        true,
    ));
    let center = [q("5"), q("-2"), q("3")];
    let two = q("2");
    let sixteen = q("16");
    for start in [0, 2, 4, 6] {
        for column in 0..2 {
            let h: [[Rational; 4]; 3] = std::array::from_fn(|row| {
                let index = (start + row) * 2 + column;
                let weight = &surface.geometry.weights()[index];
                let point = &surface.geometry.controls()[index];
                std::array::from_fn(|i| {
                    if i == 3 {
                        weight.clone()
                    } else {
                        &(&point[i] - &center[i]) * weight
                    }
                })
            });
            let coefficients: [[Rational; 4]; 3] = std::array::from_fn(|power| {
                std::array::from_fn(|i| match power {
                    0 => h[0][i].clone(),
                    1 => &two * &(&h[1][i] - &h[0][i]),
                    2 => &(&h[0][i] - &(&two * &h[1][i])) + &h[2][i],
                    _ => unreachable!(),
                })
            });
            let axial = Rational::from_i64(if column == 0 { 0 } else { 6 });
            for coefficient in &coefficients {
                assert_eq!(&coefficient[0] + &coefficient[1], &axial * &coefficient[3]);
            }
            let mut equation: [Rational; 5] = std::array::from_fn(|_| Rational::zero());
            for a in 0..3 {
                for b in 0..3 {
                    let x = &coefficients[a][0] - &coefficients[a][1];
                    let y = &coefficients[b][0] - &coefficients[b][1];
                    let xy = &x * &y;
                    let zz = &two * &(&coefficients[a][2] * &coefficients[b][2]);
                    let ww = &sixteen * &(&coefficients[a][3] * &coefficients[b][3]);
                    let sum = &(&equation[a + b] + &xy) + &zz;
                    equation[a + b] = &sum - &ww;
                }
            }
            assert!(equation.iter().all(Rational::is_zero));
        }
    }
}

#[test]
fn a_curved_rational_meridian_keeps_its_weights_coupled_to_rotation() {
    let original = source("0.,0.,1.", "0.,0.,0.", ["2.,0.,0.", "2.,0.,3."], true);
    let source=original.replace("B_SPLINE_CURVE(1,(#40,#41)","B_SPLINE_CURVE(3,(#40,#44,#45,#41)")
        .replace("B_SPLINE_CURVE_WITH_KNOTS((2,2)","B_SPLINE_CURVE_WITH_KNOTS((4,4)")
        .replace("RATIONAL_B_SPLINE_CURVE((1.,2.))","RATIONAL_B_SPLINE_CURVE((1.,3.,2.,4.))")
        .replace("#43=AXIS1_PLACEMENT", "#44=CARTESIAN_POINT('',(4.,0.,1.));\n#45=CARTESIAN_POINT('',(3.,0.,2.));\n#43=AXIS1_PLACEMENT");
    let surface = geometry(&source);
    assert_eq!(surface.geometry.degrees(), [2, 3]);
    assert_eq!(surface.geometry.control_count(), [9, 4]);
    for (u, expected) in [
        ("0", [q("3.2"), q("0"), q("1.65")]),
        ("0.5", [q("1.92"), q("2.56"), q("1.65")]),
        ("1", [q("0"), q("3.2"), q("1.65")]),
    ] {
        for (actual, expected) in point(&surface, &q(u), &q("0.5")).iter().zip(expected) {
            contains(actual, &expected);
        }
    }
}

#[test]
fn proportional_axis_ratios_produce_the_same_bounded_control_net() {
    let make = |axis| {
        geometry(&source(
            axis,
            "5.,-2.,3.",
            ["7.,-4.,3.", "10.,-1.,3."],
            true,
        ))
    };
    let expected = make("1.,1.,0.");
    for axis in [
        "1000000000000.,1000000000000.,0.",
        "0.000000000001,0.000000000001,0.",
    ] {
        assert_eq!(make(axis), expected);
    }
}
