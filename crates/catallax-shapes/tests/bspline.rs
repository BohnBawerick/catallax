#![allow(missing_docs)]

use catallax_numbers::{Ball, Decider, Rational, Real, with_operation_budget};
use catallax_shapes::{BSplineCurve, BSplineSurface, Point3, SplineBoundary};

fn r(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn b(s: &str) -> Ball {
    Ball::from_rational(&r(s), 128)
}
fn p(x: &str, y: &str, z: &str) -> Point3 {
    Point3::new(r(x).into(), r(y).into(), r(z).into())
}
fn linear_knots() -> Vec<Rational> {
    ["0", "0", "1", "1"].map(r).to_vec()
}
fn bilinear(weights: Vec<Real>) -> BSplineSurface {
    BSplineSurface::new(
        [1, 1],
        [linear_knots(), linear_knots()],
        [2, 2],
        vec![
            p("0", "0", "0"),
            p("0", "1", "0"),
            p("1", "0", "0"),
            p("1", "1", "1"),
        ],
        weights,
        &mut Decider::default(),
    )
    .unwrap()
}

#[test]
fn bilinear_box_contains_values_and_all_chart_derivatives() {
    let surface = bilinear(vec![Real::one(); 4]);
    let u = b("1/4").hull(&b("3/4"));
    let v = b("1/4").hull(&b("1/2"));
    let jet = surface.jet(&u, &v, &Decider::default()).unwrap();
    for x in [r("1/4"), r("3/4")] {
        for y in [r("1/4"), r("1/2")] {
            assert!(jet.value[0].contains(&x));
            assert!(jet.value[1].contains(&y));
            assert!(jet.value[2].contains(&(&x * &y)));
            assert!(jet.first[0][0].contains(&r("1")));
            assert!(jet.first[0][1].contains(&r("0")));
            assert!(jet.first[0][2].contains(&y));
            assert!(jet.first[1][0].contains(&r("0")));
            assert!(jet.first[1][1].contains(&r("1")));
            assert!(jet.first[1][2].contains(&x));
        }
    }
    for (i, second) in jet.second.iter().enumerate() {
        for (axis, value) in second.iter().enumerate() {
            assert!(value.contains(&Rational::from_i64(i64::from(i == 1 && axis == 2))));
        }
    }
}

#[test]
fn rational_quotient_derivatives_include_the_weight_derivatives() {
    // x(u)=4u/(1+u), y=0, z=v; at u=1/2, x=4/3,
    // dx/du=16/9 and d²x/du²=-64/27.
    let surface = BSplineSurface::new(
        [1, 1],
        [linear_knots(), linear_knots()],
        [2, 2],
        vec![
            p("0", "0", "0"),
            p("0", "0", "1"),
            p("2", "0", "0"),
            p("2", "0", "1"),
        ],
        vec![
            Real::one(),
            Real::one(),
            Real::from_i64(2),
            Real::from_i64(2),
        ],
        &mut Decider::default(),
    )
    .unwrap();
    let jet = surface
        .jet(&b("1/2"), &b("1/4"), &Decider::default())
        .unwrap();
    assert!(jet.value[0].contains(&r("4/3")));
    assert!(jet.value[2].contains(&r("1/4")));
    assert!(jet.first[0][0].contains(&r("16/9")));
    assert!(jet.second[0][0].contains(&r("-64/27")));
    assert!(jet.first[1][2].contains(&r("1")));
    assert!(jet.second[1][2].contains(&r("0")));
    assert!(jet.second[2][2].contains(&r("0")));
}

#[test]
fn repeated_knot_spans_cover_a_cubic_and_both_endpoint_limits() {
    let u_knots = ["0", "0", "0", "0", "1/2", "1/2", "1/2", "1", "1", "1", "1"]
        .map(r)
        .to_vec();
    // Two Bezier spans represent the same global x=u³, z=v polynomial.
    let xs = ["0", "0", "0", "1/8", "1/4", "1/2", "1"];
    let controls = xs
        .into_iter()
        .flat_map(|x| [p(x, "0", "0"), p(x, "0", "1")])
        .collect();
    let surface = BSplineSurface::new(
        [3, 1],
        [u_knots, linear_knots()],
        [7, 2],
        controls,
        vec![Real::one(); 14],
        &mut Decider::default(),
    )
    .unwrap();
    for t in ["0", "1/4", "1/2", "3/4", "1"] {
        let t = r(t);
        let jet = surface
            .jet(
                &Ball::from_rational(&t, 128),
                &b("1/2"),
                &Decider::default(),
            )
            .unwrap();
        assert!(jet.value[0].contains(&(&(&t * &t) * &t)));
        assert!(jet.first[0][0].contains(&(&r("3") * &(&t * &t))));
        assert!(jet.second[0][0].contains(&(&r("6") * &t)));
    }
    let cover = surface
        .jet(&b("1/4").hull(&b("3/4")), &b("1/2"), &Decider::default())
        .unwrap();
    assert!(cover.value[0].contains(&r("1/64")));
    assert!(cover.value[0].contains(&r("27/64")));
    assert!(cover.first[0][0].contains(&r("3/16")));
    assert!(cover.first[0][0].contains(&r("27/16")));
}

#[test]
fn invalid_nets_domains_and_interrupted_construction_fail_closed() {
    let mut dec = Decider::default();
    for weights in [
        vec![Real::one(); 3],
        vec![Real::zero(); 4],
        vec![Real::from_i64(-1); 4],
    ] {
        assert!(
            BSplineSurface::new(
                [1, 1],
                [linear_knots(), linear_knots()],
                [2, 2],
                vec![p("0", "0", "0"); 4],
                weights,
                &mut dec
            )
            .is_err()
        );
    }
    for knots in [vec![r("0"); 4], ["0", "1", "0", "1"].map(r).to_vec()] {
        assert!(
            BSplineSurface::new(
                [1, 1],
                [knots, linear_knots()],
                [2, 2],
                vec![p("0", "0", "0"); 4],
                vec![Real::one(); 4],
                &mut dec
            )
            .is_err()
        );
    }
    let surface = bilinear(vec![Real::one(); 4]);
    assert!(surface.jet(&b("-1/4"), &b("1/2"), &dec).is_err());
    assert!(
        surface
            .jet(&b("3/4").hull(&b("5/4")), &b("1/2"), &dec)
            .is_err()
    );
    let (result, work) = with_operation_budget(0, || {
        BSplineSurface::new(
            [1, 1],
            [linear_knots(), linear_knots()],
            [2, 2],
            vec![p("0", "0", "0"); 4],
            vec![Real::one(); 4],
            &mut dec,
        )
    });
    assert!(result.is_err());
    assert!(work.exhausted);
}

#[test]
fn non_dyadic_finite_endpoints_remain_usable_without_clamping_external_boxes() {
    let u = [r("1/3"), r("2/3")];
    let v = [r("1/5"), r("2/5")];
    let surface = BSplineSurface::new(
        [1, 1],
        [
            vec![u[0].clone(), u[0].clone(), u[1].clone(), u[1].clone()],
            vec![v[0].clone(), v[0].clone(), v[1].clone(), v[1].clone()],
        ],
        [2, 2],
        vec![
            p("1/3", "0", "1/5"),
            p("1/3", "0", "2/5"),
            p("2/3", "0", "1/5"),
            p("2/3", "0", "2/5"),
        ],
        vec![Real::one(); 4],
        &mut Decider::default(),
    )
    .unwrap();
    let dec = Decider::default();
    for (x, z) in [(&u[0], &v[0]), (&u[1], &v[1])] {
        let jet = surface.jet_box([x, x], [z, z], 128, &dec).unwrap();
        assert!(jet.value[0].contains(x) && jet.value[2].contains(z));
        assert!(jet.first[0][0].contains(&Rational::one()));
        assert!(jet.first[1][2].contains(&Rational::one()));
        assert!(jet.value[0].hi() - jet.value[0].lo() < r("1e-20"));
    }
    let edge = BSplineCurve::boundary(&surface, SplineBoundary::VMin);
    for x in &u {
        let point = edge.evaluate_box([x, x], 128, &dec).unwrap();
        assert!(point[0].contains(x) && point[2].contains(&v[0]));
    }
    assert!(
        surface
            .jet(
                &Ball::from_rational(&u[0], 128),
                &Ball::from_rational(&v[0], 128),
                &dec
            )
            .is_err()
    );
    assert!(
        surface
            .jet_box([&r("0"), &u[1]], v.each_ref(), 128, &dec)
            .is_err()
    );
}

#[test]
fn shared_boundary_curves_keep_rational_position_and_parameter_derivatives() {
    let surface = BSplineSurface::new(
        [1, 1],
        [linear_knots(), linear_knots()],
        [2, 2],
        vec![
            p("0", "0", "0"),
            p("0", "0", "1"),
            p("2", "0", "0"),
            p("2", "0", "1"),
        ],
        vec![
            Real::one(),
            Real::one(),
            Real::from_i64(2),
            Real::from_i64(2),
        ],
        &mut Decider::default(),
    )
    .unwrap();
    let t = r("1/2");
    for (side, z) in [
        (SplineBoundary::VMin, r("0")),
        (SplineBoundary::VMax, r("1")),
    ] {
        let edge = BSplineCurve::boundary(&surface, side);
        let jet = edge.jet_box([&t, &t], 128, &Decider::default()).unwrap();
        assert!(jet.value[0].contains(&r("4/3")) && jet.value[2].contains(&z));
        assert!(jet.first[0].contains(&r("16/9")));
        assert!(jet.second[0].contains(&r("-64/27")));
        assert!(jet.first[2].contains(&r("0")));
    }
    let edge = BSplineCurve::boundary(&surface, SplineBoundary::UMax);
    let jet = edge.jet_box([&t, &t], 128, &Decider::default()).unwrap();
    assert!(jet.value[0].contains(&r("2")) && jet.value[2].contains(&t));
    assert!(jet.first[0].contains(&r("0")) && jet.first[2].contains(&r("1")));
}

#[test]
fn exact_curve_preimages_separate_incidence_from_nearby_points() {
    let mut dec = Decider::default();
    let curve = BSplineCurve::new(
        2,
        ["0", "0", "0", "1", "1", "1"].map(r).to_vec(),
        vec![p("0", "0", "0"), p("1/2", "0", "0"), p("1", "1", "0")],
        vec![Real::one(); 3],
        &mut dec,
    )
    .unwrap();
    let point = p("1/3", "1/9", "0");
    let parameter = curve.parameter_of(&point, &mut dec).unwrap().unwrap();
    assert!(
        dec.compare(parameter.value(), &r("1/3").into())
            .is_coincident()
    );
    let tangent = parameter.first(&dec).unwrap();
    for (component, expected) in [tangent.x, tangent.y, tangent.z]
        .iter()
        .zip(["1", "2/3", "0"])
    {
        assert!(dec.compare(component, &r(expected).into()).is_coincident());
    }
    assert!(
        curve
            .parameter_of(&p("1/3", "2/9", "0"), &mut dec)
            .unwrap()
            .is_none()
    );
    assert!(
        curve
            .parameter_of(&p("4/3", "16/9", "0"), &mut dec)
            .unwrap()
            .is_none()
    );
}
