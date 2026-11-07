use super::*;
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn p(x: &str, y: &str, z: &str) -> RationalPoint {
    [q(x), q(y), q(z)]
}

#[test]
fn nonclamped_weighted_curve_has_exact_endpoint_controls() {
    let dec = Decider::default();
    let curve = RationalSplineCurve::new(
        2,
        [0, 1, 2, 3, 4, 5].map(Rational::from_i64).to_vec(),
        vec![p("0", "0", "0"), p("2", "4", "0"), p("4", "0", "2")],
        vec![q("1"), q("2"), q("1")],
        &dec,
    )
    .unwrap();
    let extracted = curve.extracted(curve.domain(), &dec).unwrap();
    assert_eq!(
        extracted.knots(),
        &[q("2"), q("2"), q("2"), q("3"), q("3"), q("3")]
    );
    // Independent homogeneous endpoint insertion for this uniform quadratic.
    let expected = vec![
        [q("2"), q("4"), q("0"), q("1.5")],
        [q("4"), q("8"), q("0"), q("2")],
        [q("4"), q("4"), q("1"), q("1.5")],
    ];
    assert_eq!(
        weighted(extracted.controls(), extracted.weights(), &dec).unwrap(),
        expected
    );
    let mut dec = Decider::default();
    let evaluator = curve.to_curve(&mut dec).unwrap();
    let t = q("2.5");
    let value = evaluator.evaluate_box([&t, &t], 128, &dec).unwrap();
    for (actual, expected) in value
        .iter()
        .zip([q("2"), q("24") / q("7"), q("1") / q("7")])
    {
        assert!(actual.lo() <= expected && actual.hi() >= expected);
    }
}

#[test]
fn extraction_preserves_a_sub_picometre_parameter_span() {
    let dec = Decider::default();
    let narrow = q("1.000000000000000000000000000001");
    let curve = RationalSplineCurve::new(
        1,
        vec![q("0"), q("0"), q("1"), narrow.clone(), q("2"), q("2")],
        vec![
            p("0", "0", "0"),
            p("1", "0", "0"),
            p("2", "0", "0"),
            p("3", "0", "0"),
        ],
        vec![q("1"); 4],
        &dec,
    )
    .unwrap();
    let piece = curve.extracted([&q("1"), &narrow], &dec).unwrap();
    assert_eq!(piece.knots(), &[q("1"), q("1"), narrow.clone(), narrow]);
    assert_eq!(piece.controls(), &[p("1", "0", "0"), p("2", "0", "0")]);
}

#[test]
fn tensor_extraction_preserves_weighting_and_uv_order() {
    let mut dec = Decider::default();
    let mut controls = Vec::new();
    let mut weights = Vec::new();
    for (x, y, w) in [("0", "0", "1"), ("2", "4", "2"), ("4", "0", "1")] {
        for z in ["0", "2"] {
            controls.push(p(x, y, z));
            weights.push(q(w));
        }
    }
    let surface = RationalSplineSurface::new(
        [2, 1],
        [
            [0, 1, 2, 3, 4, 5].map(Rational::from_i64).to_vec(),
            vec![q("2"), q("2"), q("4"), q("4")],
        ],
        [3, 2],
        controls,
        weights,
        &dec,
    )
    .unwrap();
    let evaluator = surface.to_surface(&mut dec).unwrap();
    let (u, v) = (q("2.5"), q("3"));
    let value = evaluator
        .jet_box([&u, &u], [&v, &v], 128, &dec)
        .unwrap()
        .value;
    for (actual, expected) in value.iter().zip([q("2"), q("24") / q("7"), q("1")]) {
        assert!(actual.lo() <= expected && actual.hi() >= expected);
    }
}

#[test]
fn malformed_knots_weights_and_outside_extraction_refuse() {
    let dec = Decider::default();
    let controls = vec![p("0", "0", "0"), p("1", "0", "0")];
    for weights in [vec![q("1"), q("0")], vec![q("1"), q("-1")], vec![q("1")]] {
        assert!(
            RationalSplineCurve::new(
                1,
                vec![q("0"), q("0"), q("1"), q("1")],
                controls.clone(),
                weights,
                &dec
            )
            .is_err()
        );
    }
    let curve = RationalSplineCurve::new(
        1,
        vec![q("0"), q("0"), q("1"), q("1")],
        controls,
        vec![q("1"); 2],
        &dec,
    )
    .unwrap();
    assert!(curve.extracted([&q("-0.1"), &q("1")], &dec).is_err());
    assert!(curve.extracted([&q("1"), &q("1")], &dec).is_err());
}
