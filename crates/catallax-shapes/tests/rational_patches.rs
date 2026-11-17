//! Exact rational chart proofs and their complete finite consumers.

use catallax_numbers::{Decider, Rational, Real};
use catallax_shapes::{BSplineSurface, Plane, Point3, RationalPatch, SectionSurface, Vector3};

fn r(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn p(x: &str, y: &str, z: &str) -> Point3 {
    Point3::new(r(x).into(), r(y).into(), r(z).into())
}
fn linear_knots() -> Vec<Rational> {
    [r("0"), r("0"), r("1"), r("1")].to_vec()
}

#[test]
fn display_stations_keep_exact_order_and_deduplicate_algebraic_values() {
    let mut dec = Decider::default();
    let patch = RationalPatch::new(
        BSplineSurface::new(
            [1, 1],
            [linear_knots(), linear_knots()],
            [2, 2],
            vec![
                p("0", "0", "0"),
                p("0", "1", "0"),
                p("1", "0", "0"),
                p("1", "1", "0"),
            ],
            vec![Real::one(); 4],
            &mut dec,
        )
        .unwrap(),
        &mut dec,
    )
    .unwrap();
    let algebraic = Real::from_rational(r("1/2")).sqrt();
    let required = [
        r("3/4").into(),
        algebraic.clone(),
        r("1/4").into(),
        algebraic.clone(),
        Real::one(),
        Real::zero(),
    ];
    let grid = patch
        .display_grid(&r("1e-6"), None, None, 16, [&required, &[]], &mut dec)
        .unwrap();
    let expected = [
        Real::zero(),
        r("1/4").into(),
        algebraic,
        r("3/4").into(),
        Real::one(),
    ];
    assert_eq!(grid.parameters[0].len(), expected.len());
    for (actual, expected) in grid.parameters[0].iter().zip(expected) {
        assert!(matches!(
            dec.sign(&(actual - &expected)),
            catallax_numbers::Outcome::ProvenCoincident(_)
        ));
    }
    assert!(grid.grid.chord <= r("1e-6"));
    assert!(grid.embedding.accepts_rounding(&r("1e-12")));
}

#[test]
fn a_g1_speed_jump_keeps_angular_mesh_control() {
    let mut dec = Decider::default();
    let carrier = BSplineSurface::new(
        [1, 1],
        [["0", "0", "1/2", "1", "1"].map(r).to_vec(), linear_knots()],
        [3, 2],
        vec![
            p("0", "0", "0"),
            p("0", "1", "0"),
            p("1/4", "0", "0"),
            p("1/4", "1", "0"),
            p("1", "0", "0"),
            p("1", "1", "0"),
        ],
        vec![Real::one(); 6],
        &mut dec,
    )
    .unwrap();
    let patch = RationalPatch::new(carrier, &mut dec).unwrap();
    let grid = patch
        .display_grid(&r("1e-6"), Some(&r("1/10")), None, 64, [&[], &[]], &mut dec)
        .unwrap();
    assert!(grid.grid.chord <= r("1e-30"));
    assert!(grid.embedding.accepts_rounding(&r("1e-12")));
    for normal in grid.normal_boxes.as_ref().expect("requested angular proof") {
        assert!(normal[0].contains(&Rational::zero()));
        assert!(normal[1].contains(&Rational::zero()));
        assert!(normal[2].lo() > Rational::zero());
    }
}

fn saddle(weighted: bool, dec: &mut Decider) -> RationalPatch {
    RationalPatch::new(
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
            [
                "1",
                if weighted { "2" } else { "1" },
                "1",
                if weighted { "2" } else { "1" },
            ]
            .map(|q| Real::from_rational(r(q)))
            .to_vec(),
            dec,
        )
        .unwrap(),
        dec,
    )
    .unwrap()
}

fn plane(z: &str) -> Plane {
    Plane {
        origin: p("0", "0", z),
        normal: Vector3::new(Real::zero(), Real::zero(), Real::one()),
        ref_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
    }
}

#[test]
fn positive_rational_weights_keep_exact_positions_and_finite_incidence() {
    let mut dec = Decider::default();
    let patch = saddle(true, &mut dec);
    let point = patch
        .point_at(&[r("1/2").into(), r("1/2").into()], &mut dec)
        .unwrap();
    assert_eq!(
        point.rational_components().unwrap().map(Clone::clone),
        [r("1/2"), r("2/3"), r("1/3")]
    );
    assert!(patch.parameters_at(&point, &mut dec).unwrap().is_some());
    assert!(
        patch
            .parameters_at(
                &p("1/2", "2/3", "0.333333333333333333333333333333"),
                &mut dec
            )
            .unwrap()
            .is_none()
    );
    assert!(
        patch
            .parameters_at(&p("2", "1/2", "1"), &mut dec)
            .unwrap()
            .is_none()
    );
    assert!(
        patch
            .point_at(&[r("-1/100").into(), r("1/2").into()], &mut dec)
            .is_err()
    );
}

#[test]
fn transverse_line_roots_are_complete_and_tangency_is_not_a_false_miss() {
    let mut dec = Decider::default();
    let patch = saddle(false, &mut dec);
    let z = Vector3::new(Real::zero(), Real::zero(), Real::one());
    let hits = patch
        .line_hits(&p("1/2", "1/2", "-1"), &z, &mut dec)
        .unwrap();
    assert_eq!(
        hits.iter()
            .map(|p| p.rational_components().unwrap().map(Clone::clone))
            .collect::<Vec<_>>(),
        vec![[r("1/2"), r("1/2"), r("1/4")]]
    );
    assert!(
        patch
            .line_hits(&p("2", "2", "-1"), &z, &mut dec)
            .unwrap()
            .is_empty()
    );
    let tangent = Vector3::new(Real::one(), Real::one(), Real::zero());
    assert_eq!(
        patch
            .line_hits(&p("0", "0", "0"), &tangent, &mut dec)
            .unwrap_err()
            .code,
        "shapes.rational-patch-unproven"
    );
    let contained = Vector3::new(Real::zero(), Real::one(), Real::zero());
    assert_eq!(
        patch
            .line_hits(&p("0", "0", "0"), &contained, &mut dec)
            .unwrap_err()
            .code,
        "shapes.rational-patch-unproven"
    );
}

#[test]
fn hyperbolic_plane_sections_retain_the_exact_curve_and_definition() {
    let mut dec = Decider::default();
    let patch = saddle(false, &mut dec);
    let section = patch
        .sections(&SectionSurface::Plane(plane("1/4")), &mut dec)
        .unwrap();
    assert!(section.patch().shares_definition_nodes(&patch));
    let pieces: Vec<_> = section.curves().collect();
    assert_eq!(pieces.len(), 1);
    let (curve, crossing) = pieces[0];
    assert!(crossing);
    assert!(section.points().is_empty());
    for parameter in curve.domain() {
        let point = curve
            .evaluate_box([parameter, parameter], 160, &dec)
            .unwrap();
        assert!(point[2].contains(&r("1/4")));
        assert!(
            point
                .iter()
                .all(|coordinate| coordinate.hi() - coordinate.lo() < r("1e-30"))
        );
        assert!(point[0].lo() >= r("1/4") - r("1e-30"));
        assert!(point[1].hi() <= r("1") + r("1e-30"));
    }
}

#[test]
fn folded_and_creased_definitions_cannot_acquire_an_embedding_certificate() {
    let mut dec = Decider::default();
    let folded = BSplineSurface::new(
        [2, 1],
        [
            ["0", "0", "0", "1", "1", "1"].map(r).to_vec(),
            linear_knots(),
        ],
        [3, 2],
        vec![
            p("0", "0", "0"),
            p("0", "1", "0"),
            p("1", "0", "0"),
            p("1", "1", "0"),
            p("0", "0", "0"),
            p("0", "1", "0"),
        ],
        vec![Real::one(); 6],
        &mut dec,
    )
    .unwrap();
    assert_eq!(
        RationalPatch::new(folded, &mut dec).unwrap_err().code,
        "shapes.rational-patch-unproven"
    );
    let creased = BSplineSurface::new(
        [1, 1],
        [["0", "0", "1/2", "1", "1"].map(r).to_vec(), linear_knots()],
        [3, 2],
        vec![
            p("0", "0", "0"),
            p("0", "1", "0"),
            p("1/2", "0", "0"),
            p("1/2", "1", "0"),
            p("1", "0", "1/2"),
            p("1", "1", "1/2"),
        ],
        vec![Real::one(); 6],
        &mut dec,
    )
    .unwrap();
    assert_eq!(
        RationalPatch::new(creased, &mut dec).unwrap_err().code,
        "shapes.rational-patch-unproven"
    );
}

fn cubic_sheet(cubic: [[&str; 2]; 2], dec: &mut Decider) -> BSplineSurface {
    let cubic = cubic.map(|row| row.map(r));
    let linear = ["-1", "-1/3", "1/3", "1"].map(r);
    let cube = ["-1", "1", "-1", "1"].map(r);
    let square = ["1", "-1/3", "-1/3", "1"].map(r);
    let mut points = Vec::with_capacity(16);
    for u in 0..4 {
        for v in 0..4 {
            let x = &cubic[0][0] * &cube[u] + &cubic[0][1] * &cube[v];
            let y = &cubic[1][0] * &cube[u] + &cubic[1][1] * &cube[v];
            points.push(Point3::new(
                (&linear[u] + &x).into(),
                (&linear[v] + &y).into(),
                (&square[u] * &square[v] * r("1/10")).into(),
            ));
        }
    }
    BSplineSurface::new(
        [3, 3],
        std::array::from_fn(|_| ["0", "0", "0", "0", "1", "1", "1", "1"].map(r).to_vec()),
        [4, 4],
        points,
        vec![Real::one(); 16],
        dec,
    )
    .unwrap()
}

#[test]
fn a_positive_projected_jacobian_need_not_be_diagonally_dominant() {
    let mut dec = Decider::default();
    // In coordinates u,v in [-1,1], the projection is
    // (u + 2v^3, v + 16v^3/3). Its symmetric Jacobian has determinant
    // 1 + 16v^2 - 9v^4 >= 1, but Gershgorin fails near |v| = 1.
    // The nonzero u^2 v^2 height prevents the ruled-chart special case.
    let patch =
        RationalPatch::new(cubic_sheet([["0", "2"], ["0", "16/3"]], &mut dec), &mut dec).unwrap();
    assert!(patch.embedding_lower_bound() > &Rational::zero());
    let point = patch
        .point_at(&[r("3/4").into(), r("3/4").into()], &mut dec)
        .unwrap();
    assert_eq!(
        point.rational_components().unwrap().map(Clone::clone),
        [r("3/4"), r("7/6"), r("1/160")]
    );
}

#[test]
fn indefinite_or_negative_projected_diagonals_do_not_certify_a_chart() {
    for cubic in [[["0", "2"], ["2", "16/3"]], [["-1", "0"], ["0", "-1"]]] {
        let mut dec = Decider::default();
        assert_eq!(
            RationalPatch::new(cubic_sheet(cubic, &mut dec), &mut dec)
                .unwrap_err()
                .code,
            "shapes.rational-patch-unproven"
        );
    }
}

#[test]
fn huge_translated_non_dyadic_knots_do_not_lose_the_chart() {
    let mut dec = Decider::default();
    let lo = r("1e100");
    let hi = &lo + &r("1/3");
    let patch = RationalPatch::new(
        BSplineSurface::new(
            [1, 1],
            [vec![lo.clone(), lo, hi.clone(), hi], linear_knots()],
            [2, 2],
            vec![
                p("0", "0", "0"),
                p("0", "1", "0"),
                p("1", "0", "0"),
                p("1", "1", "0"),
            ],
            vec![Real::one(); 4],
            &mut dec,
        )
        .unwrap(),
        &mut dec,
    )
    .unwrap();
    let grid = patch
        .display_grid(
            &r("1e-20"),
            Some(&r("1/10")),
            None,
            16,
            [&[], &[]],
            &mut dec,
        )
        .unwrap();
    assert!(grid.grid.chord < r("1e-30"));
    let corners = [
        [r("0"), r("0"), r("0")],
        [r("0"), r("1"), r("0")],
        [r("1"), r("0"), r("0")],
        [r("1"), r("1"), r("0")],
    ];
    for (actual, expected) in grid.grid.points.iter().zip(corners) {
        for axis in 0..3 {
            assert!(actual[axis].contains(&expected[axis]));
        }
    }
}

#[test]
fn exact_ranges_enclose_geometry_rather_than_the_control_polygon() {
    let mut dec = Decider::default();
    let patch = RationalPatch::new(
        BSplineSurface::new(
            [2, 1],
            [
                ["0", "0", "0", "1", "1", "1"].map(r).to_vec(),
                linear_knots(),
            ],
            [3, 2],
            vec![
                p("0", "0", "1"),
                p("0", "1", "1"),
                p("1/2", "0", "3/2"),
                p("1/2", "1", "3/2"),
                p("1", "0", "1"),
                p("1", "1", "1"),
            ],
            vec![Real::one(); 6],
            &mut dec,
        )
        .unwrap(),
        &mut dec,
    )
    .unwrap();
    let ranges = patch.ranges(&r("1e-7"), &mut dec).unwrap();
    for (axis, expected) in [[r("0"), r("1")], [r("0"), r("1")], [r("1"), r("5/4")]]
        .iter()
        .enumerate()
    {
        for end in 0..2 {
            assert!(ranges[axis][end][0] <= expected[end] && ranges[axis][end][1] >= expected[end]);
            assert!(&ranges[axis][end][1] - &ranges[axis][end][0] <= r("1e-7"));
        }
    }
}
