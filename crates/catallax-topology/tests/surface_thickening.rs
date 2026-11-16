use catallax_numbers::{Decider, Rational, Real, with_operation_budget};
use catallax_shapes::{SplineSupport, Surface, spline_tools::RationalSplineSurface};
use catallax_topology::numerical::{NumericalBody, NumericalOptions, surface::SplineSurfaceBody};

fn q(text: &str) -> Rational {
    Rational::parse(text).unwrap()
}

fn cap_point(body: &NumericalBody, key: &str, dec: &mut Decider) -> [Rational; 3] {
    let face = body
        .definition()
        .faces
        .iter()
        .find(|face| face.key == key)
        .unwrap();
    let Surface::BSpline(SplineSupport::RationalPatch(patch)) = &face.surface else {
        panic!("expected the retained spline cap");
    };
    let point = patch
        .point_at(&[Real::from(q("34/125")), Real::from(q("1/2"))], dec)
        .unwrap();
    [&point.x, &point.y, &point.z]
        .map(|coordinate| dec.rational_value(coordinate).unwrap().unwrap())
}

#[test]
fn thickening_retains_a_bump_between_global_sample_nodes() {
    with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let knots = [
            "0", "0", "0", "0", "270/1000", "271/1000", "272/1000", "273/1000", "274/1000", "1",
            "1", "1", "1",
        ]
        .map(q)
        .to_vec();
        let mut controls = Vec::new();
        for i in 0..9 {
            let x = &(&knots[i + 1] + &knots[i + 2]) + &knots[i + 3];
            let x = &x / &q("3");
            let z = if i == 4 { q("1/10000") } else { q("0") };
            for y in [q("0"), q("1")] {
                controls.push([x.clone(), y, z.clone()]);
            }
        }
        let source = RationalSplineSurface::new(
            [3, 1],
            [knots, ["0", "0", "1", "1"].map(q).to_vec()],
            [9, 2],
            controls,
            vec![q("1"); 18],
            &dec,
        )
        .unwrap();
        let tolerance = q("1/1000000");
        let sheet = SplineSurfaceBody::new(
            "local-bump",
            source,
            NumericalOptions {
                tolerance_mm: tolerance.clone(),
            },
            &mut dec,
        )
        .unwrap();
        let solid = sheet.thicken(q("1/10000"), &mut dec).unwrap();
        // The cubic basis peaks at 2/3 here, with zero slope and an upward normal.
        let expected = [q("34/125"), q("1/2"), q("1/15000")];
        assert_eq!(cap_point(&solid, "local-bump", &mut dec), expected);
        let offset = cap_point(&solid, "local-bump/offset", &mut dec);
        let expected = [
            expected[0].clone(),
            expected[1].clone(),
            &expected[2] + &q("1/10000"),
        ];
        for (actual, target) in offset.iter().zip(&expected) {
            assert!(
                (actual - target).abs() <= tolerance,
                "offset cap erased the knot-local bump"
            );
        }
    });
}
