//! Finite graph queries retain chart bounds and continuous overlap scope.

use catallax_numbers::{Decider, Outcome, Rational, Real, ZeroProof, with_operation_budget};
use catallax_shapes::{BSplineCurve, BSplineSurface, Point3, RationalAtlas, SectionHits};

fn q(n: i64, d: i64) -> Real {
    Rational::from_ratio(n, d).unwrap().into()
}

fn graph(dec: &mut Decider) -> RationalAtlas {
    let controls = (0..3)
        .flat_map(|u| {
            (0..2).map(move |v| {
                Point3::new(
                    q(u, 2),
                    Real::from_i64(i64::from(u == 2) * v),
                    Real::from_i64(v),
                )
            })
        })
        .collect();
    let source = BSplineSurface::new(
        [2, 1],
        [
            [0, 0, 0, 1, 1, 1].map(Rational::from_i64).to_vec(),
            [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
        ],
        [3, 2],
        controls,
        vec![Real::one(); 6],
        dec,
    )
    .unwrap();
    RationalAtlas::partition(source, [false, false], dec)
        .unwrap()
        .remove(0)
}

fn line(from: Point3, to: Point3, dec: &mut Decider) -> BSplineCurve {
    BSplineCurve::new(
        1,
        [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
        vec![from, to],
        vec![Real::one(); 2],
        dec,
    )
    .unwrap()
}

#[test]
fn spline_graph_queries_keep_crossings_touches_and_off_domain_controls() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let atlas = graph(&mut dec);
        for (height, z, expected) in [
            (q(1, 8), q(1, 2), 1),
            (Real::zero(), q(1, 2), 1),
            (q(-1, 8), q(1, 2), 0),
            (q(1, 8), q(3, 2), 0),
        ] {
            let spline = line(
                Point3::new(Real::zero(), height.clone(), z.clone()),
                Point3::new(Real::one(), height, z),
                &mut dec,
            );
            let SectionHits::Points(points) =
                atlas.spline_hits(&spline, &mut dec).unwrap().unwrap()
            else {
                panic!("horizontal curve has no contained interval")
            };
            assert_eq!(points.len(), expected);
            for point in points {
                let at = atlas.charts()[0]
                    .patch()
                    .parameters_at(&point, &mut dec)
                    .unwrap()
                    .unwrap();
                let actual = atlas.point_at(&at, &mut dec).unwrap();
                assert!(
                    matches!(dec.sign(&actual.minus(&point).norm_squared()),Outcome::ProvenCoincident(proof) if proof.proof!=ZeroProof::SeparationBound)
                );
                assert!(spline.parameter_of(&point, &mut dec).unwrap().is_some());
            }
        }
    });
    assert!(!work.exhausted, "{} work units", work.used);
}

#[test]
fn containment_needs_the_entire_spline_domain() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let atlas = graph(&mut dec);
        for (end, contained) in [(1, true), (2, false)] {
            let spline = line(
                Point3::new(Real::zero(), Real::zero(), Real::zero()),
                Point3::new(Real::from_i64(end), Real::zero(), Real::zero()),
                &mut dec,
            );
            let answer = atlas.spline_hits(&spline, &mut dec);
            if contained {
                assert!(matches!(answer.unwrap(), Some(SectionHits::Contained)));
            } else {
                assert!(
                    answer
                        .unwrap_err()
                        .detail
                        .contains("partial spline containment")
                );
            }
        }
    });
    assert!(!work.exhausted, "{} work units", work.used);
}
