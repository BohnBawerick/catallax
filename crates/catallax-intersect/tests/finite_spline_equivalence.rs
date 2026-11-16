//! Finite spline equality retains whole spans, domains and opposed senses.

use catallax_intersect::same_curve;
use catallax_numbers::{Decider, Outcome, Rational, Real, ZeroProof, with_operation_budget};
use catallax_shapes::{Apart, BSplineCurve, Curve, GeomStore, Point3};

fn q(value: &str) -> Real {
    Rational::parse(value).unwrap().into()
}

fn bezier(
    controls: &[[&str; 2]],
    domain: [&str; 2],
    weight: i64,
    dec: &mut Decider,
) -> BSplineCurve {
    let degree = controls.len() - 1;
    BSplineCurve::new(
        degree,
        domain
            .into_iter()
            .flat_map(|value| std::iter::repeat_n(Rational::parse(value).unwrap(), degree + 1))
            .collect(),
        controls
            .iter()
            .map(|[x, y]| Point3::new(q(x), q(y), Real::zero()))
            .collect(),
        vec![Real::from_i64(weight); controls.len()],
        dec,
    )
    .unwrap()
}

#[test]
fn full_span_identities_cover_inserted_knots_degree_weights_and_reversal() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let mut geom = GeomStore::new();
        let base = bezier(
            &[["0", "0"], ["1", "2"], ["2", "0"]],
            ["3", "7"],
            1,
            &mut dec,
        );
        let a = geom.add_curve(Curve::BSpline(base));
        let left = bezier(
            &[["0", "0"], ["1/2", "1"], ["1", "1"]],
            ["-9", "-2"],
            3,
            &mut dec,
        );
        let right = bezier(
            &[["1", "1"], ["3/2", "1"], ["2", "0"]],
            ["2", "10"],
            5,
            &mut dec,
        );
        let joined =
            BSplineCurve::joined(&[(left.clone(), true), (right.clone(), true)], &mut dec).unwrap();
        let reversed =
            BSplineCurve::joined(&[(right, false), (left.clone(), false)], &mut dec).unwrap();
        let elevated = bezier(
            &[["0", "0"], ["2/3", "4/3"], ["4/3", "4/3"], ["2", "0"]],
            ["-2", "6"],
            7,
            &mut dec,
        );
        for other in [joined, reversed, elevated] {
            let b = geom.add_curve(Curve::BSpline(other));
            for (a, b) in [(a, b), (b, a)] {
                assert!(
                    matches!(same_curve(&geom, a, b, &mut dec), Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound)
                );
            }
        }
        let partial = geom.add_curve(Curve::BSpline(left));
        for (a, b) in [(a, partial), (partial, a)] {
            assert!(matches!(
                same_curve(&geom, a, b, &mut dec),
                Outcome::Proven(Apart)
            ));
        }
    });
    assert!(!work.exhausted, "{} work units", work.used);
}

#[test]
fn off_curve_witnesses_prove_difference_but_matching_samples_never_prove_equality() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let mut geom = GeomStore::new();
        let a = geom.add_curve(Curve::BSpline(bezier(
            &[["0", "0"], ["1", "2"], ["2", "0"]],
            ["0", "1"],
            1,
            &mut dec,
        )));
        let off = geom.add_curve(Curve::BSpline(bezier(
            &[["0", "0"], ["1", "2000001/1000000"], ["2", "0"]],
            ["0", "1"],
            1,
            &mut dec,
        )));
        for (a, b) in [(a, off), (off, a)] {
            assert!(matches!(
                same_curve(&geom, a, b, &mut dec),
                Outcome::Proven(Apart)
            ));
        }
        let mut hidden = bezier(
            &[["0", "0"], ["2/3", "4/3"], ["4/3", "4/3"], ["2", "0"]],
            ["0", "1"],
            1,
            &mut dec,
        )
        .controls()
        .cloned()
        .collect::<Vec<_>>();
        hidden[1].y = &hidden[1].y - &q("1/6000000");
        hidden[2].y = &hidden[2].y + &q("1/6000000");
        let hidden = BSplineCurve::new(
            3,
            [0, 0, 0, 0, 1, 1, 1, 1].map(Rational::from_i64).to_vec(),
            hidden,
            vec![Real::one(); 4],
            &mut dec,
        )
        .unwrap();
        let b = geom.add_curve(Curve::BSpline(hidden));
        for (a, b) in [(a, b), (b, a)] {
            assert!(
                matches!(same_curve(&geom, a, b, &mut dec), Outcome::Refused(error) if error.code == "intersect.parametric-coincidence")
            );
        }
    });
    assert!(!work.exhausted, "{} work units", work.used);
}
