use super::*;
use catallax_numbers::{BudgetPolicy, Real};
use catallax_shapes::Point3;

fn q(text: &str) -> Rational {
    Rational::parse(text).unwrap()
}

fn curve(degree: usize, knots: &[&str], x: &[&str], weights: &[&str]) -> BSplineCurve {
    BSplineCurve::new(
        degree,
        knots.iter().map(|s| q(s)).collect(),
        x.iter()
            .map(|s| Point3::new(q(s).into(), (-7).into(), 0.into()))
            .collect(),
        weights.iter().map(|s| Real::from(q(s))).collect(),
        &mut Decider::default(),
    )
    .unwrap()
}

fn check(curve: &BSplineCurve, domain: [&Rational; 2], expected: [&str; 2]) {
    let accuracy = Accuracy {
        relative: q("0"),
        absolute: q("1e-12"),
    };
    let dec = Decider::default();
    let extrema = curve_extrema(
        curve,
        domain,
        &accuracy,
        &mut Budget::new(&dec).unwrap(),
        &dec,
    )
    .unwrap();
    assert!(extrema.reached(&accuracy));
    for (enclosure, exact) in [
        (&extrema.min[0], q(expected[0])),
        (&extrema.max[0], q(expected[1])),
    ] {
        assert!(
            enclosure.lo <= exact && exact <= enclosure.hi,
            "{enclosure:?}, {exact}"
        );
    }
    for extrema in [&extrema.min, &extrema.max] {
        assert!(extrema[1].lo <= q("-7") && extrema[1].hi >= q("-7"));
        assert!(extrema[2].lo <= q("0") && extrema[2].hi >= q("0"));
    }
}

#[test]
fn cubic_interior_extremum_converges_beyond_its_control_hull() {
    // x=9t(1-t), max=9/4 at t=1/2; the control hull would give 3.
    let curve = curve(
        3,
        &["0", "0", "0", "0", "1", "1", "1", "1"],
        &["0", "3", "3", "0"],
        &["1"; 4],
    );
    check(&curve, curve.domain(), ["0", "9/4"]);
    // The admitted finite subinterval excludes the carrier's interior maximum.
    check(&curve, [&q("0"), &q("1/4")], ["0", "27/16"]);
}

#[test]
fn positive_rational_weights_bound_the_quotient_extremum() {
    // x=12t(1-t)/(1+2t(1-t)); max=2 at t=1/2, not hull maximum 3.
    let curve = curve(
        2,
        &["0", "0", "0", "1", "1", "1"],
        &["0", "3", "0"],
        &["1", "2", "1"],
    );
    check(&curve, curve.domain(), ["0", "2"]);
}

#[test]
fn repeated_knots_split_both_one_sided_extrema_and_negative_coordinates() {
    let curve = curve(
        2,
        &["0", "0", "0", "1/3", "1/3", "1", "1", "1"],
        &["0", "2", "0", "-4", "0"],
        &["1"; 5],
    );
    check(&curve, curve.domain(), ["-2", "1"]);
    check(&curve, [&q("1/3"), &q("1")], ["-2", "0"]);
}

#[test]
fn insufficient_measurement_work_never_publishes_an_unconverged_hull() {
    let curve = curve(
        3,
        &["0", "0", "0", "0", "1", "1", "1", "1"],
        &["0", "3", "3", "0"],
        &["1"; 4],
    );
    let dec = Decider::new(BudgetPolicy {
        max_precision_bits: 256,
        max_work: 32,
    });
    let accuracy = Accuracy {
        relative: q("0"),
        absolute: q("1e-30"),
    };
    let error = curve_extrema(
        &curve,
        curve.domain(),
        &accuracy,
        &mut Budget::new(&dec).unwrap(),
        &dec,
    )
    .err()
    .unwrap();
    assert_eq!(error.code, "measure.budget-exhausted");
}
