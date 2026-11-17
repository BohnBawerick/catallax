use super::{BSplineCurve, BSplineSurface, Point3, SplineBoundary};
use catallax_numbers::{Ball, Decider, Rational, Real};

fn r(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn joined_curves(kink: bool, reversed: bool) -> Vec<(BSplineCurve, Rational)> {
    let mut dec = Decider::default();
    let ys = if kink {
        ["1000", "1000", "1000", "1001", "1002"]
    } else {
        ["1000"; 5]
    };
    let mut controls: Vec<_> = ["0", "1/2", "1", "3/2", "2"]
        .into_iter()
        .zip(ys)
        .map(|(x, y)| Point3::new(r(x).into(), r(y).into(), Real::zero()))
        .collect();
    if reversed {
        controls.reverse();
    }
    let curve = BSplineCurve::new(
        2,
        ["0", "0", "0", "1/2", "1/2", "1", "1", "1"]
            .map(r)
            .to_vec(),
        controls,
        ["1", "1/2", "1", "1/2", "1"]
            .map(|value| r(value).into())
            .to_vec(),
        &mut dec,
    )
    .unwrap();
    let mut cases = vec![(curve.clone(), Rational::zero())];
    for axis in 0..2 {
        let mut degrees = [1, 1];
        degrees[axis] = curve.degree();
        let mut knots: [Vec<Rational>; 2] =
            std::array::from_fn(|_| ["0", "0", "1", "1"].map(r).to_vec());
        knots[axis] = curve.knots().to_vec();
        let mut count = [2, 2];
        count[axis] = curve.control_count();
        let mut controls = Vec::with_capacity(count[0] * count[1]);
        let mut weights = Vec::with_capacity(controls.capacity());
        for u in 0..count[0] {
            for v in 0..count[1] {
                let index = if axis == 0 { u } else { v };
                let fixed = if axis == 0 { v } else { u };
                let mut point = curve.control(index).clone();
                point.z = Real::from_i64(fixed as i64);
                controls.push(point);
                weights.push(curve.weight(index).clone());
            }
        }
        let surface = BSplineSurface::new(degrees, knots, count, controls, weights, &mut dec).unwrap();
        let sides = if axis == 0 {
            [SplineBoundary::VMin, SplineBoundary::VMax]
        } else {
            [SplineBoundary::UMin, SplineBoundary::UMax]
        };
        for (fixed, side) in sides.into_iter().enumerate() {
            cases.push((BSplineCurve::boundary(&surface, side), Rational::from_i64(fixed as i64)));
        }
    }
    cases
}

fn contains_tightly(value: &Ball, lower: &Rational, upper: &Rational) {
    assert!(value.contains(lower) && value.contains(upper), "{value:?}");
    assert!(
        value.hi() - value.lo() <= upper - lower + r("1/1000000000"),
        "{value:?}"
    );
}

fn check_join(kink: bool) {
    let dec = Decider::default();
    let at = r("1/2");
    for reversed in [false, true] {
        let x_first = r(if reversed { "-1" } else { "1" });
        let y_first = if kink {
            if reversed { [r("-2"), r("0")] } else { [r("0"), r("2")] }
        } else {
            [r("0"), r("0")]
        };
        for (curve, z) in joined_curves(kink, reversed) {
            for epsilon in [r("0"), r("1/1000000000000000000000000")] {
                let lo = &at - &epsilon;
                let hi = &at + &epsilon;
                let jet = curve.jet_box([&lo, &hi], 128, &dec).unwrap();
                let first = curve.first_box_reusing([&lo, &hi], 128, &dec, &mut Vec::new()).unwrap();
                let position = curve.evaluate_box([&lo, &hi], 128, &dec).unwrap();
                for value in [&jet.value, &first[0], &position] {
                    for (actual, expected) in value.iter().zip([r("1"), r("1000"), z.clone()]) {
                        contains_tightly(actual, &expected, &expected);
                    }
                }
                for derivative in [&jet.first, &first[1]] {
                    contains_tightly(&derivative[0], &x_first, &x_first);
                    contains_tightly(&derivative[1], &y_first[0], &y_first[1]);
                    contains_tightly(&derivative[2], &r("0"), &r("0"));
                }
                contains_tightly(&jet.second[0], &r("-8"), &r("8"));
                contains_tightly(&jet.second[1], &r("0"), &r(if kink { "16" } else { "0" }));
                contains_tightly(&jet.second[2], &r("0"), &r("0"));
            }
        }
    }
}

#[test]
fn rational_knot_jets_preserve_cartesian_cancellation() {
    check_join(false);
}

#[test]
fn rational_knot_jets_keep_both_genuine_corner_limits() {
    check_join(true);
}
