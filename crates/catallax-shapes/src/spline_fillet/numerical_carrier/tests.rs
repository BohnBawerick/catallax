use super::*;

fn polynomial(coefficients: &[i64], t: &Rational) -> Rational {
    coefficients.iter().fold(Rational::zero(), |value, coefficient| {
        &value * t + Rational::from_i64(*coefficient)
    })
}

#[test]
fn nonuniform_station_jets_preserve_quartics_and_constant_homogeneous_weights() {
    let parameters = ["-1", "-1/3", "0", "1/7", "2/3", "1", "2"]
        .map(|value| Rational::parse(value).unwrap());
    for offset in [0, 1_000_000_000_000] {
        let offset = Rational::from_i64(offset);
        let values: Vec<_> = parameters.iter().map(|t| [
            polynomial(&[-5, 4, -3, 2, 1], t) + offset.clone(),
            Rational::from_ratio(7, 5).unwrap() + offset.clone(),
            t * t - offset.clone(),
        ]).collect();
        let slopes = station_slopes(parameters.len(), |i| &parameters[i], |i| &values[i]);
        for (t, slope) in parameters.iter().zip(slopes) {
            assert_eq!(slope, [
                polynomial(&[-20, 12, -6, 2], t),
                Rational::zero(),
                t.mul_pow2(1),
            ]);
        }
    }
}

#[test]
fn initial_three_and_four_station_fits_keep_their_exact_lower_order_jets() {
    let parameters = ["-3/4", "-1/5", "2/7", "3/2"]
        .map(|value| Rational::parse(value).unwrap());
    for count in [3, 4] {
        let degree = count - 1;
        let values: Vec<_> = parameters[..count].iter().map(|t| [
            (0..degree).fold(Rational::one(), |value, _| &value * t),
        ]).collect();
        let slopes = station_slopes(count, |i| &parameters[i], |i| &values[i]);
        for (t, slope) in parameters[..count].iter().zip(slopes) {
            let expected = (0..degree - 1).fold(Rational::from_i64(degree as i64), |value, _| &value * t);
            assert_eq!(slope, [expected]);
        }
    }
}

fn analytic_row(t: &Rational) -> Packed {
    let square = t * t;
    let x = &Rational::from_i64(1000) + &square;
    let y = &Rational::from_i64(-1000) + t;
    let mut result: Packed = std::array::from_fn(|_| Rational::zero());
    for row in 0..3 {
        let weight = if row == 1 {
            (&Rational::from_i64(5) + t).mul_pow2(-3)
        } else {
            Rational::one()
        };
        let point = [
            &x + &Rational::from_i64([0, 1, 3][row]),
            &y + &Rational::from_i64((2 * row) as i64),
            &square + &(t * &Rational::from_i64(row as i64)),
        ];
        for axis in 0..3 {
            result[4 * row + axis] = &point[axis] * &weight;
        }
        result[4 * row + 3] = weight;
    }
    for side in 0..2 {
        result[12 + 2 * side] =
            &square + &Rational::from_ratio(side as i64, 3).unwrap();
        result[13 + 2 * side] = &Rational::from_i64(side as i64 + 1) - t;
    }
    result
}

#[test]
fn batched_sections_and_native_rails_keep_the_analytic_rational_patch() {
    let parameters = ["-2/3", "0", "1/7", "4/5", "3/2"]
        .map(|value| Rational::parse(value).unwrap());
    let values: Vec<_> = parameters.iter().map(analytic_row).collect();
    let nodes: Vec<_> = parameters.iter().zip(&values).map(|(tau, row)| {
        NumericalSpineNode {
            tau: tau.clone(),
            parameters: [
                [row[12].clone(), row[13].clone()],
                [row[14].clone(), row[15].clone()],
            ],
            charts: [0, 0],
            centre: [row[0].clone(), row[1].clone(), row[2].clone()],
        }
    }).collect();
    let mut dec = Decider::default();
    let fitted = candidate(&nodes, &values, &mut dec).unwrap();
    let basis = validation::Basis::new();
    let bound = Rational::one().mul_pow2(-80);
    for (index, pair) in nodes.windows(2).enumerate() {
        let span = validation::Span::new(
            fitted.data[3 * index..3 * index + 4].try_into().unwrap(),
        );
        for (quarter, share) in basis.shares.iter().enumerate() {
            let tau = &pair[0].tau + &(&(&pair[1].tau - &pair[0].tau) * share);
            let truth = analytic_row(&tau);
            let actual = span.at(quarter, &basis).unwrap();
            for (j, across) in basis.section_shares.iter().enumerate() {
                let expected = section(&truth, across);
                assert!(distance(&actual.points[j], &expected) < bound);
                let source = fitted.surface
                    .evaluate_box([&tau, &tau], [across, across], 128, &dec)
                    .unwrap();
                assert!(distance(&source, &expected) < bound);
            }
            for side in 0..2 {
                let expected = &truth[4 * side * 2..4 * side * 2 + 3];
                let expected: [Rational; 3] = std::array::from_fn(|i| expected[i].clone());
                let rail = fitted.rails[side].evaluate_box([&tau, &tau], 128, &dec).unwrap();
                assert!(distance(&rail, &expected) < bound);
                for axis in 0..2 {
                    assert!(
                        (&actual.parameters[side][axis] - &truth[12 + 2 * side + axis]).abs() < bound
                    );
                }
            }
        }
    }
}
