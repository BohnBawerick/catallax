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
