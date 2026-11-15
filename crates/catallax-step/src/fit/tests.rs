use super::*;

type Written<const N: usize> = ([number::Number; N], [[Rational; 2]; N]);

fn matrix<const N: usize>(
    coefficients: &[[Ball; N]],
    width: &Rational,
    digits: u32,
) -> Vec<Written<N>> {
    (0..=DEGREE)
        .map(|j| {
            let mut control: [[Rational; 2]; N] =
                std::array::from_fn(|_| [Rational::zero(), Rational::zero()]);
            for (k, coefficient) in coefficients.iter().enumerate() {
                let weight = (0..=k.min(j)).fold(Rational::zero(), |sum, l| {
                    sum + Rational::from_ratio(binomial(k, l) * binomial(j, l), binomial(DEGREE, l))
                        .unwrap()
                        * Rational::from_ratio(-1, 2)
                            .unwrap()
                            .pow((k - l) as i32)
                            .unwrap()
                }) * width.pow(k as i32).unwrap();
                for i in 0..N {
                    let a = &coefficient[i].lo() * &weight;
                    let b = &coefficient[i].hi() * &weight;
                    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
                    control[i][0] = &control[i][0] + &lo;
                    control[i][1] = &control[i][1] + &hi;
                }
            }
            (
                control
                    .each_ref()
                    .map(|b| number::rounded(&b[0], &b[1], digits)),
                control,
            )
        })
        .collect()
}

fn compare<const N: usize>() {
    for seed in 0..6 {
        let coefficients: [[Ball; N]; DEGREE + 1] = std::array::from_fn(|k| {
            std::array::from_fn(|i| {
                let lo =
                    Rational::from_ratio(seed * 11 + k as i64 * 7 + i as i64 * 3 - 29, 11).unwrap();
                let hi = &lo + &Rational::from_ratio((k + i + 1) as i64, 37).unwrap();
                if seed == 0 {
                    Ball::from_rational(&Rational::zero(), 192)
                } else {
                    Ball::from_rational(&lo, 192).hull(&Ball::from_rational(&hi, 192))
                }
            })
        });
        for width in ["0", "1", "-1", "3/7", "1e-20", "1e30"] {
            let width = Rational::parse(width).unwrap();
            for digits in [1, 17, 32] {
                let mut actual = Vec::new();
                taylor_controls(&coefficients, &width, digits, |point, enclosure| {
                    actual.push((point, enclosure));
                });
                let expected = matrix(&coefficients, &width, digits);
                assert_eq!(actual.len(), expected.len());
                for ((actual, bounds), (expected, expected_bounds)) in actual.iter().zip(&expected)
                {
                    assert_eq!(bounds, expected_bounds);
                    for (a, b) in actual.iter().zip(expected) {
                        assert_eq!(a.text, b.text);
                        assert_eq!(a.value, b.value);
                        assert_eq!(a.error, b.error);
                    }
                }
            }
        }
    }
}

#[test]
fn reflected_controls_match_exact_matrix_enclosures_and_written_controls() {
    compare::<2>();
    compare::<3>();
}
