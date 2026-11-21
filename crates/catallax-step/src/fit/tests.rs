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

fn work(run: impl FnOnce()) -> catallax_numbers::OperationWork {
    catallax_numbers::with_operation_budget(1_000_000, run).1
}

// Limbs are 64 bits. A rational counts its numerator's plus its denominator's.
// Every count adds the unit a charge pays for its deadline check.
#[test]
fn exact_arithmetic_is_counted_at_the_deciders_limb_rate_and_not_charged() {
    let third = Rational::from_ratio(1, 3).unwrap();
    let wide = Rational::from_i64(1).mul_pow2(64);
    let counted = |run: &dyn Fn()| {
        let work = work(run);
        assert_eq!(work.used, 0);
        assert!(!work.exhausted);
        work.uncharged
    };
    assert_eq!(counted(&|| drop(counted::mul(&third, &wide))), 2 * 3 + 1);
    assert_eq!(counted(&|| drop(counted::div(&wide, &third))), 3 * 2 + 1);
    assert_eq!(counted(&|| drop(counted::add(&third, &wide))), 2 + 3 + 1);
    assert_eq!(counted(&|| drop(counted::sub(&wide, &wide))), 3 + 3 + 1);
    // Three products: 1 by 1/3, 1/3 by 1/3 and 1/9 by 1/3, two limbs each side.
    assert_eq!(counted(&|| drop(counted::pow(&third, 3))), 3 * (2 * 2 + 1));
    // The midpoint sum, the value times itself for its decimal conversion, the
    // two differences from the neighbouring decimals 333/1000 and 334/1000,
    // and the two differences of the written value from the ends.
    assert_eq!(
        counted(&|| drop(counted::rounded(&third, &third, 3))),
        (2 + 2 + 1) + (2 * 2 + 1) + 4 * (2 + 2 + 1)
    );
    // The same value as the uncounted rounding.
    let (plain, priced) = (
        number::rounded(&third, &wide, 17),
        counted::rounded(&third, &wide, 17),
    );
    assert_eq!(plain.text, priced.text);
    assert_eq!(plain.value, priced.value);
    assert_eq!(plain.error, priced.error);
}

fn enclosures<const N: usize>(precision: u32) -> [[Ball; N]; DEGREE + 1] {
    std::array::from_fn(|k| {
        std::array::from_fn(|i| {
            let lo = Rational::from_ratio(k as i64 * 7 + i as i64 * 3 - 18, 11).unwrap();
            let hi = &lo + &Rational::from_ratio((k + i + 1) as i64, 37).unwrap();
            Ball::from_rational(&lo, precision).hull(&Ball::from_rational(&hi, precision))
        })
    })
}

#[test]
fn a_fit_on_wider_numbers_counts_more() {
    let width = Rational::from_ratio(3, 7).unwrap();
    let fit = |precision| {
        let coefficients = enclosures::<3>(precision);
        let work = work(|| taylor_controls(&coefficients, &width, 17, |_, _| {}));
        assert_eq!(work.used, 0);
        work.uncharged
    };
    let (narrow, wide) = (fit(64), fit(256));
    // Four times the limbs in each enclosure end. Most products have one narrow
    // operand, a weight or a power of the width, so the count grows about with
    // the limbs. No term grows faster than their square.
    assert!(narrow < wide && wide < 16 * narrow, "{narrow} then {wide}");
}

// The count is no cap: a fit far beyond the allowance in force writes every
// control and leaves the allowance whole.
#[test]
fn a_fit_beyond_the_allowance_still_writes_every_control() {
    let coefficients = enclosures::<3>(192);
    let width = Rational::from_ratio(3, 7).unwrap();
    let mut written = 0;
    let ((), work) = catallax_numbers::with_operation_budget(100, || {
        taylor_controls(&coefficients, &width, 17, |_, _| written += 1);
    });
    assert_eq!(written, DEGREE + 1);
    assert_eq!(work.used, 0);
    assert!(!work.exhausted);
    assert!(work.uncharged > 100, "{}", work.uncharged);
}

#[test]
fn the_power_basis_of_a_span_is_counted() {
    let controls: [[Rational; 2]; DEGREE + 1] = std::array::from_fn(|k| {
        std::array::from_fn(|i| Rational::from_ratio(k as i64 * 5 - 9, i as i64 + 3).unwrap())
    });
    let mut powers = None;
    let work = work(|| powers = Some(bezier_power(controls.clone())));
    // One product for each power and coordinate, and one difference for each
    // entry of the difference table below its first row. Every operand here
    // has two limbs.
    let products = (DEGREE as u64 + 1) * 2;
    let differences = (DEGREE as u64 * (DEGREE as u64 + 1) / 2) * 2;
    assert_eq!(
        work.uncharged,
        products * (2 * 2 + 1) + differences * (2 + 2 + 1)
    );
    assert_eq!(work.used, 0);
    // The leading power is the first control and the next is seven times the
    // first difference.
    let powers = powers.unwrap();
    assert_eq!(powers[0], controls[0]);
    assert_eq!(
        powers[1][0],
        (&controls[1][0] - &controls[0][0]) * Rational::from_i64(7)
    );
}

// A knot after the first is one sum and one writing in decimal. The integers
// here have two limbs, so the sum counts 2 + 2 and the conversion 2 x 2, each
// with the unit of a charge's deadline check. The knots are the ones the fit
// wrote before it counted them.
#[test]
fn each_knot_counts_its_sum_and_its_decimal_conversion() {
    for (digits, last) in [(1, "300"), (3, "21"), (17, "21")] {
        let mut made = Vec::new();
        let work = work(|| made = knots(21, digits));
        assert_eq!(made.len(), 22);
        assert_eq!(made[0], Rational::zero());
        assert_eq!(made[9], Rational::from_i64(9));
        // At one digit each knot past 10 is the next one-digit decimal: 20,
        // 30 and on through 100 and 200 to 300.
        assert_eq!(made[21], Rational::parse(last).unwrap(), "{digits} digits");
        assert_eq!(work.uncharged, 21 * ((2 + 2 + 1) + (2 * 2 + 1)));
        assert_eq!(work.used, 0);
    }
}
