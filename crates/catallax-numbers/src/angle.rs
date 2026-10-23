//! Rational-degree angles selected by exact real cyclotomic root order.

use crate::{Decider, Integer, Rational, Real, Refusal, UPoly, algebraic, decider::Work};

/// Exact sine and cosine of a rational angle in degrees, in that order.
///
/// Reduction and sector selection use integer arithmetic. Non-quadratic values
/// retain a certified root of the real cyclotomic polynomial for the primitive
/// angle, rather than a larger Chebyshev equation containing unrelated roots.
/// Construction and isolation share the caller's deterministic work allowance.
pub fn sin_cos_degrees(angle: &Rational, dec: &mut Decider) -> Result<[Real; 2], Refusal> {
    dec.angle_degrees(angle)
}

pub(crate) fn construct(angle: &Rational, work: &mut Work) -> Result<[Real; 2], &'static str> {
    work.exact(angle, &Rational::from_i64(360), true)?;
    let period = angle.denominator() * &Integer::from_i64(360);
    let (_, numerator) = angle
        .numerator()
        .div_rem_floor(&period)
        .expect("positive period");
    let mut acute = Rational::new(numerator, angle.denominator().clone()).unwrap();
    let negative_sine = acute > Rational::from_i64(180);
    if negative_sine {
        acute = &Rational::from_i64(360) - &acute;
    }
    let negative_cosine = acute > Rational::from_i64(90);
    if negative_cosine {
        acute = &Rational::from_i64(180) - &acute;
    }
    work.charge(1)?;
    let half = Real::from(Rational::from_ratio(1, 2).unwrap());
    let [sine, cosine] = if acute.is_zero() {
        [Real::zero(), Real::one()]
    } else if acute == Rational::from_i64(90) {
        [Real::one(), Real::zero()]
    } else if acute == Rational::from_i64(30) {
        [half.clone(), &Real::from_i64(3).sqrt() * &half]
    } else if acute == Rational::from_i64(45) {
        let value = &Real::from_i64(2).sqrt() * &half;
        [value.clone(), value]
    } else if acute == Rational::from_i64(60) {
        [&Real::from_i64(3).sqrt() * &half, half]
    } else if acute == Rational::from_ratio(45, 2).unwrap() {
        let root_two = Real::from_i64(2).sqrt();
        [
            (&Real::from_i64(2) - &root_two).sqrt() * half.clone(),
            (&Real::from_i64(2) + &root_two).sqrt() * half,
        ]
    } else if acute == Rational::from_ratio(135, 2).unwrap() {
        let root_two = Real::from_i64(2).sqrt();
        [
            (&Real::from_i64(2) + &root_two).sqrt() * half.clone(),
            (&Real::from_i64(2) - &root_two).sqrt() * half,
        ]
    } else if acute == Rational::from_i64(15) || acute == Rational::from_i64(75) {
        let root_six = Real::from_i64(6).sqrt();
        let root_two = Real::from_i64(2).sqrt();
        let quarter = Real::from(Rational::from_ratio(1, 4).unwrap());
        let small = (&root_six - &root_two) * quarter.clone();
        let big = (&root_six + &root_two) * quarter;
        if acute == Rational::from_i64(15) {
            [small, big]
        } else {
            [big, small]
        }
    } else {
        let turn = &acute / &Rational::from_i64(180);
        let denominator = turn
            .denominator()
            .to_i64()
            .and_then(|n| usize::try_from(n).ok())
            .ok_or("numbers.budget-exhausted")?;
        let numerator = turn
            .numerator()
            .to_i64()
            .and_then(|n| usize::try_from(n).ok())
            .ok_or("numbers.budget-exhausted")?;
        let (order, numerator) = if numerator % 2 == 0 {
            (denominator, numerator / 2)
        } else {
            (
                denominator
                    .checked_mul(2)
                    .ok_or("numbers.budget-exhausted")?,
                numerator,
            )
        };
        // Reserve before allocating degree-sized polynomials or visiting indices.
        work.charge((order as u64).saturating_mul(order as u64))?;
        let polynomial = real_cyclotomic(order, work)?;
        let roots = algebraic::isolate(&polynomial, work)?;
        let mut count = 0;
        let mut index = 0;
        for j in 1..=order / 2 {
            if coprime(j, order) {
                count += 1;
                if j > numerator {
                    index += 1;
                }
            }
        }
        if roots.len() != count || index >= count {
            return Err("numbers.root-not-certified");
        }
        // Primitive j in (0,n/2) give all conjugates. Cosine decreases with j,
        // so the chosen primitive numerator has this exact ascending root index.
        let cosine = Real::from_algebraic(roots[index].value.clone());
        let sine = (&Real::one() - &cosine.square()).sqrt();
        [sine, cosine]
    };
    work.check_deadline()?;
    Ok([
        if negative_sine { -sine } else { sine },
        if negative_cosine { -cosine } else { cosine },
    ])
}

fn coprime(mut a: usize, mut b: usize) -> bool {
    while b != 0 {
        let remainder = a % b;
        a = b;
        b = remainder;
    }
    a == 1
}
fn mobius(mut n: usize) -> i64 {
    let mut sign = 1;
    let mut p = 2;
    while p <= n / p {
        if n % p == 0 {
            n /= p;
            if n % p == 0 {
                return 0;
            }
            sign = -sign;
        }
        p += 1;
    }
    if n > 1 { -sign } else { sign }
}
fn real_cyclotomic(order: usize, work: &mut Work) -> Result<UPoly, &'static str> {
    // Phi_n(z) = product_(d|n) (z^d-1)^mu(n/d).
    let mut numerator = UPoly::constant(Rational::one());
    let mut denominator = numerator.clone();
    for d in 1..=order {
        if order % d != 0 {
            continue;
        }
        let power = mobius(order / d);
        if power == 0 {
            continue;
        }
        work.charge(d as u64 + 1)?;
        let mut coefficients = vec![Rational::zero(); d + 1];
        coefficients[0] = -Rational::one();
        coefficients[d] = Rational::one();
        let factor = UPoly::new(coefficients);
        let product = if power > 0 {
            &mut numerator
        } else {
            &mut denominator
        };
        work.poly(product, &factor)?;
        *product = &*product * &factor;
    }
    work.poly(&numerator, &denominator)?;
    let (cyclotomic, remainder) = numerator
        .div_rem_checked(&denominator, &|| work.check_deadline())?
        .ok_or("numbers.root-not-certified")?;
    let degree = cyclotomic.degree().ok_or("numbers.root-not-certified")?;
    if !remainder.is_zero() || degree % 2 != 0 {
        return Err("numbers.root-not-certified");
    }
    let half = degree / 2;
    let coefficients = cyclotomic.coefficients();
    // Phi_n(z)/z^half is a_k + sum a_(k+j)(z^j+z^-j).
    // D_0=2, D_1=2x, D_(j+1)=2x D_j-D_(j-1) therefore gives
    // its polynomial in x=cos(2*pi*j/n), of degree phi(n)/2.
    let two_x = UPoly::new(vec![Rational::zero(), Rational::from_i64(2)]);
    let mut previous = UPoly::constant(Rational::from_i64(2));
    let mut current = two_x.clone();
    let mut result = UPoly::constant(coefficients[half].clone());
    for j in 1..=half {
        let coefficient = &coefficients[half + j];
        if !coefficient.is_zero() {
            work.poly(&current, &UPoly::constant(coefficient.clone()))?;
            let term = current.scale(coefficient);
            work.poly(&result, &term)?;
            result = &result + &term;
        }
        if j < half {
            work.poly(&two_x, &current)?;
            let product = &two_x * &current;
            work.poly(&product, &previous)?;
            let next = &product - &previous;
            previous = current;
            current = next;
        }
    }
    Ok(result)
}
