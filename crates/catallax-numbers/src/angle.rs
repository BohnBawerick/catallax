//! Rational-degree angles selected by exact Chebyshev root order.

use crate::{Decider, Integer, Rational, Real, Refusal, UPoly, algebraic, decider::Work};

/// Exact sine and cosine of a rational angle in degrees, in that order.
///
/// Reduction and sector selection use integer arithmetic. Non-quadratic values
/// retain a certified root of a rational Chebyshev polynomial. Construction and
/// isolation share the caller's deterministic per-decision work allowance.
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
        let degree = turn
            .denominator()
            .to_i64()
            .and_then(|n| usize::try_from(n).ok())
            .ok_or("numbers.budget-exhausted")?;
        let numerator = turn
            .numerator()
            .to_i64()
            .and_then(|n| usize::try_from(n).ok())
            .ok_or("numbers.budget-exhausted")?;
        // Reserve before allocating or visiting a degree-sized polynomial.
        work.charge((degree as u64).saturating_mul(degree as u64))?;
        let two_x = UPoly::new(vec![Rational::zero(), Rational::from_i64(2)]);
        let mut previous = UPoly::constant(Rational::one());
        let mut current = UPoly::x();
        for _ in 1..degree {
            work.poly(&two_x, &current)?;
            let product = &two_x * &current;
            work.poly(&product, &previous)?;
            let next = &product - &previous;
            previous = current;
            current = next;
        }
        let parity = numerator % 2;
        let constant = UPoly::constant(Rational::from_i64(if parity == 0 { 1 } else { -1 }));
        work.poly(&current, &constant)?;
        let polynomial = &current - &constant;
        let roots = algebraic::isolate(&polynomial, work)?;
        // T_n(cos(k*pi/n)) = (-1)^k. Cosine decreases on [0, pi],
        // so the k with the chosen parity run in reverse root order.
        let largest = degree - (degree % 2 != parity) as usize;
        let index = (largest - numerator) / 2;
        if roots.len() != (largest - parity) / 2 + 1 {
            return Err("numbers.root-not-certified");
        }
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
