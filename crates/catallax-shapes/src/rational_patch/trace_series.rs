//! Formal series shared by the certified implicit graph producers.

use catallax_numbers::{Ball, Decider, Rational, RealPoly, Refusal};

fn product(
    a: &[Ball],
    b: &[Ball],
    order: usize,
    bits: u32,
    dec: &Decider,
) -> Result<Vec<Ball>, Refusal> {
    let mut out = vec![Ball::from_rational(&Rational::zero(), bits); order + 1];
    for (i, a) in a.iter().enumerate().take(order + 1) {
        for (j, b) in b.iter().enumerate().take(order + 1 - i) {
            dec.check_deadline()?;
            out[i + j] = &out[i + j] + &(a * b);
        }
    }
    Ok(out)
}

pub(super) fn series(
    poly: &RealPoly,
    at: &[Vec<Ball>],
    order: usize,
    bits: u32,
    dec: &Decider,
) -> Result<Vec<Ball>, Refusal> {
    let zero = || Ball::from_rational(&Rational::zero(), bits);
    let mut out = vec![zero(); order + 1];
    for (monomial, coefficient) in poly.terms() {
        let mut term = vec![zero(); order + 1];
        term[0] = dec.enclose_at(&coefficient, bits)?;
        for (axis, &power) in monomial.exponents().iter().enumerate() {
            for _ in 0..power {
                term = product(&term, &at[axis], order, bits, dec)?;
            }
        }
        for (sum, coefficient) in out.iter_mut().zip(term) {
            *sum = &*sum + &coefficient;
        }
    }
    Ok(out)
}
