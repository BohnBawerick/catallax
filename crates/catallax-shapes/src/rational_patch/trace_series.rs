//! Formal series shared by the certified implicit graph producers.

use catallax_numbers::{Ball, Decider, Rational, RealPoly, Refusal};

fn nonzero_len(coefficients: &[Ball]) -> usize {
    coefficients
        .iter()
        .rposition(|ball| !ball.is_exact_zero())
        .map_or(0, |index| index + 1)
}

fn product(
    a: &[Ball],
    b: &[Ball],
    order: usize,
    zero: &Ball,
    dec: &Decider,
    out: &mut Vec<Ball>,
) -> Result<(), Refusal> {
    // The free chart coordinate is affine, and the implicit recurrence
    // supplies its not-yet-solved coefficients as exact zero.
    let a = &a[..nonzero_len(a)];
    let b = &b[..nonzero_len(&b[..b.len().min(order + 1)])];
    out.clear();
    if a.is_empty() || b.is_empty() {
        return Ok(());
    }
    out.resize((a.len() + b.len() - 1).min(order + 1), zero.clone());
    for (i, a) in a.iter().enumerate() {
        if a.is_exact_zero() {
            continue;
        }
        for (j, b) in b.iter().enumerate().take(order + 1 - i) {
            if !b.is_exact_zero() {
                dec.check_deadline()?;
                out[i + j] = &out[i + j] + &(a * b);
            }
        }
    }
    Ok(())
}

pub(super) fn series(
    poly: &RealPoly,
    at: &[Vec<Ball>],
    order: usize,
    bits: u32,
    dec: &Decider,
) -> Result<Vec<Ball>, Refusal> {
    let zero = Ball::from_rational(&Rational::zero(), bits);
    let mut out = vec![zero.clone(); order + 1];
    let mut term = Vec::with_capacity(order + 1);
    let mut scratch = Vec::with_capacity(order + 1);
    for (monomial, coefficient) in poly.terms() {
        dec.check_deadline()?;
        term.clear();
        term.push(dec.enclose_at(&coefficient, bits)?);
        for (axis, &power) in monomial.exponents().iter().enumerate() {
            for _ in 0..power {
                product(&term, &at[axis], order, &zero, dec, &mut scratch)?;
                std::mem::swap(&mut term, &mut scratch);
            }
        }
        for (sum, coefficient) in out.iter_mut().zip(&term) {
            if !coefficient.is_exact_zero() {
                *sum = &*sum + coefficient;
            }
        }
    }
    Ok(out)
}
