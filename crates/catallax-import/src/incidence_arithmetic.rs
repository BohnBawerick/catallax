//! Charged exact arithmetic for import incidence certificates.

use crate::Vec3;
use catallax_numbers::{Decider, Integer, Rational, Refusal, Sign};
use std::cmp::Ordering;

pub(super) fn failure(error: Refusal) -> String {
    format!(
        "exact incidence arithmetic: {}: {}",
        error.code, error.detail
    )
}

pub(super) fn checkpoint(dec: &mut Decider) -> Result<(), String> {
    dec.check_deadline().map_err(failure)
}

// Reserve the numeric layer's limb-based cost before each rational operation.
pub(super) fn arithmetic(
    a: &Rational,
    b: &Rational,
    op: char,
    dec: &mut Decider,
) -> Result<Rational, String> {
    dec.reserve_exact_arithmetic(a, b, matches!(op, '*' | '/'))
        .map_err(failure)?;
    Ok(match op {
        '+' => a + b,
        '-' => a - b,
        '*' => a * b,
        '/' => a / b,
        _ => unreachable!("exact incidence arithmetic"),
    })
}

pub(super) fn compare(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Ordering, String> {
    // Fraction ordering can cross-multiply twice, as in the patch proofs.
    dec.reserve_exact_arithmetic(a, b, true).map_err(failure)?;
    dec.reserve_exact_arithmetic(a, b, true).map_err(failure)?;
    Ok(a.cmp(b))
}

pub(super) fn same(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<bool, String> {
    for i in 0..3 {
        dec.reserve_exact_arithmetic(&a[i], &b[i], false)
            .map_err(failure)?;
        if a[i] != b[i] {
            return Ok(false);
        }
    }
    Ok(true)
}

pub(super) fn minimum(a: Rational, b: Rational, dec: &mut Decider) -> Result<Rational, String> {
    Ok(if compare(&a, &b, dec)?.is_gt() { b } else { a })
}

pub(super) fn sub(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<Vec3, String> {
    Ok([
        arithmetic(&a[0], &b[0], '-', dec)?,
        arithmetic(&a[1], &b[1], '-', dec)?,
        arithmetic(&a[2], &b[2], '-', dec)?,
    ])
}

pub(super) fn dot(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<Rational, String> {
    let mut result = Rational::zero();
    for i in 0..3 {
        let product = arithmetic(&a[i], &b[i], '*', dec)?;
        result = arithmetic(&result, &product, '+', dec)?;
    }
    Ok(result)
}

pub(super) fn determinant(
    a: &Rational,
    b: &Rational,
    c: &Rational,
    d: &Rational,
    dec: &mut Decider,
) -> Result<Rational, String> {
    let ab = arithmetic(a, b, '*', dec)?;
    let cd = arithmetic(c, d, '*', dec)?;
    arithmetic(&ab, &cd, '-', dec)
}

pub(super) fn cross(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<Vec3, String> {
    Ok([
        determinant(&a[1], &b[2], &a[2], &b[1], dec)?,
        determinant(&a[2], &b[0], &a[0], &b[2], dec)?,
        determinant(&a[0], &b[1], &a[1], &b[0], dec)?,
    ])
}

pub(super) fn maximum(a: Rational, b: Rational, dec: &mut Decider) -> Result<Rational, String> {
    Ok(if compare(&a, &b, dec)?.is_lt() { b } else { a })
}

pub(super) fn equal(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<bool, String> {
    dec.reserve_exact_arithmetic(a, b, false).map_err(failure)?;
    Ok(a == b)
}

pub(super) fn same_sequence(
    a: &[Rational],
    b: &[Rational],
    dec: &mut Decider,
) -> Result<bool, String> {
    if a.len() != b.len() {
        return Ok(false);
    }
    for (x, y) in a.iter().zip(b) {
        checkpoint(dec)?;
        if !equal(x, y, dec)? {
            return Ok(false);
        }
    }
    Ok(true)
}

pub(super) fn integer_arithmetic(
    a: &Integer,
    b: &Integer,
    op: char,
    dec: &mut Decider,
) -> Result<Integer, String> {
    let ra = Rational::from_integer(a.clone());
    let rb = Rational::from_integer(b.clone());
    dec.reserve_exact_arithmetic(&ra, &rb, matches!(op, '*' | '/'))
        .map_err(failure)?;
    Ok(match op {
        '+' => a + b,
        '*' => a * b,
        '/' => {
            a.div_rem_floor(b)
                .ok_or("exact incidence division by zero")?
                .0
        }
        _ => unreachable!("exact incidence integer arithmetic"),
    })
}

// Preserve the import policy's existing 1e-12 mm outward enclosure.
pub(super) fn sqrt_enclosure(
    q: &Rational,
    dec: &mut Decider,
) -> Result<(Rational, Rational), String> {
    checkpoint(dec)?;
    let scale = Rational::from_i64(1_000_000_000_000);
    let scaled = arithmetic(q, &scale, '*', dec)?;
    let scaled = arithmetic(&scaled, &scale, '*', dec)?;
    let n = integer_arithmetic(scaled.numerator(), scaled.denominator(), '/', dec)?;
    let k = if n.sign() != Sign::Positive {
        Integer::zero()
    } else {
        let mut x = Integer::one().mul_pow2(n.bit_length().div_ceil(2));
        loop {
            checkpoint(dec)?;
            let quotient = integer_arithmetic(&n, &x, '/', dec)?;
            let sum = integer_arithmetic(&x, &quotient, '+', dec)?;
            let y = integer_arithmetic(&sum, &Integer::from_i64(2), '/', dec)?;
            let ry = Rational::from_integer(y.clone());
            let rx = Rational::from_integer(x.clone());
            dec.reserve_exact_arithmetic(&ry, &rx, false)
                .map_err(failure)?;
            if y >= x {
                break x;
            }
            x = y;
        }
    };
    let lo = arithmetic(&Rational::from_integer(k.clone()), &scale, '/', dec)?;
    let square = integer_arithmetic(&k, &k, '*', dec)?;
    let exact = equal(&Rational::from_integer(square), &scaled, dec)?;
    let hi = if exact {
        lo.clone()
    } else {
        let next = integer_arithmetic(&k, &Integer::one(), '+', dec)?;
        arithmetic(&Rational::from_integer(next), &scale, '/', dec)?
    };
    Ok((lo, hi))
}

pub(super) fn binomial(n: usize, k: usize, dec: &mut Decider) -> Result<Rational, String> {
    let mut value = Rational::one();
    for i in 0..k.min(n - k) {
        checkpoint(dec)?;
        value = arithmetic(&value, &Rational::from_i64((n - i) as i64), '*', dec)?;
        value = arithmetic(&value, &Rational::from_i64((i + 1) as i64), '/', dec)?;
    }
    Ok(value)
}

pub(super) fn vector_compare(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<Ordering, String> {
    for i in 0..3 {
        let order = compare(&a[i], &b[i], dec)?;
        if order != Ordering::Equal {
            return Ok(order);
        }
    }
    Ok(Ordering::Equal)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    #[test]
    fn charged_roots_preserve_the_existing_decimal_enclosures() {
        for text in [
            "0",
            "1/4",
            "2",
            "1/1000000000000000000000000000000",
            "123456789123456789",
        ] {
            let q = Rational::parse(text).unwrap();
            let (result, work) =
                with_operation_budget(20_000_000, || sqrt_enclosure(&q, &mut Decider::default()));
            assert_eq!(result.unwrap(), crate::gap::sqrt_enclosure(&q));
            assert!(work.used > 0);
            assert!(!work.exhausted);
        }
        let (result, work) = with_operation_budget(0, || {
            sqrt_enclosure(&Rational::one(), &mut Decider::default())
        });
        assert!(result.is_err());
        assert!(work.exhausted);
    }
}
