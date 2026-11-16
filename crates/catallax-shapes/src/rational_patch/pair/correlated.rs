//! Field bounds under a proved common affine projection of polynomial graphs.

use super::*;

fn quotient(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, true)?;
    Ok(a / b)
}

fn scale(p: &MPoly, coefficient: Rational, dec: &mut Decider) -> Result<MPoly, Refusal> {
    product(p, &MPoly::constant(2, coefficient), dec)
}

fn sum(a: &MPoly, b: &MPoly, dec: &mut Decider) -> Result<MPoly, Refusal> {
    difference(a, &scale(b, -Rational::one(), dec)?, dec)
}

fn constant(p: &MPoly, dec: &mut Decider) -> Result<Option<Rational>, Refusal> {
    dec.charge_work(1)?;
    if p.total_degree() != Some(0) {
        return Ok(None);
    }
    Ok(Some(p.terms()[0].1.clone()))
}

fn affine(p: &MPoly, dec: &mut Decider) -> Result<Option<[Rational; 3]>, Refusal> {
    dec.charge_work(1)?;
    if p.total_degree().is_some_and(|d| d > 1) {
        return Ok(None);
    }
    let mut coefficients = std::array::from_fn(|_| Rational::zero());
    for (m, c) in p.terms() {
        dec.charge_work(1)?;
        let index = m
            .exponents()
            .iter()
            .position(|e| *e == 1)
            .map_or(0, |i| i + 1);
        coefficients[index] = c;
    }
    Ok(Some(coefficients))
}

fn unturned(
    h: &[MPoly; 4],
    rotation: &[[Rational; 3]; 3],
    dec: &mut Decider,
) -> Result<[MPoly; 3], Refusal> {
    let mut coordinates = Vec::new();
    for axis in 0..3 {
        let mut value = MPoly::zero(2);
        for row in 0..3 {
            dec.charge_work(1)?;
            value = sum(
                &value,
                &scale(&h[row], rotation[row][axis].clone(), dec)?,
                dec,
            )?;
        }
        coordinates.push(value);
    }
    // The integer quaternion rotation has R^T R = c I with c > 0.
    // This common scale cancels when recovering the second parameters.
    Ok(coordinates.try_into().expect("three projected coordinates"))
}

fn compose(p: &MPoly, map: &[MPoly; 2], dec: &mut Decider) -> Result<MPoly, Refusal> {
    let mut result = MPoly::zero(2);
    for (m, c) in p.terms() {
        dec.charge_work(1)?;
        let mut term = MPoly::constant(2, c);
        for (axis, &power) in m.exponents().iter().enumerate() {
            for _ in 0..power {
                dec.charge_work(1)?;
                term = product(&term, &map[axis], dec)?;
            }
        }
        result = sum(&result, &term, dec)?;
    }
    Ok(result)
}

fn normalize(p: &MPoly, dec: &mut Decider) -> Result<MPoly, Refusal> {
    let mut maximum = Rational::zero();
    for (_, coefficient) in p.terms() {
        dec.charge_work(1)?;
        let value = coefficient.abs();
        if net::less(&maximum, &value, dec)? {
            maximum = value;
        }
    }
    if maximum.is_zero() {
        return Ok(p.clone());
    }
    scale(p, quotient(&Rational::one(), &maximum, dec)?, dec)
}

fn fields(
    first: &SpanPolynomials,
    second: &SpanPolynomials,
    rotation: &[[Rational; 3]; 3],
    dec: &mut Decider,
) -> Result<Option<[MPoly; 2]>, Refusal> {
    for source in [first, second] {
        for p in &source.homogeneous[..3] {
            dec.charge_work(1)?;
            if p.total_degree().is_some_and(|d| d > 6) {
                return Ok(None);
            }
        }
        if source.homogeneous[..3]
            .iter()
            .all(|p| p.total_degree().is_none_or(|d| d <= 1))
        {
            return Ok(None);
        }
    }
    let Some(wa) = constant(&first.homogeneous[3], dec)? else {
        return Ok(None);
    };
    let Some(wb) = constant(&second.homogeneous[3], dec)? else {
        return Ok(None);
    };
    if wa.sign() != Sign::Positive || wb.sign() != Sign::Positive {
        return Ok(None);
    }
    let a = unturned(&first.homogeneous, rotation, dec)?;
    let b = unturned(&second.homogeneous, rotation, dec)?;
    for [x, y] in [[0, 1], [0, 2], [1, 2]] {
        if affine(&a[x], dec)?.is_none() || affine(&a[y], dec)?.is_none() {
            continue;
        }
        let (Some(bx), Some(by)) = (affine(&b[x], dec)?, affine(&b[y], dec)?) else {
            continue;
        };
        let determinant = net::sub(
            &net::mul(&bx[1], &by[2], dec)?,
            &net::mul(&bx[2], &by[1], dec)?,
            dec,
        )?;
        if determinant.is_zero() {
            continue;
        }
        let weight = quotient(&wb, &wa, dec)?;
        let rhsx = difference(
            &scale(&a[x], weight.clone(), dec)?,
            &MPoly::constant(2, bx[0].clone()),
            dec,
        )?;
        let rhsy = difference(
            &scale(&a[y], weight, dec)?,
            &MPoly::constant(2, by[0].clone()),
            dec,
        )?;
        let inverse = quotient(&Rational::one(), &determinant, dec)?;
        let map = [
            scale(
                &difference(
                    &scale(&rhsx, by[2].clone(), dec)?,
                    &scale(&rhsy, bx[2].clone(), dec)?,
                    dec,
                )?,
                inverse.clone(),
                dec,
            )?,
            scale(
                &difference(
                    &scale(&rhsy, bx[1].clone(), dec)?,
                    &scale(&rhsx, by[1].clone(), dec)?,
                    dec,
                )?,
                inverse,
                dec,
            )?,
        ];
        let normal = [
            compose(&second.normal[0], &map, dec)?,
            compose(&second.normal[1], &map, dec)?,
            compose(&second.normal[2], &map, dec)?,
        ];
        let mut result = Vec::new();
        for tangent in &first.tangents {
            let mut field = MPoly::zero(2);
            for (t, n) in tangent.iter().zip(&normal) {
                dec.charge_work(1)?;
                field = sum(&field, &product(t, n, dec)?, dec)?;
            }
            // A positive constant scale preserves the directional certificate.
            result.push(normalize(&field, dec)?);
        }
        return Ok(Some(
            result.try_into().expect("two native field components"),
        ));
    }
    Ok(None)
}

pub(super) fn prepare(
    ctx: &Context,
    spans: [usize; 2],
    dec: &mut Decider,
) -> Result<Option<[MPoly; 2]>, Refusal> {
    if ctx.lengths() != 0 {
        return Ok(None);
    }
    fields(
        &ctx.sides[0].spans[spans[0]],
        &ctx.sides[1].spans[spans[1]],
        &ctx.rotation,
        dec,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn affine_correlation_keeps_near_crossings_and_never_certifies_coincidence() {
        let mut dec = Decider::default();
        let u = MPoly::var(2, 0);
        let v = MPoly::var(2, 1);
        let height = MPoly::from_terms(
            2,
            vec![
                (vec![0, 0], Rational::from_i64(2)),
                (vec![0, 1], Rational::one()),
            ],
        );
        let roof = &(&u - &(&u * &u)) * &height;
        let a = [u.clone(), roof, v, MPoly::constant(2, Rational::one())];
        let rotation = quaternion([1, 1, 0, 0]).map(|row| row.map(Rational::from_i64));
        for gap in [
            Rational::zero(),
            Rational::from_ratio(1, 1000).unwrap(),
            Rational::from_ratio(-1, 1000).unwrap(),
        ] {
            let mut b = a.clone();
            b[0] = &b[0] + &MPoly::constant(2, gap.clone());
            let first = SpanPolynomials::new(&a, &rotation, &mut dec).unwrap();
            let second = SpanPolynomials::new(&b, &rotation, &mut dec).unwrap();
            let field = fields(&first, &second, &rotation, &mut dec)
                .unwrap()
                .unwrap();
            let domain = [
                [Rational::zero(), Rational::one()],
                [Rational::zero(), Rational::one()],
            ];
            let range = net::Dense::local(&field[0], &domain, &mut dec)
                .unwrap()
                .bernstein(&mut dec)
                .unwrap()
                .range(&mut dec)
                .unwrap();
            assert_eq!(range[0].sign(), gap.sign());
            assert_eq!(range[1].sign(), gap.sign());
        }
    }
}
