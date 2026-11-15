//! Complete boundary components when the remaining cut equation stays nonzero.

use super::*;

fn divide(poly: &RealPoly, axis: usize, bound: &Rational) -> RealPoly {
    let mut terms = Vec::new();
    for (monomial, coefficient) in poly.terms() {
        let degree = monomial.exponents()[axis];
        let mut power = Real::one();
        for index in 0..degree {
            let mut exponents = monomial.exponents().to_vec();
            exponents[axis] = degree - 1 - index;
            terms.push((exponents, &coefficient * &power));
            power = &power * &Real::from_rational(bound.clone());
        }
    }
    RealPoly::from_terms(2, terms)
}

fn nonzero(poly: &RealPoly, domain: &Domain, dec: &mut Decider) -> Result<bool, Refusal> {
    let bits = dec.policy().max_precision_bits.min(128);
    if bits == 0 {
        return Ok(false);
    }
    let mut pending = vec![(domain.clone(), 0)];
    while let Some((cell, depth)) = pending.pop() {
        if !poly
            .eval_ball(&box_balls(&cell, bits), dec)?
            .contains_zero()
        {
            continue;
        }
        if depth == 8 {
            return Ok(false);
        }
        let axis = depth % 2;
        let middle = (&cell[axis][0] + &cell[axis][1]).mul_pow2(-1);
        let mut low = cell.clone();
        let mut high = cell;
        low[axis][1] = middle.clone();
        high[axis][0] = middle;
        pending.push((high, depth + 1));
        pending.push((low, depth + 1));
    }
    Ok(true)
}

pub(super) fn section(
    patch: &RationalPatch,
    surface: &SectionSurface,
    dec: &mut Decider,
) -> Result<Option<RationalPatchSection>, Refusal> {
    surface.validate(dec)?;
    let mut curves = Vec::new();
    for span in &patch.data.spans {
        let mut residual = surface.equation(span);
        if residual.is_zero(dec)? {
            return Ok(None);
        }
        let mut row_axis = None;
        for axis in 0..2 {
            for bound in &span.domain[axis] {
                let value = Real::from_rational(bound.clone());
                let mut order = 0;
                while residual
                    .substitute(axis, &RealPoly::constant(2, value.clone()))
                    .is_zero(dec)?
                {
                    if row_axis.is_some_and(|previous| previous != axis) {
                        return Ok(None);
                    }
                    row_axis = Some(axis);
                    let quotient = divide(&residual, axis, bound);
                    let factor = RealPoly::var(2, axis) - RealPoly::constant(2, value.clone());
                    if !(&residual - &(&factor * &quotient)).is_zero(dec)? {
                        return Ok(None);
                    }
                    residual = quotient;
                    order += 1;
                }
                if order > 0 {
                    curves.push((iso(span, axis, &value, dec)?, order % 2 == 1));
                }
            }
        }
        // A covered nonzero remainder proves there are no missing interior branches.
        if !nonzero(&residual, &span.domain, dec)? {
            return Ok(None);
        }
    }
    Ok(Some(RationalPatchSection {
        patch: patch.clone(),
        surface: surface.clone(),
        curves,
        points: Vec::new(),
        branches: Vec::new(),
        search: None,
    }))
}
