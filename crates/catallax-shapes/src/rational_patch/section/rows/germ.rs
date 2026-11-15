//! Inward signs require a nonzero leading coefficient along the complete row.

use super::*;

pub(in crate::rational_patch::section) fn inward_sign(
    patch: &RationalPatch,
    surface: &SectionSurface,
    ends: &[[Real; 2]; 2],
    dec: &mut Decider,
) -> Result<Option<Sign>, Refusal> {
    for axis in 0..2 {
        if exact_sign(&(&ends[0][axis] - &ends[1][axis]), dec)? != Sign::Zero {
            continue;
        }
        for upper in [false, true] {
            let bound = &patch.domain()[axis][usize::from(upper)];
            let value = Real::from_rational(bound.clone());
            if exact_sign(&(&ends[0][axis] - &value), dec)? != Sign::Zero {
                continue;
            }
            let other = 1 - axis;
            let [low, high] =
                if exact_sign(&(&ends[1][other] - &ends[0][other]), dec)? == Sign::Positive {
                    [&ends[0][other], &ends[1][other]]
                } else {
                    [&ends[1][other], &ends[0][other]]
                };
            let mut found = None;
            for span in &patch.data.spans {
                if span.domain[axis][usize::from(upper)] != *bound
                    || exact_sign(&(&Real::from(span.domain[other][1].clone()) - low), dec)?
                        != Sign::Positive
                    || exact_sign(&(high - &Real::from(span.domain[other][0].clone())), dec)?
                        != Sign::Positive
                {
                    continue;
                }
                let mut residual = surface.equation(span);
                if residual.is_zero(dec)? {
                    return Ok(None);
                }
                let mut order = 0;
                while residual
                    .substitute(axis, &RealPoly::constant(2, value.clone()))
                    .is_zero(dec)?
                {
                    let quotient = divide(&residual, axis, bound);
                    let factor = RealPoly::var(2, axis) - RealPoly::constant(2, value.clone());
                    if !(&residual - &(&factor * &quotient)).is_zero(dec)? {
                        return Ok(None);
                    }
                    residual = quotient;
                    order += 1;
                }
                if order == 0 {
                    return Ok(None);
                }
                let leading = residual.substitute(axis, &RealPoly::constant(2, value.clone()));
                if !nonzero(&leading, &span.domain, dec)? {
                    return Ok(None);
                }
                let middle = span
                    .domain
                    .each_ref()
                    .map(|[lo, hi]| Real::from((lo + hi).mul_pow2(-1)));
                let mut sign = exact_sign(&leading.eval_real(&middle), dec)?;
                if upper && order % 2 == 1 {
                    sign = match sign {
                        Sign::Negative => Sign::Positive,
                        Sign::Positive => Sign::Negative,
                        Sign::Zero => return Ok(None),
                    };
                }
                if found.is_some_and(|previous| previous != sign) {
                    return Ok(None);
                }
                found = Some(sign);
            }
            return Ok(found);
        }
    }
    Ok(None)
}
