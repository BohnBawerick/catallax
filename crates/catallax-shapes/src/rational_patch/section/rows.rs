//! Complete boundary components when the remaining cut equation stays nonzero.

use super::*;

mod germ;
pub(super) use germ::inward_sign;
pub(crate) struct BoundaryRows(Vec<(usize, Rational, bool)>);

impl BoundaryRows {
    pub(crate) fn rows(&self) -> impl Iterator<Item = (usize, &Rational, bool)> {
        self.0
            .iter()
            .map(|(axis, bound, crossing)| (*axis, bound, *crossing))
    }
}

fn ratio(index: usize, power: usize, degree: usize) -> Rational {
    (0..power).fold(Rational::one(), |value, step| {
        &value
            * &Rational::from_ratio((index - step) as i64, (degree - step) as i64)
                .expect("power does not exceed degree")
    })
}

fn bernstein_nonzero(poly: &RealPoly, domain: &Domain, dec: &mut Decider) -> Result<bool, Refusal> {
    let mut local = poly.clone();
    for (axis, [lo, hi]) in domain.iter().enumerate() {
        local = local.substitute(
            axis,
            &(RealPoly::var(2, axis).scale(&Real::from_rational(hi - lo))
                + RealPoly::constant(2, Real::from_rational(lo.clone()))),
        );
    }
    let terms = local.terms();
    let degree: [usize; 2] = std::array::from_fn(|axis| {
        terms
            .iter()
            .map(|(m, _)| m.exponents()[axis] as usize)
            .max()
            .unwrap_or(0)
    });
    let mut wanted = None;
    for u in 0..=degree[0] {
        for v in 0..=degree[1] {
            let mut coefficient = Real::zero();
            for (monomial, value) in &terms {
                let [a, b] = [
                    monomial.exponents()[0] as usize,
                    monomial.exponents()[1] as usize,
                ];
                if a <= u && b <= v {
                    let factor = ratio(u, a, degree[0]) * ratio(v, b, degree[1]);
                    coefficient = &coefficient + &(value * &Real::from_rational(factor));
                }
            }
            let sign = exact_sign(&coefficient, dec)?;
            if sign == Sign::Zero || wanted.is_some_and(|old| old != sign) {
                return Ok(false);
            }
            wanted = Some(sign);
        }
    }
    Ok(true)
}

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
        let mut local = poly.clone();
        let mut centered = cell.clone();
        for axis in 0..2 {
            let middle = (&cell[axis][0] + &cell[axis][1]).mul_pow2(-1);
            local = local.substitute(
                axis,
                &(RealPoly::var(2, axis)
                    + RealPoly::constant(2, Real::from_rational(middle.clone()))),
            );
            centered[axis] = cell[axis].each_ref().map(|bound| bound - &middle);
        }
        if !local
            .eval_ball(&box_balls(&centered, bits), dec)?
            .contains_zero()
        {
            continue;
        }
        if bernstein_nonzero(poly, &cell, dec)? {
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

/// Complete zero localization to rectangle sides; an empty list excludes all zeros.
pub(crate) fn boundary_rows(
    equation: &RealPoly,
    domain: &Domain,
    dec: &mut Decider,
) -> Result<Option<BoundaryRows>, Refusal> {
    let mut residual = equation.clone();
    if residual.is_zero(dec)? {
        return Ok(None);
    }
    let mut rows = Vec::new();
    for (axis, bounds) in domain.iter().enumerate() {
        for bound in bounds {
            let value = Real::from_rational(bound.clone());
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
            if order > 0 {
                rows.push((axis, bound.clone(), order % 2 == 1));
            }
        }
    }
    if !nonzero(&residual, domain, dec)? {
        return Ok(None);
    }
    Ok(Some(BoundaryRows(rows)))
}

pub(super) fn section(
    patch: &RationalPatch,
    surface: &SectionSurface,
    dec: &mut Decider,
) -> Result<Option<RationalPatchSection>, Refusal> {
    if !matches!(surface, SectionSurface::Plane(_)) {
        return Ok(None);
    }
    surface.validate(dec)?;
    let mut curves = Vec::new();
    for span in &patch.data.spans {
        let Some(rows) = boundary_rows(&surface.equation(span), &span.domain, dec)? else {
            return Ok(None);
        };
        let mut row_axis = None;
        for (axis, bound, crossing) in rows.rows() {
            if row_axis.is_some_and(|previous| previous != axis) {
                return Ok(None);
            }
            row_axis = Some(axis);
            curves.push((
                iso(span, axis, &Real::from_rational(bound.clone()), dec)?,
                crossing,
            ));
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
