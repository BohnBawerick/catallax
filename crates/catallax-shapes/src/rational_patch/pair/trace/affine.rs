//! Eliminate affine source constraints before enclosing the same finite roots.

use super::*;

fn substitute(
    poly: &MPoly,
    axis: usize,
    value: &MPoly,
    dec: &mut Decider,
) -> Result<MPoly, Refusal> {
    dec.check_deadline()?;
    let terms = poly.terms();
    dec.charge_work(terms.len() as u64)?;
    let degree = terms
        .iter()
        .map(|(m, _)| m.exponents()[axis])
        .max()
        .unwrap_or(0);
    let mut powers = vec![MPoly::constant(poly.nvars(), Rational::one())];
    for _ in 0..degree {
        dec.charge_work(1)?;
        powers.push(product(powers.last().unwrap(), value, dec)?);
    }
    let mut result = MPoly::zero(poly.nvars());
    for (monomial, coefficient) in terms {
        dec.charge_work(1)?;
        let mut exponents = monomial.exponents().to_vec();
        let power = exponents[axis] as usize;
        exponents[axis] = 0;
        let negative = MPoly::from_terms(poly.nvars(), vec![(exponents, -coefficient)]);
        result = difference(&result, &product(&negative, &powers[power], dec)?, dec)?;
    }
    Ok(result)
}

pub(super) fn reduce(mut equations: [MPoly; 3], dec: &mut Decider) -> Result<[MPoly; 3], Refusal> {
    let mut used_rows = [false; 3];
    let mut used_axes = [false; 4];
    for _ in 0..3 {
        let mut selected = None;
        for (row, equation) in equations.iter().enumerate() {
            dec.charge_work(1)?;
            if used_rows[row] || equation.total_degree() != Some(1) {
                continue;
            }
            let terms = equation.terms();
            dec.charge_work(terms.len() as u64)?;
            for axis in (0..4).rev() {
                dec.charge_work(terms.len() as u64 + 1)?;
                if !used_axes[axis]
                    && let Some((_, coefficient)) =
                        terms.iter().find(|(m, _)| m.exponents()[axis] == 1)
                {
                    selected = Some((row, axis, coefficient.clone(), terms));
                    break;
                }
            }
            if selected.is_some() {
                break;
            }
        }
        let Some((row, axis, coefficient, terms)) = selected else {
            break;
        };
        let mut replacement = Vec::new();
        for (monomial, value) in terms {
            dec.charge_work(1)?;
            if monomial.exponents()[axis] == 0 {
                dec.reserve_exact_arithmetic(&value, &coefficient, true)?;
                replacement.push((monomial.exponents().to_vec(), -value / coefficient.clone()));
            }
        }
        let replacement = MPoly::from_terms(equations[row].nvars(), replacement);
        used_rows[row] = true;
        used_axes[axis] = true;
        for (other, equation) in equations.iter_mut().enumerate() {
            dec.charge_work(1)?;
            if other != row {
                *equation = substitute(equation, axis, &replacement, dec)?;
            }
        }
    }
    Ok(equations)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn affine_source_constraints_remove_a_false_interval_singularity() {
        let mut dec = Decider::default();
        let u = MPoly::var(5, 0);
        let v = MPoly::var(5, 2);
        let z = MPoly::var(5, 1);
        let w = MPoly::var(5, 3);
        let equations = [
            &(&u - &v) - &MPoly::constant(5, Rational::from_ratio(1, 4).unwrap()),
            &z - &w,
            &(&u - &(&u * &u)) - &(&v - &(&v * &v)),
        ];
        let parameter = RealPoly::var(5, 1) - RealPoly::var(5, 4);
        let at: [Ball; 5] = std::array::from_fn(|_| {
            Ball::from_rational(&Rational::zero(), 128)
                .hull(&Ball::from_rational(&Rational::one(), 128))
        });
        let matrix = |polys: &[MPoly; 3], dec: &Decider| {
            std::array::from_fn(|row| {
                let p = if row < 3 {
                    RealPoly::from(polys[row].clone())
                } else {
                    parameter.clone()
                };
                std::array::from_fn(|axis| p.partial(axis).eval_ball(&at, dec).unwrap())
            })
        };
        assert!(
            Factor::new(matrix(&equations, &dec), &dec)
                .unwrap()
                .is_none()
        );
        let reduced = reduce(equations.clone(), &mut dec).unwrap();
        assert!(Factor::new(matrix(&reduced, &dec), &dec).unwrap().is_some());
        for u in [
            Rational::from_ratio(5, 8).unwrap(),
            Rational::from_ratio(1, 2).unwrap(),
        ] {
            let point = [
                u.clone(),
                Rational::from_ratio(1, 2).unwrap(),
                &u - &Rational::from_ratio(1, 4).unwrap(),
                Rational::from_ratio(1, 2).unwrap(),
                Rational::from_ratio(1, 2).unwrap(),
            ];
            let before = equations.iter().all(|p| p.eval(&point).is_zero());
            let after = reduced.iter().all(|p| p.eval(&point).is_zero());
            assert_eq!(before, u == Rational::from_ratio(5, 8).unwrap());
            assert_eq!(before, after);
        }
    }
}
