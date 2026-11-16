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

pub(super) fn reduce<const N: usize>(
    mut equations: [MPoly; N],
    dec: &mut Decider,
) -> Result<[MPoly; N], Refusal> {
    let mut used_rows = [false; N];
    let mut used_axes = [false; 4];
    for _ in 0..N {
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

pub(super) fn rational_point(
    equations: &[RealPoly],
    bounds: &[[Rational; 2]; 4],
    dec: &mut Decider,
) -> Result<Option<[Real; 4]>, Refusal> {
    let mut rational = Vec::new();
    for equation in equations {
        dec.charge_work(1)?;
        let Some(poly) = equation.rational() else {
            return Ok(None);
        };
        rational.push(poly.clone());
    }
    let original: [MPoly; 4] = rational.try_into().expect("four pinned graph equations");
    let reduced = reduce(original.clone(), dec)?;
    let mut point: [Option<Rational>; 4] = [None, None, None, None];
    for equation in &reduced {
        dec.charge_work(1)?;
        if equation.total_degree() != Some(1) {
            return Ok(None);
        }
        let terms = equation.terms();
        dec.charge_work(terms.len() as u64)?;
        let mut variable = None;
        let mut constant = Rational::zero();
        for (monomial, coefficient) in terms {
            if monomial.degree() == 0 {
                constant = coefficient;
                continue;
            }
            let axis = monomial.exponents().iter().position(|e| *e == 1).unwrap();
            if variable.is_some() || axis >= 4 {
                return Ok(None);
            }
            variable = Some((axis, coefficient));
        }
        let Some((axis, coefficient)) = variable else {
            return Ok(None);
        };
        if point[axis].is_some() {
            return Ok(None);
        }
        point[axis] = Some(-constant / coefficient);
    }
    let [Some(a), Some(b), Some(c), Some(d)] = point else {
        return Ok(None);
    };
    let point = [a, b, c, d];
    for (coordinate, [lo, hi]) in point.iter().zip(bounds) {
        dec.charge_work(1)?;
        if coordinate < lo || coordinate > hi {
            return Ok(None);
        }
    }
    for equation in &original {
        dec.check_deadline()?;
        let work = equation
            .terms()
            .iter()
            .try_fold(0_u64, |total, (m, _)| {
                total.checked_add(2 * u64::from(m.degree()) + 5)
            })
            .ok_or_else(|| failure("rational graph evaluation work overflowed"))?;
        dec.charge_work(work)?;
        if !equation.eval(&point).is_zero() {
            return Ok(None);
        }
    }
    Ok(Some(point.map(Real::from)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rational_graph_points_keep_the_pin_and_finite_native_bounds() {
        let mut dec = Decider::default();
        let u = MPoly::var(4, 0);
        let v = MPoly::var(4, 2);
        let z = MPoly::var(4, 1);
        let w = MPoly::var(4, 3);
        let height = MPoly::from_terms(
            4,
            vec![
                (vec![0, 0, 0, 0], Rational::from_i64(2)),
                (vec![0, 1, 0, 0], Rational::one()),
                (vec![0, 2, 0, 0], Rational::one()),
            ],
        );
        let curved = &(&u - &(&u * &u)) - &(&v - &(&v * &v));
        let equations = [
            &(&u - &v) - &MPoly::constant(4, Rational::from_ratio(1, 4).unwrap()),
            &z - &w,
            &curved * &height,
            &z - &MPoly::constant(4, Rational::from_ratio(1, 2).unwrap()),
        ]
        .map(RealPoly::from);
        let mut bounds = std::array::from_fn(|_| [Rational::zero(), Rational::one()]);
        let point = rational_point(&equations, &bounds, &mut dec)
            .unwrap()
            .unwrap();
        for (value, expected) in point.iter().zip([(5, 8), (1, 2), (3, 8), (1, 2)]) {
            assert_eq!(
                value.rational_ref(),
                Some(&Rational::from_ratio(expected.0, expected.1).unwrap())
            );
        }
        bounds[0][1] = Rational::from_ratio(1, 2).unwrap();
        assert!(
            rational_point(&equations, &bounds, &mut dec)
                .unwrap()
                .is_none()
        );
        let (answer, work) = catallax_numbers::with_operation_budget(0, || {
            rational_point(&equations, &bounds, &mut Decider::default())
        });
        assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
        assert!(work.exhausted);
    }

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
