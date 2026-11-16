//! Coordinate interval Newton before the coupled finite graph solve.

use super::*;

pub(super) fn unchanged(a: &[Ball; 4], b: &[Ball; 4]) -> bool {
    a.iter()
        .zip(b)
        .all(|(a, b)| a.lo() == b.lo() && a.hi() == b.hi())
}

pub(super) fn coordinate(
    equation: &RealPoly,
    derivative: &RealPoly,
    axis: usize,
    at: &mut [Ball; 4],
    parameter: &Ball,
    dec: &Decider,
) -> Result<(), Refusal> {
    dec.check_deadline()?;
    let mut values: [Ball; 5] = std::array::from_fn(|i| {
        if i < 4 {
            at[i].clone()
        } else {
            parameter.clone()
        }
    });
    let slope = derivative.eval_ball(&values, dec)?;
    if slope.contains_zero() {
        return Ok(());
    }
    let middle = Ball::from_rational(&at[axis].mid(), parameter.precision());
    values[axis] = middle.clone();
    let residual = equation.eval_ball(&values, dec)?;
    let Some(delta) = residual.checked_div(&slope) else {
        return Ok(());
    };
    // Holding the other coordinates in their complete intervals makes this
    // image contain every root in the current box, including its boundary.
    let image = &middle - &delta;
    let lo = at[axis].lo().max(image.lo());
    let hi = at[axis].hi().min(image.hi());
    if lo > hi {
        return Err(failure(
            "coordinate Newton contraction lost its certified graph",
        ));
    }
    at[axis] = Ball::from_rational(&lo, parameter.precision())
        .hull(&Ball::from_rational(&hi, parameter.precision()));
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_pinned_coordinate_contracts_before_the_coupled_jacobian_is_regular() {
        let dec = Decider::default();
        let x = RealPoly::var(5, 0);
        let y = RealPoly::var(5, 1);
        let z = RealPoly::var(5, 2);
        let w = RealPoly::var(5, 3);
        let t = RealPoly::var(5, 4);
        let equations = [
            &y - &x,
            &z - &y,
            &(&w * &w) + &(&x * &w)
                - RealPoly::constant(5, Real::from(Rational::from_ratio(3, 2).unwrap())),
            &x - &t,
        ];
        let derivative = equations
            .each_ref()
            .map(|p| std::array::from_fn::<_, 4, _>(|axis| p.partial(axis)));
        let mut at = std::array::from_fn(|_| {
            Ball::from_rational(&Rational::zero(), 128)
                .hull(&Ball::from_rational(&Rational::from_i64(2), 128))
        });
        let parameter = Ball::from_rational(&Rational::from_ratio(1, 2).unwrap(), 128);
        let initial: [Ball; 5] = std::array::from_fn(|i| {
            if i < 4 {
                at[i].clone()
            } else {
                parameter.clone()
            }
        });
        let jacobian = derivative
            .each_ref()
            .map(|row| row.each_ref().map(|p| p.eval_ball(&initial, &dec).unwrap()));
        assert!(Factor::new(jacobian, &dec).unwrap().is_none());
        for row in [3, 0, 1, 2] {
            for axis in 0..4 {
                coordinate(
                    &equations[row],
                    &derivative[row][axis],
                    axis,
                    &mut at,
                    &parameter,
                    &dec,
                )
                .unwrap();
            }
        }
        for (axis, expected) in [
            Rational::from_ratio(1, 2).unwrap(),
            Rational::from_ratio(1, 2).unwrap(),
            Rational::from_ratio(1, 2).unwrap(),
            Rational::one(),
        ]
        .iter()
        .enumerate()
        {
            assert!(at[axis].lo() <= *expected && at[axis].hi() >= *expected);
        }
        assert!(at[3].hi() - at[3].lo() < Rational::from_i64(2));
        let narrowed: [Ball; 5] = std::array::from_fn(|i| {
            if i < 4 {
                at[i].clone()
            } else {
                parameter.clone()
            }
        });
        let jacobian = derivative.each_ref().map(|row| {
            row.each_ref()
                .map(|p| p.eval_ball(&narrowed, &dec).unwrap())
        });
        assert!(Factor::new(jacobian, &dec).unwrap().is_some());
    }
}
