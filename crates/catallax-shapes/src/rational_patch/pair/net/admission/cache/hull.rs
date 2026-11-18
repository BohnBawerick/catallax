use super::super::{Decider, Grid, Rational, Refusal, Value, less};

fn ordered(grid: &Grid, axis: usize, dec: &mut Decider) -> Result<bool, Refusal> {
    let degree = grid.deg[axis];
    if degree == 0 {
        return Ok(true);
    }
    let lines = grid.deg[1 - axis] + 1;
    let mut ascending = None;
    for line in 0..lines {
        let index = |k| {
            if axis == 0 {
                k * (grid.deg[1] + 1) + line
            } else {
                line * (grid.deg[1] + 1) + k
            }
        };
        for k in 0..degree {
            let a = &grid.v[index(k)];
            let b = &grid.v[index(k + 1)];
            match ascending {
                None => {
                    if !<Rational as Value>::same(a, b, dec)? {
                        ascending = Some(less(a, b, dec)?);
                    }
                }
                Some(true) if less(b, a, dec)? => return Ok(false),
                Some(false) if less(a, b, dec)? => return Ok(false),
                _ => (),
            }
        }
    }
    Ok(true)
}

pub(super) fn monotone(grid: &Grid, dec: &mut Decider) -> Result<[bool; 2], Refusal> {
    Ok([ordered(grid, 0, dec)?, ordered(grid, 1, dec)?])
}

// Positive tensor restrictions preserve a common control-row order. Its
// extrema therefore stay on the boundary rows, with the original hull values.
pub(super) fn range(
    grid: &Grid,
    monotone: [bool; 2],
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    let [m, n] = grid.deg;
    let rows = (0..=m).filter(|&i| !monotone[0] || i == 0 || i == m);
    let mut values = rows.flat_map(|i| {
        (0..=n)
            .filter(move |&j| !monotone[1] || j == 0 || j == n)
            .map(move |j| &grid.v[i * (n + 1) + j])
    });
    let first = values.next().expect("nonempty Bernstein net");
    let Some(second) = values.next() else {
        return Ok([first.clone(), first.clone()]);
    };
    let (mut low, mut high) = if less(first, second, dec)? {
        (first, second)
    } else {
        (second, first)
    };
    while let Some(a) = values.next() {
        if let Some(b) = values.next() {
            let (next_low, next_high) = if less(a, b, dec)? { (a, b) } else { (b, a) };
            if less(next_low, low, dec)? {
                low = next_low;
            }
            if less(high, next_high, dec)? {
                high = next_high;
            }
        } else {
            if less(a, low, dec)? {
                low = a;
            }
            if less(high, a, dec)? {
                high = a;
            }
        }
    }
    Ok([low.clone(), high.clone()])
}

#[cfg(test)]
mod tests {
    use super::super::super::restricted_net;
    use super::*;
    use catallax_numbers::with_operation_budget;

    #[test]
    fn opposite_row_orders_keep_the_interior_value_after_restriction() {
        let ((), work) = with_operation_budget(1_000_000, || {
            let mut dec = Decider::default();
            let grid = Grid {
                deg: [2, 1],
                v: [0, 11, 10, 10, 11, 0].map(Rational::from_i64).to_vec(),
            };
            let order = monotone(&grid, &mut dec).unwrap();
            let cuts = [
                [Rational::zero(), Rational::one()],
                [
                    Rational::from_ratio(49, 51).unwrap(),
                    Rational::from_ratio(51, 100).unwrap(),
                ],
            ];
            let clipped = restricted_net(&grid, &cuts, &mut dec).unwrap();
            let actual = range(clipped.as_ref(), order, &mut dec).unwrap();
            assert_eq!(actual, clipped.range(&mut dec).unwrap());
            let interior = Rational::from_ratio(31, 4).unwrap();
            assert!(actual[0] <= interior && interior <= actual[1]);
        });
        assert!(!work.exhausted);
    }
}
