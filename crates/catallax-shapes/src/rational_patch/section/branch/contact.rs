//! Exact corner contacts in a one-sided Bernstein cell.

use super::*;

pub(super) fn corners(net: &Net, dec: &mut Decider) -> Result<Option<Vec<[Real; 2]>>, Refusal> {
    let mut sign = None;
    for value in &net.values {
        let next = exact_sign(value, dec)?;
        if next == Sign::Zero {
            continue;
        }
        if sign.is_some_and(|previous| previous != next) {
            return Ok(None);
        }
        sign = Some(next);
    }
    if sign.is_none() {
        return Ok(None);
    }
    // A nonzero coefficient makes the cell interior strict. Each edge also
    // needs a nonzero coefficient, so an entire contact interval cannot pass.
    for edge in net.edges() {
        let mut strict = false;
        for value in edge {
            strict |= exact_sign(&value, dec)? != Sign::Zero;
        }
        if !strict {
            return Ok(None);
        }
    }
    let [nu, nv] = net.count;
    let mut points = Vec::new();
    for (i, u) in [(0, &net.domain[0][0]), (nu - 1, &net.domain[0][1])] {
        for (j, v) in [(0, &net.domain[1][0]), (nv - 1, &net.domain[1][1])] {
            if exact_sign(&net.values[i * nv + j], dec)? == Sign::Zero {
                points.push([
                    Real::from_rational(u.clone()),
                    Real::from_rational(v.clone()),
                ]);
            }
        }
    }
    Ok(Some(points))
}

// This local classifier only permits subdivision. The complete cell checks
// must still account for the contact and every other zero of the equation.
pub(super) fn definite_hessian(
    equation: &RealPoly,
    at: &[Real; 2],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let uu = equation.partial(0).partial(0).eval_real(at);
    let uv = equation.partial(0).partial(1).eval_real(at);
    let vv = equation.partial(1).partial(1).eval_real(at);
    reserve_values(&uu, &vv, true, dec)?;
    reserve_values(&uv, &uv, true, dec)?;
    let a = &uu * &vv;
    let b = &uv * &uv;
    reserve_values(&a, &b, false, dec)?;
    Ok(exact_sign(&(&a - &b), dec)? == Sign::Positive)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn only_one_sided_cells_with_nonvanishing_edges_admit_corner_contacts() {
        for (values, expected) in [
            ([0, 1, 1, 1, 1, 1, 1, 1, 1], Some(1)),
            ([0, -1, -1, -1, -1, -1, -1, -1, -1], Some(1)),
            ([1, 1, 1, 1, 0, 1, 1, 1, 1], Some(0)),
            ([0, 0, 0, 1, 1, 1, 1, 1, 1], None),
            ([0, 1, 1, 1, -1, 1, 1, 1, 1], None),
            ([0; 9], None),
        ] {
            let net = Net {
                count: [3, 3],
                values: values.map(Real::from_i64).to_vec(),
                domain: [
                    [Rational::zero(), Rational::one()],
                    [Rational::zero(), Rational::one()],
                ],
                depth: 0,
            };
            assert_eq!(
                corners(&net, &mut Decider::default())
                    .unwrap()
                    .map(|points| points.len()),
                expected,
                "{values:?}"
            );
        }
    }
}
