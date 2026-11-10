//! Common positive denominators and uniform natural interpolation.

use super::*;
use catallax_recipe::LoftShape;
use catallax_shapes::BSplineSurface;

#[cfg(test)]
mod tests;

pub(super) type H = [Real; 4];

pub(super) struct Interval {
    pub bounds: [Rational; 2],
    pub sections: Vec<Vec<H>>,
}

fn choose(n: usize) -> Vec<Rational> {
    let mut values = vec![Rational::one()];
    for k in 1..=n {
        values.push(
            &values[k - 1]
                * &Rational::from_ratio((n + 1 - k) as i64, k as i64).expect("positive index"),
        );
    }
    values
}

pub(super) fn product(a: &[Real], b: &[Real]) -> Vec<Real> {
    if b.len() == 1 {
        return a.iter().map(|x| x * &b[0]).collect();
    }
    if a.len() == 1 {
        return b.iter().map(|x| x * &a[0]).collect();
    }
    let n = a.len() - 1;
    let m = b.len() - 1;
    let ca = choose(n);
    let cb = choose(m);
    let cc = choose(n + m);
    let mut result = vec![Real::zero(); n + m + 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            let factor = Real::from_rational(&(&ca[i] * &cb[j]) / &cc[i + j]);
            result[i + j] = &result[i + j] + &(&(x * y) * &factor);
        }
    }
    result
}

pub(super) fn elevate<const N: usize>(mut net: Vec<[Real; N]>, degree: usize) -> Vec<[Real; N]> {
    while net.len() <= degree {
        let n = net.len();
        let mut next = Vec::with_capacity(n + 1);
        next.push(net[0].clone());
        for k in 1..n {
            let alpha = Real::from_rational(
                Rational::from_ratio(k as i64, n as i64).expect("positive degree"),
            );
            let beta = &Real::one() - &alpha;
            next.push(std::array::from_fn(|i| {
                &net[k - 1][i] * &alpha + &net[k][i] * &beta
            }));
        }
        next.push(net[n - 1].clone());
        net = next;
    }
    net
}

fn scalar_degree(values: Vec<Real>, degree: usize) -> Vec<Real> {
    elevate(values.into_iter().map(|x| [x]).collect(), degree)
        .into_iter()
        .map(|[x]| x)
        .collect()
}

fn equal(a: &[Real], b: &[Real], dec: &mut Decider) -> Result<bool, Refusal> {
    if a.len() != b.len() {
        return Ok(false);
    }
    for (a, b) in a.iter().zip(b) {
        if !a.same_definition(b) && sign(&(a - b), dec)? != Sign::Zero {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Replaces every section's quadratic chart on one common interval by the same
/// exact Mobius reparameterization when that makes all of them rational.
///
/// Scaling the Bernstein weights of a rational Bezier chart by `lambda^i` keeps
/// its Euclidean controls and traces the same curve: the new chart at local
/// parameter `t` is the old chart at `lambda t / ((1 - t) + lambda t)`. That
/// map fixes both interval ends and is strictly increasing for positive
/// `lambda`. Applying one map to every section keeps each correspondence pair
/// `C_i(u)` together, so the loft surface is the same point set with the same
/// station curves; only the carrier's u speed changes.
///
/// With weights normalized to `(1, a, b)`, `lambda = a / b` gives
/// `(1, k, k)`, where `k = a^2 / b` is the chart's reparameterization
/// invariant. An equal-angle circular span has `k = cos^2(angle / 2)`, which is
/// rational for full circles, so their algebraic charts become rational. Every
/// replaced coordinate is proved rational by exact field reduction and every
/// replaced weight positive, which makes `lambda` positive because the source
/// weights are; anything else leaves the interval unchanged.
pub(super) fn rational_interval(inputs: &mut [Vec<H>], dec: &mut Decider) -> Result<(), Refusal> {
    if inputs.iter().any(|input| input.len() != 3)
        || inputs
            .iter()
            .flatten()
            .flatten()
            .all(|x| x.rational_ref().is_some())
    {
        return Ok(());
    }
    let lambda = &inputs[0][1][3] / &inputs[0][2][3];
    let powers = [Real::one(), lambda.clone(), lambda.square()];
    let mut result = Vec::with_capacity(inputs.len());
    for input in inputs.iter() {
        let mut net = Vec::with_capacity(3);
        for (control, power) in input.iter().zip(&powers) {
            let factor = power / &input[0][3];
            let mut exact = Vec::with_capacity(4);
            for x in control {
                let Some(value) = dec.rational_value(&(x * &factor))? else {
                    return Ok(());
                };
                exact.push(Real::from_rational(value));
            }
            if exact[3]
                .rational_ref()
                .is_none_or(|weight| weight.sign() != Sign::Positive)
            {
                return Ok(());
            }
            net.push(<[Real; 4]>::try_from(exact).expect("homogeneous control"));
        }
        result.push(net);
    }
    for (input, net) in inputs.iter_mut().zip(result) {
        *input = net;
    }
    Ok(())
}

pub(super) fn common(
    bounds: [Rational; 2],
    inputs: Vec<Vec<H>>,
    dec: &mut Decider,
) -> Result<Interval, Refusal> {
    let mut numerators = Vec::new();
    let mut denominators: Vec<Vec<Real>> = Vec::new();
    let mut which = Vec::new();
    for input in inputs {
        let scale = &input[0][3];
        let normalized: Vec<H> = input
            .iter()
            .map(|p| p.each_ref().map(|x| x / scale))
            .collect();
        let weights: Vec<_> = normalized.iter().map(|p| p[3].clone()).collect();
        if weights.iter().all(|w| w.same_definition(&weights[0]))
            || equal(&weights, &vec![Real::one(); weights.len()], dec)?
        {
            which.push(None);
        } else {
            let mut found = None;
            for (index, old) in denominators.iter().enumerate() {
                if equal(old, &weights, dec)? {
                    found = Some(index);
                    break;
                }
            }
            let index = found.unwrap_or(denominators.len());
            if found.is_none() {
                denominators.push(weights);
            }
            which.push(Some(index));
        }
        numerators.push(normalized);
    }
    let all = denominators
        .iter()
        .fold(vec![Real::one()], |p, d| product(&p, d));
    let others: Vec<_> = (0..denominators.len())
        .map(|index| {
            denominators
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != index)
                .fold(vec![Real::one()], |p, (_, d)| product(&p, d))
        })
        .collect();
    let mut values: Vec<[Vec<Real>; 3]> = numerators
        .iter()
        .zip(which)
        .map(|(net, which)| {
            let multiplier = which.map_or(&all, |index| &others[index]);
            std::array::from_fn(|axis| {
                product(
                    &net.iter().map(|p| p[axis].clone()).collect::<Vec<_>>(),
                    multiplier,
                )
            })
        })
        .collect();
    let degree = values
        .iter()
        .flat_map(|v| v.iter().map(Vec::len))
        .max()
        .expect("loft sections")
        - 1;
    let denominator = scalar_degree(all, degree);
    let sections = values
        .iter_mut()
        .map(|xyz| {
            for coordinates in xyz.iter_mut() {
                *coordinates = scalar_degree(std::mem::take(coordinates), degree);
            }
            (0..=degree)
                .map(|i| {
                    [
                        xyz[0][i].clone(),
                        xyz[1][i].clone(),
                        xyz[2][i].clone(),
                        denominator[i].clone(),
                    ]
                })
                .collect()
        })
        .collect();
    Ok(Interval { bounds, sections })
}

/// Cubic controls at integer stations; the tridiagonal solve fixes both end
/// second derivatives to zero. No chord length enters this linear system.
pub(super) fn stations<const N: usize>(
    values: &[[Real; N]],
    shape: LoftShape,
) -> Vec<Vec<[Real; N]>> {
    if shape == LoftShape::Ruled || values.len() == 2 {
        return values
            .windows(2)
            .map(|p| vec![p[0].clone(), p[1].clone()])
            .collect();
    }
    let n = values.len();
    let mut second: Vec<[Real; N]> = vec![std::array::from_fn(|_| Real::zero()); n];
    let mut upper = vec![Rational::zero(); n];
    for i in 1..n - 1 {
        let denominator = &Rational::from_i64(4) - &upper[i - 1];
        upper[i] = &Rational::one() / &denominator;
        second[i] = std::array::from_fn(|axis| {
            let change = &(&values[i + 1][axis] - &(&values[i][axis] * &Real::from_i64(2)))
                + &values[i - 1][axis];
            (&(change * Real::from_i64(6)) - &second[i - 1][axis])
                / Real::from_rational(denominator.clone())
        });
    }
    for i in (1..n - 1).rev() {
        second[i] = std::array::from_fn(|axis| {
            &second[i][axis] - &(&second[i + 1][axis] * &Real::from_rational(upper[i].clone()))
        });
    }
    (0..n - 1)
        .map(|i| {
            vec![
                values[i].clone(),
                std::array::from_fn(|axis| {
                    &(&values[i][axis]
                        + &((&values[i + 1][axis] - &values[i][axis]) / Real::from_i64(3)))
                        - &((&(&second[i][axis] * &Real::from_i64(2)) + &second[i + 1][axis])
                            / Real::from_i64(18))
                }),
                std::array::from_fn(|axis| {
                    &(&values[i + 1][axis]
                        - &((&values[i + 1][axis] - &values[i][axis]) / Real::from_i64(3)))
                        - &((&second[i][axis] + &(&second[i + 1][axis] * &Real::from_i64(2)))
                            / Real::from_i64(18))
                }),
                values[i + 1].clone(),
            ]
        })
        .collect()
}

fn normalized(net: &[H], adjust: bool) -> std::borrow::Cow<'_, [H]> {
    let first = &net[0][3];
    if !adjust {
        if first
            .rational_ref()
            .is_some_and(|value| *value == Rational::one())
        {
            return std::borrow::Cow::Borrowed(net);
        }
        return std::borrow::Cow::Owned(
            net.iter()
                .map(|p| p.each_ref().map(|x| x / first))
                .collect(),
        );
    }
    let last = &net[net.len() - 1][3];
    let n = net.len();
    std::borrow::Cow::Owned(
        (0..=n)
            .map(|i| {
                let alpha = Real::from_rational(
                    Rational::from_ratio(i as i64, n as i64).expect("positive degree"),
                );
                let beta = &Real::one() - &alpha;
                std::array::from_fn(|axis| {
                    let low = if i < n {
                        &(&net[i][axis] / first) * &beta
                    } else {
                        Real::zero()
                    };
                    let high = if i > 0 {
                        &(&net[i - 1][axis] / last) * &alpha
                    } else {
                        Real::zero()
                    };
                    low + high
                })
            })
            .collect(),
    )
}

pub(super) fn surface(
    intervals: &[Interval],
    shape: LoftShape,
    dec: &mut Decider,
) -> Result<BSplineSurface, Refusal> {
    // A loop whose common charts are all rational carries a running rational
    // homogeneous scale: each interval starts at the previous interval's end
    // weight, so adjacent charts share their station row without raising the
    // u degree. Other loops normalize both endpoint weights to one with a
    // positive linear factor below.
    let rational = intervals
        .iter()
        .flat_map(|part| part.sections.iter().flatten().flatten())
        .all(|x| x.rational_ref().is_some());
    let mut factors = Vec::new();
    if rational {
        let mut scale = Rational::one();
        for part in intervals {
            let section = &part.sections[0];
            let weight = |i: usize| section[i][3].rational_ref().expect("rational loop");
            let factor = &scale / weight(0);
            scale = weight(section.len() - 1) * &factor;
            factors.push(factor);
        }
    }
    let adjustments: Vec<_> = if rational {
        vec![false; intervals.len()]
    } else {
        intervals
            .iter()
            .map(|part| {
                let section = &part.sections[0];
                sign(&(&section[0][3] - &section[section.len() - 1][3]), dec)
                    .map(|s| s != Sign::Zero)
            })
            .collect::<Result<_, _>>()?
    };
    let p = intervals
        .iter()
        .zip(&adjustments)
        .map(|(part, adjust)| part.sections[0].len() - 1 + usize::from(*adjust))
        .max()
        .expect("closed loop charts");
    let n = intervals[0].sections.len();
    let q = if shape == LoftShape::Ruled || n == 2 {
        1
    } else {
        3
    };
    let nv = (n - 1) * q + 1;
    let mut u_knots = vec![intervals[0].bounds[0].clone(); p + 1];
    let nu = intervals.len() * p + 1;
    let mut net = Vec::with_capacity(nu * nv);
    let mut previous: Vec<H> = Vec::new();
    for (index, (part, adjust)) in intervals.iter().zip(adjustments).enumerate() {
        let sections: Vec<_> = part
            .sections
            .iter()
            .map(|row| {
                let row = match factors.get(index) {
                    Some(factor) if *factor != Rational::one() => std::borrow::Cow::Owned(
                        row.iter()
                            .map(|point| {
                                point.each_ref().map(|x| {
                                    Real::from_rational(
                                        x.rational_ref().expect("rational loop") * factor,
                                    )
                                })
                            })
                            .collect(),
                    ),
                    Some(_) => std::borrow::Cow::Borrowed(row.as_slice()),
                    None => normalized(row, adjust),
                };
                if row.len() <= p {
                    std::borrow::Cow::Owned(elevate(row.into_owned(), p))
                } else {
                    row
                }
            })
            .collect();
        if !previous.is_empty() {
            for (section, old) in sections.iter().zip(&previous) {
                for (a, b) in section[0].iter().zip(old) {
                    if sign(&(a - b), dec)? != Sign::Zero {
                        return Err(unproven(
                            "adjacent exact loft charts have no proved common boundary",
                        ));
                    }
                }
            }
        }
        previous.clear();
        previous.extend(sections.iter().map(|section| section[p].clone()));
        for u in usize::from(index != 0)..=p {
            let values: Vec<[Real; 3]> = sections
                .iter()
                .map(|row| std::array::from_fn(|axis| row[u][axis].clone()))
                .collect();
            for (band, controls) in stations(&values, shape).into_iter().enumerate() {
                for xyz in controls.into_iter().skip(usize::from(band != 0)) {
                    net.push([
                        xyz[0].clone(),
                        xyz[1].clone(),
                        xyz[2].clone(),
                        sections[0][u][3].clone(),
                    ]);
                }
            }
        }
        // A running scale or a positive linear homogeneous factor makes
        // adjacent endpoint weights agree. After exact endpoint identity,
        // adjacent charts share their complete station row without a
        // discontinuous knot.
        u_knots.extend(std::iter::repeat_n(
            part.bounds[1].clone(),
            p + usize::from(index + 1 == intervals.len()),
        ));
    }
    let mut v_knots = vec![Rational::zero(); q + 1];
    for station in 1..n {
        v_knots.extend(std::iter::repeat_n(
            Rational::from_i64(station as i64),
            q + usize::from(station + 1 == n),
        ));
    }
    BSplineSurface::from_homogeneous([p, q], [u_knots, v_knots], [nu, nv], net, dec)
}
