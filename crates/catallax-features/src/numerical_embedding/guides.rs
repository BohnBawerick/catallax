//! A regular clock and shared rational jets across retained trace charts.

use super::{exact, fit_failure, nurbs, rational_curve, rational_curve_definition};
use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::{BSplineCurve, Vector3};

type Fit = (Vec<f64>, Vec<[f64; 3]>);

fn endpoint_jets((knots, controls): &Fit) -> [[f64; 3]; 2] {
    let degree = knots.len() - controls.len() - 1;
    let last = controls.len() - 1;
    let width = knots[controls.len()] - knots[degree];
    let first_scale = degree as f64 * width / (knots[degree + 1] - knots[degree]);
    let last_scale = degree as f64 * width / (knots[last + 1] - knots[last]);
    [
        std::array::from_fn(|i| (controls[1][i] - controls[0][i]) * first_scale),
        std::array::from_fn(|i| (controls[last][i] - controls[last - 1][i]) * last_scale),
    ]
}

fn speed(value: &[f64; 3]) -> Result<f64, Refusal> {
    let length = value
        .iter()
        .fold(0.0_f64, |sum, coordinate| libm::hypot(sum, *coordinate));
    if !length.is_finite() || length <= 0.0 {
        return Err(fit_failure(
            "a retained guide chart has no finite nonzero endpoint speed",
        ));
    }
    Ok(length)
}

pub(super) fn join(
    mut fits: Vec<Fit>,
    ends: [[Rational; 3]; 2],
    dec: &mut Decider,
) -> Result<(BSplineCurve, Vec<Rational>), Refusal> {
    if fits.len() == 1 {
        let (knots, controls) = fits.pop().expect("one guide chart");
        let curve = rational_curve(
            &knots,
            &controls,
            [&Rational::zero(), &Rational::one()],
            ends,
            None,
            dec,
        )?;
        return Ok((curve, vec![Rational::zero(), Rational::one()]));
    }
    let jets: Vec<_> = fits.iter().map(endpoint_jets).collect();
    let mut widths = vec![1.0; fits.len()];
    for i in 1..widths.len() {
        dec.charge_work(1)?;
        widths[i] = widths[i - 1] * (speed(&jets[i][0])? / speed(&jets[i - 1][1])?);
    }
    let total: f64 = widths.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(fit_failure(
            "the retained guide clock is outside binary64 range",
        ));
    }
    let mut clock = Vec::with_capacity(fits.len() + 1);
    clock.push(0.0);
    for width in widths {
        clock.push(clock.last().unwrap() + width / total);
    }
    *clock.last_mut().unwrap() = 1.0;
    if clock
        .windows(2)
        .any(|pair| !pair[1].is_finite() || pair[0] >= pair[1])
    {
        return Err(fit_failure(
            "the retained guide clock has a collapsed chart",
        ));
    }
    let mut joins = Vec::with_capacity(fits.len() - 1);
    let mut index = fits[0].1.len() - 1;
    for i in 1..fits.len() {
        let left_width = clock[i] - clock[i - 1];
        let right_width = clock[i + 1] - clock[i];
        let direction: [f64; 3] = std::array::from_fn(|axis| {
            (jets[i - 1][1][axis] / left_width + jets[i][0][axis] / right_width) * 0.5
        });
        let [x, y, z] = direction;
        joins.push((
            index,
            Vector3::new(exact(x)?.into(), exact(y)?.into(), exact(z)?.into()),
        ));
        index += fits[i].1.len() - 1;
    }
    for (i, (knots, controls)) in fits.iter_mut().enumerate() {
        let degree = knots.len() - controls.len() - 1;
        let lo = knots[degree];
        let hi = knots[controls.len()];
        let width = clock[i + 1] - clock[i];
        let last = knots.len() - degree - 1;
        for (k, knot) in knots.iter_mut().enumerate() {
            *knot = if k <= degree {
                clock[i]
            } else if k >= last {
                clock[i + 1]
            } else {
                clock[i] + ((*knot - lo) / (hi - lo)) * width
            };
        }
    }
    let (knots, controls) = nurbs::concatenate(fits);
    let (degree, knots, mut points) = rational_curve_definition(
        &knots,
        &controls,
        [&Rational::zero(), &Rational::one()],
        ends,
        None,
        dec,
    )?;
    let degree_scalar = Rational::from_i64(degree as i64);
    for (index, direction) in joins {
        dec.charge_work(1)?;
        let at = &knots[index + 1];
        let before = Real::from_rational(-(&(at - &knots[index]) / &degree_scalar));
        let after = Real::from_rational(&(&knots[index + degree + 1] - at) / &degree_scalar);
        // One exact jet, not independently rounded Cartesian control handles.
        let point = &points[index];
        let left = point.plus(&direction.scaled(&before));
        let right = point.plus(&direction.scaled(&after));
        points[index - 1] = left;
        points[index + 1] = right;
    }
    let weights = vec![Real::one(); points.len()];
    let curve = BSplineCurve::new(degree, knots, points, weights, dec)?;
    Ok((
        curve,
        clock.into_iter().map(exact).collect::<Result<_, _>>()?,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unequal_chart_speeds_share_one_exact_cartesian_join_jet() {
        for reversed in [false, true] {
            let mut xs = [[0, 1, 2, 3], [3, 5, 7, 9]];
            if reversed {
                xs.reverse();
                for piece in &mut xs {
                    piece.reverse();
                }
            }
            let fits = xs
                .into_iter()
                .enumerate()
                .map(|(index, xs)| {
                    let bounds = [index as f64 * 0.5, (index + 1) as f64 * 0.5];
                    let knots = [vec![bounds[0]; 4], vec![bounds[1]; 4]].concat();
                    let controls = xs
                        .map(|x| [f64::from(x), 1_000_000.0, -1_000_000.0])
                        .to_vec();
                    (knots, controls)
                })
                .collect();
            let ends = if reversed { [9, 0] } else { [0, 9] }.map(|x| {
                [
                    Rational::from_i64(x),
                    Rational::from_i64(1_000_000),
                    Rational::from_i64(-1_000_000),
                ]
            });
            let mut dec = Decider::default();
            let (curve, clock) = join(fits, ends, &mut dec).unwrap();
            let at = &clock[1];
            let jet = curve.jet_box([at, at], 128, &dec).unwrap();
            for (value, expected) in jet.value.iter().zip([
                Rational::from_i64(3),
                Rational::from_i64(1_000_000),
                Rational::from_i64(-1_000_000),
            ]) {
                assert!(value.contains(&expected));
            }
            for (value, expected) in jet.first.iter().zip([
                Rational::from_i64(if reversed { -9 } else { 9 }),
                Rational::zero(),
                Rational::zero(),
            ]) {
                assert!(value.contains(&expected));
                assert!(value.hi() - value.lo() < Rational::one().mul_pow2(-64));
            }
        }
    }
}
