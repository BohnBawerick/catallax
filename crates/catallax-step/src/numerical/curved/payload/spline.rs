//! Native spline definitions with the actual written coordinates and parameters.

use super::*;

fn knots(
    source: &[Rational],
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<Vec<Rational>, Refusal> {
    let values = source
        .iter()
        .map(|q| parameter(q, errors, options, dec))
        .collect::<Result<Vec<_>, _>>()?;
    for (old, new) in source.windows(2).zip(values.windows(2)) {
        if old[0] < old[1] && new[0] >= new[1] {
            return Err(bad("rounding-topology", "distinct spline knots collapsed"));
        }
    }
    Ok(values)
}

fn weight(
    value: &Real,
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    let result = rounded(
        value,
        &mut errors.spline.max_weight_error,
        options,
        false,
        dec,
    )?;
    if result.sign() != Sign::Positive {
        return Err(bad("rounding-topology", "positive spline weight vanished"));
    }
    Ok(result.into())
}

pub(super) fn curve(
    source: &BSplineCurve,
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let knots = knots(source.knots(), errors, options, dec)?;
    let controls = source
        .controls()
        .map(|p| point(p, errors, options, dec).map(|p| point_from(&p)))
        .collect::<Result<Vec<_>, _>>()?;
    let weights = source
        .weights()
        .map(|w| weight(w, errors, options, dec))
        .collect::<Result<Vec<_>, _>>()?;
    BSplineCurve::new(source.degree(), knots, controls, weights, dec)
}

fn finite_uv_curve(
    source: &BSplineCurve,
    parameters: &[Rational; 2],
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let mut interval = parameters.clone();
    interval.sort();
    if source.domain() == interval.each_ref() {
        return Ok(source.clone());
    }
    let degree = source.degree();
    let mut knots = Vec::new();
    let mut controls = Vec::new();
    let mut weights = Vec::new();
    for (span, domain) in source.spans() {
        dec.charge_work(1)?;
        let lo = domain[0].clone().max(interval[0].clone());
        let hi = domain[1].clone().min(interval[1].clone());
        if lo >= hi {
            continue;
        }
        let net = source.span_bezier_controls_on(span, [&lo, &hi], dec)?;
        let first = controls.is_empty();
        knots.extend(std::iter::repeat_n(
            lo,
            if first { degree + 1 } else { degree },
        ));
        for [x, y, _, w] in net.into_iter().skip(usize::from(!first)) {
            dec.charge_work(1)?;
            controls.push(Point3::new(&x / &w, &y / &w, 0.into()));
            weights.push(w);
        }
    }
    if controls.is_empty() {
        return Err(bad(
            "rounding-topology",
            "finite UV interval has no retained span",
        ));
    }
    knots.extend(std::iter::repeat_n(interval[1].clone(), degree + 1));
    BSplineCurve::new(degree, knots, controls, weights, dec)
}

pub(super) fn uv_curve(
    source: &BSplineCurve,
    parameters: &[Rational; 2],
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let source = finite_uv_curve(source, parameters, dec)?;
    let knots = knots(source.knots(), errors, options, dec)?;
    let controls = source
        .controls()
        .map(|p| {
            let x = rounded(
                &p.x,
                &mut errors.spline.max_parameter_error,
                options,
                false,
                dec,
            )?;
            let y = rounded(
                &p.y,
                &mut errors.spline.max_parameter_error,
                options,
                false,
                dec,
            )?;
            Ok(Point3::new(x.into(), y.into(), 0.into()))
        })
        .collect::<Result<Vec<_>, Refusal>>()?;
    let weights = source
        .weights()
        .map(|w| weight(w, errors, options, dec))
        .collect::<Result<Vec<_>, _>>()?;
    BSplineCurve::new(source.degree(), knots, controls, weights, dec)
}

pub(super) fn surface(
    source: &BSplineSurface,
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<BSplineSurface, Refusal> {
    let knots = source
        .knots()
        .map(|row| knots(row, errors, options, dec))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .expect("two knot axes");
    let controls = source
        .controls()
        .iter()
        .map(|p| point(p, errors, options, dec).map(|p| point_from(&p)))
        .collect::<Result<Vec<_>, _>>()?;
    let weights = source
        .weights()
        .iter()
        .map(|w| weight(w, errors, options, dec))
        .collect::<Result<Vec<_>, _>>()?;
    BSplineSurface::new(
        source.degrees(),
        knots,
        source.control_count(),
        controls,
        weights,
        dec,
    )
}

/// The written atlas face: the rounded source, partitioned again, keeping
/// the component whose charts have the original domains after rounding.
pub(super) fn atlas(
    source: &RationalAtlas,
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<RationalAtlas, Refusal> {
    let carrier = surface(source.source(), errors, options, dec)?;
    let mut wanted = Vec::with_capacity(source.charts().len());
    for chart in source.charts() {
        let mut domain = chart.patch().domain().clone();
        for value in domain.iter_mut().flatten() {
            dec.charge_work(1)?;
            *value = written(value, options)?;
        }
        wanted.push(domain);
    }
    let periodic = [source.periodic_u(), source.periodic_v()];
    for candidate in RationalAtlas::partition(carrier, periodic, dec)? {
        dec.charge_work(1)?;
        if candidate
            .charts()
            .iter()
            .map(|chart| chart.patch().domain())
            .eq(wanted.iter())
        {
            return Ok(candidate);
        }
    }
    Err(bad(
        "rounding-topology",
        "the written atlas source no longer has the original charts",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    #[test]
    fn a_finite_uv_basis_keeps_native_parameters_across_partial_spans() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let source = BSplineCurve::new(
                2,
                [0, 0, 0, 1, 1, 2, 2, 2].map(Rational::from_i64).to_vec(),
                [[0, 0], [1, 1], [2, -1], [3, 2], [4, 0]]
                    .map(|[x, y]| Point3::new(x.into(), y.into(), Real::from(2).sqrt()))
                    .to_vec(),
                [1, 2, 3, 2, 1].map(Real::from).to_vec(),
                &mut dec,
            )
            .unwrap();
            let interval = [
                Rational::from_ratio(1, 4).unwrap(),
                Rational::from_ratio(7, 4).unwrap(),
            ];
            for reversed in [false, true] {
                let mut authored = interval.clone();
                if reversed {
                    authored.reverse();
                }
                let finite = finite_uv_curve(&source, &authored, &mut dec).unwrap();
                assert_eq!(finite.domain(), interval.each_ref());
                assert_eq!(finite.degree(), source.degree());
                assert_eq!(
                    finite
                        .knots()
                        .iter()
                        .filter(|t| **t == Rational::one())
                        .count(),
                    2
                );
                for i in 0..=8 {
                    let share = Rational::from_ratio(i, 8).unwrap();
                    let t = Real::from(&interval[0] + &(&share * &(&interval[1] - &interval[0])));
                    let a = source.point_at(&t, &mut dec).unwrap();
                    let b = finite.point_at(&t, &mut dec).unwrap();
                    assert_eq!(a.x.as_rational(), b.x.as_rational());
                    assert_eq!(a.y.as_rational(), b.y.as_rational());
                    assert_eq!(b.z.as_rational(), Some(Rational::zero()));
                }
                assert!(finite.point_at(&Real::zero(), &mut dec).is_err());
                assert!(finite.point_at(&Real::from(2), &mut dec).is_err());
            }
        });
        assert!(!work.exhausted, "{work:?}");
    }
}
