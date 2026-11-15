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
