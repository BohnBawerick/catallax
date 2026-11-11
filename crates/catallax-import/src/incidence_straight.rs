//! Exact regular line images of repaired rational splines.

use catallax_numbers::{Decider, Rational, Refusal, Sign};

use crate::exact::is_zero;
use crate::incidence_arithmetic::{
    arithmetic, binomial, checkpoint, compare, cross, dot, same, sub,
};
use crate::{Curve, OrientedEdge, Vec3};

fn failure(error: Refusal) -> String {
    format!(
        "straight image certificate: {}: {}",
        error.code, error.detail
    )
}

/// Nonnegative Bernstein derivative coefficients, with a positive coefficient
/// on each span, prove a strictly ordered line image. Parameter derivatives
/// may vanish at endpoints; constant spans and reversals remain refused.
pub(crate) fn reduce(edge: &OrientedEdge, dec: &mut Decider) -> Result<Option<Curve>, String> {
    checkpoint(dec)?;
    let Curve::BSpline(source) = &edge.curve else {
        return Ok(None);
    };
    let controls = source.geometry.controls();
    let origin = &controls[0];
    let direction = sub(&controls[controls.len() - 1], origin, dec)?;
    if is_zero(&direction) {
        return Ok(None);
    }
    for point in controls {
        checkpoint(dec)?;
        let offset = sub(point, origin, dec)?;
        if !is_zero(&cross(&offset, &direction, dec)?) {
            return Ok(None);
        }
    }
    let curve = source.geometry.to_curve(dec).map_err(failure)?;
    let mut first = None;
    let mut previous: Option<Vec3> = None;
    for (span, _) in curve.spans() {
        checkpoint(dec)?;
        let mut net: Vec<[Rational; 4]> = Vec::new();
        for control in curve.span_bezier_controls(span, dec).map_err(failure)? {
            checkpoint(dec)?;
            let mut point = std::array::from_fn(|_| Rational::zero());
            for i in 0..4 {
                point[i] = control[i]
                    .as_rational()
                    .ok_or("straight image needs rational homogeneous controls")?;
            }
            net.push(point);
        }
        for point in &net {
            checkpoint(dec)?;
            if point[3].sign() != Sign::Positive {
                return Ok(None);
            }
        }
        let a = [
            arithmetic(&net[0][0], &net[0][3], '/', dec)?,
            arithmetic(&net[0][1], &net[0][3], '/', dec)?,
            arithmetic(&net[0][2], &net[0][3], '/', dec)?,
        ];
        let end = &net[net.len() - 1];
        let b = [
            arithmetic(&end[0], &end[3], '/', dec)?,
            arithmetic(&end[1], &end[3], '/', dec)?,
            arithmetic(&end[2], &end[3], '/', dec)?,
        ];
        if let Some(previous) = &previous
            && !same(previous, &a, dec)?
        {
            return Ok(None);
        }
        first.get_or_insert(a);
        previous = Some(b);
        let degree = net.len() - 1;
        if degree == 0 {
            return Ok(None);
        }
        let mut scalar = Vec::new();
        for point in &net {
            checkpoint(dec)?;
            let mut v = std::array::from_fn(|_| Rational::zero());
            for i in 0..3 {
                let origin = arithmetic(&origin[i], &point[3], '*', dec)?;
                v[i] = arithmetic(&point[i], &origin, '-', dec)?;
            }
            scalar.push(dot(&v, &direction, dec)?);
        }
        let scale = Rational::from_i64(degree as i64);
        let mut dn = Vec::new();
        let mut dw = Vec::new();
        for i in 0..degree {
            checkpoint(dec)?;
            let n = arithmetic(&scalar[i + 1], &scalar[i], '-', dec)?;
            let w = arithmetic(&net[i + 1][3], &net[i][3], '-', dec)?;
            dn.push(arithmetic(&n, &scale, '*', dec)?);
            dw.push(arithmetic(&w, &scale, '*', dec)?);
        }
        let mut positive = false;
        for k in 0..=2 * degree - 1 {
            checkpoint(dec)?;
            let mut coefficient = Rational::zero();
            for (i, (dn, dw)) in dn.iter().zip(&dw).enumerate() {
                checkpoint(dec)?;
                if i <= k && k - i <= degree {
                    let j = k - i;
                    let a = binomial(degree - 1, i, dec)?;
                    let b = binomial(degree, j, dec)?;
                    let c = binomial(2 * degree - 1, k, dec)?;
                    let factor = arithmetic(&a, &b, '*', dec)?;
                    let factor = arithmetic(&factor, &c, '/', dec)?;
                    let a = arithmetic(&net[j][3], dn, '*', dec)?;
                    let b = arithmetic(&scalar[j], dw, '*', dec)?;
                    let delta = arithmetic(&a, &b, '-', dec)?;
                    let term = arithmetic(&delta, &factor, '*', dec)?;
                    coefficient = arithmetic(&coefficient, &term, '+', dec)?;
                }
            }
            if coefficient.sign() == Sign::Negative {
                return Ok(None);
            }
            positive |= coefficient.sign() == Sign::Positive;
        }
        if !positive {
            return Ok(None);
        }
    }
    let (Some(first), Some(last)) = (first, previous) else {
        return Ok(None);
    };
    let a = sub(&first, origin, dec)?;
    let b = sub(&last, origin, dec)?;
    let lo = dot(&a, &direction, dec)?;
    let hi = dot(&b, &direction, dec)?;
    if compare(&lo, &hi, dec)?.is_ge() {
        return Ok(None);
    }
    let mut parameters = Vec::new();
    for endpoint in [&edge.start.point, &edge.end.point] {
        checkpoint(dec)?;
        let relative = sub(endpoint, origin, dec)?;
        let parameter = dot(&relative, &direction, dec)?;
        if !is_zero(&cross(&relative, &direction, dec)?)
            || compare(&parameter, &lo, dec)?.is_lt()
            || compare(&parameter, &hi, dec)?.is_gt()
        {
            return Ok(None);
        }
        parameters.push(parameter);
    }
    let ordering = compare(&parameters[0], &parameters[1], dec)?;
    let ordered = if edge.curve_same_sense {
        ordering.is_lt()
    } else {
        ordering.is_gt()
    };
    Ok(ordered.then(|| Curve::Line(origin.clone(), direction)))
}
