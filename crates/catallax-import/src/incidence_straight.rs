//! Exact regular line images of repaired rational splines.

use catallax_numbers::{Decider, Rational, Refusal};

use crate::exact::{cross, dot, is_zero, sub};
use crate::{Curve, OrientedEdge, Vec3};

fn failure(error: Refusal) -> String {
    format!(
        "straight image certificate: {}: {}",
        error.code, error.detail
    )
}

/// A nonnegative Bernstein derivative numerator, with a positive coefficient
/// on every span, is strictly positive inside each span. Exact joins make the
/// complete curve a continuous, strictly ordered line image even when a raw
/// parameter derivative vanishes at an endpoint.
pub(crate) fn reduce(edge: &OrientedEdge, dec: &mut Decider) -> Result<Option<Curve>, String> {
    let Curve::BSpline(source) = &edge.curve else {
        return Ok(None);
    };
    let controls = source.geometry.controls();
    let origin = &controls[0];
    let direction = sub(&controls[controls.len() - 1], origin);
    if is_zero(&direction)
        || controls
            .iter()
            .any(|p| !is_zero(&cross(&sub(p, origin), &direction)))
    {
        return Ok(None);
    }
    let curve = source.geometry.to_curve(dec).map_err(failure)?;
    let mut first = None;
    let mut previous: Option<Vec3> = None;
    for (span, _) in curve.spans() {
        dec.check_deadline().map_err(failure)?;
        let net = curve
            .span_bezier_controls(span, dec)
            .map_err(failure)?
            .into_iter()
            .map(|control| {
                let mut point = std::array::from_fn(|_| Rational::zero());
                for i in 0..4 {
                    point[i] = control[i]
                        .as_rational()
                        .ok_or("straight image needs rational homogeneous controls")?;
                }
                Ok(point)
            })
            .collect::<Result<Vec<[Rational; 4]>, &str>>()
            .map_err(str::to_owned)?;
        if net.iter().any(|p| p[3] <= Rational::zero()) {
            return Ok(None);
        }
        let a: Vec3 = std::array::from_fn(|i| &net[0][i] / &net[0][3]);
        let end = &net[net.len() - 1];
        let b: Vec3 = std::array::from_fn(|i| &end[i] / &end[3]);
        if previous.as_ref().is_some_and(|p| p != &a) {
            return Ok(None);
        }
        first.get_or_insert(a);
        previous = Some(b);
        let degree = net.len() - 1;
        if degree == 0 {
            return Ok(None);
        }
        let scalar = net
            .iter()
            .map(|p| {
                let v = std::array::from_fn(|i| &p[i] - &(&origin[i] * &p[3]));
                dot(&v, &direction)
            })
            .collect::<Vec<_>>();
        let scale = Rational::from_i64(degree as i64);
        let dn = scalar
            .windows(2)
            .map(|p| (&p[1] - &p[0]) * scale.clone())
            .collect::<Vec<_>>();
        let dw = net
            .windows(2)
            .map(|p| (&p[1][3] - &p[0][3]) * scale.clone())
            .collect::<Vec<_>>();
        let mut positive = false;
        for k in 0..=2 * degree - 1 {
            dec.check_deadline().map_err(failure)?;
            let mut coefficient = Rational::zero();
            for (i, (dn, dw)) in dn.iter().zip(&dw).enumerate() {
                dec.check_deadline().map_err(failure)?;
                if i <= k && k - i <= degree {
                    let j = k - i;
                    let factor = crate::incidence::choose(degree - 1, i)
                        * crate::incidence::choose(degree, j)
                        / crate::incidence::choose(2 * degree - 1, k);
                    coefficient = coefficient + (&net[j][3] * dn - &scalar[j] * dw) * factor;
                }
            }
            if coefficient < Rational::zero() {
                return Ok(None);
            }
            positive |= coefficient > Rational::zero();
        }
        if !positive {
            return Ok(None);
        }
    }
    let (Some(first), Some(last)) = (first, previous) else {
        return Ok(None);
    };
    let lo = dot(&sub(&first, origin), &direction);
    let hi = dot(&sub(&last, origin), &direction);
    if lo >= hi {
        return Ok(None);
    }
    let mut parameters = Vec::new();
    for endpoint in [&edge.start.point, &edge.end.point] {
        let relative = sub(endpoint, origin);
        let parameter = dot(&relative, &direction);
        if !is_zero(&cross(&relative, &direction)) || parameter < lo || parameter > hi {
            return Ok(None);
        }
        parameters.push(parameter);
    }
    let ordered = if edge.curve_same_sense {
        parameters[0] < parameters[1]
    } else {
        parameters[0] > parameters[1]
    };
    Ok(ordered.then(|| Curve::Line(origin.clone(), direction)))
}
