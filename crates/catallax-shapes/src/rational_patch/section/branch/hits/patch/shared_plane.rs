//! Strict rational control-halfspace ownership of one shared Bezier side.

use super::*;

fn coordinates(point: &Point3) -> Option<[&Rational; 3]> {
    Some([
        point.x.rational_ref()?,
        point.y.rational_ref()?,
        point.z.rational_ref()?,
    ])
}

fn sub(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, false)?;
    Ok(a - b)
}

fn add(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, false)?;
    Ok(a + b)
}

fn mul(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, true)?;
    Ok(a * b)
}

fn offset(
    point: &[&Rational; 3],
    origin: &[&Rational; 3],
    dec: &mut Decider,
) -> Result<[Rational; 3], Refusal> {
    Ok([
        sub(point[0], origin[0], dec)?,
        sub(point[1], origin[1], dec)?,
        sub(point[2], origin[2], dec)?,
    ])
}

fn product_difference(
    a: &Rational,
    b: &Rational,
    c: &Rational,
    d: &Rational,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    sub(&mul(a, b, dec)?, &mul(c, d, dec)?, dec)
}

fn boundary_normal(
    boundary: &BSplineCurve,
    origin: &[&Rational; 3],
    dec: &mut Decider,
) -> Result<Option<[Rational; 3]>, Refusal> {
    for i in 1..boundary.control_count() {
        let Some(point) = coordinates(boundary.control(i)) else {
            return Ok(None);
        };
        let a = offset(&point, origin, dec)?;
        for j in i + 1..boundary.control_count() {
            let Some(point) = coordinates(boundary.control(j)) else {
                return Ok(None);
            };
            let b = offset(&point, origin, dec)?;
            let normal = [
                product_difference(&a[1], &b[2], &a[2], &b[1], dec)?,
                product_difference(&a[2], &b[0], &a[0], &b[2], dec)?,
                product_difference(&a[0], &b[1], &a[1], &b[0], dec)?,
            ];
            if normal.iter().any(|value| !value.is_zero()) {
                return Ok(Some(normal));
            }
        }
    }
    Ok(None)
}

fn distance(
    point: &[&Rational; 3],
    origin: &[&Rational; 3],
    normal: &[Rational; 3],
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let delta = offset(point, origin, dec)?;
    let x = mul(&normal[0], &delta[0], dec)?;
    let y = mul(&normal[1], &delta[1], dec)?;
    let z = mul(&normal[2], &delta[2], dec)?;
    add(&add(&x, &y, dec)?, &z, dec)
}

fn halfspace(
    patch: &RationalPatch,
    axis: usize,
    high: bool,
    origin: &[&Rational; 3],
    normal: &[Rational; 3],
    dec: &mut Decider,
) -> Result<Option<Sign>, Refusal> {
    let count = patch.carrier().control_count();
    let fixed = if high { count[axis] - 1 } else { 0 };
    let mut sign = None;
    for (index, point) in patch.carrier().controls().iter().enumerate() {
        dec.charge_work(1)?;
        let Some(point) = coordinates(point) else {
            return Ok(None);
        };
        let next = distance(&point, origin, normal, dec)?.sign();
        let coordinate = if axis == 0 {
            index / count[1]
        } else {
            index % count[1]
        };
        if coordinate == fixed {
            if next != Sign::Zero {
                return Ok(None);
            }
        } else {
            if next == Sign::Zero || sign.is_some_and(|old| old != next) {
                return Ok(None);
            }
            sign = Some(next);
        }
    }
    Ok(sign)
}

/// The caller has already matched the two positive-weight Bezier side definitions.
pub(super) fn separates(
    source: &RationalPatch,
    target: &RationalPatch,
    side: SplineBoundary,
    axis: usize,
    high: bool,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let boundary = source.boundary(side);
    let Some(origin) = coordinates(boundary.control(0)) else {
        return Ok(false);
    };
    let Some(normal) = boundary_normal(&boundary, &origin, dec)? else {
        return Ok(false);
    };
    let Some(a) = halfspace(source, axis, high, &origin, &normal, dec)? else {
        return Ok(false);
    };
    let Some(b) = halfspace(target, axis, !high, &origin, &normal, dec)? else {
        return Ok(false);
    };
    // Positive Bezier bases make the plane's zero set exactly the named side.
    // Opposite strict halfspaces therefore meet only on their matched side.
    Ok(a != b)
}
