//! Shared certified arithmetic for frames and polynomial definitions.

use crate::{Point3, Vector3};
use catallax_numbers::{
    Ball, Decider, Integer, MPoly, Outcome, Rational, Real, RealPoly, Refusal, Sign, StrictSign,
};

pub(crate) fn exact_sqrt(value: &Real) -> Real {
    fn integer_sqrt(n: &Integer) -> Option<Integer> {
        if n.sign() == Sign::Negative {
            return None;
        }
        if n.is_zero() {
            return Some(Integer::zero());
        }
        let mut x = Integer::one().mul_pow2(n.bit_length().div_ceil(2));
        loop {
            let y = (&x + &n.div_rem_floor(&x)?.0)
                .div_rem_floor(&Integer::from_i64(2))?
                .0;
            if y >= x {
                return (&x * &x == *n).then_some(x);
            }
            x = y;
        }
    }
    let rational = value
        .rational_ref()
        .and_then(|q| Rational::new(integer_sqrt(q.numerator())?, integer_sqrt(q.denominator())?));
    rational.map_or_else(|| value.sqrt(), Real::from)
}

pub(crate) fn validate(value: &Real, dec: &mut Decider) -> Result<(), Refusal> {
    if value.as_rational().is_none() {
        dec.enclose(value, &Rational::one())?;
    }
    Ok(())
}

pub(crate) fn validate_point(point: &Point3, dec: &mut Decider) -> Result<(), Refusal> {
    for value in [&point.x, &point.y, &point.z] {
        validate(value, dec)?;
    }
    Ok(())
}

pub(crate) fn validate_vector(vector: &Vector3, dec: &mut Decider) -> Result<(), Refusal> {
    for value in [&vector.x, &vector.y, &vector.z] {
        validate(value, dec)?;
    }
    Ok(())
}

pub(crate) fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    if let Some(value) = value.as_rational() {
        return Ok(value.sign());
    }
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(r) => Err(r),
    }
}

pub(crate) fn directions(
    axis: &Vector3,
    reference: &Vector3,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for vector in [axis, reference] {
        if sign(&vector.norm_squared(), dec)? == Sign::Zero {
            return Err(Refusal::new("shapes.zero-vector", "direction is zero"));
        }
    }
    Ok(())
}

pub(crate) fn perpendicular(
    axis: &Vector3,
    reference: &Vector3,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if sign(&axis.dot(reference), dec)? != Sign::Zero {
        return Err(Refusal::new(
            "shapes.not-perpendicular",
            "reference direction is not perpendicular",
        ));
    }
    Ok(())
}

pub(crate) fn positive(value: &Real, code: &'static str, dec: &mut Decider) -> Result<(), Refusal> {
    if sign(value, dec)? != Sign::Positive {
        return Err(Refusal::new(code, "parameter must be positive"));
    }
    Ok(())
}

pub(crate) fn displacement(origin: &Point3) -> [RealPoly; 3] {
    let coordinates = [&origin.x, &origin.y, &origin.z];
    if let Some(origin) = origin.as_rational() {
        let coordinates = [origin.x, origin.y, origin.z];
        return std::array::from_fn(|i| {
            RealPoly::from(MPoly::var(3, i) - MPoly::constant(3, coordinates[i].clone()))
        });
    }
    std::array::from_fn(|i| RealPoly::var(3, i) - RealPoly::constant(3, coordinates[i].clone()))
}

pub(crate) fn dot(p: &[RealPoly; 3], v: &Vector3) -> RealPoly {
    p[0].scale(&v.x) + p[1].scale(&v.y) + p[2].scale(&v.z)
}

pub(crate) fn norm(p: &[RealPoly; 3]) -> RealPoly {
    &p[0] * &p[0] + &p[1] * &p[1] + &p[2] * &p[2]
}

pub(crate) fn cross(p: &[RealPoly; 3], v: &Vector3) -> [RealPoly; 3] {
    [
        p[1].scale(&v.z) - p[2].scale(&v.y),
        p[2].scale(&v.x) - p[0].scale(&v.z),
        p[0].scale(&v.y) - p[1].scale(&v.x),
    ]
}

pub(crate) fn ball(r: &Rational, precision: u32) -> Ball {
    Ball::from_rational(r, precision)
}

pub(crate) fn point(p: &Point3, precision: u32, dec: &Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        dec.enclose_at(&p.x, precision)?,
        dec.enclose_at(&p.y, precision)?,
        dec.enclose_at(&p.z, precision)?,
    ])
}

fn rational_squared_norm(v: &[&Rational; 3]) -> Rational {
    v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
}

fn rational_unit(v: [&Rational; 3], precision: u32) -> [Ball; 3] {
    let n = rational_squared_norm(&v);
    v.map(|c| {
        let magnitude = ball(&(&(c * c) / &n), precision)
            .sqrt()
            .expect("nonnegative squared component");
        if c.sign() == Sign::Negative {
            -magnitude
        } else {
            magnitude
        }
    })
}

fn unit(v: &Vector3, precision: u32, dec: &Decider) -> Result<[Ball; 3], Refusal> {
    if let Some(rational) = v.rational_components() {
        return Ok(rational_unit(rational, precision));
    }
    let length = v.norm_squared().sqrt();
    Ok([
        dec.enclose_at(&(&v.x / &length), precision)?,
        dec.enclose_at(&(&v.y / &length), precision)?,
        dec.enclose_at(&(&v.z / &length), precision)?,
    ])
}

pub(crate) fn frame(
    axis: &Vector3,
    reference: &Vector3,
    precision: u32,
    dec: &Decider,
) -> Result<[[Ball; 3]; 3], Refusal> {
    Ok([
        unit(reference, precision, dec)?,
        unit(&axis.cross(reference), precision, dec)?,
        unit(axis, precision, dec)?,
    ])
}

pub(crate) fn add(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| &a[i] + &b[i])
}
pub(crate) fn scale(a: &[Ball; 3], c: &Ball) -> [Ball; 3] {
    std::array::from_fn(|i| &a[i] * c)
}

// Exact affine bounds preserve tiny local coordinates at large rational offsets.
fn rational_project(
    p: &[Ball; 3],
    origin: &[&Rational; 3],
    dir: &[&Rational; 3],
    precision: u32,
) -> Ball {
    let mut lo = Rational::zero();
    let mut hi = Rational::zero();
    for ((b, o), c) in p
        .iter()
        .zip(origin.iter().copied())
        .zip(dir.iter().copied())
    {
        let a = &(&b.lo() - o) * c;
        let z = &(&b.hi() - o) * c;
        let (lower, upper) = if a <= z { (a, z) } else { (z, a) };
        lo = lo + lower;
        hi = hi + upper;
    }
    let numerator = ball(&lo, precision).hull(&ball(&hi, precision));
    let length = ball(&rational_squared_norm(dir), precision)
        .sqrt()
        .expect("positive squared length");
    numerator
        .checked_div(&length)
        .expect("positive length enclosure")
}

pub(crate) fn project(
    p: &[Ball; 3],
    origin: &Point3,
    dir: &Vector3,
    precision: u32,
    dec: &Decider,
) -> Result<Ball, Refusal> {
    if let (Some(origin), Some(dir)) = (origin.rational_components(), dir.rational_components()) {
        return Ok(rational_project(p, &origin, &dir, precision));
    }
    let origin = point(origin, precision, dec)?;
    let dir = unit(dir, precision, dec)?;
    Ok((0..3).fold(ball(&Rational::zero(), precision), |sum, i| {
        &sum + &(&(&p[i] - &origin[i]) * &dir[i])
    }))
}

pub(crate) fn sqrt_upper(r: &Rational) -> Rational {
    ball(r, 64).sqrt().expect("nonnegative bound").hi()
}

pub(crate) fn rational_norm_squared(v: &Vector3) -> Option<Rational> {
    v.rational_components().map(|v| rational_squared_norm(&v))
}

pub(crate) fn inverse_length_upper(v: &Vector3, dec: &Decider) -> Result<Rational, Refusal> {
    if let Some(norm) = rational_norm_squared(v) {
        return Ok(sqrt_upper(&norm.recip().expect("nonzero direction")));
    }
    let inverse = &Real::one() / &v.norm_squared().sqrt();
    Ok(dec.enclose_at(&inverse, 64)?.hi())
}
