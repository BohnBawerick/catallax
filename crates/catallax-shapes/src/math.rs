//! Shared arithmetic for frames and polynomial definitions.

use crate::{RPoint3, RVector3};
use catallax_numbers::{Ball, MPoly, Rational, Refusal, Sign};

pub(crate) fn directions(axis: &RVector3, reference: &RVector3) -> Result<(), Refusal> {
    if axis.is_zero() || reference.is_zero() {
        Err(Refusal::new("shapes.zero-vector", "direction is zero"))
    } else {
        Ok(())
    }
}

pub(crate) fn perpendicular(axis: &RVector3, reference: &RVector3) -> Result<(), Refusal> {
    if !axis.dot(reference).is_zero() {
        Err(Refusal::new(
            "shapes.not-perpendicular",
            "reference direction is not perpendicular",
        ))
    } else {
        Ok(())
    }
}

pub(crate) fn positive(value: &Rational, code: &'static str) -> Result<(), Refusal> {
    if value.sign() != Sign::Positive {
        Err(Refusal::new(code, "parameter must be positive"))
    } else {
        Ok(())
    }
}

pub(crate) fn displacement(origin: &RPoint3) -> [MPoly; 3] {
    let coordinates = [&origin.x, &origin.y, &origin.z];
    std::array::from_fn(|i| MPoly::var(3, i) - MPoly::constant(3, coordinates[i].clone()))
}

pub(crate) fn dot(p: &[MPoly; 3], v: &RVector3) -> MPoly {
    p[0].scale(&v.x) + p[1].scale(&v.y) + p[2].scale(&v.z)
}

pub(crate) fn norm(p: &[MPoly; 3]) -> MPoly {
    &p[0] * &p[0] + &p[1] * &p[1] + &p[2] * &p[2]
}

pub(crate) fn cross(p: &[MPoly; 3], v: &RVector3) -> [MPoly; 3] {
    [
        p[1].scale(&v.z) - p[2].scale(&v.y),
        p[2].scale(&v.x) - p[0].scale(&v.z),
        p[0].scale(&v.y) - p[1].scale(&v.x),
    ]
}

pub(crate) fn ball(r: &Rational, precision: u32) -> Ball {
    Ball::from_rational(r, precision)
}

pub(crate) fn point(p: &RPoint3, precision: u32) -> [Ball; 3] {
    [&p.x, &p.y, &p.z].map(|r| ball(r, precision))
}

fn unit(v: &RVector3, precision: u32) -> [Ball; 3] {
    let n = v.norm_squared();
    [&v.x, &v.y, &v.z].map(|c| {
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

pub(crate) fn frame(axis: &RVector3, reference: &RVector3, precision: u32) -> [[Ball; 3]; 3] {
    [
        unit(reference, precision),
        unit(&axis.cross(reference), precision),
        unit(axis, precision),
    ]
}

pub(crate) fn add(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| &a[i] + &b[i])
}

pub(crate) fn scale(a: &[Ball; 3], c: &Ball) -> [Ball; 3] {
    std::array::from_fn(|i| &a[i] * c)
}

// Exact affine bounds before rounding avoid losing small local coordinates
// when a surface is far from the origin.
pub(crate) fn project(p: &[Ball; 3], origin: &RPoint3, dir: &RVector3, precision: u32) -> Ball {
    let mut lo = Rational::zero();
    let mut hi = Rational::zero();
    for ((b, o), c) in p
        .iter()
        .zip([&origin.x, &origin.y, &origin.z])
        .zip([&dir.x, &dir.y, &dir.z])
    {
        let a = &(&b.lo() - o) * c;
        let z = &(&b.hi() - o) * c;
        lo = lo + a.clone().min(z.clone());
        hi = hi + a.max(z);
    }
    let numerator = ball(&lo, precision).hull(&ball(&hi, precision));
    let length = ball(&dir.norm_squared(), precision)
        .sqrt()
        .expect("positive squared length");
    numerator
        .checked_div(&length)
        .expect("positive length enclosure")
}

pub(crate) fn sqrt_upper(r: &Rational) -> Rational {
    // Outward ball rounding proves this bound. No decision or policy is made.
    ball(r, 64).sqrt().expect("nonnegative bound").hi()
}

pub(crate) fn inverse_length_upper(v: &RVector3) -> Rational {
    sqrt_upper(&v.norm_squared().recip().expect("nonzero direction"))
}
