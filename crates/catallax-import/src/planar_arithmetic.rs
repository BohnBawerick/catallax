//! Planar repair uses the same reserved exact arithmetic as spline repair.

use crate::incidence_arithmetic as a;
use crate::{ImportRefusal, Vec3};
use catallax_numbers::{Decider, Rational};

fn checked<T>(result: Result<T, String>) -> Result<T, ImportRefusal> {
    result.map_err(|_| ImportRefusal::ResourceLimit)
}

pub(crate) fn checkpoint(d: &mut Decider) -> Result<(), ImportRefusal> {
    checked(a::checkpoint(d))
}
pub(crate) fn arithmetic(
    x: &Rational,
    y: &Rational,
    op: char,
    d: &mut Decider,
) -> Result<Rational, ImportRefusal> {
    checked(a::arithmetic(x, y, op, d))
}
pub(crate) fn compare(
    x: &Rational,
    y: &Rational,
    d: &mut Decider,
) -> Result<std::cmp::Ordering, ImportRefusal> {
    checked(a::compare(x, y, d))
}
pub(crate) fn same(x: &Vec3, y: &Vec3, d: &mut Decider) -> Result<bool, ImportRefusal> {
    checked(a::same(x, y, d))
}
pub(crate) fn sub(x: &Vec3, y: &Vec3, d: &mut Decider) -> Result<Vec3, ImportRefusal> {
    checked(a::sub(x, y, d))
}
pub(crate) fn dot(x: &Vec3, y: &Vec3, d: &mut Decider) -> Result<Rational, ImportRefusal> {
    checked(a::dot(x, y, d))
}
pub(crate) fn cross(x: &Vec3, y: &Vec3, d: &mut Decider) -> Result<Vec3, ImportRefusal> {
    checked(a::cross(x, y, d))
}
pub(crate) fn maximum(
    x: Rational,
    y: Rational,
    d: &mut Decider,
) -> Result<Rational, ImportRefusal> {
    checked(a::maximum(x, y, d))
}
pub(crate) fn minimum(
    x: Rational,
    y: Rational,
    d: &mut Decider,
) -> Result<Rational, ImportRefusal> {
    checked(a::minimum(x, y, d))
}
pub(crate) fn sqrt_enclosure(
    x: &Rational,
    d: &mut Decider,
) -> Result<(Rational, Rational), ImportRefusal> {
    checked(a::sqrt_enclosure(x, d))
}
pub(crate) fn squared(x: &Vec3, y: &Vec3, d: &mut Decider) -> Result<Rational, ImportRefusal> {
    let delta = sub(x, y, d)?;
    dot(&delta, &delta, d)
}
pub(crate) fn offset(
    x: &Vec3,
    direction: &Vec3,
    scale: &Rational,
    op: char,
    d: &mut Decider,
) -> Result<Vec3, ImportRefusal> {
    let mut result = x.clone();
    for i in 0..3 {
        let product = arithmetic(scale, &direction[i], '*', d)?;
        result[i] = arithmetic(&x[i], &product, op, d)?;
    }
    Ok(result)
}
