//! Binary64 helpers shared by the fitted contour generators.

use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::Vector3;

/// A binary64 value within rule T1's floor of an exact one: the nearest
/// double of a rational, otherwise the middle of an enclosure that narrow.
pub(super) fn double(value: &Real, floor: &Rational, dec: &mut Decider) -> Result<f64, Refusal> {
    if let Some(rational) = value.rational_ref() {
        return Ok(rational.to_f64_nearest());
    }
    Ok(dec.enclose(value, floor)?.mid().to_f64_nearest())
}

pub(super) fn doubles(
    v: &Vector3,
    floor: &Rational,
    dec: &mut Decider,
) -> Result<[f64; 3], Refusal> {
    Ok([
        double(&v.x, floor, dec)?,
        double(&v.y, floor, dec)?,
        double(&v.z, floor, dec)?,
    ])
}

/// The exact value of a binary64 result, refusing a nonfinite one.
pub(super) fn rational(value: f64) -> Result<Rational, Refusal> {
    Rational::from_f64_exact(value).ok_or_else(|| {
        Refusal::new(
            "app.contour-fit",
            "a fitted contour generator has a nonfinite value",
        )
    })
}

/// `a > b`; false when either is NaN, so a negated test also catches NaN.
pub(super) fn above(a: f64, b: f64) -> bool {
    a > b
}

pub(super) fn unit(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}

pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// The unit binary64 direction of a nonzero exact vector at any scale.
///
/// A rational vector is first multiplied by an exact power of two that
/// brings its largest component into `[1/2, 2)`, so that neither rounding
/// to binary64 nor squaring underflows or overflows: `1e-200 (1, -1, 1)`
/// would otherwise square to zero. IEEE 754 rounds every basic operation and
/// the square root correctly, and a power-of-two factor changes no
/// significand while results stay normal, so the scale commutes with each
/// step: a vector whose components and squares are normal binary64 numbers
/// gives the same bits as before. Other components keep [`double`].
pub(super) fn unit_direction(
    v: &Vector3,
    floor: &Rational,
    dec: &mut Decider,
) -> Result<[f64; 3], Refusal> {
    let parts = [&v.x, &v.y, &v.z].map(Real::rational_ref);
    let [Some(x), Some(y), Some(z)] = parts else {
        return Ok(unit(doubles(v, floor, dec)?));
    };
    // The bit-length difference of a nonzero rational's numerator and
    // denominator is its binary exponent, or one more.
    let exponent = [x, y, z]
        .into_iter()
        .filter(|q| !q.is_zero())
        .map(|q| q.numerator().bit_length() as i64 - q.denominator().bit_length() as i64)
        .max()
        .ok_or_else(|| Refusal::new("app.contour-direction", "view direction is zero"))?;
    Ok(unit(
        [x, y, z].map(|q| q.mul_pow2(-exponent).to_f64_nearest()),
    ))
}
