//! Binary64 helpers shared by the fitted contour generators, the
//! hidden-line views and the drawing curves.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use catallax_numbers::{Decider, Integer, Rational, Real, Refusal, Sign};
use catallax_shapes::{BSplineCurve, Point3, Vector3};

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

/// The square root of an exact square rational stays rational; anything
/// else becomes a lazy root.
pub(super) fn root(value: &Real) -> Real {
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
    value
        .rational_ref()
        .and_then(|q| Rational::new(integer_sqrt(q.numerator())?, integer_sqrt(q.denominator())?))
        .map_or_else(|| value.sqrt(), Real::from)
}

/// The exact point at binary64 angle `t` of the conic
/// `center + cos t u + sin t w`. The half-angle form keeps it exactly on the
/// conic with short numbers: `cos t = (1 - s^2) / (1 + s^2)` and
/// `sin t = 2 s / (1 + s^2)` for the rational `s` nearest `tan(t / 2)`
/// (Piegl and Tiller, "The NURBS Book", 2nd ed., section 7.2, the rational
/// parametrisation of the circle). A half turn first brings `t` within a
/// quarter turn of zero, so that `|s| <= 1`.
pub(super) fn conic_point(
    center: &Point3,
    [u, w]: [&Vector3; 2],
    t: f64,
) -> Result<Point3, Refusal> {
    let mut angle = t.rem_euclid(TAU);
    let mut flip = Real::one();
    if angle > PI {
        angle -= TAU;
    }
    if angle.abs() > FRAC_PI_2 {
        angle -= PI.copysign(angle);
        flip = -Real::one();
    }
    let s = Real::from(rational((angle / 2.0).tan())?);
    let one = Real::one();
    let den = &one + &s.square();
    let cos = &(&(&one - &s.square()) / &den) * &flip;
    let sin = &(&(&Real::from_i64(2) * &s) / &den) * &flip;
    Ok(center.plus(&u.scaled(&cos).plus(&w.scaled(&sin))))
}

/// Binary64 knots of an exact clamped B-spline, moved exactly so that its
/// domain starts at zero; `None` when two distinct knots still round to one
/// binary64 value, or a knot overflows.
///
/// A binary64 number keeps 53 significant bits of its own size, not of the
/// knots' spacing: on the exact domain `[10^16, 10^16 + 1]` both ends round
/// to `10^16`, and every span would be lost. The exact shift leaves only the
/// knots' distances from the domain's start to round. The exact parameter
/// at a binary64 one is the domain's start plus it ([`native_parameter`]),
/// and an exact affine shift of a curve's parameter changes nothing here.
pub(super) fn shifted_knots(spline: &BSplineCurve) -> Option<Vec<f64>> {
    let exact = spline.knots();
    let start = spline.domain()[0];
    let knots: Vec<f64> = exact.iter().map(|k| (k - start).to_f64_nearest()).collect();
    let merged = (1..knots.len()).any(|i| exact[i] > exact[i - 1] && knots[i] <= knots[i - 1]);
    (knots.iter().all(|k| k.is_finite()) && !merged).then_some(knots)
}

/// The exact parameter of `spline` at the binary64 parameter `t` of its
/// [`shifted_knots`], kept inside the exact domain, which rounded knots can
/// reach past.
pub(super) fn native_parameter(spline: &BSplineCurve, t: f64) -> Result<Rational, Refusal> {
    let [lo, hi] = spline.domain();
    Ok((lo + &rational(t)?).max(lo.clone()).min(hi.clone()))
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
