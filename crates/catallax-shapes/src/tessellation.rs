//! Whole-triangle interpolation bounds for parametric display charts.

use catallax_numbers::{Interval, Rational};

mod rational;
pub use rational::{Projection, RationalSpan, rational_spans};

/// Filip, Magedson and Markot's second-derivative bound on a parameter box.
///
/// `hessian` bounds the magnitudes of `S_uu`, `S_uv` and `S_vv` throughout
/// the box. Each coordinate variance is at most its squared range / 4.
/// Cauchy-Schwarz gives the mixed term. A caller adds the largest vertex
/// residual because the triangle interpolates vertex errors convexly.
pub fn chord_bound([u, v]: [&Rational; 2], hessian: &[Rational; 3]) -> Rational {
    (&(&hessian[0] * u) * u + (&(&hessian[1] * u) * v).mul_pow2(1) + &(&hessian[2] * v) * v)
        .mul_pow2(-3)
}

/// Outward-rounded binary64 evaluation of [`chord_bound`].
pub fn chord_bound_interval([u, v]: [Interval; 2], hessian: &[Interval; 3]) -> Interval {
    (hessian[0] * u.square()
        + hessian[1] * u * v * Interval::point(2.0).unwrap()
        + hessian[2] * v.square())
        * Interval::point(0.125).unwrap()
}
