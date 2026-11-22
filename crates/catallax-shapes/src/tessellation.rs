//! Whole-triangle interpolation bounds for parametric display charts.

use catallax_numbers::{Interval, Rational};

mod rational;
pub(crate) use rational::rational_span_range;
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

/// Whole-triangle chord bound using the smaller of the box and metric bounds.
///
/// The three oriented edge vectors and the box widths use the same chart
/// coordinates as the derivative bounds. For widths u, v > 0, the diagonal
/// metric (Huu + Huv v/u, Hvv + Huv u/v) dominates every directional second
/// derivative by Young's inequality. Shewchuk's longest-edge bound in this
/// metric is the largest metric length squared divided by six.
pub fn triangle_chord_bound_interval(
    edges: [[Interval; 2]; 3],
    [u, v]: [Interval; 2],
    hessian: &[Interval; 3],
) -> Interval {
    let box_bound = chord_bound_interval([u, v], hessian);
    if u.lo() <= 0.0 || v.lo() <= 0.0 {
        return box_bound;
    }
    let u_metric = hessian[0] + hessian[1] * v / u;
    let v_metric = hessian[2] + hessian[1] * u / v;
    let longest = edges
        .into_iter()
        .map(|[du, dv]| u_metric * du.square() + v_metric * dv.square())
        .map(|length| length.hi())
        .fold(0.0_f64, f64::max);
    let metric_bound =
        Interval::point(longest).map(|length| length / Interval::point(6.0).unwrap());
    metric_bound
        .filter(|metric| metric.hi() < box_bound.hi())
        .unwrap_or(box_bound)
}

/// Checks the two-sided triangle bound after replacing chart vertices.
///
/// The triangle uses convex interpolation of the vertices, so the largest
/// replacement distance bounds its displacement at every barycentric point.
/// An already-enclosed combined bound can use a zero replacement distance.
pub fn checked_triangle_bound(
    interpolation: &Rational,
    replacement: &Rational,
    tolerance: &Rational,
) -> Option<Rational> {
    let bound = interpolation + replacement;
    (bound <= *tolerance).then_some(bound)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn metric_bound_uses_the_actual_triangle_edges() {
        // For S(u,v) = (u² + v²)/2 the triangle with vertices
        // (0,0), (1,1/2), (1/2,1) has its worst error at (5/12,5/12).
        let point = |x| Interval::point(x).unwrap();
        let hessian = [1.0, 0.0, 1.0].map(point);
        let widths = [point(1.0); 2];
        let edges = [[1.0, 0.5], [-0.5, 0.5], [-0.5, -1.0]].map(|e| e.map(point));
        let bound = triangle_chord_bound_interval(edges, widths, &hessian);
        assert!(bound.hi() >= 25.0 / 144.0);
        assert!(bound.hi() < chord_bound_interval(widths, &hessian).lo());
    }
}
