//! Binary64 images of view strokes: points, derivatives and parameters of
//! their carriers, bounds of their spans along a direction, and rational
//! Bezier pieces for the crossing search.
//!
//! Lines are `o + t v`. Circles and ellipses are `c + cos t u + sin t w`
//! with `u` perpendicular to `w`, the parametrisation of their exact
//! carriers. B-splines are cut into their Bezier spans by knot insertion
//! (Piegl and Tiller, "The NURBS Book", 2nd ed., algorithm A5.6) and keep
//! their knot parameter.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use super::super::float::{above, dot, unit};

pub(super) type V3 = [f64; 3];

fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

pub(super) fn distance(a: V3, b: V3) -> f64 {
    let d = sub(a, b);
    dot(d, d).sqrt()
}

/// One Bezier span of a B-spline, over the knot interval `[lo, hi]`, with
/// homogeneous controls `(w x, w y, w z, w)`.
#[derive(Clone, Debug)]
pub(super) struct Span {
    pub lo: f64,
    pub hi: f64,
    pub ctrl: Vec<[f64; 4]>,
}

/// De Casteljau's algorithm on homogeneous controls (Farin, "Curves and
/// Surfaces for CAGD", 5th ed., section 3.1): the value at `s` and the
/// derivative in `s`.
pub(super) fn casteljau(ctrl: &[[f64; 4]], s: f64) -> ([f64; 4], [f64; 4]) {
    let n = ctrl.len() - 1;
    if n == 0 {
        return (ctrl[0], [0.0; 4]);
    }
    let mut level = ctrl.to_vec();
    for round in 1..n {
        for i in 0..=n - round {
            for k in 0..4 {
                level[i][k] += s * (level[i + 1][k] - level[i][k]);
            }
        }
    }
    let (a, b) = (level[0], level[1]);
    let value = std::array::from_fn(|k| a[k] + s * (b[k] - a[k]));
    let derivative = std::array::from_fn(|k| n as f64 * (b[k] - a[k]));
    (value, derivative)
}

/// The two halves of homogeneous Bezier controls at `s`.
pub(super) fn split<const D: usize>(ctrl: &[[f64; D]], s: f64) -> (Vec<[f64; D]>, Vec<[f64; D]>) {
    let n = ctrl.len();
    let mut level = ctrl.to_vec();
    let mut left = Vec::with_capacity(n);
    let mut right = vec![[0.0; D]; n];
    left.push(level[0]);
    right[n - 1] = level[n - 1];
    for round in 1..n {
        for i in 0..n - round {
            for k in 0..D {
                level[i][k] += s * (level[i + 1][k] - level[i][k]);
            }
        }
        left.push(level[0]);
        right[n - 1 - round] = level[n - 1 - round];
    }
    (left, right)
}

/// Controls of the part of a Bezier span over local `[a, b]`.
fn restrict(ctrl: &[[f64; 4]], a: f64, b: f64) -> Vec<[f64; 4]> {
    let (left, _) = split(ctrl, b);
    if b > 0.0 { split(&left, a / b).1 } else { left }
}

/// Decomposes a B-spline of degree `p` with knots `u` and homogeneous
/// controls into its Bezier spans (Piegl and Tiller, algorithm A5.6).
pub(super) fn decompose(p: usize, u: &[f64], pw: &[[f64; 4]]) -> Vec<Span> {
    let m = u.len() - 1;
    let mut spans = Vec::new();
    let mut a = p;
    let mut b = p + 1;
    let mut q: Vec<[f64; 4]> = pw[..=p].to_vec();
    let mut next = vec![[0.0; 4]; p + 1];
    let mut alphas = vec![0.0; p.max(1)];
    while b < m {
        let i = b;
        // Knots never decrease: equal is not above.
        while b < m && !above(u[b + 1], u[b]) {
            b += 1;
        }
        let mult = b - i + 1;
        if mult < p {
            let numer = u[b] - u[a];
            for j in (mult + 1..=p).rev() {
                alphas[j - mult - 1] = numer / (u[a + j] - u[a]);
            }
            let r = p - mult;
            for j in 1..=r {
                let save = r - j;
                let s = mult + j;
                for k in (s..=p).rev() {
                    let alpha = alphas[k - s];
                    for c in 0..4 {
                        q[k][c] = alpha * q[k][c] + (1.0 - alpha) * q[k - 1][c];
                    }
                }
                if b < m {
                    next[save] = q[p];
                }
            }
        }
        if u[b] > u[a] {
            spans.push(Span {
                lo: u[a],
                hi: u[b],
                ctrl: q.clone(),
            });
        }
        if b < m {
            for i in p - mult.min(p)..=p {
                next[i] = pw[b - p + i];
            }
            q = next.clone();
            a = b;
            b += 1;
        }
    }
    spans
}

/// How the native parameter follows a piece's Bezier parameter `s`.
#[derive(Clone, Copy, Debug)]
pub(super) enum Map {
    /// `lo + s (hi - lo)`.
    Linear { lo: f64, hi: f64 },
    /// A rational quadratic arc of half angle `half` about angle `mid`:
    /// `tan((t - mid) / 2) = (2 s - 1) tan(half / 2)`, the half-angle form of
    /// the circular arc (Piegl and Tiller, section 7.3).
    Arc { mid: f64, half: f64 },
}

impl Map {
    pub fn native(self, s: f64) -> f64 {
        match self {
            Map::Linear { lo, hi } => lo + s * (hi - lo),
            Map::Arc { mid, half } => mid + 2.0 * ((2.0 * s - 1.0) * (half / 2.0).tan()).atan(),
        }
    }
}

/// A rational Bezier piece of a stroke with its parameter map.
#[derive(Clone, Debug)]
pub(super) struct Piece {
    pub ctrl: Vec<[f64; 4]>,
    pub map: Map,
}

#[derive(Clone, Debug)]
pub(super) enum Geo {
    /// `o + t v`.
    Line { o: V3, v: V3 },
    /// `c + cos t u + sin t w`, with `u` perpendicular to `w`.
    Conic { c: V3, u: V3, w: V3 },
    /// Bezier spans in knot order.
    Spline { spans: Vec<Span> },
}

impl Geo {
    fn span(spans: &[Span], t: f64) -> (&Span, f64) {
        let i = spans.partition_point(|s| s.hi < t).min(spans.len() - 1);
        let span = &spans[i];
        (span, (t - span.lo) / (span.hi - span.lo))
    }

    /// Point and derivative in the native parameter.
    pub fn jet(&self, t: f64) -> (V3, V3) {
        match self {
            Geo::Line { o, v } => (add(*o, scale(*v, t)), *v),
            Geo::Conic { c, u, w } => {
                let (sin, cos) = t.sin_cos();
                (
                    add(*c, add(scale(*u, cos), scale(*w, sin))),
                    add(scale(*u, -sin), scale(*w, cos)),
                )
            }
            Geo::Spline { spans } => {
                let (span, s) = Self::span(spans, t);
                let (h, dh) = casteljau(&span.ctrl, s);
                let point = [h[0] / h[3], h[1] / h[3], h[2] / h[3]];
                let speed = 1.0 / (span.hi - span.lo);
                let derivative =
                    std::array::from_fn(|k| (dh[k] * h[3] - h[k] * dh[3]) / (h[3] * h[3]) * speed);
                (point, derivative)
            }
        }
    }

    pub fn point(&self, t: f64) -> V3 {
        self.jet(t).0
    }

    /// The parameter of a point on the carrier, or of its nearest point.
    pub fn param(&self, p: V3) -> f64 {
        match self {
            Geo::Line { o, v } => dot(sub(p, *o), *v) / dot(*v, *v),
            Geo::Conic { c, u, w } => {
                let r = sub(p, *c);
                (dot(r, *w) / dot(*w, *w)).atan2(dot(r, *u) / dot(*u, *u))
            }
            Geo::Spline { spans } => {
                // The nearest of each span's ends and middle, then
                // Gauss-Newton steps on |P(t) - p|^2 (Piegl and Tiller,
                // section 6.1, point inversion), which converge
                // quadratically for a point on the curve.
                let mut best = (f64::INFINITY, spans[0].lo);
                for span in spans {
                    for t in [span.lo, (span.lo + span.hi) / 2.0, span.hi] {
                        let d = distance(self.point(t), p);
                        if d < best.0 {
                            best = (d, t);
                        }
                    }
                }
                let (lo, hi) = (spans[0].lo, spans[spans.len() - 1].hi);
                let mut t = best.1;
                for _ in 0..super::cross::SOLVE_STEPS {
                    let (q, dq) = self.jet(t);
                    let f = dot(sub(q, p), dq);
                    let df = dot(dq, dq);
                    if df == 0.0 || !df.is_finite() {
                        break;
                    }
                    let next = (t - f / df).clamp(lo, hi);
                    if (next - t).abs() <= f64::EPSILON * hi.abs().max(lo.abs()) {
                        t = next;
                        break;
                    }
                    t = next;
                }
                t
            }
        }
    }

    /// The parameter domain of a closed carrier's whole loop.
    pub fn whole(&self) -> [f64; 2] {
        match self {
            Geo::Spline { spans } => [spans[0].lo, spans[spans.len() - 1].hi],
            _ => [0.0, TAU],
        }
    }

    /// Lower and upper bounds of `P(t) . e` over `[t0, t1]`: exact extremes
    /// for lines and conics, the control hull for splines (Farin,
    /// section 15.3).
    pub fn range(&self, t0: f64, t1: f64, e: V3) -> (f64, f64) {
        match self {
            Geo::Line { .. } => {
                let (a, b) = (dot(self.point(t0), e), dot(self.point(t1), e));
                (a.min(b), a.max(b))
            }
            Geo::Conic { c, u, w } => {
                let (a, b) = (dot(self.point(t0), e), dot(self.point(t1), e));
                let (mut lo, mut hi) = (a.min(b), a.max(b));
                let (ue, we) = (dot(*u, e), dot(*w, e));
                let amplitude = (ue * ue + we * we).sqrt();
                let peak = we.atan2(ue);
                let middle = dot(*c, e);
                for (angle, value) in [(peak, middle + amplitude), (peak + PI, middle - amplitude)]
                {
                    let k = ((t0 - angle) / TAU).ceil();
                    if angle + k * TAU <= t1 {
                        lo = lo.min(value);
                        hi = hi.max(value);
                    }
                }
                (lo, hi)
            }
            Geo::Spline { spans } => {
                let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
                for span in spans.iter().filter(|s| s.hi >= t0 && s.lo <= t1) {
                    for h in &span.ctrl {
                        let x = (h[0] * e[0] + h[1] * e[1] + h[2] * e[2]) / h[3];
                        lo = lo.min(x);
                        hi = hi.max(x);
                    }
                }
                (lo, hi)
            }
        }
    }

    /// Rational Bezier pieces covering `[t0, t1]`. A conic is cut into arcs
    /// of at most a quarter turn, each exactly a rational quadratic.
    pub fn pieces(&self, t0: f64, t1: f64) -> Vec<Piece> {
        let homogeneous = |p: V3| [p[0], p[1], p[2], 1.0];
        match self {
            Geo::Line { .. } => vec![Piece {
                ctrl: vec![homogeneous(self.point(t0)), homogeneous(self.point(t1))],
                map: Map::Linear { lo: t0, hi: t1 },
            }],
            Geo::Conic { c, u, w } => {
                let count = ((t1 - t0) / FRAC_PI_2).ceil().max(1.0) as usize;
                let step = (t1 - t0) / count as f64;
                (0..count)
                    .map(|k| {
                        let a = t0 + k as f64 * step;
                        let b = if k + 1 == count { t1 } else { a + step };
                        let (mid, half) = ((a + b) / 2.0, (b - a) / 2.0);
                        let weight = half.cos();
                        let (sin, cos) = mid.sin_cos();
                        let apex =
                            add(*c, scale(add(scale(*u, cos), scale(*w, sin)), 1.0 / weight));
                        Piece {
                            ctrl: vec![
                                homogeneous(self.point(a)),
                                [apex[0] * weight, apex[1] * weight, apex[2] * weight, weight],
                                homogeneous(self.point(b)),
                            ],
                            map: Map::Arc { mid, half },
                        }
                    })
                    .collect()
            }
            Geo::Spline { spans } => spans
                .iter()
                .filter(|s| s.hi > t0 && s.lo < t1)
                .map(|span| {
                    let width = span.hi - span.lo;
                    let a = ((t0 - span.lo) / width).max(0.0);
                    let b = ((t1 - span.lo) / width).min(1.0);
                    Piece {
                        ctrl: restrict(&span.ctrl, a, b),
                        map: Map::Linear {
                            lo: span.lo + a * width,
                            hi: span.lo + b * width,
                        },
                    }
                })
                .collect(),
        }
    }
}

/// An orthonormal binary64 frame of the sheet: sheet coordinates are
/// `(p . e1, p . e2)`.
#[derive(Clone, Copy, Debug)]
pub(super) struct Sheet {
    pub e1: V3,
    pub e2: V3,
}

impl Sheet {
    pub fn new(d: V3) -> Self {
        // The axis least aligned with d gives a well-conditioned first
        // sheet direction.
        let k = (0..3)
            .min_by(|&i, &j| d[i].abs().total_cmp(&d[j].abs()))
            .expect("three axes");
        let mut axis = [0.0; 3];
        axis[k] = 1.0;
        let e1 = unit(super::super::float::cross(axis, d));
        let e2 = super::super::float::cross(d, e1);
        Self { e1, e2 }
    }

    pub fn project(&self, p: V3) -> [f64; 2] {
        [dot(p, self.e1), dot(p, self.e2)]
    }

    /// Sheet bounds `[min x, min y, max x, max y]` of a span.
    pub fn bounds(&self, geo: &Geo, t0: f64, t1: f64) -> [f64; 4] {
        let (x0, x1) = geo.range(t0, t1, self.e1);
        let (y0, y1) = geo.range(t0, t1, self.e2);
        [x0, y0, x1, y1]
    }
}

/// A box in model space, `[min, max]`.
pub(super) type Box3 = [V3; 2];

pub(super) fn box3(geo: &Geo, t0: f64, t1: f64) -> Box3 {
    let mut b = [[0.0; 3]; 2];
    for k in 0..3 {
        let mut e = [0.0; 3];
        e[k] = 1.0;
        let (lo, hi) = geo.range(t0, t1, e);
        b[0][k] = lo;
        b[1][k] = hi;
    }
    b
}

pub(super) fn union(a: Box3, b: Box3) -> Box3 {
    [
        std::array::from_fn(|k| a[0][k].min(b[0][k])),
        std::array::from_fn(|k| a[1][k].max(b[1][k])),
    ]
}

pub(super) fn point_box(p: V3) -> Box3 {
    [p, p]
}

/// Whether the closed ray `p + t d`, `t >= 0`, meets the box widened by
/// `margin` (the slab test of Kay and Kajiya, "Ray tracing complex scenes",
/// SIGGRAPH 1986).
pub(super) fn ray_meets(b: &Box3, p: V3, d: V3, margin: f64) -> bool {
    let (mut near, mut far) = (0.0_f64, f64::INFINITY);
    for k in 0..3 {
        let (lo, hi) = (b[0][k] - margin, b[1][k] + margin);
        if d[k] == 0.0 {
            if p[k] < lo || p[k] > hi {
                return false;
            }
            continue;
        }
        let (a, c) = ((lo - p[k]) / d[k], (hi - p[k]) / d[k]);
        near = near.max(a.min(c));
        far = far.min(a.max(c));
        if near > far {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decomposed_spans_reproduce_a_cubic_with_an_inner_knot() {
        // A uniform clamped cubic with one inner knot: two Bezier spans
        // whose ends meet and whose joint matches de Boor's value.
        let u = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
        let pw: Vec<[f64; 4]> = [[0.0, 0.0], [1.0, 2.0], [2.0, -1.0], [3.0, 1.0], [4.0, 0.0]]
            .iter()
            .map(|p| [p[0], p[1], 0.0, 1.0])
            .collect();
        let spans = decompose(3, &u, &pw);
        assert_eq!(spans.len(), 2);
        assert_eq!(
            (spans[0].lo, spans[0].hi, spans[1].lo, spans[1].hi),
            (0.0, 0.5, 0.5, 1.0)
        );
        assert_eq!(spans[0].ctrl[3], spans[1].ctrl[0]);
        assert_eq!(spans[0].ctrl[0], pw[0]);
        assert_eq!(spans[1].ctrl[3], pw[4]);
        // The joint: N_{i,3}(0.5) for this knot vector is (0, 1/4, 1/2, 1/4, 0).
        let joint = spans[0].ctrl[3];
        assert!((joint[0] - 2.0).abs() < 1e-15 && (joint[1] - 0.25).abs() < 1e-15);
    }

    #[test]
    fn conic_ranges_and_quarter_arcs_bound_and_follow_the_circle() {
        let geo = Geo::Conic {
            c: [1.0, 2.0, 0.0],
            u: [2.0, 0.0, 0.0],
            w: [0.0, 2.0, 0.0],
        };
        let (lo, hi) = geo.range(0.0, PI, [1.0, 0.0, 0.0]);
        assert!((lo + 1.0).abs() < 1e-15 && (hi - 3.0).abs() < 1e-15);
        let (lo, hi) = geo.range(0.0, PI, [0.0, 1.0, 0.0]);
        assert!((lo - 2.0).abs() < 1e-15 && (hi - 4.0).abs() < 1e-15);
        for piece in geo.pieces(0.3, 2.9) {
            for s in [0.0, 0.25, 0.5, 1.0] {
                let (h, _) = casteljau(&piece.ctrl, s);
                let p = [h[0] / h[3], h[1] / h[3], h[2] / h[3]];
                let q = geo.point(piece.map.native(s));
                assert!(distance(p, q) < 1e-14, "{s}: {p:?} {q:?}");
            }
        }
    }

    #[test]
    fn rays_meet_boxes_only_ahead() {
        let b = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]];
        assert!(ray_meets(&b, [0.5, 0.5, -1.0], [0.0, 0.0, 1.0], 0.0));
        assert!(!ray_meets(&b, [0.5, 0.5, 2.0], [0.0, 0.0, 1.0], 0.0));
        assert!(!ray_meets(&b, [2.0, 0.5, -1.0], [0.0, 0.0, 1.0], 0.5));
        assert!(ray_meets(&b, [1.4, 0.5, -1.0], [0.0, 0.0, 1.0], 0.5));
    }
}
