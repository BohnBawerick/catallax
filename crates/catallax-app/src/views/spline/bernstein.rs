//! Bernstein polynomials in binary64, and the Bezier patches of a B-spline
//! carrier with the Bernstein form of their contour function.
//!
//! On a rational patch `S = A / w`, with `A` the weighted controls' part and
//! `w` the weights' part, the normal is
//!
//! ```text
//! S_u x S_v = (w A_u x A_v - w_v A_u x A - w_u A x A_v) / w^3
//! ```
//!
//! (Farin, "Curves and Surfaces for CAGD", 5th ed., section 13.1, the
//! derivative of a rational patch, crossed out; `A x A` vanishes). With
//! positive weights the contour function `N . d` has the sign of the
//! polynomial `K . d`, the bracket dotted with the view. Its Bernstein
//! coefficients follow from products of Bernstein polynomials (Farouki and
//! Rajan, "Algorithms for polynomials in Bernstein form", CAGD 5, 1988,
//! section 2.4), and they bound it on the patch by the convex hull property
//! (Farin, section 15.3). SMLib's silhouette tracer asks the same question of
//! its global solver through a normal cone test.

use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::BSplineSurface;

use super::super::float::double;

fn binomial(n: usize, k: usize) -> f64 {
    let mut value = 1.0;
    for i in 0..k {
        value = value * (n - i) as f64 / (i + 1) as f64;
    }
    value
}

/// A univariate Bernstein polynomial over `[0, 1]`.
#[derive(Clone, Debug)]
pub(super) struct Poly1(pub Vec<f64>);

impl Poly1 {
    /// De Casteljau's algorithm (Farin, section 3.1); also extrapolates.
    pub fn eval(&self, t: f64) -> f64 {
        let mut level = self.0.clone();
        for round in 1..level.len() {
            for i in 0..level.len() - round {
                level[i] = level[i] + t * (level[i + 1] - level[i]);
            }
        }
        level[0]
    }

    /// The two halves at `t` by de Casteljau's algorithm, each over `[0, 1]`.
    pub fn split(&self, t: f64) -> (Self, Self) {
        let n = self.0.len();
        let mut level = self.0.clone();
        let mut left = Vec::with_capacity(n);
        let mut right = vec![0.0; n];
        left.push(level[0]);
        right[n - 1] = level[n - 1];
        for round in 1..n {
            for i in 0..n - round {
                level[i] = level[i] + t * (level[i + 1] - level[i]);
            }
            left.push(level[0]);
            right[n - 1 - round] = level[n - 1 - round];
        }
        (Self(left), Self(right))
    }

    /// The polynomial over `[a, b]`, reparametrised to `[0, 1]`.
    pub fn restrict(&self, a: f64, b: f64) -> Self {
        let (left, _) = self.split(b);
        if b == 0.0 {
            return left;
        }
        left.split(a / b).1
    }
}

/// A bivariate Bernstein polynomial of degrees `(m, n)` over `[0, 1]^2`,
/// coefficients in `s`-major order.
#[derive(Clone, Debug)]
pub(super) struct Poly2 {
    pub m: usize,
    pub n: usize,
    pub c: Vec<f64>,
}

impl Poly2 {
    fn zero(m: usize, n: usize) -> Self {
        Self {
            m,
            n,
            c: vec![0.0; (m + 1) * (n + 1)],
        }
    }

    fn at(&self, i: usize, j: usize) -> f64 {
        self.c[i * (self.n + 1) + j]
    }

    /// The product, by Farouki and Rajan's formula for Bernstein
    /// coefficients of a product.
    pub fn times(&self, other: &Self) -> Self {
        let mut out = Self::zero(self.m + other.m, self.n + other.n);
        for i in 0..=self.m {
            for k in 0..=other.m {
                let su = binomial(self.m, i) * binomial(other.m, k) / binomial(out.m, i + k);
                for j in 0..=self.n {
                    for l in 0..=other.n {
                        let sv =
                            binomial(self.n, j) * binomial(other.n, l) / binomial(out.n, j + l);
                        out.c[(i + k) * (out.n + 1) + j + l] +=
                            su * sv * self.at(i, j) * other.at(k, l);
                    }
                }
            }
        }
        out
    }

    fn plus(mut self, other: &Self, scale: f64) -> Self {
        debug_assert_eq!((self.m, self.n), (other.m, other.n));
        for (a, b) in self.c.iter_mut().zip(&other.c) {
            *a += scale * b;
        }
        self
    }

    /// The restriction to `s = at`, a polynomial in `t`.
    pub fn at_s(&self, at: f64) -> Poly1 {
        Poly1(
            (0..=self.n)
                .map(|j| Poly1((0..=self.m).map(|i| self.at(i, j)).collect()).eval(at))
                .collect(),
        )
    }

    /// The restriction to `t = at`, a polynomial in `s`.
    pub fn at_t(&self, at: f64) -> Poly1 {
        Poly1(
            (0..=self.m)
                .map(|i| Poly1(self.c[i * (self.n + 1)..(i + 1) * (self.n + 1)].to_vec()).eval(at))
                .collect(),
        )
    }

    pub fn eval(&self, s: f64, t: f64) -> f64 {
        self.at_t(t).eval(s)
    }

    /// The partial derivative in `s` (degree `m - 1`), or in `t`.
    pub fn derivative(&self, axis: usize) -> Self {
        let (m, n) = (self.m, self.n);
        if axis == 0 {
            let mut out = Self::zero(m.saturating_sub(1), n);
            if m == 0 {
                return out;
            }
            for i in 0..m {
                for j in 0..=n {
                    out.c[i * (n + 1) + j] = m as f64 * (self.at(i + 1, j) - self.at(i, j));
                }
            }
            out
        } else {
            let mut out = Self::zero(m, n.saturating_sub(1));
            if n == 0 {
                return out;
            }
            for i in 0..=m {
                for j in 0..n {
                    out.c[i * n + j] = n as f64 * (self.at(i, j + 1) - self.at(i, j));
                }
            }
            out
        }
    }

    /// The two halves at the middle of one axis.
    pub fn halves(&self, axis: usize) -> (Self, Self) {
        let mut left = Self::zero(self.m, self.n);
        let mut right = Self::zero(self.m, self.n);
        if axis == 0 {
            for j in 0..=self.n {
                let (a, b) = Poly1((0..=self.m).map(|i| self.at(i, j)).collect()).split(0.5);
                for i in 0..=self.m {
                    left.c[i * (self.n + 1) + j] = a.0[i];
                    right.c[i * (self.n + 1) + j] = b.0[i];
                }
            }
        } else {
            for i in 0..=self.m {
                let row = &self.c[i * (self.n + 1)..(i + 1) * (self.n + 1)];
                let (a, b) = Poly1(row.to_vec()).split(0.5);
                left.c[i * (self.n + 1)..(i + 1) * (self.n + 1)].copy_from_slice(&a.0);
                right.c[i * (self.n + 1)..(i + 1) * (self.n + 1)].copy_from_slice(&b.0);
            }
        }
        (left, right)
    }

    /// Whether every coefficient is nonzero with one sign: the polynomial
    /// then has that sign everywhere on the cell.
    pub fn one_sign(&self) -> bool {
        self.c.iter().all(|&x| x > 0.0) || self.c.iter().all(|&x| x < 0.0)
    }
}

/// Homogeneous control `(w x, w y, w z, w)` in the frame of the fit.
pub(super) type Homogeneous = [f64; 4];

/// One nonzero span of the carrier as a rational Bezier patch, by
/// homogeneous coordinate, with its first partial derivatives.
#[derive(Clone, Debug)]
pub(super) struct Patch {
    net: [Poly2; 4],
    net_s: [Poly2; 4],
    net_t: [Poly2; 4],
    /// `K . d` over the patch in its local parameters.
    pub contour: Poly2,
    contour_s: Poly2,
    contour_t: Poly2,
}

/// The carrier in binary64: its Bezier patches and their span breaks.
#[derive(Clone, Debug)]
pub(super) struct Carrier {
    /// Span ends along each axis, ascending, domain ends included.
    pub breaks: [Vec<f64>; 2],
    /// Patches, `u`-major.
    pub patches: Vec<Patch>,
}

/// Position and first derivatives in the carrier's own parameters.
#[derive(Clone, Copy, Debug)]
pub(super) struct Jet2 {
    pub point: [f64; 3],
    pub du: [f64; 3],
    pub dv: [f64; 3],
}

/// The contour function and its gradient in the carrier's parameters. On
/// each patch it is `K . d` divided by the two span widths, which keeps its
/// sign and makes the value continuous where the carrier is `C1`.
#[derive(Clone, Copy, Debug)]
pub(super) struct Level {
    pub value: f64,
    pub grad: [f64; 2],
}

/// Blossom of one B-spline piece by de Boor's algorithm with a different
/// argument at each level (Ramshaw, "Blossoms are polar forms", CAGD 6,
/// 1989): the Bezier controls of span `k` are the blossom at `p - a` copies
/// of the span's start and `a` copies of its end.
fn bezier_row(degree: usize, knots: &[f64], row: &[Homogeneous], k: usize) -> Vec<Homogeneous> {
    let p = degree;
    (0..=p)
        .map(|a| {
            let mut d: Vec<Homogeneous> = row[k - p..=k].to_vec();
            for r in 1..=p {
                let t = if r <= p - a { knots[k] } else { knots[k + 1] };
                for j in (r..=p).rev() {
                    let i = k - p + j;
                    let alpha = (t - knots[i]) / (knots[i + p - r + 1] - knots[i]);
                    for c in 0..4 {
                        d[j][c] = (1.0 - alpha) * d[j - 1][c] + alpha * d[j][c];
                    }
                }
            }
            d[p]
        })
        .collect()
}

fn scalar(net: &[Homogeneous], m: usize, n: usize, c: usize) -> Poly2 {
    Poly2 {
        m,
        n,
        c: net.iter().map(|h| h[c]).collect(),
    }
}

/// `a f + b g` for polynomials of equal degrees.
fn combine(a: f64, f: &Poly2, b: f64, g: &Poly2) -> Poly2 {
    Poly2 {
        m: f.m,
        n: f.n,
        c: f.c.iter().zip(&g.c).map(|(x, y)| a * x + b * y).collect(),
    }
}

/// `K . d` of one patch from its homogeneous coordinate polynomials.
fn contour(net: &[Poly2; 4], d: [f64; 3]) -> Poly2 {
    let a = [net[0].clone(), net[1].clone(), net[2].clone()];
    let w = &net[3];
    let a_s = a.clone().map(|x| x.derivative(0));
    let a_t = a.clone().map(|x| x.derivative(1));
    let (w_s, w_t) = (w.derivative(0), w.derivative(1));
    // X = A_t x d and Y = A x d, componentwise in the constant d.
    let cross = |v: &[Poly2; 3]| -> [Poly2; 3] {
        [
            combine(d[2], &v[1], -d[1], &v[2]),
            combine(d[0], &v[2], -d[2], &v[0]),
            combine(d[1], &v[0], -d[0], &v[1]),
        ]
    };
    let x = cross(&a_t);
    let y = cross(&a);
    let dot = |f: &[Poly2; 3], g: &[Poly2; 3]| {
        f[0].times(&g[0])
            .plus(&f[1].times(&g[1]), 1.0)
            .plus(&f[2].times(&g[2]), 1.0)
    };
    // w (A_s x A_t).d = w A_s.(A_t x d); (A_s x A).d = A_s.(A x d);
    // (A x A_t).d = A.(A_t x d).
    let first = w.times(&dot(&a_s, &x));
    let second = w_t.times(&dot(&a_s, &y));
    let third = w_s.times(&dot(&a, &x));
    first.plus(&second, -1.0).plus(&third, -1.0)
}

impl Carrier {
    /// The carrier's spans as Bezier patches in a frame at `origin`, with the
    /// Bernstein form of `K . d` for the unit view `d`. Returns the work
    /// spent: one unit per coefficient formed.
    pub fn new(
        surface: &BSplineSurface,
        origin: [&Rational; 3],
        d: [f64; 3],
        floor: &Rational,
        dec: &mut Decider,
    ) -> Result<(Self, u64), Refusal> {
        let [p, q] = surface.degrees();
        let [nu, nv] = surface.control_count();
        let mut work = 0_u64;
        let mut homogeneous = Vec::with_capacity(nu * nv);
        for (point, weight) in surface.controls().iter().zip(surface.weights()) {
            let local = [&point.x, &point.y, &point.z]
                .into_iter()
                .zip(origin)
                .map(|(x, o)| double(&(&(x - &Real::from(o.clone())) * weight), floor, dec))
                .collect::<Result<Vec<_>, _>>()?;
            homogeneous.push([local[0], local[1], local[2], double(weight, floor, dec)?]);
            work += 4;
        }
        let knots = surface.knots();
        let mut spans: [Vec<usize>; 2] = [Vec::new(), Vec::new()];
        let mut breaks: [Vec<f64>; 2] = [Vec::new(), Vec::new()];
        let float_knots: [Vec<f64>; 2] =
            std::array::from_fn(|axis| knots[axis].iter().map(Rational::to_f64_nearest).collect());
        for axis in 0..2 {
            let (degree, count) = ([p, q][axis], [nu, nv][axis]);
            for k in degree..count {
                if knots[axis][k] < knots[axis][k + 1] {
                    let (a, b) = (float_knots[axis][k], float_knots[axis][k + 1]);
                    if a >= b {
                        return Err(Refusal::new(
                            "app.contour-fit",
                            "two spline knots are closer than binary64 separates",
                        ));
                    }
                    if breaks[axis].is_empty() {
                        breaks[axis].push(a);
                    }
                    breaks[axis].push(b);
                    spans[axis].push(k);
                }
            }
        }
        let mut patches = Vec::with_capacity(spans[0].len() * spans[1].len());
        for &k in &spans[0] {
            // Bezier rows in u for every control column.
            let columns: Vec<Vec<Homogeneous>> = (0..nv)
                .map(|j| {
                    let row: Vec<Homogeneous> = (0..nu).map(|i| homogeneous[i * nv + j]).collect();
                    bezier_row(p, &float_knots[0], &row, k)
                })
                .collect();
            for &l in &spans[1] {
                dec.check_deadline()?;
                let mut net = Vec::with_capacity((p + 1) * (q + 1));
                for a in 0..=p {
                    let row: Vec<Homogeneous> = (0..nv).map(|j| columns[j][a]).collect();
                    net.extend(bezier_row(q, &float_knots[1], &row, l));
                }
                let net: [Poly2; 4] = std::array::from_fn(|c| scalar(&net, p, q, c));
                let contour = contour(&net, d);
                work += (4 * net[0].c.len() + contour.c.len()) as u64;
                patches.push(Patch {
                    net_s: net.clone().map(|x| x.derivative(0)),
                    net_t: net.clone().map(|x| x.derivative(1)),
                    net,
                    contour_s: contour.derivative(0),
                    contour_t: contour.derivative(1),
                    contour,
                });
            }
        }
        Ok((Self { breaks, patches }, work))
    }

    pub fn spans(&self, axis: usize) -> usize {
        self.breaks[axis].len() - 1
    }

    pub fn patch(&self, k: usize, l: usize) -> &Patch {
        &self.patches[k * self.spans(1) + l]
    }

    /// The span holding `x` on one axis, the nearest end span outside the
    /// domain, and the local parameter there.
    pub fn locate(&self, axis: usize, x: f64) -> (usize, f64) {
        let b = &self.breaks[axis];
        let last = b.len() - 2;
        let k = b[1..=last].partition_point(|&e| e <= x);
        (k, (x - b[k]) / (b[k + 1] - b[k]))
    }

    pub fn width(&self, axis: usize, k: usize) -> f64 {
        self.breaks[axis][k + 1] - self.breaks[axis][k]
    }

    /// Position and first derivatives at `(u, v)`.
    pub fn jet(&self, u: f64, v: f64) -> Jet2 {
        let (k, s) = self.locate(0, u);
        let (l, t) = self.locate(1, v);
        let patch = self.patch(k, l);
        let at =
            |polys: &[Poly2; 4]| -> Homogeneous { std::array::from_fn(|c| polys[c].eval(s, t)) };
        let (h, hs, ht) = (at(&patch.net), at(&patch.net_s), at(&patch.net_t));
        let point = [h[0] / h[3], h[1] / h[3], h[2] / h[3]];
        let tangent = |dh: Homogeneous, width: f64| -> [f64; 3] {
            std::array::from_fn(|c| (dh[c] - point[c] * dh[3]) / (h[3] * width))
        };
        Jet2 {
            point,
            du: tangent(hs, self.width(0, k)),
            dv: tangent(ht, self.width(1, l)),
        }
    }

    /// The contour function and its gradient at `(u, v)`.
    pub fn level(&self, u: f64, v: f64) -> Level {
        let (k, s) = self.locate(0, u);
        let (l, t) = self.locate(1, v);
        let patch = self.patch(k, l);
        let (wu, wv) = (self.width(0, k), self.width(1, l));
        let scale = wu * wv;
        Level {
            value: patch.contour.eval(s, t) / scale,
            grad: [
                patch.contour_s.eval(s, t) / (scale * wu),
                patch.contour_t.eval(s, t) / (scale * wv),
            ],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn products_and_splits_agree_with_pointwise_values() {
        let f = Poly2 {
            m: 2,
            n: 1,
            c: vec![1.0, -2.0, 0.5, 3.0, -1.0, 2.0],
        };
        let g = Poly2 {
            m: 1,
            n: 2,
            c: vec![0.25, 1.0, -3.0, 2.0, 0.0, 1.5],
        };
        let fg = f.times(&g);
        let (left, right) = f.halves(0);
        let (low, high) = f.halves(1);
        for (s, t) in [(0.1, 0.7), (0.5, 0.5), (0.9, 0.2)] {
            assert!((fg.eval(s, t) - f.eval(s, t) * g.eval(s, t)).abs() < 1e-12);
            let x = if s < 0.5 {
                left.eval(2.0 * s, t)
            } else {
                right.eval(2.0 * s - 1.0, t)
            };
            assert!((x - f.eval(s, t)).abs() < 1e-12);
            let y = if t < 0.5 {
                low.eval(s, 2.0 * t)
            } else {
                high.eval(s, 2.0 * t - 1.0)
            };
            assert!((y - f.eval(s, t)).abs() < 1e-12);
            let h = 1e-6;
            let ds = (f.eval(s + h, t) - f.eval(s - h, t)) / (2.0 * h);
            assert!((f.derivative(0).eval(s, t) - ds).abs() < 1e-8);
        }
        let line = Poly1(vec![1.0, -1.0, 2.0, 0.5]);
        let piece = line.restrict(0.25, 0.75);
        for x in [0.0, 0.3, 1.0] {
            assert!((piece.eval(x) - line.eval(0.25 + 0.5 * x)).abs() < 1e-12);
        }
    }
}
