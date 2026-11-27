use std::cmp::Ordering;

use catallax_shapes::{BSplineCurve, BSplineSurface};

use super::{DisplaySpline, Xyz, evaluate::*};

// Reuse basis storage across all evaluations of a carrier. Only controls in
// the active spans contribute to the evaluated point.
struct Axis<'a> {
    knots: &'a [f64],
    degree: usize,
    count: usize,
    basis: Vec<f64>,
}

impl<'a> Axis<'a> {
    fn new(knots: &'a [f64], degree: usize, count: usize) -> Self {
        Self {
            knots,
            degree,
            count,
            basis: vec![0.0; degree + 1],
        }
    }

    fn at(&mut self, t: f64) -> (usize, &[f64]) {
        let p = self.degree;
        let span = (self.knots.partition_point(|&k| k <= t) - 1).min(self.count - 1);
        self.basis[0] = 1.0;
        for j in 1..=p {
            let mut tail = 0.0;
            for r in 0..j {
                let lo = self.knots[span + 1 + r - j];
                let hi = self.knots[span + 1 + r];
                let q = self.basis[r] / (hi - lo);
                self.basis[r] = tail + (hi - t) * q;
                tail = (t - lo) * q;
            }
            self.basis[j] = tail;
        }
        (span - p, &self.basis)
    }
}

pub(super) fn net(s: &BSplineSurface, origin: Xyz) -> Option<DisplaySpline> {
    let [nu, nv] = s.control_count();
    let mut points = Vec::with_capacity(nu * nv);
    let mut weights = Vec::with_capacity(nu * nv);
    for v in 0..nv {
        for u in 0..nu {
            let i = u * nv + v;
            points.push(sub(point(&s.controls()[i])?, origin));
            weights.push(positive(scalar(&s.weights()[i]))?);
        }
    }
    let [u, v] = s.knots();
    Some(DisplaySpline {
        degree: s.degrees(),
        knots: [
            knots(u.iter().map(|k| k.to_f64_nearest()))?,
            knots(v.iter().map(|k| k.to_f64_nearest()))?,
        ],
        size: [nu, nv],
        points,
        weights: (!weights
            .iter()
            .all(|w| w.partial_cmp(&weights[0]) == Some(Ordering::Equal)))
        .then_some(weights),
    })
}

fn knots(values: impl Iterator<Item = f64>) -> Option<Vec<f64>> {
    values.map(finite).collect()
}

pub(super) struct Patch<'a> {
    net: &'a DisplaySpline,
    axes: [Axis<'a>; 2],
    pub(super) bounds: [[f64; 2]; 2],
}

impl<'a> Patch<'a> {
    pub(super) fn new(net: &'a DisplaySpline) -> Option<Self> {
        let bounds =
            std::array::from_fn(|i| [net.knots[i][net.degree[i]], net.knots[i][net.size[i]]]);
        if bounds.iter().any(|[lo, hi]| lo >= hi) {
            return None;
        }
        Some(Self {
            net,
            axes: std::array::from_fn(|i| Axis::new(&net.knots[i], net.degree[i], net.size[i])),
            bounds,
        })
    }

    pub(super) fn closed(&self, axis: usize) -> bool {
        let n = self.net.size[axis];
        let degree = self.net.degree[axis];
        let [lo, hi] = self.bounds[axis].map(f64::to_bits);
        let knots = &self.net.knots[axis];
        if !knots[..degree + 1].iter().all(|x| x.to_bits() == lo)
            || !knots[n..].iter().all(|x| x.to_bits() == hi)
        {
            return false;
        }
        let stride = if axis == 0 { 1 } else { self.net.size[0] };
        let next = if axis == 0 { self.net.size[0] } else { 1 };
        (0..self.net.size[1 - axis]).all(|i| {
            let a = i * next;
            let b = a + (n - 1) * stride;
            self.net.points[a].map(f64::to_bits) == self.net.points[b].map(f64::to_bits)
                && self
                    .net
                    .weights
                    .as_ref()
                    .is_none_or(|w| w[a].to_bits() == w[b].to_bits())
        })
    }

    pub(super) fn at(&mut self, uv: [f64; 2]) -> Xyz {
        let [u, v] = &mut self.axes;
        let (u0, bu) = u.at(uv[0].clamp(self.bounds[0][0], self.bounds[0][1]));
        let (v0, bv) = v.at(uv[1].clamp(self.bounds[1][0], self.bounds[1][1]));
        let mut p = [0.0; 3];
        let mut sum = 0.0;
        for (j, &v) in bv.iter().enumerate() {
            for (i, &u) in bu.iter().enumerate() {
                let k = (v0 + j) * self.net.size[0] + u0 + i;
                let w = u * v * self.net.weights.as_ref().map_or(1.0, |ws| ws[k]);
                p = add(p, scale(self.net.points[k], w));
                sum += w;
            }
        }
        scale(p, 1.0 / sum)
    }

    pub(super) fn jet(&mut self, uv: [f64; 2]) -> (Xyz, [Xyz; 2]) {
        let p = self.at(uv);
        let tangents = std::array::from_fn(|i| {
            let [lo, hi] = self.bounds[i];
            let h = (hi - lo) * 0.000001;
            let mut a = uv;
            let mut b = uv;
            a[i] = lo.max(uv[i] - h);
            b[i] = hi.min(uv[i] + h);
            scale(sub(self.at(b), self.at(a)), 1.0 / (b[i] - a[i]))
        });
        (p, tangents)
    }

    pub(super) fn project(
        &mut self,
        p: Xyz,
        previous: Option<[f64; 2]>,
        periods: [Option<f64>; 2],
    ) -> [f64; 2] {
        let mut best = previous.map(|uv| self.fit(p, uv));
        if best.is_some_and(|(_, error)| error <= 1e-20) {
            return best.unwrap().0;
        }
        // A control-net seed finds narrow knot spans that a uniform UV search
        // can miss. A preceding trim sample also seeds its own chart branch.
        let n = self.net.size[0];
        let i = self
            .net
            .points
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                dot(sub(**a, p), sub(**a, p)).total_cmp(&dot(sub(**b, p), sub(**b, p)))
            })
            .unwrap()
            .0;
        let seed = [i % n, i / n];
        let uv = std::array::from_fn(|axis| {
            let d = self.net.degree[axis];
            self.net.knots[axis][seed[axis] + 1..=seed[axis] + d]
                .iter()
                .sum::<f64>()
                / d as f64
        });
        let estimate = self.fit(p, uv);
        if best.is_none_or(|(_, error)| estimate.1 < error) {
            best = Some(estimate);
        }
        if let Some(previous) = previous {
            for (axis, period) in periods.into_iter().enumerate() {
                if let Some(period) = period
                    && best.unwrap().1 > 1e-20
                {
                    let [lo, hi] = self.bounds[axis];
                    let mut seed = previous;
                    seed[axis] = if previous[axis] < lo + period / 2.0 {
                        hi
                    } else {
                        lo
                    };
                    let estimate = self.fit(p, seed);
                    if estimate.1 < best.unwrap().1 {
                        best = Some(estimate);
                    }
                }
            }
        }
        best.unwrap().0
    }

    fn fit(&mut self, p: Xyz, mut uv: [f64; 2]) -> ([f64; 2], f64) {
        for (v, [lo, hi]) in uv.iter_mut().zip(self.bounds) {
            *v = v.clamp(lo, hi);
        }
        let mut residual = sub(p, self.at(uv));
        let mut error = dot(residual, residual);
        for _ in 0..24 {
            if error <= 1e-24 {
                break;
            }
            let (_, [du, dv]) = self.jet(uv);
            let Some(delta) = solve_gram(du, dv, [dot(du, residual), dot(dv, residual)]) else {
                break;
            };
            let mut improved = false;
            let mut step = 1.0;
            for _ in 0..8 {
                let next = std::array::from_fn(|i| {
                    (uv[i] + step * delta[i]).clamp(self.bounds[i][0], self.bounds[i][1])
                });
                let r = sub(p, self.at(next));
                let e = dot(r, r);
                if e < error {
                    (uv, residual, error) = (next, r, e);
                    improved = true;
                    break;
                }
                step *= 0.5;
            }
            if !improved {
                break;
            }
        }
        (uv, error)
    }
}

pub(super) struct SplineEdge {
    degree: usize,
    knots: Vec<f64>,
    points: Vec<Xyz>,
    weights: Vec<f64>,
}

impl SplineEdge {
    pub(super) fn new(curve: &BSplineCurve, origin: Xyz) -> Option<Self> {
        Some(Self {
            degree: curve.degree(),
            knots: knots(curve.knots().iter().map(|k| k.to_f64_nearest()))?,
            points: curve
                .controls()
                .map(|p| point(p).map(|p| sub(p, origin)))
                .collect::<Option<_>>()?,
            weights: curve
                .weights()
                .map(|w| positive(scalar(w)))
                .collect::<Option<_>>()?,
        })
    }

    pub(super) fn evaluator(&self) -> impl FnMut(f64) -> Xyz + '_ {
        let count = self.points.len();
        let [lo, hi] = [self.knots[self.degree], self.knots[count]];
        let mut axis = Axis::new(&self.knots, self.degree, count);
        move |t| {
            let (offset, basis) = axis.at(t.clamp(lo, hi));
            let mut point = [0.0; 3];
            let mut weight = 0.0;
            for (i, &b) in basis.iter().enumerate() {
                let w = b * self.weights[offset + i];
                point = add(point, scale(self.points[offset + i], w));
                weight += w;
            }
            scale(point, 1.0 / weight)
        }
    }

    pub(super) fn breaks(&self) -> &[f64] {
        &self.knots[self.degree..=self.points.len()]
    }

    pub(super) fn samples(
        &self,
        ends: [Option<Xyz>; 2],
        tolerance: f64,
        max: usize,
    ) -> Option<Vec<Xyz>> {
        let count = self.points.len();
        let [lo, hi] = [self.knots[self.degree], self.knots[count]];
        let mut at = self.evaluator();
        let mut domain = [lo, hi];
        for (i, end) in ends.into_iter().enumerate() {
            if let Some(end) = end
                && norm(sub(end, at(domain[i]))) > tolerance * 0.001
            {
                domain[i] = nearest_parameter(&mut at, &self.knots[self.degree..=count], end);
            }
        }
        if domain[1] <= domain[0] {
            return None;
        }
        let first = ends[0].unwrap_or_else(|| at(domain[0]));
        let last = ends[1].unwrap_or_else(|| at(domain[1]));
        let mut result = vec![first];
        let mut a = domain[0];
        let mut start = first;
        for b in self.knots[self.degree..count]
            .iter()
            .copied()
            .filter(|&k| k > domain[0] && k < domain[1])
            .chain([domain[1]])
        {
            if b.partial_cmp(&a) == Some(Ordering::Equal) {
                continue;
            }
            let end = if b.partial_cmp(&domain[1]) == Some(Ordering::Equal) {
                last
            } else {
                at(b)
            };
            super::edges::sample_span(&mut at, [a, b], [start, end], tolerance, &mut |_, p| {
                (result.len() < max).then(|| result.push(p))
            })?;
            (a, start) = (b, end);
        }
        Some(result)
    }
}

fn nearest_parameter<F: FnMut(f64) -> Xyz>(at: &mut F, knots: &[f64], point: Xyz) -> f64 {
    let error = |at: &mut F, t| {
        let d = sub(at(t), point);
        dot(d, d)
    };
    let mut best = (knots[0], f64::INFINITY);
    for pair in knots.windows(2) {
        let [lo, hi] = [pair[0], pair[1]];
        if lo.partial_cmp(&hi) == Some(Ordering::Equal) {
            continue;
        }
        // Bracket a local minimum before narrowing, even if the whole knot
        // span is not unimodal in distance from the point.
        let width = (hi - lo) / 16.0;
        let (index, _) = (0..=16)
            .map(|i| {
                (
                    i,
                    error(at, if i == 16 { hi } else { lo + width * i as f64 }),
                )
            })
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap();
        let mut a = lo + width * (index - 1).max(0) as f64;
        let mut b = hi.min(lo + width * (index + 1) as f64);
        for _ in 0..48 {
            let l = a + (b - a) / 3.0;
            let r = b - (b - a) / 3.0;
            if error(at, l) <= error(at, r) {
                b = r;
            } else {
                a = l;
            }
        }
        let t = a.midpoint(b);
        let e = error(at, t);
        if e < best.1 {
            best = (t, e);
        }
    }
    best.0
}
